# coding: utf-8

"""
A new and generalised approach for Trigger-Object matching
"""

from typing import Optional

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

from httcp.util import filter_by_triggers, get_objs_p4, trigger_object_matching_deep

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Muon.pt"    , "Muon.eta"    , "Muon.phi"    , "Muon.mass",
        "Tau.pt",      "Tau.eta",      "Tau.phi",      "Tau.mass",
    },
    #produces={
    #"trigger_ids",
    #},
    exposed=False
)
def match_trigobj(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        etau_indices_pair: ak.Array,
        mutau_indices_pair: ak.Array,
        tautau_indices_pair: ak.Array,
        dotrigobjmatch: Optional[bool] = False,
        **kwargs
) -> ak.Array:
        
    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_names         = trigger_results.x.trigger_names
    trigger_types         = trigger_results.x.trigger_types
    trigger_ids           = trigger_results.x.trigger_ids
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs

    # check etau triggers, filter out the indices based on any etau trigger passed or not and then get the etau pair
    has_e_triggers     = ((trigger_types == "single_e") | (trigger_types == "cross_e_tau"))
    etau_indices_pair  = filter_by_triggers(etau_indices_pair, has_e_triggers)
    etau_pair = ak.concatenate([events.Electron[etau_indices_pair[:,0:1]],
                                events.Tau[etau_indices_pair[:,1:2]]],
                               axis=1)
    # same for mutau pairs
    has_mu_triggers     = ((trigger_types == "single_mu") | (trigger_types == "cross_mu_tau"))
    mutau_indices_pair  = filter_by_triggers(mutau_indices_pair, has_mu_triggers)
    mutau_pair          = ak.concatenate([events.Muon[mutau_indices_pair[:,0:1]], 
                                          events.Tau[mutau_indices_pair[:,1:2]]],
                                         axis=1)
    # same for tautau pairs
    has_tau_triggers     = (trigger_types == "cross_tau_tau")
    tautau_indices_pair  = filter_by_triggers(tautau_indices_pair, has_tau_triggers)
    tautau_pair          = ak.concatenate([events.Tau[tautau_indices_pair[:,0:1]], 
                                           events.Tau[tautau_indices_pair[:,1:2]]], 
                                          axis=1)

    # Event level masks
    # if events have electron, muon or tau
    # because of bla-tau pair, each list contains 2 elements
    has_ele = ak.num(etau_pair, axis=1)   == 2
    has_muo = ak.num(mutau_pair, axis=1)  == 2
    has_tau = ak.num(tautau_pair, axis=1) == 2


    # get the triggers separately for events that have etau, mutau or tautau pairs
    # An event can be fired by both electron and cross-tau triggers
    # So, here separating the triggers as
    # etau, mutau and tautau triggers
    # e.g.
    #        trigger_ids      = [ [111000,11151,15151], [111000,11151], [11151,15151] ]
    # then,  has_e_triggers   = [ [T,T,F], [T,T], [T,F] ]
    #        has_tau_triggers = [ [F,F,T], [F,F], [F,T] ]
    # etau_trigger_ids   = [ [111000,11151], [111000,11151], [11151] ]
    # tautau_trigger_ids = [ [15151],        [],             [15151] ]
    # also make sure, that the presence of e,mu or tau is included in the decisions
    # i.e.
    # has_ele                       : [          True           ,       True     ,     False      ]
    # has_e_triggers                : [ [True, True, True, False], [False, False], [True, True]   ]
    # mask_has_e_triggers_and_has_e : [ [True, True, True, False], [False, False], [False, False] ]
    mask_has_e_triggers_and_has_e     = has_e_triggers & has_ele
    mask_has_mu_triggers_and_has_mu   = has_mu_triggers & has_muo
    mask_has_tau_triggers_and_has_tau = has_tau_triggers & has_tau

    # filtering out the info based on the three masks defined just above
    # etau
    etau_trigger_names          = trigger_names[mask_has_e_triggers_and_has_e]
    etau_trigger_types          = trigger_types[mask_has_e_triggers_and_has_e]
    etau_trigger_ids            = trigger_ids[mask_has_e_triggers_and_has_e]
    etau_leg_1_minpt            = leg1_minpt[mask_has_e_triggers_and_has_e] 
    etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_e_triggers_and_has_e]

    # mutau
    mutau_trigger_names          = trigger_names[mask_has_mu_triggers_and_has_mu]
    mutau_trigger_types          = trigger_types[mask_has_mu_triggers_and_has_mu]
    mutau_trigger_ids            = trigger_ids[mask_has_mu_triggers_and_has_mu]
    mutau_leg_1_minpt            = leg1_minpt[mask_has_mu_triggers_and_has_mu] 
    mutau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_mu_triggers_and_has_mu]

    # tautau
    tautau_trigger_names          = trigger_names[mask_has_tau_triggers_and_has_tau]
    tautau_trigger_types          = trigger_types[mask_has_tau_triggers_and_has_tau]
    tautau_trigger_ids            = trigger_ids[mask_has_tau_triggers_and_has_tau]    
    tautau_leg_1_minpt            = leg1_minpt[mask_has_tau_triggers_and_has_tau]
    tautau_leg_2_minpt            = leg2_minpt[mask_has_tau_triggers_and_has_tau] 
    tautau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_tau_triggers_and_has_tau]
    tautau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_tau_triggers_and_has_tau]

    if dotrigobjmatch:
        # electrons, muons and taus are extracted from the pairs
        # to match with the respective trigger-objects
        #    [ele  - triggerObjects_leg1]
        #    [muo  - triggerObjects_leg1]
        #    [tau1 - triggerObjects_leg1 & tau2 - triggerObjects_leg2] | [tau1 - triggerObjects_leg2 | tau2 - triggerObjects_leg1]
        ele  = etau_pair[:,0:1]
        muo  = mutau_pair[:,0:1]
        tau1 = tautau_pair[:,0:1]
        tau2 = tautau_pair[:,1:2]
        # get the p4s
        p4_ele  = get_objs_p4(ele)
        p4_muo  = get_objs_p4(muo)
        p4_tau1 = get_objs_p4(tau1)
        p4_tau2 = get_objs_p4(tau2)

        # to convert the masks to event level
        # e.g. events with etau pair and pass electron triggers
        mask_has_e_triggers_and_has_e_evt_level     = ak.any(mask_has_e_triggers_and_has_e, axis=1)
        mask_has_mu_triggers_and_has_mu_evt_level   = ak.any(mask_has_mu_triggers_and_has_mu, axis=1)
        mask_has_tau_triggers_and_has_tau_evt_level = ak.any(mask_has_tau_triggers_and_has_tau, axis=1)
        
        # dummy bool array
        trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])
        
        # It matches the electron for etau pair to the trigger objects responsible for firing a trigger
        # The output follows the shape of the trigger_names or trigger_ids
        el_trigobj_matched_mask = ak.where(mask_has_e_triggers_and_has_e_evt_level,
                                           trigger_object_matching_deep(p4_ele, etau_leg_1_matched_trigobjs, etau_leg_1_minpt, True),
                                           trigobj_matched_mask_dummy)
        el_trigobj_matched_mask = el_trigobj_matched_mask[:,0]
        # later, we will use ak.any to make it an event level mask to filter the trigger object matched electrons / etau pairs
        
        # same for mutau
        mu_trigobj_matched_mask = ak.where(mask_has_mu_triggers_and_has_mu_evt_level,
                                           trigger_object_matching_deep(p4_muo, mutau_leg_1_matched_trigobjs, mutau_leg_1_minpt, True),
                                           trigobj_matched_mask_dummy)
        mu_trigobj_matched_mask = mu_trigobj_matched_mask[:,0]
        
        # For tau-tau, it is a bit lengthy
        # matching tau1 to the passed triggers leg1
        tau1_trigobj_matched_mask_leg1 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau1, tautau_leg_1_matched_trigobjs, tautau_leg_1_minpt, True),
                                                  trigobj_matched_mask_dummy)
        tau1_trigobj_matched_mask_leg1 = tau1_trigobj_matched_mask_leg1[:,0]
        
        # tau1 to leg2
        tau1_trigobj_matched_mask_leg2 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau1, tautau_leg_2_matched_trigobjs, tautau_leg_2_minpt, True),
                                                  trigobj_matched_mask_dummy)
        tau1_trigobj_matched_mask_leg2 = tau1_trigobj_matched_mask_leg2[:,0]
        
        # tau2 to leg1
        tau2_trigobj_matched_mask_leg1 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau2, tautau_leg_1_matched_trigobjs, tautau_leg_1_minpt, True),
                                                  trigobj_matched_mask_dummy)
        tau2_trigobj_matched_mask_leg1 = tau2_trigobj_matched_mask_leg1[:,0]
        
        # tau2 to leg2
        tau2_trigobj_matched_mask_leg2 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau2, tautau_leg_2_matched_trigobjs, tautau_leg_2_minpt, True),
                                                  trigobj_matched_mask_dummy)
        tau2_trigobj_matched_mask_leg2 = tau2_trigobj_matched_mask_leg2[:,0]
        
        # combining decisions
        # tau1 is matched with leg1 and, tau2 is matched with leg2
        # or tau1 with leg2 and, tau2 with leg1
        tau_trigobj_matched_mask = ( (tau1_trigobj_matched_mask_leg1 & tau2_trigobj_matched_mask_leg2)
                                     | (tau1_trigobj_matched_mask_leg2 & tau2_trigobj_matched_mask_leg1) )
        
        


        # ---------- For DEBUGGING
        #from IPython import embed; embed()
        #e_dr = p4_ele.metric_table(etau_leg_1_matched_trigobjs)
        #mu_dr = p4_muo.metric_table(mutau_leg_1_matched_trigobjs)
        #tau1leg1_dr = p4_tau1.metric_table(tautau_leg_1_matched_trigobjs)
        #tau1leg2_dr = p4_tau1.metric_table(tautau_leg_2_matched_trigobjs)
        #tau2leg1_dr = p4_tau2.metric_table(tautau_leg_1_matched_trigobjs)
        #tau2leg2_dr = p4_tau2.metric_table(tautau_leg_2_matched_trigobjs)
        #for i in range(10100):
        #    if not (has_ele[i] | has_tau[i] | has_muo[i]): continue
        #    print("Ele: ", has_ele[i], ele.pt[i], has_e_triggers[i], etau_trigger_ids[i], etau_leg_1_minpt[i], etau_leg_1_matched_trigobjs.pt[i], '--', e_dr[i], el_trigobj_matched_mask[i])
        #    print("Muo: ", has_muo[i], muo.pt[i], has_mu_triggers[i], mutau_trigger_ids[i], mutau_leg_1_minpt[i], mutau_leg_1_matched_trigobjs.pt[i], '--', mu_dr[i], mu_trigobj_matched_mask[i])
        #    print("Tau: ", has_tau[i], tau1.pt[i], tau2.pt[i], has_tau_triggers[i], tautau_trigger_ids[i], tautau_leg_1_minpt[i], tautau_leg_2_minpt[i], tautau_leg_1_matched_trigobjs.pt[i], tautau_leg_2_matched_trigobjs.pt[i], '--', tau1leg1_dr[i], tau1_trigobj_matched_mask_leg1[i], tau1leg2_dr[i], tau1_trigobj_matched_mask_leg2[i], tau2leg1_dr[i], tau2_trigobj_matched_mask_leg1[i], tau2leg2_dr[i], tau2_trigobj_matched_mask_leg2[i], "===>>>", tau_trigobj_matched_mask[i])
        #    print('\n')


        
        # filter out the triggers 
        etau_trigger_ids   = etau_trigger_ids[el_trigobj_matched_mask]
        etau_trigger_names = etau_trigger_names[el_trigobj_matched_mask]
        
        mutau_trigger_ids   = mutau_trigger_ids[mu_trigobj_matched_mask]
        mutau_trigger_names = mutau_trigger_names[mu_trigobj_matched_mask]
        
        tautau_trigger_ids   = tautau_trigger_ids[tau_trigobj_matched_mask]
        tautau_trigger_names = tautau_trigger_names[tau_trigobj_matched_mask]
        
        
        # Get the event level mask from trigger level
        # to see, if an electron matches to any of the trigger objects of any of the triggers
        el_trigobj_matched_mask_evt_level  = ak.any(el_trigobj_matched_mask, axis=1)
        mu_trigobj_matched_mask_evt_level  = ak.any(mu_trigobj_matched_mask, axis=1)
        tau_trigobj_matched_mask_evt_level = ak.any(tau_trigobj_matched_mask, axis=1)
        
        # apply on etau
        etau_indices_pair_dummy = etau_indices_pair[:,:0]
        etau_indices_pair = ak.where(el_trigobj_matched_mask_evt_level, etau_indices_pair, etau_indices_pair_dummy)

        # apply on mutau
        mutau_indices_pair_dummy = mutau_indices_pair[:,:0]
        mutau_indices_pair = ak.where(mu_trigobj_matched_mask_evt_level, mutau_indices_pair, mutau_indices_pair_dummy)
        
        # apply on tautau
        tautau_indices_pair_dummy = tautau_indices_pair[:,:0]
        tautau_indices_pair = ak.where(tau_trigobj_matched_mask_evt_level, tautau_indices_pair, tautau_indices_pair_dummy)
        tautau_pair         = ak.concatenate([events.Tau[tautau_indices_pair[:,0:1]], 
                                              events.Tau[tautau_indices_pair[:,1:2]]], 
                                             axis=1)
        # DO WE NEED TO SORT THE TAU-TAU PAIR CANDIDATES PT SORTED?
        # ANYWAY ... I AM DOING THAT - BABUSHCHA
        #from IPython import embed; embed()
        tautau_indices_pair = tautau_indices_pair[ak.argsort(tautau_pair.pt, axis=1, ascending=False)]



    matchedDict = {
        "etau" : {
            "pair_indices"   : etau_indices_pair,
            "trigger_ids"    : etau_trigger_ids,
            "trigger_names"  : etau_trigger_names,
        },
        "mutau" : {
            "pair_indices"   : mutau_indices_pair,
            "trigger_ids"    : mutau_trigger_ids,
            "trigger_names"  : mutau_trigger_names,
        },
        "tautau" : {
            "pair_indices"   : tautau_indices_pair,
            "trigger_ids"    : tautau_trigger_ids,
            "trigger_names"  : tautau_trigger_names,
        },
    }

    matchedResults = SelectionResult(
        aux=matchedDict,
    )

    return events, matchedResults
