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

    # check etau triggers, filter out the indices based on any etau
    # trigger passed or not and then get the etau pair
    # separately done for single and cross ele triggers
    # e.g.
    #   trigger_ids = [ [111000, 112000, 11151], [111000], [11151] ] 
    #   has_single_e_triggers =
    #                 [ [ True ,  True , False], [ True ], [False] ] 
    #   has_cross_e_triggers = 
    #                 [ [False , False ,  True], [False ], [ True] ] 
    #   has_e_triggers = has_single_e_triggers | has_cross_e_triggers
    #                 [ [ True ,  True ,  True], [ True ], [ True] ]
    # Make sure that events with etau pair are fired by single or cross ele triggers
    # and filter_by_triggers func is basically checking if an event has any of the ele triggers
    # Finally, etau_pair will exist in those events only that has any of the single/cross ele triggers 
    has_single_e_triggers = trigger_types == "single_e"
    has_cross_e_triggers  = trigger_types == "cross_e_tau"
    has_e_triggers = (has_single_e_triggers | has_cross_e_triggers)
    etau_indices_pair  = filter_by_triggers(etau_indices_pair, has_e_triggers)
    etau_pair = ak.concatenate([events.Electron[etau_indices_pair[:,0:1]],
                                events.Tau[etau_indices_pair[:,1:2]]],
                               axis=1)
    # same for mutau pairs
    has_single_mu_triggers  = trigger_types == "single_mu"
    has_cross_mu_triggers   = trigger_types == "cross_mu_tau"
    has_mu_triggers     = (has_single_mu_triggers | has_cross_mu_triggers)
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
    mask_has_single_e_triggers_and_has_e = has_single_e_triggers & has_ele
    mask_has_cross_e_triggers_and_has_e  = has_cross_e_triggers & has_ele
    mask_has_e_triggers_and_has_e        = has_e_triggers & has_ele

    mask_has_single_mu_triggers_and_has_mu = has_single_mu_triggers & has_muo
    mask_has_cross_mu_triggers_and_has_mu  = has_cross_mu_triggers & has_muo
    mask_has_mu_triggers_and_has_mu        = has_mu_triggers & has_muo

    mask_has_tau_triggers_and_has_tau = has_tau_triggers & has_tau

    
    # filtering out the info based on the masks defined just above
    # for etau and mutau type of events, separate masks are created
    # for single and cross triggered events
    # for single e triggers
    single_etau_trigger_names          = trigger_names[mask_has_single_e_triggers_and_has_e]
    single_etau_trigger_types          = trigger_types[mask_has_single_e_triggers_and_has_e]
    single_etau_trigger_ids            = trigger_ids[mask_has_single_e_triggers_and_has_e]
    single_etau_leg_1_minpt            = leg1_minpt[mask_has_single_e_triggers_and_has_e] 
    single_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_single_e_triggers_and_has_e]
    # for cross e triggers
    cross_etau_trigger_names          = trigger_names[mask_has_cross_e_triggers_and_has_e]
    cross_etau_trigger_types          = trigger_types[mask_has_cross_e_triggers_and_has_e]
    cross_etau_trigger_ids            = trigger_ids[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_1_minpt            = leg1_minpt[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_2_minpt            = leg2_minpt[mask_has_cross_e_triggers_and_has_e] 
    cross_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_cross_e_triggers_and_has_e]
    # concatenating single and cross info, so that it does not remain depend on the trigger hierarchy
    etau_trigger_names          = ak.concatenate([single_etau_trigger_names, cross_etau_trigger_names], axis=-1)
    etau_trigger_types          = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=-1)
    etau_trigger_ids            = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=-1)
    

    # same for mutau
    # for single mu triggers
    single_mutau_trigger_names          = trigger_names[mask_has_single_mu_triggers_and_has_mu]
    single_mutau_trigger_types          = trigger_types[mask_has_single_mu_triggers_and_has_mu]
    single_mutau_trigger_ids            = trigger_ids[mask_has_single_mu_triggers_and_has_mu]
    single_mutau_leg_1_minpt            = leg1_minpt[mask_has_single_mu_triggers_and_has_mu] 
    single_mutau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_single_mu_triggers_and_has_mu]
    # for cross mu triggers
    cross_mutau_trigger_names          = trigger_names[mask_has_cross_mu_triggers_and_has_mu]
    cross_mutau_trigger_types          = trigger_types[mask_has_cross_mu_triggers_and_has_mu]
    cross_mutau_trigger_ids            = trigger_ids[mask_has_cross_mu_triggers_and_has_mu]
    cross_mutau_leg_1_minpt            = leg1_minpt[mask_has_cross_mu_triggers_and_has_mu]
    cross_mutau_leg_2_minpt            = leg2_minpt[mask_has_cross_mu_triggers_and_has_mu] 
    cross_mutau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_cross_mu_triggers_and_has_mu]
    cross_mutau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_cross_mu_triggers_and_has_mu]
    # mutau
    mutau_trigger_names          = ak.concatenate([single_mutau_trigger_names, cross_mutau_trigger_names], axis=-1)
    mutau_trigger_types          = ak.concatenate([single_mutau_trigger_types, cross_mutau_trigger_types], axis=-1)
    mutau_trigger_ids            = ak.concatenate([single_mutau_trigger_ids, cross_mutau_trigger_ids], axis=-1)

    
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
        tauele = etau_pair[:,1:2]
        muo  = mutau_pair[:,0:1]
        taumuo = mutau_pair[:,1:2]
        tau1 = tautau_pair[:,0:1]
        tau2 = tautau_pair[:,1:2]
        # get the p4s
        p4_ele  = get_objs_p4(ele)
        p4_tauele = get_objs_p4(tauele)
        p4_muo  = get_objs_p4(muo)
        p4_taumuo = get_objs_p4(taumuo)
        p4_tau1 = get_objs_p4(tau1)
        p4_tau2 = get_objs_p4(tau2)

        # to convert the masks to event level
        # e.g. events with etau pair and pass electron triggers
        mask_has_single_e_triggers_and_has_e_evt_level = ak.any(mask_has_single_e_triggers_and_has_e, axis=1)
        mask_has_cross_e_triggers_and_has_e_evt_level  = ak.any(mask_has_cross_e_triggers_and_has_e, axis=1)
        mask_has_e_triggers_and_has_e_evt_level        = ak.any(mask_has_e_triggers_and_has_e, axis=1)

        mask_has_single_mu_triggers_and_has_mu_evt_level = ak.any(mask_has_single_mu_triggers_and_has_mu, axis=1)
        mask_has_cross_mu_triggers_and_has_mu_evt_level  = ak.any(mask_has_cross_mu_triggers_and_has_mu, axis=1)
        mask_has_mu_triggers_and_has_mu_evt_level        = ak.any(mask_has_mu_triggers_and_has_mu, axis=1)

        mask_has_tau_triggers_and_has_tau_evt_level = ak.any(mask_has_tau_triggers_and_has_tau, axis=1)
        
        # dummy bool array
        trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])
        
        # It matches the electron for etau pair to the trigger objects responsible for firing a trigger
        # The output follows the shape of the trigger_names or trigger_ids
        # etau
        # Algo:
        #   if matches to a single e trigger, check the matching for e only
        #   else if matches to a cross e trigger, check the matching for both e and tau
        #   same for mu
        #   for tautau, match both legs
        single_el_trigobj_matched_mask = ak.where(mask_has_single_e_triggers_and_has_e_evt_level,
                                                  trigger_object_matching_deep(p4_ele,
                                                                               single_etau_leg_1_matched_trigobjs,
                                                                               single_etau_leg_1_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        cross_el_trigobj_matched_mask_leg1 = ak.where(mask_has_cross_e_triggers_and_has_e_evt_level,
                                                      trigger_object_matching_deep(p4_ele,
                                                                                   cross_etau_leg_1_matched_trigobjs,
                                                                                   cross_etau_leg_1_minpt,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)
        cross_el_trigobj_matched_mask_leg2 = ak.where(mask_has_cross_e_triggers_and_has_e_evt_level,
                                                      trigger_object_matching_deep(p4_tauele,
                                                                                   cross_etau_leg_2_matched_trigobjs,
                                                                                   cross_etau_leg_2_minpt,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)
        # ensures that both legs are matched
        cross_el_trigobj_matched_mask = (cross_el_trigobj_matched_mask_leg1 & cross_el_trigobj_matched_mask_leg2)
        # concatenating masks for single and cross ele
        el_trigobj_matched_mask = ak.concatenate([single_el_trigobj_matched_mask,
                                                  cross_el_trigobj_matched_mask],
                                                 axis=-1)

        el_trigobj_matched_mask = el_trigobj_matched_mask[:,0]
        # later, we will use ak.any to make it an event level mask to filter the trigger object matched electrons / etau pairs
        


        # same for mutau
        single_mu_trigobj_matched_mask = ak.where(mask_has_single_mu_triggers_and_has_mu_evt_level,
                                                  trigger_object_matching_deep(p4_muo,
                                                                               single_mutau_leg_1_matched_trigobjs,
                                                                               single_mutau_leg_1_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        cross_mu_trigobj_matched_mask_leg1 = ak.where(mask_has_cross_mu_triggers_and_has_mu_evt_level,
                                                      trigger_object_matching_deep(p4_muo,
                                                                                   cross_mutau_leg_1_matched_trigobjs,
                                                                                   cross_mutau_leg_1_minpt,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)
        cross_mu_trigobj_matched_mask_leg2 = ak.where(mask_has_cross_mu_triggers_and_has_mu_evt_level,
                                                      trigger_object_matching_deep(p4_taumuo,
                                                                                   cross_mutau_leg_2_matched_trigobjs,
                                                                                   cross_mutau_leg_2_minpt,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)
        cross_mu_trigobj_matched_mask = (cross_mu_trigobj_matched_mask_leg1 & cross_mu_trigobj_matched_mask_leg2)
        mu_trigobj_matched_mask = ak.concatenate([single_mu_trigobj_matched_mask,
                                                  cross_mu_trigobj_matched_mask],
                                                 axis=-1)

        mu_trigobj_matched_mask = mu_trigobj_matched_mask[:,0]
        

        

        # For tau-tau, it is a bit lengthy
        # matching tau1 to the passed triggers leg1
        tau1_trigobj_matched_mask_leg1 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau1,
                                                                               tautau_leg_1_matched_trigobjs,
                                                                               tautau_leg_1_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        tau1_trigobj_matched_mask_leg1 = tau1_trigobj_matched_mask_leg1[:,0]
        
        # tau1 to leg2
        tau1_trigobj_matched_mask_leg2 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau1,
                                                                               tautau_leg_2_matched_trigobjs,
                                                                               tautau_leg_2_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        tau1_trigobj_matched_mask_leg2 = tau1_trigobj_matched_mask_leg2[:,0]
        
        # tau2 to leg1
        tau2_trigobj_matched_mask_leg1 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau2,
                                                                               tautau_leg_1_matched_trigobjs,
                                                                               tautau_leg_1_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        tau2_trigobj_matched_mask_leg1 = tau2_trigobj_matched_mask_leg1[:,0]
        
        # tau2 to leg2
        tau2_trigobj_matched_mask_leg2 = ak.where(mask_has_tau_triggers_and_has_tau_evt_level,
                                                  trigger_object_matching_deep(p4_tau2,
                                                                               tautau_leg_2_matched_trigobjs,
                                                                               tautau_leg_2_minpt,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        tau2_trigobj_matched_mask_leg2 = tau2_trigobj_matched_mask_leg2[:,0]
        
        # combining decisions
        # tau1 is matched with leg1 and, tau2 is matched with leg2
        # or tau1 with leg2 and, tau2 with leg1
        tau_trigobj_matched_mask = ( (tau1_trigobj_matched_mask_leg1 & tau2_trigobj_matched_mask_leg2)
                                     | (tau1_trigobj_matched_mask_leg2 & tau2_trigobj_matched_mask_leg1) )
        
        


        # ---------- For DEBUGGING
        #from IPython import embed; embed()
        """
        ## etau
        etau_e_dr      = p4_ele.metric_table(single_etau_leg_1_matched_trigobjs)
        etau_e_dr_leg1 = p4_ele.metric_table(cross_etau_leg_1_matched_trigobjs)
        etau_t_dr_leg2 = p4_tauele.metric_table(cross_etau_leg_2_matched_trigobjs)
        ## mutau
        mutau_m_dr      = p4_muo.metric_table(single_mutau_leg_1_matched_trigobjs)
        mutau_m_dr_leg1 = p4_muo.metric_table(cross_mutau_leg_1_matched_trigobjs)
        mutau_t_dr_leg2 = p4_taumuo.metric_table(cross_mutau_leg_2_matched_trigobjs)
        ## tautau
        tau1leg1_dr = p4_tau1.metric_table(tautau_leg_1_matched_trigobjs)
        tau1leg2_dr = p4_tau1.metric_table(tautau_leg_2_matched_trigobjs)
        tau2leg1_dr = p4_tau2.metric_table(tautau_leg_1_matched_trigobjs)
        tau2leg2_dr = p4_tau2.metric_table(tautau_leg_2_matched_trigobjs)
        #
        for i in range(10100):
            if not (has_ele[i] | has_tau[i] | has_muo[i]): continue
            print(f"Trigger IDs: {trigger_ids[i]}") 
            print(f"Ele?  : {has_ele[i]}")
            print(f" ele_pt : {ele.pt[i]}, tau_pt: {tauele.pt[i]}")
            print(f" etau_trigId : {etau_trigger_ids[i]}, single? {has_single_e_triggers[i]}, cross? {has_cross_e_triggers[i]}")
            print(f"   single    : trigger_ids {single_etau_trigger_ids[i]}")
            print(f"    leg1pt   : {single_etau_leg_1_minpt[i]}, l1_match_trigobjs_pt: {single_etau_leg_1_matched_trigobjs.pt[i]}")
            print(f"    eleg1dr  : {etau_e_dr[i]}, mask: {single_el_trigobj_matched_mask[i]}")
            print(f"   cross     : trigger_ids {cross_etau_trigger_ids[i]}")
            print(f"    leg1pt   : {cross_etau_leg_1_minpt[i]}, l1_match_trigobjs_pt: {cross_etau_leg_1_matched_trigobjs.pt[i]}")
            print(f"    leg2pt   : {cross_etau_leg_2_minpt[i]}, l2_match_trigobjs_pt: {cross_etau_leg_2_matched_trigobjs.pt[i]}")
            print(f"    eleg1dr  : {etau_e_dr_leg1[i]}, tleg2dr: {etau_t_dr_leg2[i]}, eleg1dr & tleg2dr : {cross_el_trigobj_matched_mask[i]}")
            print(f"   comb      : trig_matched_mask: {el_trigobj_matched_mask[i]}")
            print(f"Muo?  : {has_muo[i]}")
            print(f" muo_pt : {ele.pt[i]}, tau_pt: {tauele.pt[i]}")
            print(f" mutau_trigId : {mutau_trigger_ids[i]}, single? {has_single_mu_triggers[i]}, cross? {has_cross_mu_triggers[i]}")
            print(f"   single     : trigger_ids {single_mutau_trigger_ids[i]}")
            print(f"    leg1pt    : {single_mutau_leg_1_minpt[i]}, l1_match_trigobjs_pt: {single_mutau_leg_1_matched_trigobjs.pt[i]}")
            print(f"    mleg1dr   : {mutau_m_dr[i]}, mask: {single_mu_trigobj_matched_mask[i]}")
            print(f"   cross      : trigger_ids {cross_mutau_trigger_ids[i]}")
            print(f"    leg1pt    : {cross_mutau_leg_1_minpt[i]}, l1_match_trigobjs_pt: {cross_mutau_leg_1_matched_trigobjs.pt[i]}")
            print(f"    leg2pt    : {cross_mutau_leg_2_minpt[i]}, l2_match_trigobjs_pt: {cross_mutau_leg_2_matched_trigobjs.pt[i]}")
            print(f"    mleg1dr   : {mutau_mu_dr_leg1[i]}, tleg2dr: {mutau_t_dr_leg2[i]}, muleg1dr & tleg2dr : {cross_mu_trigobj_matched_mask[i]}")
            print(f"   comb       : trig_matched_mask: {mu_trigobj_matched_mask[i]}")
            print(f"Tau?  : {has_tau[i]}, t1pt: {tau1.pt[i]}, t2pt: {tau2.pt[i]}, trig? {has_tau_triggers[i]}, Ids: {tautau_trigger_ids[i]}")
            print(f" \t l1pt: {tautau_leg_1_minpt[i]}, l2pt: {tautau_leg_2_minpt[i]}, l1_match_tobjspt: {tautau_leg_1_matched_trigobjs.pt[i]}, l2_match_tobjspt: : {tautau_leg_2_matched_trigobjs.pt[i]}")
            print(f" \t t1l1dr     : {tau1leg1_dr[i]}")
            print(f" \t t1l1drmask : {tau1_trigobj_matched_mask_leg1[i]}")
            print(f" \t t1l2dr     : {tau1leg2_dr[i]}")
            print(f" \t t1l2drmask : {tau1_trigobj_matched_mask_leg2[i]}")
            print(f" \t t2l1dr     : {tau2leg1_dr[i]}")
            print(f" \t t2l1drmask : {tau2_trigobj_matched_mask_leg1[i]}")
            print(f" \t t2l2dr     : {tau2leg2_dr[i]}")
            print(f" \t t2l2drmask : {tau2_trigobj_matched_mask_leg2[i]}")
            print(f" \t trig_matched_mask  : {tau_trigobj_matched_mask[i]}")
            print('\n')
        """
        
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




