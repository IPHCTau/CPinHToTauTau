# coding: utf-8
"""
A new and generalised approach for Trigger-Object matching
Not very smartly written ... but algo is smart, if works ofc :(
by - GS
"""
import time

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
        "Tau.pt",      "Tau.eta",      "Tau.phi",      "Tau.mass",
    },
    exposed=False
)
def match_trigobj_etau(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        etau_pair: ak.Array,
        dotrigobjmatch: Optional[bool] = True,
        **kwargs
) -> ak.Array:

    start = time.time()

    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_names         = trigger_results.x.trigger_names
    trigger_types         = trigger_results.x.trigger_types
    trigger_ids           = trigger_results.x.trigger_ids
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
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
    from IPython import embed; embed()
    etau_pair  = filter_by_triggers(etau_pair, has_e_triggers)

    # etau pair is a cartesian object
    # so unzipping
    ele, tau = ak.unzip(etau_pair)
    
    # Event level masks
    # if events have electron, muon or tau
    # because of bla-tau pair, each list contains 2 elements
    has_ele = ak.num(ele, axis=1)   > 0


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
    
    # filtering out the info based on the masks defined just above
    # for etau and mutau type of events, separate masks are created
    # for single and cross triggered events
    # for single e triggers
    single_etau_trigger_names          = trigger_names[mask_has_single_e_triggers_and_has_e]
    single_etau_trigger_types          = trigger_types[mask_has_single_e_triggers_and_has_e]
    single_etau_trigger_ids            = trigger_ids[mask_has_single_e_triggers_and_has_e]
    single_etau_leg_1_minpt            = leg1_minpt[mask_has_single_e_triggers_and_has_e]
    single_etau_leg_1_maxeta           = leg1_maxeta[mask_has_single_e_triggers_and_has_e] 
    single_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_single_e_triggers_and_has_e]
    # for cross e triggers
    cross_etau_trigger_names          = trigger_names[mask_has_cross_e_triggers_and_has_e]
    cross_etau_trigger_types          = trigger_types[mask_has_cross_e_triggers_and_has_e]
    cross_etau_trigger_ids            = trigger_ids[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_1_minpt            = leg1_minpt[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_2_minpt            = leg2_minpt[mask_has_cross_e_triggers_and_has_e] 
    cross_etau_leg_1_maxeta           = leg1_maxeta[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_2_maxeta           = leg2_maxeta[mask_has_cross_e_triggers_and_has_e] 
    cross_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_cross_e_triggers_and_has_e]
    cross_etau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_cross_e_triggers_and_has_e]
    # concatenating single and cross info, so that it does not remain depend on the trigger hierarchy
    etau_trigger_names          = ak.concatenate([single_etau_trigger_names, cross_etau_trigger_names], axis=-1)
    etau_trigger_types          = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=-1)
    etau_trigger_ids            = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=-1)
        
    if dotrigobjmatch:
        # get the p4s
        p4_ele = get_objs_p4(ele)
        p4_tau = get_objs_p4(tau)

        # to convert the masks to event level
        # e.g. events with etau pair and pass electron triggers
        mask_has_single_e_triggers_and_has_e_evt_level = ak.any(mask_has_single_e_triggers_and_has_e, axis=1)
        mask_has_cross_e_triggers_and_has_e_evt_level  = ak.any(mask_has_cross_e_triggers_and_has_e, axis=1)
        mask_has_e_triggers_and_has_e_evt_level        = ak.any(mask_has_e_triggers_and_has_e, axis=1)

        from IPython import embed; embed()
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
                                                                               single_etau_leg_1_maxeta,
                                                                               True),
                                                  trigobj_matched_mask_dummy)
        single_el_trigobj_matched_mask = ak.fill_none(ak.firsts(single_el_trigobj_matched_mask, axis=1), False)
        
        print("aaa")
        from IPython import embed; embed()
        cross_el_trigobj_matched_mask_leg1 = ak.where(mask_has_cross_e_triggers_and_has_e_evt_level,
                                                      trigger_object_matching_deep(p4_ele,
                                                                                   cross_etau_leg_1_matched_trigobjs,
                                                                                   cross_etau_leg_1_minpt,
                                                                                   cross_etau_leg_1_maxeta,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)

        cross_el_trigobj_matched_mask_leg1 = ak.fill_none(ak.firsts(cross_el_trigobj_matched_mask_leg1, axis=1), False)

        print("bbb")
        from IPython import embed; embed()
        cross_el_trigobj_matched_mask_leg2 = ak.where(mask_has_cross_e_triggers_and_has_e_evt_level,
                                                      trigger_object_matching_deep(p4_tau,
                                                                                   cross_etau_leg_2_matched_trigobjs,
                                                                                   cross_etau_leg_2_minpt,
                                                                                   cross_etau_leg_2_maxeta,
                                                                                   True),
                                                      trigobj_matched_mask_dummy)

        cross_el_trigobj_matched_mask_leg2 = ak.fill_none(ak.firsts(cross_el_trigobj_matched_mask_leg2, axis=1), False)



        # ensures that both legs are matched

        print("ccc")
        from IPython import embed; embed()
        cross_el_trigobj_matched_mask = (cross_el_trigobj_matched_mask_leg1 & cross_el_trigobj_matched_mask_leg2)
        # concatenating masks for single and cross ele
        print("ccc2")
        from IPython import embed; embed()
        #el_trigobj_matched_mask = ak.concatenate([single_el_trigobj_matched_mask,
        #                                          cross_el_trigobj_matched_mask],
        #                                         axis=-1)
        el_trigobj_matched_mask = ak.concatenate([single_el_trigobj_matched_mask,
                                                  cross_el_trigobj_matched_mask],
                                                 axis=1)

        print("ddd")
        from IPython import embed; embed()
        #el_trigobj_matched_mask = el_trigobj_matched_mask[:,0]
        el_trigobj_matched_mask = ak.fill_none(ak.firsts(el_trigobj_matched_mask, axis=-1), False)
        # later, we will use ak.any to make it an event level mask to filter the trigger object matched electrons / etau pairs
        
        print("eee")
        from IPython import embed; embed()

        # ---------- For DEBUGGING
        if self.config_inst.x.verbose.selection.trigobject_matching:
            ## etau
            etau_e_dr      = p4_ele.metric_table(single_etau_leg_1_matched_trigobjs)
            etau_e_dr_leg1 = p4_ele.metric_table(cross_etau_leg_1_matched_trigobjs)
            etau_t_dr_leg2 = p4_tauele.metric_table(cross_etau_leg_2_matched_trigobjs)
            #
            for i in range(1000):
                if not (has_ele[i]): continue
                print(f"Event       : {events.event[i]}")
                print(f"Trigger IDs : {trigger_ids[i]}")
                print(f"-----------------------------------------------------------")
                print(f"has etau pair?  : {has_ele[i]}")
                print(f" ele_pt         : {ele.pt[i]}, tau_pt: {tauele.pt[i]}")
                print(f" trig?          : {has_e_triggers[i]}")
                print(f" etau_trigId    : {etau_trigger_ids[i]}, single? {has_single_e_triggers[i]}, cross? {has_cross_e_triggers[i]}")
                print(f"   single       : trigger_ids {single_etau_trigger_ids[i]}")
                print(f"    leg1pt      : {single_etau_leg_1_minpt[i]},  l1_match_trigobjs_pt  : {single_etau_leg_1_matched_trigobjs.pt[i]}")
                print(f"    leg1abseta  : {single_etau_leg_1_maxeta[i]}, l1_match_trigobjs_abseta : {np.abs(single_etau_leg_1_matched_trigobjs.eta[i])}")
                print(f"    eleg1dr     : {etau_e_dr[i]}, mask: {single_el_trigobj_matched_mask[i]}")
                print(f"   cross        : trigger_ids {cross_etau_trigger_ids[i]}")
                print(f"    leg1pt      : {cross_etau_leg_1_minpt[i]},  l1_match_trigobjs_pt     : {cross_etau_leg_1_matched_trigobjs.pt[i]}")
                print(f"    leg1abseta  : {cross_etau_leg_1_maxeta[i]}, l1_match_trigobjs_abseta : {np.abs(cross_etau_leg_1_matched_trigobjs.eta[i])}")
                print(f"    leg2pt      : {cross_etau_leg_2_minpt[i]},  l2_match_trigobjs_pt     : {cross_etau_leg_2_matched_trigobjs.pt[i]}")
                print(f"    leg2abseta  : {cross_etau_leg_2_maxeta[i]}, l2_match_trigobjs_abseta : {np.abs(cross_etau_leg_2_matched_trigobjs.eta[i])}")
                print(f"    eleg1dr     : {etau_e_dr_leg1[i]}, tleg2dr: {etau_t_dr_leg2[i]}, eleg1dr & tleg2dr : {cross_el_trigobj_matched_mask[i]}")
                print(f"   comb ------> : {el_trigobj_matched_mask[i]}")
                print(f"-----------------------------------------------------------")            
                print('\n')
        
        print("fff")
        from IPython import embed; embed()
            
        # filter out the triggers 
        etau_trigger_ids   = etau_trigger_ids[el_trigobj_matched_mask]
        etau_trigger_names = etau_trigger_names[el_trigobj_matched_mask]
        etau_trigger_types = etau_trigger_types[el_trigobj_matched_mask]
        
        # Get the event level mask from trigger level
        # to see, if an electron matches to any of the trigger objects of any of the triggers
        el_trigobj_matched_mask_evt_level  = ak.any(el_trigobj_matched_mask, axis=1)
        
        # apply on etau
        etau_pair_dummy = etau_pair[:,:0]
        etau_pair = ak.where(el_trigobj_matched_mask_evt_level, etau_pair, etau_pair_dummy)


    matchedDict = {
        "etau" : {
            "pairs"   : etau_pair,
            "trigger_ids"    : etau_trigger_ids,
            "trigger_names"  : etau_trigger_names,
            "trigger_types"  : etau_trigger_types,
        },
    }

    matchedResults = SelectionResult(
        aux=matchedDict,
    )

    stop = time.time()
    if self.config_inst.x.verbose.selection.trigobject_matching:
        print(f"etau triggerobj matching takes : {round((stop-start)/60.0, 3)}")
    
    return events, matchedResults




