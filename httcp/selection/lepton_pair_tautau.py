# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

from httcp.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep



def match_trigobjs(
        leps_pair: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[ak.Array, ak.Array]:

    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_ids           = trigger_results.x.trigger_ids
    trigger_types         = trigger_results.x.trigger_types
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs

    has_tau_triggers     = ((trigger_types == "cross_tau_tau") | (trigger_types == "cross_tau_tau_jet"))

    old_leps_pair = leps_pair

    leps_pair = filter_by_triggers(leps_pair, has_tau_triggers) 

    taus1, taus2 = ak.unzip(leps_pair)
    
    # Event level masks
    # if events have tau
    has_tau_pairs = ak.fill_none(ak.num(taus1, axis=1) > 0, False)

    # events must be fired by tau triggers and there is ta inside
    mask_has_tau_triggers_and_has_tau_pairs = has_tau_triggers & has_tau_pairs

    tautau_trigger_ids            = trigger_ids[mask_has_tau_triggers_and_has_tau_pairs]    
    tautau_trigger_types          = trigger_types[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_1_minpt            = leg1_minpt[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_minpt            = leg2_minpt[mask_has_tau_triggers_and_has_tau_pairs] 
    tautau_leg_1_maxeta           = leg1_maxeta[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_maxeta           = leg2_maxeta[mask_has_tau_triggers_and_has_tau_pairs] 
    tautau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]


    # tau1            : [ [    t11,        t12      ], [  t11  ] ]
    # tau2            : [ [    t21,        t22      ], [  t21  ] ]

    # trigger_ids     : [ [    1120,       1121     ], [  1120 ] ]
    # leg1/2_minpt    : [ [     40 ,        40      ], [   40  ] ]
    # leg1/2_trigobjs : [ [[o1,o2,o3], [o1,o2,o3,o4]], [[o1,o2]] ]

    # In [67]: taus1.pt[1044]
    # Out[67]: <Array [111, 111, 58.7] type='3 * float32'>

    # In [72]: ak.to_list(tautau_leg_1_matched_trigobjs.pt[1044])
    # Out[72]: 
    # [[110.203125, 53.0, 61.3828125, 210.625],
    #  [110.203125, 53.0, 61.3828125, 210.625, 35.0703125]]
    
    # In [70]: ak.to_list(ak.any(dr[1044] < 0.5, axis=-1))
    # Out[70]: [[True, True], [True, True], [True, True]]

    # same for tau1-legs2 and for tau2
    
    # final mask may look like this:
    # [[True, True], [True, True], [False, True]]

    """
    dr_taus1_leg1 = taus1.metric_table(tautau_leg_1_matched_trigobjs)
    pass_1_1_triglevel = ak.fill_none(ak.any(dr_taus1_leg1 < 0.5, axis=-1), False)
    dr_taus1_leg2 = taus1.metric_table(tautau_leg_2_matched_trigobjs)
    pass_1_2_triglevel = ak.fill_none(ak.any(dr_taus1_leg2 < 0.5, axis=-1), False)
    pass_taus1_legs = pass_1_1_triglevel | pass_1_2_triglevel

    dr_taus2_leg1 = taus2.metric_table(tautau_leg_1_matched_trigobjs)
    pass_2_1_triglevel = ak.fill_none(ak.any(dr_taus2_leg1 < 0.5, axis=-1), False)
    dr_taus2_leg2 = taus2.metric_table(tautau_leg_2_matched_trigobjs)
    pass_2_2_triglevel = ak.fill_none(ak.any(dr_taus2_leg2 < 0.5, axis=-1), False)
    pass_taus2_legs = pass_2_1_triglevel | pass_2_2_triglevel
    """

    # is tau1 matched to leg1?
    pass_tau1_leg1_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau1 matched to leg2?
    pass_tau1_leg2_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)

    # is tau2 matched to leg1?
    pass_tau2_leg1_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau2 matched to leg2?
    pass_tau2_leg2_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)


    # tau1 to leg1 & tau2 to leg2 or, tau1 to leg2 & tau2 to leg1
    pass_taus_legs = (pass_tau1_leg1_triglevel & pass_tau2_leg2_triglevel) | (pass_tau1_leg2_triglevel & pass_tau2_leg1_triglevel)

    
    mask_has_tau_triggers_and_has_tau_pairs_evt_level = ak.fill_none(ak.any(mask_has_tau_triggers_and_has_tau_pairs, axis=1), False)
    trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])
    
    pass_taus_legs = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, pass_taus_legs, trigobj_matched_mask_dummy)
    pass_taus_legs = ak.enforce_type(ak.values_astype(pass_taus_legs, "bool"), "var * var * bool") # 100000 * var * var * bool
    
    
    pass_taus = ak.fill_none(ak.any(pass_taus_legs, axis=-1), False)

    #tautau_trigger_ids_dummy = tautau_trigger_ids[:,:0]
    #tautau_trigger_ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, tautau_trigger_ids, tautau_trigger_ids_dummy)
    
    
    tautau_trigger_ids_brdcst, _ = ak.broadcast_arrays(tautau_trigger_ids[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_ids_brdcst = tautau_trigger_ids_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids
    ids = ak.fill_none(ak.firsts(tautau_trigger_ids_brdcst, axis=-1), -1)

    #ids_dummy = ak.from_regular((trigger_ids > 0)[:,:0])
    #ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, ids, ids_dummy)

    
    ids = ak.values_astype(ids, 'int64')

    new_taus1 = taus1[pass_taus]
    new_taus2 = taus2[pass_taus]
    #ids = ids[pass_taus]
    
    leps_pair = ak.zip([new_taus1, new_taus2])



    tautau_trigger_types_brdcst, _ = ak.broadcast_arrays(tautau_trigger_types[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_types_brdcst = tautau_trigger_types_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids
    types = ak.fill_none(ak.firsts(tautau_trigger_types_brdcst, axis=-1), "")

    #from IPython import embed; embed()
    
    
    return leps_pair, ids, types




def sort_pairs(dtrpairs: ak.Array)->ak.Array:
    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    where_same_iso_2 = ak.fill_none(where_same_iso_2, False)

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False
    )
    
    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)

    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    #lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    #lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    #dtrpair    = ak.concatenate([lep1, lep2], axis=1) 

    return dtrpairs



@selector(
    uses={
        optional("Tau.pt"),
        optional("Tau.pt_tautau"),
        optional("Tau.mass"),
        optional("Tau.mass_tautau"),
        "Tau.eta", "Tau.phi",
        "Tau.rawIdx", optional("Tau.genPartFlav"),
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    taus            = events.Tau[lep_indices]
    # Extra channel specific selections on tau
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["tautau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["tautau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["tautau"]

    is_good_tau     = (
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    # rename the {channel}_pt/mass to pt/mass
    if self.dataset_inst.is_mc:
        taus = ak.without_field(taus, "pt")
        taus = ak.with_field(taus, taus.pt_tautau, "pt")
        taus = ak.without_field(taus, "mass")
        taus = ak.with_field(taus, taus.mass_tautau, "mass")
    
    # Sorting leps [Tau] by deeptau [descending]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]

    leps_pair  = ak.combinations(taus, 2, axis=1)    
    lep1, lep2 = ak.unzip(leps_pair)

    preselection = {
        #"tautau_tau1_iso"      : (lep1.idDeepTau2018v2p5VSjet >= tau_tagger_wps.vs_j[vs_jet_wp]),
        "tautau_is_pt_40"      : (lep1.pt > 40) & (lep2.pt > 40),
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        #"tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        "tautau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40, # invariant_mass(lep1, lep2) > 40
    }

    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    category_selections = {}
    pair_selection_steps["tautau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair = leps_pair[good_pair_mask]

    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["tautau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs(leps_pair, trigger_results)

    
    pair_selection_steps["tautau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    lep1, lep2 = ak.unzip(leps_pair)

    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    trigId = trigIds[:,:1]
    trigTypes = trigTypes[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    sort_idx = ak.argsort(leps_pair.pt, ascending=False)
    leps_pair = leps_pair[sort_idx]


    return SelectionResult(
        aux = pair_selection_steps,
    ), leps_pair, trigIds, trigTypes
