# coding: utf-8

"""
Exemplary selection methods.
"""

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior
from columnflow.production.categories import category_ids

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.selection.physics_objects import *
from httcp.selection.trigger import trigger_selection
from httcp.selection.lepton_pair_etau import etau_selection
from httcp.selection.lepton_pair_mutau import mutau_selection
from httcp.selection.lepton_pair_tautau import tautau_selection
from httcp.selection.event_category import get_categories
#from httcp.selection.match_trigobj import match_trigobj
from httcp.selection.trigobject_matching import match_trigobj
from httcp.selection.lepton_veto import *
from httcp.selection.higgscand import higgscand, higgscandprod

#from httcp.production.main import cutflow_features
#from httcp.production.weights import scale_mc_weight
from httcp.production.extra_weights import scale_mc_weight
from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars

from httcp.util import filter_by_triggers, get_objs_p4, trigger_object_matching_deep, IF_DATASET_IS_DY_LO

from httcp.selection.debug import debug_main


np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@selector(uses={"process_id", optional("mc_weight")}) #TODO: Move it to utils
def custom_increment_stats(
    self: Selector,
    events: ak.Array,
    results: SelectionResult,
    stats: dict,
    **kwargs,
) -> ak.Array:
    """
    Unexposed selector that does not actually select objects but instead increments selection
    *stats* in-place based on all input *events* and the final selection *mask*.
    """
    # get event masks
    event_mask = results.event

    # get a list of unique process ids present in the chunk
    unique_process_ids = np.unique(events.process_id)
    # increment plain counts
    n_evt_per_file = self.dataset_inst.n_events/self.dataset_inst.n_files
    stats["num_events"] = n_evt_per_file
    stats["num_events_selected"] += float(ak.sum(event_mask, axis=0))
    #print(f"stats : {stats}")
    #from IPython import embed; embed()
    if self.dataset_inst.is_mc:
        stats[f"sum_mc_weight"] = n_evt_per_file
        stats.setdefault(f"sum_mc_weight_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_mc_weight_per_process"][int(p)] = n_evt_per_file
        
    # create a map of entry names to (weight, mask) pairs that will be written to stats
    weight_map = OrderedDict()
    if self.dataset_inst.is_mc:
        # mc weight for selected events
        weight_map["mc_weight_selected"] = (events.mc_weight, event_mask)

    # get and store the sum of weights in the stats dictionary
    for name, (weights, mask) in weight_map.items():
        joinable_mask = True if mask is Ellipsis else mask

        # sum of different weights in weight_map for all processes
        stats[f"sum_{name}"] += float(ak.sum(weights[mask]))
        # sums per process id
        stats.setdefault(f"sum_{name}_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_{name}_per_process"][int(p)] += float(ak.sum(
                weights[(events.process_id == p) & joinable_mask],
            ))

    return events, results


def get_2n_pairs(etau_indices_pair,
                 mutau_indices_pair,
                 tautau_indices_pair):

    has_one_etau_pair   = ak.num(etau_indices_pair, axis=1) == 2
    has_one_mutau_pair  = ak.num(mutau_indices_pair, axis=1) == 2
    has_one_tautau_pair = ak.num(tautau_indices_pair, axis=1) == 2

    has_at_least_one_pair = (has_one_etau_pair | has_one_mutau_pair | has_one_tautau_pair)

    return has_at_least_one_pair    

    
# exposed selectors
# (those that can be invoked from the command line)
@selector(
    uses={
        "event",
        # selectors / producers called within _this_ selector
        attach_coffea_behavior,
        json_filter, 
        met_filters, 
        scale_mc_weight, 
        process_ids,
        trigger_selection,
        IF_DATASET_IS_DY_LO(genZ_selection),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories,
        extra_lepton_veto, 
        double_lepton_veto, 
        match_trigobj,
        increment_stats, 
        custom_increment_stats,
        higgscand,
        gentau_selection,
        higgscandprod,
        rel_charge,
        category_ids,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        scale_mc_weight, 
        trigger_selection,
        IF_DATASET_IS_DY_LO(genZ_selection),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories, 
        process_ids,
        extra_lepton_veto, 
        double_lepton_veto, 
        match_trigobj,
        higgscandprod,
        gentau_selection,
        rel_charge,
        category_ids,
        "trigger_ids",
        "single_triggered",
        "cross_triggered",
        "single_e_triggered",
        "single_mu_triggered",
        "cross_e_triggered",
        "cross_mu_triggered",
        "cross_tau_triggered",
        "cross_tau_jet_triggered",
    },
    exposed=True,
)
def main(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    # ################### ensure coffea behaviors are loaded ################### #
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # ############# Sel : Good Lumi JSON for DATA only ############# #
    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    else:
        results += SelectionResult(steps={"json": np.ones(len(events), dtype=bool)})


    # -------- Sel : Met Filters -------- #
    # met filter selection
    events, met_filter_results = self[met_filters](events, **kwargs)
    results += met_filter_results


    # --------- Sel : Trigger Selection -------- #
    # trigger selection
    events, trigger_results = self[trigger_selection](events, **kwargs)
    results += trigger_results


    # -------- Get genZ collection for Zpt reweight
    if self.dataset_inst.has_tag("is_dy_LO"):
        events = self[genZ_selection](events, **kwargs)
        
    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ] 
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += ele_results

    
    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ] 
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += muon_results


    # tau selection
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ]
    good_tau_indices = None
    good_tau_indices_mutau = None
    good_tau_indices_etau  = None
    good_tau_indices_tautau= None

    if self.dataset_inst.is_mc:
        events, tau_results, good_tau_indices_mutau = self[tau_selection](events, "mutau",  call_force=True, **kwargs)
        results += tau_results
        _, _, good_tau_indices_etau           = self[tau_selection](events, "etau",   call_force=True, **kwargs)
        _, _, good_tau_indices_tautau         = self[tau_selection](events, "tautau", call_force=True, **kwargs)
    else:
        events, tau_results, good_tau_indices = self[tau_selection](events, call_force=True, **kwargs)
    
    # double lepton veto
    events, extra_double_lepton_veto_results = self[double_lepton_veto](events,
                                                                        dlveto_ele_indices,
                                                                        dlveto_muon_indices)
    results += extra_double_lepton_veto_results

    
    # check if there are at least two leptons with at least one tau

    #_lepton_indices = ak.concatenate([good_muon_indices, good_ele_indices, good_tau_indices], axis=1)
    #nlep_mask = ((ak.num(_lepton_indices, axis=1) >= 2) & (ak.num(good_tau_indices, axis=1) >= 1))
    #match_nlep = SelectionResult(
    #    steps = {
    #        "has_2_or_more_leps_with_at_least_1_tau"  : nlep_mask,
    #    },
    #)
    #results += match_nlep
    
    #from IPython import embed; embed()
    
    # e-tau pair i.e. hcand selection
    # e.g. [ [], [e1, tau1], [], [], [e1, tau2] ]
    etau_results, etau_pair = self[etau_selection](events,
                                                   good_ele_indices,
                                                   good_tau_indices_etau if self.dataset_inst.is_mc else good_tau_indices,
                                                   call_force=True,
                                                   **kwargs)
    results += etau_results

    # mu-tau pair i.e. hcand selection
    # e.g. [ [mu1, tau1], [], [mu1, tau2], [], [] ]
    mutau_results, mutau_pair = self[mutau_selection](events,
                                                      good_muon_indices,
                                                      good_tau_indices_mutau if self.dataset_inst.is_mc else good_tau_indices,
                                                      call_force=True,
                                                      **kwargs)
    results += mutau_results

    # tau-tau pair i.e. hcand selection
    # e.g. [ [], [tau1, tau2], [], [], [] ]
    tautau_results, tautau_pair = self[tautau_selection](events,
                                                         good_tau_indices_tautau if self.dataset_inst.is_mc else good_tau_indices,
                                                         call_force=True,
                                                         **kwargs)
    results += tautau_results

    pre_match_at_least_one_pair = get_2n_pairs(etau_pair.rawIdx,
                                               mutau_pair.rawIdx,
                                               tautau_pair.rawIdx)


    pre_match_pair_result = SelectionResult(
        steps = {
            "has_at_least_1_pair_before_trigobj_matching"  : pre_match_at_least_one_pair,
        },
    )
    results += pre_match_pair_result

    #from IPython import embed; embed()
    #for i in range(1000):
    #  print(etau_indices_pair[i], mutau_indices_pair[i], tautau_indices_pair[i], pre_match_at_least_one_pair[i])

    
    # Trigger objects matching
    # ################## W A R N I N G ################# #
    # Mind the SWITCH: True  -> TriggerObject matching ON 
    #                  False -> TriggerObject matching OFF
    # although, the trigger requirement per channel is automatic
    # ################################################## #
    events, matchedResults = self[match_trigobj](events,
                                                 trigger_results,
                                                 etau_pair,
                                                 mutau_pair,
                                                 tautau_pair,
                                                 True)
    results += matchedResults


    etau_pair    = matchedResults.x.etau["pairs"]
    mutau_pair   = matchedResults.x.mutau["pairs"]
    tautau_pair  = matchedResults.x.tautau["pairs"]
    
    
    post_match_at_least_one_pair = get_2n_pairs(etau_pair.rawIdx,
                                                mutau_pair.rawIdx,
                                                tautau_pair.rawIdx)
    post_match_pair_result = SelectionResult(
        steps = {
            "has_at_least_1_pair_after_trigobj_matching"  : post_match_at_least_one_pair,
        },
    )
    results += post_match_pair_result


    # channel selection
    # channel_id is now in columns
    events, channel_results = self[get_categories](events,
                                                   trigger_results,
                                                   etau_pair.rawIdx,
                                                   mutau_pair.rawIdx,
                                                   tautau_pair.rawIdx)
    results += channel_results

    etau_trigger_types = matchedResults.x.etau["trigger_types"]
    etau_trigger_ids   = matchedResults.x.etau["trigger_ids"]

    mutau_trigger_types = matchedResults.x.mutau["trigger_types"]
    mutau_trigger_ids   = matchedResults.x.mutau["trigger_ids"]
    
    tautau_trigger_types = matchedResults.x.tautau["trigger_types"]
    tautau_trigger_ids   = matchedResults.x.tautau["trigger_ids"]
    
    trigger_types = ak.concatenate([etau_trigger_types, mutau_trigger_types, tautau_trigger_types], axis=1)
    trigger_ids   = ak.concatenate([etau_trigger_ids, mutau_trigger_ids, tautau_trigger_ids], axis=1)
    trigger_ids = ak.values_astype(ak.fill_none(ak.firsts(trigger_ids, axis=1), 999), np.uint64)

    # save the trigger_ids column
    events = set_ak_column(events, "trigger_ids", trigger_ids)

    # save single_triggered and cross_triggered
    single_e_triggered = ak.any(ak.fill_none((trigger_types == 'single_e'), False), axis=1)
    single_mu_triggered = ak.any(ak.fill_none((trigger_types == 'single_mu'), False), axis=1)
    cross_e_triggered = ak.any(ak.fill_none((trigger_types == 'cross_e_tau'), False), axis=1)
    cross_mu_triggered = ak.any(ak.fill_none((trigger_types == 'cross_mu_tau'), False), axis=1)
    cross_tau_triggered = ak.any(ak.fill_none((trigger_types == 'cross_tau_tau'), False), axis=1)
    cross_tau_jet_triggered = ak.any(ak.fill_none((trigger_types == 'cross_tau_jet_tau'), False), axis=1)

    events = set_ak_column(events, "single_e_triggered", single_e_triggered)
    events = set_ak_column(events, "single_mu_triggered", single_mu_triggered)
    events = set_ak_column(events, "cross_e_triggered", cross_e_triggered)
    events = set_ak_column(events, "cross_mu_triggered", cross_mu_triggered)
    events = set_ak_column(events, "cross_tau_triggered", cross_tau_triggered)
    events = set_ak_column(events, "cross_tau_jet_triggered", cross_tau_jet_triggered)

    events = set_ak_column(events, "single_triggered", (single_e_triggered | single_mu_triggered))
    events = set_ak_column(events, "cross_triggered", (cross_e_triggered | cross_mu_triggered | cross_tau_triggered | cross_tau_jet_triggered))


    # hcand pair: [ [[mu1,tau1]], [[e1,tau1],[tau1,tau2]], [[mu1,tau2]], [], [[e1,tau2]] ]
    hcand_pairs = ak.concatenate([etau_pair[:,None], mutau_pair[:,None], tautau_pair[:,None]], axis=1)

    # extra lepton veto
    # it is only applied on the events with one higgs candidate only
    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_ele_indices,
                                                                veto_muon_indices,
                                                                hcand_pairs)
    results += extra_lepton_veto_results

    #from IPython import embed; embed()
    # hcand results
    events, hcand_array, hcand_results = self[higgscand](events, hcand_pairs)
    results += hcand_results
    
    # hcand prod results
    # -------------------------------------------- #
    # Here, inside hcandprod, the replacement of raw
    # taus by calibrated taus takes place
    # -------------------------------------------- #
    events, hcandprod_results = self[higgscandprod](events, hcand_array)
    results += hcandprod_results

    
    # -------- Sel : b-veto -------- #
    # jet selection
    # -------------------------------------------- #
    # this is moved here, because now the jets are
    # cleaned against the tau cadidates of hacnd
    # -------------------------------------------- #
    events, bjet_veto_result = self[jet_selection](events, 
                                                   call_force=True, 
                                                   **kwargs)
    results += bjet_veto_result
    
    
    # gen particles info
    # ############################################ #
    # After building the higgs candidates, one can
    # switch on the production of GenTau. Those gentaus
    # will be selected which get matched to the hcand
    # hcand-gentau match = True/False
    # ############################################ #
    if self.config_inst.x.extra_tags.genmatch:
        if "is_signal" in list(self.dataset_inst.aux.keys()):
            print(" --->>> hcand-gentau matching")
            events, gentau_results = self[gentau_selection](events, True)
            results += gentau_results


    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    
    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[scale_mc_weight](events, **kwargs)

    # rel-charge
    events = self[rel_charge](events, **kwargs)

    # add cutflow features, passing per-object masks
    #events = self[cutflow_features](events, results.objects, **kwargs)

    events = self[process_ids](events, **kwargs)
    events = self[category_ids](events, **kwargs)     
    #events = set_ak_column(events, 'category_ids', ak.ones_like(events.event, dtype=np.uint8))
    
    events, results = self[custom_increment_stats]( 
        events,
        results,
        stats,
    )

    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events, results, self.config_inst.x.triggers)
        
    return events, results
