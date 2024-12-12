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



def get_sorted_pair(
        dtrpairs: ak.Array,
)->ak.Array:
    # redundant, because taus were sorted by the deeptau before
    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = (
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = (
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = (
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1)
    )

    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)

    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]


    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    dtrpair    = ak.concatenate([lep1, lep2], axis=1) 

    return dtrpair



@selector(
    uses={
        optional("Tau.pt"),
        optional("Tau.pt_tautau"),
        optional("Tau.mass"),
        optional("Tau.mass_tautau"),
        "Tau.eta", "Tau.phi",
        "Tau.rawIdx",
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:

    taus            = events.Tau[lep_indices]
    # Extra channel specific selections on tau
    # -------------------- #
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["tautau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["tautau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["tautau"]

    is_good_tau     = (
        (taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    if self.dataset_inst.is_mc:
        taus = ak.without_field(taus, "pt")
        taus = ak.with_field(taus, taus.pt_tautau, "pt")
        taus = ak.without_field(taus, "mass")
        taus = ak.with_field(taus, taus.mass_tautau, "mass")
    
    # -------------------- # 

    #from IPython import embed; embed()
    
    # Sorting leps [Tau] by deeptau [descending]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]

    leps_pair        = ak.combinations(taus, 2, axis=1)
    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)

    preselection = {
        "tautau_is_pt_40"      : (lep1.pt > 40) & (lep2.pt > 40),
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        "tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        "tautau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40, # invariant_mass(lep1, lep2) > 40
    }

    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["tautau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask
        
    leps_pair_sel = leps_pair[good_pair_mask]

    lep1, lep2  = ak.unzip(leps_pair)
    
    leps_pair_sel_single = ak.concatenate([lep1, lep2], axis=1)

    where_many   = ak.num(leps_pair_sel_single, axis=1) > 2

    pairs = ak.where(where_many, 
                     get_sorted_pair(leps_pair_sel),
                     leps_pair_sel_single)
    
    return SelectionResult(
        aux = pair_selection_steps,
    ), pairs
