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

from httcp.util import transverse_mass
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



def get_sorted_pair(
        dtrpairs: ak.Array,
)->ak.Array:

    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso03_all, ascending=True)

    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pfRelIso03_all[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pfRelIso03_all[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False
    )

    # if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if rawDeepTau2017v2p1VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    # Extract the first object in each pair (lep1) and the second object (lep2)
    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    # Concatenate lep1 and lep2 to create the final dtrpair
    dtrpair    = ak.concatenate([lep1, lep2], axis=1)

    return dtrpair




@selector(
    uses={
        # muon
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass",
        "Muon.charge", "Muon.pfRelIso03_all", "Muon.rawIdx",
        # tau
        optional("Tau.pt"),
        optional("Tau.pt_mutau"),
        "Tau.eta", "Tau.phi",
        optional("Tau.mass"),
        optional("Tau.mass_mutau"),
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet", "Tau.rawIdx",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
        # met
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
    },
    exposed=False,
)
def mutau_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:

    # lep1 and lep2 e.g.
    # lep1: [ [m1], [m1],    [m1,m2], [],   [m1,m2] ]
    # lep2: [ [t1], [t1,t2], [t1],    [t1], [t1,t2] ]
    muons = events.Muon[lep1_indices] 
    taus  = events.Tau[lep2_indices]

    # Extra channel specific selections on m or tau
    # -------------------- #
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["mutau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["mutau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["mutau"]
    
    is_good_tau     = (
        (taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    if self.dataset_inst.is_mc:
        taus = ak.without_field(taus, "pt")
        taus = ak.with_field(taus, taus.pt_mutau, "pt")
        taus = ak.without_field(taus, "mass")
        taus = ak.with_field(taus, taus.mass_mutau, "mass")

    # -------------------- # 
    
    met = events.MET if self.config_inst.campaign.x.year < 2022 else events.PuppiMET

    # Sorting lep1 [Electron] by isolation [ascending]
    muons_sort_idxs = ak.argsort(muons.pfRelIso03_all, axis=-1, ascending=True)
    muons = muons[muons_sort_idxs]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
        
    leps_pair        = ak.cartesian([muons, taus], axis=1)
    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)

    preselection = {
        "mutau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "mutau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        "mutau_mT_50"         : transverse_mass(lep1, met) < 50,
        "mutau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40,  # invariant_mass(lep1, lep2) > 40
    }

    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["mutau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask
        
    leps_pair_sel = leps_pair[good_pair_mask]

    lep1, lep2 = ak.unzip(leps_pair_sel)
    
    leps_pair_sel_single = ak.concatenate([lep1,lep2], axis=1)

    where_many   = ak.num(leps_pair_sel_single, axis=1) > 2
    
    pairs = ak.where(where_many,
                     get_sorted_pair(leps_pair_sel),
                     leps_pair_sel_single)

    return SelectionResult(
        aux = pair_selection_steps,
    ), pairs
