# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.util import enforce_hcand_type, IF_RUN2, IF_RUN3

from httcp.calibration.tau import insert_calibrated_taus

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")



@selector(
    uses={
        "channel_id",
    },
    exposed=False,
)
def higgscand(
        self: Selector,
        events: ak.Array,
        hcand_pair: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:

    sel_hcand = ak.sum(ak.num(hcand_pair.pt, axis=-1), axis=1) == 2

    empty_hcand_pair = hcand_pair[:,:0][:,None]
    hcand_pair_concat = ak.where(events.channel_id == 1, hcand_pair[:,0][:,None], empty_hcand_pair)
    hcand_pair_concat = ak.where(events.channel_id == 2, hcand_pair[:,1][:,None], hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 4, hcand_pair[:,2][:,None], hcand_pair_concat)

    hcand_pair_concat = ak.where(events.channel_id == 3, 
                                 ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,1][:,None]], axis=1),
                                 hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 5,
                                 ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,2][:,None]], axis=1),
                                 hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 6, 
                                 ak.concatenate([hcand_pair[:,1][:,None], hcand_pair[:,2][:,None]], axis=1),
                                 hcand_pair_concat)

    hcand_array = None
    if self.config_inst.campaign.x.year >= 2022 :
        #from IPython import embed; embed()
        hcand_array = enforce_hcand_type(hcand_pair_concat, 
                                         {"pt"            : "float64",
                                          "eta"           : "float64",
                                          "phi"           : "float64",
                                          "mass"          : "float64",
                                          "charge"        : "int64",
                                          "decayMode"     : "int64",
                                          "rawIdx"        : "int64",
                                          "IPx"           : "float64",
                                          "IPy"           : "float64",
                                          "IPz"           : "float64",
                                          "IPsig"         : "float64",
                                          }
                                         )
    else:
        hcand_array = enforce_hcand_type(hcand_pair_concat, 
                                         {"pt"            : "float64",
                                          "eta"           : "float64",
                                          "phi"           : "float64",
                                          "mass"          : "float64",
                                          "charge"        : "int64",
                                          "decayMode"     : "int64",
                                          "rawIdx"        : "int64",
                                          }
                                         )
        

    sel_hcand = ak.fill_none(ak.num(ak.firsts(hcand_array.pt, axis=1), axis=1) == 2, False)

    return events, hcand_array, SelectionResult(
        steps={
            "One_higgs_cand_per_event": sel_hcand,
        },
    )


def select_tauprods(hcand_idx, tauprods):
    hcand_idx_brdcst, tauprod_tauIdx = ak.broadcast_arrays(ak.firsts(hcand_idx,axis=1), tauprods.tauIdx)
    hcandprod_mask                   = tauprod_tauIdx == hcand_idx_brdcst
    hcandprods                       = tauprods[hcandprod_mask]

    return hcandprods


@selector(
    uses={
        "TauProd.pdgId",
    },
    produces={
        "TauProd.mass", "TauProd.charge",
    },
    exposed=False,
)
def assign_tauprod_mass_charge(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    pionp  =  211
    pionm  = -211
    kaonp  =  321
    kaonm  = -321
    gamma  =  22

    mass = ak.where(np.abs(events.TauProd.pdgId) == pionp, 
                    0.13957, 
                    ak.where(np.abs(events.TauProd.pdgId) == kaonp,
                             0.493677,
                             ak.where(events.TauProd.pdgId == gamma,
                                      0.0, 0.0)) 
    )
    charge = ak.where(((events.TauProd.pdgId == pionp) | (events.TauProd.pdgId == kaonp)),
                      1.0,
                      ak.where(((events.TauProd.pdgId == pionm) | (events.TauProd.pdgId == kaonm)),
                               -1.0,
                               0.0)
                  )

    events = set_ak_column(events, "TauProd.mass", mass)
    events = set_ak_column(events, "TauProd.charge", charge)    

    return events

def build_hcand_mask(hcand, hcandprods, dummy):
    is_pion         = lambda prods : ((np.abs(prods.pdgId) == 211) | (np.abs(prods.pdgId) == 321))
    is_photon       = lambda prods : prods.pdgId == 22
    has_one_pion    = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 1)[:,None]
    has_atleast_one_pion = lambda prods : (ak.sum(is_pion(prods),   axis = 1) >= 1)[:,None] # new
    has_three_pions = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 3)[:,None]
    has_photons     = lambda prods : (ak.sum(is_photon(prods), axis = 1) >  0)[:,None]
    has_no_photons  = lambda prods : (ak.sum(is_photon(prods), axis = 1) == 0)[:,None]

    hcand_mask = ak.where(hcand.decayMode == 0,
                          has_atleast_one_pion(hcandprods),
                          ak.where(((hcand.decayMode == 1) | (hcand.decayMode == 2)),
                                   (has_atleast_one_pion(hcandprods) & has_photons(hcandprods)),
                                   ak.where(hcand.decayMode == 10,
                                            has_three_pions(hcandprods),
                                            ak.where(hcand.decayMode == 11,
                                                     (has_three_pions(hcandprods) & has_photons(hcandprods)),
                                                     dummy)
                                            )
                                   )
                          )
    return hcand_mask


@selector(
    uses={
        "Tau.*",
        "channel_id", "TauProd.*",
        assign_tauprod_mass_charge,
        insert_calibrated_taus,
    },
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.charge", "hcand.rawIdx", "hcand.decayMode",
        IF_RUN3("hcand.IPx", "hcand.IPy", "hcand.IPz"), "hcand.IPsig",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass", "hcandprod.charge", "hcandprod.pdgId", "hcandprod.tauIdx",
        assign_tauprod_mass_charge,
        insert_calibrated_taus,
    },
    exposed=False,
)
def higgscandprod(
        self: Selector,
        events: ak.Array,
        hcand_array: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id

    events   = self[assign_tauprod_mass_charge](events)

    tauprods = events.TauProd
    #hcand    = hcand_array
    hcand    = hcand_array[:,0]
    #hcand1 = ak.firsts(hcand[:,:,0:1], axis=1)
    #hcand2 = ak.firsts(hcand[:,:,1:2], axis=1)
    hcand1 = hcand[:,0:1]
    hcand2 = hcand[:,1:2]
    #hcand_concat = ak.concatenate([hcand1, hcand2], axis=1)
        
    hcand1_idx = hcand1.rawIdx
    hcand2_idx = hcand2.rawIdx

    hcand1prods = ak.where(events.channel_id == tautau_id,
                           select_tauprods(hcand1_idx, tauprods), 
                           tauprods[:,:0])
    hcand2prods = ak.where(((events.channel_id == etau_id) | (events.channel_id == mutau_id) | (events.channel_id == tautau_id)), 
                           select_tauprods(hcand2_idx, tauprods),
                           tauprods[:,:0])

    dummy = (events.event >= 0)[:,None]

    hcand1_mask = ak.fill_none(build_hcand_mask(hcand1, hcand1prods, dummy), False)
    hcand2_mask = ak.fill_none(build_hcand_mask(hcand2, hcand2prods, dummy), False)

    hcand_prod_mask = ak.concatenate([hcand1_mask, hcand2_mask], axis=1)
    
    hcand_prods = ak.concatenate([hcand1prods[:,None], hcand2prods[:,None]], axis=1)

    hcand_prods_array = enforce_hcand_type(ak.from_regular(hcand_prods),
                                           {"pt"            : "float64",
                                            "eta"           : "float64",
                                            "phi"           : "float64",
                                            "mass"          : "float64",
                                            "charge"        : "int64",
                                            "pdgId"         : "int64",
                                            "tauIdx"        : "int64"}
                                       )

    #FOR DEBUGGING
    if self.config_inst.x.verbose.selection.higgscand:
        for i in range(1000):
            if not events.channel_id[i] > 0: continue
            print(f"channel_id: {events.channel_id[i]}")
            print(f"h1 == DM : {hcand1.decayMode[i]}, PROD: {hcand1prods.pdgId[i]}, MASK: {hcand1_mask[i]}")
            print(f"h2 == DM : {hcand2.decayMode[i]}, PROD: {hcand2prods.pdgId[i]}, MASK: {hcand2_mask[i]}")
            print(f"Comb Mask: {hcand_prod_mask[i]}")
            print("\n")

    
    events = set_ak_column(events, "hcand",     hcand)
    events = set_ak_column(events, "hcandprod", hcand_prods_array)

    if self.dataset_inst.is_mc:
        events = self[insert_calibrated_taus](events)
    
    # set Muon, Electron and Tau from hcands
    hcand_idx = events.hcand.rawIdx
    raw_ele_idx_dummy = events.Electron.rawIdx[:,:0]
    raw_mu_idx_dummy  = events.Muon.rawIdx[:,:0]
    raw_tau_idx_dummy = events.Tau.rawIdx[:,:0]
    
    hcand_ele_idxs = ak.where(events.channel_id == etau_id, hcand_idx[:,0:1], raw_ele_idx_dummy)
    hcand_muo_idxs = ak.where(events.channel_id == mutau_id, hcand_idx[:,0:1], raw_mu_idx_dummy)
    hcand_tau_idxs = ak.where(events.channel_id == etau_id,
                              hcand_idx[:,1:2],
                              ak.where(events.channel_id == mutau_id,
                                       hcand_idx[:,1:2],
                                       ak.where(events.channel_id == tautau_id,
                                                hcand_idx,
                                                raw_tau_idx_dummy)
                                       )
                              )

    return events, SelectionResult(
        steps={
            "has_proper_tau_decay_products": ak.sum(hcand_prod_mask, axis=1) == 2,
        },
        objects={
            "Muon": {
                "Muon": hcand_muo_idxs,
            },
            "Electron": {
                "Electron": hcand_ele_idxs,
            },
            "Tau": {
                "Tau": hcand_tau_idxs,
            },
        },
    )
