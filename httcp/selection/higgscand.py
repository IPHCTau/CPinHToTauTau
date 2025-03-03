# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

import law
from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.util import enforce_hcand_type, IF_RUN2, IF_RUN3

from httcp.calibration.tau import insert_calibrated_taus
from httcp.production.ReArrangeHcandProds import getphotons, getpions, presel_decay_pis, presel_decay_pi0s, reconstructPi0

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

logger = law.logger.get_logger(__name__)


def get_energy_split(h, pions, pizeros):
    """
      energy split has appropriate values for DM 1 and 2
      for rest of the DMs, it would be empty
    """
    pi  = pions[:,:1]
    pi0 = pizeros[:,:1]

    has_both = (ak.num(pi, axis=1) == 1) & (ak.num(pi0, axis=1) == 1)
    
    dummy = pi.energy[:,:0]
    dm_mask = ak.fill_none(ak.firsts(((h.decayMode == 1) | (h.decayMode == 2)), axis=1), False)
    E_pi = ak.where(dm_mask & has_both, pi.energy, dummy)
    E_pi0 = ak.where(dm_mask & has_both, pi0.energy, dummy)
    
    E_split = np.abs((E_pi - E_pi0)/(E_pi + E_pi0))

    return E_split



@selector(
    uses={
        "channel_id",
        "Muon.*", "Electron.*", "Tau.*",
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
        hcand_array = enforce_hcand_type(hcand_pair_concat, 
                                         {"pt"            : "float64",
                                          "eta"           : "float64",
                                          "phi"           : "float64",
                                          "mass"          : "float64",
                                          "charge"        : "int32",
                                          "decayMode"     : "int32",
                                          "rawIdx"        : "int32",
                                          "IPx"           : "float64",
                                          "IPy"           : "float64",
                                          "IPz"           : "float64",
                                          "IPsig"         : "float64",
                                          "idVsJet"       : "int32",
                                          "genPartFlav"   : "int32",
                                          }
                                         )
    else:
        hcand_array = enforce_hcand_type(hcand_pair_concat, 
                                         {"pt"            : "float64",
                                          "eta"           : "float64",
                                          "phi"           : "float64",
                                          "mass"          : "float64",
                                          "charge"        : "int32",
                                          "decayMode"     : "int32",
                                          "rawIdx"        : "int32",
                                          "idVsJet"       : "int32",
                                          "genPartFlav"   : "int32",
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

    charge = ak.values_astype(charge, np.int32)
    mass   = ak.values_astype(mass, np.float32)
    
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
        "channel_id",
        "Electron.*","Muon.*","Tau.*","TauProd.*",
        assign_tauprod_mass_charge,
        insert_calibrated_taus,
    },
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass",
        "hcand.charge", "hcand.rawIdx", "hcand.decayMode",
        "hcand.energy_split",
        IF_RUN3("hcand.IPx", "hcand.IPy", "hcand.IPz"), "hcand.IPsig",
        "hcand.idVsJet", "hcand.genPartFlav",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass",
        "hcandprod.charge", "hcandprod.pdgId", "hcandprod.tauIdx",
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
    hcand    = hcand_array[:,0]
    hcand1 = hcand[:,0:1]
    hcand2 = hcand[:,1:2]
        
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

    # reconstruct pi-zeros here
    # save those in the hcand prods instead of photons

    hcand1prod_photons = getphotons(hcand1prods)
    hcand2prod_photons = getphotons(hcand2prods)
    
    hcand1prod_pions = getpions(hcand1prods)
    hcand2prod_pions = getpions(hcand2prods)

    # hcand1 and its decay products
    p4_hcand1     = ak.with_name(hcand1, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi  = ak.with_name(hcand1prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi  = presel_decay_pis(p4_hcand1, p4_hcand1_pi) # safe
    p4_hcand1_pi0 = reconstructPi0(p4_hcand1, hcand1prod_photons, method="simpleIC") # simpleIC, simpleMB
    p4_hcand1_pi0 = presel_decay_pi0s(p4_hcand1, p4_hcand1_pi0) # safe

    # hcand2 and its decay products
    p4_hcand2     = ak.with_name(hcand2, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi  = ak.with_name(hcand2prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi  = presel_decay_pis(p4_hcand2, p4_hcand2_pi)	# safe 
    p4_hcand2_pi0 = reconstructPi0(p4_hcand2, hcand2prod_photons, method="simpleIC") # simpleIC, simpleMB
    p4_hcand2_pi0 = presel_decay_pi0s(p4_hcand2, p4_hcand2_pi0)	# safe  

    # energy split
    hcand1_E_split = get_energy_split(p4_hcand1, p4_hcand1_pi, p4_hcand1_pi0)
    hcand1_pass_E_split_mask = ak.fill_none(ak.firsts(hcand1_E_split > 0.2, axis=1), True)
    hcand1_E_split = ak.fill_none(ak.firsts(hcand1_E_split, axis=1), -1.0)[:,None] # to wrap it in the hcand array

    hcand2_E_split = get_energy_split(p4_hcand2, p4_hcand2_pi, p4_hcand2_pi0)
    hcand2_pass_E_split_mask = ak.fill_none(ak.firsts(hcand2_E_split > 0.2, axis=1), True)
    hcand2_E_split = ak.fill_none(ak.firsts(hcand2_E_split, axis=1), -1.0)[:,None] # to wrap it in the hcand array

    hcand1prods = ak.concatenate([p4_hcand1_pi, p4_hcand1_pi0], axis=1)
    hcand2prods = ak.concatenate([p4_hcand2_pi, p4_hcand2_pi0], axis=1)

    hcand_prods = ak.concatenate([hcand1prods[:,None], hcand2prods[:,None]], axis=1)

    hcand_prods_array = enforce_hcand_type(ak.from_regular(hcand_prods),
                                           {"pt"            : "float64",
                                            "eta"           : "float64",
                                            "phi"           : "float64",
                                            "mass"          : "float64",
                                            "charge"        : "int32",
                                            "pdgId"         : "int64",
                                            "tauIdx"        : "int32"}
                                       )

    # saving hcand
    events = set_ak_column(events, "hcand", hcand)
    # saving hcand E split
    hcand_E_split = ak.from_regular(ak.concatenate([hcand1_E_split, hcand2_E_split], axis=1), axis=1)
    hcand_E_split_dummy = ak.from_regular(hcand_E_split[:,:0], axis=1)
    _mask = ak.num(events.hcand.decayMode, axis=1) == 2
    hcand_E_split = ak.where(_mask, hcand_E_split, hcand_E_split_dummy)    
    events = set_ak_column(events, "hcand.energy_split", hcand_E_split)

    # saving hcandprods
    events = set_ak_column(events, "hcandprod", hcand_prods_array)

    #FOR DEBUGGING
    if self.config_inst.x.verbose.selection.higgscand:
        for i in range(1000):
            if not events.channel_id[i] > 0: continue
            logger.info(f"channel_id: {events.channel_id[i]}")
            logger.info("hcand : [X,X], hcandprod : [ [y,z], [y,y] ]")
            logger.info(f"h DecayMode --- hprod pdgId --- hprod pt --- hcand mask")
            logger.info(f"{events.hcand.decayMode[i]} --- {events.hcandprod.pdgId[i]} --- {events.hcandprod.pt[i]} --- {hcand_prod_mask[i]}\n")

    #from IPython import embed; embed()
    
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

    result = SelectionResult(
        steps={
            "has_proper_tau_decay_products" : ak.sum(hcand_prod_mask, axis=1) == 2,
            "pass_energy_split_for_DM_1_2"  : (hcand1_pass_E_split_mask & hcand2_pass_E_split_mask),
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

    # to make channel specific
    if self.config_inst.x.is_channel_specific:
        ch_mask = events.event < 0 # all False
        for ch,mask in self.config_inst.x.channel_specific_info.items():
            if mask == True:
                logger.warning(f"Keeping events for {ch} channel only")
                ch_mask = ch_mask | (events.channel_id == self.config_inst.get_channel(ch).id) # False | (True/False)
        ch_mask_result = SelectionResult(
            steps = {
                "channel_mask" : ch_mask,
            },
        )
        result += ch_mask_result
            
    return events, result



