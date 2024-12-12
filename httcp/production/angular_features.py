# coding: utf-8
"""
A producer to create Acoplanarity angles for differennt methods
"""

import os
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from httcp.production.PhiCP_Estimator import GetPhiCP
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

#from IPython import embed


def get_rho_from_a1(p4hpi_m):
    pi0RecoM = 0.136 #approximate pi0 peak from fits in PF paper
    pi0RecoW = 0.013 #approximate width of pi0 peak from early in PF paper

    p4_pi_comb = ak.combinations(p4hpi_m, 2)

    pi1 = p4_pi_comb["0"]
    pi2 = p4_pi_comb["1"]

    mass = (pi1 + pi2).mass
    mass_sort_idx = ak.argsort(np.abs(mass-pi0RecoM), ascending=True)

    pi1 = pi1[mass_sort_idx]
    pi2 = pi2[mass_sort_idx]
    
    #mask = (pi1.charge * pi2.charge < 0) & (np.abs((pi1+pi2).mass - pi0RecoM) < 2*pi0RecoW)
    mask = (pi1.charge * pi2.charge < 0)
    
    pi1 = pi1[mask][:,:1]
    pi2 = pi2[mask][:,:1]
    
    rho = pi1+pi2
    
    return rho



@producer(
    uses={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode",
        "hcand.IPx", "hcand.IPy", "hcand.IPz", "hcand.charge",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass", "hcandprod.charge",
        "hcandprod.pdgId",
        optional("GenTau.pt"), optional("GenTau.eta"),
        optional("GenTau.phi"), optional("GenTau.mass"),
        optional("GenTau.decayMode"), optional("GenTau.charge"),
        optional("GenTau.IPx"), optional("GenTau.IPy"), optional("GenTau.IPz"),
        optional("GenTauProd.pt"), optional("GenTauProd.eta"),
        optional("GenTauProd.phi"),  optional("GenTauProd.mass"), optional("GenTauProd.charge"),
    },
)
def ProduceCosPsi(
        self: Producer,
        events: ak.Array,
        p4hcandinfodict: dict,
        **kwargs,
) -> tuple[ak.Array,ak.Array] :
    # https://github.com/DesyTau/DesyTauAnalysesUL/blob/bbHTT/BBHTT/interface/functionsCP.h#L72-L87
    # TLorentzVector momentum is the momentum of negative prong (negative muon in the case
    # of tau->mu+v+v decay and negative pion in the case of tau->pi+v or tau->rho+v decays)
    # TLorentzVector ipvec is IP vector of negative prong in the case of tau->mu+v+v or
    # tau->pi+v decay or momentum of neutral pion in the case of tau->rho+v decay 
    #    - ez is the unit vector in the direction of the proton beam (positive
    #        direction in the CMS coordinate frame)
    #    - p is the unit vector in the direction of the negative prong momentum
    #    - n is the unit vector in the direction of the IP vector of the
    #        negative prong (or pi0 momentum in the case when negative
    #        tau decays to rho+v)
    p4h1    = p4hcandinfodict["p4h1"]
    p4h1pi  = p4hcandinfodict["p4h1pi"]
    p4h1pi0 = p4hcandinfodict["p4h1pi0"]
    p4h2    = p4hcandinfodict["p4h2"]
    p4h2pi  = p4hcandinfodict["p4h2pi"]
    p4h2pi0 = p4hcandinfodict["p4h2pi0"]

    # to decide the negative prong decay
    h1_is_negative = ak.fill_none(ak.firsts(p4h1.charge < 0, axis=1), False)

    # select the e/mu/tauh, pions and pizeros for the negative leg
    p4h_m     = ak.where(h1_is_negative, p4h1, p4h2)
    p4hpi_m   = ak.where(h1_is_negative, p4h1pi, p4h2pi)
    p4hpi0_m  = ak.where(h1_is_negative, p4h1pi0, p4h2pi0)

    # create masks on the basis of decay modes
    is_dm_lep_pi = ak.fill_none(ak.firsts(p4h_m.decayMode <= 0,  axis=1), False) # e/mu/pi
    is_dm_rho    = ak.fill_none(ak.firsts(p4h_m.decayMode == 1,  axis=1), False) # rho
    is_dm_a1     = ak.fill_none(ak.firsts(p4h_m.decayMode == 10, axis=1), False) # a1
    
    # extract the info for the decay modes only
    # mandatory to make the cross and dot products successful !!!
    comb_mask = (is_dm_lep_pi | is_dm_rho | is_dm_a1)
    p4h_m = ak.where(comb_mask, p4h_m, p4h_m[:,:0])
    p4hpi_m = ak.where(comb_mask, p4hpi_m, p4hpi_m[:,:0])
    p4hpi0_m = ak.where(comb_mask, p4hpi0_m, p4hpi0_m[:,:0])
    
    # calculate the p4 of rho for the tau of decay mode 10 only
    p4hpi_m_for_a1 = ak.where(is_dm_a1, p4hpi_m, p4hpi_m[:,:0])
    rho = get_rho_from_a1(p4hpi_m_for_a1)

    dummy = rho.x[:,:0]
    
    nx = ak.where(is_dm_lep_pi, p4h_m.IPx, # for e/mu/pi : IP
                  ak.where(is_dm_rho, p4hpi0_m.x, # for rho: pi0
                           ak.where(is_dm_a1, rho.x, # for a1: rho
                                    dummy)))
    ny = ak.where(is_dm_lep_pi, p4h_m.IPy,
                  ak.where(is_dm_rho, p4hpi0_m.y,
                           ak.where(is_dm_a1, rho.y, dummy)))
    nz = ak.where(is_dm_lep_pi, p4h_m.IPz,
                  ak.where(is_dm_rho, p4hpi0_m.z,
                           ak.where(is_dm_a1, rho.z, dummy)))
    
    n = ak.zip({"x": nx, "y": ny, "z": nz, "t": ak.zeros_like(nx)},
               with_name = "LorentzVector",
               behavior = coffea.nanoevents.methods.vector.behavior)

    #zeros = ak.where(ak.num(dummy, axis=1) == 0, 1, 0)
    #ones  = ak.where(ak.num(dummy, axis=1) == 0, 0, 1)
    ez = ak.zip({"x": ak.zeros_like(p4h_m.x), "y": ak.zeros_like(p4h_m.y), "z": ak.ones_like(p4h_m.z), "t": ak.zeros_like(p4h_m.t)},
                with_name = "LorentzVector",
                behavior = coffea.nanoevents.methods.vector.behavior)

    #from IPython import embed; embed()

    # Get the unit vectors
    p  = p4h_m.pvec.unit
    n  = n.pvec.unit
    ez = ez.pvec.unit
    
    # calculate angle
    num    = (ez.cross(p)).dot(n.cross(p))
    den    = (ez.cross(p).absolute() * n.cross(p).absolute())
    # cospsi : always +ve
    cospsi = np.abs(num/den)
    psi    = np.arccos(cospsi)

    return events, psi

    

# ------------ DETECTOR LEVEL ----------- #
@producer(
    uses={
        ProduceCosPsi,
    },
    produces={
        "alpha",
    },
)
def ProduceDetCosPsi(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    print(" ===>>> ProduceDetCosPsi ===>>> ")
    events, psi = self[ProduceCosPsi](events, p4hcandinfo)
    psi = ak.nan_to_num(psi, 0.0)
    events = set_ak_column(events, "alpha", psi, value_type=np.float32)

    return events, psi



# ------------ DETECTOR LEVEL ----------- #
@producer(
    uses={
        ProduceCosPsi,
    },
    produces={
        "alphaGen",
    },
)
def ProduceGenCosPsi(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    print(" ===>>> ProduceGenCosPsi ===>>> ")
    events, psi = self[ProduceCosPsi](events, p4hcandinfo)
    psi = ak.nan_to_num(psi, 0.0)
    events = set_ak_column(events, "alphaGen", psi, value_type=np.float32)
    
    return events, psi
