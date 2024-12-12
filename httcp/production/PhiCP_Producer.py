# coding: utf-8
"""
A producer to create Acoplanarity angles for differennt methods
"""

import os
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from httcp.production.PhiCP_Estimator import GetPhiCP
from httcp.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional

import law
np = maybe_import("numpy")
ak = maybe_import("awkward")

logger = law.logger.get_logger(__name__)
#from IPython import embed



@producer(
    uses={
        "channel_id",
        "hcand.pt","hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode", "hcand.charge",
        "hcand.IPx", "hcand.IPy", "hcand.IPz",
        "hcandprod.pt","hcandprod.eta", "hcandprod.phi", "hcandprod.mass", "hcandprod.pdgId",
        optional("GenTau.pt"), optional("GenTau.eta"), optional("GenTau.phi"), optional("GenTau.mass"),
        optional("GenTau.IPx"), optional("GenTau.IPy"), optional("GenTau.IPz"),
        optional("GenTau.decayMode"), optional("GenTau.charge"), 
        optional("GenTauProd.pt"), optional("GenTauProd.eta"), optional("GenTauProd.phi"), optional("GenTauProd.mass"),
    },
)
def ProducePhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> tuple[ak.Array,ak.Array,ak.Array,ak.Array,ak.Array]:

    PrepareP4 = lambda p4dict, mask : {key: ak.where(mask, val, val[:,:0]) for key, val in p4dict.items()}

    is_e      = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -1, axis=1), False)
    is_mu     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -2, axis=1), False)
    is_pi     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  0, axis=1), False)
    is_rho    = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  1, axis=1), False)
    is_a1DM2  = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  2, axis=1), False)
    is_a1DM10 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 10, axis=1), False)
    is_a1DM11 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 11, axis=1), False)
    is_a1     = is_a1DM10 # for now

    p4h1      = p4hcandinfo["p4h1"]
    p4h2      = p4hcandinfo["p4h2"]

    # etau
    mask_e_pi      = is_e(p4h1)   &   is_pi(p4h2)
    mask_e_rho     = is_e(p4h1)   &   is_rho(p4h2)
    mask_e_a1DM2   = is_e(p4h1)   &   is_a1DM2(p4h2)
    mask_e_a1DM10  = is_e(p4h1)   &   is_a1DM10(p4h2)
    mask_e_a1DM11  = is_e(p4h1)   &   is_a1DM11(p4h2)


    # mutau
    mask_mu_pi     = is_mu(p4h1)  &   is_pi(p4h2)
    mask_mu_rho    = is_mu(p4h1)  &   is_rho(p4h2)
    mask_mu_a1DM2  = is_mu(p4h1)  &   is_a1DM2(p4h2)
    mask_mu_a1DM10 = is_mu(p4h1)  &   is_a1DM10(p4h2)
    mask_mu_a1DM11 = is_mu(p4h1)  &   is_a1DM11(p4h2)    


    # tautau
    mask_pi_pi   = is_pi(p4h1)  &   is_pi(p4h2)

    mask_pi_rho  = is_pi(p4h1)  &   is_rho(p4h2)
    mask_rho_pi  = is_rho(p4h1) &   is_pi(p4h2)

    mask_pi_a1DM2   = is_pi(p4h1)  &   is_a1DM2(p4h2)
    mask_pi_a1DM10  = is_pi(p4h1)  &   is_a1DM10(p4h2)
    mask_pi_a1DM11  = is_pi(p4h1)  &   is_a1DM11(p4h2)
    mask_a1DM2_pi   = is_a1DM2(p4h1)  &   is_pi(p4h2)
    mask_a1DM10_pi  = is_a1DM10(p4h1) &   is_pi(p4h2)
    mask_a1DM11_pi  = is_a1DM11(p4h1) &   is_pi(p4h2)

    mask_rho_rho    = is_rho(p4h1) &   is_rho(p4h2)

    mask_rho_a1DM2  = is_rho(p4h1) &   is_a1DM2(p4h2)
    mask_rho_a1DM10 = is_rho(p4h1) &   is_a1DM10(p4h2)
    mask_rho_a1DM11 = is_rho(p4h1) &   is_a1DM11(p4h2)
    mask_a1DM2_rho  = is_a1DM2(p4h1)   &   is_rho(p4h2)
    mask_a1DM10_rho = is_a1DM10(p4h1)  &   is_rho(p4h2)
    mask_a1DM11_rho = is_a1DM11(p4h1)  &   is_rho(p4h2)
    
    mask_a1DM2_a1DM2    = is_a1DM2(p4h1)  &  is_a1DM2(p4h2)
    mask_a1DM2_a1DM10   = is_a1DM2(p4h1)  &  is_a1DM10(p4h2)
    mask_a1DM10_a1DM2   = is_a1DM10(p4h1) &  is_a1DM2(p4h2)
    mask_a1DM2_a1DM11   = is_a1DM2(p4h1)  &  is_a1DM11(p4h2)
    mask_a1DM11_a1DM2   = is_a1DM11(p4h1) &  is_a1DM2(p4h2)    

    mask_a1DM10_a1DM10   = is_a1DM10(p4h1) &  is_a1DM10(p4h2)
    mask_a1DM10_a1DM11   = is_a1DM10(p4h1) &  is_a1DM11(p4h2)
    mask_a1DM11_a1DM10   = is_a1DM11(p4h1) &  is_a1DM10(p4h2)

    mask_a1DM11_a1DM11   = is_a1DM11(p4h1) &  is_a1DM11(p4h2)
    
    # GetPhiCP
    # Inputs:
    #   masked dict of p4hcand
    #   leg1_method, leg2_method
    #   leg1_mode,lmeg2_mode

    dummyPhiCP = ak.values_astype(events.event[:,None][:,:0], np.float32)



    
    ##################################################################################
    #                                      IPIP                                      #
    #                        Only possible for (e/mu/pi)-pi                          #
    ##################################################################################
    logger.info("IP-IP")
    
    PhiCP_IPIP = ak.where(mask_e_pi,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_pi   ), "IP", "IP", "e",   "pi" ), dummyPhiCP)
    PhiCP_IPIP = ak.where(mask_mu_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_pi  ), "IP", "IP", "mu",  "pi" ), PhiCP_IPIP)
    PhiCP_IPIP = ak.where(mask_pi_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_pi  ), "IP", "IP", "pi",  "pi" ), PhiCP_IPIP)



    
    ##################################################################################
    #                                       DPDP                                     #
    #                    Only possible for (rho/a1)-(rho/a1) channel                 #
    ##################################################################################
    logger.info("DP-DP")
    
    # -- rho-rho
    PhiCP_DPDP = ak.where(mask_rho_rho,    GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_rho    ), "DP", "DP", "rho", "rho"), dummyPhiCP)
    # -- rho-a1
    PhiCP_DPDP = ak.where(mask_rho_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM2  ), "DP", "DP", "rho", "rho"), PhiCP_DPDP) # DM2 -> DM1
    PhiCP_DPDP = ak.where(mask_rho_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM10 ), "DP", "DP", "rho", "a1" ), PhiCP_DPDP)
    PhiCP_DPDP = ak.where(mask_rho_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM11 ), "DP", "DP", "rho", "a1" ), PhiCP_DPDP) # DM11 -> DM10
    # -- a1-rho
    PhiCP_DPDP = ak.where(mask_a1DM2_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_rho  ), "DP", "DP", "rho", "rho"), PhiCP_DPDP) # DM2 -> DM1 
    PhiCP_DPDP = ak.where(mask_a1DM10_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_rho ), "DP", "DP", "a1",  "rho"), PhiCP_DPDP) 
    PhiCP_DPDP = ak.where(mask_a1DM11_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_rho ), "DP", "DP", "a1",  "rho"), PhiCP_DPDP) # DM11 -> DM10
    # -- a1-a1
    PhiCP_DPDP = ak.where(mask_a1DM2_a1DM2,   GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM2  ), "DP", "DP", "rho", "rho" ), PhiCP_DPDP) # DM2 -> DM1 
    PhiCP_DPDP = ak.where(mask_a1DM2_a1DM10,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM10 ), "DP", "DP", "rho", "a1"  ), PhiCP_DPDP) # DM2 -> DM1
    PhiCP_DPDP = ak.where(mask_a1DM10_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM2 ), "DP", "DP", "a1", "rho"  ), PhiCP_DPDP) # DM2 -> DM1 
    PhiCP_DPDP = ak.where(mask_a1DM2_a1DM11,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM11 ), "DP", "DP", "rho",  "a1" ), PhiCP_DPDP) # DM2 -> DM1, DM11 -> DM10
    PhiCP_DPDP = ak.where(mask_a1DM11_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM2 ), "DP", "DP", "a1",  "rho" ), PhiCP_DPDP) # DM2 -> DM1, DM11 -> DM10 
    PhiCP_DPDP = ak.where(mask_a1DM10_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM10), "DP", "DP", "a1",   "a1" ), PhiCP_DPDP) 
    PhiCP_DPDP = ak.where(mask_a1DM10_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM11), "DP", "DP", "a1",   "a1" ), PhiCP_DPDP) # DM11 -> DM10
    PhiCP_DPDP = ak.where(mask_a1DM11_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM10), "DP", "DP", "a1",   "a1" ), PhiCP_DPDP) # DM11 -> DM10
    PhiCP_DPDP = ak.where(mask_a1DM11_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM11), "DP", "DP", "a1",   "a1" ), PhiCP_DPDP) # DM11 -> DM10
    

    

    ##################################################################################
    #                                      PVPVP                                     #
    #                         Only possible for tau-tau channel                      #
    ##################################################################################
    logger.info("PV-PV")
    
    # -- pi-pi
    PhiCP_PVPV = ak.where(mask_pi_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_pi  ), "PV", "PV", "pi",  "pi" ), dummyPhiCP)
    # -- pi-rho
    PhiCP_PVPV = ak.where(mask_pi_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho ), "PV", "PV", "pi",  "rho"), PhiCP_PVPV)
    # -- rho-pi
    PhiCP_PVPV = ak.where(mask_rho_pi,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_pi ), "PV", "PV", "rho", "pi" ), PhiCP_PVPV)
    # -- rho-rho
    PhiCP_PVPV = ak.where(mask_rho_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_rho), "PV", "PV", "rho", "rho"), PhiCP_PVPV)
    # -- pi-a1
    PhiCP_PVPV = ak.where(mask_pi_a1DM2,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM2  ), "PV", "PV", "pi",  "rho"), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_pi_a1DM10,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM10 ), "PV", "PV", "pi",  "a1" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_pi_a1DM11,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM11 ), "PV", "PV", "pi",  "a1" ), PhiCP_PVPV) # DM11 -> DM10
    # -- a1-pi
    PhiCP_PVPV = ak.where(mask_a1DM2_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_pi  ), "PV", "PV", "rho", "pi" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM10_pi,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_pi ), "PV", "PV", "a1",  "pi" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_a1DM11_pi,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_pi ), "PV", "PV", "a1",  "pi" ), PhiCP_PVPV) # DM11 -> DM10
    # -- rho-a1
    PhiCP_PVPV = ak.where(mask_rho_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM2 ), "PV", "PV", "rho", "rho"), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_rho_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM10), "PV", "PV", "rho", "a1" ),  PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_rho_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1DM11), "PV", "PV", "rho", "a1" ),  PhiCP_PVPV) # DM11 -> DM10
    # -- a1-rho
    PhiCP_PVPV = ak.where(mask_a1DM2_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_rho ), "PV", "PV", "rho", "rho"), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM10_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_rho), "PV", "PV", "a1",  "rho"), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_a1DM11_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_rho), "PV", "PV", "a1",  "rho"), PhiCP_PVPV) # DM11 -> DM10
    # -- a1-a1
    PhiCP_PVPV = ak.where(mask_a1DM2_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM2 ), "PV", "PV", "rho", "rho" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM2_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM10), "PV", "PV", "rho",  "a1" ), PhiCP_PVPV) # DM11 -> DM10
    PhiCP_PVPV = ak.where(mask_a1DM10_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM2), "PV", "PV", "a1",  "rho" ), PhiCP_PVPV) # DM11 -> DM10
    PhiCP_PVPV = ak.where(mask_a1DM2_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM2_a1DM11), "PV", "PV", "rho",  "a1" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM11_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM2), "PV", "PV", "a1",  "rho" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM10_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM10),"PV", "PV",  "a1",  "a1" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM10_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM10_a1DM11),"PV", "PV",  "a1",  "a1" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM11_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM10),"PV", "PV",  "a1",  "a1" ), PhiCP_PVPV) # DM2 -> DM1
    PhiCP_PVPV = ak.where(mask_a1DM11_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_a1DM11_a1DM11),"PV", "PV",  "a1",  "a1" ), PhiCP_PVPV) # DM2 -> DM1
    



    ##################################################################################
    #                                      IPDP                                      #
    #                   Only possible for (e/mu/pi)-(rho/a1) channel                 #
    ##################################################################################
    logger.info("IP-DP")
    
    # -- e-rho
    PhiCP_IPDP = ak.where(mask_e_rho,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_rho    ), "IP", "DP", "e",   "rho"), dummyPhiCP)
    # -- e-a1
    PhiCP_IPDP = ak.where(mask_e_a1DM2,  GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM2  ), "IP", "DP", "e",   "rho"), PhiCP_IPDP) # treat DM2 as rho
    PhiCP_IPDP = ak.where(mask_e_a1DM10, GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM10 ), "IP", "DP", "e",   "a1" ), PhiCP_IPDP) 
    PhiCP_IPDP = ak.where(mask_e_a1DM11, GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM11 ), "IP", "DP", "e",   "a1" ), PhiCP_IPDP) # treat DM11 as a1
    # -- mu-rho
    PhiCP_IPDP = ak.where(mask_mu_rho,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_rho   ), "IP", "DP", "mu",  "rho"), PhiCP_IPDP)
    # -- mu-a1
    PhiCP_IPDP = ak.where(mask_mu_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM2 ), "IP", "DP", "mu",  "rho"), PhiCP_IPDP) # treat DM2 as rho
    PhiCP_IPDP = ak.where(mask_mu_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM10), "IP", "DP", "mu",  "a1" ), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_mu_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM11), "IP", "DP", "mu",  "a1" ), PhiCP_IPDP) # treat DM11 as a1
    # -- pi-rho
    PhiCP_IPDP = ak.where(mask_pi_rho,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho   ), "IP", "DP", "pi",  "rho"), PhiCP_IPDP)
    # -- pi-a1
    PhiCP_IPDP = ak.where(mask_pi_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM2 ), "IP", "DP", "pi",  "rho"), PhiCP_IPDP) # treat DM2 as rho  
    PhiCP_IPDP = ak.where(mask_pi_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM10), "IP", "DP", "pi",  "a1" ), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_pi_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM11), "IP", "DP", "pi",  "a1" ), PhiCP_IPDP) # treat DM11 as a1 



    
    ##################################################################################
    #                                      IPPV                                      #
    #                Only possible for (e/mu/pi)-(pi/rho/a1) channel                 #
    ##################################################################################
    logger.info("IP-PV")
    
    # -- e-pi
    PhiCP_IPPV = ak.where(mask_e_pi,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_pi    ), "IP", "PV", "e",  "pi"),   dummyPhiCP)    
    # -- e-rho
    PhiCP_IPPV = ak.where(mask_e_rho,   GetPhiCP(PrepareP4(p4hcandinfo, mask_e_rho   ), "IP", "PV", "e",   "rho"), PhiCP_IPPV)
    # -- e-a1
    PhiCP_IPPV = ak.where(mask_e_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM2 ), "IP", "PV", "e",   "rho"), PhiCP_IPPV) # treat DM2 as rho
    PhiCP_IPPV = ak.where(mask_e_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM10), "IP", "PV", "e",   "a1" ), PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_e_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1DM11), "IP", "PV", "e",   "a1" ), PhiCP_IPPV) # treat DM11 as a1 
    # -- mu-pi
    PhiCP_IPPV = ak.where(mask_mu_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_pi   ), "IP", "PV", "mu",  "pi"),  PhiCP_IPPV)
    # -- mu-rho
    PhiCP_IPPV = ak.where(mask_mu_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_rho ), "IP", "PV", "mu",  "rho"), PhiCP_IPPV)
    # -- mu-a1
    PhiCP_IPPV = ak.where(mask_mu_a1DM2, GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM2 ),"IP", "PV", "mu",  "rho"), PhiCP_IPPV) # treat DM2 as rho
    PhiCP_IPPV = ak.where(mask_mu_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM10),"IP", "PV", "mu",  "a1"),  PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_mu_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1DM11),"IP", "PV", "mu",  "a1"),  PhiCP_IPPV) # treat DM11 as a1
    # -- pi-pi
    PhiCP_IPPV = ak.where(mask_pi_pi,    GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_pi ), "IP", "PV", "pi",  "pi"), PhiCP_IPPV)        
    # -- pi-rho
    PhiCP_IPPV = ak.where(mask_pi_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho ), "IP", "PV", "pi",  "rho"), PhiCP_IPPV)
    # -- pi-a1
    PhiCP_IPPV = ak.where(mask_pi_a1DM2,GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM2  ), "IP", "PV", "pi",  "rho" ), PhiCP_IPPV) # treat DM2 as rho
    PhiCP_IPPV = ak.where(mask_pi_a1DM10,GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM10), "IP", "PV", "pi",  "a1" ),  PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_pi_a1DM11,GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1DM11), "IP", "PV", "pi",  "a1" ),  PhiCP_IPPV) # treat DM11 as a1

    
    
    return events, PhiCP_IPIP, PhiCP_DPDP, PhiCP_PVPV, PhiCP_IPDP, PhiCP_IPPV




# ------------ DETECTOR LEVEL ----------- #
@producer(
    uses={
        ProducePhiCP, ProduceDetCosPsi,
    },
    produces={
        ProduceDetCosPsi,
        "PhiCP_IPIP",
        "PhiCP_IPIP_alpha_lt_piby4",
        "PhiCP_IPIP_alpha_gt_piby4",
        "PhiCP_DPDP",
        "PhiCP_PVPV",
        "PhiCP_IPDP",
        "PhiCP_IPDP_alpha_lt_piby4",
        "PhiCP_IPDP_alpha_gt_piby4",
        "PhiCP_IPPV",
    },
)
def ProduceDetPhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    events, alpha = self[ProduceDetCosPsi](events, p4hcandinfo)

    print(f" ===>>> ProduceDetPhiCP ===>>> ")
    with_alpha_lt_piby4 = ak.fill_none(ak.firsts(alpha < np.pi/4.0, axis=1), False)
    with_alpha_gt_piby4 = ak.fill_none(ak.firsts(alpha > np.pi/4.0, axis=1), False)

    events, PhiCP_IPIP, PhiCP_DPDP, PhiCP_PVPV, PhiCP_IPDP, PhiCP_IPPV = self[ProducePhiCP](events, p4hcandinfo)

    dummy = PhiCP_IPIP[:,:0]
    
    #psi_IPIP = ak.where(ak.num(PhiCP_IPIP, axis=1) == 1, psi, dummy)
    PhiCP_IPIP_alpha_lt_piby4 = ak.where(with_alpha_lt_piby4, PhiCP_IPIP, dummy)
    PhiCP_IPIP_alpha_gt_piby4 = ak.where(with_alpha_gt_piby4, PhiCP_IPIP, dummy)

    #psi_IPDP = ak.where(ak.num(PhiCP_IPDP, axis=1) == 1, psi, dummy)
    PhiCP_IPDP_alpha_lt_piby4 = ak.where(with_alpha_lt_piby4, PhiCP_IPDP, dummy)
    PhiCP_IPDP_alpha_gt_piby4 = ak.where(with_alpha_gt_piby4, PhiCP_IPDP, dummy)
    
    #PhiCP_PVPV = ak.nan_to_num(PhiCP_PVPV, nan=0.0)
    PhiCP_IPIP = ak.nan_to_num(PhiCP_IPIP, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPIP_alpha_lt_piby4 = ak.nan_to_num(PhiCP_IPIP_alpha_lt_piby4, 0.0)
    PhiCP_IPIP_alpha_gt_piby4 = ak.nan_to_num(PhiCP_IPIP_alpha_gt_piby4, 0.0)
    PhiCP_PVPV = ak.nan_to_num(PhiCP_PVPV, 0.0)
    PhiCP_DPDP = ak.nan_to_num(PhiCP_DPDP, 0.0)
    PhiCP_IPPV = ak.nan_to_num(PhiCP_IPPV, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPDP = ak.nan_to_num(PhiCP_IPDP, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPDP_alpha_lt_piby4 = ak.nan_to_num(PhiCP_IPDP_alpha_lt_piby4, 0.0)
    PhiCP_IPDP_alpha_gt_piby4 = ak.nan_to_num(PhiCP_IPDP_alpha_gt_piby4, 0.0)
    
    events = set_ak_column(events, "PhiCP_IPIP", PhiCP_IPIP)
    events = set_ak_column(events, "PhiCP_IPIP_alpha_lt_piby4", PhiCP_IPIP_alpha_lt_piby4)
    events = set_ak_column(events, "PhiCP_IPIP_alpha_gt_piby4", PhiCP_IPIP_alpha_gt_piby4)
    events = set_ak_column(events, "PhiCP_DPDP", PhiCP_DPDP)
    events = set_ak_column(events, "PhiCP_PVPV", PhiCP_PVPV)
    events = set_ak_column(events, "PhiCP_IPDP", PhiCP_IPDP)
    events = set_ak_column(events, "PhiCP_IPDP_alpha_lt_piby4", PhiCP_IPDP_alpha_lt_piby4)
    events = set_ak_column(events, "PhiCP_IPDP_alpha_gt_piby4", PhiCP_IPDP_alpha_gt_piby4)
    events = set_ak_column(events, "PhiCP_IPPV", PhiCP_IPPV)

    return events



# ------------ GENERATOR LEVEL ----------- #
@producer(
    uses={
        ProducePhiCP, ProduceGenCosPsi
    },
    produces={
        ProduceGenCosPsi,
        "PhiCPGen_IPIP",
        "PhiCPGen_IPIP_alpha_lt_piby4",
        "PhiCPGen_IPIP_alpha_gt_piby4",        
        "PhiCPGen_DPDP",
        "PhiCPGen_PVPV",
        "PhiCPGen_IPDP",
        "PhiCPGen_IPDP_alpha_lt_piby4",
        "PhiCPGen_IPDP_alpha_gt_piby4",        
        "PhiCPGen_IPPV",
    },
    mc_only=True,
)
def ProduceGenPhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    events, alpha = self[ProduceGenCosPsi](events, p4hcandinfo)

    print(f" ===>>> ProduceGenPhiCP ===>>> ")
    with_alpha_lt_piby4 = ak.fill_none(ak.firsts(alpha < np.pi/4.0, axis=1), False)
    with_alpha_gt_piby4 = ak.fill_none(ak.firsts(alpha > np.pi/4.0, axis=1), False)
    
    events, PhiCP_IPIP, PhiCP_DPDP, PhiCP_PVPV, PhiCP_IPDP, PhiCP_IPPV = self[ProducePhiCP](events, p4hcandinfo)

    dummy = PhiCP_IPIP[:,:0]

    #from IPython import embed; embed()
    
    #psi_IPIP = ak.where(ak.num(PhiCP_IPIP, axis=1) == 1, psi, dummy)
    PhiCP_IPIP_alpha_lt_piby4 = ak.where(with_alpha_lt_piby4, PhiCP_IPIP, dummy)
    PhiCP_IPIP_alpha_gt_piby4 = ak.where(with_alpha_gt_piby4, PhiCP_IPIP, dummy)

    #psi_IPDP = ak.where(ak.num(PhiCP_IPDP, axis=1) == 1, psi, dummy)
    PhiCP_IPDP_alpha_lt_piby4 = ak.where(with_alpha_lt_piby4, PhiCP_IPDP, dummy)
    PhiCP_IPDP_alpha_gt_piby4 = ak.where(with_alpha_gt_piby4, PhiCP_IPDP, dummy)
    

    PhiCP_IPIP = ak.nan_to_num(PhiCP_IPIP, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPIP_alpha_lt_piby4 = ak.nan_to_num(PhiCP_IPIP_alpha_lt_piby4, 0.0)
    PhiCP_IPIP_alpha_gt_piby4 = ak.nan_to_num(PhiCP_IPIP_alpha_gt_piby4, 0.0)
    PhiCP_PVPV = ak.nan_to_num(PhiCP_PVPV, 0.0)
    PhiCP_DPDP = ak.nan_to_num(PhiCP_DPDP, 0.0)
    PhiCP_IPPV = ak.nan_to_num(PhiCP_IPPV, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPDP = ak.nan_to_num(PhiCP_IPDP, 0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPDP_alpha_lt_piby4 = ak.nan_to_num(PhiCP_IPDP_alpha_lt_piby4, 0.0)
    PhiCP_IPDP_alpha_gt_piby4 = ak.nan_to_num(PhiCP_IPDP_alpha_gt_piby4, 0.0)
    
    
    events = set_ak_column(events, "PhiCPGen_IPIP", PhiCP_IPIP)
    events = set_ak_column(events, "PhiCPGen_IPIP_alpha_lt_piby4", PhiCP_IPIP_alpha_lt_piby4)
    events = set_ak_column(events, "PhiCPGen_IPIP_alpha_gt_piby4", PhiCP_IPIP_alpha_gt_piby4)
    events = set_ak_column(events, "PhiCPGen_DPDP", PhiCP_DPDP)
    events = set_ak_column(events, "PhiCPGen_PVPV", PhiCP_PVPV)
    events = set_ak_column(events, "PhiCPGen_IPDP", PhiCP_IPDP)
    events = set_ak_column(events, "PhiCPGen_IPDP_alpha_lt_piby4", PhiCP_IPDP_alpha_lt_piby4)
    events = set_ak_column(events, "PhiCPGen_IPDP_alpha_gt_piby4", PhiCP_IPDP_alpha_gt_piby4)
    events = set_ak_column(events, "PhiCPGen_IPPV", PhiCP_IPPV)

    return events
