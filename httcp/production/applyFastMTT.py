# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

from httcp.production.FastMTT import FastMTT

np = maybe_import("numpy")
ak = maybe_import("awkward")
pd = maybe_import("pandas")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


@producer(
    uses={
        # nano columns
        "hcand.*", "channel_id",
        "PuppiMET.pt", "PuppiMET.phi", "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
    },
    produces={
        # new columns
        "hcand.pt_fastMTT", "hcand.eta_fastMTT", "hcand.phi_fastMTT", "hcand.mass_fastMTT",
        "hcand_invm_fastMTT",
    },
)
def apply_fastMTT(
        self: Producer, 
        events: ak.Array,
        run_fmtt = True,
        **kwargs
) -> ak.Array:

    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hmass  = (hcand_[:,:1] + hcand_[:,1:2]).mass

    # dummy if fastMTT does not run
    hcand_pt_fastMTT = hcand_.pt
    hcand_eta_fastMTT = hcand_.eta
    hcand_phi_fastMTT = hcand_.phi
    hcand_mass_fastMTT = hcand_.mass
    mass_h = hmass * 1.2

    
    if run_fmtt:
        
        etau_id   = self.config_inst.get_channel("etau").id
        mutau_id  = self.config_inst.get_channel("mutau").id
        tautau_id = self.config_inst.get_channel("tautau").id
        
        # decay_type:
        #1 - TauToHad
        #2 - TauToElec
        #3 - TauToMu
        type_1 = ak.ones_like(hmass)
        type_2 = 2 * type_1
        type_3 = 3 * type_1
        
        pt1    = ak.to_numpy(events.hcand.pt[:,0:1])
        eta1   = ak.to_numpy(events.hcand.eta[:,0:1])
        phi1   = ak.to_numpy(events.hcand.phi[:,0:1])
        mass1  = ak.to_numpy(events.hcand.mass[:,0:1])
        dm1    = events.hcand.decayMode[:,0:1]
        dm1_dummy = ak.values_astype(-1 * ak.ones_like(dm1), np.int32)
        dm1    = ak.to_numpy(ak.where(events.channel_id < 4, dm1_dummy, dm1))
        type1  = ak.to_numpy(ak.where(events.channel_id == etau_id,
                                      type_2,
                                      ak.where(events.channel_id == mutau_id,
                                               type_3,
                                               type_1)
                                      )
                             )
        
        pt2    = ak.to_numpy(events.hcand.pt[:,1:2])
        eta2   = ak.to_numpy(events.hcand.eta[:,1:2])
        phi2   = ak.to_numpy(events.hcand.phi[:,1:2])
        mass2  = ak.to_numpy(events.hcand.mass[:,1:2])
        dm2    = ak.to_numpy(events.hcand.decayMode[:,1:2])
        type2  = ak.to_numpy(type_1)
        
        metpt  = ak.to_numpy(events.PuppiMET.pt[:,None])
        metphi = ak.to_numpy(events.PuppiMET.phi[:,None])
        metcovxx = ak.to_numpy(events.PuppiMET.covXX[:,None])
        metcovxy = ak.to_numpy(events.PuppiMET.covXY[:,None])
        metcovyy = ak.to_numpy(events.PuppiMET.covYY[:,None])
        
        
        higgs_mass = ak.to_numpy(hmass)
        
        
        # prepare inputs properly
        metcov = np.concatenate((metcovxx,metcovxy,metcovxy,metcovyy), axis=1)
        metcov = np.reshape(metcov, (metcov.shape[0], 2, 2))
        measuredTauLeptons = np.concatenate((dm1,
                                             pt1,
                                             eta1,
                                             phi1,
                                             mass1,
                                             type1,
                                             dm2,
                                             pt2,
                                             eta2,
                                             phi2,
                                             mass2,
                                             type2), axis=1)
        measuredTauLeptons = np.reshape(measuredTauLeptons, (measuredTauLeptons.shape[0], 2, 6))
        metx   = (metpt * np.cos(metphi)).reshape(-1)
        mety   = (metpt * np.sin(metphi)).reshape(-1)
        
        
        # Launch FastMTT
        fMTT = FastMTT()
        #You can choose to plot likelihood for one of the events. -1 means no plot.
        fMTT.WhichLikelihoodPlot = -1
        #You can also choose to calculate uncertainties by:
        fMTT.CalculateUncertainties = True
        
        fMTT.run(measuredTauLeptons, metx, mety, metcov)
        
        mass_h = fMTT.mass
        mass_h = ak.from_regular(ak.Array(mass_h.reshape(mass_h.shape[0],1)))
        
        p4_h1   = fMTT.tau1P4
        px_h1   = ak.from_regular(ak.Array(p4_h1[:,0:1]))
        py_h1   = ak.from_regular(ak.Array(p4_h1[:,1:2]))
        pz_h1   = ak.from_regular(ak.Array(p4_h1[:,2:3]))
        energy_h1 = ak.from_regular(ak.Array(p4_h1[:,3:4]))
        p4_h1_reg = ak.zip(
            {
                "x": px_h1, "y": py_h1, "z": pz_h1, "t": energy_h1,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
    
        p4_h2  = fMTT.tau2P4
        px_h2   = ak.from_regular(ak.Array(p4_h2[:,0:1]))
        py_h2   = ak.from_regular(ak.Array(p4_h2[:,1:2]))
        pz_h2   = ak.from_regular(ak.Array(p4_h2[:,2:3]))
        energy_h2 = ak.from_regular(ak.Array(p4_h2[:,3:4]))
        p4_h2_reg = ak.zip(
            {
                "x": px_h2, "y": py_h2, "z": pz_h2, "t": energy_h2,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

        hcand_pt_fastMTT  = ak.concatenate([p4_h1_reg.pt, p4_h2_reg.pt], axis=1)
        hcand_eta_fastMTT = ak.concatenate([p4_h1_reg.eta, p4_h2_reg.eta], axis=1)
        hcand_phi_fastMTT = ak.concatenate([p4_h1_reg.phi, p4_h2_reg.phi], axis=1)
        h1_mass = p4_h1_reg.mass
        #from IPython import embed; embed()
        mass_check = ak.nan_to_num(ak.firsts(h1_mass, axis=1), nan=-100.0)
        nan_mass_mask = mass_check < 0.0
        nan_mass_mask_et = nan_mass_mask[events.channel_id == etau_id]
        nan_mass_mask_mt = nan_mass_mask[events.channel_id == mutau_id]
        nan_mass_mask_tt = nan_mass_mask[events.channel_id == tautau_id]
        logger.info("Checking h1 mass")
        logger.warning("in general: check h1 regressed mass for leptonic channel, many nan values found for etau channel with WW dataset though")
        logger.warning(f"Out of {len(h1_mass)} events, {ak.sum(nan_mass_mask)} events have unphysical values for fastMTT h1 mass")
        logger.warning(f"Channel wise, e-tau has {ak.sum(nan_mass_mask_et)}, mu-tau has {ak.sum(nan_mass_mask_mt)}, and tau-tau has {ak.sum(nan_mass_mask_tt)} events with unphysical values")
        h1_mass = ak.where(events.channel_id == etau_id,
                           ak.nan_to_num(h1_mass, nan=0.00051),
                           ak.where(events.channel_id == mutau_id,
                                    ak.nan_to_num(h1_mass, nan=0.10566),
                                    h1_mass))
        
        h2_mass = p4_h2_reg.mass
        #from IPython import embed; embed()
        mass_check = ak.nan_to_num(ak.firsts(h2_mass, axis=1), nan=-100.0)
        nan_mass_mask = mass_check < -99.0
        nan_mass_mask_et = nan_mass_mask[events.channel_id == etau_id]
        nan_mass_mask_mt = nan_mass_mask[events.channel_id == mutau_id]
        nan_mass_mask_tt = nan_mass_mask[events.channel_id == tautau_id]
        logger.info("Checking h2 mass")
        logger.warning(f"Out of {len(h1_mass)} events, {ak.sum(nan_mass_mask)} events have unphysical values for fastMTT h2 mass")
        logger.warning(f"Channel wise, e-tau has {ak.sum(nan_mass_mask_et)}, mu-tau has {ak.sum(nan_mass_mask_mt)}, and tau-tau has {ak.sum(nan_mass_mask_tt)} events with unphysical values")
        
        hcand_mass_fastMTT = ak.concatenate([h1_mass, h2_mass], axis=1)
        
        # any nan fastMTT mass?
        mass_check = ak.nan_to_num(ak.firsts(mass_h, axis=1), nan=-100.0)
        nan_mass_mask = mass_check < -99.0
        nan_mass_mask_et = nan_mass_mask[events.channel_id == etau_id]
        nan_mass_mask_mt = nan_mass_mask[events.channel_id == mutau_id]
        nan_mass_mask_tt = nan_mass_mask[events.channel_id == tautau_id]
        logger.info("Checking fastMTT mass")
        logger.warning(f"Out of {len(mass_h)} events, {ak.sum(nan_mass_mask)} events have unphysical values for fastMTT mass")
        logger.warning(f"Channel wise, e-tau has {ak.sum(nan_mass_mask_et)}, mu-tau has {ak.sum(nan_mass_mask_mt)}, and tau-tau has {ak.sum(nan_mass_mask_tt)} events with unphysical values")

    else:
        logger.critical("Set [cfg.x.run_fastMTT = True] in the main config to make fastMTT run ... setting dummy variables as output")


    events = set_ak_column(events, "hcand.pt_fastMTT",   hcand_pt_fastMTT)
    events = set_ak_column(events, "hcand.eta_fastMTT",  hcand_eta_fastMTT)
    events = set_ak_column(events, "hcand.phi_fastMTT",  hcand_phi_fastMTT)
    events = set_ak_column(events, "hcand.mass_fastMTT", hcand_mass_fastMTT)        
    events = set_ak_column_f32(events, "hcand_invm_fastMTT", mass_h)
        
    return events


