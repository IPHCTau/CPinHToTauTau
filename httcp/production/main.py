# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.normalization import stitched_normalization_weights

from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.pdf import pdf_weights
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.production.cms.mc_weight import mc_weight

from columnflow.production.util import attach_coffea_behavior

from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

from columnflow.config_util import get_events_from_categories
#from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP
#from httcp.production.weights import tauspinner_weight
from httcp.production.extra_weights import zpt_reweight, zpt_reweight_v2, ff_weight # ff_weight : dummy
from httcp.production.muon_weights import muon_id_weights, muon_iso_weights, muon_trigger_weights, muon_xtrigger_weights
from httcp.production.electron_weights import electron_idiso_weights, electron_trigger_weights, electron_xtrigger_weights
from httcp.production.tau_weights import tau_all_weights, tauspinner_weights


from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars
#from httcp.production.weights import pu_weight, muon_weight, tau_weight
#from httcp.production.weights import tau_weight
from httcp.production.sample_split import split_dy
from httcp.production.processes import build_abcd_masks

from httcp.production.columnvalid import make_column_valid

#from httcp.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi

from httcp.util import IF_DATASET_HAS_LHE_WEIGHTS, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL
from httcp.util import IF_RUN2, IF_RUN3, IF_ALLOW_STITCHING

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


@producer(
    uses={
        # nano columns
        "hcand.*", optional("GenTau.*"), optional("GenTauProd.*"),
        "Jet.pt",
        "PuppiMET.pt", "PuppiMET.phi",
        #reArrangeDecayProducts,
        #reArrangeGenDecayProducts,
        #ProduceGenPhiCP, ####ProduceGenCosPsi, 
        #ProduceDetPhiCP, ####ProduceDetCosPsi,
    },
    produces={
        # new columns
        "hcand_invm",
        "hcand_dr",
        "n_jet",
        #ProduceGenPhiCP, ####ProduceGenCosPsi,
        #ProduceDetPhiCP, ####ProduceDetCosPsi,
        "dphi_met_h1", "dphi_met_h2",
        "met_var_qcd_h1", "met_var_qcd_h2",
        "hT",
    },
)
def hcand_features(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]

    met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")

    mass = (hcand1 + hcand2).mass
    dr = hcand1.delta_r(hcand2)

    # deltaPhi between MET and hcand1 & 2
    dphi_met_h1 = met.delta_phi(hcand1)
    dphi_met_h2 = met.delta_phi(hcand2)
    met_var_qcd_h1 = met.pt * np.cos(dphi_met_h1)/hcand1.pt
    met_var_qcd_h2 = met.pt * np.cos(dphi_met_h2)/hcand2.pt

    # scalar sum pt  
    hT = ak.sum(events.Jet.pt, axis=1)
    #hT = ak.where(ak.num(events.Jet.pt, axis=1) > 0, hT, events.Jet.pt[:,:1])
    events = set_ak_column_f32(events, "hT",  hT)

    events = set_ak_column_f32(events, "hcand_invm",  mass)
    events = set_ak_column_f32(events, "hcand_dr",    dr)
    events = set_ak_column_f32(events, "dphi_met_h1", np.abs(dphi_met_h1))
    events = set_ak_column_f32(events, "dphi_met_h2", np.abs(dphi_met_h2))
    events = set_ak_column_f32(events, "met_var_qcd_h1", met_var_qcd_h1)
    events = set_ak_column_f32(events, "met_var_qcd_h2", met_var_qcd_h2)
    
    events = set_ak_column_i32(events, "n_jet", ak.num(events.Jet.pt, axis=1))

    # ########################### #
    # -------- For PhiCP -------- #
    # ########################### #
    """
    events, P4_dict = self[reArrangeDecayProducts](events)
    events   = self[ProduceDetPhiCP](events, P4_dict)
    #events  = self[ProduceDetCosPsi](events, P4_dict) # for CosPsi only
    
    if self.config_inst.x.extra_tags.genmatch:
        if "is_signal" in list(self.dataset_inst.aux.keys()):
            events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
            events = self[ProduceGenPhiCP](events, P4_gen_dict) 
            #events = self[ProduceGenCosPsi](events, P4_gen_dict) # for CosPsi only
    """
    # ########################### #
    
    return events


@producer(
    uses={
        make_column_valid,
        ##deterministic_seeds,
        attach_coffea_behavior,
        normalization_weights,
        "hcand.decayMode",
        IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- electron -- #
        electron_idiso_weights,
        electron_trigger_weights,
        electron_xtrigger_weights,
        # -- tau -- #
        tau_all_weights,
        IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_DY(zpt_reweight),
        #IF_DATASET_IS_DY(zpt_reweight_v2),
        hcand_features,
        hcand_mass,
        category_ids,
        build_abcd_masks,
        "channel_id",
        ff_weight,
    },
    produces={
        make_column_valid,
        ##deterministic_seeds,
        attach_coffea_behavior,
        normalization_weights,
        IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- electron -- #
        electron_idiso_weights,
        electron_trigger_weights,
        electron_xtrigger_weights,
        # -- tau -- #
        tau_all_weights,
        IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_DY(zpt_reweight),
        #IF_DATASET_IS_DY(zpt_reweight_v2),
        hcand_features,
        hcand_mass,
        #"channel_id",
        #"trigger_ids",
        category_ids,
        build_abcd_masks,
        ff_weight,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[make_column_valid](events, **kwargs)
    
    events = self[attach_coffea_behavior](events, **kwargs)
    # deterministic seeds
    ##events = self[deterministic_seeds](events, **kwargs)

    logger.warning("NO b-veto cut for tautau categories : Imperial")
    events = self[build_abcd_masks](events, **kwargs)
    # building category ids
    events, category_ids_debug_dict = self[category_ids](events, debug=False, **kwargs)

    events = self[hcand_features](events, **kwargs)       
    
    # debugging categories
    if self.config_inst.x.verbose.production.main:
        from httcp.production.debug import category_flow
        #print(category_ids_debug_dict)
        od_cats_in_config = self.config_inst.categories
        cat_etau = [cat for cat in od_cats_in_config if cat.name == "etau"][0]
        etau_events = get_events_from_categories(events, cat_etau)
        logger.info(f"Analysis : etau")
        category_flow("etau", etau_events)

        cat_mutau = [cat for cat in od_cats_in_config if cat.name == "mutau"][0]
        mutau_events = get_events_from_categories(events, cat_mutau)
        logger.info(f"Analysis : mutau")
        category_flow("mutau", mutau_events)

        cat_tautau = [cat for cat in od_cats_in_config if cat.name == "tautau"][0]
        tautau_events = get_events_from_categories(events, cat_tautau)
        logger.info(f"Analysis : tautau")
        category_flow("tautau", tautau_events)
        

    if self.dataset_inst.is_mc:
        # allow stitching is applicable only when datasets are DY or wjets, only if the stitching booleans are true in config
        allow_stitching = bool(ak.any([(self.dataset_inst.has_tag("is_dy") and self.config_inst.x.allow_dy_stitching),
                                       (self.dataset_inst.has_tag("is_w") and self.config_inst.x.allow_w_stitching)]))
        #from IPython import embed; embed()
        if allow_stitching:
            events = self[stitched_normalization_weights](events, **kwargs)
        else:
            events = self[normalization_weights](events, **kwargs)
        #events = self[stitched_normalization_weights](events, allow_stitching=allow_stitching, **kwargs)


        # TODO : pileup weight is constrained to max value 10
        # TODO : check columnflow production/pileup
        events = self[pu_weight](events, **kwargs)
        if self.has_dep(pdf_weights):
            events = self[pdf_weights](events, **kwargs)
        # ----------- Muon weights ----------- #
        events = self[muon_id_weights](events, **kwargs)
        events = self[muon_iso_weights](events, **kwargs)
        events = self[muon_trigger_weights](events, **kwargs)
        events = self[muon_xtrigger_weights](events, **kwargs)
        # ----------- Electron weights ----------- #
        events = self[electron_idiso_weights](events, **kwargs)
        events = self[electron_trigger_weights](events, **kwargs)
        events = self[electron_xtrigger_weights](events, **kwargs)
        # ----------- Tau weights ----------- #
        events = self[tau_all_weights](events, do_syst=True, **kwargs)

        #from IPython import embed; embed()
        if self.has_dep(tauspinner_weights):
            events = self[tauspinner_weights](events, **kwargs)

        # -- Z-pT reweighting with corrections from Imperial (Danny) and Kansas (Dennis : v2)
        if self.has_dep(zpt_reweight):
        #if self.has_dep(zpt_reweight_v2):
            events = self[zpt_reweight](events, **kwargs)
            #events = self[zpt_reweight_v2](events, **kwargs)

        #processes = self.dataset_inst.processes.names()
        #if ak.any(['dy_' in proc for proc in processes]):
        #    logger.info("splitting (any) Drell-Yan dataset ... ")
        #    events = self[split_dy](events,**kwargs)

    events = self[ff_weight](events, **kwargs)        

    # features
    events = self[hcand_mass](events, **kwargs)
    # events = self[mT](events, **kwargs)

    #from IPython import embed; embed()
    #events_cat = self[get_events_from_categories](events, ["tautau","real_1","hadD","has_1j","rho_1"], self.config_inst)
    
    
    return events
