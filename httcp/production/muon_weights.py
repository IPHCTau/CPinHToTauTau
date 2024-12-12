import functools

from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior
from columnflow.production.cms.muon import muon_weights

from columnflow.util import maybe_import, safe_div, InsertableDict

from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route, flat_np_view, layout_ak_array
from columnflow.columnar_util import optional_column as optional

from httcp.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)



# ------------------------------------------------- #
# Muon ID, ISO, Trigger weights
# ID, ISO, Trigger are from the same external json
# ------------------------------------------------- #

muon_id_weights = muon_weights.derive(
    "muon_id_weights",
    cls_dict={
        "weight_name": "muon_id_weight",
        "get_muon_config": (lambda self: self.config_inst.x.muon_id_sf_names),
    }
)

muon_iso_weights = muon_weights.derive(
    "muon_iso_weights",
    cls_dict={
        "weight_name": "muon_iso_weight",
        "get_muon_config": (lambda self: self.config_inst.x.muon_iso_sf_names),
    }
)

muon_IsoMu24_trigger_weights = muon_weights.derive(
    "muon_IsoMu24_trigger_weights",
    cls_dict={
        "weight_name": "muon_IsoMu24_trigger_weight",
        "get_muon_config": (lambda self: self.config_inst.x.muon_IsoMu24_trigger_sf_names),
    }
)

# ################################################## #
# Use muon_IsoMu24_trigger_weights in the func below
# to use the pt and trigger_id condition
# ################################################## #

@producer(
    uses={
        "trigger_ids",
        muon_IsoMu24_trigger_weights,
    },
    produces={
        *[f"muon_IsoMu24_trigger_weight{tag}" for tag in ["", "_up", "_down"]],
    },
)
def muon_trigger_weights(self: Producer,
                         events: ak.Array,
                         **kwargs) -> ak.Array:
    #Producer that calculates the single lepton trigger weights.

    # get trigger ids for IsoMu24
    trigger_id_map = get_trigger_id_map(self.config_inst.x.triggers)
    trigger_id = trigger_id_map["HLT_IsoMu24"]
    
    # compute muon trigger SF weights (NOTE: trigger SFs are only defined for muons with
    # pt > 26 GeV, so create a copy of the events array with with all muon pt < 26 GeV set to 26 GeV)
    trigger_sf_events = set_ak_column_f32(events, "Muon.pt", ak.where(events.Muon.pt > 26., events.Muon.pt, 26.))
    trigger_sf_events = self[muon_IsoMu24_trigger_weights](trigger_sf_events, **kwargs)
    for route in self[muon_IsoMu24_trigger_weights].produced_columns:
        events = set_ak_column_f32(events, route, ak.where(events.trigger_ids == trigger_id,
                                                           route.apply(trigger_sf_events),
                                                           1.0))
    # memory cleanup
    del trigger_sf_events

    return events


# ################################################## #
# xTrig SF is from different file
# check: events fired by cross mu trigger only and
# muon_pt < 25
# ################################################## #

@producer(
    uses={
        "Muon.pt", "Muon.eta",
        "single_mu_triggered","cross_mu_triggered",
    },
    produces={
        *[f"muon_xtrig_weight{tag}" for tag in ["", "_up", "_down"]],
    },
    mc_only=True,
    # function to determine the correction file
    get_muon_file=(lambda self, external_files: external_files.muon_xtrig_sf),
    # function to determine the muon weight config
    get_muon_config=(lambda self: self.config_inst.x.muon_xtrig_sf_names),
    supported_versions=(1, 2),
)
def muon_xtrigger_weights(self: Producer,
                          events: ak.Array,
                          **kwargs) -> ak.Array:
    # get xtrig sf
    # https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossMuTauHlt.json?ref_type=heads

    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Muon.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Muon.pt, axis=1)
    abseta = flat_np_view(abs(events.Muon.eta), axis=1)
    
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    args = lambda mask, syst : (abseta[mask], pt[mask], syst)

    xtrig_mask = events.cross_mu_triggered & ~events.single_mu_triggered & (events.Muon.pt < 25) & (np.abs(events.Muon.eta) <= 2.1)
    xtrig_mask = flat_np_view(xtrig_mask)
    
    for syst, postfix in [
            ("nominal", ""),
            ("systup", "_up"),
            ("systdown", "_down"),
    ]:
        sf_values = sf_nom.copy()
        sf_values[xtrig_mask] = self.sf_corrector(*args(xtrig_mask, syst))

        events = set_ak_column(events,
                               f"muon_xtrig_weight{postfix}",
                               reduce_mul(sf_values),
                               value_type=np.float32)

    return events

@muon_xtrigger_weights.requires
def muon_xtrigger_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@muon_xtrigger_weights.setup
def muon_xtrigger_weights_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    #json_data = bundle.files.muon_xtrig_sf.load()
    #correction_set = correctionlib.CorrectionSet.from_string(
    #    bundle.files.muon_xtrig_sf.load(),
    #)
    #correction_set = orrectionlib.CorrectionSet(json_data)
    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.muon_xtrig_sf.path,
    )
    corrector_name, self.year = self.get_muon_config()
    self.sf_corrector = correction_set[corrector_name]
