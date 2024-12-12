import functools

from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior
from columnflow.production.cms.electron import electron_weights

from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, layout_ak_array
from columnflow.columnar_util import optional_column as optional

from columnflow.util import maybe_import, safe_div, InsertableDict

from httcp.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)



# ------------------------------------------------- #
# Electron ID, ISO, Trigger weights
# ------------------------------------------------- #

electron_idiso_weights = electron_weights.derive(
    "electron_idiso_weights",
    cls_dict={
        "weight_name": "electron_idiso_weight",
        "get_electron_config": (lambda self: self.config_inst.x.electron_sf_names),
    }
)

# ################################################## #
# trig SF is from different file
# check: events fired by single ele trigger and
# check trigger_id too
# ################################################## #

@producer(
    uses={
        "Electron.pt",
        "Electron.eta",
        "trigger_ids",
        "single_e_triggered","cross_e_triggered",
    },
    produces={
        *[f"electron_Ele30_WPTight_trigger_weight{tag}" for tag in ["", "_up", "_down"]],
    },
    mc_only=True,
    # function to determine the correction file
    get_electron_file=(lambda self, external_files: external_files.electron_trig_sf),
    # function to determine the muon weight config
    get_electron_config=(lambda self: self.config_inst.x.electron_trig_sf_names),
    supported_versions=(1, 2),    
)
def electron_trigger_weights(self: Producer,
                             events: ak.Array,
                             **kwargs) -> ak.Array:
    #Producer that calculates the single lepton trigger weights.

    # get trigger ids for IsoMu24
    trigger_id_map = get_trigger_id_map(self.config_inst.x.triggers)
    trigger_id = trigger_id_map["HLT_Ele30_WPTight_Gsf"]
    
    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Electron.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt  = flat_np_view(events.Electron.pt, axis=1)
    eta = flat_np_view(events.Electron.eta, axis=1)
    
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    args = lambda mask, syst : (self.year, syst, self.wp, eta[mask], pt[mask])

    ele_mask = np.abs(events.Electron.pt) >= 31.0
    
    trig_mask = events.single_e_triggered & (events.trigger_ids == trigger_id) & ele_mask
    trig_mask = flat_np_view(trig_mask)

    for syst, postfix in [
            ("sf", ""),
            ("sfup", "_up"),
            ("sfdown", "_down"),
    ]:
        sf_values = sf_nom.copy()
        sf_values[trig_mask] = self.sf_corrector(*args(trig_mask, syst))

        events = set_ak_column(events,
                               f"electron_Ele30_WPTight_trigger_weight{postfix}",
                               reduce_mul(sf_values),
                               value_type=np.float32)

    return events

@electron_trigger_weights.requires
def electron_trigger_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@electron_trigger_weights.setup
def electron_trigger_weights_setup(
        self: Producer,
        reqs: dict,
        inputs: dict,
        reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.electron_trig_sf.load(formatter="gzip").decode("utf-8"),
    )
    corrector_name, self.year, self.wp = self.get_electron_config()
    self.sf_corrector = correction_set[corrector_name]


# ################################################## #
# xTrig SF is from different file
# check: events fired by cross ele trigger only and
# electron_pt < 30
# ################################################## #

@producer(
    uses={
        "Electron.pt", "Electron.eta",
        "single_e_triggered","cross_e_triggered",
    },
    produces={
        *[f"electron_xtrig_weight{tag}" for tag in ["", "_up", "_down"]],
    },
    mc_only=True,
    # function to determine the correction file
    get_electron_file=(lambda self, external_files: external_files.electron_xtrig_sf),
    # function to determine the muon weight config
    get_electron_config=(lambda self: self.config_inst.x.electron_xtrig_sf_names),
    supported_versions=(1, 2),
)
def electron_xtrigger_weights(self: Producer,
                              events: ak.Array,
                              **kwargs) -> ak.Array:
    # get xtrig sf
    # https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossMuTauHlt.json?ref_type=heads

    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Electron.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt  = flat_np_view(events.Electron.pt, axis=1)
    eta = flat_np_view(events.Electron.eta, axis=1)
    
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    args = lambda mask, syst : (self.year, syst, self.wp, eta[mask], pt[mask])

    xtrig_mask = events.cross_e_triggered & ~events.single_e_triggered & (events.Electron.pt < 30) & (np.abs(events.Electron.eta) <= 2.1)
    xtrig_mask = flat_np_view(xtrig_mask)
    
    for syst, postfix in [
            ("sf", ""),
            ("sfup", "_up"),
            ("sfdown", "_down"),
    ]:
        sf_values = sf_nom.copy()
        sf_values[xtrig_mask] = self.sf_corrector(*args(xtrig_mask, syst))

        events = set_ak_column(events,
                               f"electron_xtrig_weight{postfix}",
                               reduce_mul(sf_values),
                               value_type=np.float32)

    return events

@electron_xtrigger_weights.requires
def electron_xtrigger_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@electron_xtrigger_weights.setup
def electron_xtrigger_weights_setup(
        self: Producer,
        reqs: dict,
        inputs: dict,
        reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    #correction_set = correctionlib.CorrectionSet.from_string(
    #    bundle.files.electron_xtrig_sf.load(formatter="gzip").decode("utf-8"),
    #)
    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.electron_xtrig_sf.path,
    )    
    corrector_name, self.year, self.wp = self.get_electron_config()
    self.sf_corrector = correction_set[corrector_name]
