import re
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from httcp.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


# ------------------------------------------------- #
# Assign MC weight [gen Weight / LHEWeight]
# ------------------------------------------------- #
@producer(
    uses={"genWeight", optional("LHEWeight.originalXWGTUP")},
    produces={"mc_weight"},
    mc_only=True,
)
def scale_mc_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Reads the genWeight and LHEWeight columns and makes a decision about which one to save. This
    should have been configured centrally [1] and stored in genWeight, but there are some samples
    where this failed.

    Strategy:

      1. Use LHEWeight.originalXWGTUP when it exists and genWeight is always 1.
      2. In all other cases, use genWeight.

    [1] https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD?rev=99#Weigths
    """
    # determine the mc_weight
    mc_weight = np.sign(events.genWeight)
    if has_ak_column(events, "LHEWeight.originalXWGTUP") and ak.all(events.genWeight == 1.0):
        mc_weight = np.sign(events.LHEWeight.originalXWGTUP)

    # store the column
    events = set_ak_column(events, "mc_weight", mc_weight, value_type=np.float32)

    return events


# ------------------------------------------------- #
# Calculate Zpt weights
# ------------------------------------------------- #

@producer(
    uses={
        "GenZ.pt", "GenZ.mass",
    },
    produces={
        "zpt_reweight", "zpt_reweight_up", "zpt_reweight_down",
    },
    mc_only=True,
)
def zpt_reweight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # if no gen particle found, all fields of GenZ will be zero
    # Those events will have nominal zpt rewt = 1.0
    # the same will be applied also for the evnts outside the range

    is_outside_range = (
        ((events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0))
        | ((events.GenZ.pt >= 600.0) | (events.GenZ.mass >= 1000.0))
    )

    # for safety
    zm  = ak.where(events.GenZ.mass > 1000.0, 999.99, events.GenZ.mass)
    zpt = ak.where(events.GenZ.pt > 600.0, 599.99, events.GenZ.pt)

    sf_nom = ak.where(is_outside_range,
                      1.0,
                      self.zpt_corrector.evaluate(zm,zpt))

    events = set_ak_column(events, "zpt_reweight", sf_nom, value_type=np.float32)

    sf_up   = ak.where(is_outside_range, 1.0, 1.67 * sf_nom)
    events = set_ak_column(events, "zpt_reweight_up", sf_up, value_type=np.float32)
    sf_down = ak.where(is_outside_range, 1.0, 0.33 * sf_nom)
    events = set_ak_column(events, "zpt_reweight_down", sf_down, value_type=np.float32)    
    
    return events


@zpt_reweight.requires
def zpt_reweight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@zpt_reweight.setup
def zpt_reweight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.zpt_rewt_v1_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["zptreweight"]



# ------------------------------------------------- #
# Calculate Zpt weights (Dennis)
# ------------------------------------------------- #

@producer(
    uses={
        "GenZ.pt", "GenZ.mass",
    },
    produces={
        "zpt_reweight", "zpt_reweight_up", "zpt_reweight_down",
    },
    mc_only=True,
)
def zpt_reweight_v2(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # if no gen particle found, all fields of GenZ will be zero
    # Those events will have nominal zpt rewt = 1.0
    # the same will be applied also for the evnts outside the range

    is_outside_range = (
        (events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0)
    )

    zpt = events.GenZ.pt

    dataset_type = "LO" if "madgraph" in self.dataset_inst.name else "NLO"
    era = f"{self.config_inst.campaign.x.year}{self.config_inst.campaign.x.postfix}_{dataset_type}"
    
    for syst in ["nom", "up1", "down1"]:
        tag = re.match(r'([a-zA-Z]+)\d*', syst).group(1)
        name = "zpt_reweight" if syst=="nom" else f"zpt_reweight_{tag}"
        sf = ak.where(is_outside_range,
                      1.0,
                      self.zpt_corrector.evaluate(era,
                                                  zpt,
                                                  syst)
                      )
        
        events = set_ak_column(events, name, sf, value_type=np.float32)
    
    return events


@zpt_reweight_v2.requires
def zpt_reweight_v2_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@zpt_reweight_v2.setup
def zpt_reweight_v2_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.zpt_rewt_v2_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["DY_pTll_reweighting"]






# ------------------------------------------------- #
# Calculate FF weights (Dummy)
# ------------------------------------------------- #

@producer(
    uses={
        "channel_id",
        "category_ids",
        "hcand.pt", "hcand.decayMode", "n_jet",
        "met_var_qcd_h1",
    },
    produces={
        "ff_weight",
        #"closure_weight",
        #"ff_ext_corr_weight",
    },
    mc_only=False,
)
def ff_weight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # Leading candidate
    hcand1 = events.hcand[:,0] 

    is_B = (
        (events.channel_id == 4)
        & ~events.is_os
        & events.is_real_1
        & events.is_iso_2
        & ~events.is_iso_1
    )
    is_C0 = (
        (events.channel_id == 4)
        & events.is_os
        & events.is_real_1
        & ~events.is_iso_2
        & ~events.is_iso_1
    )
    is_C = (
        (events.channel_id == 4)
        & events.is_os
        & events.is_real_1
        & events.is_iso_2
        & ~events.is_iso_1
    )

    is_B_category = ak.to_numpy(is_B)
    is_C0_category = ak.to_numpy(is_C0)
    is_C_category = ak.to_numpy(is_C)

    
    # Get the pt of the leading candidate
    pt1 = flat_np_view(hcand1.pt[:,None])
    metvarqcdh1 = flat_np_view(events.met_var_qcd_h1[:,None])
    
    # Get met_var_qcd_h1
    #met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")
    #dphi_met_h1 = met.delta_phi(hcand1)
    #met_var_qcd_h1 = met.pt * np.cos(dphi_met_h1)/hcand1.pt
    #met_var_qcd_h1 = flat_np_view(met_var_qcd_h1[:,None])

    njet = ak.where(events.n_jet > 2, 2, events.n_jet)
    njet = ak.to_numpy(njet)
    dm = ak.where(hcand1.decayMode < 0, 0, hcand1.decayMode)
    dm = ak.to_numpy(dm)
        
    fake_factors_nom = self.ff_corrector.evaluate(
        pt1,
        dm,
        njet,
        "nom",
    )
    #fake_0_factors_nom = self.ff0_corrector.evaluate(
    #    pt1,
    #    dm,
    #    njet
    #)
    #ext_corr_nom = self.ext_corrector.evaluate(
    #    metvarqcdh1,
    #    "nom",
    #)
    
    # Apply the fake factor only for the C category
    ff_nom = np.where((is_C_category | is_B_category), fake_factors_nom, 1.0)
    #ff_nom = np.where(is_C0_category, fake_0_factors_nom, ff_nom)

    #ext_corr_nom = np.where((is_C_category | is_C0_category), ext_corr_nom, 1.0)
    
    #closure_nom = np.where(is_C_category, closure_nom, 1.0)

    # Add the column to the events
    events = set_ak_column(events, "ff_weight", ff_nom, value_type=np.float32)
    # events = set_ak_column(events, "closure_weight", closure_nom, value_type=np.float32)
    #events = set_ak_column(events, "ff_ext_corr_weight", ext_corr_nom, value_type=np.float32)
    
    return events


@ff_weight.requires
def ff_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@ff_weight.setup
def ff_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    #from IPython import embed; embed()
    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.tautau_ff.path,
    ) 
    self.ff_corrector = correction_set["fake_factors_fit"]

    #correction_set_0 = correctionlib.CorrectionSet.from_file(
    #    bundle.files.tautau_ff0.path,
    #)
    #self.ff0_corrector = correction_set["fake_factors_fit"]

    ## extrapolation correction on FF
    #ext_correction_set = correctionlib.CorrectionSet.from_file(
    #    bundle.files.tautau_ext_corr.path,
    #)
    #self.ext_corrector = ext_correction_set["extrapolation_correction"]
    
