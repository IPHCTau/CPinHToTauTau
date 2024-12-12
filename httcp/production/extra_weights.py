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
        "zpt_reweight",
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
        bundle.files.zpt_rewt_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["zptreweight"]





# ------------------------------------------------- #
# Calculate FF weights (Dummy)
# ------------------------------------------------- #

@producer(
    uses={
        "channel_id",
        "category_ids",
        "hcand.pt", "hcand.mass",
    },
    produces={
        "ff_weight",
    },
    mc_only=True,
)
def ff_weight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    hcand1 = events.hcand[:,0]
    target_id = self.config_inst.x.ff_apply_id_map.tautau["id"]
    print(f"target_id : {target_id}")
    is_AR_id = ak.fill_none(ak.any(events.category_ids == target_id, axis=1), False)
    is_AR_id = flat_np_view(is_AR_id[:,None])

    is_outside_range = (
        ((hcand1.pt == 0.0) & (hcand1.mass == 0.0))
        | ((hcand1.pt >= 600.0) | (hcand1.mass >= 1000.0))
    )
    #from IPython import embed; embed()
    is_outside_range = flat_np_view(is_outside_range[:,None])
    # for safety
    zm  = ak.where(hcand1.mass > 1000.0, 999.99, hcand1.mass)
    zpt = ak.where(hcand1.pt > 600.0, 599.99, hcand1.pt)

    zm = flat_np_view(zm[:,None])
    zpt = flat_np_view(zpt[:,None])

    sf_nom_temp = 0.8*self.ff_corrector.evaluate(np.abs(zm), zpt)

    #from IPython import embed; embed()
    
    #sf_nom = np.where(is_outside_range, 1.0, np.where(is_AR_id, self.ff_corrector.evaluate(zm,zpt), 1.0))
    sf_nom = np.where(is_outside_range, 1.0, np.where(is_AR_id, sf_nom_temp, 1.0))

    events = set_ak_column(events, "ff_weight", sf_nom, value_type=np.float32)
    
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
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.tautau_ff.load(formatter="gzip").decode("utf-8"),
    )
    self.ff_corrector    = correction_set["zptreweight"]
