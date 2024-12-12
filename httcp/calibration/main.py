# coding: utf-8

"""
main calibration script
"""

import functools

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.jets import jets, jec, jec_nominal, jer
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

from httcp.calibration.tau import tau_energy_scale
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

# derive calibrators to add settings
#jec_full = jec.derive("jec_full", cls_dict={"mc_only": True, "nominal_only": True})

@calibrator(
    uses={
        deterministic_seeds,
        "Jet.*",
        "RawPuppiMET.*",
        "PuppiMET.*",
        jets,
        tau_energy_scale,
    },
    produces={
        deterministic_seeds,
        "Jet.pt_no_corr", "Jet.phi_no_corr", "Jet.eta_no_corr", "Jet.mass_no_corr",
        "PuppiMET.pt_no_corr", "PuppiMET.phi_no_corr",
        jets,
        tau_energy_scale,
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    events = self[deterministic_seeds](events, **kwargs)

    events = set_ak_column_f32(events, "Jet.pt_no_corr", events.Jet.pt)
    events = set_ak_column_f32(events, "Jet.phi_no_corr", events.Jet.phi)
    events = set_ak_column_f32(events, "Jet.eta_no_corr", events.Jet.eta)
    events = set_ak_column_f32(events, "Jet.mass_no_corr", events.Jet.mass)
    
    #PuppiMET variables before applying energy corrections
    events = set_ak_column_f32(events, "PuppiMET.pt_no_corr", events.PuppiMET.pt)
    events = set_ak_column_f32(events, "PuppiMET.phi_no_corr", events.PuppiMET.phi)

    if self.config_inst.campaign.x.run == 3:
        events = ak.with_field(events, events.RawPuppiMET, "RawMET")
        events = ak.with_field(events, events.PuppiMET, "MET")
        events = ak.without_field(events, "RawPuppiMET")
        events = ak.without_field(events, "PuppiMET")

    events = self[jets](events, **kwargs)

    if self.config_inst.campaign.x.run == 3:
        events = ak.with_field(events, events.RawMET, "RawPuppiMET")
        events = ak.with_field(events, events.MET, "PuppiMET")
        events = ak.without_field(events, "RawMET")
        events = ak.without_field(events, "MET")

    if self.dataset_inst.is_mc: 
        ##Apply tau energy scale correction
        events = self[tau_energy_scale](events, **kwargs)


    return events
