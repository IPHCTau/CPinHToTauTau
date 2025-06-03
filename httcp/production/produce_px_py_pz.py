import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

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
        "hcand.*", optional("GenTau.*"),
    },
    produces={
        # new columns
        "GenTau.px", "GenTau.py", "GenTau.pz", 
        "hcand.px", "hcand.py", "hcand.pz", 
        
        # Add new resolution columns here
        "hcand.resolution_rel__gen_to_reco__px","hcand.resolution_rel__gen_to_fastMTT__px","hcand.resolution_rel__reco_to_fastMTT__px",
        "hcand.resolution_rel__gen_to_reco__py","hcand.resolution_rel__gen_to_fastMTT__py","hcand.resolution_rel__reco_to_fastMTT__py",
        "hcand.resolution_rel__gen_to_reco__pz","hcand.resolution_rel__gen_to_fastMTT__pz","hcand.resolution_rel__reco_to_fastMTT__pz",
        "hcand.resolution_rel__gen_to_reco__pt","hcand.resolution_rel__gen_to_fastMTT__pt","hcand.resolution_rel__reco_to_fastMTT__pt",
        "hcand.resolution_abs__gen_to_reco__px","hcand.resolution_abs__gen_to_fastMTT__px","hcand.resolution_abs__reco_to_fastMTT__px",
        "hcand.resolution_abs__gen_to_reco__py","hcand.resolution_abs__gen_to_fastMTT__py","hcand.resolution_abs__reco_to_fastMTT__py",
        "hcand.resolution_abs__gen_to_reco__pz","hcand.resolution_abs__gen_to_fastMTT__pz","hcand.resolution_abs__reco_to_fastMTT__pz",
        "hcand.resolution_abs__gen_to_reco__pt","hcand.resolution_abs__gen_to_fastMTT__pt","hcand.resolution_abs__reco_to_fastMTT__pt",
        "hcand.resolution_abs__gen_to_reco__eta","hcand.resolution_abs__gen_to_fastMTT__eta","hcand.resolution_abs__reco_to_fastMTT__eta",
        "hcand.resolution_abs__gen_to_reco__phi","hcand.resolution_abs__gen_to_fastMTT__phi","hcand.resolution_abs__reco_to_fastMTT__phi",
    },

)

def produce_px_py_pz(
        self: Producer,
        events: ak.Array,
        **kwargs,
):
    # Ensure we have the correct behavior
    if not hasattr(events, "__behavior__") or events.__behavior__ is not coffea.nanoevents.methods.nanoaod.behavior:
        events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)

    # Lorentz vectors
    gentau_p4 = ak.with_name(events.GenTau, "PtEtaPhiMLorentzVector")
    hcand_p4 = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")

    # Px, Py, Pz calculations
    gentau_px = gentau_p4.px
    gentau_py = gentau_p4.py
    gentau_pz = gentau_p4.pz

    hcand_px = hcand_p4.px
    hcand_py = hcand_p4.py
    hcand_pz = hcand_p4.pz

    # Set the new columns for Px, Py, Pz
    events = set_ak_column(events, "GenTau.px", gentau_px)
    events = set_ak_column(events, "GenTau.py", gentau_py)
    events = set_ak_column(events, "GenTau.pz", gentau_pz)
    events = set_ak_column(events, "hcand.px", hcand_px)
    events = set_ak_column(events, "hcand.py", hcand_py)
    events = set_ak_column(events, "hcand.pz", hcand_pz)

    # Calculate and add the resolution columns
    variables = ["px", "py", "pz", "pt", "eta", "phi"]
    comparisons = {
        "gen_to_reco": ("GenTau", "hcand"),
        "gen_to_fastMTT": ("GenTau", "hcand_fastMTT"),
        "reco_to_fastMTT": ("hcand", "hcand_fastMTT")
    }

    for var in variables:
        for name, (denom, num) in comparisons.items():

            if "fastMTT" in num:
                num_var = f"{var}_fastMTT"
                if num_var in events["hcand"].fields:
                    num_vals = events["hcand"][num_var]
                else:
                    num_vals = None
            elif num == "hcand":
                num_vals = events["hcand"][var]
            else:
                num_vals = events[num][var]

            if denom == "GenTau":
                denom_vals = events["GenTau"][var]
            elif denom == "hcand":
                denom_vals = events["hcand"][var]
            else:
                denom_vals = events[denom][var]

            ## absolute resolutions for momenta
            abs_reso = denom_vals - num_vals
            events = set_ak_column(events, f"hcand.resolution_abs__{name}__{var}", abs_reso)

            ## relatove resolutions for momenta
            if var in {"px", "py", "pz", "pt"}:
                rel_reso = (denom_vals - num_vals) / denom_vals
                events = set_ak_column(events, f"hcand.resolution_rel__{name}__{var}", rel_reso)

    return events

@producer(
    uses={"GenTau.*",},
    produces={"GenHiggs_mass",}
)

def calculate_higgs_mass_genlevel(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    # Ensure we have the correct behavior
    if not hasattr(events, "__behavior__") or events.__behavior__ is not coffea.nanoevents.methods.nanoaod.behavior:
        events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)

    # Lorentz vectors
    gentau_p4 = ak.with_name(events.GenTau, "PtEtaPhiMLorentzVector")
    hcand_p4 = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")

    
    # Umrechnung von pt, eta und Masse in Lorentz-Vektoren
    gentau1_pt = events.GenTau.pt[:, 0]
    gentau1_eta = events.GenTau.eta[:, 0]
    gentau1_phi = events.GenTau.phi[:, 0]
    gentau1_mass = events.GenTau.mass[:, 0]
    
    gentau2_pt = events.GenTau.pt[:, 1]
    gentau2_eta = events.GenTau.eta[:, 1]
    gentau2_phi = events.GenTau.phi[:, 1]
    gentau2_mass = events.GenTau.mass[:, 1]
    
    # Berechnung der Energie (E) jedes Taus unter der Annahme, dass wir pt, eta und mass kennen
    gentau1_E = np.sqrt(gentau1_pt**2 + gentau1_mass**2) * np.cosh(gentau1_eta)
    gentau2_E = np.sqrt(gentau2_pt**2 + gentau2_mass**2) * np.cosh(gentau2_eta)
    
    # Berechnung der px, py und pz Komponenten der Tau-Leptonen
    gentau1_px = gentau1_pt * np.cos(gentau1_phi)
    gentau1_py = gentau1_pt * np.sin(gentau1_phi)
    gentau1_pz = gentau1_pt * np.sinh(gentau1_eta)
    
    gentau2_px = gentau2_pt * np.cos(gentau2_phi)
    gentau2_py = gentau2_pt * np.sin(gentau2_phi)
    gentau2_pz = gentau2_pt * np.sinh(gentau2_eta)
    
    # Berechnung der Gesamtmasse (Higgs-Masse) als Vektoraddition der beiden Taus
    total_px = gentau1_px + gentau2_px
    total_py = gentau1_py + gentau2_py
    total_pz = gentau1_pz + gentau2_pz
    total_E = gentau1_E + gentau2_E
    
    # Berechnung der Masse des Higgs Bosons
    higgs_mass = np.sqrt(total_E**2 - (total_px**2 + total_py**2 + total_pz**2))
    
    events["GenHiggs_mass"] = higgs_mass
    
    return events

    