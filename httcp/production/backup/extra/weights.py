import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from columnflow.production.cms.muon import muon_weights
from columnflow.production.cms.electron import electron_weights

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


# -----> But, xTrig SFs are from different file

@producer(
    uses={
        "hcand.pt",
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
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Tau.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Muon.pt, axis=1)
    abseta = flat_np_view(abs(events.Muon.eta), axis=1)
    
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)

    xtrig_mask = cross_mu_triggered & ~single_mu_triggered & (events.hcand.pt[:,0] < 25)

    args = lambda mask, syst : (abseta[mask], pt[mask], syst)

    shifts=["nom", "up", "down"]

    sf_values = {}

    for the_shift in shifts:
        sf_values[the_shift] = sf_nom.copy()
        sf_values[the_shift][xtrig_mask] = self.sf_corrector(*args(xtrig_mask, the_shift))

        wt_name = "muon_xtrig_weight" if the_shift == "nom" else f"muon_xtrig_weight_{the_shift}"
        events = set_ak_column(events, wt_name, reduce_mul(sf_values[the_shift]), value_type=np.float32)

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
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.muon_xtrig_sf.load().decode("utf-8"),
    )
    corrector_name, self.year = self.get_electron_config()
    self.sf_corrector = correction_set[corrector_name]




# ------------------------------------------------- #
# Calculate Tau weights
# ------------------------------------------------- #

@producer(
    uses={
        f"Tau.{var}" for var in [
            "pt",
            "eta",
            #"decayMode",
            "decayModeHPS",
            "genPartFlav"
        ]
    },
    produces={
        "tau_weight",
    } | {
        f"tau_weight_{direction}"
        for direction in ["up", "down"]
    },
    mc_only=True,
)
def tau_weight(self: Producer, events: ak.Array, do_syst: bool, **kwargs) -> ak.Array:
    """
    Producer for tau scale factors derived by the TAU POG. Requires an external file in the
    config under ``tau_correction``:

        cfg.x.external_files = DotDict.wrap({
            "tau_correction": "/afs/cern.ch/user/s/stzakhar/work/higgs_cp/corrections/tau/POG/TAU/2022_preEE/tau_DeepTau2018v2p5_2022_preEE.json.gz", 
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files. A correction set named ``"tau_trigger"`` is extracted from it.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """

    #Helper function to deal with the case when two taus exist at the events. In that case one should multiply sf values to get the sf per event
    shape_sf = lambda sf: ak.prod(ak.unflatten(sf, 
                                               ak.num(events.Tau.pt, axis=1)), 
                                  axis=1, 
                                  mask_identity=False)
    
    #Make masks for each channel
    mask = {}
    for the_ch in ["etau","mutau", "tautau"]: mask[the_ch] = events.channel_id == self.config_inst.get_channel(the_ch).id
    
    
    #Prepare flat arrays of the inputs to send into the 
    pt = flat_np_view(events.Tau.pt, axis=1)
    eta = flat_np_view(abs(events.Tau.eta), axis=1)
    #dm = flat_np_view(events.Tau.decayMode, axis=1)
    dm = flat_np_view(events.Tau.decayModeHPS, axis=1)
    genmatch = flat_np_view(events.Tau.genPartFlav, axis=1)

    #from IPython import embed; embed()

    channel_id_brdcast, _ = ak.broadcast_arrays(events.channel_id[:,None],
                                                events.Tau.genPartFlav)
    channel_id_flat = flat_np_view(channel_id_brdcast, axis=1)

    
    deep_tau_tagger = self.config_inst.x.deep_tau_tagger
    args_vs_e = lambda mask, discr, syst   : (eta[mask],
                                              dm[mask],
                                              genmatch[mask],
                                              discr, #self.config_inst.x.deep_tau.vs_e, 
                                              syst)
    tes_args_vs_e = lambda mask, discr, syst : (pt[mask],
                                                eta[mask],
                                                dm[mask],
                                                genmatch[mask],
                                                deeptauid,
                                                id_wp,
                                                discr, #self.config_inst.x.deep_tau.vs_e, 
                                                syst)
    
    args_vs_mu = lambda mask, discr, syst  : (eta[mask],
                                              #dm[mask],
                                              genmatch[mask],
                                              discr, #self.config_inst.x.deep_tau.vs_mu, 
                                              syst)
    tes_args_vs_e = lambda mask, discr, syst : (pt[mask],
                                                eta[mask],
                                                dm[mask],
                                                genmatch[mask],
                                                deeptauid,
                                                id_wp,
                                                discr, #self.config_inst.x.deep_tau.vs_e, 
                                                syst)

    args_vs_jet = lambda mask, discr1, discr2, syst : (pt[mask],
                                                       dm[mask],
                                                       genmatch[mask],
                                                       discr1, #self.config_inst.x.deep_tau.vs_jet,
                                                       discr2, #self.config_inst.x.deep_tau.vs_e, 
                                                       '' if syst =="nom" else syst,
                                                       "dm")
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    
    shifts = ["nom"]
    if  do_syst:    shifts=[*shifts,"up", "down"]
    
    _sf = np.ones_like(pt, dtype=np.float32)
    sf_values = {}
    tes_sf_values = {}
    
    no_dm_5_and_6 = flat_np_view(((events.Tau.decayModeHPS != 5) & (events.Tau.decayModeHPS != 6)), axis=1)
    
    for the_shift in shifts:
        #print(f"tau_weight_shift: {the_shift}")
        #from IPython import embed; embed()
        sf_values[the_shift] = _sf.copy()
        #Calculate scale factors for tau vs electron classifier 
        e_mask = ((genmatch == tau_part_flav["prompt_e"]) | (genmatch == tau_part_flav["tau->e"])) & no_dm_5_and_6
        #sf_values[the_shift][e_mask] = self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask,the_shift))
        chid_e_mask = channel_id_flat[e_mask]
        sf_values[the_shift][e_mask] = \
            ak.where(
                chid_e_mask == 1,
                self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["etau"], the_shift)),
                ak.where(
                    chid_e_mask == 2,
                    self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["mutau"], the_shift)),
                    self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["tautau"], the_shift))))
        tes_sf_values[the_shift][e_mask] = \
            ak.where(
                chid_e_mask == 1,
                self.tes_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["etau"], the_shift)),
                ak.where(
                    chid_e_mask == 2,
                    self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["mutau"], the_shift)),
                    self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["tautau"], the_shift))))

        #Calculate scale factors for tau vs muon classifier 
        mu_mask = ((genmatch == tau_part_flav["prompt_mu"]) | (genmatch == tau_part_flav["tau->mu"])) & no_dm_5_and_6
        #sf_values[the_shift][mu_mask] = self.id_vs_e_corrector.evaluate(*args_vs_mu(mu_mask,the_shift))
        chid_mu_mask = channel_id_flat[mu_mask]
        sf_values[the_shift][mu_mask] = \
            ak.where(
                chid_mu_mask == 1,
                self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_m"]["etau"], the_shift)),
                ak.where(
                    chid_mu_mask == 2,
                    self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_m"]["mutau"], the_shift)),
                    self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_m"]["tautau"], the_shift))))
        
        
        #Calculate scale factors for tau vs jet classifier 
        tau_mask = (genmatch == tau_part_flav["tau_had"]) & no_dm_5_and_6
        #sf_values[the_shift][tau_mask] = self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,the_shift))
        chid_tau_mask = channel_id_flat[tau_mask]
        sf_values[the_shift][tau_mask] = \
            ak.where(
                chid_tau_mask == 1,
                self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                               self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_j"]["etau"],
                                                               self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["etau"],
                                                               the_shift)),
                ak.where(
                    chid_tau_mask == 2,
                    self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                   self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_j"]["mutau"],
                                                                   self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["mutau"],
                                                                   the_shift)),
                    self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                   self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_j"]["tautau"],
                                                                   self.config_inst.x.deep_tau_info[deep_tau_tagger]["vs_e"]["tautau"],
                                                                   the_shift))))
        
        #Save weights and their systematic variations to the main tree. For tautau the result is a product of two vs_jet weights of the corresponding taus
        wt_name = "tau_weight" if the_shift == "nom" else f"tau_weight_{the_shift}"
        events = set_ak_column(events,
                               wt_name,
                               shape_sf(sf_values[the_shift]),
                               value_type=np.float32)
        
    return events

@tau_weight.requires
def tau_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@tau_weight.setup
def tau_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        #bundle.files.tau_correction.load(formatter="gzip").decode("utf-8"),
        bundle.files.tau_sf.load(formatter="gzip").decode("utf-8"),
    )
    tagger_name = self.config_inst.x.deep_tau_tagger
    self.id_vs_jet_corrector    = correction_set[f"{tagger_name}VSjet"]
    self.id_vs_e_corrector      = correction_set[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector     = correction_set[f"{tagger_name}VSmu"]
    self.tes_corrector     = correction_set[f"tau_energy_scale"]


# ------------------------------------------------- #
# Calculate Tau-Spinner weights
# ------------------------------------------------- #

@producer(
    uses={
        f"TauSpinner*" 
    },
    produces={
        "tauspinner_weight_up", "tauspinner_weight", "tauspinner_weight_down"
    },
    mc_only=True,
)
def tauspinner_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    # https://github.com/hephysicist/CPinHToTauTau/blob/desy_dev/httcp/production/weights.py#L411
    """
    names = ["_up", "", "_down"]
    for the_name in names:
        #if  "TauSpinner" in list(events.fields):
        if the_name == "_up": the_weight = events.TauSpinner.weight_cp_0
        #elif the_name == "":  the_weight = ak.ones_like(events.TauSpinner.weight_cp_0p5)
        elif the_name == "":  the_weight =  (events.TauSpinner.weight_cp_0p5 + events.TauSpinner.weight_cp_0)/2.
        elif the_name == "_down": the_weight = events.TauSpinner.weight_cp_0p5
        else:  raise NotImplementedError('CP hypothesis is not known to the tauspinner weight producer!')   
        buf = ak.to_numpy(the_weight)
        if any(np.isnan(buf)):
            warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
            buf[np.isnan(buf)] = 0
            the_weight = buf

        events = set_ak_column_f32(events, f"tauspinner_weight{the_name}", the_weight)
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
