import functools

from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

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


@producer(
    uses={
        # custom columns created upstream, probably by a selector
        "single_triggered", "cross_triggered",
        "single_e_triggered", "cross_e_triggered",
        "single_mu_triggered", "cross_mu_triggered",
        "cross_tau_triggered", "cross_tau_jet_triggered",
        # nano columns
        "Tau.pt", "Tau.eta", "Tau.genPartFlav", "Tau.decayModeHPS",
    },
    produces={
        "tau_weight",
    } | {
        f"tau_weight_{unc}_{direction}"
        for direction in ["up", "down"]
        for unc in [
                "jet_dm0", "jet_dm1", "jet_dm10", "e_barrel", "e_endcap",
                "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
        ]
    } | {
        "tau_trigger_weight"
    } | {
        f"tau_trigger_weight_{direction}"
        for direction in ["up", "down"]
    },
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_tau_file=(lambda self, external_files: external_files.tau_sf),
    # function to determine the tau tagger name
    get_tau_tagger=(lambda self: self.config_inst.x.deep_tau_tagger),
)
def tau_weights(self: Producer, events: ak.Array, do_syst: bool, **kwargs) -> ak.Array:
    """
    Producer for tau ID weights. Requires an external file in the config under ``tau_sf``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "tau_sf": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/TAU/2017_UL/tau.json.gz",  # noqa
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files.

    The name of the tagger should be given as an auxiliary entry in the config:

    .. code-block:: python

        cfg.x.tau_tagger = "DeepTau2017v2p1"

    It is used to extract correction set names such as "DeepTau2017v2p1VSjet". *get_tau_tagger* can
    be adapted in a subclass in case it is stored differently in the config.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """
    wp_config = self.config_inst.x.deep_tau_info[self.config_inst.x.deep_tau_tagger]

    # get channels from the config
    ch_etau = self.config_inst.get_channel("etau")
    ch_mutau = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Tau.pt), axis=1, mask_identity=False)
    #shape_sf = lambda sf: ak.prod(ak.unflatten(sf, ak.num(events.Tau.pt, axis=1)), axis=1, mask_identity=False)
    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayModeHPS, axis=1)
    genmatch = flat_np_view(events.Tau.genPartFlav, axis=1)

    # define channel / trigger dependent masks
    single_triggered = events.single_triggered
    cross_triggered = events.cross_triggered

    channel_id_brdcast, _  = ak.broadcast_arrays(events.channel_id[:,None], events.Tau.genPartFlav)
    channel_id_flat        = flat_np_view(channel_id_brdcast, axis=1)

    single_e_triggered, _  = ak.broadcast_arrays(events.single_e_triggered[:,None], events.Tau.genPartFlav)
    single_e_triggered     = flat_np_view(single_e_triggered, axis=1)

    cross_e_triggered, _   = ak.broadcast_arrays(events.cross_e_triggered[:,None], events.Tau.genPartFlav)
    cross_e_triggered      = flat_np_view(cross_e_triggered, axis=1)

    single_mu_triggered, _ = ak.broadcast_arrays(events.single_mu_triggered[:,None], events.Tau.genPartFlav)
    single_mu_triggered    = flat_np_view(single_mu_triggered, axis=1)
    
    cross_mu_triggered, _  = ak.broadcast_arrays(events.cross_mu_triggered[:,None], events.Tau.genPartFlav)
    cross_mu_triggered     = flat_np_view(cross_mu_triggered, axis=1)
    
    cross_tau_triggered, _ = ak.broadcast_arrays(events.cross_tau_triggered[:,None], events.Tau.genPartFlav)
    cross_tau_triggered    = flat_np_view(cross_tau_triggered, axis=1)

    cross_tau_jet_triggered, _ = ak.broadcast_arrays(events.cross_tau_jet_triggered[:,None], events.Tau.genPartFlav)
    cross_tau_jet_triggered    = flat_np_view(cross_tau_jet_triggered, axis=1)
    
    
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }

    dm_mask = (
        (events.Tau.decayModeHPS == 0) |
        (events.Tau.decayModeHPS == 1) |
        (events.Tau.decayModeHPS == 10)|
        (events.Tau.decayModeHPS == 11)
    )
  
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)


    # helpers to create corrector arguments
    if self.id_vs_e_corrector.version == 0:
        args_vs_e   = lambda mask, id_vs_e_wp, syst  : (abseta[mask], genmatch[mask], id_vs_e_wp, syst)
    elif self.id_vs_e_corrector.version in (1,):
        args_vs_e   = lambda mask, id_vs_e_wp, syst  : (abseta[mask], dm[mask], genmatch[mask], id_vs_e_wp, syst)
    else:
        raise NotImplementedError
        
    args_vs_mu  = lambda mask, id_vs_m_wp, syst  : (abseta[mask], genmatch[mask], id_vs_m_wp, syst)    
    args_vs_jet = lambda mask, id_vs_j_wp, id_vs_e_wp, syst : (pt[mask], dm[mask], genmatch[mask], id_vs_j_wp, id_vs_e_wp, syst, "dm")
    
    trig_eval_args     = lambda mask, type, discr, syst: (pt[mask], dm[mask], type, discr, "sf", syst)
    trig_eff_eval_args = lambda mask, type, discr, node, syst: (pt[mask], dm[mask], type, discr, node, syst)



    shifts = ["nom"]
    if  do_syst:
        shifts=[*shifts,"up", "down"]

    sf_values = {}
    sf_values_dm = {}
    sf_values_e = {}
    sf_values_mu = {}
    trig_sf_values = {}


    
    for the_shift in shifts:
        sf_values[the_shift] = sf_nom.copy() # sf_values = {"nom": [ , , , , ...]}
        trig_sf_values[the_shift] = sf_nom.copy() # sf_values = {"nom": [ , , , , ...]}
        
        ##############################################################
        # ------------------- for ele -> tau fake ------------------ #
        ##############################################################
        e_mask = ((events.Tau.genPartFlav == tau_part_flav["prompt_e"]) | (events.Tau.genPartFlav == tau_part_flav["tau->e"]))
        if self.config_inst.campaign.x.run == 3:
            e_mask = e_mask & (events.Tau.decayModeHPS != 5) & (events.Tau.decayModeHPS != 6)
        e_mask = flat_np_view(e_mask, axis=1)
        ch_id_e_mask = channel_id_flat[e_mask]
        sf_values[the_shift][e_mask] = ak.where(ch_id_e_mask == ch_etau.id,
                                                self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, wp_config["vs_e"]["etau"], the_shift)),
                                                ak.where(ch_id_e_mask == ch_mutau.id,
                                                         self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, wp_config["vs_e"]["mutau"], the_shift)),
                                                         self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, wp_config["vs_e"]["tautau"], the_shift))))

        if do_syst:
            # --->>> electron fakes -> split into 2 eta regions [up/down only]
            for region, region_mask in [
                    ("barrel", (abseta < 1.5)),
                    ("endcap", (abseta >= 1.5)),
            ]:
                if the_shift == "nom": continue
                sf_values_e[the_shift] = sf_nom.copy()
                e_single_region_mask = e_mask & single_e_triggered & region_mask
                e_cross_region_mask = e_mask & cross_e_triggered & ~single_e_triggered & region_mask
                
                sf_values_e[the_shift][e_single_region_mask] = self.id_vs_e_corrector.evaluate(*args_vs_e(e_single_region_mask, wp_config["vs_e"]["etau"], the_shift))
                sf_values_e[the_shift][e_cross_region_mask] = self.id_vs_e_corrector.evaluate(*args_vs_e(e_cross_region_mask, wp_config["vs_e"]["etau"], the_shift))
                
                wt_name = f"tau_weight_e_{region}" if the_shift == "nom" else f"tau_weight_e_{region}_{the_shift}"
                events = set_ak_column(events, wt_name, reduce_mul(sf_values_e[the_shift]), value_type=np.float32)

            
        # trigger sf
        trig_e_mask = (channel_id_flat == ch_etau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 25.0)), axis=1) & cross_e_triggered & ~single_e_triggered
        trig_sf_values[the_shift][trig_e_mask] = self.trig_corrector.evaluate(*trig_eval_args(trig_e_mask, 'etau', wp_config["vs_j"]["etau"], the_shift))

        ##############################################################

        ##############################################################
        # ------------------ for muon -> tau fake ------------------ #
        ##############################################################
        
        mu_mask = ((events.Tau.genPartFlav == tau_part_flav["prompt_mu"]) | (events.Tau.genPartFlav == tau_part_flav["tau->mu"]))
        mu_mask = flat_np_view(mu_mask, axis=1)
        ch_id_mu_mask = channel_id_flat[mu_mask]
        sf_values[the_shift][mu_mask] = ak.where(ch_id_mu_mask == ch_etau.id,
                                                 self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, wp_config["vs_m"]["etau"], the_shift)),
                                                 ak.where(ch_id_mu_mask == ch_mutau.id,
                                                          self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, wp_config["vs_m"]["mutau"], the_shift)),
                                                          self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, wp_config["vs_m"]["tautau"], the_shift))))

        if do_syst:
            # --->>> muon fakes -> split into 5 eta regions [up/down]
            for region, region_mask in [
                    ("0p0To0p4", (abseta < 0.4)),
                    ("0p4To0p8", ((abseta >= 0.4) & (abseta < 0.8))),
                    ("0p8To1p2", ((abseta >= 0.8) & (abseta < 1.2))),
                    ("1p2To1p7", ((abseta >= 1.2) & (abseta < 1.7))),
                    ("1p7To2p3", (abseta >= 1.7)),
            ]:
                if the_shift == "nom": continue
                sf_values_mu[the_shift] = sf_nom.copy()
                mu_single_region_mask = mu_mask & single_mu_triggered & region_mask
                mu_cross_region_mask = mu_mask & cross_mu_triggered & ~single_mu_triggered & region_mask
                
                sf_values_mu[the_shift][mu_single_region_mask] = self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_single_region_mask, wp_config["vs_m"]["etau"], the_shift))
                sf_values_mu[the_shift][mu_cross_region_mask] = self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_cross_region_mask, wp_config["vs_m"]["etau"], the_shift))
                
                wt_name = f"tau_weight_mu_{region}" if the_shift == "nom" else f"tau_weight_mu_{region}_{the_shift}"
                events = set_ak_column(events, wt_name, reduce_mul(sf_values_e[the_shift]), value_type=np.float32)
                

        # trigger sf
        trig_mu_mask = (channel_id_flat == ch_mutau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 25.0)), axis=1) & cross_mu_triggered & ~single_mu_triggered
        trig_sf_values[the_shift][trig_mu_mask] = self.trig_corrector.evaluate(*trig_eval_args(trig_mu_mask, 'mutau', wp_config["vs_j"]["mutau"], the_shift))

        ##############################################################            
            
        ##############################################################
        # -------------------- for genuine taus -------------------- #
        ##############################################################
        tau_mask = dm_mask & (events.Tau.genPartFlav == tau_part_flav["tau_had"])
        tau_mask = flat_np_view(tau_mask, axis=1)
        ch_id_tau_mask = channel_id_flat[tau_mask]
        sf_values[the_shift][tau_mask] = ak.where(ch_id_tau_mask == ch_etau.id,
                                                  self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                                                 wp_config["vs_j"]["etau"],
                                                                                                 wp_config["vs_e"]["etau"],
                                                                                                 the_shift)),
                                                  ak.where(ch_id_tau_mask == ch_mutau.id,
                                                           self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                                                          wp_config["vs_j"]["mutau"],
                                                                                                          wp_config["vs_e"]["mutau"],
                                                                                                          the_shift)),
                                                           self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                                                          wp_config["vs_j"]["tautau"],
                                                                                                          wp_config["vs_e"]["tautau"],
                                                                                                          the_shift))))

        if do_syst:
            # --->>> genuine taus for DM 0, 1, 10 [only up/down variations]
            for idm in [0, 1, 10]:
                if the_shift == "nom": continue
                sf_values_dm[the_shift] = sf_nom.copy()
                tau_dmX_mask = tau_mask & (dm == idm)
                ch_id_tau_dmX_mask = channel_id_flat[tau_dmX_mask]
                sf_values_dm[the_shift][tau_dmX_mask] = ak.where(ch_id_tau_dmX_mask == ch_etau.id,
                                                                 self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_dmX_mask,
                                                                                                                wp_config["vs_j"]["etau"],
                                                                                                                wp_config["vs_e"]["etau"],
                                                                                                                the_shift)),
                                                                 ak.where(ch_id_tau_dmX_mask == ch_mutau.id,
                                                                          self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_dmX_mask,
                                                                                                                         wp_config["vs_j"]["mutau"],
                                                                                                                         wp_config["vs_e"]["mutau"],
                                                                                                                         the_shift)),
                                                                          self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_dmX_mask,
                                                                                                                         wp_config["vs_j"]["tautau"],
                                                                                                                         wp_config["vs_e"]["tautau"],
                                                                                                                         the_shift))))
                
                wt_name = f"tau_weight_jet_dm{idm}" if the_shift == "nom" else f"tau_weight_jet_dm{idm}_{the_shift}"
                events = set_ak_column(events, wt_name, reduce_mul(sf_values_dm[the_shift]), value_type=np.float32)

        # trig sf
        """
        trig_tau_mask = (channel_id_flat == ch_tautau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 40.0)), axis=1) & cross_tau_triggered & ~cross_tau_jet_triggered
        trig_sf_values[the_shift][trig_tau_mask] = self.trig_corrector.evaluate(*trig_eval_args(trig_tau_mask,
                                                                                                'ditau',
                                                                                                wp_config["vs_j"]["tautau"],
                                                                                                the_shift))
        trig_tau_jet_mask = (channel_id_flat == ch_tautau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 40.0)), axis=1) & cross_tau_jet_triggered & ~cross_tau_triggered
        trig_sf_values[the_shift][trig_tau_jet_mask] = self.trig_corrector.evaluate(*trig_eval_args(trig_tau_jet_mask,
                                                                                                    'ditaujet',
                                                                                                    wp_config["vs_j"]["tautau"],
                                                                                                    the_shift))
        """
        # Inclusive OR
        #trig_tau_mask = (channel_id_flat == ch_tautau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 40.0)), axis=1) & (cross_tau_triggered | cross_tau_jet_triggered)
        #from IPython import embed; embed()

        trig_tautau_mask = (channel_id_flat == ch_tautau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 40.0)), axis=1) & cross_tau_triggered
        trig_tautaujet_mask = (channel_id_flat == ch_tautau.id) & flat_np_view((dm_mask & (events.Tau.pt >= 40.0)), axis=1) & cross_tau_jet_triggered
        trig_mask = (trig_tautau_mask | trig_tautaujet_mask)
        

        #from IPython import embed; embed()

        
        tautau_trig_eff_data = self.trig_corrector.evaluate(*trig_eff_eval_args(trig_mask,
                                                                                'ditau',
                                                                                wp_config["vs_j"]["tautau"],
                                                                                "eff_data",
                                                                                the_shift))
        tautau_trig_eff_mc  = self.trig_corrector.evaluate(*trig_eff_eval_args(trig_mask,
                                                                               'ditau',
                                                                               wp_config["vs_j"]["tautau"],
                                                                               "eff_mc",
                                                                               the_shift))
        tautaujet_trig_eff_data = self.trig_corrector.evaluate(*trig_eff_eval_args(trig_mask,
                                                                                   'ditaujet',
                                                                                   wp_config["vs_j"]["tautau"],
                                                                                   "eff_data",
                                                                                   the_shift))
        tautaujet_trig_eff_mc  = self.trig_corrector.evaluate(*trig_eff_eval_args(trig_mask,
                                                                                  'ditaujet',
                                                                                  wp_config["vs_j"]["tautau"],
                                                                                  "eff_mc",
                                                                                  the_shift))
        trig_tautau_mask_int = trig_tautau_mask.astype(int)[trig_mask]
        trig_tautaujet_mask_int = trig_tautaujet_mask.astype(int)[trig_mask]

        eff_data = trig_tautau_mask_int*tautau_trig_eff_data + trig_tautaujet_mask_int*tautaujet_trig_eff_data \
            - (trig_tautau_mask_int * trig_tautaujet_mask_int * tautau_trig_eff_data * tautaujet_trig_eff_data)
        eff_mc = trig_tautau_mask_int*tautau_trig_eff_mc + trig_tautaujet_mask_int*tautaujet_trig_eff_mc \
            - (trig_tautau_mask_int * trig_tautaujet_mask_int * tautau_trig_eff_mc * tautaujet_trig_eff_mc)

        sf = eff_data/eff_mc

        trig_sf_values[the_shift][trig_mask] = sf

        
        trig_wt_name = "tau_trigger_weight" if the_shift == "nom" else f"tau_trigger_weight_{the_shift}"
        events = set_ak_column(events, trig_wt_name, reduce_mul(trig_sf_values[the_shift]), value_type=np.float32)
        ##############################################################
            
        wt_name = "tau_weight" if the_shift == "nom" else f"tau_weight_{the_shift}"
        events = set_ak_column(events, wt_name, reduce_mul(sf_values[the_shift]), value_type=np.float32)

        
    return events



@tau_weights.requires
def tau_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@tau_weights.setup
def tau_weights_setup(self: Producer, reqs: dict, inputs: dict, reader_targets: InsertableDict) -> None:
    bundle = reqs["external_files"]

    # create the trigger and id correctors
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_tau_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    tagger_name = self.get_tau_tagger()
    # id
    self.id_vs_jet_corrector = correction_set[f"{tagger_name}VSjet"]
    self.id_vs_e_corrector = correction_set[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector = correction_set[f"{tagger_name}VSmu"]
    # trigger
    self.trig_corrector = correction_set["tau_trigger"]

    # check versions
    assert self.id_vs_jet_corrector.version in (0, 1, 2, 3)
    assert self.id_vs_e_corrector.version in (0, 1)
    assert self.id_vs_mu_corrector.version in (0, 1)



# ------------------------------------------------- #
# Calculate Tau-Spinner weights
# ------------------------------------------------- #

#@producer(
#    uses={
#        f"TauSpinner*" 
#    },
#    produces={
#        "tauspinner_weight_up", "tauspinner_weight", "tauspinner_weight_down"
#    },
#    mc_only=True,
#)
#def tauspinner_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
#    """
#    A simple function that sets tauspinner_weight according to the cp_hypothesis
#    # https://github.com/hephysicist/CPinHToTauTau/blob/desy_dev/httcp/production/weigh#ts.py#L411
#    """
#    names = ["_up", "", "_down"]
#
#    for the_name in names:
#        #if  "TauSpinner" in list(events.fields):
#        if the_name == "_up": the_weight = events.TauSpinner.weight_cp_0
#        #elif the_name == "":  the_weight = ak.ones_like(events.TauSpinner.weight_cp_0p#5)
#        elif the_name == "":  the_weight =  (events.TauSpinner.weight_cp_0p5 + events.T#auSpinner.weight_cp_0)/2.
#        elif the_name == "_down": the_weight = events.TauSpinner.weight_cp_0p5
#        else:  raise NotImplementedError('CP hypothesis is not known to the tauspinner #weight producer!')   
#        buf = ak.to_numpy(the_weight)
#        if any(np.isnan(buf)):
#            warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
#            buf[np.isnan(buf)] = 0
#            the_weight = buf
#    
#        events = set_ak_column_f32(events, f"tauspinner_weight{the_name}", the_weight)
#    return events


@producer(
    uses={
        f"TauSpinner*" 
    },
    produces={
        # Version with _alt are duplicated weights computed with different setup of
        # neutral currents parameterization and can be used to estimate uncertainty
        # of the weighting. However it was neglected in the Run-2 analysis
        "tauspinner_weight",
        "tauspinner_weight_up",   # same as cpeven
        "tauspinner_weight_down", # same as cpodd 
        "tauspinner_weight_cpeven",    # 0
        "tauspinner_weight_cpeven_alt",# 0
        "tauspinner_weight_cpmix",     # 0.25
        "tauspinner_weight_cpmix_alt", # 0.25_alt
        "tauspinner_weight_cpmixm",    # -0.25
        "tauspinner_weight_cpmixm_alt",# -0.25_alt
        "tauspinner_weight_cpalpha0p375",     # 0.375
        "tauspinner_weight_cpalpha0p375_alt", # 0.375_alt
        "tauspinner_weight_cpodd",     # 0.5
        "tauspinner_weight_cpodd_alt", # 0.5
    },
    mc_only=True,
)
def tauspinner_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    # https://github.com/hephysicist/CPinHToTauTau/blob/desy_dev/httcp/production/weights.py#L411
    """
    wt_names     = ['nom',
                    'cpeven',
                    'cpeven_alt',
                    'cpmix',
                    'cpmix_alt',
                    'cpmixm',
                    'cpmixm_alt',
                    'cpalpha0p375',
                    'cpalpha0p375_alt',
                    'cpodd',
                    'cpodd_alt']

    branch_names = ['',
                    'weight_cp_0',
                    'weight_cp_0_alt',
                    'weight_cp_0p25',
                    'weight_cp_0p25_alt',
                    'weight_cp_minus0p25',
                    'weight_cp_minus0p25_alt',
                    'weight_cp_0p375',
                    'weight_cp_0p375_alt',
                    'weight_cp_0p5',
                    'weight_cp_0p5_alt']
    
    weight_map   = zip(wt_names, branch_names)
    
    for (wt_name, branch) in weight_map:
        _name = ""
        if wt_name == "nom":
            weight = (events.TauSpinner.weight_cp_0p5 + events.TauSpinner.weight_cp_0)/2.
        else:
            _name = f"_{wt_name}"
            weight = events.TauSpinner[branch]

        buf = ak.to_numpy(weight)
        if any(np.isnan(buf)):
            warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
            buf[np.isnan(buf)] = 0
            weight = buf
        
        events = set_ak_column_f32(events, f"tauspinner_weight{_name}", weight)
        if _name == "_cpeven": # redundant, needs to be resolved later
            events = set_ak_column_f32(events, f"tauspinner_weight_up", weight)
        elif _name == "_cpodd": # redundant, needs to be resolved later
            events = set_ak_column_f32(events, f"tauspinner_weight_down", weight)

    return events

