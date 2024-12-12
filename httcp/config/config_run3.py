# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools
import itertools

import os
import re
import law
import yaml
from glob import glob
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.columnar_util import ColumnCollection, EMPTY_FLOAT
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category,
    add_shift_aliases, 
    get_shifts_from_sources,
    verify_config_processes,
)
#from httcp.config.config_util import add_expanded_shift_aliases

logger = law.logger.get_logger(__name__)


ak = maybe_import("awkward")

#thisdir = os.path.dirname(os.path.abspath(__file__))
thisdir = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/httcp/config"
#print(f"thisdir: {thisdir}")
corrdir = os.path.join(os.path.dirname(thisdir), "data")
#print(f"corrdir: {corrdir}")

def add_config (ana: od.Analysis,
                campaign: od.Campaign,
                config_name           = None,
                config_id             = None,
                limit_dataset_files   = None
) -> od.Config :

    logger.info(f"campaign    : {campaign.name}")
    logger.info(f"config name : {config_name}")
    logger.info(f"config id   : {config_id}")

    # gather campaign data
    year = campaign.x.year
    postfix = campaign.x.postfix
    
    # some validations
    assert year in [2022, 2023, 2024]

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)

    
    # --------------------------------------------------------------------------------------------- #
    # create a config by passing the campaign, so id and name will be identical    
    # --------------------------------------------------------------------------------------------- #

    cfg = ana.add_config(campaign,
                         name  = config_name,
                         id    = config_id)

    
    # --------------------------------------------------------------------------------------------- #
    # add processes we are interested in
    # --------------------------------------------------------------------------------------------- #
    process_names = [
        ## Data
        "data",
        ## W + jets
        "w_lnu",
        ## Drell-Yan
        "dy",
        #"dy_m10to50",
        #"dy_m50toinf",
        "dy_m50toinf_lep",
        "dy_m50toinf_tau",
        "dy_m50toinf_jet",
        ## TTJets
        "tt",
        ## Single top
        "st",
        ## VV [diboson inclusive]
        "vv",
        "vvv",
        ## Signal
        "h_ggf_htt",
        "zh_htt",
        "wh_htt",
    ]

    for process_name in process_names:
        # development switch in case datasets are not _yet_ there
        if process_name not in procs:
            logger.warning(f"WARNING: {process_name} not in cmsdb processes")
            continue
        # add the process
        #if process_name == "h_ggf_tautau":
        #    procs.get(process_name).is_signal = True
        proc = cfg.add_process(procs.get(process_name))
        #print(procs.get(process_name))
        #if proc.name == "h_ggf_tautau":
        #    proc.is_signal = True

    # configuration of colors, labels, etc. can happen here
    from httcp.config.styles import stylize_processes
    stylize_processes(cfg)    

    
    # --------------------------------------------------------------------------------------------- #
    # add datasets we need to study
    # --------------------------------------------------------------------------------------------- #

    dataset_names = [
        ##W+jets
        # --- LO --- #
        "wj_incl_madgraph",
        "wj_1j_madgraph",
        "wj_2j_madgraph",
        "wj_3j_madgraph",
        "wj_4j_madgraph",
        #"wj_ht40to100_madgraph",
        #"wj_ht100to400_madgraph",
        #"wj_ht400to800_madgraph",
        #"wj_ht800to1500_madgraph",
        #"wj_ht1500to2500_madgraph",
        # --- NLO --- #
        #"wj_incl_amcatnlo",
        #"wj_0j_amcatnlo",
        #"wj_1j_amcatnlo",
        #"wj_2j_amcatnlo",

        ##Drell-Yan
        # --- LO --- #
        "dy_lep_m10to50_madgraph",
        "dy_lep_m50_madgraph",
        "dy_lep_m50_1j_madgraph",
        "dy_lep_m50_2j_madgraph",
        "dy_lep_m50_3j_madgraph",
        "dy_lep_m50_4j_madgraph",
        # --- NLO --- #
        #"dy_lep_m10to50_amcatnlo",
        #"dy_lep_m50_amcatnlo",
        #"dy_lep_m50_0j_amcatnlo",
        #"dy_lep_m50_1j_amcatnlo",
        #"dy_lep_m50_2j_amcatnlo",
        #"dy_lep_m50_1j_pt40to100_amcatnlo",
        #"dy_lep_m50_2j_pt40to100_amcatnlo",
        #"dy_lep_m50_1j_pt100to200_amcatnlo",
        #"dy_lep_m50_2j_pt100to200_amcatnlo",
        #"dy_lep_m50_1j_pt200to400_amcatnlo",
        #"dy_lep_m50_2j_pt200to400_amcatnlo",
        #"dy_lep_m50_1j_pt400to600_amcatnlo",
        #"dy_lep_m50_2j_pt400to600_amcatnlo",
        #"dy_lep_m50_1j_pt600toInf_amcatnlo",
        #"dy_lep_m50_2j_pt600toInf_amcatnlo",

        ## ttbar
        "tt_sl",
        "tt_dl",
        "tt_fh",

        ##single top
        "st_tchannel_t",
        "st_tchannel_tbar",
        "st_tw_t_sl",
        "st_tw_tb_sl",
        "st_tw_t_dl",
        "st_tw_tb_dl",
        "st_tw_t_fh",
        "st_tw_tb_fh",
        "st_schannel_t",
        "st_schannel_tbar",

        ##Diboson
        "ww",
        "wz",
        "zz",

        ##Triboson
        "www",
        "wwz",
        "wzz",
        "zzz",

        ##Signal
        "h_ggf_tautau_uncorrelated_filter",
        "h_ggf_tautau_uncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
        "h_ggf_tautau_uncorrelatedDecay_CPodd_UnFiltered_ProdAndDecay",
        "h_ggf_tautau_uncorrelatedDecay_MM_Filtered_ProdAndDecay",
        "h_ggf_tautau_uncorrelatedDecay_MM_UnFiltered_ProdAndDecay",
        "h_ggf_tautau_uncorrelatedDecay_SM_Filtered_ProdAndDecay",
        "h_ggf_tautau_uncorrelatedDecay_SM_UnFiltered_ProdAndDecay",
        "h_ggf_tautau_M125_amcatnloFXFX",
        "h_ggf_tautau_prod_cp_even_sm",
        "zh_tautau_uncorrelatedDecay_Filtered",
        "zh_tautau_uncorrelatedDecay_UnFiltered",
        "wph_tautau_uncorrelatedDecay_Filtered",
        "wph_tautau_uncorrelatedDecay_UnFiltered",
        "wmh_tautau_uncorrelatedDecay_Filtered",
        "wmh_tautau_uncorrelatedDecay_UnFiltered",
    ]
    
    datasets_data = []
    if year == 2022:
        if postfix == "preEE":
            datasets_data = ["data_e_C",   "data_e_D",
                             "data_single_mu_C",
                             "data_mu_C",  "data_mu_D", 
                             "data_tau_C", "data_tau_D"]

        elif postfix == "postEE":
            datasets_data = ["data_e_E",   "data_e_G",
                             "data_mu_E",  "data_mu_G",
                             "data_tau_E", "data_tau_F", "data_tau_G"]
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2023:
        if postfix == "preBPix":
            datasets_data = ["data_e0_C",  "data_e1_C",
                             "data_mu0_C", "data_mu1_C",
                             "data_tau_C"]
        elif postfix == "postBPix":
            datasets_data = ["data_e0_D",  "data_e1_D",
                             "data_mu0_D", "data_mu1_D",
                             "data_tau_D"]
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2024:
        raise RuntimeError("2024? Too Early ...")

    else:
        raise RuntimeError(f"Check Year in __init__.py in cmsdb campaign: {year}")


    dataset_names = datasets_data + dataset_names
    for dataset_name in dataset_names:
        # development switch in case datasets are not _yet_ there
        if dataset_name not in campaign.datasets:
            logger.warning(f"WARNING: {dataset_name} not in cmsdb campaign")
            continue
        
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if re.match(r"^(ww|wz|zz|www|wwz|wzz|zzz)$", dataset.name):
            dataset.add_tag("no_lhe_weights")
        elif re.match(r"^dy_lep_m50.*$", dataset.name):
            dataset.add_tag("is_dy")
        elif re.match(r"^wj.*$", dataset.name):
            dataset.add_tag("is_w")
        elif re.match(r"^h_ggf_tautau.*$", dataset.name):
            dataset.add_tag("is_ggf_signal")
        elif re.match(r"^(.*)h_tautau(.*)Filtered$", dataset.name):
            dataset.add_tag("is_vh_signal")
        
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    
    # --------------------------------------------------------------------------------------------- #
    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    # --------------------------------------------------------------------------------------------- #

    cfg.x.default_calibrator      = "main"
    cfg.x.default_selector        = "main"
    cfg.x.default_producer        = "main"
    cfg.x.default_ml_model        = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories      = ("incl",)
    #cfg.x.default_variables = ("n_jet", "jet1_pt")
    cfg.x.default_variables       = ("event","channel_id")
    cfg.x.default_weight_producer = "main"
    
    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "backgrounds": (backgrounds := [
            "w_lnu", #"w_lnu",
            "tt",
            "dy_m50toinf", "dy_m10to50",
            "st",
            "vv",
        ]),
        "bkg_sig"       : (bkg_sig       := [*backgrounds, "h_ggf_htt"]),
        "data_bkg"      : (data_bkg      := [*backgrounds, "data"]),
        "data_bkg_sig"  : (data_bkg_sig  := [*backgrounds, "data", "h_ggf_htt"]),
    }
    cfg.x.process_settings_groups = {
        "unstack_processes": {proc: {"unstack": True, "scale": 10.0} for proc in ("h_ggf_htt")},
    }
    cfg.x.general_settings_groups = {
        "compare_shapes": {"skip_ratio": False, "shape_norm": False, "yscale": "log"} #"cms_label": "simpw"},
    }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json",
                    "met_filter",
                    "trigger",
                    "b_veto",
                    "has_2_or_more_leps_with_at_least_1_tau",
                    "dilepton_veto",
                    "has_at_least_1_pair_before_trigobj_matching",
                    "has_at_least_1_pair_after_trigobj_matching",
                    "extra_lepton_veto",
                    "One_higgs_cand_per_event",
                    "has_proper_tau_decay_products"],
    }

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    # define inclusive datasets for the stitched process identification with corresponding leaf processes
    # drell-yan [NLO]
    cfg.x.allow_dy_stitching = False
    cfg.x.dy_stitching = {
        "dy": {
            "inclusive_dataset": cfg.datasets.n.dy_lep_m50_madgraph,
            "leaf_processes": [
                # the following processes cover the full njet phasespace
                *(
                    procs.get(f"dy_m50toinf_{nj}j")
                    for nj in [0, 1, 2, 3, 4]
                ),
            ],
        },
    }
    """
    cfg.x.dy_stitching = {
        "dy": {
            "inclusive_dataset": cfg.datasets.n.dy_lep_m50_amcatnlo,
            "leaf_processes": [
                # the following processes cover the full njet and pt phasespace
                procs.n.dy_m50toinf_0j,
                *(
                    procs.get(f"dy_m50toinf_{nj}j_pt{pt}")
                    for nj in [1, 2]
                    for pt in ["0to40", "40to100", "100to200", "200to400", "400to600", "600toinf"]
                ),
                #procs.n.dy_m50toinf_ge3j,
            ],
        },
    }
    """
    # w+jets [NLO]
    cfg.x.allow_w_stitching = True
    cfg.x.w_stitching = {
        "wj": {
            "inclusive_dataset": cfg.datasets.n.wj_incl_madgraph,
            "leaf_processes": [
                # the following processes cover the full njet phasespace
                *(
                    procs.get(f"w_lnu_{nj}j") # njet from LO samples
                    for nj in [0,1,2,3,4]
                ),
            ],
        },
    }
    """
    cfg.x.w_stitching = {
        "wj": {
            "inclusive_dataset": cfg.datasets.n.wj_incl_amcatnlo,
            "leaf_processes": [
                # the following processes cover the full njet and pt phasespace
                *(
                    procs.get(f"w_lnu_{nj}j") # njet from NLO samples
                    for nj in [0,1,2]
                ),
            ],
        },
    }
    """
    
    
    # --------------------------------------------------------------------------------------------- #
    # Luminosity and Normalization
    # lumi values in inverse fb
    # TODO later: preliminary luminosity using norm tag. Must be corrected, when more data is available
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # --------------------------------------------------------------------------------------------- #

    if year == 2022:
        if postfix == "preEE":
            cfg.x.luminosity = Number(7_980, {
                #"total": 0.014j,
                "lumi_13p6TeV_correlated": 0.014j,
            })
        elif postfix == "postEE":
            cfg.x.luminosity = Number(26_671.7, {
                #"total": 0.014j,
                "lumi_13p6TeV_correlated": 0.014j,
            })
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2023:
        if postfix == "preBPix":
            cfg.x.luminosity = Number(17_794 , {
                "lumi_13TeV_correlated": 0.0j,
            })
        elif postfix == "postBPix":
            cfg.x.luminosity = Number(9_451, {
                "lumi_13TeV_correlated": 0.0j,
            })
            
    elif year == 2024:
        cfg.x.luminosity = Number(0, {
            "lumi_13TeV_correlated": 0.0j,
        })
        
    else:
        raise RuntimeError(f"Wrong year: {year}")


    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    # TODO later: Error for run three not available yet. Using error from run 2.
    cfg.x.minbias_xs = Number(69.2, 0.046j)

    year_postfix = ""
    if year == 2022:
        year_postfix = "EE" if postfix == "postEE" else ""
    elif year == 2023:
        year_postfix = "BPix" if postfix == "postBPix" else ""

        


    # --------------------------------------------------------------------------------------------- #
    # Adding triggers
    # --------------------------------------------------------------------------------------------- #

    if year == 2022:
        from httcp.config.triggers import add_triggers_run3_2022        
        add_triggers_run3_2022(cfg, postfix)
    elif year == 2023:
        from httcp.config.triggers import add_triggers_run3_2023
        add_triggers_run3_2023(cfg, postfix)
    elif year == 2024:
        raise RuntimeError("Babushcha: too early")
    else:
        raise RuntimeError(f"Wrong year: {year}. Check __init__.py in cmsdb campaign")


    # --------------------------------------------------------------------------------------------- #
    # Adding met filters
    # --------------------------------------------------------------------------------------------- #
    
    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)


    # --------------------------------------------------------------------------------------------- #
    # Adding external files e.g. goodlumi, normtag, SFs etc
    # restructure the postfix names to build appropriate tags
    # --------------------------------------------------------------------------------------------- #

    year2 = year%100
    #external_path_parent = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/corrections")
    external_path_parent = os.path.join(corrdir, "corrections")
    external_path_tail   = f"{year}{postfix}" if postfix else f"{year}"
    external_path        = os.path.join(external_path_parent, f"{external_path_tail}")
    normtag_path = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags"

    #print(f"external_path_parent : {external_path_parent}")
    #print(f"external_path        : {external_path}")

    #json_mirror = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/jsonpog-integration")
    json_mirror = os.path.join(corrdir, "jsonpog-integration")
    #print(f"json_mirror          : {json_mirror}")

    normtagjson = None
    if year == 2022:
        normtagjson = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json"
    elif year == 2023:
        normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
    elif year == 2024:
        raise RuntimeWarning("too early")
    else:
        raise RuntimeError(f"Check year : {year}")

    goldenjson  = glob(f"{external_path}/Lumi/*.json")[0]

    #print(f"GoldenJSON           : {goldenjson}")
    #print(f"NormtagJSON          : {normtagjson}")

    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden"  : (goldenjson,  "v1"),  # noqa
            "normtag" : (normtagjson, "v1"),
        },

        "pu_sf"             : (f"{json_mirror}/POG/LUM/{year}_Summer{year2}{year_postfix}/puWeights.json.gz",          "v1"), # PU
        "jet_jerc"          : (f"{json_mirror}/POG/JME/{year}_Summer{year2}{year_postfix}/jet_jerc.json.gz",           "v1"), # JEC
        "muon_sf"           : (f"{json_mirror}/POG/MUO/{year}_Summer{year2}{year_postfix}/muon_Z.json.gz",             "v1"), # Mu POG SF
        # https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossMuTauHlt.json?ref_type=heads
        "muon_xtrig_sf"     : (f"{json_mirror}/POG/MUO/{year}_Summer{year2}{year_postfix}/CrossMuTauHlt.json",         "v1"), # Mu xTrig SF
        "electron_sf"       : (f"{json_mirror}/POG/EGM/{year}_Summer{year2}{year_postfix}/electron.json.gz",           "v1"), # Ele POG SF
        "electron_trig_sf"  : (f"{json_mirror}/POG/EGM/{year}_Summer{year2}{year_postfix}/electronHlt.json.gz",        "v1"), # Ele HLT SF
        # https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossEleTauHlt.json?ref_type=heads
        "electron_xtrig_sf" : (f"{json_mirror}/POG/EGM/{year}_Summer{year2}{year_postfix}/CrossEleTauHlt.json",        "v1"), # Ele xTrig SF
        "tau_sf"            : (f"{json_mirror}/POG/TAU/{year}_{postfix}/tau_DeepTau2018v2p5_{year}_{postfix}.json.gz", "v1"), # TEC and ID SF
        # https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/ditaujet_jetleg_SFs_preEE.json?ref_type=heads
        "ditau_jet_trig_sf" : (f"{json_mirror}/POG/TAU/{year}_{postfix}/ditaujet_jetleg_SFs_preEE.json",               "v1"),
        "jet_veto_map"      : (f"{json_mirror}/POG/JME/{year}_Summer{year2}{year_postfix}/jetvetomaps.json.gz",        "v1"), # JetVeto
        "zpt_rewt_sf"       : (f"{external_path}/Zpt/myZptCorrections.json.gz",                                        "v1"), # Zpt Rewt
        "tautau_ff"         : (f"{external_path}/Zpt/myZptCorrections.json.gz",                                        "v1"), # DUMMY !!!
        #"btag_sf_corr": (f"{json_mirror}/POG/BTV/{year}_Summer{year2}{year_postfix}/btagging.json.gz",                "v1"),
        #"met_phi_corr": (f"{json_mirror}/POG/JME/2018_UL/met.json.gz",                                                "v1"), #met phi, unavailable Run3
    })


    # --------------------------------------------------------------------------------------------- #
    # electron settings
    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    # --------------------------------------------------------------------------------------------- #
    
    electron_sf_tag = ""
    if year == 2022:
        electron_sf_tag = "2022Re-recoE+PromptFG" if year_postfix else "2022Re-recoBCD"
    elif year == 2023:
        electron_sf_tag = "2023PromptC" if year_postfix else "2023PromptD"
    elif year == 2024:
        raise RuntimeWarning("too early")
    else:
        raise RuntimeError("wrong year")
    
    cfg.x.electron_sf_names = (
        "Electron-ID-SF",
        electron_sf_tag,
        "wp80iso",
    )
    cfg.x.electron_trig_sf_names = (
        "Electron-HLT-SF",
        electron_sf_tag,
        "HLT_SF_Ele30_TightID",
    )
    cfg.x.electron_xtrig_sf_names = (
        "Electron-HLT-SF",
        electron_sf_tag,
        "HLT_SF_Ele24_TightID",
    )
    
    
    # --------------------------------------------------------------------------------------------- #
    # muon settings
    # names of muon correction sets and working points
    # (used in the muon producer)
    # --------------------------------------------------------------------------------------------- #

    cfg.x.muon_id_sf_names  = (
        "NUM_MediumID_DEN_TrackerMuons",
        f"{year}_{postfix}"
    )
    cfg.x.muon_iso_sf_names = (
        "NUM_TightPFIso_DEN_MediumID",
        f"{year}_{postfix}"
    )
    cfg.x.muon_IsoMu24_trigger_sf_names = (
        "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium",
        f"{year}_{postfix}"
    )
    cfg.x.muon_xtrig_sf_names = (
        "NUM_IsoMu20_DEN_CutBasedIdMedium_and_PFIsoMedium",
        f"{year}_{postfix}"
    )
    

    # --------------------------------------------------------------------------------------------- #
    # jet settings
    # TODO: keep a single table somewhere that configures all settings: btag correlation, year
    #       dependence, usage in calibrator, etc
    # common jec/jer settings configuration
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107
    # --------------------------------------------------------------------------------------------- #

    if year == 2022:
        # TODO: check versions
        jec_campaign = f"Summer{year2}{year_postfix}_22Sep2023"
        jec_version = {2022: "V2"}[year]
        jer_campaign = f"Summer{year2}{year_postfix}_22Sep2023"
        jer_version = "JR" + {2022: "V1"}[year]
        jet_type = "AK4PFPuppi"
    else:
        assert False

    cfg.x.jec = DotDict.wrap({
        "campaign": jec_campaign,     # campaign name for this JEC corrections
        "version": jec_version,       # version of the corrections
        "jet_type": jet_type,         # Type of jets that the corrections should be applied on
        "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"], # relevant levels in the derivation process of the JEC 
        "levels_for_type1_met": ["L1FastJet"], # relevant levels in the derivation process of the Type 1 MET JEC
        "uncertainty_sources": [      # names of the uncertainties to be applied
            # "AbsoluteStat",
            # "AbsoluteScale",
            # "AbsoluteSample",
            # "AbsoluteFlavMap",
            # "AbsoluteMPFBias",
            # "Fragmentation",
            # "SinglePionECAL",
            # "SinglePionHCAL",
            # "FlavorQCD",
            # "TimePtEta",
            # "RelativeJEREC1",
            # "RelativeJEREC2",
            # "RelativeJERHF",
            # "RelativePtBB",
            # "RelativePtEC1",
            # "RelativePtEC2",
            # "RelativePtHF",
            # "RelativeBal",
            # "RelativeSample",
            # "RelativeFSR",
            # "RelativeStatFSR",
            # "RelativeStatEC",
            # "RelativeStatHF",
            # "PileUpDataMC",
            # "PileUpPtRef",
            # "PileUpPtBB",
            # "PileUpPtEC1",
            # "PileUpPtEC2",
            # "PileUpPtHF",
            # "PileUpMuZero",
            # "PileUpEnvelope",
            # "SubTotalPileUp",
            # "SubTotalRelative",
            # "SubTotalPt",
            # "SubTotalScale",
            # "SubTotalAbsolute",
            # "SubTotalMC",
            "Total",
            # "TotalNoFlavor",
            # "TotalNoTime",
            # "TotalNoFlavorNoTime",
            # "FlavorZJet",
            # "FlavorPhotonJet",
            # "FlavorPureGluon",
            # "FlavorPureQuark",
            # "FlavorPureCharm",
            # "FlavorPureBottom",
            #"CorrelationGroupMPFInSitu",
            #"CorrelationGroupIntercalibration",
            #"CorrelationGroupbJES",
            #"CorrelationGroupFlavor",
            #"CorrelationGroupUncorrelated",
        ],
    })

    # JER
    cfg.x.jer = DotDict.wrap({
        "campaign": jer_campaign,
        "version": jer_version,
        "jet_type": jet_type,
    })


    # --------------------------------------------------------------------------------------------- #
    # bjet settings
    # b-tag working points
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
    # TODO later: complete WP when data becomes available
    # --------------------------------------------------------------------------------------------- #

    btag_key = f"{year}{year_postfix}"
    cfg.x.btag_working_points = DotDict.wrap({
        "deepjet": {
            "loose"  : {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0583, "2023BPix": 0.0583, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.3086, "2023BPix": 0.3086, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.7183, "2022EE": 0.73,   "2023": 0.7183, "2023BPix": 0.7183, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.8111, "2023BPix": 0.8111, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.9512, "2023BPix": 0.9512, "2024": 0.0}[btag_key],
        },
        "robustParticleTransformer": {
            "loose"  : {"2022": 0.0849, "2022EE": 0.0897, "2023": 0.0849, "2023BPix": 0.0849, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.4319, "2022EE": 0.451,  "2023": 0.4319, "2023BPix": 0.4319, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.8482, "2022EE": 0.8604, "2023": 0.8482, "2023BPix": 0.8482, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.9151, "2022EE": 0.9234, "2023": 0.9151, "2023BPix": 0.9151, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9874, "2022EE": 0.9893, "2023": 0.9874, "2023BPix": 0.9874, "2024": 0.0}[btag_key],
        },
        "particleNet": {
            "loose"  : {"2022": 0.047,  "2022EE": 0.0499, "2023": 0.0499, "2023BPix": 0.0499, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.245,  "2022EE": 0.2605, "2023": 0.2605, "2023BPix": 0.2605, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.6915, "2023BPix": 0.6915, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.8033, "2023BPix": 0.8033, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.961,  "2022EE": 0.9664, "2023": 0.9664, "2023BPix": 0.9664, "2024": 0.0}[btag_key],
        },
    })
    # 2023 is dummy ... CORRECT IT LATER ===>> TODO


    # --------------------------------------------------------------------------------------------- #
    # tau settings
    # tau-id working points
    # --------------------------------------------------------------------------------------------- #

    cfg.x.deep_tau_tagger = "DeepTau2018v2p5"
    cfg.x.deep_tau_info = DotDict.wrap({
        "DeepTau2018v2p5": {
            "wp": {
                "vs_e": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
                "vs_m": {"VVVLoose": 1, "VVLoose": 1, "VLoose": 1, "Loose": 2, "Medium": 3, "Tight": 4, "VTight": 4, "VVTight": 4},
                "vs_j": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
            },
            "vs_e": {
                "etau"   : "Tight",
                "mutau"  : "VVLoose",
                "tautau" : "VVLoose",
            },
            "vs_m": {
                "etau"   : "Loose",
                "mutau"  : "Tight",
                "tautau" : "VLoose",
            },
            "vs_j": {
                "etau"   : "Tight",
                "mutau"  : "Medium",
                "tautau" : "Medium",
            },
        },
    })
    

    # --------------------------------------------------------------------------------------------- #
    # met settings
    # --------------------------------------------------------------------------------------------- #

    
    
    # --------------------------------------------------------------------------------------------- #
    # register shifts
    # --------------------------------------------------------------------------------------------- #

    cfg.add_shift(name="nominal", id=0)


    # --- >>> PU weight <<< --- #
    cfg.add_shift(name="minbias_xs_up", id=7, type="shape")
    cfg.add_shift(name="minbias_xs_down", id=8, type="shape")
    add_shift_aliases(
        cfg,
        "minbias_xs",
        {
            "pu_weight": "pu_weight_{name}",
            #"normalized_pu_weight": "normalized_pu_weight_{name}",
        },
    )
    
    # load jec sources
    with open(os.path.join(thisdir, "jec_sources.yaml"), "r") as f:
        all_jec_sources = yaml.load(f, yaml.Loader)["names"]

    for jec_source in cfg.x.jec.uncertainty_sources:
        idx = all_jec_sources.index(jec_source)
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5000 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5001 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        add_shift_aliases(
            cfg,
            f"jec_{jec_source}",
            {
                "Jet.pt": "Jet.pt_{name}",
                "Jet.mass": "Jet.mass_{name}",
                "MET.pt": "MET.pt_{name}",
                "MET.phi": "MET.phi_{name}",
            },
        )

    cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    add_shift_aliases(
        cfg,
        "jer",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            "MET.pt": "MET.pt_{name}",
            "MET.phi": "MET.phi_{name}",
        },
    )
    
    # --- >>> tau weight <<< --- #
    cfg.add_shift(name=f"tau_up", id=50, type="shape")
    cfg.add_shift(name=f"tau_down", id=51, type="shape")
    add_shift_aliases(
        cfg,
        "tau",
        {
            "tau_weight": "tau_weight_{direction}",
        },
    )
    #cfg.add_shift(name=f"tes_up", id=52, type="shape")
    #cfg.add_shift(name=f"tes_down", id=53, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "tes",
    #    {
    #        "tes_weight": "tes_weight_{direction}",
    #    },
    #)

    # --- >>> ele weight <<< --- #
    cfg.add_shift(name="e_up", id=90, type="shape")
    cfg.add_shift(name="e_down", id=91, type="shape")
    add_shift_aliases(
        cfg,
        "e",
        {
            "electron_idiso_weight": "electron_idiso_weight_{direction}",
        },
    )
    cfg.add_shift(name="e_trig_up", id=92, type="shape")
    cfg.add_shift(name="e_trig_down", id=93, type="shape")
    add_shift_aliases(
        cfg,
        "e_trig",
        {
            "electron_Ele30_WPTight_trigger_weight": "electron_Ele30_WPTight_trigger_weight_{direction}",
        },
    )
    cfg.add_shift(name="e_xtrig_up", id=94, type="shape")
    cfg.add_shift(name="e_xtrig_down", id=95, type="shape")
    add_shift_aliases(
        cfg,
        "e_xtrig",
        {
            "electron_xtrig_weight": "electron_xtrig_weight_{direction}",
        },
    )

    # --- >>> mu weight <<< --- #
    cfg.add_shift(name="mu_id_up", id=100, type="shape")
    cfg.add_shift(name="mu_id_down", id=101, type="shape")
    add_shift_aliases(
        cfg,
        "mu_id",
        {
            "muon_id_weight": "muon_id_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_iso_up", id=102, type="shape")
    cfg.add_shift(name="mu_iso_down", id=103, type="shape")
    add_shift_aliases(
        cfg,
        "mu_iso",
        {
            "muon_iso_weight": "muon_iso_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_trig_up", id=104, type="shape")
    cfg.add_shift(name="mu_trig_down", id=105, type="shape")
    add_shift_aliases(
        cfg,
        "mu_trig",
        {
            "muon_IsoMu24_trigger_weight": "muon_IsoMu24_trigger_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_xtrig_up", id=106, type="shape")
    cfg.add_shift(name="mu_xtrig_down", id=107, type="shape")
    add_shift_aliases(
        cfg,
        "mu_xtrig",
        {
            "muon_xtrig_weight": "muon_xtrig_weight_{direction}",
        },
    )

    # --- >>> pdf weight <<< --- #    
    cfg.add_shift(name="pdf_up", id=130, type="shape")
    cfg.add_shift(name="pdf_down", id=131, type="shape")
    add_shift_aliases(
        cfg,
        "pdf",
        {
            "pdf_weight": "pdf_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )

    # --- >>> tau spinner weight <<< --- #
    cfg.add_shift(name="tauspinner_up",   id=150, type="shape") # cp-even
    cfg.add_shift(name="tauspinner_down", id=151, type="shape") # cp-odd
    add_shift_aliases(
        cfg,
        "tauspinner",
        {
            "tauspinner_weight": "tauspinner_weight_{direction}",
        },
    )

    # --- >>> zpt weight <<< --- #    
    cfg.add_shift(name="zpt_up", id=160, type="shape")
    cfg.add_shift(name="zpt_down", id=161, type="shape")
    add_shift_aliases(
        cfg,
        "zpt",
        {
            "zpt_reweight": "zpt_reweight_{direction}",
        },
    )


    # --- >>> fake factor weight <<< --- # D U M M Y !!!!!!!!   
    cfg.add_shift(name="ff_up", id=170, type="shape")
    cfg.add_shift(name="ff_down", id=171, type="shape")
    add_shift_aliases(
        cfg,
        "ff",
        {
            "ff_weight": "ff_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )



    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    

    #---------------------------------------------------------------------------------------------#
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    # get_shifts = functools.partial(get_shifts_from_sources, cfg)
    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    #---------------------------------------------------------------------------------------------#

    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight"                  : [],
        "pu_weight"                             : [], #get_shifts("minbias_xs"),
        "electron_idiso_weight"                 : [], #get_shifts("e"),
        "electron_Ele30_WPTight_trigger_weight" : [], #get_shifts("e_xtrig"),
        "electron_xtrig_weight"                 : [], #get_shifts("e_xtrig"),
        "muon_id_weight"                        : [], #get_shifts("mu_id"),
        "muon_iso_weight"                       : [], #get_shifts("mu_iso"),
        "muon_IsoMu24_trigger_weight"           : [], #get_shifts("mu_trig"),
        "muon_xtrig_weight"                     : [], #get_shifts("mu_xtrig"),
        "tau_weight"                            : [], #get_shifts("tau"),
        #"ff_weight"                             : [],
        #"tes_weight"                           : [], #get_shifts("tes"),
        "tauspinner_weight"                     : get_shifts("tauspinner"),
        "pdf_weight"                            : [],
        "zpt_reweight"                          : [],
    })

    #---------------------------------------------------------------------------------------------#
    # No Idea
    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    #---------------------------------------------------------------------------------------------#
    
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev1'

    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }

    
    #---------------------------------------------------------------------------------------------#
    # Channels [possible to avoid as category id does the same]
    #---------------------------------------------------------------------------------------------#

    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    cfg.add_channel(name="tautau", id=4)


    #---------------------------------------------------------------------------------------------#
    # LFN settings
    # Campaigns must have creator names as desy or IPHC
    # the main path is in the campaign __init__, the basepath is mentioned in the dataset info
    #---------------------------------------------------------------------------------------------#

    campaign_tag = cfg.campaign.x("custom").get("creator")
    if campaign_tag == "desy" or campaign_tag == "IPHC":
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            print(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                print("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eoscms_redirector",
            )
            print(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eoscms_redirector"
        
    #---------------------------------------------------------------------------------------------#
    # Add categories described in categorization.py
    #---------------------------------------------------------------------------------------------#
        
    from httcp.config.categories import add_categories
    add_categories(cfg)

    cfg.x.ff_apply_id_map = DotDict.wrap({
        "etau"  : {},
        "mutau" : {},
        "tautau" : {
            "id"  : [cfg.get_category("tautau").id, cfg.get_category("real_1").id, cfg.get_category("hadC").id],  # category_id for AR C 
            "id0" : [cfg.get_category("tautau").id, cfg.get_category("real_1").id, cfg.get_category("hadC0").id], # category_id for AR C0
        },
    })
    
    #---------------------------------------------------------------------------------------------#
    # Add variables described in variables.py
    #---------------------------------------------------------------------------------------------#
    
    from httcp.config.variables import add_variables
    add_variables(cfg)

    #---------------------------------------------------------------------------------------------#
    # columns to keep after reduce events, MergeSelectionMasks and UniteColumns tasks
    #---------------------------------------------------------------------------------------------#
    
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # TauProds
            "TauProd.*",
            # general event info
            "run", "luminosityBlock", "event", "LHEPdfWeight",
            "PV.npvs","Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP",
            "trigger_ids",
            "single_triggered", "cross_triggered",
            "single_e_triggered", "cross_e_triggered",
            "single_mu_triggered", "cross_mu_triggered",
            "cross_tau_triggered", "cross_tau_jet_triggered",
            # --- new for categorization --- #
            "is_os", "is_iso_1", "is_iso_2", "is_low_mt", "is_b_veto",
            "is_real_1", "is_real_2", "is_fake_1", "is_fake_2",
            "is_lep_1", "is_pi_1", "is_pi_2", "is_rho_1", "is_rho_2",
            "is_a1_1pr_2pi0_1", "is_a1_1pr_2pi0_2",
            "is_a1_3pr_0pi0_1", "is_a1_3pr_0pi0_2",
            "is_a1_3pr_1pi0_1", "is_a1_3pr_1pi0_2",
            "is_ipsig_0to1_1",
            "has_0jet", "has_1jet", "has_2jet",
            #"LHE.NpLO", "LHE.NpNLO", "LHE.Vpt", "LHE.Njets",
        } | {
            f"GenPart.{var}" for var in [
                "pt", "eta", "phi", "mass",
                "status", "pdgId", "statusFlags",
            ]
        } | {
            f"GenZ.{var}" for var in [
                "pt", "eta", "phi", "mass",
            ]
        } | {
            f"PuppiMET.{var}" for var in [
                "pt_no_corr", "phi_no_corr",
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]
            #} | {
            #f"MET.{var}" for var in [
            #    "pt", "phi", "significance",
            #    "covXX", "covXY", "covYY",
            #]
        } | {
            f"Jet.{var}" for var in [
                "pt_no_corr", "eta_no_corr",
                "phi_no_corr", "mass_no_corr",
                "pt", "eta", "phi", "mass",
                "btagDeepFlavB", "hadronFlavour"
            ]
        } | {
            f"bJet.{var}" for var in [
                "pt", "eta", "phi", "mass",
                "btagDeepFlavB", "hadronFlavour"
            ]
        } | { # raw tau in events before any selection
            f"RawTau.{var}" for var in [
                "pt","eta","phi","mass",
                "decayMode",
                "rawIdx",
            ]
        } | { # taus after good selection
            f"GoodTau.{var}" for var in [
                "pt","eta","phi","mass",
                "decayMode",
                "rawIdx",
            ]
        } | { # taus in hcand 
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", "IPx", "IPy", "IPz",
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu",
                "decayMode",
                "decayModeHPS",
                "decayModePNet",
                "genPartFlav",
                "rawIdx",
                "pt_no_tes", "mass_no_tes"
            ]
        } | {
            f"TauSpinner.weight_cp_{var}" for var in [
                "0", "0_alt", "0p25", "0p25_alt", "0p375",
                "0p375_alt", "0p5", "0p5_alt", "minus0p25", "minus0p25_alt"
            ]
        } | { # muon before any selection 
            f"RawMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # good muons after selection
            f"GoodMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # to check if extra lep veto work correctly
            f"VetoMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # to check if dilep veto work correctly
            f"DoubleVetoMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # muons from hcand
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", "IPx", "IPy", "IPz",
                "decayMode", "pfRelIso04_all","mT", "rawIdx"
            ]
        } | { # eles before any selection
            f"RawElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # eles after good selection
            f"GoodElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # to check if extra lep veto work correctly
            f"VetoElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # to check if dilep veto work correctly
            f"DoubleVetoElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # electrons from hcand
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", "IPx", "IPy", "IPz",
                "decayMode", "pfRelIso03_all", "mT", "rawIdx",
                "deltaEtaSC",
            ]
        } | {
            f"TrigObj.{var}" for var in [
                "id", "pt", "eta", "phi", "filterBits",
            ]
        } | {
            f"hcand.{var}" for var in [
                "pt","eta","phi","mass", "charge",
                "decayMode", "rawIdx", "idVsJet",
                "IPx", "IPy", "IPz", "IPsig",
            ]
        } | {
            "GenTau.*", "GenTauProd.*",
        } | {
            f"hcandprod.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
                "pdgId", "tauIdx",
            ]
        } | {ColumnCollection.ALL_FROM_SELECTOR},
        "cf.MergeSelectionMasks": {
            "normalization_weight",
            "cutflow.*", "process_id", "category_ids",
        },
        "cf.UniteColumns": {
            "run","luminosityBlock","event","LHEPdfWeight",
            "PV.npvs","Pileup.nTrueInt","Pileup.nPU","genWeight",
            "single_triggered", "cross_triggered",
            "single_e_triggered", "cross_e_triggered",
            "single_mu_triggered", "cross_mu_triggered",
            "cross_tau_triggered", "cross_tau_jet_triggered",
            # --- new for categorization --- #
            "is_os", "is_iso_1", "is_iso_2", "is_low_mt", "is_b_veto",
            "is_real_1", "is_real_2", "is_fake_1", "is_fake_2",
            "is_lep_1", "is_pi_1", "is_pi_2", "is_rho_1", "is_rho_2",
            "is_a1_1pr_2pi0_1", "is_a1_1pr_2pi0_2",
            "is_a1_3pr_0pi0_1", "is_a1_3pr_0pi0_2",
            "is_a1_3pr_1pi0_1", "is_a1_3pr_1pi0_2",
            "is_ipsig_0to1_1",
            "has_0jet", "has_1jet", "has_2jet",
            "PuppiMET.*", "Jet.*",
            "TauSpinner.*", "hcand.*", "hcandprod.*"
            "GenTau.*", "GenTauProd.*",
            "process_id", "category_ids"
        },
    })

    
    #---------------------------------------------------------------------------------------------#
    # Helper switch for debugging
    #---------------------------------------------------------------------------------------------#

    cfg.x.verbose = DotDict.wrap({
        "calibration": {
            "main"                    : False,
            "tau"                     : False,
        },
        "selection": {
            "main"                    : False,
            "trigobject_matching"     : False,
            "extra_lep_veto"          : False,
            "dilep_veto"              : False,
            "higgscand"               : False,
        },        
    })


    cfg.x.extra_tags = DotDict.wrap({
        "genmatch"       : False,
    })


    cfg.x.is_channel_specific = False
    cfg.x.channel_specific_info = DotDict.wrap({
        "etau"   : False,
        "mutau"  : False,
        "tautau" : True,
    })
