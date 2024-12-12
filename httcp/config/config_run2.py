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

    print(f"\ncampaign: {campaign.name}")
    print(f"config name: {config_name}, id: {config_id}")
    # gather campaign data
    year = campaign.x.year
    postfix = campaign.x.postfix
    
    # some validations
    assert year in [2016, 2017, 2018]

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    
    # create a config by passing the campaign, so id and name will be identical
    cfg = ana.add_config(campaign,
                         name  = config_name,
                         id    = config_id)

    # add processes we are interested in
    process_names = [
        ## Data
        "data",
        ## Drell-Yan
        "dy_lep_m10to50",
        "dy_lep_m50",
        ## W + jets
        "w_lnu",
        ## TTJets
        "tt",
        ## Single top
        "st_tchannel",
        "st_twchannel",
        ## VV [diboson inclusive]
        "vv",
        ## Signal
        "h_ggf_tautau",
    ]

    for process_name in process_names:
        # development switch in case datasets are not _yet_ there
        if process_name not in procs:
            print(f"WARNING: {process_name} not in cmsdb processes")
            continue
        # add the process
        #if process_name == "h_ggf_tautau":
        #    procs.get(process_name).is_signal = True
        proc = cfg.add_process(procs.get(process_name))
        #if proc.name == "h_ggf_tautau":
        #    proc.is_signal = True
    # configuration of colors, labels, etc. can happen here
    from httcp.config.styles import stylize_processes
    stylize_processes(cfg)    

    # add datasets we need to study
    dataset_names = [
        ##Drell-Yan
        "dy_lep_m10to50_madgraph",
        "dy_lep_m50_madgraph",
        "dy_lep_m50_madgraph_ext1",
        ##W+jets
        "w_lnu_madgraph",
        ## ttbar
        "tt_sl",
        "tt_dl",
        "tt_fh",
        ##single top
        "st_tchannel_t",
        "st_tchannel_tbar",
        "st_tw_t",
        "st_tw_tb",
        ##Diboson
        "ww_incl",
        "wz_incl",
        "zz_incl",
        ## single top tW channel
        #"st_t_wminus_to_2l2nu",
        #"st_tbar_wplus_to_lnu2q",
        #"st_tbar_wplus_to_2l2nu",
        "h_ggf_tautau_prod_cp_even_sm",
        "h_ggf_tautau_prod_cp_odd_flat",
        "h_ggf_tautau_prod_cp_even_flat",
        "h_ggf_tautau_prod_max_mix_flat",
    ]
    
    datasets_data = []
    if year == 2016:
        if postfix == "preVFP":
            datasets_data = ["data_e_B",   "data_e_C",   "data_e_D",   "data_e_E",   "data_e_F",
                             "data_mu_B",  "data_mu_C",  "data_mu_D",  "data_mu_E",  "data_mu_F",
                             "data_tau_B", "data_tau_C", "data_tau_D", "data_tau_E", "data_tau_F"]

        elif postfix == "postVFP":
            datasets_data = ["data_e_F",   "data_e_G",   "data_e_H",
                             "data_mu_F",  "data_mu_G",  "data_mu_H",
                             "data_tau_F", "data_tau_G", "data_tau_H"]
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2017:
        datasets_data = ["data_e_B",   "data_e_C",   "data_e_D",   "data_e_E",   "data_e_F",
                         "data_mu_B",  "data_mu_C",  "data_mu_D",  "data_mu_E",  "data_mu_F",
                         "data_tau_B", "data_tau_C", "data_tau_D", "data_tau_E", "data_tau_F"]

    elif year == 2018:
        datasets_data = ["data_e_A",   "data_e_B",   "data_e_C",   "data_e_D",
                         "data_mu_A",  "data_mu_B",  "data_mu_C",  "data_mu_D",
                         "data_tau_A", "data_tau_B", "data_tau_C", "data_tau_D"]


    else:
        raise RuntimeError(f"Check Year in __init__.py in cmsdb campaign: {year}")


    dataset_names = datasets_data + dataset_names
    for dataset_name in dataset_names:
        # development switch in case datasets are not _yet_ there
        if dataset_name not in campaign.datasets:
            print(f"WARNING: {dataset_name} not in cmsdb campaign")
            continue
        
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if re.match(r"^(ww|wz|zz)_.*", dataset.name):
            dataset.add_tag("no_lhe_weights")
        
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator      = "main"
    cfg.x.default_selector        = "main"
    cfg.x.default_producer        = "main"
    cfg.x.default_ml_model        = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories      = ("incl",)
    #cfg.x.default_variables = ("n_jet", "jet1_pt")
    cfg.x.default_variables       = ("event","channel_id")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {}
    """
    cfg.x.process_groups = {
        "backgrounds": (backgrounds := [
            "w_lnu",
            "tt",
            "dy_lep_m50",
            "st",
            "vv",
        ]),
        "bkg_sig"       : (bkg_sig       := [*backgrounds, "h_ggf_tautau"]),
        "data_bkg"      : (data_bkg      := [*backgrounds, "data"]),
        "data_bkg_sig"  : (data_bkg_sig  := [*backgrounds, "data", "h_ggf_tautau"]),
    }
    cfg.x.process_settings_groups = {
        "unstack_processes": {proc: {"unstack": True, "scale": 10.0} for proc in ("h_ggf_tautau")},
    }
    cfg.x.general_settings_groups = {
        "compare_shapes": {"skip_ratio": False, "shape_norm": False, "yscale": "log"} #"cms_label": "simpw"},
    }
    """
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
                    "trigger",
                    "met_filter",
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
    
    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # difference pre-post VFP: https://cds.cern.ch/record/2854610/files/DP2023_006.pdf
    if year == 2016:
        if postfix == "preVFP":
            cfg.x.luminosity = Number(19_500, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
        elif postfix == "postVFP":
            cfg.x.luminosity = Number(16_800, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2017:
        cfg.x.luminosity = Number(41_480, {
            "lumi_13TeV_2017": 0.02j,
            "lumi_13TeV_1718": 0.006j,
            "lumi_13TeV_correlated": 0.009j,
        })
        
    elif year == 2018: # 59_830
        cfg.x.luminosity = Number(59_830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })

    else:
        raise RuntimeError(f"Wrong year: {year}")


    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    # TODO later: Error for run three not available yet. Using error from run 2.
    cfg.x.minbias_xs = Number(69.2, 0.046j)

    year_postfix = ""
    if year == 2016:
        year_postfix = "APV" if postfix == "preVFP" else ""
        
    # b-tag working points
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
    # TODO later: complete WP when data becomes available
    btag_key = f"{year}{year_postfix}"
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP?rev=6
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP?rev=8
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17?rev=15
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18?rev=18
    cfg.x.btag_working_points = DotDict.wrap({
        "deepjet": {
            "loose" : {"2016APV": 0.0508, "2016": 0.0480, "2017": 0.0532, "2018": 0.0490}[btag_key],
            "medium": {"2016APV": 0.2598, "2016": 0.2489, "2017": 0.3040, "2018": 0.2783}[btag_key],
            "tight" : {"2016APV": 0.6502, "2016": 0.6377, "2017": 0.7476, "2018": 0.7100}[btag_key],
        },
        "deepcsv": {
            "loose" : {"2016APV": 0.2027, "2016": 0.1918, "2017": 0.1355, "2018": 0.1208}[btag_key],
            "medium": {"2016APV": 0.6001, "2016": 0.5847, "2017": 0.4506, "2018": 0.4168}[btag_key],
            "tight" : {"2016APV": 0.8819, "2016": 0.8767, "2017": 0.7738, "2018": 0.7665}[btag_key],
        },
    })

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


  
    # Adding triggers
    if year == 2016:
        from httcp.config.triggers import add_triggers_UL2017
        add_triggers_UL2016(cfg, postfix)
    elif year == 2017:
        from httcp.config.triggers import add_triggers_UL2017
        add_triggers_UL2017(cfg)
    elif year == 2018:
        from httcp.config.triggers import add_triggers_UL2018
        add_triggers_UL2018(cfg)
    else:
        raise RuntimeError(f"Wrong year: {year}. Check __init__.py in cmsdb campaign")


    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    year2 = year%100
    #external_path_parent = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/corrections")
    external_path_parent = os.path.join(corrdir, "corrections")
    external_path_tail   = f"{year}{postfix}" if postfix else f"{year}"
    external_path        = os.path.join(external_path_parent, f"{external_path_tail}")

    #print(f"external_path_parent : {external_path_parent}")
    #print(f"external_path        : {external_path}")

    #json_mirror = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/jsonpog-integration")
    json_mirror = os.path.join(corrdir, "jsonpog-integration")
    #print(f"json_mirror          : {json_mirror}")

    normtagjson = None
    goldenjson  = None
    if year == 2016:
        normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
        goldenjson  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
    elif year == 2017:
        normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
        goldenjson  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
    elif year == 2018:
        normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
        goldenjson  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"        
    else:
        raise RuntimeError(f"Check year : {year}")

    #print(f"GoldenJSON           : {goldenjson}")
    #print(f"NormtagJSON          : {normtagjson}")

    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden"  : (goldenjson,  "v1"),  # noqa
            "normtag" : (normtagjson, "v1"),
        },

        # pileup weight corrections
        "pu_sf": (f"{json_mirror}/POG/LUM/{year}{postfix}_UL/puWeights.json.gz", "v1"),

        # jet energy correction
        "jet_jerc": (f"{json_mirror}/POG/JME/{year}{postfix}_UL/jet_jerc.json.gz", "v1"),

        # Add Muon POG scale factors
        "muon_sf": (f"{json_mirror}/POG/MUO/{year}{postfix}_UL/muon_Z.json.gz", "v1"),
        
        # electron scale factors
        "electron_sf": (f"{json_mirror}/POG/EGM/{year}{postfix}_UL/electron.json.gz", "v1"),
        
        # tau energy correction and scale factors
        #"tau_sf": (f"{external_path_parent}/tau_DeepTau2018v2p5_2022_preEE.json.gz", "v1"),  # noqa
        "tau_correction": (f"{json_mirror}/POG/TAU/{year}{postfix}_UL/tau.json.gz", "v1"),  # noqa

        # met phi corrector 
        # unavailable for Run3
        # "met_phi_corr": (f"{json_mirror}/POG/JME/2018_UL/met.json.gz", "v1"),
    })

    # register shifts
    cfg.add_shift(name="nominal", id=0)

    """
    # ----------------------------------------------------------------------- #
    #                       jec & jer configuration                           #
    # ----------------------------------------------------------------------- #
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    # TODO later: check this corrections summary correction_file (jet_jerc.json.gz) after setting sandbox_dev
    # https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME?ref_type=heads
    cfg.x.jec = DotDict.wrap({
        "campaign": f"Summer{year2}{year_postfix}_22Sep2023",                  # campaign name for this JEC corrections
        "version": {2022: "V2", 2023: "V1"}[year],                             # version of the corrections 
        "jet_type": "AK4PFPuppi",                                              # Type of jets that the corrections should be applied on
        "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],   # relevant levels in the derivation process of the JEC
        "levels_for_type1_met": ["L1FastJet"],                                 # relevant levels in the derivation process of the Type 1 MET JEC
        "uncertainty_sources": [                                               # names of the uncertainties to be applied
            "Total",
            "CorrelationGroupMPFInSitu",
            "CorrelationGroupIntercalibration",
            "CorrelationGroupbJES",
            "CorrelationGroupFlavor",
            "CorrelationGroupUncorrelated",
        ],
    })

    # JER
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107 # TODO later: check this
    cfg.x.jer = DotDict.wrap({
        "campaign": f"Summer{year2}{year_postfix}_22Sep2023",
        "version": "JR" + {2022: "V1", 2023: "V1"}[year],
        "jet_type": "AK4PFPuppi",
    })
    """

    # ----------------------------------------------------------------------- #
    #                                electron id SF                           #
    # ----------------------------------------------------------------------- #
    # names of electron correction sets and working points
    # (used in the electron_sf producer)    
    electron_sf_tag = ""
    if year == 2022:
        electron_sf_tag = "2022Re-recoE+PromptFG" if year_postfix else "2022Re-recoBCD"
    elif year == 2023:
        electron_sf_tag = "2023PromptC" if year_postfix else "2023PromptD"
    
    cfg.x.electron_sf_names = (
        "UL-Electron-ID-SF",
        f"{year}{postfix}",
        "wp80iso",
    )

    # ----------------------------------------------------------------------- #
    #                                muon id SF                               #
    # ----------------------------------------------------------------------- #
    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = (
        "NUM_TightRelIso_DEN_MediumPromptID", 
        f"{year}{postfix}_UL",
    )

    """
    # register shifts
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="minbias_xs_up", id=7, type="shape")
    cfg.add_shift(name="minbias_xs_down", id=8, type="shape")
    add_shift_aliases(
        cfg,
        "minbias_xs",
        {
            "pu_weight": "pu_weight_{name}",
            "normalized_pu_weight": "normalized_pu_weight_{name}",
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
    """
    """
    for i, (match, dm) in enumerate(itertools.product(["jet", "e"], [0, 1, 10, 11])):
        cfg.add_shift(name=f"tec_{match}_dm{dm}_up", id=20 + 2 * i, type="shape", tags={"tec"})
        cfg.add_shift(name=f"tec_{match}_dm{dm}_down", id=21 + 2 * i, type="shape", tags={"tec"})
        add_shift_aliases(
            cfg,
            f"tec_{match}_dm{dm}",
            {
                "Tau.pt": "Tau.pt_{name}",
                "Tau.mass": "Tau.mass_{name}",
                "MET.pt": "MET.pt_{name}",
                "MET.phi": "MET.phi_{name}",
            },
        )
    """
    """
    # start at id=50
    cfg.x.tau_unc_names = [
        "jet_dm0", "jet_dm1", "jet_dm10",
        "e_barrel", "e_endcap",
        "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
    ]
    for i, unc in enumerate(cfg.x.tau_unc_names):
        cfg.add_shift(name=f"tau_{unc}_up", id=50 + 2 * i, type="shape")
        cfg.add_shift(name=f"tau_{unc}_down", id=51 + 2 * i, type="shape")
        add_shift_aliases(cfg, f"tau_{unc}", {"tau_weight": f"tau_weight_{unc}_" + "{direction}"})
    """
    cfg.add_shift(name="e_up", id=90, type="shape")
    cfg.add_shift(name="e_down", id=91, type="shape")
    add_shift_aliases(cfg, "e", {"electron_weight": "electron_weight_{direction}"})

    cfg.add_shift(name="mu_up", id=100, type="shape")
    cfg.add_shift(name="mu_down", id=101, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

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


    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    #from httcp.config.variables import keep_columns
    #keep_columns(cfg)


    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    # get_shifts = functools.partial(get_shifts_from_sources, cfg)
    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        #"pu_weight"           : [], #get_shifts("minbias_xs"),
        #"electron_weight"     : [], #get_shifts("e"),
        #"muon_weight"         : [], #get_shifts("mu"),
        #"tau_weight"          : [], #get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
    })
    # define per-dataset event weights
    for dataset in cfg.datasets:
        if dataset.x("no_lhe_weights", False):
            dataset.x.event_weights = {
                #"pdf_weight": [], #get_shifts("pdf"),
            }
    cfg.x.default_weight_producer = "all_weights"

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
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
    # channels
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    cfg.add_channel(name="tautau", id=4)
    
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
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from httcp.config.categories import add_categories
    add_categories(cfg)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # TauProds                                                                                                            
            "TauProd.*",
            # general event info                                                                                                  
            "run", "luminosityBlock", "event", "LHEPdfWeight",
            "PV.npvs","Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP",
            "trigger_ids",
        } | {
            f"PuppiMET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]
        } | {
            f"MET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]
        } | {
            f"Jet.{var}" for var in [
                "pt", "eta", "phi", "mass",
                "btagDeepFlavB", "hadronFlavour"
            ]
        } | {
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", #"IPx", "IPy", "IPz",
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu",
                "decayMode",
                "decayModeHPS",
                "decayModePNet",
                "genPartFlav",
                "rawIdx",
                "pt_no_tes", "mass_no_tes"
            ]
            #} | {
            #f"TauSpinner.weight_cp_{var}" for var in [
            #    "0", "0_alt", "0p25", "0p25_alt", "0p375",
            #    "0p375_alt", "0p5", "0p5_alt", "minus0p25", "minus0p25_alt"
            #]
        } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", #"IPx", "IPy", "IPz",
                "decayMode", "pfRelIso04_all","mT", "rawIdx"
            ]
        } | {
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", #"IPx", "IPy", "IPz",
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
                "decayMode", "rawIdx"
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
            "*",
        },
    })

    # For debugging
    cfg.x.verbose = DotDict.wrap({
        "selection": {
            "main"                    : True,
            "trigobject_matching"     : False,
            "extra_lep_veto"          : False,
            "dilep_veto"              : False,
        },
    })
