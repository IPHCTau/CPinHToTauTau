import os
import re
import sys
import pickle
import numpy as np
import awkward as ak
from tabulate import tabulate
from IPython import embed

lumi = 7.98
categories_to_check = [1,101,201,301]
isMC = True
campaign      = "run2_UL2018_nano_cp_tau_v09"
version       = "Run2_2018_full_20240830_181952"
cutflowpath   = "/eos/user/g/gsaha/cf_store/analysis_httcp/cf.CreateCutflowHistograms"
mergehistpath = "/eos/user/g/gsaha/cf_store/analysis_httcp/cf.MergeHistograms"

filesMC = {
    "signal"     : [f"{cutflowpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    3.274821, #0.0142, #3.274821,
                    175963.0],
    "wj_incl"    : [f"{cutflowpath}/{campaign}/w_lnu_madgraph/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/w_lnu_madgraph/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1, #44.7, #63425.1,
                    82442496.0],
    "dy_lep_m50" : [f"{cutflowpath}/{campaign}/dy_lep_m50_madgraph_ext1/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/dy_lep_m50_madgraph_ext1/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    6282.6, #3.59, #6282.6,
                    101415750.0],
    "ww"         : [f"{cutflowpath}/{campaign}/ww_incl/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/ww_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    122.3, #0.453, #122.3,
                    15679000.0],
    "wz"         : [f"{cutflowpath}/{campaign}/wz_incl/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wz_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    41.1, #0.352, #41.1,
                    7940000.0],
    "zz"         : [f"{cutflowpath}/{campaign}/zz_incl/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/zz_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    19.4, #0.271, #19.4,
                    3526000.0],
    "tt_dl"      : [f"{cutflowpath}/{campaign}/tt_dl/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_dl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    98.04,
                    144831008.0],
    "tt_sl"      : [f"{cutflowpath}/{campaign}/tt_sl/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_sl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    405.75,
                    475111154.0],
    "tt_fh"      : [f"{cutflowpath}/{campaign}/tt_fh/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_fh/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    419.80,
                    340474926.0],
}

filesDATA = {
    "data_e_A"   : [f"{cutflowpath}/{campaign}/data_e_A/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_e_A/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    1000, 175963.0],
    "data_e_B"   : [f"{cutflowpath}/{campaign}/data_e_B/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_e_B/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1, 59463375.0],
    "data_e_C"   : [f"{cutflowpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1, 59463375.0],
    #"data_e_D"   : [f"{cutflowpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
    #                f"{mergehistpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
    #                63425.1, 59463375.0],
    "data_mu_A"  : [f"{cutflowpath}/{campaign}/data_mu_A/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_mu_A/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    1000, 175963.0],
    "data_mu_B"  : [f"{cutflowpath}/{campaign}/data_mu_B/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_mu_B/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1, 59463375.0],
    "data_mu_C"  : [f"{cutflowpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    1000, 175963.0],
    "data_mu_D"  : [f"{cutflowpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1, 59463375.0],
    "data_tau_A" : [f"{cutflowpath}/{campaign}/data_tau_A/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_tau_A/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    1000, 175963.0],
    "data_tau_B": [f"{cutflowpath}/{campaign}/data_tau_B/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                   f"{mergehistpath}/{campaign}/data_tau_B/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                   63425.1, 59463375.0],    
    "data_tau_C" : [f"{cutflowpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    1000, 175963.0],
    "data_tau_D": [f"{cutflowpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                   f"{mergehistpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                   63425.1, 59463375.0],    
}

files = filesMC if isMC else filesDATA


selections = [
    'Initial',
    'json',
    'trigger',
    'met_filter',
    'b_veto',
    'has_2_or_more_leps_with_at_least_1_tau',
    'dilepton_veto',
    'has_at_least_1_pair_before_trigobj_matching',
    'has_at_least_1_pair_after_trigobj_matching',
    'extra_lepton_veto',
    'One_higgs_cand_per_event',
    'has_proper_tau_decay_products'
]

file_names = list(files.keys())
files_count_dict = {}
for name, filelist in files.items():
    file = filelist[0]
    print(f"File: {file}")
    # load pickle file
    file_obj = open(file, 'rb')
    data = pickle.load(file_obj)
    #from IPython import embed; embed()
    axes = data.axes
    category_axis  = axes['category']
    
    category_idxs  = [category_axis.index(catid) for catid in categories_to_check]
    selection_axis = axes['step']
    
    count_dict = {key:[] for key in categories_to_check}
    for i,catidx in enumerate(category_idxs):
        catid = categories_to_check[i]
        print(f"categoryID : {catid}")
        for j in range(len(selections)):
            count = data[catidx, :, j, :, :].values()[0,0,0]
            count_dict[catid].append(count)
            print(f"\t {selections[j]} : {count}")

    files_count_dict[name] = count_dict

embed()
#print(files_count_dict)

rev_dict = {key:{} for key in categories_to_check}
counts_from_merged_hist = []
for file_name, count_dict in files_count_dict.items():
    print(f"file : {file_name}")
    mergehistfile = files[file_name][1]
    fptr = open(mergehistfile, 'rb')
    merge_data = pickle.load(fptr)
    
    xsec = files[file_name][2]
    nev  = files[file_name][-1]
    normwt = (xsec * lumi * 1000)/nev
    print(normwt)
    temp_merged = []
    for key, val in count_dict.items():
        cat_idx = merge_data.axes['category'].index(key)
        count_from_merged_hist = merge_data[cat_idx,:,:,:].values()[0,0,0]
        count_from_merged_hist = count_from_merged_hist * normwt if not re.match(r"^data_*", file_name) else count_from_merged_hist
        temp_merged.append(count_from_merged_hist)
        #print(f"cat id: {key}")
        counts = ak.Array(val)
        # abs eff
        abs_eff = ak.values_astype(counts/(val[0]*ak.ones_like(counts)), np.float32)
        # rel eff
        den = ak.Array([val[0]]+val[:-1])
        rel_eff = ak.values_astype(counts/den, np.float32)

        #count_1 = (counts[0] * xsec * lumi * 1000)/nev if not re.match(r"^data_*", file_name) else counts[0]
        #count_2 = (counts[-1] * xsec * lumi * 1000)/nev if not re.match(r"^data_*", file_name) else counts[-1]
        count_1 = counts[0] * normwt if not re.match(r"^data_*", file_name) else counts[0]
        count_2 = counts[-1] * normwt if not re.match(r"^data_*", file_name) else counts[-1]

        abs_eff = ak.to_list(abs_eff)
        rel_eff = ak.to_list(rel_eff)
        #ext_abs_eff = np.round(np.array([count_1] + abs_eff + [count_2]), decimals=3)[:,None]
        #rev_dict[key][file_name] = ext_abs_eff

        ext_rel_eff = np.round(np.array([count_1] + rel_eff + [count_2]), decimals=3)[:,None]
        rev_dict[key][file_name] = ext_rel_eff
        
    counts_from_merged_hist.append(temp_merged)

#print(rev_dict)

#print(counts_from_merged_hist)
counts_from_merged_hist = np.round(np.array(counts_from_merged_hist), decimals=3)
#print(counts_from_merged_hist)
#from IPython import embed; embed()
counts_from_merged_hist = counts_from_merged_hist.T

headers = ["selections"] + file_names

for j, (key, val) in enumerate(rev_dict.items()):
    print(f"cat id: {key}")
    sels = ak.to_numpy(ak.Array(["Begins with (LumiScaled)"] + selections + ["Ends with (LumiScaled)"]))[:,None]
    for file, arr in val.items():
        sels = np.concatenate([sels, arr], axis=1)
    #from IPython import embed; embed()
    y = sels[:,1:][-1]
    y = y.astype(np.float64)
    total_y = round(np.sum(y), 2)
    sels = sels.tolist()
    sels.append(["From Merged Hist"]+counts_from_merged_hist[j].tolist())
    z = round(np.sum(counts_from_merged_hist[j]), 2)
    table = tabulate(sels, headers, tablefmt="pretty")
    #from IPython import embed; embed()
    print(table)
    print(f"Total Yield      : {total_y}")
    print(f"From Merged Hist : {z}\n")




"""
files = {
    "signal"     : [f"{cutflowpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    3.274821,
                    175963.0],
    "wj_incl"    : [f"{cutflowpath}/{campaign}/wj_incl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wj_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1,
                    10626959.0],
    "dy_lep_m50" : [f"{cutflowpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    6282.6,
                    145286646.0],
    "ww"         : [f"{cutflowpath}/{campaign}/ww/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/ww/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    122.3,
                    15405496.0],
    "wz"         : [f"{cutflowpath}/{campaign}/wz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    41.1,
                    7527043.0],
    "zz"         : [f"{cutflowpath}/{campaign}/zz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/zz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    19.4,
                    1181750.0],
    "tt_dl"      : [f"{cutflowpath}/{campaign}/tt_dl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_dl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    98.04,
                    21890633.0],
    "tt_sl"      : [f"{cutflowpath}/{campaign}/tt_sl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_sl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    405.75,
                    133416781.0],
    "tt_fh"      : [f"{cutflowpath}/{campaign}/tt_fh/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_fh/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    419.80,
                    42646240.0],
}

"""
"""
files = {
    "data_e_C"       : [f"{cutflowpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_e_D"       : [f"{cutflowpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_mu_C"      : [f"{cutflowpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_mu_D"      : [f"{cutflowpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_single_mu_C" : [f"{cutflowpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                          f"{mergehistpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                          63425.1, 59463375.0],
    "data_tau_C"     : [f"{cutflowpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_tau_D"     : [f"{cutflowpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
}
"""
