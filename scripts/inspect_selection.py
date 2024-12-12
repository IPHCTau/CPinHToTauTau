import os
import re
import sys
import pickle
import awkward as ak
import numpy as np
from tabulate import tabulate
from IPython import embed

basepath = "/eos/user/g/gsaha/cf_store/analysis_httcp/cf.SelectEvents"
#dataset = "h_ggf_tautau_prod_cp_even_sm"
datasets = ["h_ggf_tautau_prod_cp_even_sm", "dy_lep_m50_madgraph"]
version = "v3"
campaign = "run2_UL2018_nano_cp_tau_v09_limited"

#file = f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/results_0.parquet"
#cols = f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/columns_0.parquet"


#files = [f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/results_0.parquet" for dataset in datasets]
#cols  = [f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/columns_0.parquet" for dataset in datasets]

selections = [
    'json',
    'trigger',
    'met_filter',
    'b_veto',
    'dilepton_veto',
    'has_2_or_more_leps_with_at_least_1_tau',
    'has_at_least_1_pair_before_trigobj_matching',
    'has_at_least_1_pair_after_trigobj_matching',
    'extra_lepton_veto',
    'One_higgs_cand_per_event',
    'has_proper_tau_decay_products',
]




for dataset in datasets:
    print(f"dataset: {dataset}")
    file = f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/results_0.parquet"
    cols = f"{basepath}/{campaign}/{dataset}/nominal/calib__main/sel__main/{version}/columns_0.parquet"

    results = ak.from_parquet(file)
    columns = ak.from_parquet(cols)
    
    steps = results.steps
    events = results.event
    objects = results.objects
    
    nselevnents = ak.sum(events)
    
    initial = len(objects.Muon)
    count = ak.Array(np.ones(initial, dtype=bool))
    initial = ak.sum(count)
    headers = ["selections", "nevents", "abs eff"]
    rows = []
    rows.append(["Initial", initial, 1.0])
    for i in range(len(selections)):
        sel = selections[i]
        count = count & steps[sel]
        nsel = ak.sum(count)
        abseff = round(nsel/initial, 3)
        rows.append([sel, nsel, abseff])
        
    rows = np.array(rows)
    evts = rows[:,1]
    evts_den = np.array([initial] + evts[:-1].tolist())
    rel_eff = np.round(np.asarray(evts, float)/np.asarray(evts_den, float), decimals=3)
    
    rows = np.concatenate([rows, rel_eff[:,None]], axis=1).tolist()
    headers.append("rel eff")
    
    table = tabulate(rows, headers, tablefmt="pretty")
    print(f"CutFlow :\n{table}")
    print(f"Events selected at the end (all categories) : {nselevnents} [ {round(100*nselevnents/initial,6)}% ]")
    
    etau_mask = columns.channel_id == 1
    mutau_mask = columns.channel_id == 2
    tautau_mask = columns.channel_id == 4
    other_ch_mask = ~(etau_mask | mutau_mask | tautau_mask)
    
    etau_ev = ak.sum(events[etau_mask])
    mutau_ev = ak.sum(events[mutau_mask])
    tautau_ev = ak.sum(events[tautau_mask])
    total = etau_ev + mutau_ev + tautau_ev
    print(f"Channel wise nEvents ===>")
    print(f"etau   : {etau_ev}\t[ {round(100*etau_ev/initial,6)}% ]")
    print(f"mutau  : {mutau_ev}\t[ {round(100*mutau_ev/initial,6)}% ]")
    print(f"tautau : {tautau_ev}\t[ {round(100*tautau_ev/initial,6)}% ]")
    print(f"etau + mutau + tautau : {etau_ev + mutau_ev + tautau_ev}")
    print(f"In other channels  : {ak.sum(events[other_ch_mask])}\n\n")
    
