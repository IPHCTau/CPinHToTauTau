from typing import Optional

import law
#import matplotlib.pyplot as plt
#import mplhep #as hep
from columnflow.util import maybe_import
from columnflow.selection import SelectionResult

mpl = maybe_import("matplotlib")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")
np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

logger = law.logger.get_logger(__name__)

def get_event_level_eff(events, results, dataset_name):
    from tabulate import tabulate
    steps_ = results.steps.keys()
    indiv_selections_ = []
    comb_selections_ = []
    indiv_headers_ = ["selections", "nevents", "abs eff"]
    comb_headers_ = ["selections", "nevents", "abs eff", "rel eff"]
    init = len(events)
    comb_count_array = ak.Array(np.ones(init, dtype=bool))
    for step_ in steps_:
        count_array      = results.steps[step_]
        comb_count_array = comb_count_array & count_array
        count            = ak.sum(count_array)
        comb_count       = ak.sum(comb_count_array)
        indiv_selections_.append([step_, count, round(count/init,3)])
        comb_selections_.append([step_, comb_count, round(comb_count/init,3)])
    indiv_table_ = tabulate(indiv_selections_,
                            indiv_headers_,
                            tablefmt="pretty")
    logger.info(f"---> Efficiencies of individual selections: \n{indiv_table_}")
    # save the table as a .txt
    save_efficiency_table(indiv_table_, dataset_name, "Individual Selection Efficiencies")
    # plot the table as stair plot
    plot_eff(indiv_selections_, indiv_headers_, dataset_name, title="Individual Selection Efficiencies")

    comb_selections_ = np.array(comb_selections_)
    comb_selections_counts_ = comb_selections_[:,1]
    comb_den_ = np.array([init] + comb_selections_counts_[:-1].tolist())
    rel_eff_ = np.round(np.asarray(comb_selections_counts_, float)/np.asarray(comb_den_, float), decimals=3)
    comb_selections_ = np.concatenate([comb_selections_, rel_eff_[:,None]], axis=1).tolist()
    comb_table_ = tabulate(comb_selections_,
                           comb_headers_,
                           tablefmt="pretty")
    logger.info(f"---> Efficiencies of combined selections: \n{comb_table_}")
    

def get_object_eff(results,
                   tag,
                   key : Optional[str]=None,
                   final_mask : Optional[ak.Array]=None):
    from tabulate import tabulate
    logger.info(f"{tag}")
    aux = results.aux
    if key:
        logger.info(f" --- {key} --- ")
        aux = aux[f"{tag}_{key}"]
    #keys = [key for key in results.aux.keys() if key.startswith(f"{tag}_")]
    keys = [key for key in aux.keys() if key.startswith(f"{tag}_")]
    rows = []
    rows_evt_level = []

    n0 = 0
    nevt0 = 0

    for i, key in enumerate(keys):
        #mask = results.aux[key]
        mask = aux[key]
        mask_dummy = mask[:,:0]
        #from IPython import embed; embed()
        if final_mask:
            mask = ak.where(final_mask, mask, mask_dummy)
            #mask = mask & final_mask
        n = ak.sum(ak.sum(mask, axis=1))
        nevt = ak.sum(ak.any(mask, axis=1))
        if i == 0:
            n0 = n
            nevt0 = nevt
        rows.append([key, n, round(n/n0, 3)])
        rows_evt_level.append([key, nevt, round(nevt/nevt0, 3)])
    table = tabulate(rows, ["selection", f"n_{tag}", "abseff"], tablefmt="pretty")
    evt_table = tabulate(rows_evt_level, ["selection", f"n_{tag}", "abseff"], tablefmt="pretty")
    
    logger.info(f"object level : \n{table}")
    logger.info(f"event level  : \n{evt_table}")
    # save the table as a .txt file
    save_efficiency_table(table, dataset_name, title=f"Object-Level Efficiency: {tag} {key}")
    save_efficiency_table(evt_table, dataset_name, title=f"Event-Level Efficiency: {tag} {key}")
    # plot the table as stair plot
    plot_eff(rows, ["Selection", f"n_{tag}", "Abs Eff"], dataset_name, title=f"Object-Level Efficiency: {tag} {key}")
    plot_eff(rows_evt_level, ["Selection", f"n_{tag}", "Abs Eff"], dataset_name, title=f"Event-Level Efficiency: {tag} {key}")

    


def plot_eff(selections, headers, dataset_name, title="Selection Yield"): #self: Selector,
    """
    Plot efficiencies in CMS style.
    
    selections : list of lists,     Contains selection names, event counts, and efficiencies.
    headers : list of str,          Column headers where the first element is selection names.
    title : str,                    Title of the plot.
    """
    #dataset_name = self.dataset_inst.name
    mplhep.style.use("CMS")
    
    # Extract data
    step_labels = [row[0] for row in selections]  # Selection names
    event_counts = [row[1] for row in selections]  # Event counts
    efficiencies = [float(row[2]) for row in selections]
    #efficiencies = [row[2] for row in selections]  # Absolute efficiencies
    
    x_positions = np.arange(len(step_labels))
    
    # Create figure
    fig, ax1 = plt.subplots(figsize=(10.7, 10.7)) 
    
    # Plot event counts as stair plot
    edges = np.arange(len(event_counts) + 1) - 0.5
    ax1.stairs(event_counts, edges = edges , color="blue", label="Number of Events")
    ax1.set_ylabel("Event Count")
    ax1.tick_params(axis="y")
    
    # Add a second axis for efficiency
    ax2 = ax1.twinx()
    ax2.plot(x_positions, efficiencies, marker="o", color="red", linestyle="-", label="Efficiency")
    ax2.set_ylabel("Efficiency")
    ax2.tick_params(axis="y")
    ax2.set_ylim(0, 1.1)  # Efficiency should be between 0 and 1
    
    # Labeling
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(step_labels, rotation=45, ha="right")
    ax1.set_title(title, pad=50)
    
    # CMS label
    mplhep.cms.label("Private Work", loc=0)

    # Get handles and labels from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    # Combine them into one legend
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    
    # **Add efficiency values on top of each step**
    for x, eff in zip(x_positions, efficiencies):
        ax2.text(x, eff + 0.05, f"{eff:.2%}", ha="center", va="bottom", fontsize=10, color="red")


    fig.tight_layout()
    fig.savefig(f"/eos/user/m/mwitt/CPinHToTauTauOutput/debug/{dataset_name}_{title}.pdf", dpi=300)
    
    return fig

def save_efficiency_table(table, dataset_name, title="Selection Yield"):
    file_path = f"/eos/user/m/mwitt/CPinHToTauTauOutput/debug/{dataset_name}_{title}.txt"

    with open(file_path, "w") as f: # w = write
        f.write(table)
        logger.info(f"Efficiencies saved to {file_path}")


def debug_main(events, results, triggers, dataset_name, **kwargs):
    
    logger.info(f"---> ################### Inspecting event selections ################### <---\n")
    get_event_level_eff(events, results)

    _steps = list(results.steps.keys())
    final_mask = results.steps[_steps[0]]
    for key in _steps:
        final_mask = final_mask & results.steps[key]
    
    logger.info(f"---> ################### Inspecting trigger selections ################### <---\n")
    
    from tabulate import tabulate        
    # trigger details
    trigger_names = results.aux["trigger_names"]
    trigger_ids   = results.aux["trigger_ids"]
    
    #HLT_names = [trigger.name for trigger in self.config_inst.x.triggers]
    #HLT_ids   = [trigger.id for trigger in self.config_inst.x.triggers]
    HLT_names = [trigger.name for trigger in triggers]
    HLT_ids   = [trigger.id for trigger in triggers]
    
    trig_info = []
    for i in range(len(HLT_ids)):
        HLT = HLT_names[i]
        ID  = HLT_ids[i]
        nPassed = ak.sum(ak.any(trigger_ids == ID, axis=1))
        trig_info.append([HLT, ID, nPassed, round(nPassed/len(events), 3)])
    trig_table = tabulate(trig_info, ["HLT", "ID", "nEvents Passed", "Efficiency"], tablefmt="pretty")
    logger.info(trig_table)
    
    logger.info(f"---> ################### Inspecting object selections ################### <---")
    # muon
    get_object_eff(results, "muon", dataset_name, "good_selection")
    get_object_eff(results, "muon", dataset_name, "single_veto_selection")
    get_object_eff(results, "muon", dataset_name, "double_veto_selection")
    get_object_eff(results, "electron", dataset_name, "good_selection")
    get_object_eff(results, "electron", dataset_name, "single_veto_selection")
    get_object_eff(results, "electron", dataset_name, "double_veto_selection")
    get_object_eff(results, "tau", dataset_name)
    get_object_eff(results, "jet", dataset_name)
    
    logger.info(f"---> ################### Inspecting pair selections ################### <---")
    
    # pairs
    #print(f"\n---> Before trigobj matching <---")
    get_object_eff(results, "etau", dataset_name)
    get_object_eff(results, "mutau", dataset_name)
    get_object_eff(results, "tautau", dataset_name)
    
    #print(f"\n---> After trigobj matching <---")
    #pairs = []
    #pairs.append(["etau", ak.sum(ak.num(results.aux["etau"]["pairs"], axis=1) == 2)])
    #pairs.append(["mutau", ak.sum(ak.num(results.aux["mutau"]["pairs"], axis=1) == 2)])
    #pairs.append(["tautau", ak.sum(ak.num(results.aux["tautau"]["pairs"], axis=1) == 2)])
    ##pairs.append(["etau", ak.sum(ak.num(results.aux["etau"]["pair_indices"], axis=1) == 2)])
    ##pairs.append(["mutau", ak.sum(ak.num(results.aux["mutau"]["pair_indices"], axis=1) == 2)])
    ##pairs.append(["tautau", ak.sum(ak.num(results.aux["tautau"]["pair_indices"], axis=1) == 2)])
    #pair_table = tabulate(pairs, ["pair", "nEvents with pair"], tablefmt="pretty")
    
    #print(pair_table)

    logger.info(f"---> ################### Categorization ################### <---")
    cats = []
    cats.append(["is_etau", ak.sum(results.aux["cat_is_etau"])])
    cats.append(["is_mutau", ak.sum(results.aux["cat_is_mutau"])])
    cats.append(["is_tautau", ak.sum(results.aux["cat_is_tautau"])])
    cats.append(["is_etau_mutau", ak.sum(results.aux["cat_is_etau_mutau"])])
    cats.append(["is_etau_tautau", ak.sum(results.aux["cat_is_etau_tautau"])])
    cats.append(["is_mutau_tautau", ak.sum(results.aux["cat_is_mutau_tautau"])])
    cats_table = tabulate(cats, ["category", "nEvents"], tablefmt="pretty")
    
    logger.info(cats_table)
    
    logger.info(f"---> Events selected per channel <---")
    events = events[final_mask]
    sel_ev = ak.sum(events.channel_id > 0)
    logger.info(f"nSelectedEvents : {sel_ev}")
    channels = []
    etau_ev = ak.sum(events.channel_id == 1)
    mtau_ev = ak.sum(events.channel_id == 2)
    ttau_ev = ak.sum(events.channel_id == 4)
    mixed   = ak.sum(~((events.channel_id == 0)
                       | (events.channel_id == 1)
                       | (events.channel_id == 2)
                       | (events.channel_id == 4)))
    channels.append(["etau", etau_ev, round(etau_ev/sel_ev, 3)])
    channels.append(["mutau", mtau_ev, round(mtau_ev/sel_ev, 3)])
    channels.append(["tautau", ttau_ev, round(ttau_ev/sel_ev, 3)])
    channels.append(["other", mixed, round(mixed/sel_ev, 3)])
    channel_table = tabulate(channels, ["channel", "nEvents", "eff"], tablefmt="pretty")
    
    logger.info(channel_table)
    logger.info(f" ---> Total selected events in etau, mutau and tautau chennels : {etau_ev+mtau_ev+ttau_ev}\n\n")
    #from IPython import embed; embed()

    

def debug_extra_lepton_veto(nevts, *args):
    for i in range(nevts):
        if args[0].channel_id[i] < 1: continue
        logger.info(f"event : {args[0].event[i]}")
        logger.info(f"hcand_pairs_pt         : {args[1].pt[i]}")
        logger.info(f"extra leps pt          : {args[2].pt[i]}")
        logger.info(f"h1 & h2 pt             : {args[3].pt[i]}, {args[4].pt[i]}")
        logger.info(f"dr_hlep1_extraleps     : {args[5][i]}")
        logger.info(f"dr_hlep2_extraleps     : {args[6][i]}")
        logger.info(f"dr_mask                : {args[7][i]}")
        logger.info(f"has_extra_lepton       : {args[8][i]}")
        logger.info(f"has_no_extra_lepton    : {args[9][i]}\n")
        

def debug_double_lepton_veto(nevts, *args):
    for i in range(nevts):
        logger.info(f"event              : {args[0].event[i]}")
        logger.info(f"dl_veto_mu_pair_pt : {args[1].pt[i]}, {args[2].pt[i]}")
        logger.info(f"is Z like pair ?   : {args[3][i]}")
        logger.info(f"dl_veto_el_pair_pt : {args[4].pt[i]}, {args[5].pt[i]}")
        logger.info(f"is Z like pair ?   : {args[6][i]}")
        logger.info(f"concat masks       : {args[7][i]}")
        logger.info(f" --->> Any True means Z like pair exists --->> ")
        logger.info(f"has no Z like pair : {args[8][i]}\n")
        
