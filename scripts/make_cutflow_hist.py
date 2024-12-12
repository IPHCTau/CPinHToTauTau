import matplotlib.pyplot as plt
import mplhep
import pickle
import numpy as np


def create_cutflow_histogram(cuts, data, xlabel="Selections", ylabel="Selection efficiency", title="", log=False, rel=False, save_path=None):
    """
    Create a cutflow histogram using Matplotlib with CMS style and save it to a PDF file.

    Parameters:
    - cuts: List of strings representing the names of the cuts.
    - data: List of integers representing the corresponding event counts for each cut.
    - xlabel: Label for the x-axis (default is "Cuts").
    - ylabel: Label for the y-axis (default is "Events").
    - title: Title of the plot (default is "Cutflow Histogram").
    - save_path: Path to save the PDF file. If None, the plot will be displayed but not saved.

    Returns:
    - fig: Matplotlib figure object.
    - ax: Matplotlib axis object.
    """

    # Set CMS style
    plt.style.use(mplhep.style.CMS)
   
    # Create Matplotlib figure and axis
    fig, ax = plt.subplots()
    if log: plt.yscale('log')
    # Plotting the cutflow histogram
    color = ['black','red']
    
    for i, (name, n_evt) in enumerate(data.items()):
        
        n_evt = np.array(n_evt)
        for cut_name, the_n_evt in zip(cuts,n_evt):
            print(f"{cut_name}: {the_n_evt}")
        if rel:
            x = cuts[1:]
            y = n_evt[1:]/n_evt[:-1]
        else:
            x = cuts
            y = n_evt/n_evt[0]
        print(f'Event numbers:')
        ax.scatter(x, y , color=color[i], marker='o', alpha=0.8, label=f"Data {name}")
        n_evt_prev = n_evt

    if log: ax.set_ylim((1*10**-2,2*10**0))
    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticklabels(cuts[1:] if rel else cuts, rotation=45, ha='right')
    
    # Add legend
    ax.legend()
    ax.grid(True, which='both', axis='y')
    ax.grid(True, which='major', axis='x')
    label_options = {
        "wip": "Work in progress",
        "pre": "Preliminary",
        "pw": "Private work",
        "sim": "Simulation",
        "simwip": "Simulation work in progress",
        "simpre": "Simulation preliminary",
        "simpw": "Simulation private work",
        "od": "OpenData",
        "odwip": "OpenData work in progress",
        "odpw": "OpenData private work",
        "public": "",
    }
    cms_label_kwargs = {
        "ax": ax,
        "llabel": label_options.get("pw"),
        "fontsize": 22,
        "data": True,
        'rlabel': "Data to Data"
    }
    mplhep.cms.label(**cms_label_kwargs)
   
    plt.tight_layout()
    
     # Save to PDF if save_path is provided
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    else:
        # Show the plot if save_path is not provided
        plt.show()


    return fig, ax

def get_hist_values(pickle_file):
    file_obj = open(pickle_file, 'rb')
    data = pickle.load(file_obj)
    hist = data.profile(axis=0)
    cuts = []
    values = []
    print(hist)
    for cut_name in hist.axes[1]:
        cuts.append(f'{cut_name}')
        values.append(hist[0,f'{cut_name}',0,0].count)
    return cuts, values

 



#path22 = "/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.CreateCutflowHistograms/run3_2022_postEE_nano_tau_v12/data_mu_g/nominal/calib__example/sel__default__steps_trigger_muon_pt_26_muon_eta_2p4_mediumID_muon_dxy_0p045_muon_dz_0p2_muon_iso_0p15_DeepTauVSjet_DeepTauVSe_DeepTauVSmu_tau_eta_2p3_tau_dz_0p2_tau_pt_20_single_pair_extra_lep_veto_dilep_veto/condor_production/cutflow_hist__event.pickle"

path22 = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/data/cf_store/analysis_httcp/cf.CreateCutflowHistograms/run3_2022_preEE_nano_cp_tau_v12_limited/dy_lep_m50/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/v1/cutflow_hist__event.pickle"

#path18 = "/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.CreateCutflowHistograms/run2_UL2018_nano_tau_v10_limited/data_ul2018_a_single_mu/nominal/calib__example/sel__default__steps_trigger_muon_pt_26_muon_eta_2p4_mediumID_muon_dxy_0p045_muon_dz_0p2_muon_iso_0p15_DeepTauVSjet_DeepTauVSe_DeepTauVSmu_tau_eta_2p3_tau_dz_0p2_tau_pt_20_mutau_os_mutau_dr_0p5_mutau_mt_50_single_pair/dev1/cutflow_hist__event.pickle"
#path22 = "/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.CreateCutflowHistograms/run3_2022_postEE_nano_tau_v12_limited/data_mu_f/nominal/calib__example/sel__default__steps_trigger_muon_pt_26_muon_eta_2p4_mediumID_muon_dxy_0p045_muon_dz_0p2_muon_iso_0p15_DeepTauVSjet_DeepTauVSe_DeepTauVSmu_tau_eta_2p3_tau_dz_0p2_tau_pt_20_single_pair_extra_lep_veto_dilep_veto/multi_pair/cutflow_hist__event.pickle"
cuts, values22 = get_hist_values(path22)
#cuts, values18 = get_hist_values(path18)

create_cutflow_histogram(cuts, 
                         data={"2022_preEE": values22},
                         save_path="cutflow_histogram_l.pdf",
                         log=False,
                         rel=False)
