# CPinHToTauTau Analysis

A Columnflow based analysis framework from IPHC, forked from [desyTau](https://github.com/DesyTau/CPinHToTauTau). Although the basic workflow is similar to the desyTau CP analysis setup, but some differences are there mainly in array structures, categorizations and implementations of some methods. However, synchronization study is being performed recently.

### Link to the AN : [AN-25-008](https://gitlab.cern.ch/tdr/notes/AN-25-008)

<!-- marker-before-logo -->

<div style="text-align: center;">
    <img src="assets/logo.png" alt="Logo" style="width: 400px; height: 220px; display: block; margin: 0 auto;">
</div>

<!-- marker-after-logo -->

### Resources

- [luigi](https://github.com/spotify/luigi)
- [law](https://github.com/riga/law)
- [columnflow](https://github.com/IPHCTau/columnflow)
- [order](https://github.com/riga/order)
- [cmsdb](https://github.com/IPHCTau/cmsdb)

### Installation

```sh
git clone --recursive https://github.com/IPHCTau/CPinHToTauTau.git
cd CPinHToTauTau

# set environment
source setup.sh hcp
```

Then follow the details:

<span style="color:red">In LxPlus, in order to use HTCondor, make sure of the location of `CF_DATA` and `CF_JOB_BASE` must be in `afs` as HTCOndor does not support job submission from `eos` yet. </span>
e.g.

```sh
CERN username (CF_CERN_USER, default gsaha):
Local data directory (CF_DATA, default ./data):
Local directory for installing software (CF_SOFTWARE_BASE, default $CF_DATA/software):  /eos/user/g/gsaha/CPinHToTauTauData
Local directory for storing job files (CF_JOB_BASE, default $CF_DATA/jobs):             
Relative path used in store paths (see next queries) (CF_STORE_NAME, default cf_store):  
Default local output store (CF_STORE_LOCAL, default $CF_DATA/$CF_STORE_NAME):  
Local directory for caching remote files (CF_WLCG_CACHE_ROOT, default ''):  
Automatically update virtual envs if needed (CF_VENV_SETUP_MODE_UPDATE, default False):  
Use a local scheduler for law tasks (CF_LOCAL_SCHEDULER, default True):  
Flavor of the columnflow setup ('', 'cms') (CF_FLAVOR, default cms):  
storage element for crab specific job outputs (e.g. T2_DE_DESY) (CF_CRAB_STORAGE_ELEMENT, default ''):  
base directory on storage element for crab specific job outputs (CF_CRAB_BASE_DIRECTORY, default /store/user/$CF_CERN_USER/cf_crab_outputs):
```

This will install the `softwares` in `eos`, but the `data` and `jobs` directories are in `afs`.
In the next step, the big output files have to be stored in `eos`.
So, `law.cfg` needs to be set as used in this repo. 

For any computing system, please make sure that this hard-coded (stupid) line is changed to the proper one:

`thisdir` in the main config e.g. [here](https://github.com/IPHCTau/CPinHToTauTau/blob/main/httcp/config/config_run3.py#L35)

### Best Practices :+1:

To ensure a smooth test of the workflow, here are some of the "best(?)" practices:

- **Testing Changes:**
  Use a configuration file (e.g. `httcp.config.config_run3.py`) with the `_limited` suffix to the campaign name to speed up testing. It will only process one `ROOT` file for faster execution. One has to modify the `campaign_dict` inside `httcp.config.analysis_httcp.py` before using the `limited` or `full` campaign. 

- **For Production Runs:**  
  Use `cf.ReduceEventsWrapper` to pre-process your events before running `cf.ProduceColumns`. This ensures that all object and event level selections are done, and the workflow is ready for production. The typical process involves:
  - **Plotting**
  - **Saving skimmed/corrected events (as parquet files or flat ROOT ntuples)**
  - **Creating DataCards**
  - If necessary, you can include **ML training** or **evaluation tasks**.  
  Follow the default task graph in ColumnFlow [task-graph](https://github.com/columnflow/columnflow/wiki#default-task-graph) is recommended.

- **Handling Large Datasets:**  
  If you have multiple datasets or processes, this [script](https://github.com/IPHCTau/CPinHToTauTau/blob/main/cf_run.py) can help generate the necessary commands for batch execution.  
  Steps:
  1. Modify the `yaml` (e.g. [2022PreEE.yml](https://github.com/IPHCTau/CPinHToTauTau/blob/main/yamls/2022PreEE.yml)) based on your requirements.
  2. Run the `cf_run` script:  
     ```sh
     python cf_run.py -i <yaml/2022PreEE_full.yml> -f <ReduceEvents or other> -l <if you want to use limited config>
     ```
     The script will generate the necessary commands for execution.

- **Plotting Tips:**  
  When using `cf.PlotVariables1/2D`, avoid specifying all categories in one command as it can take a long time to generate. Break the plotting into several smaller steps for faster execution.


### Useful links :paperclip: 
These are useful links that can save you time and effort.  Because let's face it, even though you know these websites exist, you always manage to misplace the link.

- [NanoAOD doc](https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html#GenPart)

### Tips and trick :hear_no_evil:
Welcome to the Tips and Tricks section! Here, you'll find recipes to make your journey with ColumnFlow smoother and more efficient. Let's unlock the secrets to mastering ColumnFlow with ease! :rocket:


#### How to use ROOT in columnflow :deciduous_tree:
Before each session and before configuring your ColumnFlow setup, ensure that you import the necessary library. For ROOT, if you're utilizing lxplus, you can execute the following command:

```
CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh 
```

If you want to use ROOT in one of your script, at the beginning of your file, import the ROOT module as `import ROOT`.

#### Plot quickly ak.array :bar_chart:	
Sometimess it can be nice to quickly have a look a the distribution of your favorite variable. For this you can include the following modules: 
```
mpl = maybe_import("matplotlib")
plt = maybe_import("matplotlib.pyplot")
```

And save your distribution in few lines:

```
# Make the plot
fig, ax = plt.subplots()
plt.hist(ak.to_numpy(leps1["mass"]).flatten(), bins=20, range=(0, 2), alpha = 0.5)
plt.savefig("mass.pdf") 
```