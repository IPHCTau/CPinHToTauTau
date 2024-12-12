A Columnflow based analysis framework from IPHC and DESY

<!-- marker-before-logo -->

<div style="text-align: center;">
    <img src="assets/logo.png" alt="Logo" style="width: 400px; height: 220px; display: block; margin: 0 auto;">
</div>

<!-- marker-after-logo -->

### Resources

- [columnflow](https://github.com/columnflow/columnflow/tree/master)
- [law](https://github.com/riga/law)
- [order](https://github.com/riga/order)
- [luigi](https://github.com/spotify/luigi)

### Setting up the framework (on the main branch)
Clone from the main github repo:
```
git clone --recurse-submodules git@github.com:DesyTau/CPinHToTauTau.git
```
ATTENTION: If you are setting up the framework on lxplus (i.e. on you /afs/cern.ch/user/user_first_char/cern_user_name/ you will encounter problems with micromaba.
A work around to this is selecting the /eos/user/user_first_char/cern_user_name/software as the path for the vens, cmssw and conda folders. As far as I know only afs is giving issues, columnflow people are aware of this but they can not do anything (Jacopo 11/12/24)
This can be achieved by specifying the path that you want to use in the iterative setting up procedure 

With this is mind, now you can set up the framework by:
```
cd CPinHToTauTau
source setup.sh your_setup_name
```
The first time this can take ~10min, but then you are good to go! 
### Setting up the framework used at DESY
Before the previous last step, there are few extra things to do.

But lets do the things step by step starting for cloning the repository:
REMEMBER: if you are planning to develop the code, forking the repo and then cloninf the fork is mandatory!
```
git clone --recurse-submodules git@github.com:DesyTau/CPinHToTauTau.git
cd CPinHToTauTau
```
Then you need to use the following setup:
```
https://github.com/DesyTau/CPinHToTauTau/blob/main/.gitmodules
```
This mean that essentially these 3 additional steps need to be performed:
1) CPinHToTauTau: https://github.com/DesyTau/CPinHToTauTau/tree/desy_dev; branch desy_dev
```
git checkout origin/desy_dev
```
This is enough to get the core of the analysis scripts and to have the modules/jsonpog-integration repository.
This is necessary for any kind of correction based on file.json

2) cmsdb: https://github.com/DesyTau/cmsdb/tree/desy_dev; branch desy_dev 
```
cd modules/cmsdb
```
```
git chechout origin/desy_dev
```
This is necessary for the definition of campaigns, processes and datasets.

3) columnflow: https://github.com/DesyTau/columnflow; branch production 
```
cd ../columnflow/
```
```
git chechout origin/production
```
This is necessary to have a stable version of columnflow that it is known to work with the current analysis set up.

Now you are finally ready to set up the framework!
```
cd ~
cd CPinHToTauTau
```
ATTENTION: If you are setting up the framework on lxplus (or on any other afs path) see above.
```
source setup.sh your_set_up_name
```
Enjoy the potential of columnflow and do not hesitate to reach us for any issue you may encounter!


### Installation

```sh
git clone --recursive https://github.com/gsaha009/CPinHToTauTau.git
cd CPinHToTauTau

# set environment
source setup.sh hcp
```

Then follow the details:

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
So, `law.cfg` needs to be set like here in this repo.

For any computing system, please make sure that this hard-coded (stupid) line is changed to the proper one:

`thisdir` in the config e.g. here: https://github.com/gsaha009/CPinHToTauTau/blob/main/httcp/config/config_run3.py#L33

### Best Practices

 - If you want to check whether your changes are working properly, use the config with the "_limited" tag in the end of your actual config. It will take one ROOT file as input to make the execution faster.
 - For production, if you need to send jobs to any batch system configured in `law`, it is better to run `cd.ReduceEventsWrapper` first. Then it is safe to use `cf.ProduceColumns`. By the end of this task, all possible event and object level corrections have already been applied. So, now usually there are three targets: Plotting, Saving the skimmed and corrected events array (porbably in flat ROOT ntuple) or to produce DataCard. In between, ML training or evaluation task can be used if needed. It is better to follow this `columnflow` [task-graph](https://github.com/columnflow/columnflow/wiki#default-task-graph), always.
 - Now, you may have a ton of datasets and processes. This [script](https://github.com/gsaha009/CPinHToTauTau/blob/main/cf_run.py) might help you to get the commands, which you can copy and paste in a `tmux/screen` session.
   - To start with, you can modify this [`yaml`](https://github.com/gsaha009/CPinHToTauTau/blob/main/yamls/2022PreEE_full.yml) depending on your need. It is recommended to make separate `yaml` for different eras.
   - Then you can run the `cf_run` script like `python cf_run.py -i <yaml/2022PreEE_full.yml or other> -f <ReduceEvents or other>`
   - Notice the version name. If it is dummy, the script will produce some version on it's own, otherwise of you want to use some older version, just specify the name there.
   - And, at the end, you will get the command you want to run.
 - While using `cf.PlotVariables1/2D`, better not to mention all the categories in the command, because it will take ages as it creates the categories in the runtime. So, it would be better to produce plots in several steps. 