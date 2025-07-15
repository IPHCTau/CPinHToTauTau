# Fake Factor Measurement Project

## Overview

This project provides a set of scripts for deriving **fake factors (FF)**, **closure corrections (CR)**, and **extrapolation corrections (ER)**. The output includes fitted correction factors stored in **ROOT** files and **CorrectionLib-compatible JSON files**, which can be used in high-energy physics analyses.

## Features

- 🏗 **Fake factor derivation** for systematic studies
- 🛠 **Closure and extrapolation corrections** for improved accuracy
- 📈 **Fitting functions** (Landau, polynomial, etc.) for corrections
- 📁 **JSON output generation** for compatibility with CorrectionLib
- ⚙️ **Configurable and modular pipeline**

## Project Structure

```
fake_factor_derivation/
│
├── outputs/                    # Store the output JSON and ROOT files
│   ├── 2022_preEE              # Dataset before the End-of-Year Electron Run (pre-EE)
│   │   ├── fake_factors/       # Fake factor calculations
│   │   │   ├── json/fake_factor_2022_preEE.json   # Combined FF JSON file
│   │   │   ├── pi_1/fake_factor_pi_1_has0j_2022_preEE.json  # Specific JSON files for decay modes (DMs)
│   │   │   ├── rho_1/...
│   │   ├── closure_corrections/   # Same structure as fake_factors/
│   │   ├── extrapolation_corrections/   # Same structure as fake_factors/
│   ├── 2022_postEE             # Dataset after the End-of-Year Electron Run (post-EE), same structure as preEE
│
├── src/                         # Source code directory
│   ├── fit_and_save_ratios.py        # Script to plot ratios and fit FFs with Landau/polynomial functions
│   ├── config.json              # Configuration file specifying options (categories, corrections, input paths, etc.)
│   ├── save_correctionlib.py    # Script to convert fit results into a CorrectionLib-compatible JSON file
│   ├── fit_functions.py         # Defines fitting functions for FF, CR, ER, etc.
│   ├── combine_jsonFF.py        # Merges JSON files from different DMs and NJets into a single JSON file
│
├── README.md                    # Documentation file
└── requirements.txt              # Required python3 dependencies
```

## Installation

Set usaul CF envrionment:
```sh
source setup httcp_env
```

Ensure you have **ROOT** installed for handling histograms and fits. If not run:
```sh
CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh 
```

## Usage

### 1. <span style="color:green">Run CF</span> to produce the FF ratio in a pkl file:
Produce ratio of data_minus_mc_catA/data_minus_mc_catB = FF and data_minus_mc_catA0/data_minus_mc_catB0 = FF0. The FF and FF0 are pt dependent and can be produce by DM and Njet or inclusively. This is done with the following line:
```sh
law run Plot1dVariable ... 
```
And save the output path of the ratio .pkl file.
Use this path in the config.json, add the variable (pt), the dm, njet (-1 if inclusive), year, type of correction. See example `config_FF_2022_preEE.json`.

### 2. Fit and produce json of the Fake Factors
```sh
python3 src/fit_and_save_ratios.py --config config_FF_2022_preEE.json
```
This will generate FF histograms and save fitted functions in a ROOT file and json format for each category. For example `fake_factors_2022_preEE.json`.

You can combine dm and njet multiplicity files with: 

```sh
python3 src/combine_jsonFF.py --config config_FF_2022_preEE.json
```


### 3. <span style="color:green">Run CF</span> to produce Closure plots, control plots:
Now you can plug the FF file `fake_factors_2022_preEE.json` in columnflow. Change the tag `cfg.x.regions_to_extrapolate_fake = "AB"` in `httcp/config/config_run3.py`. You can produce closure plots by applying FFxB = A, like:
```sh
law run Plot1dVariable ... 
```
Or change the tag `cfg.x.regions_to_extrapolate_fake = "CD"` in `httcp/config/config_run3.py` and you can produce control plots by applying FFxC = D, like: 
```sh
law run Plot1dVariable ... 
```

Note: here you run the same version as step 1, but you have to run from cf.ProduceColumn

### 4. If needed, fit and produce json of the closure correction
If required, closure correction can be produced like fake foctors with the ratio pkl files. This time the ratio is data/MC of A category. Similarly to step 2 you can run:
```sh
python3 src/fit_and_save_ratios.py --config config_CC_2022_preEE.json
```
This time you need to modify the config json to specify the "closure _correction" as correction type and the varaibles used.

### 5. If needed, <span style="color:green">Run CF</span> to produce control plots and new A0/BO = FF0 with closure correction applied


### 6. <span style="color:green">Run CF</span> to produce extrapolation plots FF0xC0 = D0
```sh
law run Plot1dVariable ... 
```
And save the ratio pkl .file.

### 6. If needed, fit and produce json of the extrapolation corrections
If required, extrapolation correction can be produced like fake factors with the ratio pkl files. This time the ratio is data/MC of D0 category. Similarly to step 2 you can run:
```sh
python3 src/fit_and_save_ratios.py --config config_EC_2022_preEE.json
```
This time you need to modify the config json to specify the "extrapolation_correction" as correction type and the variables used.

### 7. <span style="color:green">Run CF</span> to produce final control plots. 
FF, CR and EC should be apply, and you can run:
```sh
python3 src/fit_and_save_ratios.py --config config_CC_2022_preEE.json
```

## Configuration
The `config.json` file specifies options such as:
```json
{
    "era": "2022_postEE",
    "correction_type": "fake_factors",
    "variable": "hcand_1_pt",
    "input_base_path": "/eos/project/i/iphctau/public/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.PlotVariables1D/run3_2022_postEE_nano_cp_tau_v14/FF",
    "hist_name": "ratio",
    "dm":["tau1a1DM10", "tau1a1DM11", "tau1a1DM2", "tau1pi", "tau1rho"],
    "n_jets": ["has0j", "has1j", "has2j"]
}
```

## Output Format
The final CorrectionLib JSON format follows:
```json
{
    "schema_version": 2,
    "corrections": [
        {
            "name": "fake_factor",
            "inputs": [{ "name": "pt", "type": "real" }],
            "output": { "name": "ff", "type": "real" },
            "data": {
                "nodetype": "formula",
                "expression": "a + b*pt"
            }
        }
    ]
}
```

