# LenPrEP

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Zenodo](https://img.shields.io/badge/Zonodo-doi:10.5281/zenodo.14154555-green)](https://doi.org/10.5281/zenodo.16612336) 


This repository provides a framework for pharmacokinetic (PK) and pharmacodynamic (PD) model fitting of Lenacapavir using clinical data for various dosing routes (oral, subcutaneous, intramuscular). It supports both local and global optimizers via [lmfit](https://lmfit.github.io/lmfit-py/), customizable parameter initialization via .ini files, and automatic plotting and result export.

## Table of Contents

- [System requirements](#system-requirements)
    - [Dependencies](#dependencies)
- [Folder Structure](#folder-structure)
- [Clinical Data](#clinical-data)
- [Installation](#installation)
- [Code](#code)
  - [Running from Terminal](#running-from-terminal)
  - [Arguments](#arguments)
  - [Show CLI Help](#show-cli-help)

## System requirements

The code was tested on Ubuntu 24.04.2 LTS (Noble Numbat) with Linux kernel version 6.14.0-24-generic.

### Dependencies

| Package     | Version  |
|-------------|----------|
| Python      | ≥ 3.12   |
| numpy       | 2.1.3    |
| pandas      | 2.2.3    |
| scipy       | 1.14.1   |
| matplotlib  | 3.9.2    |
| lmfit       | 1.3.2    |

They are included in `requirements.txt` and can be installed as explained in [Installation](#installation).

## Folder Structure

```
project-root/
├── Code/
|   ├── requirements.txt        # python dependencies
│   ├── cli.py                  # CLI argument parser
│   ├── run_model.py            # Execute fitting
│   ├── pk_fitting.py           # PK model fitting
│   ├── pd_fitting.py           # PD model fitting
│   ├── utils.py                # Shared utilities (ODEs, solving, plotting)
│   └── ini/                    # Parameter .ini files per mode
│       ├── PK_PO.ini
│       ├── PK_IM.ini
│       ├── PK_SC_PEG.ini
│       ├── PK_SC_AQ1.ini
│       ├── PK_SC_AQ2.ini
│       └── PD_SC.ini
├── Data/
│   ├── PK/
│   │   ├── oral/               # PK oral clinical data (PO, Study 2)
│   │   ├── IM/                 # PK intramuscular clinical data (IM, Study 5)
│   │   └── SC/                 # PK subcutaneous formulations (SC)
│   │       ├── SC_PEG/         # PEG/water solution (Study 4)
│   │       ├── SC_AQ1/         # SC aqueous suspension (Study 1)
│   │       └── SC_AQ2/         # SC aqueous suspension (Study 3)
│   └── PD/                     # PD viral load clinical data (VL, Study 3)
|
├── Results/                    # Output plots and files
└── README.md                   # Documentation


```

## Clinical Data

| StudyID | Focus / Population                                   | Formulation             | Dataset(s)                                   |
| -------- | ---------------------------------------------------- | ----------------------- | -------------------------------------------- |
| Study 1  | Phase I, healthy volunteers                          | SC aqueous suspension   | 30, 100, 300, 450 mg                         |
| Study 2  | Phase I, healthy volunteers                          | Oral tablets            | 300, 900, 1,800 mg                           |
| Study 3  | Phase Ib, HIV-1 positive subjects (antiviral effect) | SC aqueous suspension   | 20, 50, 150, 450 mg                          |
| Study 4  | Phase I, healthy volunteers                          | SC PEG/water suspension | 309, 927 mg                                  |
| Study 5  | Phase I, healthy volunteers                          | IM ethanol/water        | 5,000 mg with 5% (F1) and 10% ethanol (F2)   |


## Installation

```bash
git clone https://github.com/KleistLab/LenPrEP.git
cd LenPrEP/Code            # go to folder
python3 -m venv venv       # create isolated environment
source venv/bin/activate   # activate it
pip install -r requirements.txt  # install dependencies
```

## Code

### Running from Terminal

```bash
# Minimal example: run default fitting for a specific mode
python run_model.py --mode PK_SC_PEG

# Run fitting and generate plot comparing model to clinical data
python run_model.py --mode PK_IM --plot True

# Use a global optimizer (automatically sets --max_count=1)
python run_model.py --mode PD_SC --method differential_evolution

# Track optimizer progress per iteration (only for local methods)
python run_model.py --mode PK_PO --track_iterations True

# Specify a different output directory
python run_model.py --mode PK_SC_AQ1 --output /user/output

```

### Arguments

| Argument (Flags)        | Type | Default         | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
| ----------------------- | ---- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--mode`             | str  | **(required)**  | Fitting mode: choose one of the following combinations:<br>• `PK_PO` – PK, oral dosing (Study 2)<br>• `PK_SC_PEG` – PK, SC PEG formulation (Study 4)<br>• `PK_SC_AQ1` – PK, SC aqueous formulation 1 (Study 1)<br>• `PK_SC_AQ2` – PK, SC aqueous formulation 2 (Sudy 3)<br>• `PK_IM` – PK, intramuscular dosing (Study 5)<br>• `PD_SC` – pharmacodynamics, SC only (Study 3)                                                                             |
| `--method`           | str  | `least_squares` | Optimization method for `lmfit.minimize()`.<br><br>**Local methods** :<br>• `least_squares`, `nelder`, `lbfgsb`, `powell`, `cg`, `newton`, `cobyla`, `tnc`, `trust-constr`<br><br>**Global methods** :<br>• `differential_evolution`, `brute`, `basinhopping`, `ampgo`<br><br>Note: global methods automatically set `--max_count=1`. |
| `--max_count`        | int  | `10`            | Maximum number of fitting iterations. For local optimizers, multiple randomized initializations are used. Global optimizers force this to `1`.                                                                                                                                                                                                                                                               |
| `--plot`             | bool | `False`         | Set to `True` to plot best fit vs. clinical data after fitting.                                                                                                                                                                                                                                                                                                                                              |
| `--output`           | str  | `Results/`      | Folder to save output files (default: `Results/`). Relative to project root.                                                                                                                                                                                                                                                                                                                                 |
| `--track_iterations` | bool | `False`         | If `True`, prints RSS and parameter values at every optimizer iteration (only applicable to local methods like `least_squares`).                                                                                                                                                                                                                                                                             |

### Show CLI Help

To see all available options and descriptions directly from the command line:

```bash
python run_model.py --help

usage: run_model.py [-h] --mode {PK_PO,PK_SC_PEG,PK_SC_AQ1,PK_SC_AQ2,PK_IM,PD_SC}
                    [--method {least_squares,nelder,lbfgsb,powell,cg,newton,cobyla,tnc,trust-constr,differential_evolution,brute,basinhopping,ampgo}]
                    [--max_count MAX_COUNT] [--plot True|False] [--output OUTPUT]
                    [--track_iterations True|False]

Run PK/PD model fitting based on selected mode. Available modes combine the model type (PK or PD) with the dosing route.

options:
  -h, --help            Show this help message and exit
  --mode                Fitting mode: choose one of the following dosings:
                          PK_PO       - Oral dosing (Study 2)
                          PK_SC_PEG   - SC PEG formulation (Study 4)
                          PK_SC_AQ1   - SC aqueous formulation 1 (Study 1)
                          PK_SC_AQ2   - SC aqueous formulation 2 (Study 3)
                          PK_IM       - Intramuscular dosing (Study 5)
                          PD_SC       - Subcutaneous dosing only (Study 3)
  --method              Optimization method for lmfit.minimize(). Default: least_squares.
                          Local methods (better for large parameter ranges):
                            least_squares, nelder, lbfgsb, powell, cg, newton, cobyla, tnc, trust-constr
                          Global methods (better for smaller parameter ranges):
                            differential_evolution, brute, basinhopping, ampgo
                          Note: global methods automatically set --max_count=1.
  --max_count           Maximum number of fitting iterations (default: 10)
  --plot                Set to True to plot best fit vs clinical data after fitting (default: False)
  --output              Folder to save output files (default: Results/)
  --track_iterations    If True, prints RSS and parameter values at every optimizer iteration (only for local methods)

```

