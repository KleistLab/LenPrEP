import pandas as pd
from pathlib import Path
import configparser
from lmfit import Parameters
from scipy.integrate import solve_ivp
from types import SimpleNamespace
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import math

def get_doses(path):
    """
    Extract base dose values from filenames in the given directory.

    Args:
        path (str): Path to the directory containing dose CSV files

    Returns:
        tuple:
            - doses (list of int): Unique sorted list of dose values (e.g., [5000])
            - dose_counts (dict): Number of datasets per dose (e.g., {5000: 2})
    """
    pattern = re.compile(r'^(\d+)mg(?:_(\d))?\.csv$', re.IGNORECASE)
    dose_counts = {}

    for f in os.listdir(path):
        m = pattern.match(f)
        if m:
            dose = int(m.group(1))
            dose_counts[dose] = dose_counts.get(dose, 0) + 1

    doses = sorted(dose_counts.keys())
    return doses, dose_counts


def read_data(mode, path, doses, dose_counts):
    """
    Reads PK or PD data files from 'path' for given 'doses'.

    Args:
        mode (str): 'PK_PO', 'PK_SC', 'PK_IM', 'PD_SC', etc.
        path (str or Path): directory containing the CSV files
        doses (list of int): list of base doses (e.g. [dose_1, dose_2, ...])
        dose_counts (dict): how many CSV files per dose

    Returns:
        data (list of dict): [{'dose': 5000, 'time': [...], 'conc': [...]}, ...]
                             One dict per file, even if same dose
    """
    data = []
    pattern = re.compile(rf'^(\d+)mg(?:_(\d))?\.csv$', re.IGNORECASE)

    for dose in doses:
        expected_count = dose_counts.get(dose, 1)
        matched = []

        for f in sorted(os.listdir(path)):
            m = pattern.match(f)
            if m and int(m.group(1)) == dose:
                matched.append(f)

        if len(matched) != expected_count:
            raise FileNotFoundError(f"Expected {expected_count} files for {dose}mg but found {len(matched)}")

        for filename in matched:
            filepath = os.path.join(path, filename)
            df = pd.read_csv(filepath).dropna()

            if mode.startswith('PK'):
                df.columns = ['time', 'conc']
                df['conc'] = np.log10(df['conc'])
                entry = {
                    'dose': dose,
                    'time': df['time'].tolist(),
                    'conc': df['conc'].tolist()
                }

            elif mode.startswith('PD'):
                df.columns = ['time', 'change']
                entry = {
                    'dose': dose,
                    'time': df['time'].tolist(),
                    'change': df['change'].tolist()
                }

            else:
                raise ValueError(f"Unsupported mode: {mode}")

            data.append(entry)

    return data




def create_initial_conditions(mode, doses, dose_counts=None):
    """
    Create initial condition vectors z0 for a given mode and list of  (PK only).

    Args:
        mode (str): Mode string, e.g., 'PK_PO', 'PK_SC_PEG', 'PK_IM'
        doses (list of int): List of base doses
        dose_counts (dict, optional): Number of datasets per dose. If None, assume 1 per dose

    Returns:
        list of list: Initial condition vectors for each dataset
    """
    z0 = []

    # Fallback: 1 file per dose if dose_counts not provided
    if dose_counts is None:
        dose_counts = {d: 1 for d in doses}

    for dose in doses:
        for _ in range(dose_counts[dose]):
            if mode == 'PK_PO':
                z0.append([dose, 0, 0])
            elif mode.startswith('PK_SC') or mode == 'PK_IM':
                z0.append([0, dose, 0])
            else:
                raise ValueError(f"Unsupported mode for z0 initialization: {mode}")

    return z0


def read_ini_file(path):
    """
    Read an .ini file containing parameter definitions and convert to dict.

    Args:
        path (str or Path): Path to the INI file

    Returns:
        dict: Parsed parameter dictionary for model initialization
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(path)

    params = {}

    for section in config.sections():
        section_dict = {}
        for key, val in config[section].items():
            try:
                parsed_val = float(val)
            except ValueError:
                if val.lower() == 'true':
                    parsed_val = True
                elif val.lower() == 'false':
                    parsed_val = False
                else:
                    parsed_val = val
            section_dict[key] = parsed_val

        if 'value' in section_dict:
            params[section] = {
                'value': float(section_dict['value']),
                'min': float(section_dict.get('min', '-inf')),
                'max': float(section_dict.get('max', 'inf')),
                'vary': bool(section_dict.get('vary', True)),
            }
        else:
            params[section] = section_dict

    return params



def initialize_from_ini(ini_params, doses, mode):
    """
    Initialize lmfit.Parameters object from ini_params dictionary.

    Args:
        ini_params (dict): Parsed ini file parameters with structure:
                           {'Section': {param: value, ...}, ...}
        doses (list): Used to expand 'Vd_default' into Vd_{dose}
        mode (str): e.g. 'PK_PO', 'PK_SC_AQ1', 'PD_SC', ...

    Returns:
        lmfit.Parameters object
    """
    params = Parameters()

    # Add all lmfit-style parameters (sections with value/min/max/vary)
    for section, entry in ini_params.items():
        if isinstance(entry, dict) and 'value' in entry:
            params.add(section, value=entry['value'],
                       min=entry.get('min', -np.inf),
                       max=entry.get('max', np.inf),
                       vary=entry.get('vary', True))

    # Add viral dynamics parameters (fixed)
    for k, v in ini_params.get('VD parameters', {}).items():
        if isinstance(v, str):
            v = float(v.split('#')[0].strip())
        params.add(k, value=v, vary=False)

    # Add PK model parameters (fixed)
    for k, v in ini_params.get('Model parameters', {}).items():
        if isinstance(v, str):
            v = float(v.split('#')[0].strip())
        elif not isinstance(v, (int, float)):
            continue  # skip non-numeric entries
        params.add(k, value=v, vary=False)

    # Dose-specific Vd values from Vd_default if present (for PK PO)
    if 'Vd_default' in ini_params and mode.startswith('PK'):
        vd_entry = ini_params['Vd_default']
        for dose in doses:
            params.add(f"Vd_{dose}", value=vd_entry['value'],
                       min=vd_entry.get('min', -np.inf),
                       max=vd_entry.get('max', np.inf),
                       vary=vd_entry.get('vary', True))

    return params


def randomize_params(params):
    """
    Randomize parameter values for fitting within their allowed range.

    Args:
        params (lmfit.Parameters): Parameter object to randomize

    Returns:
        lmfit.Parameters: Randomized parameter object
    """
    randomized = params.copy()

    for name, p in randomized.items():
        if not p.vary:
            continue  # skip fixed parameters

        center = p.value

        # Skip if center is zero or too small
        if center == 0 or np.isclose(center, 0):
            continue

        # One order of magnitude around the center
        lower = max(p.min, center / 10)
        upper = min(p.max, center * 10)

        if lower >= upper:
            continue  # invalid range

        randomized[name].set(value=np.random.uniform(lower, upper))

    return randomized

def per_iteration(pars, iteration, resid, *args, **kws):
    """
    Callback function for lmfit that prints parameter info of each iteration.

    Args:
        pars (lmfit.Parameters): Current parameter object
        iteration (int): Current iteration number
        resid (ndarray): Current residuals
    """
    print(f"\nIteration {iteration + 1}")
    print(f"RSS = {np.sum(np.square(resid)):.6f}")
    print("Parameters (vary=True):")
    for p in pars.values():
        if p.vary:
            print(f"  {p.name:10s} = {p.value:.5g}")


def ode_system(t, y, ka, ke, kd, kid, Vd, n, Frac):
    """
    ODE system for PK model.

    Args:
        t (float): Time
        y (list): State vector
        ka, ke, kd, kid, Vd, n, Frac: Model parameters

    Returns:
        list: Derivatives
    """
    Oral, SC, C = y

    dOral = -ka * Oral
    dSC = -(kd + kid) * SC

    dC = (
        ka / Vd * Oral +
        kd / Vd * SC * Frac +
        kid / Vd * SC * (1 - Frac) * ((kid * t) ** n / math.factorial(n)) * np.exp(-kid * t) -
        ke * C
    )

    return [dOral, dSC, dC]


def eta(D, IC50, m):
    """Compute drug effect using Emax Hill function."""
    return D ** m / (IC50 ** m + D ** m)

def NT_T(D, IC50, m, NT):
    return (1 - eta(D, IC50, m)) * NT

def NT_M(D, IC50, m, NM):
    return (1 - eta(D, IC50, m)) * NM


def get_model_params(params):
    """
    Convert lmfit.Parameters to SimpleNamespace for easier access.

    Args:
        params (lmfit.Parameters): Parameters object.

    Returns:
        SimpleNamespace: Attribute-style access to parameters.
    """
    d = {k: v.value if hasattr(v, 'value') else v for k, v in params.items()}
    return SimpleNamespace(**d)


def ode_complex(t, z, params, IC50):
    """
    Full ODE system incorporating PD and viral dynamics models.

    This system includes:
    - Uninfected and infected T-cells (Tu, T1, T2)
    - Uninfected and infected macrophages (Mu, M1, M2)
    - Virus populations (V, VN)
    - Drug compartments (SC, C)

    The drug effect is modeled as a concentration-dependent inhibition using an
    Emax/Hill-type function (`eta`) with IC50 and m.

    Args:
        t (float): Current time
        z (list of float): State vector containing all compartments
        params (lmfit.Parameters): Model parameters
        IC50 (float): 50% inhibitory concentration for the drug

    Returns:
        list of float: Derivatives for the full state vector at time t
    """
    p = get_model_params(params)
    m = int(p.m)
    n = int(p.n)

    Tu, T1, T2, V, Mu, M1, M2, VN, SC, C = z

    dTu = p.lambdaT - p.deltaT * Tu - p.betaT * V * Tu + p.deltaPICT * T1
    dMu = p.lambdaM - p.deltaM * Mu - p.betaM0 * V * Mu + p.deltaPICM * M1
    dT1 = p.betaT * V * Tu - (p.deltaT1 + p.kT + p.deltaPICT) * T1
    dM1 = p.betaM0 * V * Mu - (p.deltaM1 + p.kM + p.deltaPICM) * M1
    dT2 = p.kT * T1 - p.deltaT2 * T2
    dM2 = p.kM * M1 - p.deltaM2 * M2

    dV = NT_M(C, IC50, m, p.NM) * M2 + NT_T(C, IC50, m, p.NT) * T2 \
         - V * (p.CL + 2 * p.betaT * Tu + 2 * p.betaM0 * Mu)

    dVN = ((p.NThat - NT_T(C, IC50, m, p.NT)) * T2 +
           (p.NMhat - NT_M(C, IC50, m, p.NM)) * M2) - p.CL * VN

    dSC = -(p.kd + p.kid) * SC

    dC = (
        p.kd / p.Vd * SC * p.Frac +
        p.kid / p.Vd * SC * (1 - p.Frac) *
        ((p.kid * t) ** n / math.factorial(n)) * np.exp(-p.kid * t) -
        p.ke * C
    )

    return [dTu, dT1, dT2, dV, dMu, dM1, dM2, dVN, dSC, dC]


def solve_pk_ode(z0, times, params, doses, mode, conversion_factor_ng_ml=1000):
    """
    Solve the PK ODE system for each dataset.

    Args:
        z0 (list of list): Initial conditions
        times (list of list): Timepoints for each dataset
        params (lmfit.Parameters): Fitting parameters
        doses (list of int): Doses administered
        mode (str): Dosing route for PK data
        conversion_factor_ng_ml (float): Unit conversion factor (L -> mL)

    Returns:
        list of list: Log10 concentrations over time
    """
    solutions = []

    for i in range(len(z0)):
        ke = params['ke'].value

        if mode == 'PK_PO':
            ka = params['ka'].value
            kd, kid, n, Frac = 0, 0, 0, 0
            Vd = params[f'Vd_{doses[i]}'].value

        elif mode == 'PK_IM':
            ka = 0
            kd = params['kd'].value
            kid = params['kid'].value
            n = int(params['n'].value)
            Frac = params[f'Frac{i+1}'].value
            Vd = params['Vd'].value

        else:  # PK_SC, PEG, AQ, ...
            ka = 0
            kd = params['kd'].value
            kid = params['kid'].value
            n = int(params['n'].value)
            Frac = params['Frac'].value
            Vd = params['Vd'].value

        '''# Debug-Ausgabe pro Datensatz
        print(f"[PK-Fit] Dataset {i+1}")
        print(f"         Mode = {mode}")
        print(f"         Dose = {doses[i] if i < len(doses) else doses[0]} mg")
        print(f"         z0   = {z0[i]}")
        print(f"         ka   = {ka}")
        print(f"         ke   = {ke}")
        print(f"         kd   = {kd}")
        print(f"         kid  = {kid}")
        print(f"         n    = {n}")
        print(f"         Frac = {Frac}")
        print(f"         Vd   = {Vd}\n")'''

        args = (ka, ke, kd, kid, Vd, n, Frac)

        sol = solve_ivp(
            ode_system,
            t_span=(0, times[i][-1]),
            y0=z0[i],
            t_eval=times[i],
            args=args,
            method='LSODA'
        )

        conc = sol.y[2] * conversion_factor_ng_ml
        log_conc = [np.log10(c) if c > 0 else -8 for c in conc]
        solutions.append(log_conc)

    return solutions



def solve_pd_ode(z0, t, params):
    """
    Solve PD ODE model for multiple datasets.

    Args:
        z0 (list of list): Initial conditions per dataset
        t (list of list): Timepoints per dataset (days)
        params (lmfit.Parameters): Model parameters

    Returns:
        list of list: Simulated log10 VL fold changes
    """
    solution = []
    for i in range(len(z0)):
        sol = solve_ivp(
            ode_complex,
            (0, t[i][-1]),
            z0[i],
            t_eval=t[i],
            args=(params, params['IC50'].value,)
        )
        V_total = 2 * (sol.y[3] + sol.y[7]) / (50 * 9600 + 3100) # VL in plasma
        V_change = V_total / V_total[0]
        solution.append(np.log10(V_change))
    return solution

def solve_vl(z, t, params):
    """
    Compute viral load (VL) change from initial patient-specific baseline VL.

    Args:
        z (list): Initial conditions
        t (list): Time points
        params: Model parameters

    Returns:
        list of float: Change from baseline viral load
    """
    VL_profiles = []
    
    for i in range(len(z)):
        sol = solve_ivp(
            ode_complex,
            (0, t[i][-1]),
            z[i],
            t_eval=t[i],
            args=(params, params['IC50'].value)
        )
        V_total = sol.y[3] + sol.y[7]
        solution_vl = np.log10(2 * V_total / (50 * 9600 + 3100)) # In plasma
        VL_profiles.append((4.5 - solution_vl[0]))
    return VL_profiles


def residuals(params, initial_conditions, times, concentrations, mode, doses):
    """
    Compute residuals for model vs. observed PK data for fitting.

    Args:
        params (lmfit.Parameters): Model parameters
        initial_conditions (list of list): Initial values for ODE
        times (list of list): Timepoints
        concentrations (list of list): Observed concentrations (log10)
        mode (str): Dosing route
        doses (list of int): Doses used

    Returns:
        list of float: Flattened residual vector
    """
    solutions = solve_pk_ode(initial_conditions, times, params, doses, mode)
    residual = []

    for j in range(len(solutions)):
        model = np.array(solutions[j])
        obs = np.array(concentrations[j])
        residual += list(obs - model)

    return residual

def steady_states(params):
    """
    Simulate PD system until steady state to initialize state variables.

    Args:
        params (lmfit.Parameters): Model parameters

    Returns:
        list: Steady-state values of viral dynamics system
    """
    lambdaT = params['lambdaT']
    deltaT = params['deltaT']
    lambdaM = params['lambdaM']

    v0 = [lambdaT/deltaT, 0, 0, 10, lambdaM/params['deltaM'], 0, 0, 0, 0, 0]
    t_obs = np.linspace(0, 200, num=1000)

    res = solve_ivp(
        ode_complex,
        (0, t_obs[-1]),
        v0,
        t_eval=t_obs,
        args=(params, 5),
        method='LSODA'
    )

    return [res.y[i][-1] for i in range(8)]  # Tu0, T10, T20, V0, Mu0, M10, M20, VNI0


def residual_CV(params, data_time, data_mat, doses):
    """
   Compute residuals for LEN viral decay in plasma.

   Args:
       params (lmfit.Parameters): Model parameters
       data_time (list of list): Timepoints per dataset
       data_mat (list of list): Log10 fold change per dataset
       doses (list of int): Administered doses

   Returns:
       ndarray: Residual vector
   """

    Tu0, T10, T20, V0, Mu0, M10, M20, VNI0 = steady_states(params)

    z0 = []
    for k in range(len(doses)):
        initial_SC = doses[k] * 1e6
        z0.append([Tu0, T10, T20, V0, Mu0, M10, M20, VNI0, initial_SC, 0])

    deltaCV = solve_pd_ode(z0, data_time, params)
    VL = solve_vl(z0, data_time, params)

    resid = []
    for i in range(len(deltaCV)):
        resid.append(np.subtract(deltaCV[i], data_mat[i]))

    for vl in VL:
        resid.append(np.array([vl]))

    return np.concatenate(resid)

def plot_results(route, results_df, data_conc, data_time, z0, doses, output_dir):
    """
    Plot the best fitting result vs. clinical data (by minimal AIC).

    Args:
        route (str): Dosing mode (e.g., 'PK_PO', 'PK_SC_PEG', ...)
        results_df (pd.DataFrame): DataFrame with AIC, RSS, and parameter columns
        data_conc (list of list): Observed concentrations or changes (log10)
        data_time (list of list): Observed timepoints for each dose (hours for PK, days for PD)
        z0 (list of list): Initial conditions for the ODE system
        doses (list of int): Sorted list of unique base dose values
        output_dir (str or Path): Location to save the plot
    """

    best_result = results_df.iloc[0]
    best_params = Parameters()

    for param, value in best_result.items():
        if param not in ['AIC', 'RSS']:
            best_params.add(param, value=value)

    # Solve ODE & convert hours to days for plotting
    if route.startswith("PK"):
        solutions_best = solve_pk_ode(z0, data_time, best_params, doses, route)
        data_time_plot = [[t / 24 for t in l] for l in data_time]
    else:
        solutions_best = solve_pd_ode(z0, data_time, best_params)
        data_time_plot = data_time  # Already in days

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    print(len(solutions_best))

    for i in range(len(solutions_best)):
        if route == "PK_IM":
            dose_str = f"{doses[0]} mg"
            label_model = f"{dose_str} (model {i + 1})"
        else:
            dose_str = f"{doses[i]} mg"
            label_model = f"{dose_str} (model)"

        ax.scatter(data_time_plot[i], data_conc[i], alpha=0.7)
        ax.plot(data_time_plot[i], solutions_best[i], linestyle='-', linewidth=2, label=label_model)

    ax.set_xlabel("Time post-dose (days)")

    if route.startswith("PK"):
        ax.set_ylabel(r"LEN plasma concentration (log$_{10}$ ng/mL)")
    else:
        ax.set_ylabel(r"HIV-1 RNA (log$_{10}$ fold change)")

    ax.legend()
    plt.tight_layout()

    output_file = os.path.join(output_dir, f"Fitting_{route}.png")
    plt.savefig(output_file)
    print(f"Plot saved to: {output_file}")
    plt.show()
