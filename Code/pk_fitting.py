from utils import create_initial_conditions, read_ini_file, initialize_from_ini, residuals, randomize_params, per_iteration
from lmfit import minimize
import numpy as np
import pandas as pd


def run_fitting(doses, data, mode, param_file, max_count, dose_counts, method, track_iterations):

    # Read .ini file
    ini_params = read_ini_file(param_file)

    # Initialize Parameters as lmfit object
    base_params = initialize_from_ini(ini_params, doses, mode)

    #for name, p in base_params.items():
    #    print(f"{name}: value={p.value}, vary={p.vary}")

    # Initial condition for ODE
    z0 = create_initial_conditions(mode, doses, dose_counts)

    times = [entry['time'] for entry in data]
    concentrations = [entry['conc'] if mode.startswith('PK') else entry['change'] for entry in data]

    results_list = []

    for iteration in range(max_count):
        print(f"\nIteration {iteration + 1} / {max_count}")

        # Randomize parameter start values
        params = randomize_params(base_params)

        # Fit
        result = minimize(
            residuals,
            params,
            args=(z0, times, concentrations, mode, doses),
            iter_cb=per_iteration if track_iterations and method == 'least_squares' else None,
            method=method
        )

        rss = np.sum(result.residual ** 2)
        result_dict = {
            'AIC': result.aic,
            'RSS': rss,
        }

        for name, p in result.params.items():
            if name == "n":
                result_dict[name] = int(round(p.value))
            else:
                result_dict[name] = p.value

        results_list.append(result_dict)

    df = pd.DataFrame(results_list).sort_values(by="AIC")
    #print(df.head())

    return df, z0
