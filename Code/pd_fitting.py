from utils import read_ini_file, initialize_from_ini, residual_CV, steady_states, randomize_params, per_iteration
from lmfit import minimize
import pandas as pd


def run_fitting(doses, data, mode, param_file, max_count, method, track_iterations):
    data_time = [entry['time'] for entry in data]
    data_mat = [entry['change'] for entry in data]

    # Read .ini file
    ini_params = read_ini_file(param_file)

    # Initialize Parameters as lmfit object
    base_params = initialize_from_ini(ini_params, doses, mode)

    #for name, p in base_params.items():
    #    print(f"{name}: value={p.value}, vary={p.vary}")

    results = []

    for iteration in range(max_count):
        print(f"\nIteration {iteration + 1} / {max_count}")

        # Randomize parameter start values
        params = randomize_params(base_params)

        # Fit
        result = minimize(
            residual_CV,
            params,
            args=(data_time, data_mat, doses),
            iter_cb=per_iteration if track_iterations and method == 'least_squares' else None,
            method=method
        )

        rss = sum(result.residual[:-4]**2)
        result_dict = {
            'AIC': result.aic,
            'RSS': rss,
        }

        result_dict.update({k: p.value for k, p in result.params.items()})
        results.append(result_dict)

    df = pd.DataFrame(results).sort_values(by="AIC")

    # Save current parameters for plotting
    z0_base = steady_states(base_params)
    z0 = []
    for k in range(len(doses)):
        z0.append(z0_base + [doses[k] * 1e6, 0])

    return df, z0
