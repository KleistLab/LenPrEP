from cli import parse_arguments
from utils import get_doses, read_data, read_ini_file, initialize_from_ini, plot_results
from pk_fitting import run_fitting as run_pk_fitting
from pd_fitting import run_fitting as run_pd_fitting
import sys

def run_model():
    args = parse_arguments()

    doses, dose_counts = get_doses(args.data_path)

    data = read_data(args.mode, args.data_path, doses, dose_counts)

    # Run fitting and get results as DataFrame
    if args.mode.startswith('PK'):
        df, z0 = run_pk_fitting(doses, data, args.mode, args.param_file, args.max_count, dose_counts, args.method, args.track_iterations)
    elif args.mode.startswith('PD'):
        df, z0 = run_pd_fitting(doses, data, args.mode, args.param_file, args.max_count, args.method, args.track_iterations)
        # Only save fitted parameters
        df_reduced = df[['AIC', 'RSS', 'IC50', 'm', 'betaT']]

    if args.mode.startswith('PD'):
        print('\n', df_reduced.head())
    else:
        print('\n', df.head())

    # Plot optional
    if args.plot:
        # Use the best fit parameter values (minimal AIC)
        best_params = df.iloc[0].to_dict()

        ini_params = read_ini_file(args.param_file)
        param_obj = initialize_from_ini(ini_params, doses, args.mode)
        #print(param_obj)

        for key, val in best_params.items():
            if key in param_obj:
                param_obj[key].value = val

        times = [entry['time'] for entry in data]
        concentrations = [entry['conc'] if args.mode.startswith('PK') else entry['change'] for entry in data]

        # Plotting function
        plot_results(
            route=args.mode,
            results_df=df,
            data_conc=concentrations,
            data_time=times,
            z0=z0,
            doses=doses,
            output_dir=args.output
        )

        # Create output folder if it doesn’t exist
        if not args.output.exists():
            print(f"Creating output folder: {args.output}")
            args.output.mkdir(parents=True, exist_ok=True)

        # Save results to CSV
        result_file = args.output / f"{args.mode}_results.csv"

        # Check if output file already exist
        if result_file.exists():
            answer = input(f"File {result_file} already exists. Overwrite? (y/n): ").strip().lower()
            if answer != 'y':
                print("Aborted by user. Exiting.")
                sys.exit(0)

        if args.mode.startswith('PD'):
            df_reduced.to_csv(result_file, index=False)
        else:
            df.to_csv(result_file, index=False)
        print(f"Results saved to: {result_file}")


if __name__ == '__main__':
    run_model()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)


