import argparse
from pathlib import Path
import sys

def str2bool(v):
    return str(v).lower() in ('true', 't', 'yes', '1')

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Run PK/PD model fitting using clinical data and ODE-based simulations.\n\n"
            "You can select from several pharmacokinetic (PK) and pharmacodynamic (PD) models,\n"
            "and fit parameters using local or global optimization methods from lmfit (https://lmfit.github.io/lmfit-py/)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_model.py --mode PK_IM --method least_squares --max_count 50 --plot True\n"
            "  python run_model.py --mode PD_SC --method differential_evolution\n\n"
            "Notes:\n"
            "  - Global optimizers like 'differential_evolution' automatically set --max_count=1.\n"
            "  - Clinical data must be located in Data/ according to the mode.\n"
            "  - Parameter files must exist under Code/ini/MODE.ini\n"
        )
    )

    parser.add_argument(
        '--mode',
        type=str,
        choices=[
            'PK_PO',        # oral (Study 2)
            'PK_SC_PEG',    # SC PEG (Study 4)
            'PK_SC_AQ1',    # SC AQ1 (Study 1)
            'PK_SC_AQ2',    # SC AQ2 (Study 3)
            'PK_IM',        # IM (Study 5)
            'PD_SC'         # PD SC (Study 3)
        ],
        required=True,
        help=(
            "Fitting mode (required):\n"
            "  PK_PO        - PK, oral dosing\n"
            "  PK_SC_PEG    - PK, subcutaneous PEG formulation\n"
            "  PK_SC_AQ1    - PK, subcutaneous aqueous formulation 1\n"
            "  PK_SC_AQ2    - PK, subcutaneous aqueous formulation 2\n"
            "  PK_IM        - PK, intramuscular (IM) with Frac1 & Frac2\n"
            "  PD_SC        - PD, subcutaneous only"
        )
    )

    parser.add_argument(
        '--method',
        type=str,
        default='least_squares',
        choices=[
            'least_squares',         # local
            'nelder',                # local
            'lbfgsb',                # local
            'powell',                # local
            'cg',                    # local
            'newton',                # local
            'cobyla',                # local
            'tnc',                   # local
            'trust-constr',          # local
            'differential_evolution',# global
            'brute',                 # global
            'basinhopping',          # global
            'ampgo',                 # global
        ],
        help=(
            "Optimization method for parameter fitting (default: least_squares).\n\n"
            "Local methods (better for large parameter ranges):\n"
            "  least_squares, nelder, lbfgsb, powell, cg, newton, cobyla, tnc, trust-constr\n\n"
            "Global methods (better for smaller parameter ranges):\n"
            "  differential_evolution, brute, basinhopping, ampgo\n\n"
            "Note: global methods force --max_count=1."
        )
    )

    parser.add_argument(
        '--track_iterations',
        type=str2bool,
        default=False,
        help="If True, prints RSS and parameter values at every optimizer iteration (only for local methods)"
    )


    parser.add_argument(
        '--output',
        type=str,
        default='Results',
        help="Folder to save output files (default: Results/)"
    )

    parser.add_argument(
        '--max_count',
        type=int,
        default=10,
        help="Number of randomized fits to run (ignored for global methods). Default: 10"
    )

    parser.add_argument(
        '--plot',
        type=str2bool,
        default=False,
        help="Whether to plot best fit vs. data after fitting (default: False)"
    )

    # Show help if no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Force max_count = 1 for global optimizers
    global_methods = {'differential_evolution', 'brute', 'basinhopping', 'ampgo'}
    if args.method in global_methods:
        print(f"[Info] Global optimizer '{args.method}' selected — forcing max_count = 1.")
        args.max_count = 1

    # Resolve paths
    project_root = Path(__file__).resolve().parent.parent
    data_base = project_root / 'Data'
    mode_to_data_path = {
        'PK_PO': data_base / 'PK' / 'oral',
        'PK_SC_PEG': data_base / 'PK' / 'SC' / 'SC_PEG',
        'PK_SC_AQ1': data_base / 'PK' / 'SC' / 'SC_AQ1',
        'PK_SC_AQ2': data_base / 'PK' / 'SC' / 'SC_AQ2',
        'PK_IM': data_base / 'PK' / 'IM',
        'PD_SC': data_base / 'PD',
    }

    args.data_path = mode_to_data_path[args.mode]
    args.param_file = project_root / 'Code' / 'ini' / f'{args.mode}.ini'
    args.output = project_root / args.output

    # Show which defaults are used
    provided_args = set(arg.split('=')[0] for arg in sys.argv[1:] if arg.startswith('--'))

    def show_default(name, value):
        if f'--{name}' not in provided_args:
            print(f"[Default] --{name} = {value}")

    show_default('method', args.method)
    show_default('max_count', args.max_count)
    show_default('output', args.output)
    show_default('plot', args.plot)


    return args
