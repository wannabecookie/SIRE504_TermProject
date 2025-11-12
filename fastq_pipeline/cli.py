import argparse

# Import the execution functions from statistics module
from fastq_pipeline.analysis.statistics import run_stats, run_plot, run_filter

def argparser() -> argparse.ArgumentParser:
    """Defines the command-line interface structure using subparsers."""
    parser = argparse.ArgumentParser(
        prog='seq_analysis.py', 
        description='Nanopore QC Tool: Calculate stats, plot, and filter reads.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    subparsers = parser.add_subparsers(
        title='Available Commands', 
        dest='command', 
        required=True, 
        help="The action to perform: 'stats', 'plot', or 'filter'."
    )

    # 1. Setup for the 'stats' command 
    parser_stats = subparsers.add_parser(
        'stats', 
        help='Calculate read statistics (length, quality) from FASTQ.GZ files and save to CSV.'
    )
    parser_stats.add_argument(
        '-i', '--input', 
        nargs='+', 
        required=True,
        help="One or more paths to input FASTQ.GZ files (e.g., barcode*.fastq.gz)"
    )
    parser_stats.add_argument(
        '-o', '--output_dir', 
        type=str,
        default='read_stats_output',
        help="Directory to save output CSV files (default: read_stats_output)"
    )
    parser_stats.set_defaults(func=run_stats) 

    # 2. Setup for the 'plot' command 
    parser_plot = subparsers.add_parser(
        'plot', 
        help='Generate histograms and boxplots from pre-calculated CSV statistics files.'
    )
    parser_plot.add_argument(
        '-i', '--input', 
        nargs='+', 
        required=True,
        help="One or more paths to input CSV files (e.g., read_stats_output/*.csv)"
    )
    parser_plot.add_argument(
        '-o', '--output_dir', 
        type=str,
        default='read_plots_output',
        help="Directory to save output PNG plots (default: read_plots_output)"
    )
    parser_plot.set_defaults(func=run_plot)
    
    # 3. Setup for the 'filter' command 
    parser_filter = subparsers.add_parser(
        'filter', 
        help='Filter reads based on cut-off criteria (length/quality) and save a new FASTQ file.'
    )
    parser_filter.add_argument(
        '-i', '--input', 
        required=True,
        help="Path to the input FASTQ.GZ file to be filtered."
    )
    parser_filter.add_argument(
        '-o', '--output_file', 
        required=True,
        help="Path for the filtered output FASTQ.GZ file."
    )
    # parser_filter.add_argument(
    #     '--min_len', 
    #     type=int, 
    #     default=500,
    #     help="Minimum read length cut-off (default: 500 bases)"
    # )
    # parser_filter.add_argument(
    #     '--max_len', 
    #     type=float, 
    #     default=100000000.0,
    #     help="Maximum read length cut-off (default: 10000000 bases)"
    # )
    # parser_filter.add_argument(
    #     '--min_qual', 
    #     type=float, 
    #     default=10.0,
    #     help="Minimum mean Phred quality score cut-off (default: 7.0)"
    # )
    # Link the command to the run_filter function
    parser_filter.set_defaults(func=run_filter) 
    
    return parser


def main(argv=None):
    """Main function to parse arguments and execute the selected command."""
    parser = argparser()
    args = parser.parse_args(argv)
    
    # args.funcion are run_stats, run_plot, or run_filter, as set above.
    args.func(args)


if __name__ == "__main__":
    # If cli.py is run directly
    main()