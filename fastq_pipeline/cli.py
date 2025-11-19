import argparse
from barcode.barcode_utils import *
from filtering.filter import *
from analysis.statistics import *


def argparser():
    parser = argparse.ArgumentParser(
        prog='fastq_pipeline', 
        description='Extract barcode, calculate stats, show stats visualization, and filter reads to create new fastq file.'
    )

    subparsers = parser.add_subparsers(
        title='Available Commands', 
        description="Please choose command below:",
        dest='command',  
        required=True,
    )

    fastq_pipeline = subparsers.add_parser('filteredProg', help='A utility that extracts barcodes from compressed FASTQ files and generates a filtered output file based on a statistically derived quality cutoff.')
    fastq_pipeline.add_argument('-i', '--input-file', type=str, default=None, required=True, help="Input fastq.gz file")
    fastq_pipeline.add_argument('-o', '--output-dir', type=str, default='../data/processed_data', required=False, help="Output directory")
    
    barcode_extract = subparsers.add_parser('extractBQ', help='Extract barcode from fastq.gz to new file')
    barcode_extract.add_argument('-i', '--input-file', type=str, default=None, required=True, help="Input fastq.gz file")
    barcode_extract.add_argument('-o', '--output-dir', type=str, default='../data/processed_data', required=False, help="Output directory")
    
    cal_stat = subparsers.add_parser('statCal', help='Calculate statistics to design cut-off')
    cal_stat.add_argument('-d', '--input-dir', type=str, default=None, required=True, help="Input processed fastq directory")
    #cal_stat.add_argument('-o', '--output-dir', type=str, default='../../data/processed_data', required=False, help="Output directory")
    
    filter_file = subparsers.add_parser('filterRead', help='Filter read to new fastq file')
    filter_file.add_argument('-d', '--input-dir', type=str, default=None, required=True, help="Input stat directory")
    filter_file.add_argument('-b', '--barcode-dir', type=str, default='../data/processed_data', required=False, help="Input extracted barcode directory")
    
    return parser
