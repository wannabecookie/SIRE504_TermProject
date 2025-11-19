from barcode.barcode_utils import *
from filtering.filter import *
from analysis.statistics import *
from cli import argparser


parser = argparser()
args = parser.parse_args()

if args.command == 'filteredProg': 
    input_file = args.input_file
    output_dir = args.output_dir
    if not input_file.endswith('.fastq.gz'):
        print("wrong file format, please input fastq.gz format")
    else:
        extractBQ_dir = grouped_fastq(input_file, output_dir)
        
    stat_dir = stat(extractBQ_dir)
    pass_id = filter_csv(stat_dir)
    filtered_fastq(pass_id, extractBQ_dir)
    print(f"\n============Program completed============")
    

elif args.command == 'extractBQ':
    input_file = args.input_file
    output_dir = args.output_dir
    if not input_file.endswith('.fastq.gz'):
        print("wrong file format, please input fastq.gz format")
    else:
        grouped_fastq(input_file, output_dir)
        
elif args.command == 'statCal':
    input_dir = args.input_dir
    stat(input_dir)
    
elif args.command == 'filterRead':
    input_dir = args.input_dir
    bq_dir = args.barcode_dir
    pass_id = filter_csv(input_dir)
    filtered_fastq(pass_id, bq_dir)