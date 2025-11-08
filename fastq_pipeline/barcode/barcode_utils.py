import gzip
import re
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def extract_barcode(record):
    match = re.search(r'barcode=(\S+)', record.description)
    if match:
        return match.group(1) #barcode
    return None 

def grouped_fastq(input_file, output_dir='.'): #output -> current directory
    
    print(f"Starting grouping...", flush=True)
    print(f"Input file: {input_file}", flush=True)
    print(f"Output directory: {output_dir}", flush=True)
    
    output_dir = Path(output_dir)
    
    barcode_files = {} #to store barcode
    barcode_counts = defaultdict(int) #set default new key to be 0
    no_barcode_count = 0
    total_reads = 0
    
    try:
        print("Reading FASTQ file...", flush=True)
        with gzip.open(input_file, 'rt') as file: #open input file
            for index, record in enumerate(SeqIO.parse(file, 'fastq')):
                total_reads += 1
        
                if total_reads % 10000 == 0: #show progress every 10000 reads
                    print(f"Processed {total_reads} reads...", end='\r', flush=True)
                
                barcode = extract_barcode(record)
                
                if barcode:
                    if barcode not in barcode_files: #create file handle if this is a new barcode
                        output_file = output_dir / f"{barcode}.fastq.gz"
                        barcode_files[barcode] = gzip.open(output_file, 'wt')
                        print(f"\nCreating output file: {output_file}", flush=True)
                    
                    SeqIO.write(record, barcode_files[barcode], 'fastq') #write record
                    barcode_counts[barcode] += 1
                else:
                    no_barcode_count += 1
        
        print(f"\nProcessed {total_reads} total reads", flush=True)
    
    finally:
        #close all output files
        print("Closing output files...", flush=True)
        for i in barcode_files.values():
            i.close()
    
    #Print summary
    print("\n=== Grouping Barcode Summary ===", flush=True)
    for barcode, count in sorted(barcode_counts.items()):
        print(f"{barcode}: {count} reads", flush=True)
    if no_barcode_count > 0:
        print(f"No barcode found: {no_barcode_count} reads", flush=True)
    print(f"Total barcodes: {len(barcode_counts)}", flush=True)
    print(f"Total reads processed: {sum(barcode_counts.values())}", flush=True)
    
    return {
        'barcode_counts': dict(barcode_counts),
        'no_barcode_count': no_barcode_count,
        'total_barcodes': len(barcode_counts),
        'total_reads': sum(barcode_counts.values())
    }

if __name__ == "__main__":
    import sys
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    grouped_fastq(input_file, output_dir)