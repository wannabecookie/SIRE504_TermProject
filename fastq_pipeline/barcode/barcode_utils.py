import gzip
import re
from Bio import SeqIO
import os


def extract_barcode(record):
    match = re.search(r'barcode=(\S+)', record.description)
    return match.group(1) #barcode


def grouped_fastq(input_file, output_dir): 
    os.makedirs(output_dir, exist_ok=True)
    barcode_files = {} #to store barcode as key and value as file handle
    barcode_counts = {} #store barcode counts
    total_reads = 0
    
    print(f"===Starting grouping barcode===")
    with gzip.open(input_file, 'rt') as file: #open input gz file
        for record in SeqIO.parse(file, 'fastq'):
            total_reads += 1
                    
            if total_reads % 10000 == 0: #show progress every 10000 reads
                print(f"Processed {total_reads} reads..", end='\r')

            barcode = extract_barcode(record)
                                     
            if barcode not in barcode_files: #create file handle if this is a new barcode file
                output_file = os.path.join(output_dir, f"{barcode}.fastq.gz")
                barcode_files[barcode] = gzip.open(output_file, 'wt') #open compress file to store as file handle of barcode key
                print(f"Creating output file: {output_file}")
            
            SeqIO.write(record, barcode_files[barcode], 'fastq') #write record
            
            barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
            
        for file_handle in barcode_files.values(): #close file
            file_handle.close()
                    
        print(f"\nProcessed {total_reads} total reads")
    
    #Print summary
    print("\n=== Grouping Barcode Summary ===")
    for barcode, count in barcode_counts.items():
        print(f"{barcode}: {count} reads")
    print(f"Total barcodes: {len(barcode_counts)}")
    print(f"Total reads processed: {sum(barcode_counts.values())}")
    
    return output_dir
    

if __name__ == "__main__":
    import sys
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    grouped_fastq(input_file, output_dir)