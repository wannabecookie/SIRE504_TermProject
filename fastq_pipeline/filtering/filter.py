import pandas as pd
import os
import glob
from Bio import SeqIO
import gzip


def filter_csv(input_dir):
    csv_files = glob.glob(os.path.join(input_dir, "*.csv")) #find all csv files in input_dir -> list of csv files
    output_dir = os.path.join(input_dir, "filtered") #create subdirectory inside input_dir
    os.makedirs(output_dir, exist_ok=True) #create output directory 
    pass_id = {} #to store key:barcode id and value:list of seqid
    for files in csv_files: 
        file_name = os.path.basename(files) #get only filename
        print(f"\nWorking on: {file_name}")
        
        try:    
            seq_len = int(input(f"Enter minimum cut-off sequence length for {file_name}: "))
            quality = float(input(f"Enter minimum Q-score for {file_name}: "))
        except ValueError:
            print("Invalid input, please enter number (cut-off: integer only, q-score: decimals allowed)")
            seq_len = int(input(f"Enter minimum cut-off sequence length for {file_name}: "))
            quality = float(input(f"Enter minimum Q-score for {file_name}: "))
            
        
        df = pd.read_csv(files)
        filtered_csv = df.query('`Seq_Length` > @seq_len and `Mean_Phred_Score` > @quality')
        total_read = len(df)
        total_filtered = len(filtered_csv)
        print(f"Total reads: {total_read}, Passed reads: {total_filtered}")
        
        output_file = os.path.join(output_dir, f"filtered_{file_name}")
        filtered_csv.to_csv(output_file, index=False) #write filtered to csv file
        print(f"Saved filtered csv to: {os.path.join(output_file)}")

        barcode = file_name.replace('_stats.csv','') #to get barcode id as key value, extract from file name
        pass_id_list = filtered_csv['Read_ID'].tolist()
        pass_id[barcode] = pass_id_list #store value of each barcode key as list of pass read id
    
    return pass_id

    
def filtered_fastq(pass_id, barcode_dir):
    
    fastq_files = glob.glob(os.path.join(barcode_dir, "*.fastq.gz")) #get list of files path in barcode_dir, end with fastq.gz
    output_dir = os.path.join("../data", "output_filtered_fastq")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n===Creating filtered fastq.gz file===")
    
    for files in fastq_files: 
        file_name = os.path.basename(files).replace('.fastq.gz', '') #get the file name to get only barcode id
        pass_ids = pass_id[file_name] #pass_ids store value as list of barcode that match from barcode file, and filtered csv barcode
        pass_ids_set = set(pass_ids) #set will make the program faster
        count = 0

        output_file_path = os.path.join(output_dir, f"filtered_{file_name}.fastq.gz")
        
        input_file = gzip.open(files, 'rt')
        output_file = gzip.open(output_file_path, 'wt')
        
        for record in SeqIO.parse(input_file, "fastq"):
            if record.id in pass_ids_set: #direct look up for record id in pass ids
                SeqIO.write(record, output_file, "fastq")
                count +=1
                if count % 10000 == 0:
                    print(f"{file_name}: Processed {count} reads..", end='\r', flush=True)
        
        input_file.close()
        output_file.close()
        
        print(f"Write {count} read to {os.path.basename(output_file_path)}")
        print(f"Files available at folder: {os.path.join(output_dir)}")
    
      
if __name__ == "__main__":
    import sys
    
    csv_dir = sys.argv[1]
    barcode_fqdir = sys.argv[2]
    
    pass_id = filter_csv(csv_dir)
    filtered_fastq(pass_id, barcode_fqdir)