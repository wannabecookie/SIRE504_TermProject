import pandas as pd
from Bio import SeqIO
import gzip
import os
import matplotlib.pyplot as plt
import sys 
import glob 


# 1. Core Functions: Calculation and Summary 

def calculate_read_stats(fastq_file): 
    stats_data = []
    # retrive file's name 
    base_name = os.path.basename(fastq_file)
    file_barcode_id = base_name.split('.')[0] 
    
    try:
        with gzip.open(fastq_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                
                read_length = len(record.seq)
                
                if record.letter_annotations.get("phred_quality"):
                    qual_scores = record.letter_annotations["phred_quality"]
                    mean_quality = sum(qual_scores) / len(qual_scores) if len(qual_scores) > 0 else 0
                else:
                    mean_quality = 0 
                
                stats_data.append({
                    'Read_ID': record.id,
                    'Barcode': file_barcode_id, 
                    'Seq_Length': read_length,
                    'Mean_Phred_Score': mean_quality,
                })
    except Exception as e:
        print(f"Error processing {fastq_file}: {e}", file=sys.stderr)
        return pd.DataFrame(), "" 
        
    stats_df = pd.DataFrame(stats_data)
    return stats_df, file_barcode_id

def calculate_summary_stats(stats_df):
    if stats_df.empty:
        return None
    length_summary = stats_df['Seq_Length'].describe(percentiles=[.25, .5, .75])
    quality_summary = stats_df['Mean_Phred_Score'].describe(percentiles=[.25, .5, .75])

    summary_data = {
        'Total_Reads': int(length_summary['count']),
        'Length_Mean': length_summary['mean'],
        'Length_Q25': length_summary['25%'],
        'Length_Median': length_summary['50%'],
        'Length_Q75': length_summary['75%'],
        'Quality_Mean': quality_summary['mean'],
        'Quality_Q25': quality_summary['25%'],
        'Quality_Median': quality_summary['50%'],
        'Quality_Q75': quality_summary['75%'],
    }
    return summary_data

def save_stats_to_csv(stats_df, output_dir, barcode):
    os.makedirs(output_dir, exist_ok=True) 
    output_csv = os.path.join(output_dir, f"{barcode}_stats.csv")
    stats_df.to_csv(output_csv, index=False)
    print(f"-> Saved read statistics to: {output_csv}")
    return output_dir

# 2. Plotting Functions 

def plot_histograms(stats_df, output_dir, barcode):
    os.makedirs(output_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    stats_df['Seq_Length'].plot(kind='hist', bins=50, title=f'{barcode} Read Length Distribution', ax=axes[0])
    axes[0].set_xlabel("Read Length (bases)")
    
    stats_df['Mean_Phred_Score'].plot(kind='hist', bins=50, title=f'{barcode} Mean Phred Quality Score Distribution', ax=axes[1])
    axes[1].set_xlabel("Mean Quality Score (Phred)")
    
    plt.tight_layout()
    output_plot = os.path.join(output_dir, f'{barcode}_histograms.png')
    plt.savefig(output_plot)
    plt.close(fig) 
    print(f"-> Saved histograms to: {output_plot}")
    

def plot_boxplots(all_stats_dfs, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    if not all_stats_dfs:
        print("Warning: No data to generate plots.")
        return
    
    final_df_combined = pd.concat(all_stats_dfs, ignore_index=True)

    unique_barcodes = final_df_combined['Barcode'].unique()
    
    # Create boxPlots for Each Barcode ---
    for barcode in unique_barcodes:
        
        df_barcode = final_df_combined[final_df_combined['Barcode'] == barcode]
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        # Box Plot of Read Length 
        df_barcode.boxplot(column='Seq_Length', grid=False, 
                           medianprops=dict(color='red'), ax=axes[0])
        axes[0].set_title('Read Length Distribution') 
        axes[0].set_ylabel('Read Length (bases)')
        axes[0].set_xlabel(f'Barcode {barcode}')
        
        # Box Plot of Mean Quality Score 
        df_barcode.boxplot(column='Mean_Phred_Score', grid=False, 
                           medianprops=dict(color='red'), ax=axes[1])
        
        axes[1].set_title('Mean Quality Score Distribution')
        axes[1].set_ylabel('Mean Quality Score (Phred)')
        axes[1].set_xlabel(f'Barcode {barcode}') 
        
        plt.tight_layout()
        output_plot = os.path.join(output_dir, f'{barcode}_boxplots.png')
        plt.savefig(output_plot)
        plt.close(fig)
        print(f"-> Saved boxplots for {barcode} to: {output_plot}")


def stat(input_dir):
    print(f"\n===Start calculating statistic and visualization===")
    if os.path.isdir(input_dir):
    # If input is a directory, find all gzipped FASTQ files
        fastq_files = glob.glob(os.path.join(input_dir, '*.fastq.gz'))
        if not fastq_files:
            print(f"Error: No *.fastq.gz files found in directory: {input_dir}")
            sys.exit(1)
        print(f"Found {len(fastq_files)} FASTQ files in directory.")
    else:
        # If input is a single file
        fastq_files = [input_dir]
        print(f"Processing single file: {input_dir}")
    # -------------------------------
    
    output_test_dir =  os.path.join("../data", "stat_result")
    os.makedirs(output_test_dir, exist_ok=True)
    all_stats_dfs = []

    for i, fastq_file in enumerate(fastq_files):
        
        print(f"\n--- Processing File {i+1}/{len(fastq_files)}: {os.path.basename(fastq_file)} ---")

        stats_df, barcode_id = calculate_read_stats(fastq_file)

        if stats_df.empty:
            print(f"Warning: Skipping {os.path.basename(fastq_file)} due to empty data.")
            continue 
        
        
        print(f"Summary Stats on {barcode_id} data:")
        summary_data = calculate_summary_stats(stats_df)
        print(pd.Series(summary_data).to_string())

        print("\nSaving CSV and Visualization:")
        returned_dir = save_stats_to_csv(stats_df, output_test_dir, barcode_id)
        plot_histograms(stats_df, returned_dir, barcode_id)
        
        all_stats_dfs.append(stats_df)


    # if all_stats_dfs:
    #     print("\n--- Running Combined Plot Test ---")
        plot_boxplots(all_stats_dfs, output_test_dir)
    # else:
    #     print("\nNo valid data found to run combined plots.")

    print("\nTest analysis completed. Check the 'stat_result' folder.")
    
    return output_test_dir


# TEST STANDALONE MODE


if __name__ == "__main__":
    # 0. Argument Check
    if len(sys.argv) < 2:
        print("Usage: python statistics.py <path_to_fastq.gz or path_to_directory_with_fastqs>")
        sys.exit(1)
        
    input_path = sys.argv[1]
    
    print("--- Running Standalone Test ---")
    
    
    # --- File/Directory Handling ---
    if os.path.isdir(input_path):
        # If input is a directory, find all gzipped FASTQ files
        fastq_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))
        if not fastq_files:
            print(f"Error: No *.fastq.gz files found in directory: {input_path}")
            sys.exit(1)
        print(f"Found {len(fastq_files)} FASTQ files in directory.")
    else:
        # If input is a single file
        fastq_files = [input_path]
        print(f"Processing single file: {input_path}")
    # -------------------------------
    
    output_test_dir =  os.path.join("../data", "stat_result")
    os.makedirs(output_test_dir, exist_ok=True)
    all_stats_dfs = [] 
    
    for i, fastq_file in enumerate(fastq_files):
        
        print(f"\n--- Processing File {i+1}/{len(fastq_files)}: {os.path.basename(fastq_file)} ---")

        stats_df, barcode_id = calculate_read_stats(fastq_file)

        if stats_df.empty:
            print(f"Warning: Skipping {os.path.basename(fastq_file)} due to empty data.")
            continue 
        
        
        print(f"2. Summary Stats on {barcode_id} data:")
        summary_data = calculate_summary_stats(stats_df)
        print(pd.Series(summary_data).to_string())

        print("3. Saving CSV and Individual Plots:")
        returned_dir = save_stats_to_csv(stats_df, output_test_dir, barcode_id)
        plot_histograms(stats_df, returned_dir, barcode_id)
        
        all_stats_dfs.append(stats_df)

    if all_stats_dfs:
        print("\n--- Running Plot Test ---")
        plot_boxplots(all_stats_dfs, output_test_dir)
    else:
        print("\nNo valid data found to run plots.")

    print("\nTest analysis complete. Check the 'stat_result' folder.")
    
    pass
