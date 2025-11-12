import pandas as pd
from Bio import SeqIO
import gzip
import os
import matplotlib.pyplot as plt
import sys 

# 1. Core Functions: Calculation and Summary 

def calculate_read_stats(fastq_file, barcode_id):
    """
    Parses a gzipped FASTQ file and calculates read length and mean Phred score
    for each read, returning the results as a Pandas DataFrame.
    """
    stats_data = []
    try:
        with gzip.open(fastq_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_length = len(record.seq)
                
                if record.letter_annotations.get("phred_quality"):
                    qual_scores = record.letter_annotations["phred_quality"]
                    mean_quality = sum(qual_scores) / len(qual_scores)
                else:
                    mean_quality = 0 
                
                stats_data.append({
                    'Read_ID': record.id,
                    'Barcode': barcode_id,
                    'Seq_Length': read_length,
                    'Mean_Phred_Score': mean_quality,
                })
    except Exception as e:
        print(f"Error processing {fastq_file}: {e}", file=sys.stderr)
        return pd.DataFrame()
        
    return pd.DataFrame(stats_data)

# calculate 25th 50th 75th of sequence data
def calculate_summary_stats(df: pd.DataFrame):
    """
    Calculates 25th, 50th (median), and 75th percentiles (Q25, Q50, Q75) for key metrics.
    """
    if df.empty:
        return None
    length_summary = df['Seq_Length'].describe(percentiles=[.25, .5, .75])
    quality_summary = df['Mean_Phred_Score'].describe(percentiles=[.25, .5, .75])

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


def save_stats_to_csv(df, output_dir, barcode):
    """Saves a DataFrame of per-read statistics to a CSV file."""
    output_csv = os.path.join(output_dir, f"{barcode}_stats.csv")
    df.to_csv(output_csv, index=False)
    print(f"-> Saved read statistics to: {output_csv}")
    
    
# 2. Plotting Functions 

def plot_individual_histograms(df, output_dir, barcode):
    """Generates and saves separate histograms for length and quality for one barcode."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram of Read length
    df['Seq_Length'].plot(kind='hist', bins=50, title=f'{barcode} Read Length Distribution', ax=axes[0])
    axes[0].set_xlabel("Read Length (bases)")
    
    # Histogram of Mean Quality Score  
    df['Mean_Phred_Score'].plot(kind='hist', bins=50, title=f'{barcode} Mean Phred Quality Score Distribution', ax=axes[1])
    axes[1].set_xlabel("Mean Quality Score (Phred)")
    
    plt.tight_layout()
    output_plot = os.path.join(output_dir, f'{barcode}_histograms.png')
    plt.savefig(output_plot)
    plt.close(fig) 
    print(f"-> Saved histograms to: {output_plot}")
    

def plot_combined_boxplots(all_stats_dfs, output_dir):
    """Generates and saves combined box plots for all barcodes/DFs."""
    if not all_stats_dfs:
        print("Warning: No data to generate combined plots.")
        return

    print("Generating combined Box Plots for barcode comparison...")
    final_df_combined = pd.concat(all_stats_dfs, ignore_index=True)
    
    # Box Plot of Read Length 
    plt.figure(figsize=(10, 6))
    final_df_combined.boxplot(column='Seq_Length', by='Barcode', grid=False, medianprops=dict(color='red'), ax=plt.gca())
    plt.suptitle('Read Length Comparison Across Barcodes')
    plt.title('') 
    plt.ylabel('Read Length (bases)')
    plt.xlabel('Barcode')
    plt.tight_layout()
    output_plot_length = os.path.join(output_dir, 'combined_length_boxplot.png')
    plt.savefig(output_plot_length)
    plt.close()
    print(f"-> Saved combined length boxplot to: {output_plot_length}")

    # Box Plot of Mean Quality Score ---
    plt.figure(figsize=(10, 6))
    final_df_combined.boxplot(column='Mean_Phred_Score', by='Barcode', grid=False, medianprops=dict(color='red'), ax=plt.gca())
    plt.suptitle('Mean Quality Score Comparison Across Barcodes')
    plt.title('')
    plt.ylabel('Mean Quality Score (Phred)')
    plt.xlabel('Barcode')
    plt.tight_layout()
    output_plot_quality = os.path.join(output_dir, 'combined_quality_boxplot.png')
    plt.savefig(output_plot_quality)
    plt.close()
    print(f"-> Saved combined quality boxplot to: {output_plot_quality}")


# 3. Execution Functions 

def run_stats(args):
    """
    Entry point for the 'stats' command. Calculates per-read and summary statistics.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    print("--- Running STATS Command: Calculate Statistics and Save CSV ---")
    
    for file_path in args.input:
        if os.path.exists(file_path) and file_path.endswith('.fastq.gz'):
            base_name = os.path.basename(file_path)
            if not base_name.endswith('.fastq.gz'):
                 print(f"Warning: Skipping {file_path}. Input file name should end with '.fastq.gz'.", file=sys.stderr)
                 continue

            barcode = base_name.split('.')[0] 
            print(f"Processing {barcode} from file: {file_path}")
            
            df = calculate_read_stats(file_path, barcode)
            
            if not df.empty:
                # 1. Save the per-read data
                save_stats_to_csv(df, args.output_dir, barcode)
                
                # 2. Calculate and display the statistic summary 
                summary = calculate_summary_stats(df)
                if summary:
                    print("\n  Summary Statistics:")
                    print(f"  Total Reads: {summary['Total_Reads']:,}")
                    print("  --------------------------------------")
                    print("  Metric         | Q25 | Median | Q75")
                    print("  --------------------------------------")
                    print(f"  Length (bases) | {summary['Length_Q25']:.0f} | {summary['Length_Median']:.0f} | {summary['Length_Q75']:.0f}")
                    print(f"  Quality (Phred)| {summary['Quality_Q25']:.2f} | {summary['Quality_Median']:.2f} | {summary['Quality_Q75']:.2f}")
                    print("  --------------------------------------\n")
        else:
            print(f"Error: Invalid or non-gzipped file path skipped: {file_path}", file=sys.stderr)
            
    print("\n✅ Stats Command Complete.")


def run_plot(args):
    """Entry point for the 'plot' command."""
    os.makedirs(args.output_dir, exist_ok=True)
    all_stats_dfs = [] 
    
    print("--- Running PLOT Command: Generating Visualizations ---")
    
    for csv_file in args.input:
        try:
            df = pd.read_csv(csv_file)
            if df.empty:
                print(f"Warning: CSV file {csv_file} is empty. Skipping.")
                continue

            base_name = os.path.basename(csv_file)
            barcode = base_name.split('_stats.csv')[0] 
            df['Barcode'] = barcode 
            
            print(f"Processing {barcode} from file: {csv_file}")
            
            plot_individual_histograms(df, args.output_dir, barcode)
            all_stats_dfs.append(df)
            
        except Exception as e:
            print(f"Error reading or plotting {csv_file}: {e}", file=sys.stderr)

    if all_stats_dfs:
        plot_combined_boxplots(all_stats_dfs, args.output_dir)

    print("\n Plot Command Complete.")
    
    
def run_filter(args):
    """Placeholder for the filtering function (Step 3)."""
    
    print("\n Filtering command (Step 3) not yet implemented!")
    print(f"Criteria: Min Length={args.min_len}, Min Quality={args.min_qual}")
    print("Implement the logic here to filter reads from the input FASTQ.GZ and write to the output FASTQ.GZ.")
