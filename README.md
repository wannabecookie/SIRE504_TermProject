# 🧬 Fastq Data Analysis Pipeline

## Table of Contents
1.  [🚀 Overview](#1.-Overview)
2.  [📦 Project Structure](#-project-structure)
3.  [✨ Features](#-features)
4.  [🛠️ Installation](#%EF%B8%8F-installation)
5.  [⚙️ Usage](#%EF%B8%8F-usage)
    * [Full Pipeline: `filteredProg`](#full-pipeline-filteredprog)
    * [Step-by-Step Commands](#step-by-step-commands)
6.  [💡 Contacts](#Contacts)

---

## 🚀 1. Overview
The **`fastq_pipeline`** is a command-line utility for processing large gzipped FASTQ files, primarily designed for data where barcode information is embedded in the sequence description (e.g., Nanopore sequencing data).

The pipeline performs three main functions:
1.  **Barcode Extraction:** Splits a single input FASTQ file into multiple gzipped FASTQ files, one for each unique barcode.
2.  **Statistical Analysis & Visualization:** Calculates per-read statistics (sequence length, mean Phred quality score) for each barcode, saves them to CSV, and generates insightful plots (histograms, boxplots).
3.  **Quality Filtering:** Allows the user to apply custom cut-offs (minimum sequence length and mean Q-score) based on the statistical analysis, generating a final, filtered set of high-quality FASTQ files.

---

## 📦 2. Project Structure
The core functionality is organized into three modules: `barcode`, `filtering`, and `analysis`.

```bash
fastq_pipeline/
├── main.py             # Main entry point, command-line arguments.
├── cli.py              # Defines the command-line interface.
├── barcode/
│   └── barcode_utils.py  # Functions for barcode extraction and file grouping.
├── filtering/
│   └── filter.py         # Functions for applying statistical filters and creating filtered FASTQ.
└── analysis/
    └── statistics.py     # Functions for calculating read stats, summarized, and plotting.

```
---

## ✨ 3. Features
* **Handle Compressed Data:** Directly reads and writes gzipped FASTQ files (`.fastq.gz`).
* **Barcode-Specific Processing:** Groups reads by embedded barcode for isolated analysis.
* **Statistical Reporting:** Generates per-read CSV files containing **Read ID**, **Sequence Length**, and **Mean Phred Score**.
* **Visualization:** Creates **Histograms** and **Boxplots** for length and quality distributions, aiding in manual quality-cutoff selection. * **Custom Filtering:** Filters reads based on user-defined minimum sequence length and mean quality score.

---

## 🛠️ 4. Installation

### Prerequisites
The project requires **Python 3.x** and the following libraries:

* **Biopython (`Bio`)**: For handling sequence data (FASTQ parsing).
* **Pandas**: For statistical data manipulation (DataFrames, CSV).
* **Matplotlib (`matplotlib`)**: For generating statistical plots.

Or you can install the required packages using conda:

```bash
conda env create -y -f environment.yml
```
---

## ⚙️ 5. Usage

The pipeline is executed via **main.py** using subcommands defined in **cli.py**.
## Full Pipeline: `filteredProg`
This command runs all steps sequentially: **Barcode Extraction** $\rightarrow$ **Statistics Calculation** $\rightarrow$ **Sequences Filtering**, outputting a filtered fastq.gz file.

---

### Bash Command Example
```bash
python main.py filteredProg -i <input_file.fastq.gz> [-o] <output_directory>
```
### ⚙️ Additional Arguments for `filteredProg`

| Argument | Description | Default | Required |
| :--- | :--- | :--- | :--- |
| **-i, --input-file** | The single input FASTQ file (must be .fastq.gz). | None | Yes |
| **-o, --output-dir** | Directory where intermediate (grouped FASTQ) files are saved. | ../data/processed_data | No |

## Step-by-step commands
you can run each step individually following:

## 5.1. Extract Barcodes: `extractBQ`
Splits the input .fastq.gz file into multiple .fastq.gz files based on the barcode= tag in the read description.

```bash
python main.py extractBQ -i <input_file.fastq.gz> [-o] <output_directory>
```

### ⚙️ Additional Arguments for `extractBQ`

| Argument | Description | Default | Required |
| :--- | :--- | :--- | :--- |
| **-i, --input-file** | The single input FASTQ file (must be .fastq.gz). | None | Yes |
| **-o, --output-dir** | Directory where intermediate (grouped FASTQ) files are saved. | ../data/processed_data | No |


## 5.2. Calculate Statistics & Plot: `statCal`
Calculates stat. for all of the FASTQ files in the input directory and create CSV files and plots(histogram and boxplot).

```bash
python main.py statCal -d <directory_with_grouped_fastqs>
```

### ⚙️ Additional Arguments for `statCal`

| Argument | Description | Default | Required |
| :--- | :--- | :--- | :--- |
| **-d, --input-dir** | Directory containing the grouped .fastq.gz files from extractBQ. | None | Yes |

## 5.3. Filter Reads: `filterRead`
Reads the statistics CSV files, ask the user for length and quality cut-offs per barcode, and creates new filtered .fastq.gz files based on that values

```bash
python main.py filterRead -d <directory_with_stat_csvs> [-b] <directory_with_grouped_fastqs>
```

### ⚙️ Additional Arguments for `filterRead`

| Argument | Description | Default | Required |
| :--- | :--- | :--- | :--- |
| **-d, --input-dir** | Directory contain the *_stats.csv files from statCal. | None | Yes |
| **-b, --barcode-dir** | Directory contain the original grouped .fastq.gz files. | ../data/processed_data | No |


## 💡 6. Contacts
Saranyaporn Kunawongkrit - [@wannabecookie](https://github.com/wannabecookie) - Saranyaporn.kun@student.mahidol.edu

Wongsatorn Thummawong - [@WONGSATORN-SIMB](https://github.com/WONGSATORN-SIMB) - Wongsatorn.thu@student.mahidol.edu
