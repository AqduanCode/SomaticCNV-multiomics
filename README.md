# CNVDetect pipeline

<p align="center">
  <b>A Julia-based toolkit for Copy Number Variation detection from sequencing depth data</b>
</p>

# 

## Overview

**CNVDetect pipeline** is a command-line pipeline base on Julia for genome-wide copy number variation (CNV) detection from read depth data. It uses Gaussian Mixture Models (GMM) for parameter estimation and statistical testing for identifying copy number differences between samples.

### Key Features

- **GMM-based estimation**: Robust parameter estimation using Expectation-Maximization algorithm
- **Multi-sample comparison**: Pairwise CNV detection across multiple samples
- **FDR correction**: Benjamini-Hochberg procedure for controlling false discovery rate
- **False negative reduction**: Iterative extension algorithm to minimize missed CNVs
- **Visualization**: Publication-ready plots for CNV regions

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Step 1: Data Loading](#step-1-data-loading)
  - [Step 2: Filtering and Normalization](#step-2-filtering-and-normalization)
  - [Step 3: GMM Parameter Estimation](#step-3-gmm-parameter-estimation)
  - [Step 4: Genome Segmentation](#step-4-genome-segmentation)
  - [Step 5: Copy Number Determination](#step-5-copy-number-determination)
  - [Step 6: Pairwise Comparison](#step-6-pairwise-comparison)
  - [Step 7: Visualization](#step-7-visualization)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [Parameters](#parameters)
- [Examples](#examples)
- [Citation](#citation)
- [Contact](#contact)

## Installation

### Requirements

- Julia 1.6 or higher
- Required packages (auto-installed):
  - CSV
  - DataFrames
  - Plots
  - StatsBase
  - Distributions
  - DelimitedFiles
  - MultipleTesting

### Setup

```bash
# Clone the repository
git clone https://github.com/AqduanCode/SomaticCNV-multiomics
cd SomaticCNV-multiomics

# Start Julia and load the modules
julia
```

```julia
# Load required functions
include("MyFunctions.jl")
include("MyProcess.jl")
```

## Quick Start

```julia
# 1. Load modules
include("MyFunctions.jl")
include("MyProcess.jl")

# 2. Set data path and files
temp_path = "/path/to/your/data/"
file_names = ["sample1.csv", "sample2.csv", "sample3.csv"]

# 3. Process data
combineData = process_genomic_data(temp_path, file_names)
metaDataDF = filtering_normalizing_seqDepth(combineData, [100, 2500], [15, 250])

# 4. Estimate parameters
metaData = Matrix(metaDataDF[:, 4:end])
esti_param = fitting_gmm(metaData, [35, 65, 95, 125, 165, 195], 10)

# 5. Detect CNVs
CNVscaffold = genome_segment(metaDataDF, esti_param, 0.001, 10, 10)
copynumber_info = copynumber_determination(metaDataDF, CNVscaffold, esti_param, 10)
pairedCNV = copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, esti_param, 10)

# 6. Visualize results
plotCNVregionMulti(metaDataDF, CNVscaffold, pairedCNV, [1 2]; chr="chr1")
```

## Usage

### Step 1: Data Loading

Load and combine multiple sample files into a single DataFrame.

```julia
current_dir = pwd()
temp_path = joinpath(current_dir,"sample data")
file_names = [
    "H.1_his_rd_p_1000_chr1.csv",
    "H.1.1_his_rd_p_1000_chr1.csv",
    "H.1.2_his_rd_p_1000_chr1.csv"
]

combineData = process_genomic_data(temp_path, file_names)
```

**Parameters:**

| Parameter    | Type           | Description                           |
| ------------ | -------------- | ------------------------------------- |
| `temp_path`  | String         | Directory path containing input files |
| `file_names` | Vector{String} | List of CSV file names                |

### Step 2: Filtering and Normalization

Filter low-quality fragments and normalize sequencing depth.

```julia
metaDataDF = filtering_normalizing_seqDepth(
    combineData,    # Input DataFrame
    [100, 2500],    # Raw depth filter [min, max]
    [15, 250]       # Normalized depth filter [min, max])
```

**Parameters:**

| Parameter       | Type   | Default     | Description                                |
| --------------- | ------ | ----------- | ------------------------------------------ |
| `depth_filter1` | Vector | [100, 2500] | Total depth threshold before normalization |
| `depth_filter2` | Vector | [15, 250]   | Depth threshold after normalization        |

### Step 3: GMM Parameter Estimation

Estimate Gaussian Mixture Model parameters for each sample.

```julia
metaData = Matrix(metaDataDF[:, 4:end])
temp_mean = [35, 65, 95, 125, 165, 195]  # Initial means for N components(According to your data)
esti_param = fitting_gmm(metaData, temp_mean, 10)

# Visualize estimation results
moving_average_histogram(metaData[:, 1], esti_param[1].means, 50)
```

**Parameters:**

| Parameter    | Type   | Example                | Description                               |
| ------------ | ------ | ---------------------- | ----------------------------------------- |
| `inputData`  | Matrix | -                      | M×N depth matrix (M fragments, N samples) |
| `pre_mean`   | Vector | [35,65,95,125,165,195] | Initial component means                   |
| `windowSize` | Int    | 10                     | Sliding window size                       |

**Optional: Manual Parameter Tuning**

```julia
# Adjust component mean directly
esti_param[1].means[1] = esti_param[1].means[1] - 2

# Use reference sample for adjustment
esti_param = paramFineTuning(esti_param, 1, 2, [2,3,4], 1)
# Args: (params, ref_sample, target_sample, ref_components, target_component)
```

### Step 4: Genome Segmentation

Identify breakpoints and construct CNV examination intervals.

```julia
CNVscaffold = genome_segment(
    metaDataDF,     # Normalized data
    esti_param,     # GMM parameters
    0.001,          # Breakpoint threshold
    10,             # Flanking length
    10              # Minimum fragment length
)
```

**Parameters:**

| Parameter          | Type  | Example | Description                                    |
| ------------------ | ----- | ------- | ---------------------------------------------- |
| `breaks_threshold` | Float | 0.001   | Probability threshold for breakpoint detection |
| `flanking_length`  | Int   | 10      | Flanking region size for comparison            |
| `minimum_length`   | Int   | 10      | Minimum CNV fragment length                    |

### Step 5: Copy Number Determination

Determine copy number for each genomic fragment.

```julia
copynumber_info = copynumber_determination(
    metaDataDF,
    CNVscaffold,
    esti_param,
    10;                    # flanking_length
    thresholdBFactor = 2.0 # Bayes factor threshold
)
```

**Output:** 3D array `[samples × fragments × (copy_number + probabilities)]`

### Step 6: Pairwise Comparison

Compare copy numbers between sample pairs.

```julia
# Basic comparison
pairedCNV = copyNumber_comparison(
    metaDataDF,
    CNVscaffold,
    copynumber_info,
    esti_param,
    10;                          # windowSize
    pvalue_theSame = 0.05e-3,    # p-value threshold
    fdr_threshold = 0.1e-3       # FDR threshold
)

# Reduce false negatives
pairedCNV_FN = copyNumber_comparison_filter_fn(
    metaDataDF,
    CNVscaffold,
    copynumber_info,
    pairedCNV;
    highRDthreshold = 20         # RD difference threshold
)
```

**Parameters:**

| Parameter         | Type  | Default | Description                        |
| ----------------- | ----- | ------- | ---------------------------------- |
| `pvalue_theSame`  | Float | 0.05    | p-value for same copy number test  |
| `fdr_threshold`   | Float | 0.1     | FDR threshold for region extension |
| `highRDthreshold` | Float | 100     | RD difference for FN reduction     |

### Step 7: Visualization

Generate plots for CNV regions.

```julia
# Define sample pairs
pairs = [1 2; 1 3; 2 3; 4 5; 4 6; 5 6]

# Multi-panel plot
result = plotCNVregionMulti(
    metaDataDF,
    CNVscaffold,
    pairedCNV_FN,
    pairs;
    layout = (2, 3),           # 2 rows, 3 columns
    positive_cnv_content = 0.8, # CNV purity threshold
    ylims_pre = [0 0],         # Auto y-axis
    xlims_pre = [0 0],         # Auto x-axis
    chr = "chr20"              # Specific chromosome,ddefault for all chromosomes
)

# Export results
CSV.write("cnv_regions.csv", result[1])
CSV.write("cnv_regions_filtered.csv", result[2])
```

**Plot Parameters:**

| Parameter              | Type   | Default | Description                   |
| ---------------------- | ------ | ------- | ----------------------------- |
| `layout`               | Tuple  | (1,1)   | Subplot layout (rows, cols)   |
| `positive_cnv_content` | Float  | 0.8     | Min positive fragment ratio   |
| `chr`                  | String | "all"   | Chromosome to plot            |
| `xlims_pre`            | Vector | [0 0]   | X-axis limits (auto if [0 0]) |
| `ylims_pre`            | Vector | [0 0]   | Y-axis limits (auto if [0 0]) |

## Input Format

Each input CSV file should contain:

| Column | Name       | Description                        |
d| ------ | ---------- | ---------------------------------- |
| 1      | Chromosome | Chromosome identifier (chr1-chr22) |
| 2      | BinLowEdge | Bin start position                 |
| 3      | BinUpEdge  | Bin end position                   |
| 4      | Content    | Read depth value                   |

**Example:**

```csv
Chromosome,BinLowEdge,BinUpEdge,Content
chr1,0,1000,125
chr1,1000,2000,132
chr1,2000,3000,118
```

## Output Format

### CNV Regions DataFrame

| Column       | Description            |
| ------------ | ---------------------- |
| pair_sampleA | First sample index     |
| pair_sampleB | Second sample index    |
| chr          | Chromosome             |
| start_pos    | CNV start position     |
| end_pos      | CNV end position       |
| length       | CNV region length      |
| mean_depth_A | Mean depth in sample A |
| mean_depth_B | Mean depth in sample B |
| depth_diff   | Depth difference       |

## Parameters

### Recommended Parameter Ranges

| Parameter          | Range               | Description                           |
| ------------------ | ------------------- | ------------------------------------- |
| `depth_filter1`    | [50-200, 2000-5000] | Adjust based on sequencing depth      |
| `depth_filter2`    | [10-20, 200-300]    | Post-normalization filter             |
| `temp_mean`        | Data-dependent      | Expected depth for each copy number   |
| `windowSize`       | 5-20                | Larger = smoother but less resolution |
| `breaks_threshold` | 0.0001-0.01         | Smaller = stricter                    |
| `pvalue_theSame`   | 10⁻⁵ - 10⁻²         | Controls false positives              |
| `fdr_threshold`    | 10⁻⁴ - 10⁻¹         | Controls false negatives              |
| `highRDthreshold`  | 10-50               | RD difference for extension           |

## Examples

### Example 1: Basic CNV Detection

```julia
include("MyFunctions.jl")
include("MyProcess.jl")

# Load data
combineData = process_genomic_data("./data/", ["ctrl.csv", "tumor.csv"])
metaDataDF = filtering_normalizing_seqDepth(combineData, [100, 2500], [15, 250])

# Analyze
metaData = Matrix(metaDataDF[:, 4:end])
esti_param = fitting_gmm(metaData, [35, 65, 95, 125, 165, 195], 10)
CNVscaffold = genome_segment(metaDataDF, esti_param, 0.001, 10, 10)
copynumber_info = copynumber_determination(metaDataDF, CNVscaffold, esti_param, 10)
pairedCNV = copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, esti_param, 10)

# Visualize
plotCNVregionMulti(metaDataDF, CNVscaffold, pairedCNV, [1 2]; chr="all")
```

### Example 2: Specific Chromosome Analysis

```julia
# Plot chr1 with custom range
plotCNVregionMulti(
    metaDataDF, CNVscaffold, pairedCNV, [1 2];
    layout = (1, 1),
    chr = "chr1",
    xlims_pre = [1e7 1.3e7],
    positive_cnv_content = 0.8
)
```

## Contact

For questions and bug reports, please  please send an e-mail to Helab (duananqi@fudan.edu.cn).

---

