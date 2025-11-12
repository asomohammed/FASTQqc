# Enhanced 10x Genomics FASTQ QC Analyzer v2.0

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/yourusername/fastq-qc-analyzer)

A comprehensive, production-ready quality control analysis tool for 10x Genomics single-cell RNA-seq FASTQ files. Provides detailed QC metrics, interactive visualizations, and automated quality assessment with both **sampled** and **full dataset** analysis modes.

---

##  Key Features

###  Dual Analysis Modes
- **Sampled Mode** (Default): Fast QC using 100k reads per file
- **Full Dataset Mode**: Comprehensive analysis of ALL reads with memory optimization

###  10x Genomics-Specific Analysis
- **Barcode Analysis**: 16bp cell barcode QC and diversity metrics
- **UMI Analysis**: 12bp UMI quality and complexity assessment
- **Knee Plot**: Cell calling estimation and barcode rank visualization
- **Library Complexity**: Duplication rate and unique molecule detection

###  Comprehensive QC Metrics
- Quality score distribution and per-base quality
- GC content analysis
- Base composition (A/T/G/C/N)
- Read length distribution
- Adapter contamination detection
- K-mer analysis

###  Advanced Features
- **Pass/Warn/Fail Indicators**: Automated quality assessment
- **Interactive Visualizations**: Plotly-based interactive plots
- **Multi-Sample Comparison**: Batch analysis and comparison reports
- **Automated Recommendations**: Smart suggestions based on QC results
- **Progress Tracking**: Real-time progress bars with tqdm
- **Memory Efficient**: Optimized for large datasets (>100M reads)

---

##  Table of Contents

- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage Guide](#-usage-guide)
- [Analysis Modes](#-analysis-modes)
- [Output Structure](#-output-structure)
- [Quality Thresholds](#-quality-thresholds)
- [Examples](#-examples)
- [Performance](#-performance)
- [Troubleshooting](#-troubleshooting)

---

##  Installation

### Prerequisites
- Python 3.9 or higher
- Conda or Miniconda (recommended)
- 4GB RAM minimum (8GB+ recommended for full analysis)

### Method 1: Automated Setup (Recommended)

```bash
# Download the setup script
wget https://github.com/asomohammed/FASTQqc

# Make it executable
chmod +x setup_enhanced.sh

# Run setup
./setup_enhanced.sh
```

### Method 2: Manual Installation

```bash
# Create conda environment
conda create -n fastq_qc python=3.9 -y

# Activate environment
conda activate fastq_qc

# Install dependencies
conda install -c conda-forge -c bioconda \
    biopython \
    matplotlib \
    seaborn \
    pandas \
    numpy \
    plotly \
    tqdm -y

# Install additional packages
pip install kaleido
```

### Method 3: Using pip

```bash
# Create virtual environment
python3 -m venv fastq_qc_env
source fastq_qc_env/bin/activate

# Install requirements
pip install -r requirements_enhanced.txt
```

### Verify Installation

```bash
conda activate fastq_qc
python fastq_qc.py --help
```

---

## Quick Start

### Basic Usage (Sampled Mode)

```bash
# Activate environment
conda activate fastq_qc

# Run analysis
python fastq_qc.py \
  -i /path/to/fastq/directory \
  -o /path/to/output/directory
```

### Full Dataset Analysis

```bash
python fastq_qc.py \
  -i /path/to/fastq/directory \
  -o /path/to/output/directory \
  --all
```

### View Results

```bash
# Open HTML report in browser
firefox output_directory/qc_report.html

# Or use default browser
xdg-open output_directory/qc_report.html
```

---

##  Usage Guide

### Command-Line Options

```bash
python fastq_qc.py [OPTIONS]

Required Arguments:
  -i, --input DIR          Input directory containing FASTQ files
  -o, --output DIR         Output directory for results

Optional Arguments:
  -s, --sample-size N      Number of reads to sample (default: 100000)
  -t, --threads N          Number of threads to use (default: 4)
  --all, --use-all-reads   Analyze ALL reads instead of sampling
  -h, --help              Show help message and exit
```

### Input File Requirements

Your FASTQ files must follow the standard 10x Genomics naming convention:

 **Correct Format:**
```
SampleName_S1_L001_R1_001.fastq.gz
SampleName_S1_L001_R2_001.fastq.gz
```

 **Also Supported:**
```
SampleName_R1.fastq.gz
SampleName_R2.fastq
Sample_S2_R1_001.fastq.gz
```

**Not Supported:**
```
Sample_1.fastq.gz          # Missing R1/R2 identifier
Sample_read1.fastq.gz      # Wrong identifier format
Sample.fastq.gz            # No read type specified
```

### File Organization

```
your_fastq_directory/
├── Sample1_S1_L001_R1_001.fastq.gz
├── Sample1_S1_L001_R2_001.fastq.gz
├── Sample2_S2_L001_R1_001.fastq.gz
├── Sample2_S2_L001_R2_001.fastq.gz
└── ...
```

---

##  Analysis Modes

### Sampled Mode (Default)

**Best for:** Quick QC checks, large datasets, routine analysis

```bash
python fastq_qc.py -i INPUT -o OUTPUT -s 100000
```

**Features:**
- Fast analysis (minutes)
- Low memory usage
- Representative statistics
- Suitable for >100M reads

**What's Analyzed:**
- Quality metrics: 100k reads sampled
- Barcode/UMI: All reads counted
- Base composition: All reads counted

### Full Dataset Mode

**Best for:** Final QC, publication, accurate cell counts

```bash
python fastq_qc.py -i INPUT -o OUTPUT --all
```

**Features:**
-  Comprehensive analysis
-  Maximum accuracy
-  Complete statistics
-  Memory-optimized for large files

**What's Analyzed:**
- Quality metrics: All reads
- Barcode/UMI: All reads
- Base composition: All reads
- Duplication: All reads

### Comparison Table

| Feature | Sampled Mode | Full Mode |
|---------|--------------|-----------|
| **Speed** | ⚡⚡⚡ Fast |  Slower |
| **Memory** |  Low (2-4GB) |  Higher (4-16GB) |
| **Accuracy** | ✓ Good | ✓✓✓ Excellent |
| **Best For** | Quick checks | Final analysis |
| **1M reads** | ~2 min | ~5 min |
| **10M reads** | ~2 min | ~30 min |
| **100M reads** | ~2 min | ~4 hours |

---

##  Output Structure

```
output_directory/
│
├── qc_report.html                      #  Main HTML report (open this!)
│
├── plots/                              #  Static PNG plots
│   ├── Sample_R1_quality_dist.png
│   ├── Sample_R1_quality_per_pos.png
│   ├── Sample_R1_gc_content.png
│   ├── Sample_R1_base_comp.png
│   ├── Sample_R1_read_length.png
│   ├── Sample_barcode_analysis.png
│   ├── Sample_R2_*.png
│   └── ...
│
├── interactive_plots/                  #  Interactive HTML plots
│   ├── Sample_R1_interactive.html
│   ├── Sample_R2_interactive.html
│   ├── Sample_knee_plot.html
│   ├── Sample_umi_analysis.html
│   └── multi_sample_comparison.html    # (if multiple samples)
│
├── statistics/                         #  JSON statistics
│   ├── Sample1_statistics.json
│   ├── Sample2_statistics.json
│   └── ...
│
└── reports/                            #  CSV reports
    └── comparison_table.csv            # (if multiple samples)
```

### Key Output Files

#### 1. Main HTML Report (`qc_report.html`)
- **Overview**: Summary of all samples
- **Quality Assessment**: Pass/Warn/Fail indicators
- **Detailed Metrics**: Comprehensive statistics tables
- **Visualizations**: All plots in one place
- **Recommendations**: Automated quality suggestions
- **Interactive Elements**: Tabbed interface for easy navigation

#### 2. Interactive Plots
- **Zoomable**: Click and drag to zoom
- **Hover Information**: Detailed values on hover
- **Downloadable**: Save as PNG from browser
- **Responsive**: Works on all screen sizes

#### 3. Statistics JSON
```json
{
  "R1": {
    "sample": "Sample1",
    "total_reads": 1000000,
    "analyzed_reads": 100000,
    "is_full_analysis": false,
    "mean_quality": 35.2,
    "q30_percentage": 89.5,
    "unique_barcodes": 8543,
    "barcode_diversity": 0.8543,
    "mean_umi_quality": 34.8
  },
  "R2": { ... }
}
```

---

##  Quality Thresholds

### Automated Quality Assessment

Each metric is automatically assessed as **PASS**, **WARN**, or **FAIL**:

| Metric | Pass (✅) | Warn (⚠️) | Fail (❌) |
|--------|----------|----------|----------|
| **Q30 Percentage** | ≥ 75% | 65-75% | < 65% |
| **Mean Quality** | ≥ 30 | 25-30 | < 25 |
| **GC Content** | 40-60% | Outside range | - |
| **Duplication Rate** | < 20% | 20-40% | > 40% |

### Quality Indicators

🟢 **PASS**: All metrics meet quality standards
- Proceed with downstream analysis
- Data quality is excellent

🟡 **WARN**: Some metrics need attention
- Review specific warnings
- May proceed with caution
- Consider re-sequencing if critical

🔴 **FAIL**: Critical quality issues detected
- Do not proceed with analysis
- Re-sequencing recommended
- Check library preparation

---

##  Examples

### Example 1: Single Sample Quick Check

```bash
# Fast QC check of one sample
python fastq_qc.py \
  -i ./fastq_data/Sample1 \
  -o ./qc_results/Sample1_QC

# View results
firefox ./qc_results/Sample1_QC/qc_report.html
```

### Example 2: Full Analysis with Multiple Threads

```bash
# Comprehensive analysis using 8 threads
python fastq_qc.py \
  -i ./fastq_data/Sample1 \
  -o ./qc_results/Sample1_Full_QC \
  --all \
  -t 8
```

### Example 3: Custom Sample Size

```bash
# Analyze 200k reads per file
python fastq_qc.py \
  -i ./fastq_data/Sample1 \
  -o ./qc_results/Sample1_QC \
  -s 200000 \
  -t 4
```

### Example 4: Batch Processing Multiple Samples

```bash
# Process all samples in a directory
for sample_dir in ./fastq_data/S*; do
  sample_name=$(basename $sample_dir)
  echo "Processing $sample_name..."
  
  python fastq_qc.py \
    -i "$sample_dir" \
    -o "./qc_results/${sample_name}_QC" \
    -t 8
done

# Generate comparison report
python generate_comparison_report.py \
  -i ./qc_results \
  -o ./qc_results/Combined_Report
```

### Example 5: Using Batch Script

```bash
# Edit batch_analyze.sh
nano batch_analyze.sh

# Set your parameters:
# BASE_INPUT_DIR="/path/to/fastq"
# BASE_OUTPUT_DIR="./QC_Results"
# USE_ALL_READS=false  # or true for full analysis
# THREADS=8

# Run batch analysis
chmod +x batch_analyze.sh
./batch_analyze.sh
```

### Example 6: High-Throughput Pipeline

```bash
#!/bin/bash
# Pipeline for multiple experiments

EXPERIMENTS=("Exp1" "Exp2" "Exp3")
BASE_DIR="/data/sequencing"

for exp in "${EXPERIMENTS[@]}"; do
  echo "=== Processing $exp ==="
  
  # Run QC
  python fastq_qc.py \
    -i "$BASE_DIR/$exp/fastq" \
    -o "$BASE_DIR/$exp/qc" \
    --all \
    -t 16
  
  # Generate report
  echo "QC completed for $exp"
done

# Generate combined comparison
python generate_comparison_report.py \
  -i "$BASE_DIR" \
  -o "$BASE_DIR/Combined_QC_Report"
```

---

## ⚡ Performance

### Benchmark Results

Tested on: Intel Xeon E5-2680 v4 @ 2.40GHz, 64GB RAM, SSD storage

| Dataset Size | Sampled Mode | Full Mode | Memory Usage |
|--------------|--------------|-----------|--------------|
| 1M reads | 1-2 min | 4-6 min | 2-3 GB |
| 5M reads | 1-2 min | 15-20 min | 3-4 GB |
| 10M reads | 2-3 min | 30-40 min | 4-6 GB |
| 50M reads | 2-3 min | 2-3 hours | 6-10 GB |
| 100M reads | 2-3 min | 4-5 hours | 8-16 GB |

### Performance Tips

1. **Use Sampled Mode** for routine QC
   ```bash
   python fastq_qc.py -i INPUT -o OUTPUT
   ```

2. **Increase Threads** for faster processing
   ```bash
   python fastq_qc.py -i INPUT -o OUTPUT -t 16
   ```

3. **Use SSD Storage** for I/O-intensive operations

4. **Close Other Applications** when running full analysis

5. **Process Samples in Parallel** for multiple samples
   ```bash
   # Use GNU parallel
   parallel -j 4 python fastq_qc.py -i {} -o {}_QC ::: Sample*
   ```

---

##  Troubleshooting

### Common Issues and Solutions

#### 1. "No FASTQ files found"

**Problem:** Script can't find your FASTQ files

**Solutions:**
```bash
# Check file naming
ls -lh /path/to/fastq/

# Files must contain _R1_ or _R2_
# Correct: Sample_S1_R1_001.fastq.gz
# Wrong: Sample_1.fastq.gz

# Check file extensions
# Supported: .fastq, .fastq.gz, .fq, .fq.gz
```

#### 2. Out of Memory Error

**Problem:** System runs out of RAM during analysis

**Solutions:**
```bash
# Solution 1: Use sampled mode
python fastq_qc.py -i INPUT -o OUTPUT

# Solution 2: Reduce sample size
python fastq_qc.py -i INPUT -o OUTPUT -s 50000

# Solution 3: Close other applications
# Solution 4: Use a machine with more RAM
```

#### 3. Analysis Too Slow

**Problem:** Full analysis taking too long

**Solutions:**
```bash
# Solution 1: Use sampled mode
python fastq_qc.py -i INPUT -o OUTPUT

# Solution 2: Increase threads
python fastq_qc.py -i INPUT -o OUTPUT -t 16

# Solution 3: Use faster storage (SSD)
# Solution 4: Reduce sample size
python fastq_qc.py -i INPUT -o OUTPUT -s 50000
```

#### 4. Import Errors

**Problem:** Missing Python packages

**Solutions:**
```bash
# Reinstall dependencies
conda activate fastq_qc
conda install -c conda-forge -c bioconda biopython matplotlib seaborn pandas numpy plotly tqdm -y
pip install kaleido

# Or recreate environment
conda env remove -n fastq_qc
./setup_enhanced.sh
```

#### 5. Permission Denied

**Problem:** Can't write to output directory

**Solutions:**
```bash
# Check permissions
ls -ld /path/to/output

# Create directory with correct permissions
mkdir -p /path/to/output
chmod 755 /path/to/output

# Or use a different output directory
python fastq_qc.py -i INPUT -o ~/qc_results
```

#### 6. Corrupted FASTQ Files

**Problem:** Error reading FASTQ files

**Solutions:**
```bash
# Check file integrity
gunzip -t file.fastq.gz

# Check file format
zcat file.fastq.gz | head -n 4

# Re-download or regenerate files if corrupted
```

---


##  Technical Details

### Algorithm Overview

1. **File Discovery**: Automatically pairs R1/R2 files
2. **Read Sampling**: Intelligent sampling for large files
3. **Quality Analysis**: Per-base and overall quality metrics
4. **Barcode/UMI Extraction**: 10x-specific sequence analysis
5. **Statistical Analysis**: Comprehensive metric calculation
6. **Visualization**: Static and interactive plot generation
7. **Report Generation**: HTML report with all results

### Memory Optimization

For large files (>10M reads), the script uses:
- Streaming analysis to reduce memory footprint
- Intelligent sampling for quality metrics (every Nth read)
- Full counting for barcode/UMI statistics
- Efficient data structures (Counters, defaultdicts)

### Quality Control Pipeline

```
FASTQ Files
    ↓
File Validation
    ↓
Read Sampling/Streaming
    ↓
Quality Metrics ← Base Composition
    ↓              ↓
GC Content → Statistical Analysis
    ↓              ↓
Barcode/UMI → Diversity Metrics
    ↓              ↓
Duplication → Assessment
    ↓
Visualization
    ↓
HTML Report
```

---

