#!/bin/bash
# Batch processing script for multiple samples with full/sampled analysis option

# Configuration
BASE_INPUT_DIR="/path/to/your/fastq/directory"
BASE_OUTPUT_DIR="./QC_Results"
SAMPLE_SIZE=100000
THREADS=8
USE_ALL_READS=false  # Set to true for full analysis

# Create output directory
mkdir -p "$BASE_OUTPUT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=========================================="
echo "Batch FASTQ QC Analysis"
echo "=========================================="
echo ""

# Display analysis mode
if [ "$USE_ALL_READS" = true ]; then
    echo -e "${BLUE}Analysis Mode: FULL DATASET (all reads)${NC}"
else
    echo -e "${BLUE}Analysis Mode: SAMPLED ($SAMPLE_SIZE reads per file)${NC}"
fi
echo ""

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastq_qc_enhanced

# Counter for samples
total_samples=0
successful=0
failed=0

# Find all sample directories (assuming S1, S2, S3, etc.)
for sample_dir in "$BASE_INPUT_DIR"/S*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        output_dir="$BASE_OUTPUT_DIR/${sample_name}_QC"
        
        total_samples=$((total_samples + 1))
        
        echo ""
        echo "${YELLOW}Processing: $sample_name${NC}"
        echo "Input:  $sample_dir"
        echo "Output: $output_dir"
        echo ""
        
        # Build command
        cmd="python fastq_qc_analyzer_enhanced.py -i \"$sample_dir\" -o \"$output_dir\" -t $THREADS"
        
        if [ "$USE_ALL_READS" = true ]; then
            cmd="$cmd --all"
            echo -e "${BLUE}Using full dataset analysis${NC}"
        else
            cmd="$cmd -s $SAMPLE_SIZE"
            echo -e "${BLUE}Using sampled analysis ($SAMPLE_SIZE reads)${NC}"
        fi
        
        # Run analysis
        if eval $cmd; then
            echo "${GREEN}✅ Success: $sample_name${NC}"
            successful=$((successful + 1))
        else
            echo "${RED}❌ Failed: $sample_name${NC}"
            failed=$((failed + 1))
        fi
    fi
done

# Summary
echo ""
echo "=========================================="
echo "Batch Analysis Summary"
echo "=========================================="
echo "Total samples processed: $total_samples"
echo -e "${GREEN}Successful: $successful${NC}"
echo -e "${RED}Failed: $failed${NC}"
echo ""
echo "Results saved in: $BASE_OUTPUT_DIR"
echo ""

# Generate combined report
if [ $successful -gt 1 ]; then
    echo "Generating combined comparison report..."
    python generate_comparison_report.py -i "$BASE_OUTPUT_DIR" -o "$BASE_OUTPUT_DIR/Combined_Report"
fi
