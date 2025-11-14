#!/usr/bin/env python3
"""
FASTQ QC Analysis Pipeline
"""

import os
import sys
import gzip
import argparse
from pathlib import Path
from collections import defaultdict, Counter
import json
from datetime import datetime
import multiprocessing as mp
from functools import partial
import Carefulings
warnings.filterwarnings('ignore')

# Progress bars
from tqdm import tqdm

# Data processing
import numpy as np
import pandas as pd
from Bio import SeqIO

# Plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

class EnhancedFastqAnalyzer:
    def __init__(self, fastq_dir, output_dir, sample_size=None, threads=4, use_all_reads=False):
        self.fastq_dir = Path(fastq_dir)
        self.output_dir = Path(output_dir)
        self.sample_size = sample_size
        self.threads = threads
        self.use_all_reads = use_all_reads
        self.results = {}
        
        # 10x Genomics specifications
        self.barcode_length = 16
        self.umi_length = 12
        self.umi_start = 16
        self.umi_end = 28
        
        # Quality thresholds
        self.thresholds = {
            'q30_pass': 75,
            'q30_Careful': 65,
            'mean_quality_pass': 30,
            'mean_quality_Careful': 25,
            'gc_content_min': 40,
            'gc_content_max': 60
        }
        
        # Create output directories
        self.plots_dir = self.output_dir / "plots"
        self.interactive_dir = self.output_dir / "interactive_plots"
        self.stats_dir = self.output_dir / "statistics"
        self.reports_dir = self.output_dir / "reports"
        
        for dir_path in [self.plots_dir, self.interactive_dir, 
                         self.stats_dir, self.reports_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def find_fastq_files(self):
        """Find and pair R1 and R2 FASTQ files"""
        fastq_files = defaultdict(dict)
        
        for file in self.fastq_dir.glob("*.fastq*"):
            if "_R1_" in file.name or "_R1." in file.name:
                sample_name = file.name.split("_R1")[0]
                fastq_files[sample_name]['R1'] = file
            elif "_R2_" in file.name or "_R2." in file.name:
                sample_name = file.name.split("_R2")[0]
                fastq_files[sample_name]['R2'] = file
        
        return fastq_files
    
    def open_fastq(self, filepath):
        """Open FASTQ file (gzipped or not)"""
        if str(filepath).endswith('.gz'):
            return gzip.open(filepath, 'rt')
        else:
            return open(filepath, 'r')
    
    def count_reads_in_fastq(self, fastq_path):
        """Count total number of reads in FASTQ file"""
        print(f"    Counting reads in {fastq_path.name}...")
        count = 0
        with self.open_fastq(fastq_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                count += 1
        return count
    
    def analyze_fastq_enhanced(self, fastq_path, read_type):
        """ FASTQ analysis"""
        
        # Determine how many reads to analyze
        if self.use_all_reads:
            total_reads = self.count_reads_in_fastq(fastq_path)
            reads_to_analyze = total_reads
            print(f"  📊 Analyzing ALL {total_reads:,} reads from {read_type}: {fastq_path.name}")
        else:
            reads_to_analyze = self.sample_size
            print(f"  📊 Analyzing {reads_to_analyze:,} sampled reads from {read_type}: {fastq_path.name}")
        
        stats = {
            'total_reads': 0,
            'total_bases': 0,
            'quality_scores': [],
            'read_lengths': [],
            'gc_content': [],
            'base_composition': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0},
            'quality_per_position': defaultdict(list),
            'quality_distribution': defaultdict(int),
            'barcode_data': {} if read_type == 'R1' else None,
            'umi_data': {} if read_type == 'R1' else None,
            'duplication_estimate': {},
            'adapter_content': defaultdict(int),
            'kmer_content': Counter(),
            'analyzed_reads': reads_to_analyze,
            'is_full_analysis': self.use_all_reads
        }
        
        sequences = []
        
        # For memory efficiency with large files, we'll use streaming approach
        if self.use_all_reads and reads_to_analyze > 1000000:
            # Use sampling for quality metrics but count everything for other stats
            stats = self._analyze_large_fastq(fastq_path, read_type, reads_to_analyze)
        else:
            # Standard analysis
            with self.open_fastq(fastq_path) as handle:
                for i, record in enumerate(tqdm(SeqIO.parse(handle, "fastq"), 
                                               total=reads_to_analyze,
                                               desc=f"    Processing {read_type}",
                                               disable=False)):
                    if not self.use_all_reads and i >= reads_to_analyze:
                        break
                    
                    stats['total_reads'] += 1
                    seq = str(record.seq)
                    stats['total_bases'] += len(seq)
                    stats['read_lengths'].append(len(seq))
                    sequences.append(seq)
                    
                    # Quality scores
                    qualities = record.letter_annotations["phred_quality"]
                    stats['quality_scores'].extend(qualities)
                    
                    # Quality distribution
                    for q in qualities:
                        stats['quality_distribution'][q] += 1
                    
                    # Quality per position
                    for pos, qual in enumerate(qualities):
                        stats['quality_per_position'][pos].append(qual)
                    
                    # GC content
                    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                    stats['gc_content'].append(gc)
                    
                    # Base composition
                    for base in seq:
                        if base in stats['base_composition']:
                            stats['base_composition'][base] += 1
                    
                    # K-mer analysis (5-mers) - sample only for large datasets
                    if stats['total_reads'] <= 100000 or stats['total_reads'] % 10 == 0:
                        for j in range(len(seq) - 4):
                            kmer = seq[j:j+5]
                            stats['kmer_content'][kmer] += 1
                    
                    # R1-specific analysis
                    if read_type == 'R1':
                        # Barcode analysis (first 16bp)
                        if len(seq) >= self.barcode_length:
                            barcode = seq[:self.barcode_length]
                            if 'barcodes' not in stats['barcode_data']:
                                stats['barcode_data']['barcodes'] = []
                                stats['barcode_data']['barcode_quality'] = []
                            stats['barcode_data']['barcodes'].append(barcode)
                            bc_qual = np.mean(qualities[:self.barcode_length])
                            stats['barcode_data']['barcode_quality'].append(bc_qual)
                        
                        # UMI analysis (bases 17-28 for 10x v3)
                        if len(seq) >= self.umi_end:
                            umi = seq[self.umi_start:self.umi_end]
                            if 'umis' not in stats['umi_data']:
                                stats['umi_data']['umis'] = []
                                stats['umi_data']['umi_quality'] = []
                            stats['umi_data']['umis'].append(umi)
                            umi_qual = np.mean(qualities[self.umi_start:self.umi_end])
                            stats['umi_data']['umi_quality'].append(umi_qual)
        
        # Duplication analysis
        seq_counts = Counter(sequences)
        stats['duplication_estimate']['unique_reads'] = len(seq_counts)
        stats['duplication_estimate']['total_reads'] = len(sequences)
        stats['duplication_estimate']['duplication_rate'] = \
            1 - (len(seq_counts) / len(sequences)) if sequences else 0
        
        return stats
    
    def _analyze_large_fastq(self, fastq_path, read_type, total_reads):
        """Memory-efficient analysis for very large FASTQ files"""
        print(f"    Using memory-efficient mode for large file...")
        
        stats = {
            'total_reads': 0,
            'total_bases': 0,
            'quality_scores': [],
            'read_lengths': [],
            'gc_content': [],
            'base_composition': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0},
            'quality_per_position': defaultdict(list),
            'quality_distribution': defaultdict(int),
            'barcode_data': {} if read_type == 'R1' else None,
            'umi_data': {} if read_type == 'R1' else None,
            'duplication_estimate': {},
            'adapter_content': defaultdict(int),
            'kmer_content': Counter(),
            'analyzed_reads': total_reads,
            'is_full_analysis': True
        }
        
        # Sample every Nth read for detailed quality analysis
        sample_interval = max(1, total_reads // 100000)  # Sample ~100k reads
        sequences_sample = []
        
        with self.open_fastq(fastq_path) as handle:
            for i, record in enumerate(tqdm(SeqIO.parse(handle, "fastq"), 
                                           total=total_reads,
                                           desc=f"    Processing {read_type} (efficient mode)",
                                           disable=False)):
                
                stats['total_reads'] += 1
                seq = str(record.seq)
                stats['total_bases'] += len(seq)
                stats['read_lengths'].append(len(seq))
                
                # Sample for detailed analysis
                if i % sample_interval == 0:
                    sequences_sample.append(seq)
                    qualities = record.letter_annotations["phred_quality"]
                    stats['quality_scores'].extend(qualities)
                    
                    for q in qualities:
                        stats['quality_distribution'][q] += 1
                    
                    for pos, qual in enumerate(qualities):
                        stats['quality_per_position'][pos].append(qual)
                    
                    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                    stats['gc_content'].append(gc)
                
                # Always count base composition
                for base in seq:
                    if base in stats['base_composition']:
                        stats['base_composition'][base] += 1
                
                # R1-specific analysis
                if read_type == 'R1':
                    if len(seq) >= self.barcode_length:
                        barcode = seq[:self.barcode_length]
                        if 'barcodes' not in stats['barcode_data']:
                            stats['barcode_data']['barcodes'] = []
                            stats['barcode_data']['barcode_quality'] = []
                        stats['barcode_data']['barcodes'].append(barcode)
                        
                        if i % sample_interval == 0:
                            qualities = record.letter_annotations["phred_quality"]
                            bc_qual = np.mean(qualities[:self.barcode_length])
                            stats['barcode_data']['barcode_quality'].append(bc_qual)
                    
                    if len(seq) >= self.umi_end:
                        umi = seq[self.umi_start:self.umi_end]
                        if 'umis' not in stats['umi_data']:
                            stats['umi_data']['umis'] = []
                            stats['umi_data']['umi_quality'] = []
                        stats['umi_data']['umis'].append(umi)
                        
                        if i % sample_interval == 0:
                            qualities = record.letter_annotations["phred_quality"]
                            umi_qual = np.mean(qualities[self.umi_start:self.umi_end])
                            stats['umi_data']['umi_quality'].append(umi_qual)
        
        # Duplication analysis on sample
        seq_counts = Counter(sequences_sample)
        stats['duplication_estimate']['unique_reads'] = len(seq_counts)
        stats['duplication_estimate']['total_reads'] = len(sequences_sample)
        stats['duplication_estimate']['duplication_rate'] = \
            1 - (len(seq_counts) / len(sequences_sample)) if sequences_sample else 0
        stats['duplication_estimate']['note'] = f"Estimated from {len(sequences_sample):,} sampled reads"
        
        return stats
    
    def assess_quality(self, stats, read_type):
        """Assess quality and return pass/Careful/warn status"""
        assessment = {
            'overall': 'PASS',
            'metrics': {}
        }
        
        # Q30 percentage
        q30_bases = sum(1 for q in stats['quality_scores'] if q >= 30)
        q30_pct = (q30_bases / len(stats['quality_scores']) * 100) if stats['quality_scores'] else 0
        
        if q30_pct >= self.thresholds['q30_pass']:
            assessment['metrics']['q30'] = 'PASS'
        elif q30_pct >= self.thresholds['q30_Careful']:
            assessment['metrics']['q30'] = 'Careful'
            assessment['overall'] = 'Careful'
        else:
            assessment['metrics']['q30'] = 'warn'
            assessment['overall'] = 'warn'
        
        # Mean quality
        mean_qual = np.mean(stats['quality_scores'])
        if mean_qual >= self.thresholds['mean_quality_pass']:
            assessment['metrics']['mean_quality'] = 'PASS'
        elif mean_qual >= self.thresholds['mean_quality_Careful']:
            assessment['metrics']['mean_quality'] = 'Careful'
            if assessment['overall'] == 'PASS':
                assessment['overall'] = 'Careful'
        else:
            assessment['metrics']['mean_quality'] = 'warn'
            assessment['overall'] = 'warn'
        
        # GC content
        mean_gc = np.mean(stats['gc_content'])
        if self.thresholds['gc_content_min'] <= mean_gc <= self.thresholds['gc_content_max']:
            assessment['metrics']['gc_content'] = 'PASS'
        else:
            assessment['metrics']['gc_content'] = 'Careful'
            if assessment['overall'] == 'PASS':
                assessment['overall'] = 'Careful'
        
        # Duplication rate
        dup_rate = stats['duplication_estimate']['duplication_rate']
        if dup_rate < 0.2:
            assessment['metrics']['duplication'] = 'PASS'
        elif dup_rate < 0.4:
            assessment['metrics']['duplication'] = 'Careful'
            if assessment['overall'] == 'PASS':
                assessment['overall'] = 'Careful'
        else:
            assessment['metrics']['duplication'] = 'warn'
            assessment['overall'] = 'warn'
        
        return assessment
    
    def create_knee_plot(self, stats, sample_name):
        """Create knee plot for cell calling estimation"""
        if 'barcode_data' not in stats or not stats['barcode_data']:
            return None
        
        barcodes = stats['barcode_data']['barcodes']
        barcode_counts = Counter(barcodes)
        
        # Sort by count (descending)
        sorted_counts = sorted(barcode_counts.values(), reverse=True)
        ranks = list(range(1, len(sorted_counts) + 1))
        
        # Create interactive plot
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=ranks,
            y=sorted_counts,
            mode='lines',
            name='Barcode Counts',
            line=dict(color='steelblue', width=2)
        ))
        
        # Estimate knee point (simple method)
        if len(sorted_counts) > 100:
            # Find inflection point
            log_counts = np.log10(sorted_counts[:min(10000, len(sorted_counts))])
            log_ranks = np.log10(ranks[:len(log_counts)])
            
            # Simple derivative method
            derivatives = np.diff(log_counts) / np.diff(log_ranks)
            knee_idx = np.argmax(derivatives < np.percentile(derivatives, 10))
            
            fig.add_vline(x=knee_idx, line_dash="dash", line_color="red",
                         annotation_text=f"Estimated cells: ~{knee_idx}")
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        
        fig.update_layout(
            title=f'{sample_name}: Barcode Rank Plot (Knee Plot)<br><sub>{analysis_note}</sub>',
            xaxis_title='Barcode Rank',
            yaxis_title='UMI Counts',
            xaxis_type='log',
            yaxis_type='log',
            hovermode='closest',
            template='plotly_white'
        )
        
        output_path = self.interactive_dir / f"{sample_name}_knee_plot.html"
        fig.write_html(str(output_path))
        
        return output_path.name
    
    def create_interactive_quality_plot(self, stats, sample_name, read_type):
        """Create interactive quality distribution plot"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Quality Score Distribution', 
                          'Quality per Position',
                          'GC Content Distribution',
                          'Read Length Distribution')
        )
        
        # Quality distribution
        fig.add_trace(
            go.Histogram(x=stats['quality_scores'], name='Quality Scores',
                        marker_color='steelblue', nbinsx=50),
            row=1, col=1
        )
        
        # Quality per position
        positions = sorted(stats['quality_per_position'].keys())
        mean_qual = [np.mean(stats['quality_per_position'][pos]) for pos in positions]
        q25 = [np.percentile(stats['quality_per_position'][pos], 25) for pos in positions]
        q75 = [np.percentile(stats['quality_per_position'][pos], 75) for pos in positions]
        
        fig.add_trace(
            go.Scatter(x=positions, y=mean_qual, name='Mean Quality',
                      line=dict(color='blue')),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=positions, y=q75, fill=None, mode='lines',
                      line_color='lightblue', showlegend=False),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=positions, y=q25, fill='tonexty', mode='lines',
                      line_color='lightblue', name='25-75 percentile'),
            row=1, col=2
        )
        
        # GC content
        fig.add_trace(
            go.Histogram(x=stats['gc_content'], name='GC Content',
                        marker_color='green', nbinsx=50),
            row=2, col=1
        )
        
        # Read length
        fig.add_trace(
            go.Histogram(x=stats['read_lengths'], name='Read Length',
                        marker_color='orange', nbinsx=30),
            row=2, col=2
        )
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        
        fig.update_layout(
            title_text=f'{sample_name} - {read_type}: Quality Metrics<br><sub>{analysis_note}</sub>',
            showlegend=True,
            height=800,
            template='plotly_white'
        )
        
        output_path = self.interactive_dir / f"{sample_name}_{read_type}_interactive.html"
        fig.write_html(str(output_path))
        
        return output_path.name
    
    def create_umi_analysis_plot(self, stats, sample_name):
        """Create UMI analysis plots"""
        if 'umi_data' not in stats or not stats['umi_data']:
            return None
        
        umis = stats['umi_data']['umis']
        umi_counts = Counter(umis)
        
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('UMI Diversity', 'UMI Quality Distribution')
        )
        
        # UMI diversity
        unique_umis = len(umi_counts)
        total_umis = len(umis)
        diversity = unique_umis / total_umis
        
        fig.add_trace(
            go.Bar(x=['Total UMIs', 'Unique UMIs'],
                  y=[total_umis, unique_umis],
                  marker_color=['steelblue', 'orange']),
            row=1, col=1
        )
        
        # UMI quality
        fig.add_trace(
            go.Histogram(x=stats['umi_data']['umi_quality'],
                        name='UMI Quality',
                        marker_color='purple', nbinsx=30),
            row=1, col=2
        )
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled reads"
        
        fig.update_layout(
            title_text=f'{sample_name}: UMI Analysis (Diversity: {diversity:.4f})<br><sub>{analysis_note}</sub>',
            showlegend=False,
            height=400,
            template='plotly_white'
        )
        
        output_path = self.interactive_dir / f"{sample_name}_umi_analysis.html"
        fig.write_html(str(output_path))
        
        return output_path.name
    
    def generate_summary_stats(self, stats, sample_name, read_type):
        """Generate comprehensive summary statistics"""
        q30_bases = sum(1 for q in stats['quality_scores'] if q >= 30)
        q30_percentage = (q30_bases / len(stats['quality_scores']) * 100) if stats['quality_scores'] else 0
        
        summary = {
            'sample': sample_name,
            'read_type': read_type,
            'total_reads': stats['total_reads'],
            'analyzed_reads': stats.get('analyzed_reads', stats['total_reads']),
            'is_full_analysis': stats.get('is_full_analysis', False),
            'total_bases': stats['total_bases'],
            'mean_read_length': float(np.mean(stats['read_lengths'])),
            'median_read_length': float(np.median(stats['read_lengths'])),
            'mean_quality': float(np.mean(stats['quality_scores'])),
            'median_quality': float(np.median(stats['quality_scores'])),
            'q30_percentage': float(q30_percentage),
            'q20_percentage': float(sum(1 for q in stats['quality_scores'] if q >= 20) / len(stats['quality_scores']) * 100),
            'mean_gc_content': float(np.mean(stats['gc_content'])),
            'base_composition': stats['base_composition'],
            'duplication_rate': float(stats['duplication_estimate']['duplication_rate']),
            'unique_reads': stats['duplication_estimate']['unique_reads']
        }
        
        # Add R1-specific metrics
        if read_type == 'R1' and stats['barcode_data']:
            barcodes = stats['barcode_data']['barcodes']
            barcode_counts = Counter(barcodes)
            summary['unique_barcodes'] = len(barcode_counts)
            summary['barcode_diversity'] = len(barcode_counts) / len(barcodes)
            summary['mean_barcode_quality'] = float(np.mean(stats['barcode_data']['barcode_quality']))
            
            if stats['umi_data']:
                umis = stats['umi_data']['umis']
                umi_counts = Counter(umis)
                summary['unique_umis'] = len(umi_counts)
                summary['umi_diversity'] = len(umi_counts) / len(umis)
                summary['mean_umi_quality'] = float(np.mean(stats['umi_data']['umi_quality']))
        
        return summary
    
    def create_comparison_report(self, all_results):
        """Create multi-sample comparison report"""
        if len(all_results) < 2:
            return None
        
        # Prepare comparison data
        comparison_data = []
        for sample_name, sample_data in all_results.items():
            for read_type in ['R1', 'R2']:
                if read_type in sample_data['summary']:
                    row = sample_data['summary'][read_type].copy()
                    row['sample'] = sample_name
                    row['read_type'] = read_type
                    comparison_data.append(row)
        
        df = pd.DataFrame(comparison_data)
        
        # Create comparison plots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Mean Quality Comparison',
                          'Q30 Percentage Comparison',
                          'GC Content Comparison',
                          'Duplication Rate Comparison')
        )
        
        # Group by sample and read type
        for read_type in ['R1', 'R2']:
            df_read = df[df['read_type'] == read_type]
            
            # Mean quality
            fig.add_trace(
                go.Bar(name=read_type, x=df_read['sample'], 
                      y=df_read['mean_quality']),
                row=1, col=1
            )
            
            # Q30 percentage
            fig.add_trace(
                go.Bar(name=read_type, x=df_read['sample'],
                      y=df_read['q30_percentage']),
                row=1, col=2
            )
            
            # GC content
            fig.add_trace(
                go.Bar(name=read_type, x=df_read['sample'],
                      y=df_read['mean_gc_content']),
                row=2, col=1
            )
            
            # Duplication rate
            fig.add_trace(
                go.Bar(name=read_type, x=df_read['sample'],
                      y=df_read['duplication_rate']),
                row=2, col=2
            )
        
        fig.update_layout(
            title_text='Multi-Sample Comparison',
            showlegend=True,
            height=800,
            template='plotly_white'
        )
        
        output_path = self.interactive_dir / "multi_sample_comparison.html"
        fig.write_html(str(output_path))
        
        # Save comparison table
        csv_path = self.reports_dir / "comparison_table.csv"
        df.to_csv(csv_path, index=False)
        
        return output_path.name, csv_path.name
    
    def create_enhanced_html_report(self, all_results):
        """Create  HTML reports"""
        
        analysis_mode = "FULL DATASET" if self.use_all_reads else f"SAMPLED ({self.sample_size:,} reads)"
        
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title> Genomics FASTQ QC Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            border-radius: 15px;
        }}
        
        .header {{
            text-align: center;
            padding: 30px 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .analysis-mode {{
            display: inline-block;
            padding: 10px 25px;
            background-color: rgba(255,255,255,0.2);
            border-radius: 20px;
            margin-top: 10px;
            font-weight: bold;
            font-size: 1.1em;
        }}
        
        .status-badge {{
            display: inline-block;
            padding: 8px 20px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 14px;
            margin: 5px;
        }}
        
        .status-pass {{
            background-color: #2ecc71;
            color: white;
        }}
        
        .status-Careful {{
            background-color: #f39c12;
            color: white;
        }}
        
        .status-warn {{
            background-color: #e74c3c;
            color: white;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }}
        
        .summary-card:hover {{
            transform: translateY(-5px);
        }}
        
        .summary-card h3 {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 10px;
        }}
        
        .summary-card .value {{
            font-size: 32px;
            font-weight: bold;
        }}
        
        .sample-section {{
            margin: 40px 0;
            padding: 30px;
            background-color: #f8f9fa;
            border-radius: 10px;
            border-left: 5px solid #667eea;
        }}
        
        .sample-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
        }}
        
        .sample-header h2 {{
            color: #2c3e50;
            font-size: 1.8em;
        }}
        
        .stats-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background-color: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .stats-table th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        
        .stats-table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #ecf0f1;
        }}
        
        .stats-table tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 25px;
            margin: 30px 0;
        }}
        
        .plot-item {{
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }}
        
        .plot-item:hover {{
            transform: translateY(-5px);
            box-shadow: 0 5px 20px rgba(0,0,0,0.15);
        }}
        
        .plot-item img {{
            width: 100%;
            height: auto;
            border-radius: 5px;
        }}
        
        .plot-item h3 {{
            margin-top: 15px;
            color: #2c3e50;
            font-size: 16px;
            text-align: center;
        }}
        
        .interactive-link {{
            display: inline-block;
            margin: 10px 0;
            padding: 10px 20px;
            background-color: #3498db;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            transition: background-color 0.3s;
        }}
        
        .interactive-link:hover {{
            background-color: #2980b9;
        }}
        
        .alert {{
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            border-left: 5px solid;
        }}
        
        .alert-info {{
            background-color: #d1ecf1;
            border-color: #0c5460;
            color: #0c5460;
        }}
        
        .alert-success {{
            background-color: #d4edda;
            border-color: #155724;
            color: #155724;
        }}
        
        .alert-warning {{
            background-color: #fff3cd;
            border-color: #856404;
            color: #856404;
        }}
        
        .alert-danger {{
            background-color: #f8d7da;
            border-color: #721c24;
            color: #721c24;
        }}
        
        .recommendations {{
            background-color: #fff3cd;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            border-left: 5px solid #ffc107;
        }}
        
        .recommendations h3 {{
            color: #856404;
            margin-bottom: 15px;
        }}
        
        .recommendations ul {{
            list-style-position: inside;
            color: #856404;
        }}
        
        .recommendations li {{
            margin: 8px 0;
        }}
        
        .footer {{
            margin-top: 50px;
            padding: 20px;
            text-align: center;
            color: #7f8c8d;
            border-top: 2px solid #ecf0f1;
        }}
        
        .tabs {{
            display: flex;
            border-bottom: 2px solid #667eea;
            margin: 20px 0;
        }}
        
        .tab {{
            padding: 15px 30px;
            cursor: pointer;
            background-color: #f8f9fa;
            border: none;
            border-bottom: 3px solid transparent;
            transition: all 0.3s;
            font-size: 16px;
            font-weight: 500;
        }}
        
        .tab:hover {{
            background-color: #e9ecef;
        }}
        
        .tab.active {{
            background-color: white;
            border-bottom-color: #667eea;
            color: #667eea;
        }}
        
        .tab-content {{
            display: none;
            padding: 20px;
            animation: fadeIn 0.5s;
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; }}
            to {{ opacity: 1; }}
        }}
        
        @media print {{
            body {{
                background: white;
            }}
            .container {{
                box-shadow: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Enhanced 10x Genomics FASTQ QC Report</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <div class="analysis-mode">📊 Analysis Mode: {analysis_mode}</div>
        </div>
        
        <div class="alert alert-info">
            <strong>📊 Analysis Summary:</strong> Analyzed {len(all_results)} sample(s) with dual-read FASTQ files.
            {'<br><strong>Mode:</strong> Full dataset analysis - all reads were processed.' if self.use_all_reads else f'<br><strong>Mode:</strong> Sampled analysis - {self.sample_size:,} reads per file for QC.'}
        </div>
"""
        
        # Overall quality assessment
        overall_pass = sum(1 for s in all_results.values() 
                          if all(a['overall'] == 'PASS' 
                          for a in s.get('assessment', {}).values()))
        overall_Careful = sum(1 for s in all_results.values() 
                          if any(a['overall'] == 'Careful' 
                          for a in s.get('assessment', {}).values()))
        overall_warn = sum(1 for s in all_results.values() 
                          if any(a['overall'] == 'warn' 
                          for a in s.get('assessment', {}).values()))
        
        # Calculate total reads analyzed
        total_reads_analyzed = sum(
            s['summary'].get('R1', {}).get('total_reads', 0) + 
            s['summary'].get('R2', {}).get('total_reads', 0)
            for s in all_results.values()
        )
        
        html_content += f"""
        <div class="summary-grid">
            <div class="summary-card">
                <h3>Total Samples</h3>
                <div class="value">{len(all_results)}</div>
            </div>
            <div class="summary-card">
                <h3>Total Reads Analyzed</h3>
                <div class="value">{total_reads_analyzed:,}</div>
            </div>
            <div class="summary-card">
                <h3>Samples Passed</h3>
                <div class="value">{overall_pass}</div>
            </div>
            <div class="summary-card">
                <h3>Samples with Warnings</h3>
                <div class="value">{overall_warn}</div>
            </div>
            <div class="summary-card">
                <h3>Samples warned</h3>
                <div class="value">{overall_warn}</div>
            </div>
        </div>
"""
        
        # Overall statistics table
        html_content += """
        <h2 style="color: #2c3e50; margin-top: 40px;">📈 Overall Statistics</h2>
        <table class="stats-table">
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Read Type</th>
                    <th>Status</th>
                    <th>Total Reads</th>
                    <th>Analyzed Reads</th>
                    <th>Mean Quality</th>
                    <th>Q30 (%)</th>
                    <th>GC (%)</th>
                    <th>Duplication (%)</th>
                </tr>
            </thead>
            <tbody>
"""
        
        for sample_name, sample_data in all_results.items():
            for read_type in ['R1', 'R2']:
                if read_type in sample_data['summary']:
                    s = sample_data['summary'][read_type]
                    assessment = sample_data['assessment'][read_type]
                    status_class = f"status-{assessment['overall'].lower()}"
                    
                    analysis_badge = "✓ Full" if s.get('is_full_analysis') else f"~ {s.get('analyzed_reads', 'N/A'):,}"
                    
                    html_content += f"""
                <tr>
                    <td><strong>{s['sample']}</strong></td>
                    <td>{s['read_type']}</td>
                    <td><span class="status-badge {status_class}">{assessment['overall']}</span></td>
                    <td>{s['total_reads']:,}</td>
                    <td>{analysis_badge}</td>
                    <td>{s['mean_quality']:.2f}</td>
                    <td>{s['q30_percentage']:.2f}</td>
                    <td>{s['mean_gc_content']:.2f}</td>
                    <td>{s['duplication_rate']*100:.2f}</td>
                </tr>
"""
        
        html_content += """
            </tbody>
        </table>
"""
        
        # Multi-sample comparison link
        if len(all_results) > 1:
            html_content += """
        <div class="alert alert-success">
            <strong>📊 Multi-Sample Comparison Available!</strong><br>
            <a href="interactive_plots/multi_sample_comparison.html" class="interactive-link" target="_blank">
                View Interactive Comparison →
            </a>
        </div>
"""
        
        # Individual sample sections
        for sample_name, sample_data in all_results.items():
            html_content += f"""
        <div class="sample-section">
            <div class="sample-header">
                <h2>📁 Sample: {sample_name}</h2>
"""
            
            # Overall status badges
            for read_type in ['R1', 'R2']:
                if read_type in sample_data['assessment']:
                    assessment = sample_data['assessment'][read_type]
                    status_class = f"status-{assessment['overall'].lower()}"
                    html_content += f"""
                <span class="status-badge {status_class}">{read_type}: {assessment['overall']}</span>
"""
            
            html_content += """
            </div>
"""
            
            # Recommendations
            recommendations = []
            for read_type in ['R1', 'R2']:
                if read_type in sample_data['assessment']:
                    assessment = sample_data['assessment'][read_type]
                    summary = sample_data['summary'][read_type]
                    
                    if assessment['metrics'].get('q30') == 'warn':
                        recommendations.append(f"{read_type}: Low Q30 percentage ({summary['q30_percentage']:.1f}%). Consider re-sequencing.")
                    if assessment['metrics'].get('duplication') == 'warn':
                        recommendations.append(f"{read_type}: High duplication rate ({summary['duplication_rate']*100:.1f}%). May indicate low library complexity.")
                    if assessment['metrics'].get('gc_content') == 'Careful':
                        recommendations.append(f"{read_type}: GC content ({summary['mean_gc_content']:.1f}%) outside expected range. Check for contamination.")
            
            if recommendations:
                html_content += """
            <div class="recommendations">
                <h3>⚠️ Recommendations</h3>
                <ul>
"""
                for rec in recommendations:
                    html_content += f"                    <li>{rec}</li>\n"
                html_content += """
                </ul>
            </div>
"""
            
            # Tabs for different views
            html_content += """
            <div class="tabs">
                <button class="tab active" onclick="openTab(event, 'metrics-""" + sample_name + """')">Metrics</button>
                <button class="tab" onclick="openTab(event, 'plots-""" + sample_name + """')">Static Plots</button>
                <button class="tab" onclick="openTab(event, 'interactive-""" + sample_name + """')">Interactive Plots</button>
            </div>
"""
            
            # Metrics tab
            html_content += f"""
            <div id="metrics-{sample_name}" class="tab-content active">
                <h3>📊 Detailed Metrics</h3>
                <table class="stats-table">
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>R1</th>
                            <th>R2</th>
                        </tr>
                    </thead>
                    <tbody>
"""
            
            # Add metrics rows
            metrics_to_show = [
                ('Total Reads', 'total_reads', ','),
                ('Analyzed Reads', 'analyzed_reads', ','),
                ('Mean Read Length', 'mean_read_length', '.1f'),
                ('Mean Quality', 'mean_quality', '.2f'),
                ('Q30 Percentage', 'q30_percentage', '.2f'),
                ('Q20 Percentage', 'q20_percentage', '.2f'),
                ('GC Content (%)', 'mean_gc_content', '.2f'),
                ('Duplication Rate (%)', 'duplication_rate', '.2%'),
            ]
            
            for metric_name, metric_key, fmt in metrics_to_show:
                r1_val = sample_data['summary'].get('R1', {}).get(metric_key, 'N/A')
                r2_val = sample_data['summary'].get('R2', {}).get(metric_key, 'N/A')
                
                if r1_val != 'N/A':
                    if fmt == ',':
                        r1_val = f"{r1_val:,}"
                    elif fmt == '.2%':
                        r1_val = f"{r1_val*100:.2f}%"
                    else:
                        r1_val = f"{r1_val:{fmt}}"
                
                if r2_val != 'N/A':
                    if fmt == ',':
                        r2_val = f"{r2_val:,}"
                    elif fmt == '.2%':
                        r2_val = f"{r2_val*100:.2f}%"
                    else:
                        r2_val = f"{r2_val:{fmt}}"
                
                html_content += f"""
                        <tr>
                            <td><strong>{metric_name}</strong></td>
                            <td>{r1_val}</td>
                            <td>{r2_val}</td>
                        </tr>
"""
            
            # Add R1-specific metrics
            if 'R1' in sample_data['summary'] and 'unique_barcodes' in sample_data['summary']['R1']:
                r1_summary = sample_data['summary']['R1']
                html_content += f"""
                        <tr>
                            <td><strong>Unique Barcodes</strong></td>
                            <td>{r1_summary['unique_barcodes']:,}</td>
                            <td>-</td>
                        </tr>
                        <tr>
                            <td><strong>Barcode Diversity</strong></td>
                            <td>{r1_summary['barcode_diversity']:.4f}</td>
                            <td>-</td>
                        </tr>
                        <tr>
                            <td><strong>Mean Barcode Quality</strong></td>
                            <td>{r1_summary['mean_barcode_quality']:.2f}</td>
                            <td>-</td>
                        </tr>
"""
                
                if 'unique_umis' in r1_summary:
                    html_content += f"""
                        <tr>
                            <td><strong>Unique UMIs</strong></td>
                            <td>{r1_summary['unique_umis']:,}</td>
                            <td>-</td>
                        </tr>
                        <tr>
                            <td><strong>UMI Diversity</strong></td>
                            <td>{r1_summary['umi_diversity']:.4f}</td>
                            <td>-</td>
                        </tr>
                        <tr>
                            <td><strong>Mean UMI Quality</strong></td>
                            <td>{r1_summary['mean_umi_quality']:.2f}</td>
                            <td>-</td>
                        </tr>
"""
            
            html_content += """
                    </tbody>
                </table>
            </div>
"""
            
            # Static plots tab
            html_content += f"""
            <div id="plots-{sample_name}" class="tab-content">
                <h3>📊 Quality Control Plots</h3>
                <div class="plot-grid">
"""
            
            for read_type in ['R1', 'R2']:
                if read_type in sample_data['plots']:
                    plots = sample_data['plots'][read_type]
                    for plot_name, plot_file in plots.items():
                        if plot_file:  # Check if plot exists
                            html_content += f"""
                    <div class="plot-item">
                        <img src="plots/{plot_file}" alt="{plot_name}">
                        <h3>{read_type}: {plot_name.replace('_', ' ').title()}</h3>
                    </div>
"""
            
            html_content += """
                </div>
            </div>
"""
            
            # Interactive plots tab
            html_content += f"""
            <div id="interactive-{sample_name}" class="tab-content">
                <h3>🎯 Interactive Visualizations</h3>
                <div style="margin: 20px 0;">
"""
            
            # Add links to interactive plots
            if 'interactive_plots' in sample_data:
                for plot_name, plot_file in sample_data['interactive_plots'].items():
                    if plot_file:
                        html_content += f"""
                    <a href="interactive_plots/{plot_file}" class="interactive-link" target="_blank">
                        {plot_name.replace('_', ' ').title()} →
                    </a>
"""
            
            html_content += """
                </div>
            </div>
        </div>
"""
        
        # Footer
        analysis_note = "Full dataset analysis" if self.use_all_reads else f"Sampled analysis ({self.sample_size:,} reads per file)"
        
        html_content += f"""
        <div class="footer">
            <p><strong>Enhanced 10x Genomics FASTQ QC Analyzer</strong></p>
            <p>Version 2.0 | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>Analysis mode: {analysis_note}</p>
            <p>For production analysis, consider using Cell Ranger for full processing.</p>
        </div>
    </div>
    
    <script>
        function openTab(evt, tabName) {{
            var i, tabcontent, tablinks;
            
            // Get all tab content and hide them
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {{
                tabcontent[i].classList.remove("active");
            }}
            
            // Get all tab buttons and remove active class
            tablinks = document.getElementsByClassName("tab");
            for (i = 0; i < tablinks.length; i++) {{
                tablinks[i].classList.remove("active");
            }}
            
            // Show current tab and mark button as active
            document.getElementById(tabName).classList.add("active");
            evt.currentTarget.classList.add("active");
        }}
    </script>
</body>
</html>
"""
        
        # Write HTML file
        html_path = self.output_dir / "qc_report.html"
        with open(html_path, 'w') as f:
            f.write(html_content)
        
        print(f"\n✅ Enhanced HTML report generated: {html_path}")
        return html_path
    
    def run_analysis(self):
        """Run complete enhanced analysis pipeline"""
        print(f"\n{'='*70}")
        print("🧬 Enhanced 10x Genomics FASTQ QC Analysis Pipeline")
        print(f"{'='*70}\n")
        
        if self.use_all_reads:
            print("📊 Analysis Mode: FULL DATASET (all reads will be analyzed)")
        else:
            print(f"📊 Analysis Mode: SAMPLED ({self.sample_size:,} reads per file)")
        
        # Find FASTQ files
        print("\n🔍 Searching for FASTQ files...")
        fastq_files = self.find_fastq_files()
        
        if not fastq_files:
            print("❌ No FASTQ files found!")
            return
        
        print(f"✅ Found {len(fastq_files)} sample(s)")
        print(f"⚙️  Using {self.threads} threads for analysis\n")
        
        all_results = {}
        
        # Analyze each sample
        for sample_name, files in fastq_files.items():
            print(f"\n{'='*70}")
            print(f"📊 Analyzing sample: {sample_name}")
            print(f"{'='*70}")
            
            sample_results = {
                'plots': {},
                'summary': {},
                'stats': {},
                'assessment': {},
                'interactive_plots': {}
            }
            
            for read_type in ['R1', 'R2']:
                if read_type not in files:
                    print(f"  ⚠️  {read_type} file not found, skipping...")
                    continue
                
                # Analyze FASTQ
                stats = self.analyze_fastq_enhanced(files[read_type], read_type)
                sample_results['stats'][read_type] = stats
                
                # Assess quality
                assessment = self.assess_quality(stats, read_type)
                sample_results['assessment'][read_type] = assessment
                print(f"  📋 Quality Assessment: {assessment['overall']}")
                
                # Generate plots
                print(f"  📈 Generating plots for {read_type}...")
                plots = {}
                
                # Static plots (keeping originals)
                plots['quality_distribution'] = self.plot_quality_distribution(
                    stats, sample_name, read_type)
                plots['quality_per_position'] = self.plot_quality_per_position(
                    stats, sample_name, read_type)
                plots['gc_content'] = self.plot_gc_content(
                    stats, sample_name, read_type)
                plots['base_composition'] = self.plot_base_composition(
                    stats, sample_name, read_type)
                plots['read_length'] = self.plot_read_length_distribution(
                    stats, sample_name, read_type)
                
                sample_results['plots'][read_type] = plots
                
                # Interactive plots
                print(f"  🎯 Creating interactive visualizations...")
                interactive_plot = self.create_interactive_quality_plot(
                    stats, sample_name, read_type)
                sample_results['interactive_plots'][f'{read_type}_interactive'] = interactive_plot
                
                # Generate summary statistics
                summary = self.generate_summary_stats(stats, sample_name, read_type)
                sample_results['summary'][read_type] = summary
                
                # R1-specific analysis
                if read_type == 'R1':
                    print(f"  🧬 Analyzing barcodes and UMIs...")
                    
                    # Knee plot
                    knee_plot = self.create_knee_plot(stats, sample_name)
                    if knee_plot:
                        sample_results['interactive_plots']['knee_plot'] = knee_plot
                    
                    # UMI analysis
                    umi_plot = self.create_umi_analysis_plot(stats, sample_name)
                    if umi_plot:
                        sample_results['interactive_plots']['umi_analysis'] = umi_plot
                    
                    # Barcode analysis (static)
                    barcode_result = self.analyze_barcode_diversity(stats, sample_name)
                    if barcode_result:
                        plot_name, unique_bc, diversity = barcode_result
                        plots['barcode_analysis'] = plot_name
            
            all_results[sample_name] = sample_results
            
            # Save individual sample statistics
            stats_file = self.stats_dir / f"{sample_name}_statistics.json"
            with open(stats_file, 'w') as f:
                json_summary = {}
                for read_type, summ in sample_results['summary'].items():
                    json_summary[read_type] = {
                        k: (int(v) if isinstance(v, (np.integer, np.int64)) else 
                            float(v) if isinstance(v, (np.floating, np.float64)) else v)
                        for k, v in summ.items() if k != 'base_composition'
                    }
                json.dump(json_summary, f, indent=2)
            
            print(f"  ✅ Statistics saved: {stats_file}")
        
        # Multi-sample comparison
        if len(all_results) > 1:
            print(f"\n{'='*70}")
            print("📊 Generating multi-sample comparison...")
            print(f"{'='*70}")
            self.create_comparison_report(all_results)
        
        # Generate enhanced HTML report
        print(f"\n{'='*70}")
        print("📝 Generating enhanced HTML report...")
        print(f"{'='*70}")
        self.create_enhanced_html_report(all_results)
        
        print(f"\n{'='*70}")
        print("✅ Analysis Complete!")
        print(f"{'='*70}")
        print(f"\n📂 Output directory: {self.output_dir}")
        print(f"   📄 HTML Report: {self.output_dir}/qc_report.html")
        print(f"   📊 Static Plots: {self.plots_dir}")
        print(f"   🎯 Interactive Plots: {self.interactive_dir}")
        print(f"   📈 Statistics: {self.stats_dir}")
        print(f"   📋 Reports: {self.reports_dir}")
        
        if self.use_all_reads:
            print(f"\n   ✅ Full dataset analysis completed")
        else:
            print(f"\n   ℹ️  Sampled {self.sample_size:,} reads per file")
        print()
    
    # Keep original plotting methods from previous version
    
    def plot_quality_distribution(self, stats, sample_name, read_type):
        """Plot quality score distribution"""
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(stats['quality_scores'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Quality Score (Phred)')
        ax.set_ylabel('Frequency')
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        ax.set_title(f'{sample_name} - {read_type}: Quality Score Distribution\n({analysis_note})')
        ax.axvline(x=30, color='r', linestyle='--', label='Q30', linewidth=2)
        ax.axvline(x=20, color='orange', linestyle='--', label='Q20', linewidth=2)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        output_path = self.plots_dir / f"{sample_name}_{read_type}_quality_dist.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name
    
    def plot_quality_per_position(self, stats, sample_name, read_type):
        """Plot quality scores across read positions"""
        positions = sorted(stats['quality_per_position'].keys())
        mean_qual = [np.mean(stats['quality_per_position'][pos]) for pos in positions]
        q25 = [np.percentile(stats['quality_per_position'][pos], 25) for pos in positions]
        q75 = [np.percentile(stats['quality_per_position'][pos], 75) for pos in positions]
        
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(positions, mean_qual, 'b-', label='Mean Quality', linewidth=2)
        ax.fill_between(positions, q25, q75, alpha=0.3, label='25-75 percentile', color='lightblue')
        ax.axhline(y=30, color='r', linestyle='--', label='Q30', linewidth=2)
        ax.axhline(y=20, color='orange', linestyle='--', label='Q20', linewidth=2)
        ax.set_xlabel('Position in Read (bp)')
        ax.set_ylabel('Quality Score')
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        ax.set_title(f'{sample_name} - {read_type}: Quality Scores Across Read Positions\n({analysis_note})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        output_path = self.plots_dir / f"{sample_name}_{read_type}_quality_per_pos.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name
    
    def plot_gc_content(self, stats, sample_name, read_type):
        """Plot GC content distribution"""
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(stats['gc_content'], bins=50, edgecolor='black', alpha=0.7, color='green')
        ax.set_xlabel('GC Content (%)')
        ax.set_ylabel('Frequency')
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        ax.set_title(f'{sample_name} - {read_type}: GC Content Distribution\n({analysis_note})')
        mean_gc = np.mean(stats['gc_content'])
        ax.axvline(x=mean_gc, color='r', linestyle='--', 
                   label=f"Mean: {mean_gc:.2f}%", linewidth=2)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        output_path = self.plots_dir / f"{sample_name}_{read_type}_gc_content.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name
    
    def plot_base_composition(self, stats, sample_name, read_type):
        """Plot base composition"""
        bases = list(stats['base_composition'].keys())
        counts = list(stats['base_composition'].values())
        total = sum(counts)
        percentages = [c/total*100 for c in counts]
        
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
        bars = ax.bar(bases, percentages, color=colors, edgecolor='black', alpha=0.8)
        ax.set_ylabel('Percentage (%)')
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled reads"
        ax.set_title(f'{sample_name} - {read_type}: Base Composition\n({analysis_note})')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add percentage labels on bars
        for bar, pct in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{pct:.2f}%', ha='center', va='bottom', fontweight='bold')
        
        output_path = self.plots_dir / f"{sample_name}_{read_type}_base_comp.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name
    
    def plot_read_length_distribution(self, stats, sample_name, read_type):
        """Plot read length distribution"""
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(stats['read_lengths'], bins=30, edgecolor='black', alpha=0.7, color='orange')
        ax.set_xlabel('Read Length (bp)')
        ax.set_ylabel('Frequency')
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled {stats.get('analyzed_reads', 'N/A'):,} reads"
        ax.set_title(f'{sample_name} - {read_type}: Read Length Distribution\n({analysis_note})')
        mean_len = np.mean(stats['read_lengths'])
        ax.axvline(x=mean_len, color='r', linestyle='--', 
                   label=f"Mean: {mean_len:.1f} bp", linewidth=2)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        output_path = self.plots_dir / f"{sample_name}_{read_type}_read_length.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name
    
    def analyze_barcode_diversity(self, stats, sample_name):
        """Analyze barcode diversity for R1"""
        if 'barcode_data' not in stats or not stats['barcode_data']:
            return None
        
        barcodes = stats['barcode_data']['barcodes']
        barcode_counts = Counter(barcodes)
        
        # Plot top barcodes
        top_n = 20
        top_barcodes = barcode_counts.most_common(top_n)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Top barcodes
        bc_names = [f"BC{i+1}" for i in range(len(top_barcodes))]
        bc_counts = [count for _, count in top_barcodes]
        
        ax1.bar(range(len(bc_counts)), bc_counts, color='steelblue', edgecolor='black', alpha=0.8)
        ax1.set_xlabel('Barcode Rank')
        ax1.set_ylabel('Count')
        ax1.set_title(f'{sample_name}: Top {top_n} Barcodes')
        ax1.set_xticks(range(len(bc_counts)))
        ax1.set_xticklabels(bc_names, rotation=45, ha='right')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Barcode diversity
        unique_barcodes = len(barcode_counts)
        total_barcodes = len(barcodes)
        diversity_ratio = unique_barcodes / total_barcodes
        
        analysis_note = "Full dataset" if stats.get('is_full_analysis') else f"Sampled reads"
        
        ax2.text(0.5, 0.75, f'Analysis: {analysis_note}', 
                ha='center', fontsize=12, transform=ax2.transAxes, style='italic')
        ax2.text(0.5, 0.6, f'Total Barcodes Analyzed: {total_barcodes:,}', 
                ha='center', fontsize=14, transform=ax2.transAxes, fontweight='bold')
        ax2.text(0.5, 0.45, f'Unique Barcodes: {unique_barcodes:,}', 
                ha='center', fontsize=14, transform=ax2.transAxes, fontweight='bold')
        ax2.text(0.5, 0.3, f'Diversity Ratio: {diversity_ratio:.4f}', 
                ha='center', fontsize=14, transform=ax2.transAxes, fontweight='bold')
        ax2.axis('off')
        ax2.set_title(f'{sample_name}: Barcode Diversity Statistics')
        
        output_path = self.plots_dir / f"{sample_name}_barcode_analysis.png"
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path.name, unique_barcodes, diversity_ratio


def main():
    parser = argparse.ArgumentParser(
        description='Enhanced 10x Genomics FASTQ QC Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with sampling (default 100k reads)
  python fastq_qc_analyzer_enhanced.py -i /path/to/fastq_files -o /path/to/output
  
  # Analyze ALL reads (full dataset)
  python fastq_qc_analyzer_enhanced.py -i /path/to/fastq_files -o /path/to/output --all
  
  # Custom sample size with multiple threads
  python fastq_qc_analyzer_enhanced.py -i ./fastq -o ./qc_results -s 200000 -t 8
  
  # Full analysis with 8 threads
  python fastq_qc_analyzer_enhanced.py -i ./fastq -o ./qc_results --all -t 8
  
  # Analyze specific directory
  python fastq_qc_analyzer_enhanced.py -i ./data/sample1 -o ./results/sample1_qc
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input directory containing FASTQ files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for results')
    parser.add_argument('-s', '--sample-size', type=int, default=100000,
                       help='Number of reads to sample per file (default: 100000). Ignored if --all is used.')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads to use (default: 4)')
    parser.add_argument('--all', '--use-all-reads', action='store_true', dest='use_all_reads',
                       help='Analyze ALL reads instead of sampling (WARNING: may take longer and use more memory)')
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.exists(args.input):
        print(f"❌ Error: Input directory '{args.input}' does not exist!")
        sys.exit(1)
    
    # Warning for full analysis
    if args.use_all_reads:
        print("\n" + "="*70)
        print("⚠️  WARNING: Full dataset analysis selected")
        print("="*70)
        print("This will analyze ALL reads in the FASTQ files.")
        print("This may take significantly longer and use more memory.")
        print("For very large files (>10M reads), consider using sampling mode.")
        print("="*70)
        
        response = input("\nDo you want to continue? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("Analysis cancelled.")
            sys.exit(0)
    
    # Run analysis
    analyzer = EnhancedFastqAnalyzer(
        args.input, 
        args.output, 
        args.sample_size,
        args.threads,
        args.use_all_reads
    )
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
