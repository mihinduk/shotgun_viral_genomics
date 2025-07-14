#!/usr/bin/env python3
"""
Viral Culture QC Visualizations
Creates publication-ready plots for viral culture quality assessment
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

def parse_diagnostic_report(report_file):
    """Parse viral diagnostic report to extract key metrics"""
    data = {
        'sample': '',
        'total_reads': 0,
        'duplicate_reads': 0,
        'unique_reads': 0,
        'raw_mapping_reads': 0,
        'raw_mapping_percent': 0,
        'dedup_mapping_reads': 0,
        'dedup_mapping_percent': 0,
        'duplication_rate': 0,
        'total_contigs': 0,
        'contigs_gt1000': 0
    }
    
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            
        # Extract sample name
        for line in content.split('\n'):
            if line.startswith('Sample:'):
                data['sample'] = line.split(':', 1)[1].strip()
                break
                
        # Extract mapping statistics
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if 'Total Reads:' in line:
                data['total_reads'] = int(line.split(':')[1].strip().replace(',', ''))
            elif 'Duplicate Reads:' in line:
                parts = line.split(':')[1].strip().split()
                data['duplicate_reads'] = int(parts[0].replace(',', ''))
                if '(' in line:
                    data['duplication_rate'] = float(line.split('(')[1].split('%')[0])
            elif 'Unique Reads:' in line:
                data['unique_reads'] = int(line.split(':')[1].strip().replace(',', ''))
            elif 'Raw Mapping:' in line:
                parts = line.split(':')[1].strip().split()
                data['raw_mapping_reads'] = int(parts[0].replace(',', ''))
                if '(' in line:
                    data['raw_mapping_percent'] = float(line.split('(')[1].split('%')[0])
            elif 'Deduplicated Mapping:' in line:
                parts = line.split(':')[1].strip().split()
                data['dedup_mapping_reads'] = int(parts[0].replace(',', ''))
                if '(' in line:
                    data['dedup_mapping_percent'] = float(line.split('(')[1].split('%')[0])
            elif 'Total Contigs:' in line:
                data['total_contigs'] = int(line.split(':')[1].strip())
            elif 'Contigs >1000bp:' in line:
                data['contigs_gt1000'] = int(line.split(':')[1].strip())
                
    except Exception as e:
        print(f"Warning: Could not parse {report_file}: {e}")
        
    return data

def parse_blast_results(blast_file):
    """Parse BLAST results to extract contamination information"""
    contamination_data = []
    
    try:
        if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
            df = pd.read_csv(blast_file, sep='\t')
            
            if len(df) > 0:
                for _, row in df.iterrows():
                    # Categorize organisms
                    title = row['Subject_Title'].lower()
                    if 'virus' in title:
                        category = 'Virus'
                    elif any(x in title for x in ['mycoplasma', 'mesomycoplasma']):
                        category = 'Mycoplasma'
                    elif any(x in title for x in ['escherichia', 'e. coli']):
                        category = 'E. coli'
                    elif any(x in title for x in ['staphylococcus', 'staph']):
                        category = 'Staphylococcus'
                    elif any(x in title for x in ['pseudomonas']):
                        category = 'Pseudomonas'
                    elif any(x in title for x in ['candida']):
                        category = 'Candida'
                    elif any(x in title for x in ['saccharomyces']):
                        category = 'Yeast'
                    else:
                        category = 'Other'
                    
                    contamination_data.append({
                        'contig_id': row['Query_ID'],
                        'category': category,
                        'organism': row['Subject_Title'][:50] + '...' if len(row['Subject_Title']) > 50 else row['Subject_Title'],
                        'identity': row['Percent_Identity'],
                        'length': row.get('Contig_Length', 0),
                        'coverage': row.get('Query_Coverage', 0)
                    })
    except Exception as e:
        print(f"Warning: Could not parse {blast_file}: {e}")
        
    return contamination_data

def create_mapping_quality_plot(samples_data, output_dir):
    """Create mapping quality assessment plot"""
    if not samples_data:
        return
        
    df = pd.DataFrame(samples_data)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Viral Culture Quality Assessment', fontsize=16, fontweight='bold')
    
    # 1. Mapping percentage (deduplicated)
    colors = ['green' if x >= 70 else 'orange' if x >= 30 else 'red' for x in df['dedup_mapping_percent']]
    bars1 = ax1.bar(df['sample'], df['dedup_mapping_percent'], color=colors, alpha=0.7)
    ax1.set_title('Mapping to Reference Genome\n(Deduplicated Reads)', fontweight='bold')
    ax1.set_ylabel('Mapping Percentage (%)')
    ax1.set_xlabel('Sample')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add threshold lines
    ax1.axhline(y=70, color='green', linestyle='--', alpha=0.7, label='Good (>70%)')
    ax1.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Moderate (30-70%)')
    ax1.legend()
    
    # Add values on bars
    for bar, val in zip(bars1, df['dedup_mapping_percent']):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # 2. Duplication rate
    colors2 = ['red' if x >= 80 else 'orange' if x >= 60 else 'green' for x in df['duplication_rate']]
    bars2 = ax2.bar(df['sample'], df['duplication_rate'], color=colors2, alpha=0.7)
    ax2.set_title('PCR Duplication Rate', fontweight='bold')
    ax2.set_ylabel('Duplication Rate (%)')
    ax2.set_xlabel('Sample')
    ax2.tick_params(axis='x', rotation=45)
    
    # Add threshold lines  
    ax2.axhline(y=80, color='red', linestyle='--', alpha=0.7, label='High (>80%)')
    ax2.axhline(y=60, color='orange', linestyle='--', alpha=0.7, label='Moderate (60-80%)')
    ax2.legend()
    
    # Add values on bars
    for bar, val in zip(bars2, df['duplication_rate']):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # 3. Assembly quality (contigs >1000bp)
    bars3 = ax3.bar(df['sample'], df['contigs_gt1000'], color='skyblue', alpha=0.7)
    ax3.set_title('Assembly Quality\n(Contigs >1000bp)', fontweight='bold')
    ax3.set_ylabel('Number of Contigs')
    ax3.set_xlabel('Sample')
    ax3.tick_params(axis='x', rotation=45)
    
    # Add values on bars
    for bar, val in zip(bars3, df['contigs_gt1000']):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{val}', ha='center', va='bottom', fontweight='bold')
    
    # 4. Read depth (millions)
    read_depth = df['total_reads'] / 1000000
    bars4 = ax4.bar(df['sample'], read_depth, color='lightcoral', alpha=0.7)
    ax4.set_title('Sequencing Depth', fontweight='bold')
    ax4.set_ylabel('Total Reads (Millions)')
    ax4.set_xlabel('Sample')
    ax4.tick_params(axis='x', rotation=45)
    
    # Add values on bars
    for bar, val in zip(bars4, read_depth):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{val:.1f}M', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'culture_quality_overview.png'), 
                bbox_inches='tight', dpi=300)
    plt.close()

def create_contamination_plot(all_contamination, output_dir):
    """Create contamination detection summary plot"""
    if not all_contamination:
        return
        
    df = pd.DataFrame(all_contamination)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Contamination Detection Summary', fontsize=16, fontweight='bold')
    
    # 1. Contamination by category
    category_counts = df['category'].value_counts()
    colors = plt.cm.Set3(np.linspace(0, 1, len(category_counts)))
    
    wedges, texts, autotexts = ax1.pie(category_counts.values, 
                                      labels=category_counts.index,
                                      autopct='%1.1f%%',
                                      colors=colors,
                                      startangle=90)
    ax1.set_title('Contamination Types Detected', fontweight='bold')
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    # 2. Sample contamination heatmap
    sample_contamination = df.groupby(['sample', 'category']).size().unstack(fill_value=0)
    
    if len(sample_contamination) > 0:
        sns.heatmap(sample_contamination, annot=True, fmt='d', cmap='Reds', 
                   ax=ax2, cbar_kws={'label': 'Number of Contigs'})
        ax2.set_title('Contamination by Sample', fontweight='bold')
        ax2.set_xlabel('Contamination Type')
        ax2.set_ylabel('Sample')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'contamination_summary.png'), 
                bbox_inches='tight', dpi=300)
    plt.close()

def create_individual_sample_plot(sample_data, contamination_data, output_dir):
    """Create detailed plot for individual sample"""
    sample_name = sample_data['sample']
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Detailed QC Report: {sample_name}', fontsize=16, fontweight='bold')
    
    # 1. Read statistics pie chart
    read_categories = ['Mapped (Unique)', 'Unmapped (Unique)', 'Duplicates']
    mapped_unique = sample_data['dedup_mapping_reads']
    unmapped_unique = sample_data['unique_reads'] - mapped_unique
    duplicates = sample_data['duplicate_reads']
    
    read_values = [mapped_unique, unmapped_unique, duplicates]
    read_colors = ['green', 'lightcoral', 'orange']
    
    wedges, texts, autotexts = ax1.pie(read_values, labels=read_categories, 
                                      autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(read_values)):,})',
                                      colors=read_colors, startangle=90)
    ax1.set_title('Read Distribution', fontweight='bold')
    
    for autotext in autotexts:
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)
    
    # 2. Mapping quality interpretation
    mapping_pct = sample_data['dedup_mapping_percent']
    
    if mapping_pct >= 70:
        quality = 'EXCELLENT'
        color = 'green'
        interpretation = 'High confidence in organism identity'
    elif mapping_pct >= 30:
        quality = 'MODERATE'
        color = 'orange'
        interpretation = 'Possible mixed infection or variant'
    else:
        quality = 'POOR'
        color = 'red'
        interpretation = 'Wrong reference or heavy contamination'
    
    ax2.text(0.5, 0.7, f'Mapping Quality', ha='center', va='center', 
             fontsize=16, fontweight='bold', transform=ax2.transAxes)
    ax2.text(0.5, 0.5, f'{quality}', ha='center', va='center',
             fontsize=24, fontweight='bold', color=color, transform=ax2.transAxes)
    ax2.text(0.5, 0.3, f'{mapping_pct:.1f}% mapping', ha='center', va='center',
             fontsize=14, transform=ax2.transAxes)
    ax2.text(0.5, 0.1, interpretation, ha='center', va='center',
             fontsize=12, style='italic', transform=ax2.transAxes, wrap=True)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis('off')
    
    # 3. Contamination breakdown
    if contamination_data:
        cont_df = pd.DataFrame([c for c in contamination_data if c.get('sample', sample_name) == sample_name])
        
        if len(cont_df) > 0:
            category_counts = cont_df['category'].value_counts()
            bars = ax3.bar(range(len(category_counts)), category_counts.values, 
                          color=plt.cm.Set3(np.linspace(0, 1, len(category_counts))))
            ax3.set_xticks(range(len(category_counts)))
            ax3.set_xticklabels(category_counts.index, rotation=45)
            ax3.set_title('Detected Contaminants', fontweight='bold')
            ax3.set_ylabel('Number of Contigs')
            
            # Add values on bars
            for bar, val in zip(bars, category_counts.values):
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                        f'{val}', ha='center', va='bottom', fontweight='bold')
        else:
            ax3.text(0.5, 0.5, 'No contamination\ndetected', ha='center', va='center',
                    fontsize=16, fontweight='bold', color='green', transform=ax3.transAxes)
            ax3.axis('off')
    else:
        ax3.text(0.5, 0.5, 'No contamination\ndata available', ha='center', va='center',
                fontsize=14, transform=ax3.transAxes)
        ax3.axis('off')
    
    # 4. Assembly statistics
    stats_text = f"""
Assembly Statistics:
• Total contigs: {sample_data['total_contigs']:,}
• Contigs >1000bp: {sample_data['contigs_gt1000']:,}
• Assembly efficiency: {(sample_data['contigs_gt1000']/max(sample_data['total_contigs'],1)*100):.1f}%

Sequencing Statistics:
• Total reads: {sample_data['total_reads']:,}
• Duplication rate: {sample_data['duplication_rate']:.1f}%
• Unique reads: {sample_data['unique_reads']:,}
    """
    
    ax4.text(0.05, 0.95, stats_text, ha='left', va='top', fontsize=11,
             transform=ax4.transAxes, fontfamily='monospace')
    ax4.set_title('Detailed Statistics', fontweight='bold')
    ax4.axis('off')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_detailed_qc.png'), 
                bbox_inches='tight', dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create QC visualizations for viral culture analysis')
    parser.add_argument('diagnostic_dirs', nargs='+', help='Diagnostic output directories')
    parser.add_argument('-o', '--output', default='qc_plots', help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Process all samples
    samples_data = []
    all_contamination = []
    
    for diagnostic_dir in args.diagnostic_dirs:
        if not os.path.exists(diagnostic_dir):
            print(f"Warning: Directory not found: {diagnostic_dir}")
            continue
            
        # Find sample name
        sample_name = os.path.basename(diagnostic_dir)
        if sample_name.startswith('diagnostic_'):
            sample_name = sample_name[11:]  # Remove 'diagnostic_' prefix
            
        # Parse diagnostic report
        report_file = os.path.join(diagnostic_dir, f"{sample_name}_diagnostic_report.txt")
        sample_data = parse_diagnostic_report(report_file)
        
        if sample_data['sample']:
            samples_data.append(sample_data)
            
            # Parse contamination data
            blast_file = os.path.join(diagnostic_dir, f"{sample_name}_top_hits.tsv")
            contamination_data = parse_blast_results(blast_file)
            
            # Add sample name to contamination data
            for item in contamination_data:
                item['sample'] = sample_data['sample']
            all_contamination.extend(contamination_data)
            
            # Create individual sample plot
            create_individual_sample_plot(sample_data, contamination_data, args.output)
            print(f"Created detailed QC plot for {sample_data['sample']}")
    
    # Create summary plots
    if samples_data:
        create_mapping_quality_plot(samples_data, args.output)
        print("Created culture quality overview plot")
        
        if all_contamination:
            create_contamination_plot(all_contamination, args.output)
            print("Created contamination summary plot")
        
        print(f"\nQC visualizations saved to: {args.output}/")
        print("Files created:")
        print("  - culture_quality_overview.png (overall quality metrics)")
        print("  - contamination_summary.png (contamination detection)")
        print("  - [sample]_detailed_qc.png (individual sample reports)")
    else:
        print("No valid diagnostic data found")

if __name__ == "__main__":
    main()