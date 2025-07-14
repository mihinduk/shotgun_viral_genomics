#!/usr/bin/env python3
"""
Simple QC Report Generator (no external dependencies)
Creates basic text-based QC report for viral culture analysis
"""

import os
import sys
import argparse
from pathlib import Path

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
            with open(blast_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) > 1:  # Skip header
                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        title = parts[8].lower() if len(parts) > 8 else ""
                        
                        # Categorize organisms
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
                            'contig_id': parts[0],
                            'category': category,
                            'organism': parts[8][:60] + '...' if len(parts) > 8 and len(parts[8]) > 60 else parts[8] if len(parts) > 8 else 'Unknown',
                            'identity': float(parts[2]) if len(parts) > 2 else 0
                        })
    except Exception as e:
        print(f"Warning: Could not parse {blast_file}: {e}")
        
    return contamination_data

def create_text_report(samples_data, all_contamination, output_dir):
    """Create a comprehensive text-based QC report"""
    
    report_file = os.path.join(output_dir, 'qc_summary_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VIRAL CULTURE QUALITY CONTROL REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Samples analyzed: {len(samples_data)}\n\n")
        
        # Overall summary
        f.write("OVERALL QUALITY SUMMARY\n")
        f.write("-" * 40 + "\n")
        
        for sample in samples_data:
            mapping_pct = sample['dedup_mapping_percent']
            dup_rate = sample['duplication_rate']
            
            # Quality assessment
            if mapping_pct >= 70:
                quality = "EXCELLENT"
                quality_symbol = "✓✓✓"
            elif mapping_pct >= 30:
                quality = "MODERATE"
                quality_symbol = "✓✓ "
            else:
                quality = "POOR"
                quality_symbol = "✗  "
                
            f.write(f"Sample: {sample['sample']:<15} | Quality: {quality:<10} {quality_symbol}\n")
            f.write(f"  Mapping: {mapping_pct:5.1f}% | Duplication: {dup_rate:5.1f}% | Contigs >1kb: {sample['contigs_gt1000']:4d}\n")
            f.write("\n")
        
        # Detailed sample analysis
        f.write("\nDETAILED SAMPLE ANALYSIS\n")
        f.write("=" * 50 + "\n")
        
        for sample in samples_data:
            f.write(f"\nSAMPLE: {sample['sample']}\n")
            f.write("-" * 30 + "\n")
            
            # Mapping analysis
            mapping_pct = sample['dedup_mapping_percent']
            f.write(f"Mapping Quality: {mapping_pct:.1f}%\n")
            
            if mapping_pct >= 70:
                f.write("  → EXCELLENT: High confidence in organism identity\n")
                f.write("  → Culture appears pure and well-identified\n")
            elif mapping_pct >= 30:
                f.write("  → MODERATE: Possible mixed infection or sequence variant\n")
                f.write("  → Consider checking for contamination or using different reference\n")
            else:
                f.write("  → POOR: Low mapping suggests wrong reference or heavy contamination\n")
                f.write("  → URGENT: Review culture purity and reference genome choice\n")
            
            # Read statistics
            f.write(f"\nSequencing Statistics:\n")
            f.write(f"  Total reads: {sample['total_reads']:,}\n")
            f.write(f"  Unique reads: {sample['unique_reads']:,}\n")
            f.write(f"  Duplication rate: {sample['duplication_rate']:.1f}%\n")
            f.write(f"  Mapped reads (dedup): {sample['dedup_mapping_reads']:,}\n")
            
            # Assembly statistics
            f.write(f"\nAssembly Statistics:\n")
            f.write(f"  Total contigs: {sample['total_contigs']:,}\n")
            f.write(f"  Contigs >1000bp: {sample['contigs_gt1000']:,}\n")
            
            # Assembly efficiency
            if sample['total_contigs'] > 0:
                efficiency = (sample['contigs_gt1000'] / sample['total_contigs']) * 100
                f.write(f"  Assembly efficiency: {efficiency:.1f}%\n")
            
            f.write("\n" + "=" * 50)
        
        # Contamination analysis
        if all_contamination:
            f.write("\n\nCONTAMINATION ANALYSIS\n")
            f.write("=" * 40 + "\n")
            
            # Count contamination types
            contamination_counts = {}
            sample_contamination = {}
            
            for item in all_contamination:
                category = item['category']
                sample = item.get('sample', 'Unknown')
                
                if category not in contamination_counts:
                    contamination_counts[category] = 0
                contamination_counts[category] += 1
                
                if sample not in sample_contamination:
                    sample_contamination[sample] = {}
                if category not in sample_contamination[sample]:
                    sample_contamination[sample][category] = 0
                sample_contamination[sample][category] += 1
            
            # Overall contamination summary
            f.write("Contamination Types Detected:\n")
            for category, count in sorted(contamination_counts.items()):
                f.write(f"  {category:<15}: {count:3d} contigs\n")
            
            f.write(f"\nDetailed Contamination by Sample:\n")
            f.write("-" * 40 + "\n")
            
            for sample in samples_data:
                sample_name = sample['sample']
                f.write(f"\n{sample_name}:\n")
                
                if sample_name in sample_contamination:
                    for category, count in sorted(sample_contamination[sample_name].items()):
                        f.write(f"  {category:<15}: {count:3d} contigs\n")
                        
                        # Special warnings
                        if category == 'Mycoplasma' and count > 0:
                            f.write("    ⚠️  MYCOPLASMA DETECTED - Critical cell culture contaminant!\n")
                        elif category in ['E. coli', 'Staphylococcus', 'Pseudomonas'] and count > 0:
                            f.write("    ⚠️  Bacterial contamination detected\n")
                        elif category in ['Candida', 'Yeast'] and count > 0:
                            f.write("    ⚠️  Fungal contamination detected\n")
                else:
                    f.write("  No contamination detected ✓\n")
        
        # Recommendations
        f.write("\n\nRECOMMENDATIONS\n")
        f.write("=" * 30 + "\n")
        
        for sample in samples_data:
            sample_name = sample['sample']
            mapping_pct = sample['dedup_mapping_percent']
            dup_rate = sample['duplication_rate']
            
            f.write(f"\n{sample_name}:\n")
            
            # Mapping recommendations
            if mapping_pct < 30:
                f.write("  🔴 URGENT: Very low mapping - verify reference genome\n")
                f.write("     Consider BLAST analysis to identify correct organism\n")
            elif mapping_pct < 70:
                f.write("  🟡 CAUTION: Moderate mapping - check for contamination\n")
                f.write("     Consider culture purification if mapping is critical\n")
            else:
                f.write("  🟢 GOOD: High mapping confidence\n")
            
            # Duplication recommendations
            if dup_rate > 80:
                f.write("  🟡 High duplication rate - consider PCR optimization\n")
            
            # Contamination recommendations
            sample_contaminants = [item for item in all_contamination if item.get('sample', sample_name) == sample_name]
            
            mycoplasma_found = any(item['category'] == 'Mycoplasma' for item in sample_contaminants)
            bacterial_found = any(item['category'] in ['E. coli', 'Staphylococcus', 'Pseudomonas'] for item in sample_contaminants)
            
            if mycoplasma_found:
                f.write("  🔴 CRITICAL: Mycoplasma contamination detected\n")
                f.write("     Immediate culture treatment or disposal recommended\n")
            elif bacterial_found:
                f.write("  🟡 Bacterial contamination present - consider culture cleaning\n")
            elif not sample_contaminants:
                f.write("  🟢 No obvious contamination detected\n")
        
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Create simple QC report for viral culture analysis')
    parser.add_argument('diagnostic_dirs', nargs='+', help='Diagnostic output directories')
    parser.add_argument('-o', '--output', default='qc_report', help='Output directory for report')
    
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
            
            print(f"Processed sample: {sample_data['sample']}")
    
    # Create report
    if samples_data:
        create_text_report(samples_data, all_contamination, args.output)
        print(f"\nQC report saved to: {args.output}/qc_summary_report.txt")
        print("This text-based report is ready for your presentation!")
    else:
        print("No valid diagnostic data found")

if __name__ == "__main__":
    main()