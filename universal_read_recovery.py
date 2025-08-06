#!/usr/bin/env python3
import gzip
from collections import defaultdict
import os
import re
import csv
import argparse
from pathlib import Path
import shutil
from datetime import datetime


def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings of equal length"""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def fuzzy_match_barcode(query_barcode, reference_patterns, max_mismatches=2):
    """
    Match barcode against reference patterns with fuzzy matching
    Handles 'N' bases and allows for sequencing errors
    """
    best_match = None
    best_score = float('inf')
    
    for pattern, sample_name in reference_patterns.items():
        # Handle the '+' separator in patterns
        if '+' in pattern and '+' in query_barcode:
            i1_pattern, i2_pattern = pattern.split('+')
            i1_query, i2_query = query_barcode.split('+')
        else:
            # Single barcode comparison
            i1_pattern, i2_pattern = pattern, ""
            i1_query, i2_query = query_barcode, ""
        
        # Calculate mismatches for each index
        i1_mismatches = calculate_mismatches_with_n(i1_query, i1_pattern)
        i2_mismatches = calculate_mismatches_with_n(i2_query, i2_pattern) if i2_pattern else 0
        
        total_mismatches = i1_mismatches + i2_mismatches
        
        if total_mismatches <= max_mismatches and total_mismatches < best_score:
            best_score = total_mismatches
            best_match = sample_name
    
    return best_match, best_score

def calculate_mismatches_with_n(query, pattern):
    """Calculate mismatches, treating 'N' as wildcards"""
    if len(query) != len(pattern):
        return float('inf')
    
    mismatches = 0
    for q, p in zip(query, pattern):
        if q != 'N' and p != 'N' and q != p:
            mismatches += 1
    return mismatches

def extract_barcode_from_header(header):
    """Extract barcode from FASTQ header"""
    parts = header.split()
    if len(parts) >= 2:
        info_parts = parts[1].split(':')
        if len(info_parts) >= 4:
            return info_parts[3]
    return None

def discover_assigned_samples(directory=".", sample_reads=10000):
    """
    Automatically discover all assigned sample files in the directory
    Returns dictionary of sample_name -> barcode_patterns
    """
    print(f"Discovering assigned samples in: {os.path.abspath(directory)}")
    print(f"Sampling {sample_reads:,} reads per file for barcode discovery")
    
    # Find all FASTQ files that are not Unassigned/Undetermined
    fastq_files = []
    for file in os.listdir(directory):
        if (file.endswith('.fastq.gz') and 
            'R1' in file and 
            not file.startswith('Unassigned') and 
            not file.startswith('Undetermined')):
            fastq_files.append(file)
    
    print(f"Found {len(fastq_files)} assigned sample files:")
    for f in sorted(fastq_files):
        print(f"  {f}")
    
    # Extract barcodes from each sample file
    sample_patterns = {}
    
    for fastq_file in fastq_files:
        sample_name = fastq_file.replace('_R1_001.fastq.gz', '')
        print(f"\nAnalyzing {sample_name}...")
        
        barcode_counts = defaultdict(int)
        read_count = 0
        
        try:
            with gzip.open(os.path.join(directory, fastq_file), 'rt') as f:
                while read_count < sample_reads:
                    header = f.readline().strip()
                    if not header:
                        break
                    
                    # Skip the other 3 lines
                    for _ in range(3):
                        f.readline()
                    
                    barcode = extract_barcode_from_header(header)
                    if barcode:
                        barcode_counts[barcode] += 1
                    read_count += 1
            
            # Get the most common barcode as the representative pattern
            if barcode_counts:
                most_common_barcode = max(barcode_counts.items(), key=lambda x: x[1])
                sample_patterns[most_common_barcode[0]] = sample_name
                print(f"  Representative barcode: {most_common_barcode[0]} ({most_common_barcode[1]} reads)")
            
        except Exception as e:
            print(f"  Error reading {fastq_file}: {e}")
    
    print(f"\nDiscovered {len(sample_patterns)} sample patterns:")
    for pattern, sample in sample_patterns.items():
        print(f"  {pattern} -> {sample}")
    
    return sample_patterns

def find_unassigned_files(directory="."):
    """Find unassigned/undetermined FASTQ files"""
    unassigned_files = {}
    
    # Look for both "Unassigned" and "Undetermined" prefixes
    prefixes = ['Unassigned', 'Undetermined']
    
    for file in os.listdir(directory):
        if file.endswith('.fastq.gz'):
            for prefix in prefixes:
                if file.startswith(prefix):
                    if 'R1' in file:
                        unassigned_files['R1'] = os.path.join(directory, file)
                    elif 'R2' in file:
                        unassigned_files['R2'] = os.path.join(directory, file)
    
    return unassigned_files

def analyze_unassigned_barcodes(unassigned_r1_file, sample_reads=1000000):
    """Analyze barcodes in unassigned file"""
    print(f"Analyzing unassigned reads (sampling up to {sample_reads:,} reads)...")
    
    barcode_counts = defaultdict(int)
    total_reads = 0
    
    if not os.path.exists(unassigned_r1_file):
        print(f"Error: {unassigned_r1_file} not found!")
        return {}, 0
    
    with gzip.open(unassigned_r1_file, 'rt') as f:
        while total_reads < sample_reads:
            header = f.readline().strip()
            if not header:
                break
            
            # Skip the other 3 lines
            for _ in range(3):
                f.readline()
            
            barcode = extract_barcode_from_header(header)
            if barcode:
                barcode_counts[barcode] += 1
            total_reads += 1
            
            if total_reads % 1000000 == 0:
                print(f"  Analyzed {total_reads:,} unassigned reads...")
    
    print(f"Sampled {total_reads:,} unassigned reads with {len(barcode_counts)} unique barcodes")
    return barcode_counts, total_reads

def copy_and_append_files(original_files, recovered_reads, output_dir):
    """Copy original files to output directory and append recovered reads"""
    print(f"\nCopying and combining files to {output_dir}/...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    combined_counts = {}
    
    for sample_name, files in original_files.items():
        print(f"Processing {sample_name}...")
        
        # Copy and append R1
        original_r1 = files['R1']
        output_r1 = os.path.join(output_dir, os.path.basename(original_r1))
        
        # Copy original file
        print(f"  Copying {os.path.basename(original_r1)}...")
        shutil.copy2(original_r1, output_r1)
        
        # Append recovered reads
        if sample_name in recovered_reads:
            print(f"  Appending {len(recovered_reads[sample_name]):,} recovered reads to R1...")
            with gzip.open(output_r1, 'at') as out_f:
                for r1_read in recovered_reads[sample_name]['R1']:
                    out_f.write('\n'.join(r1_read) + '\n')
        
        # Copy and append R2
        original_r2 = files['R2']
        output_r2 = os.path.join(output_dir, os.path.basename(original_r2))
        
        print(f"  Copying {os.path.basename(original_r2)}...")
        shutil.copy2(original_r2, output_r2)
        
        if sample_name in recovered_reads:
            print(f"  Appending {len(recovered_reads[sample_name]):,} recovered reads to R2...")
            with gzip.open(output_r2, 'at') as out_f:
                for r2_read in recovered_reads[sample_name]['R2']:
                    out_f.write('\n'.join(r2_read) + '\n')
        
        # Count total reads in final file
        total_reads = count_reads_in_file(output_r1)
        recovered_count = len(recovered_reads.get(sample_name, {}).get('R1', []))
        original_count = total_reads - recovered_count
        
        combined_counts[sample_name] = {
            'original': original_count,
            'recovered': recovered_count,
            'total': total_reads
        }
        
        print(f"  Final counts - Original: {original_count:,}, Recovered: {recovered_count:,}, Total: {total_reads:,}")
    
    return combined_counts

def count_reads_in_file(filepath):
    """Count total reads in a FASTQ file"""
    count = 0
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count

def recover_reads(max_mismatches=2, directory=".", output_csv="recovery_report.csv", 
                 output_html="recovery_report.html", sample_reads=1000000, discovery_reads=10000):
    """Main recovery function"""
    
    print("Universal Read Recovery Tool")
    print("===========================")
    print(f"Working directory: {os.path.abspath(directory)}")
    print(f"Max mismatches allowed: {max_mismatches}")
    print(f"Sample reads for analysis: {sample_reads:,}")
    print(f"Discovery reads per sample: {discovery_reads:,}")
    
    # First, let's see what files we have
    print(f"\nAll FASTQ files in directory:")
    fastq_files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
    for f in sorted(fastq_files):
        size_mb = os.path.getsize(os.path.join(directory, f)) / (1024**2)
        print(f"  {f} ({size_mb:.1f} MB)")
    
    # Discover assigned samples automatically
    reference_patterns = discover_assigned_samples(directory, discovery_reads)
    
    if not reference_patterns:
        print("No assigned samples found! Make sure you have demultiplexed FASTQ files in the directory.")
        return
    
    # Find original files for each sample
    original_files = {}
    for pattern, sample_name in reference_patterns.items():
        r1_file = None
        r2_file = None
        
        for file in os.listdir(directory):
            if sample_name in file and file.endswith('.fastq.gz'):
                if 'R1' in file:
                    r1_file = os.path.join(directory, file)
                elif 'R2' in file:
                    r2_file = os.path.join(directory, file)
        
        if r1_file and r2_file:
            original_files[sample_name] = {'R1': r1_file, 'R2': r2_file}
    
    # Find unassigned files
    unassigned_files = find_unassigned_files(directory)
    if not unassigned_files.get('R1') or not unassigned_files.get('R2'):
        print("Unassigned/Undetermined files not found!")
        print("Looking for files starting with 'Unassigned' or 'Undetermined'")
        return
    
    print(f"\nUnassigned files found:")
    for read, path in unassigned_files.items():
        size_gb = os.path.getsize(path) / (1024**3)
        print(f"  {read}: {os.path.basename(path)} ({size_gb:.2f} GB)")
    
    # Analyze unassigned barcodes
    unassigned_barcodes, total_sampled = analyze_unassigned_barcodes(unassigned_files['R1'], sample_reads)
    
    if total_sampled == 0:
        print("No reads found in unassigned file!")
        return
    
    # Show top unassigned barcodes
    print(f"\nTop 20 unassigned barcodes:")
    for barcode, count in sorted(unassigned_barcodes.items(), key=lambda x: x[1], reverse=True)[:20]:
        pct = (count/total_sampled)*100
        print(f"  {barcode}: {count:,} reads ({pct:.1f}%)")
    
    # Test fuzzy matching
    print(f"\nTesting fuzzy matching (max {max_mismatches} mismatches):")
    recoverable_barcodes = {}
    
    for barcode, count in unassigned_barcodes.items():
        match, score = fuzzy_match_barcode(barcode, reference_patterns, max_mismatches)
        if match:
            recoverable_barcodes[barcode] = (match, score, count)
    
    # Sort by read count and show top matches
    sorted_recoverable = sorted(recoverable_barcodes.items(), key=lambda x: x[1][2], reverse=True)
    
    if sorted_recoverable:
        print(f"Top 10 recoverable patterns:")
        for barcode, (sample, score, count) in sorted_recoverable[:10]:
            pct = (count/total_sampled)*100
            print(f"  {barcode} -> {sample} (mismatches: {score}, reads: {count:,}, {pct:.1f}%)")
    else:
        print("No recoverable patterns found!")
    
    total_recoverable = sum(count for _, _, count in recoverable_barcodes.values())
    recovery_rate = (total_recoverable/total_sampled)*100 if total_sampled > 0 else 0
    
    print(f"\nRecovery Summary:")
    print(f"  Total sampled reads: {total_sampled:,}")
    print(f"  Recoverable reads: {total_recoverable:,} ({recovery_rate:.1f}%)")
    print(f"  Unique recoverable barcodes: {len(recoverable_barcodes)}")
    
    if total_recoverable == 0:
        print("No reads can be recovered with current parameters.")
        print(f"Try increasing max_mismatches (currently {max_mismatches})")
        return
    
    # Ask user if they want to proceed
    response = input(f"\nProceed with recovery? (y/n): ")
    if response.lower() != 'y':
        print("Recovery cancelled.")
        return
    
    # Perform recovery
    print("\nPerforming read recovery...")
    recovered_reads = defaultdict(lambda: {'R1': [], 'R2': []})
    recovery_counts = defaultdict(int)
    processed_reads = 0
    
    with gzip.open(unassigned_files['R1'], 'rt') as r1, \
         gzip.open(unassigned_files['R2'], 'rt') as r2:
        
        while True:
            try:
                # Read one record from each file
                r1_rec = [r1.readline().strip() for _ in range(4)]
                r2_rec = [r2.readline().strip() for _ in range(4)]
                
                if not r1_rec[0]:
                    break
                
                barcode = extract_barcode_from_header(r1_rec[0])
                
                if barcode in recoverable_barcodes:
                    sample, score, _ = recoverable_barcodes[barcode]
                    
                    # Store recovered reads
                    recovered_reads[sample]['R1'].append(r1_rec)
                    recovered_reads[sample]['R2'].append(r2_rec)
                    
                    recovery_counts[sample] += 1
                
                processed_reads += 1
                if processed_reads % 1000000 == 0:
                    print(f"  Processed {processed_reads:,} reads...")
                
            except Exception as e:
                print(f"Error at read {processed_reads}: {e}")
                break
    
    print(f"Processed {processed_reads:,} unassigned reads")
    print("Recovered reads by sample:")
    total_recovered = sum(recovery_counts.values())
    for sample, count in sorted(recovery_counts.items()):
        pct = (count/total_recovered)*100 if total_recovered > 0 else 0
        print(f"  {sample}: {count:,} reads ({pct:.1f}%)")
    
    # Create recovered directory and combine files
    output_dir = os.path.join(directory, "recovered")
    combined_counts = copy_and_append_files(original_files, recovered_reads, output_dir)
    
    # Create simple HTML report
    html_path = os.path.join(directory, output_html)
    with open(html_path, 'w') as f:
        f.write(f"""
<!DOCTYPE html>
<html>
<head>
    <title>Read Recovery Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .number {{ font-family: monospace; }}
    </style>
</head>
<body>
    <h1>Read Recovery Report</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    
    <h2>Summary</h2>
    <table>
        <tr><th>Parameter</th><th>Value</th></tr>
        <tr><td>Max Mismatches</td><td class="number">{max_mismatches}</td></tr>
        <tr><td>Reads Sampled</td><td class="number">{total_sampled:,}</td></tr>
        <tr><td>Recoverable Reads</td><td class="number">{total_recoverable:,}</td></tr>
        <tr><td>Recovery Rate</td><td class="number">{recovery_rate:.2f}%</td></tr>
    </table>
    
    <h2>Sample Results</h2>
    <table>
        <tr><th>Sample</th><th>Original</th><th>Recovered</th><th>Total</th><th>Recovery %</th></tr>
        """)
        
        for sample, counts in sorted(combined_counts.items()):
            recovery_pct = (counts['recovered']/counts['total'])*100 if counts['total'] > 0 else 0
            f.write(f"""
        <tr>
            <td>{sample}</td>
            <td class="number">{counts['original']:,}</td>
            <td class="number">{counts['recovered']:,}</td>
            <td class="number">{counts['total']:,}</td>
            <td class="number">{recovery_pct:.2f}%</td>
        </tr>""")
        
        f.write("""
    </table>
</body>
</html>""")
    
    # Create CSV report
    csv_path = os.path.join(directory, output_csv)
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sample', 'Original_Reads', 'Recovered_Reads', 'Total_Reads', 'Recovery_Percent'])
        for sample, counts in sorted(combined_counts.items()):
            recovery_pct = (counts['recovered']/counts['total'])*100 if counts['total'] > 0 else 0
            writer.writerow([sample, counts['original'], counts['recovered'], counts['total'], f"{recovery_pct:.2f}"])
    
    print(f"\nRecovery Complete!")
    print(f"Total recovered: {total_recovered:,} reads")
    print(f"Combined files saved to: {output_dir}/")
    print(f"HTML report: {html_path}")
    print(f"CSV report: {csv_path}")

def main():
    parser = argparse.ArgumentParser(description='Universal Read Recovery Tool')
    parser.add_argument('-m', '--max-mismatches', type=int, default=2,
                       help='Maximum number of mismatches allowed (default: 2)')
    parser.add_argument('-d', '--directory', type=str, default='.',
                       help='Directory containing FASTQ files (default: current directory)')
    parser.add_argument('-o', '--output', type=str, default='recovery_report.csv',
                       help='Output CSV report filename (default: recovery_report.csv)')
    parser.add_argument('-s', '--sample-reads', type=int, default=1000000,
                       help='Number of unassigned reads to sample for analysis (default: 1,000,000)')
    parser.add_argument('-r', '--discovery-reads', type=int, default=10000,
                       help='Number of reads per sample for barcode discovery (default: 10,000)')
    
    args = parser.parse_args()
    
    recover_reads(max_mismatches=args.max_mismatches, 
                 directory=args.directory, 
                 output_csv=args.output,
                 sample_reads=args.sample_reads,
                 discovery_reads=args.discovery_reads)

if __name__ == "__main__":
    main()
