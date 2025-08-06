import os
import gzip
from collections import defaultdict
from .barcode_utils import extract_barcode_from_header

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
    
    return sample_patterns

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
    
    return barcode_counts, total_reads