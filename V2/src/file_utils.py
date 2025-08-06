import os
import gzip
import shutil
from collections import defaultdict
from .barcode_utils import extract_barcode_from_header

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

def copy_and_append_files(original_files, recovered_reads, output_dir):
    """Copy original files to output directory and append recovered reads"""
    print(f"\nCopying and combining files to {output_dir}/...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    combined_counts = {}
    
    for sample_name, files in original_files.items():
        print(f"Processing {sample_name}...")
        
        # Process R1 and R2 files
        for read_type in ['R1', 'R2']:
            original_file = files[read_type]
            output_file = os.path.join(output_dir, os.path.basename(original_file))
            
            print(f"  Copying {os.path.basename(original_file)}...")
            shutil.copy2(original_file, output_file)
            
            if sample_name in recovered_reads:
                print(f"  Appending {len(recovered_reads[sample_name]):,} recovered reads to {read_type}...")
                with gzip.open(output_file, 'at') as out_f:
                    for read in recovered_reads[sample_name][read_type]:
                        out_f.write('\n'.join(read) + '\n')
        
        # Count total reads in final file
        total_reads = count_reads_in_file(os.path.join(output_dir, os.path.basename(files['R1'])))
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