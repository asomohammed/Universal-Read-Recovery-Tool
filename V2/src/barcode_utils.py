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