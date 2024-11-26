import pandas as pd
import numpy as np
import subprocess
import re
from collections import defaultdict

# Load CpG bed file into a DataFrame
def preload_cpg_bed(cpg_bed_file):
    """
    Load CpG bed file into memory with optimized data structure.
    
    Args:
        cpg_bed_file: Path to the CpG BED file (gzipped).
    
    Returns:
        cpg_map_by_chrom: Nested dictionary mapping chrom -> pos -> index
    """
    dtype = {
        'chrom': str,
        'pos': np.int32,
        'index': np.int32
    }
    
    # Read the gzipped BED file into a DataFrame
    cpg_df = pd.read_csv(
        cpg_bed_file,
        sep='\t',
        header=None,
        names=['chrom', 'pos', 'index'],
        compression='gzip',
        dtype=dtype
    )

    # Create nested dictionary for faster lookups
    cpg_map_by_chrom = defaultdict(dict)
    for _, row in cpg_df.iterrows():
        cpg_map_by_chrom[row['chrom']][row['pos']] = row['index']

    print(f"CpG map created with {len(cpg_df)} entries across {len(cpg_map_by_chrom)} chromosomes")
    return cpg_map_by_chrom

# Function to calculate the positions of methylated Cs
def calculate_positions(row, cpg_map_by_chrom):
    """
    Calculate methylation status for CpGs in a read
    
    Args:
        row: DataFrame row containing read info
        cpg_map_by_chrom: Dictionary mapping "chrom:pos" to CpG index
    """
    chrom = row['chrom']
    start = row['start']
    cigar = row['cigar']
    mm_tag = row['mm_tag']
    sequence = row['sequence']
    is_reverse = row['is_reverse']

    # Find C/G positions based on strand
    if is_reverse:
        # For reverse strand, look for Gs from right to left
        c_positions = [i for i, base in enumerate(sequence[::-1]) if base == 'G']
    else:
        # For forward strand, look for Cs from left to right
        c_positions = [i for i, base in enumerate(sequence) if base == 'C']

    # Parse MM tag and get methylated positions
    methylated_read_positions = []
    if mm_tag.startswith('MM:Z:'):
        values = mm_tag.replace('MM:Z:C+C,', '').rstrip(';').split(',')
        c_index = 0
        for v in values:
            if v.isdigit():
                c_index += int(v)
                if c_index < len(c_positions):
                    if is_reverse:
                        methylated_read_positions.append(len(sequence) - 1 - c_positions[c_index])
                    else:
                        methylated_read_positions.append(c_positions[c_index])
                c_index += 1

    # Parse CIGAR string and map read positions to reference positions
    ref_pos = start
    read_pos = 0
    read_to_ref_pos = {}
    
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    for length, op in cigar_ops:
        length = int(length)
        if op in ['M', '=', 'X']:
            for i in range(length):
                read_to_ref_pos[read_pos + i] = ref_pos + i
            read_pos += length
            ref_pos += length
        elif op in ['I', 'S']:
            read_pos += length
        elif op in ['D', 'N']:
            ref_pos += length

    # Map methylated read positions to reference positions
    methylated_ref_positions = []
    for read_pos in methylated_read_positions:
        if read_pos in read_to_ref_pos:
            methylated_ref_positions.append(read_to_ref_pos[read_pos])

    # Get methylation status for each reference position
    methylation_by_pos = {}
    min_ref_pos = min(read_to_ref_pos.values())
    max_ref_pos = max(read_to_ref_pos.values())
    
    # Find all CpGs in the read's span
    for ref_pos in range(min_ref_pos, max_ref_pos + 1):
        # For reverse reads, we need to look up the C position (one base before G)
        cpg_lookup_pos = ref_pos - 1 if is_reverse else ref_pos
        
        if chrom in cpg_map_by_chrom and cpg_lookup_pos in cpg_map_by_chrom[chrom]:
            read_pos = next((rp for rp, rp_ref in read_to_ref_pos.items() if rp_ref == ref_pos), None)
            if read_pos is None:
                methylation_by_pos[cpg_lookup_pos] = '.'  # Position in deletion or gap
            else:
                base = sequence[read_pos]
                expected_base = 'G' if is_reverse else 'C'
                
                if base == expected_base:
                    if ref_pos in methylated_ref_positions:
                        methylation_by_pos[cpg_lookup_pos] = 'C'  # Methylated
                    else:
                        methylation_by_pos[cpg_lookup_pos] = 'T'  # Unmethylated
                else:
                    methylation_by_pos[cpg_lookup_pos] = '.'  # Variant or error

    # Convert to final string (already in order by position)
    cpg_indices = []
    methylation_status = []
    
    for pos in sorted(methylation_by_pos.keys()):
        cpg_indices.append(cpg_map_by_chrom[chrom][pos])
        methylation_status.append(methylation_by_pos[pos])

    if not cpg_indices:
        return pd.Series([chrom, np.nan, '', 0])

    methylation_string = ''.join(methylation_status)
    return pd.Series([chrom, cpg_indices[0], methylation_string, 1])

# Process chunk of data
def process_chunk(file, cpg_map, chunk_size):
    """Process input file in chunks."""
    dtype = {
        'chrom': str,
        'start': np.int32,
        'cigar': str,
        'mm_tag': str,
        'sequence': str,
        'is_reverse': np.int32
    }
    
    df = pd.read_csv(file, sep='\t', header=None, 
                     names=['chrom', 'start', 'cigar', 'mm_tag', 'sequence', 'is_reverse'],
                     dtype=dtype, chunksize=chunk_size)
    
    all_results = []
    
    for i, chunk in enumerate(df):
        # Convert is_reverse to boolean (16 = reverse strand)
        chunk['is_reverse'] = chunk['is_reverse'] == 16
        
        # Filter out rows where MM tag is missing or malformed
        chunk = chunk[chunk['mm_tag'].str.startswith('MM:Z:C+C', na=False)]
        
        results = chunk.apply(calculate_positions, axis=1, cpg_map_by_chrom=cpg_map)
        results.columns = ['chrom', 'first_cpg_index', 'methylation_pattern', 'count']
        all_results.append(results)
        print(f"Processed chunk {i+1}")
    
    # Concatenate all results and group by to recalculate counts
    final_results = pd.concat(all_results)
    final_results = final_results.groupby(['chrom', 'first_cpg_index', 'methylation_pattern'])['count'].sum().reset_index()
    
    # Convert 'first_cpg_index' to int and other columns to string
    final_results['first_cpg_index'] = final_results['first_cpg_index'].astype(int)
    final_results = final_results.astype(str)

    # Write output with proper format
    output_file = f"{file}.out"
    final_results.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    chunk_size = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000
    cpg_bed_file = "/g/data/pq08/projects/biomodal/patformm/CpG.bed.hg38.gz"
    cpg_map = preload_cpg_bed(cpg_bed_file)
    print(f"CpG bed preloaded")
    print(f"Processing chunk using chunksize={chunk_size}")
    process_chunk(input_file, cpg_map, chunk_size)

