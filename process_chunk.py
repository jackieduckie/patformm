import pandas as pd
import numpy as np
import subprocess
import re

# Load CpG bed file into a DataFrame
def preload_cpg_bed(cpg_bed_file):
    dtype = {
        'chrom': str,
        'pos': np.int32,
        'index': np.int32
    }
    cpg_df = pd.read_csv(cpg_bed_file, sep='\t', header=None, names=['chrom', 'pos', 'index'], compression='gzip')
    print(f"cpg_df loaded")
    cpg_df['pos_key'] = cpg_df['chrom'] + ':' + cpg_df['pos'].astype(str)
    cpg_map = cpg_df.set_index('pos_key')['index'].to_dict()
    return cpg_map

# Function to run grep and get the CpG index
# Turns out this is much slower than reading in cpg_bed into memory
# def get_cpg_index(chrom, pos, CpG_bed):
#     try:
#         result = subprocess.check_output(['grep', f"{chrom}\t{pos}\t", CpG_bed], text=True)
        # result = subprocess.check_output(['grep', f"{chrom}\t{pos}\t", 'CpG.bed.hg38.gz'], text=True)
#         cpg_index = result.split('\t')[3]
#         return int(cpg_index)
#     except subprocess.CalledProcessError:
#         return np.nan

# Function to calculate the positions of methylated Cs
def calculate_positions(row, cpg_map):
    """
    Calculate methylation status for CpGs in a read
    
    Args:
        row: DataFrame row containing read info (chrom, start, cigar, mm_tag, sequence)
        cpg_map: Dictionary mapping "chrom:pos" to CpG index
    """
    chrom = row['chrom']
    start = row['start']
    cigar = row['cigar']
    mm_tag = row['mm_tag']
    sequence = row['sequence']

    # Parse MM tag to get methylated C positions (relative to read)
    methylated_read_positions = []
    if mm_tag.startswith('MM:Z:'):
        # Find positions of all Cs in the sequence
        c_positions = [i for i, base in enumerate(sequence) if base == 'C']
        values = mm_tag.replace('MM:Z:C+C,', '').rstrip(';').split(',')
        c_index = 0  # Index to track which C we're at
        
        for v in values:
            if v.isdigit():
                c_index += int(v)  # Skip this many Cs
                if c_index < len(c_positions):
                    methylated_read_positions.append(c_positions[c_index])  # Get the actual position of this C
                c_index += 1  # Move to next C

    # Parse CIGAR string to map read positions to reference positions
    ref_pos = start
    read_pos = 0
    read_to_ref_pos = {}  # Maps read positions to reference positions
    
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    for length, op in cigar_ops:
        length = int(length)
        if op == 'M' or op == '=' or op == 'X':  # Match or mismatch
            for i in range(length):
                read_to_ref_pos[read_pos + i] = ref_pos + i
            read_pos += length
            ref_pos += length
        elif op == 'I' or op == 'S':  # Insertion or soft clipping
            read_pos += length
        elif op == 'D' or op == 'N':  # Deletion or skip
            ref_pos += length
        # Ignore H and P operations as they don't affect coordinates

    # Map methylated read positions to reference positions
    methylated_ref_positions = []
    for read_pos in methylated_read_positions:
        if read_pos in read_to_ref_pos:
            methylated_ref_positions.append(read_to_ref_pos[read_pos])

    # Get CpG indices for methylated positions
    cpg_indices = []
    methylation_status = []
    
    # Find all possible CpGs in the read's genomic range
    min_ref_pos = min(read_to_ref_pos.values())
    max_ref_pos = max(read_to_ref_pos.values())
    
    # Check each reference position if it's in the CpG map
    for ref_pos in range(min_ref_pos, max_ref_pos + 1):
        cpg_key = f"{chrom}:{ref_pos}"
        if cpg_key in cpg_map:
            cpg_idx = cpg_map[cpg_key]
            cpg_indices.append(cpg_idx)
            # If this CpG position was methylated, mark as 'C', otherwise 'T'
            methylation_status.append('C' if ref_pos in methylated_ref_positions else 'T')

    if not cpg_indices:
        return pd.Series([chrom, np.nan, '', 0])

    methylation_string = ''.join(methylation_status)
    return pd.Series([chrom, cpg_indices[0], methylation_string, 1])

# Process chunk of data
def process_chunk(file, cpg_map, chunk_size):
    dtype = {
        'chrom': str,
        'start': np.int32,
        'cigar': str,
        'mm_tag': str,
        'sequence': str
    }
    
    df = pd.read_csv(file, sep='\t', header=None, 
                     names=['chrom', 'start', 'cigar', 'mm_tag', 'sequence'],
                     dtype=dtype, chunksize=chunk_size)
    
    all_results = []
    
    for i, chunk in enumerate(df):
        # Filter out rows where MM tag is missing or malformed
        chunk = chunk[chunk['mm_tag'].str.startswith('MM:Z:C+C', na=False)]
        
        results = chunk.apply(calculate_positions, axis=1, cpg_map=cpg_map)
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
    # process_chunk(input_file, cpg_bed_file)

