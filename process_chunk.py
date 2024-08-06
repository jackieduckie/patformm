import pandas as pd
import numpy as np
import subprocess

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
def get_cpg_index(chrom, pos, CpG_bed):
    try:
        result = subprocess.check_output(['grep', f"{chrom}\t{pos}\t", CpG_bed], text=True)
        # result = subprocess.check_output(['grep', f"{chrom}\t{pos}\t", 'CpG.bed.hg38.gz'], text=True)
        cpg_index = result.split('\t')[3]
        return int(cpg_index)
    except subprocess.CalledProcessError:
        return np.nan

# Function to calculate the positions of methylated Cs
def calculate_positions(row, cpg_map):
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    cpg_start = row['cpg_start']
    cpg_end = row['cpg_end']
    mm_tag = row['mm_tag']
    sequence = row['sequence']
    print(f"calculating position: {chrom}:{start}-{end}")
    t_values = [int(x) for x in mm_tag.rstrip(';').split(',')[1:] if x.isdigit()]

    if not t_values:
        methylation_string = 'T' * (cpg_end - cpg_start)
        return pd.Series([chrom, start, end, mm_tag, cpg_start, cpg_end, np.nan, methylation_string])

    c_positions = [start + i for i, base in enumerate(sequence) if base == 'C']
    methylated_c_positions = []
    sum_t = 0
    
    for t in t_values:
        sum_t += t + 1
        if sum_t <= len(c_positions):
            methylated_c_positions.append(c_positions[sum_t - 1])
        else:
            break

    cpg_indexes = [cpg_map.get(f"{chrom}:{pos}", np.nan) for pos in methylated_c_positions]
    # cpg_indexes = [get_cpg_index(chrom, pos, cpg_map) for pos in methylated_c_positions]
    methylation_string = ''.join(['C' if i in cpg_indexes else 'T' for i in range(cpg_start, cpg_end)])
    return pd.Series([chrom, start, end, mm_tag, cpg_start, cpg_end, cpg_indexes, methylation_string])

# Process chunk of data
def process_chunk(file, cpg_map, chunk_size):
    dtype = {
        'chrom': str,
        'start': np.int32,
        'end': np.int32,
        'cpg_start': np.int32,
        'cpg_end': np.int32,
        'mm_tag': str,
        'sequence': str
    }
    df = pd.read_csv(file, sep='\t', header=None, names=['chrom', 'start', 'end', 'cpg_start', 'cpg_end', 'mm_tag', 'sequence'], dtype=dtype, chunksize=chunk_size)
    print("bed chunk loaded")
    for i, chunk in enumerate(df):
        results = chunk.apply(calculate_positions, axis=1, cpg_map=cpg_map)
        results.to_csv(f"{file}.out", sep='\t', index=False, header=False, mode='a')
        print(f"Processed chunk {i+1}")
    # results = df.apply(calculate_positions, axis=1, cpg_map=cpg_map)
    # results.to_csv(f"{file}.out", sep='\t', index=False, header=False)
    print(f"Positions calculated to {file}.out")

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    chunk_size = int(sys.argv[2]) if len(sys.argv) > 2 else 500000
    cpg_bed_file = "/g/data/pq08/projects/biomodal/patformm/CpG.bed.hg38.gz"
    cpg_map = preload_cpg_bed(cpg_bed_file)
    print(f"CpG bed preloaded")
    print(f"Processing chunk using chunksize={chunk_size}")
    process_chunk(input_file, cpg_map, chunk_size)
    # process_chunk(input_file, cpg_bed_file)

