# patformm

A tool for converting methylation information from Biomodal BAM files (MM tags) to PAT format for methylation analysis.

## Description

`patformm` processes BAM files containing methylation information in MM tags (e.g., `MM:Z:C+C.,1,23;`) and converts them to PAT format through a two-step process:

1. **parse_mm_tags**: Extracts methylation information from BAM files
2. **calculate_cpos**: Processes methylation patterns and converts to CpG indices

### Core Components

#### parse_mm_tags
Extracts essential information from BAM files:
- Chromosome
- Start position
- CIGAR string
- MM tag
- Sequence
- Strand information (reverse/forward)

```bash
patformm parse_mm_tags --threads 8 -o output.bed input.bam
```

#### calculate_cpos
Processes the extracted information to generate PAT format:
- Maps read positions to reference positions using CIGAR strings
- Identifies methylated C positions from MM tags
- Converts to CpG indices using a reference CpG map
- Handles both forward and reverse strands
- Supports chunked processing for memory efficiency

```bash
patformm calculate_cpos --threads 8 --chunk-size 1000000 -o output.pat input.bed
```

## Implementation Details

### MM Tag Processing
- Parses MM tags in format `MM:Z:C+C,<positions>`
- Tracks methylated C positions in reads
- Handles both forward (C) and reverse (G) strand methylation

### CpG Index Mapping
- Uses a preloaded CpG reference map (`CpG.bed.hg38.gz`)
- Maps genomic positions to CpG indices
- Handles deletions and gaps with '.' notation

### Output Format
PAT format with:
```
chromosome  first_cpg_index  methylation_pattern  count
chr1        1234            CT.C                 5
```
Where:
- `methylation_pattern`: C=methylated, T=unmethylated, .=missing/invalid

## Performance Features

- Parallel processing support
- Chunked file processing
- Memory-efficient design
- Temporary file handling in scratch space

## Requirements

- Python 3.8+
- samtools
- wgbstools
- Reference files:
  - CpG.bed.hg38.gz (CpG position index)

## Usage Example

```bash
# Step 1: Parse MM tags
patformm parse_mm_tags --threads 8 -o sample.bed sample.bam

# Step 2: Calculate CpG positions
patformm calculate_cpos --threads 8 --chunk-size 1000000 -o sample.pat sample.bed
```

## Notes

- Large BAM files are processed in chunks to manage memory usage
- Supports multi-threaded processing for improved performance
