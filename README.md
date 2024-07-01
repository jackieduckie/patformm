# patformm
A little tool that takes MM tags of Biomodal bams and convert to pat format.

---

## Description

`calculate_methylated_c_positions.sh` is designed to process an input BAM file containing MM tags from Biomodal and convert them into PAT format. The script calculates the positions of methylated 'C's, converts them to CpG index values, handles invalid genomic regions, and outputs the results in a specified format.

## Usage

### Prerequisites

Ensure you have the following tools installed on your system: `samtools`, `wgbstools`, `awk`, `sort`, `uniq`

### Running the Script
1. Execute the script with the input BAM file:
   ```bash
   patformm parse_mm_tags [--threads <threads>] [-o <output_file>] <input_bam>
   ```
2. The output BED file should have the following columns:
   1. Chromosome
   2. Start position
   3. End position
   4. CpG start index
   5. CpG end index
   6. MM tag (e.g., `MM:Z:C+C.,1,23;`)
3. Execute the script with the input BED file:
   ```bash
   patformm calculate_cpos [--threads <threads>] [-o <output_file>] <input_bed>
   ```

### Output

The script produces an output file with the following format (more info [here](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md)):
- Chromosome
- CpG start index
- Methylation string
- Count of identical rows

