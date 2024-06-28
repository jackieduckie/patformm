# patformm
A little tool that takes MM tags of Biomodal bams and convert to pat format.

---

## Description

`calculate_methylated_c_positions.sh` is a Bash script designed to process an input BED file containing MM tags from Biomodal BAM files and convert them into PAT format. The script calculates the positions of methylated 'C's, converts them to CpG index values, handles invalid genomic regions, and outputs the results in a specified format.

## The scripts do the following:

1. **Filter Rows**: Removes rows where the methylation string column contains the value `NA`.
2. **Extract Columns**: Extracts specific columns (chromosome, CpG start index, and methylation string) from the input file.
3. **Sort Data**: Sorts the extracted data by the CpG start index.
4. **Count Unique Rows**: Counts unique rows and appends the count as the last column.
5. **Metadata Printing**: Optionally prints metadata about the script.

## Usage

### Prerequisites

Ensure you have the following tools installed on your system: `samtools`, `wgbstools`, `awk`, `sort`, `uniq`

### Input File Format

The input BED file should have the following columns:
1. Chromosome
2. Start position
3. End position
4. CpG start index
5. CpG end index
6. MM tag (e.g., `MM:Z:C+C.,1,23;`)

### Running the Script

1. Save the script `calculate_methylated_c_positions.sh` to your desired directory.
2. Make the script executable:
   ```bash
   chmod +x calculate_methylated_c_positions.sh
   ```
3. Execute the script with the input BED file:
   ```bash
   ./calculate_methylated_c_positions.sh input.bed
   ```
   Optionally, you can specify the number of threads and output file:
   ```bash
   ./calculate_methylated_c_positions.sh --threads 4 -o output.txt input.bed
   ```

### Output

The script produces an output file with the following format (more info [here](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md)):
- Chromosome
- CpG start index
- Methylation string
- Count of identical rows

