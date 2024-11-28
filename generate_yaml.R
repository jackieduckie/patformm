library(tidyverse)
library(yaml)

# Function to generate config.yml from BAM files
generate_config_yml <- function(base_path, output_file = "config.yml") {
  # Get all BAM files
  bam_files <- list.files(
    base_path, 
    pattern = "bam$", 
    recursive = TRUE, 
    full.names = FALSE
  ) %>% 
    str_subset("ctrl|gcp|synthetic|test|cellline", negate = TRUE)
  
  # Create a dataframe with parsed information
  bam_df <- tibble(full_path = bam_files) %>%
    dplyr::mutate(
      # Extract date from path
      date = str_extract(full_path, "\\d{4}-\\d{2}-\\d{2}"),
      # Extract sample name
      sample = str_extract(full_path, "[^/]+(?=\\.genome)"),
      # Get the directory path
      input_path = file.path(base_path, dirname(full_path)),
      # Get just the BAM filename
      bam = basename(full_path)
    ) %>%
    # Group by sample and take the most recent version
    dplyr::group_by(sample) %>%
    dplyr::slice_max(date, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, input_path, bam)
  
  # Convert to nested list format
  config_list <- list(
    bam_files = setNames(
      lapply(1:nrow(bam_df), function(i) {
        list(
          input_path = bam_df$input_path[i],
          bam = bam_df$bam[i]
        )
      }),
      bam_df$sample
    )
  )
  
  # Write to YAML file
  write_yaml(config_list, output_file)
  
  message(sprintf("Config file written to %s", output_file))
}

# Usage
base_path <- "/g/data/pq08/projects/biomodal/data_bucket"
output_path <- "/g/data/pq08/projects/biomodal/patformm"
output_yaml <- "config_271124.yaml"
generate_config_yml(base_path, glue::glue("{output_path}/{output_yaml}"))
