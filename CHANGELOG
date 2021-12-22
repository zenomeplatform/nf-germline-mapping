## Changelog

PR centric changelog with description of notable changes implemented in a PR.

## v0.5 - 22/12/2021 - Adds post-alignment QC steps
  Added:
  - Process `apply_BQSR`
  - Process for post-alignment QC:
    - `qc_samtools_flagstat`
    - `qc_insert_size`
    - `qc_alignment_summary`
    - `qc_sequencing_artifact`
    - `qc_collect_hs_metrics`
    - `fastqc_mapped`
  - Post-alignment MultiQC processes:
    - `multiqc_postalignment_report_by_sample`
    - `multiqc_postalignment_report_all`
  - optional `--bait_regions` parameter

  Changes:
  - updates docker container to include R v4.1.0 (v0.7)
  - renames `--regions` to `--target_regions`
  - test data for public test
  - V7 used instead of V6 sure select target regions for full-genome tests


## v0.4 - 18/10/2021 - Adds GATK calculate BQSR step and ref data
  Added:
  - GATK calculate BQSR process step, with options to add up to three optional `--known_sites` vcf variants, and an exome regions file
  - adds new reference files: fasta.dict files and known variants vcf files for all genomes
  - SureSelect V6 and V7 exome regions files to the current repo into `assets` folder (and s3 too)
  - Adds S_cer_build3.1 genome for better CI tests
  - Adds modified TESTX csv to test merging by LB values
  - Adds test profiles to run tests with public (test) and private (test2) data not with small yeast genome but with full human genomes reference data: test_hg37, test_hg38, test2_hg37, test2_hg38.

  Changes:
  - New parameters are now required if using a custom genome: `--known_sites` and `--known_sites_index`
  - Docker container has gatk4 now installed
  - profile test uses now S_cer_build3.1 instead of sacCer3 because the latter doesn't have all refence files needed for GATK
  - Refactors test profiles to have full genome versions

## v0.3 - 11/10/2021 - Adds mapping steps
  Added:
  - Processes:
    - map_reads
    - merge_bams_by_sample
    - sort_bams_by_name
    - mark_duplicates
    - sort_bams_by_coord
  - tools bwa, samtools and picard to the main pipeline docker container (upd to v0.5)
  - igenomes config that contains paths to all public reference data
  - paths to equivalent reference data located in private yandex bucket
  - Saccharomyces cerevisiae reference data to be used in tests

  Changes:
  - the input csv file has now one more required column - LB, which should specify the read group library ID, that can be known only by researcher who has been preparing the libraries for sequencing.
  - CI tests are now run by using full set of S. cereviseae reference data (genome and bwa index) which makes them much faster and easier
  - From now only fixed container versions will be used, not `latest` to avoid future spoils of earlier pipeline versions. Current container version: v0.5.

## v0.2 - 29/02/2021 - MultiQC changes, metadata from file name
  Added:
    - Metadata retrieval from fq file names in the input channel (bio-type, seq-type, seq-machine, flowcell-ID, lane, barcode). This is optional and can be switched off by using --metadata_from_file_name false *true by default)
    - Check for correct metadata naming format for fq file (if --metadata_from_file_name is true, as it is by default)
    - Check for matching smaple_id specified in csv file and in fq file name (can be turned off by --double_check_sample_id false, true buy default)
    - Process multiqc_prealignment_report_by_sample to create smaller realignment multiqc reports for each batch of files that corresponds to same sample_id.

  Fixed:
    - Multiqc not showing fastqc report for trimmed reads (reason - "_trimmed" substrings are automatically removed from file names by multiqc before plotting, so the trimmed pair overwrote the untrimmed pair and showed up as untrimmed pair in the report)

  Changes:
    - Required input name format now needs to have also the "bio-type" field, so the example data for test2 has been changed on s3 bucket (old files retained for old tests)



## v0.1 - 29/02/2021 - Adds first pipeline template and 4 steps
  Added:
    - the basic piplene infrastrucure with multiple auxiliary elements;
    - reading input from input csv file with sample name and paths to fastq files;
    - first 4 pipeline steps:
        1. initial fastqc on raw data,
        2. read quality filtering and trimming with flexbar,
        3. secondary fastqc,
        4. collective multiqc;
    - test profiles:
        1. with public tiny s3 data,
        2. with private tiny s3 data;
    - two ci tests for the above.
