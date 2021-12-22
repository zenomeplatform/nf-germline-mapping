# nf-germline-snv

Nextflow pipeline for germline variant calling.

## Introduction

Nf-germline-snv is a genomics pipeline designed with Nextflow and Docker to use common bioinformatic tools to perform germline variant calling.

Pipeline uses the following tools:
 - fastqc
 - flexbar
 - multiqc
 - bwa
 - samtools
 - bcftools
 - picard
 - gatk
 - manta
 - strelka

## Quick Start
1. Install [`Nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) (Needs root permissions. Don't forget to activate docker service after installation)

3. Install [`Graphviz`](https://graphviz.org/download/) (optional)

4. Download the pipeline and run it with test data with a single command:
```
nextflow run . -profile test
```

Test data but with full reference genome:
```
nextflow run . -profile test_hg37
nextflow run . -profile test_hg38
```

With external yandex s3 bucket:
```
nextflow run . -profile test2,yandex --accessKey <accessKey> --secretKey <secretKey>
nextflow run . -profile test2_hg37,yandex --accessKey <accessKey> --secretKey <secretKey>
nextflow run . -profile test2_hg38,yandex --accessKey <accessKey> --secretKey <secretKey>
```


5. Using your own data:
```
nextflow run . --input test_data_paths.csv --metadata_from_file_name false
```

## Pipeline Inputs

Pipeline takes as an input raw paired-end fastq[.gz] files, provided as paths via an input csv file:

`input.csv`:
```
sample,fastq_1,fastq_2
sample_S1_L001,path/to/data/sample_S1_L001_R1_001.fastq.gz,path/to/data/sample_S1_L001_R2_001.fastq.gz
sample_S1_L002,path/to/data/sample_S1_L001_R1_002.fastq.gz,path/to/data/sample_S1_L001_R2_002.fastq.gz
```
First line of file must be a header. First CSV column must contain sample names, other two columns must specify paths to forward and reverse fastq files. It is possible to specify s3 paths to AWS s3 buckets. If you wish to specify s3 paths to non-amazon s3 buckets, create and use a config similar to [conf/yandex.config](https://github.com/zenomeplatform/nf-germline-snv/blob/main/conf/yandex.config).

### Naming convention
By default the pipeline relies on obtaining certain metadata information directly from the fastq file names:
- sample_ID
- bio-type
- seq-type
- seq-machine
- flowcell-ID
- lane
- barcode

To provide that kind of information fastq file names must follow the naming convention:
```
{sampleID}-{bio-type}-{seq-type}-{seq-machine}-{flowcell-ID}-{lane}-{barcode}-{read-direction}.fq[.gz]
```
Where individual items are separated by single dashes `-` (can be changed with `--sample_name_format_delimeter` parameter)

Example:
```
ZD210122-GED-E080AS6-MG-300056277-3-66-F.fq.gz
```

If you do not want to provide metadata in such way and want to treat all your samples independently, you can use pipeline parameter `--metadata_from_file_name false`

### Bio-type
Possible values:
`GED` - germline DNA
`SOD` - somatic DNA
`GER` - germline RNA
`SOR` - somatic RNA
`CFD` - cell free DNA
`SCD` - single cell DNA

## Optional inputs

### Genome reference data

By default pipeline uses reference data from igenomes s3 bucket for pre-defined genomes hg37 and hg38 (GRCh38). Genome hg38 is used by default, but that can be changed with `--genome` parameter.

Alternatively, user may provide his own reference data bundle with corresponding arguments. Note that full list of reference data will have to be provided in this case (http and s3 links are also accepted). Check the `conf/igenomes.config` file to see which files are already in use and which reference data files are needed for the pipeline in general.

### Targeted sequencing regions files (exome sequening)
For targeted sequencing projects such as exome sequencing it is required to provide target and bait regions files, which should be in gatk interval_list format. The pipeline already comes with prepared SureSelect V6 S07604514 and V7 S31285117 files converted to iterval_list format for both hg37 and hg38. You can find these files and instructions how they were obtained in `assets/regions_files` folder.

Pipeline requires at least one targeted sequencing regions file:
```
--target_regions SS_V7_hg19_regions.bed.interval_list
```

If available, also the "probes" bait regions file should be provided:
```
--bait_regions SS_V7_hg19_probes.bed.interval_list
```
If only target regions file is provided, it will be used both as target and bait file in target seqencing specificity QC step. However this is less accurate then if providing both files.


## Test profiles

Using public test data:
```
nextflow run . -profile test
```

Using private test data:
```
nextflow run . -profile test2,yandex --accessKey '<accessKey>' --secretKey '<secretKey>'
```
where `<accessKey>` and `<secretKey>` are credentials to private s3 bucket, for example s3://zenome-ngs-data .

