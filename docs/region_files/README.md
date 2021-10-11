# Targeted sequencing region files
Targeted sequencing uses special kits that are designed to target only particluar parts of genome, most often - exome regions.
Current pipeline is designed to accomodate for such targeted sequencing by taking as an input (optional) regions files, that specifies the target regions. This file is used in calculate_BQSR step of the pipeline. Detailed docs on how the files below were obtained can be found in corresponding subfolders at assets/regions_files folder.

## Available Files
Here you can find a list of region files that can be used with the pipeline.

### Available Files in the repository
#### SureSelect V6
```
# hg19/hg37 (Ensemble chromosome nomenclature: "1,2,3,X,Y")
assets/regions_files/V6/hg37/SS_V6_hg19_regions.sorted.merged.no_chr.bed
assets/regions_files/V6/hg37/SS_V6_hg19_regions.bed.interval_list

# hg38 (UCSC chromosome nomenclature: "chr1,chr2,chr3,chrX,chrY")
assets/regions_files/V6/hg38/SS_V6_hg38_regions.sorted.merged.bed
assets/regions_files/V6/hg38/SS_V6_hg38_regions.bed.interval_list
```

#### SureSelect V7
```
# hg19/hg37 (Ensemble chromosome nomenclature: "1,2,3,X,Y")
assets/regions_files/V7/hg37/SS_V7_hg19_regions.sorted.merged.no_chr.bed
assets/regions_files/V7/hg37/SS_V7_hg19_regions.bed.interval_list
assets/regions_files/V7/hg37/SS_V7_hg19_probes.sorted.merged.no_chr.bed
assets/regions_files/V7/hg37/SS_V7_hg19_probes.bed.interval_list

# hg38 (UCSC chromosome nomenclature: "chr1,chr2,chr3,chrX,chrY")
assets/regions_files/V7/hg38/SS_V7_hg38_regions.sorted.merged.bed
assets/regions_files/V7/hg38/SS_V7_hg38_regions.bed.interval_list
assets/regions_files/V7/hg38/SS_V7_hg38_probes.sorted.merged.bed
assets/regions_files/V7/hg38/SS_V7_hg38_probes.bed.interval_list
```

### Available Files on S3
#### SureSelect V6
```
# hg19/hg37 (Ensemble chromosome nomenclature: "1,2,3,X,Y")
s3://zenome-ngs-data/reagents/sureselect/V6/hg37/SS_V6_hg19_regions.sorted.merged.no_chr.bed
s3://zenome-ngs-data/reagents/sureselect/V6/hg37/SS_V6_hg19_regions.bed.interval_list

# hg38 (UCSC chromosome nomenclature: "chr1,chr2,chr3,chrX,chrY")
s3://zenome-ngs-data/reagents/sureselect/V6/hg38/SS_V6_hg38_regions.sorted.merged.bed
s3://zenome-ngs-data/reagents/sureselect/V6/hg38/SS_V6_hg38_regions.bed.interval_list
```

#### SureSelect V7
```
# hg19/hg37 (Ensemble chromosome nomenclature: "1,2,3,X,Y")
s3://zenome-ngs-data/reagents/sureselect/V7/hg37/SS_V7_hg19_regions.sorted.merged.no_chr.bed
s3://zenome-ngs-data/reagents/sureselect/V7/hg37/SS_V7_hg19_regions.bed.interval_list
s3://zenome-ngs-data/reagents/sureselect/V7/hg37/SS_V7_hg19_probes.sorted.merged.no_chr.bed
s3://zenome-ngs-data/reagents/sureselect/V7/hg37/SS_V7_hg19_probes.bed.interval_list

# hg38 (UCSC chromosome nomenclature: "chr1,chr2,chr3,chrX,chrY")
s3://zenome-ngs-data/reagents/sureselect/V7/hg38/SS_V7_hg38_regions.sorted.merged.bed
s3://zenome-ngs-data/reagents/sureselect/V7/hg38/SS_V7_hg38_regions.bed.interval_list
s3://zenome-ngs-data/reagents/sureselect/V7/hg38/SS_V7_hg38_probes.sorted.merged.bed
s3://zenome-ngs-data/reagents/sureselect/V7/hg38/SS_V7_hg38_probes.bed.interval_list
```

## File description
Each file has its bed counterpart, both files should be euqally suitible for the pipeline and have identical information, although in slightly different format. Bed files are simple bed files, while interval_list files are picard-specifc files. Both files should be interchangeable, with one exception - V6 hg37 interval_list file doesn't work currently with the pipeline and gives an error. All other files work fine.

V7 files have two kinds of files - regions and probes, while V6 files have only regions. This is because V6 files only have Regions.bed file provided and no MergedProbes.bed file.

## How to use
The region files can be provided as optional input to the pipeline with `--regions` parameter:
```
nextflow run . -profile test_hg37 --regions assets/regions_files/V6/hg37/SS_V6_hg19_regions.sorted.merged.no_chr.bed
nextflow run . -profile test_hg38,yandex --accessKey <accessKey> --secretKey <secretKey> --regions s3://zenome-ngs-data/reagents/sureselect/V7/hg38/SS_V7_hg38_probes.bed.interval_list
```