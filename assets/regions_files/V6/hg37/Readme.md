## Commands to obtain data

The initial file v6 S07604514_Regions.bed must be obrtained from SureSelect website.

Then, run these commands:
```
# Remove first two rows
cat S07604514_Regions.bed | tail -n +3 > SS_V6_hg19_raw.bed
# Sort
bedtools sort -i SS_V6_hg19_raw.bed > SS_V6_hg19_regions.sorted.bed
# Merge overlapping regions
bedtools merge -i SS_V6_hg19_regions.sorted.bed > SS_V6_hg19_regions.sorted.merged.bed
# Remove "chr" from chromsome names to make bed compatible with Ensembl-type chromosome notation "1,2,3,X,Y"
sed 's/chr//g' SS_V6_hg19_regions.sorted.merged.bed > SS_V6_hg19_regions.sorted.merged.no_chr.bed

# Get the corresponding genome dict file
aws s3 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict . --no-sign-request

# Create a picard-style interval list file
picard BedToIntervalList I=SS_V6_hg19_regions.sorted.merged.no_chr.bed O=SS_V6_hg19_regions.bed.interval_list SD=human_g1k_v37_decoy.dict
```