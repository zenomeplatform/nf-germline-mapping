## Commands to obtain data

The initial V7 files S31285117_Regions.bed and S31285117_MergedProbes.bed must be obrtained from SureSelect website.

Then, run these commands:

```
# Remove first two rows
cat S31285117_Regions.bed | tail -n +3 > SS_V7_hg19_raw.bed
#Sort
bedtools sort -i SS_V7_hg19_raw.bed > SS_V7_hg19_regions.sorted.bed
# Merge overlapping regions
bedtools merge -i SS_V7_hg19_regions.sorted.bed > SS_V7_hg19_regions.sorted.merged.bed
# Remove "chr" from chromsome names to make bed compatible with Ensembl-type chromosome notation "1,2,3,X,Y"
sed 's/chr//g' SS_V7_hg19_regions.sorted.merged.bed > SS_V7_hg19_regions.sorted.merged.no_chr.bed

# Get the corresponding genome dict file
aws s3 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict . --no-sign-request

# Create a picard-style interval list file
picard BedToIntervalList I=SS_V7_hg19_regions.sorted.merged.no_chr.bed O=SS_V7_hg19_regions.bed.interval_list SD=human_g1k_v37_decoy.dict
```

Similar thing for MergedProbes bed file, with additional step of selecting only first three columns
```
cat S31285117_MergedProbes.bed | tail -n +3 > SS_V7_hg19_raw_probes.bed
cut -f 1-3 SS_V7_hg19_raw_probes.bed > SS_V7_hg19_probes.bed

bedtools sort -i  SS_V7_hg19_probes.bed > SS_V7_hg19_probes.sorted.bed
bedtools merge -i SS_V7_hg19_probes.sorted.bed > SS_V7_hg19_probes.sorted.merged.bed
sed 's/chr//g' SS_V7_hg19_probes.sorted.merged.bed > SS_V7_hg19_probes.sorted.merged.no_chr.bed

picard BedToIntervalList I=SS_V7_hg19_probes.sorted.merged.no_chr.bed O=SS_V7_hg19_probes.bed.interval_list SD=human_g1k_v37_decoy.dict
```