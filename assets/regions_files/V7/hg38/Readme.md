## Commands to obtain data

The initial file v7 S31285117_Regions.bed must be obrtained from SureSelect website.

Then, run these commands:
# Remove first two rows
cat S31285117_Regions.bed | tail -n +3 > SS_V7_hg38_raw.bed
#Sort
bedtools sort -i SS_V7_hg38_raw.bed > SS_V7_hg38_regions.sorted.bed
#Merge overlapping regions
bedtools merge -i SS_V7_hg38_regions.sorted.bed > SS_V7_hg38_regions.sorted.merged.bed

#Get the corresponding genome dict file
aws s3 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict . --no-sign-request

# Create a picard-style interval list file
picard BedToIntervalList I=SS_V7_hg38_regions.sorted.merged.bed O=SS_V7_hg38_regions.bed.interval_list SD=Homo_sapiens_assembly38.dict


## Similar thing for MergedProbes bed file, with additional step of selecting only first three columns
cat S31285117_MergedProbes.bed | tail -n +3 > SS_V7_hg38_raw_probes.bed
cut -f 1-3 SS_V7_hg38_raw_probes.bed > SS_V7_hg38_probes.bed

bedtools sort -i  SS_V7_hg38_probes.bed > SS_V7_hg38_probes.sorted.bed
bedtools merge -i SS_V7_hg38_probes.sorted.bed > SS_V7_hg38_probes.sorted.merged.bed

picard BedToIntervalList I=SS_V7_hg38_probes.sorted.merged.bed O=SS_V7_hg38_probes.bed.interval_list SD=Homo_sapiens_assembly38.dict