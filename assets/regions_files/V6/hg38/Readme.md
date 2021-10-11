## Commands to obtain data

The initial file v6 S07604514_Regions.bed must be obrtained from SureSelect website.

Then, run these commands:
# Remove first two rows
cat S07604514_Regions.bed | tail -n +3 > SS_V6_hg38_raw.bed
#Sort
bedtools sort -i SS_V6_hg38_raw.bed > SS_V6_hg38_regions.sorted.bed
#Merge overlapping regions
bedtools merge -i SS_V6_hg38_regions.sorted.bed > SS_V6_hg38_regions.sorted.merged.bed

#Get the corresponding genome dict file
aws s3 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict . --no-sign-request

# Create a picard-style interval list file
picard BedToIntervalList I=SS_V6_hg38_regions.sorted.merged.bed O=SS_V6_hg38_regions.bed.interval_list SD=Homo_sapiens_assembly38.dict


