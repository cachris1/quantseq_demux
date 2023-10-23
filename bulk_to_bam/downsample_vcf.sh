# Corrected and proofread commands:

# Save header lines to a file
bcftools view --h "$1" > "${1}_subsample.vcf"

# Select variants (without header) and shuffle 1200 lines, then append to the VCF file
bcftools view -H "$1" | shuf -n "$2" >> "${1}_subsample.vcf"

# Sort the VCF file and save the sorted version
bcftools sort -o "${1}_subsample_sort.vcf" "${1}_subsample.vcf"

# Create a compressed VCF file (gzip)
bcftools view -o "${1}_subsample_sort.vcf.gz" "${1}_subsample_sort.vcf"

# Create an index for the compressed VCF file
bcftools index "${1}_subsample_sort.vcf.gz"

# remove intermediate files
rm "${1}_subsample.vcf"
rm "${1}_subsample_sort.vcf"





