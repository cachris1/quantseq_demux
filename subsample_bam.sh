#save header lines to a file
samtools view -H "$1" > "${1}_subsample.sam"

# Select variants (without header) and shuffle 1200 lines, then append to the VCF file
samtools view "$1" | shuf -n "$2" >> "${1}_subsample.sam"

# Sort the VCF file and save the sorted version
samtools sort -o "${1}_subsample_sort.sam" "${1}_subsample.sam"

# Create a compressed VCF file (gzip)
samtools view -o "${1}_subsample_sort.bam" "${1}_subsample_sort.sam"

# Create an index for the compressed VCF file
samtools index "${1}_subsample_sort.bam"

# remove intermediate files
rm "${1}_subsample.sam"
rm "${1}_subsample_sort.sam"

