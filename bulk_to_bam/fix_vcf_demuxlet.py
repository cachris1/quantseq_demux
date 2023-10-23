import os
import subprocess
import pandas as pd

subprocess.call("bcftools query -l merged_vcf.vcf.gz | while read sample; do
bcftools query -f '%CHROM\t%POS\n' -i "GT[\"${sample}\"] != \"./.\"" -o positions_${sample}.txt merger_vcf.vcf.gz
done")
for sample in samples:
    subprocess.call("bcftools query -f 'CHROM\tPOS\n' -i "GT[\"${sample}\"] != \"./.\"" -o positions_$" + sample + ".txt merged_vcf.vcf.gz")
locs = pd.DataFrame()
for sample in samples:
    pd.read_csv()
    