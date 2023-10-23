directory = $1
find "$directory" -name "Aligned.sortedByCoord.out.bam" | parallel 'samtools view -F 0x40 {} | cut -f1 | sort | uniq' | sort | uniq -c >> size.txt
