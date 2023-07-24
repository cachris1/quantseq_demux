import subprocess 
import os
import system

# a file to align a set of fastqs using STAR.  For paired-end reads.  Takes in 

def call_star(fastqs):
    cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3'
    cmd += ' --twopassMode Basic'
    cmd += ' --runThreadN ' + \
        "15" + ' --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate '
    cmd += '--alignSJDBoverhangMin ' + "1" + \
        ' --alignIntronMax 299999 --genomeDir ' \
        + "/home/genomes/hg38/star_hg38/" + ' --sjdbGTFfile ' + "/home/genomes/hg38/star_hg38/Homo_sapiens.GRCh38.99.gtf"  
    cmd += ' --outFileNamePrefix ' + "./star_out/" + fastqs[0].split("/")[-2]  + \
        " --readFilesIn " + fastqs[0] + " " + fastqs[1] + ' --readFilesCommand zcat'
    print(cmd)
    subprocess.run(cmd, shell=True)


def main(fastq_dir):
   for sample in os.listdir(dir1):
        if sample.endswith(".fastq"):
          fastqs = []
          for fastq in os.listdir(dir1 + sample):
              fastqs.append(dir1+sample+"/"+ fastq)
          call_star(fastqs)

# __name__
if __name__=="__main__":
    main(fastq_dir, out_dir)


