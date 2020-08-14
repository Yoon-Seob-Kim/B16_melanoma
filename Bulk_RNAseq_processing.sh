#!/bin/bash
DATA=/home/kysbbubbu/B16
TRIM=$DATA/02_TRIM
BAM=$DATA/03_BAM
STRING=$DATA/04_STRINGTIE
TRIMOMATIC=/data/MRC1_data4/kysbbubbu/tools/Trimmomatic-0.39
STRINGTIE=/data/MRC1_data4/kysbbubbu/tools/stringtie-2.1.4.Linux_x86_64
for T1 in B16F0 B16F10
do
java -Xmx4g -jar $TRIMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 $DATA/01_FASTQ/$T1\_1.fastq.gz $DATA/01_FASTQ/$T1\_2.fastq.gz $TRIM/$T1\_forward_paired.fastq.gz $TRIM/$T1\_forward_unpaired.fastq.gz $TRIM/$T1\_reverse_paired.fastq.gz $TRIM/$T1\_reverse_unpaired.fastq.gz ILLUMINACLIP:$TRIMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
hisat2 --dta -p 4 -x /data/MRC1_data4/kysbbubbu/genomicdb/grcm38_snp_tran/genome_snp_tran -1 $TRIM/$T1\_forward_paired.fastq.gz -2 $TRIM/$T1\_reverse_paired.fastq.gz -S $BAM/$T1\.sam
rm $TRIM/$T1\_forward_paired.fastq.gz $TRIM/$T1\_reverse_paired.fastq.gz
samtools sort -o $BAM/$T1\.sorted.bam $BAM/$T1\.sam
rm $BAM/$T1\.sam
$STRINGTIE/stringtie $BAM/$T1\.sorted.bam -p 4 -e -B -G /data/MRC1_data4/kysbbubbu/genomicdb/cellranger/refdata-cellranger-mm10-3.0.0/genes/genes.gtf -o $STRING/$T1/$T1\.gtf -A $STRING/$T1/$T1\.txt
done


