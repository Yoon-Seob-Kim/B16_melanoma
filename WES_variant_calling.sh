#!/bin/bash
BWA=/data4/kysbbubbu/tools/bwa-0.7.17
GATK=/data4/kysbbubbu/tools/gatk-4.1.6.0
DB=/data4/kysbbubbu/genomicdb/mouse
BAM_DIR=./02_BAM
MT_DIR=./04_MT
item1=C57BL6J
item2=B16F0 
$GATK/gatk --java-options "-Xmx4g" Mutect2 -R $DB/GRCm38.primary_assembly.genome.fa -I $BAM_DIR/$item1\_mm10.bam -I $BAM_DIR/$item2\_mm10.bam -normal $item1 -tumor $item2 --intervals $DB/SureSelece_Mouse_Exom_MM10.list -O $MT_DIR/$item2\.vcf.gz
$GATK/gatk --java-options "-Xmx4g" FilterMutectCalls -V $MT_DIR/$item2\.vcf.gz -R $DB/GRCm38.primary_assembly.genome.fa -O $MT_DIR/$item2\_filt.vcf
awk -F '\t' '{if($1== "#CHROM") print; else if($7 == "PASS") print}' $MT_DIR/$item2\_filt.vcf > $item2\.vcf
