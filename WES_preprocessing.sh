#!/bin/bash

#tool_path
BWA=/data4/kysbbubbu/tools/bwa-0.7.17
GATK=/data4/kysbbubbu/tools/gatk-4.1.6.0
DB=/data4/kysbbubbu/genomicdb/mouse
QUALIMAP=/data4/kysbbubbu/tools/qualimap_v2.2.1

##OUTPUT
BAM_DIR=./02_BAM
BAMQC_DIR=./03_BAMQC
MT_DIR=./04_MT
ANNO_DIR=./05_Annovar
for i in B16F0 B16F10 C57BL6J
do
$BWA/bwa mem -t 4 $DB/GRCm38.primary_assembly.genome.fa $i\_1.fastq.gz $i\_2.fastq.gz > $i\.sam
$GATK/gatk --java-options "-Xmx4g" SortSam -I=$i\.sam -O=$i\.bam -SO=coordinate --VALIDATION_STRINGENCY=SILENT --CREATE_INDEX=true
rm $i\.sam
$GATK/gatk --java-options "-Xmx4g" MarkDuplicates -I=$i\.bam -O=$i\_dup.bam -M=$i\_dup.met --VALIDATION_STRINGENCY=SILENT --REMOVE_DUPLICATES=true
rm $i\.bam $i\.bai
$GATK/gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I=$i\_dup.bam -O=$i\_dupRemoved.bam -SO=coordinate --VALIDATION_STRINGENCY=SILENT -ID=$i -LB=$i -PL=Illumina -PU=$i -SM=$i --CREATE_INDEX TRUE
rm $i\_dup.bam $i\_dup.met
$GATK/gatk --java-options "-Xmx4g" BaseRecalibrator -R $DB/GRCm38.primary_assembly.genome.fa -I $i\_dupRemoved.bam --known-sites $DB/C57BL_6NJ.mgp.v5.snps.dbSNP142_sorted.vcf --known-sites $DB/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed_sorted.vcf -O $i\_recal_data.table
$GATK/gatk --java-options "-Xmx4g" ApplyBQSR -R $DB/GRCm38.primary_assembly.genome.fa -I $i\_dupRemoved.bam -bqsr $i\_recal_data.table -O $BAM_DIR/$i\_mm10.bam
rm $i\_dupRemoved.bam $i\_dupRemoved.bai $i\_recal_data.table
$QUALIMAP/qualimap bamqc -bam $BAM_DIR/$i\_mm10.bam -gff $DB/SureSelece_Mouse_Exom_MM10.bed -outdir $BAMQC_DIR/$i -c --java-mem-size=4G
done



