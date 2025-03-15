#!/bin/bash

GENOMEIDX1="/home/ngs-linux/workshop_CUP/reference_genome/NC_045512.2.fasta"

BASEDIRDATA="$(pwd)"
reads1=(${BASEDIRDATA}/*_P1.fastq)
reads1=("${reads1[@]##*/}")
#echo ${reads1[@]}
reads2=("${reads1[@]/_P1./_P2.}")

echo "***Following SAMPLES are submitted for Germline variant calling***"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
{
    #sample="${reads1[$i]%%.*}"
    sample="${reads1[$i]%_P1*}"
    sample="${sample%P1*}"
    echo "$sample" 
    }
    done

for ((i=0; i<=${#reads1[@]}-1; i++ )); do
{
    sample="${reads1[$i]%_P1*}"
    sample="${sample%_P1*}"
    stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "$stime Processing sample: $sample"

bwa mem -t 4 $GENOMEIDX1 ${BASEDIRDATA}/${sample}_P1.fastq ${BASEDIRDATA}/${sample}_P2.fastq -o ${sample}.sam

gatk SortSam --INPUT ${BASEDIRDATA}/${sample}.sam --OUTPUT ${sample}.bam --SORT_ORDER coordinate

gatk BuildBamIndex --INPUT ${BASEDIRDATA}/${sample}.bam 

gatk CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE $GENOMEIDX1 --INPUT ${BASEDIRDATA}/${sample}.bam --OUTPUT ${sample}_alignment_metrics.txt --VALIDATION_STRINGENCY LENIENT

gatk CollectInsertSizeMetrics --INPUT ${BASEDIRDATA}/${sample}.bam --OUTPUT ${sample}_metrics.txt --Histogram_FILE ${sample}_bwa_histogram.pdf

gatk AddOrReplaceReadGroups --INPUT ${BASEDIRDATA}/${sample}.bam --OUTPUT ${sample}_group.bam --RGLB lib1 --RGPL illumina --RGPU NONE --RGSM ${sample}

gatk MarkDuplicates --INPUT ${BASEDIRDATA}/${sample}_group.bam --OUTPUT ${sample}_dedup_reads.bam --METRICS_FILE ${BASEDIRDATA}/${sample}_metrics.txt -AS true --VALIDATION_STRINGENCY LENIENT

gatk BuildBamIndex --INPUT ${BASEDIRDATA}/${sample}_dedup_reads.bam

gatk HaplotypeCaller -R $GENOMEIDX1 -I ${sample}_dedup_reads.bam -ploidy 1 -O ${sample}.vcf

java -jar /home/ngs-linux/workshop_CUP/software/snpEff/snpEff.jar -v NC_045512.2  ${sample}.vcf >  ${sample}_ann.vcf
java -Xmx64G -jar /home/ngs-linux/workshop_CUP/software/snpEff/SnpSift.jar extractFields ${sample}_ann.vcf "CHROM" "POS" "REF" "ALT" "QUAL" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].FEATURE" "ANN[0].BIOTYPE" "ANN[0].AA" "GEN[*].GT" > ${sample}_ann.txt

}
done

