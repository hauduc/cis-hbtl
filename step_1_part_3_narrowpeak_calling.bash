#!/bin/bash

######################################################################
# Safe Mode Parameters
######################################################################
set -e
set -u
set -o pipefail

######################################################################
# Set user-inputted arguments
######################################################################
ID=$1
OUTPUT_DIR=$2
THREADS=$3
REF=$4
BAM=$5
EPIGENETIC_MARK_INPUT_BAM=$6 # CEMT_155 = /projects/epigenomics_assembly/jsteif/ChIP-seq/CEMT/Breast_normal/bams/Input/E00534.CEMT_155.Input.hg38.merged.E00534_2_lanes_dupsFlagged.bam
EPIGENETIC_MARK_TREAT_BAM=$7 # CEMT_155 = /projects/epigenomics_assembly/jsteif/ChIP-seq/CEMT/Breast_normal/bams/H3K4me3/E00547.CEMT_155.H3K4me3.hg38.merged.E00547_2_lanes_dupsFlagged.bam

######################################################################
# Additional resources for hg38_no_alt.fa
######################################################################
CHROMINFO="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/refs/chrom_info/hg38_no_alt.chromsizes.txt"
BLACKLIST="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/blacklists/hg38-blacklist.v2.bed"
DBSNP="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/dbSNP/all/00-All_chr.vcf.gz"
DBSNP_SNPS="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/dbSNP/snps/00-All_chr.snps.vcf.gz"
ALL_1KG_PHASE3_INDELS="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/callsets/1kg/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels_fixed.vcf.gz"
MILLS_AND_1KG_GOLD_INDELS="/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/callsets/mills_indels/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

######################################################################
# Set software packages employed & define dependencies
######################################################################
# Main tools
BWA="/gsc/software/linux-x86_64-centos7/bwa-0.7.6a/bwa"
SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.10/bin/samtools"
GATK="/home/ahauduc/software/gatk-4.1.4.1/gatk"
PICARD="/gsc/software/linux-x86_64-centos7/picard-tools-2.4.1/picard.jar"
BEDTOOLS="/gsc/software/linux-x86_64-centos7/bedtools-2.27.1/bin/bedtools"
BCFTOOLS="/gsc/software/linux-x86_64-centos7/bcftools-1.11/bin/bcftools"
TABIX="/home/ahauduc/anaconda3/bin/tabix"
PYTHON="/home/ahauduc/software/Python-3.8.6/python"
MACS2="/home/ahauduc/software/MACS2-2.2.7.1/bin/macs2"
bedGraphToBigWig="/home/ahauduc/software/kentUtils/bedGraphToBigWig"
HOMER="/home/ahauduc/software/homer"

# Quality control
FASTQC="/home/ahauduc/anaconda3/bin/fastqc"
FASTP="/gsc/software/linux-x86_64-centos7/fastp-0.21.0/bin/fastp"

# Virtualization
JAVA="/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/bin/java"

# Phasing
WHATSHAP="/gsc/software/linux-x86_64-centos7/whatshap-1.0/bin/whatshap"
G2GTOOLS="/home/ahauduc/anaconda3/envs/g2gtools/bin/g2gtools"

# WASP
WASP="/home/ahauduc/software/WASP"

# Need R 3.6.3 with certain packages for GATK AnalyzeCovariates installed:
# R> install.packages(c("gsalib", "ggplot2", "reshape", "gplots"))

######################################################################
# Making the directories to do the project work in
######################################################################

if [ ! -d ${OUTPUT_DIR}/${ID} ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}
fi

if [ ! -d ${OUTPUT_DIR}/${ID}/peak_calls ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}/peak_calls
fi

######################################################################
# Peak calling
######################################################################
# hg38
$MACS2 callpeak --control   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.bam \
                --treatment ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.bam \
                --format BAMPE \
                --gsize hs \
                --keep-dup all \
                --bdg \
                --qvalue 0.01 \
                --outdir ${OUTPUT_DIR}/${ID}/peak_calls \
                --name ${ID}.hg38

# Convert then remove bedgraphs
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_treat_pileup.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_treat_pileup.bw
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_control_lambda.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_control_lambda.bw
rm \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_treat_pileup.bdg \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.hg38_control_lambda.bdg

# personalized
$MACS2 callpeak --control   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.input.bam \
                --treatment ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.treatment.bam \
                --format BAMPE \
                --gsize hs \
                --keep-dup all \
                --bdg \
                --qvalue 0.01 \
                --outdir ${OUTPUT_DIR}/${ID}/peak_calls \
                --name ${ID}.personalized

# Convert then remove bedgraphs
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_treat_pileup.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_treat_pileup.bw
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_control_lambda.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_control_lambda.bw
rm \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_treat_pileup.bdg \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.personalized_control_lambda.bdg

######################################################################
# WASP peak calling
######################################################################
$MACS2 callpeak --control   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.rmdup.sorted.bam \
                --treatment ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.rmdup.sorted.bam \
                --format BAMPE \
                --gsize hs \
                --keep-dup all \
                --bdg \
                --qvalue 0.01 \
                --outdir ${OUTPUT_DIR}/${ID}/peak_calls \
                --name ${ID}.wasp

# Convert then remove bedgraphs
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_treat_pileup.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_treat_pileup.bw
$bedGraphToBigWig ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_control_lambda.bdg \
                  $CHROMINFO \
                  ${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_control_lambda.bw
rm \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_treat_pileup.bdg \
${OUTPUT_DIR}/${ID}/peak_calls/${ID}.wasp_control_lambda.bdg

######################################################################
# Removing regions that appear in ENCODE blacklist
######################################################################
for MACS_OUTPUT_FILE in ${OUTPUT_DIR}/${ID}/peak_calls/*{narrowPeak,summits.bed}
do
$BEDTOOLS intersect -v -wa -a $MACS_OUTPUT_FILE -b $BLACKLIST > ${MACS_OUTPUT_FILE}.filtered.bed
done

######################################################################
# Creating set of HOMER peak differentiation outputs
######################################################################
# Change dirs to get misc outputs in right place
pushd ${OUTPUT_DIR}/${ID}/peak_calls
# Merge peaks within the average distance given
# hg38 vs personalized
for DISTANCE in $(seq 0 100 200)
do
$HOMER/bin/mergePeaks -d $DISTANCE \
                      ${ID}.hg38_peaks.narrowPeak.filtered.bed \
                      ${ID}.personalized_peaks.narrowPeak.filtered.bed \
                      -prefix hg38_vs_pers_distance_${DISTANCE} \
                      -matrix hg38_vs_pers_distance_${DISTANCE} \
                      2> hg38_vs_pers_distance_${DISTANCE}.stderr
done

# hg38 vs WASP
for DISTANCE in $(seq 0 100 200)
do
$HOMER/bin/mergePeaks -d $DISTANCE \
                      ${ID}.hg38_peaks.narrowPeak.filtered.bed \
                      ${ID}.wasp_peaks.narrowPeak.filtered.bed \
                      -prefix hg38_vs_wasp_${DISTANCE} \
                      -matrix hg38_vs_wasp_${DISTANCE} \
                      2> hg38_vs_wasp_distance_${DISTANCE}.stderr
done

# personalized vs WASP
for DISTANCE in $(seq 0 100 200)
do
$HOMER/bin/mergePeaks -d $DISTANCE \
                      ${ID}.personalized_peaks.narrowPeak.filtered.bed \
                      ${ID}.wasp_peaks.narrowPeak.filtered.bed \
                      -prefix pers_vs_wasp_${DISTANCE} \
                      -matrix pers_vs_wasp_${DISTANCE} \
                      2> pers_vs_wasp_distance_${DISTANCE}.stderr
done
popd

# Move all peak overlap files to their own folder
mkdir ${OUTPUT_DIR}/${ID}/peak_calls/comparisons
mv ${OUTPUT_DIR}/${ID}/peak_calls/{hg38_,pers_}* ${OUTPUT_DIR}/${ID}/peak_calls/comparisons

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Optional! Remove epigenetic alignments after calling
rm ${OUTPUT_DIR}/${ID}/epi_mark_alignment/*
rmdir ${OUTPUT_DIR}/${ID}/epi_mark_alignment

# Stats for epigenetic alignments
rm \
${OUTPUT_DIR}/${ID}/stats/${ID}.personalized.treatment.sam.log \
${OUTPUT_DIR}/${ID}/stats/${ID}.personalized.input.sam.log \
${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.remap.sam.log \
${OUTPUT_DIR}/${ID}/stats/${ID}.input.remap.sam.log \
${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.sam.log \
${OUTPUT_DIR}/${ID}/stats/${ID}.input.sam.log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######################################################################
# End message
######################################################################

echo "Done with $(basename ${OUTPUT_DIR})/${ID} Stage 3: NarrowPeak Calling"
