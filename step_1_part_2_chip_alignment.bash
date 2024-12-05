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

if [ ! -d ${OUTPUT_DIR}/${ID}/epi_mark_alignment ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}/epi_mark_alignment \
                 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/vcf_subsets
fi

######################################################################
# New samtools fastq reconverting the epigenetic INPUT FASTQ reads from the BAM files
######################################################################
# Grab the input files form BAM and sort by name in prep for fastq conversion
$SAMTOOLS sort -n -@ ${THREADS} ${EPIGENETIC_MARK_INPUT_BAM} -O bam -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.raw_name_sorted.bam

# Converting sorted BAM temporary file into FASTQs, including those that aren't paired since I need all of them to map and see if mapping is improved
$SAMTOOLS fastq -n -1 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
                   -2 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
                   -0 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_supplementary.secondary.fastq.gz \
                   -s ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_singleton.fastq.gz \
                      ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.raw_name_sorted.bam

# Quality control
$FASTQC ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
        -o ${OUTPUT_DIR}/${ID}/stats

$FASTP --in1 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
       --in2 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
       --html ${OUTPUT_DIR}/${ID}/stats/${ID}.input.fastp.html \
       --json ${OUTPUT_DIR}/${ID}/stats/${ID}.input.fastp.json

######################################################################
# New samtools fastq reconverting the epigenetic TREATMENT FASTQ reads from the BAM files
######################################################################
# Grab the treatment files form BAM and sort by name in prep for fastq conversion
$SAMTOOLS sort -n -@ ${THREADS} ${EPIGENETIC_MARK_TREAT_BAM} -O bam -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.raw_name_sorted.bam

# Converting sorted BAM temporary file into FASTQs, including those that aren't paired since I need all of them to map and see if mapping is improved
$SAMTOOLS fastq -n -1 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
                   -2 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
                   -0 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_supplementary.secondary.fastq.gz \
                   -s ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_singleton.fastq.gz \
                      ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.raw_name_sorted.bam

# Quality control
$FASTQC ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
        -o ${OUTPUT_DIR}/${ID}/stats

$FASTP --in1 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
       --in2 ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
       --html ${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.fastp.html \
       --json ${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.fastp.json

######################################################################
# WASP PIPELINE
######################################################################
# Split the VCF by chromosome for WASP
$TABIX --list-chroms ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz | while read LINE
do
$BCFTOOLS view \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz \
          --threads $THREADS \
          --regions $LINE \
          --output-type z \
          --output ${OUTPUT_DIR}/${ID}/epi_mark_alignment/vcf_subsets/${ID}.recal_hard_filtered_phased_snvs.${LINE}.vcf.gz
done

# Create h5 files from individual phased split VCF
$WASP/snp2h5/snp2h5 --chrom     $CHROMINFO \
                    --format    vcf \
                    --haplotype ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.haplotypes.h5 \
                    --snp_index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_index.h5 \
                    --snp_tab   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_tab.h5 \
                                ${OUTPUT_DIR}/${ID}/epi_mark_alignment/vcf_subsets/*.vcf.gz

################# Alignment 1
# hg38_input
$BWA mem -M -t $THREADS \
               $REF \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.input.sam.log | \
$SAMTOOLS view -@ $THREADS -b -h -q 10 - | \
$SAMTOOLS sort -@ $THREADS - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.bam

# hg38_treatment
$BWA mem -M -t $THREADS \
               $REF \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.sam.log | \
$SAMTOOLS view -@ $THREADS -b -h -q 10 - | \
$SAMTOOLS sort -@ $THREADS - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.bam

# Pull out reads that need to be remapped to check for bias
# WASP Find Intersecting SNPs
# Input
$PYTHON $WASP/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir ${OUTPUT_DIR}/${ID}/epi_mark_alignment \
        --haplotype ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.haplotypes.h5 \
        --snp_index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_index.h5 \
        --snp_tab   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_tab.h5 \
        --samples "$($BCFTOOLS query -l ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz)" \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.bam

# Treatment
# Input
$PYTHON $WASP/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir ${OUTPUT_DIR}/${ID}/epi_mark_alignment \
        --haplotype ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.haplotypes.h5 \
        --snp_index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_index.h5 \
        --snp_tab   ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_tab.h5 \
        --samples "$($BCFTOOLS query -l ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz)" \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.bam

# Remap the reads, using same the program and options as before (Alignment 2)
# hg38_input
$BWA mem -M -t $THREADS \
               $REF \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.fq1.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.fq2.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.input.remap.sam.log | \
$SAMTOOLS view -@ $THREADS -b -h -q 10 - | \
$SAMTOOLS sort -@ $THREADS - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.bam

# hg38_treatment
$BWA mem -M -t $THREADS \
               $REF \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.fq1.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.fq2.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.treatment.remap.sam.log | \
$SAMTOOLS view -@ $THREADS -b -h -q 10 - | \
$SAMTOOLS sort -@ $THREADS - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.bam

# Use filter_remapped_reads.py to create filtered list of reads that correctly
# remap to same position
# hg38_input
$PYTHON $WASP/mapping/filter_remapped_reads.py \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.to.remap.bam \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.bam \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.keep.bam

# hg38_treatment
$PYTHON $WASP/mapping/filter_remapped_reads.py \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.to.remap.bam \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.bam \
        ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.keep.bam

# Create a merged BAM containing [1] reads that did
# not need remapping [2] filtered remapped reads
# Then sort, index, remove dups, then re-index
# hg38_input
$SAMTOOLS merge ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.bam \
                ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.keep.bam \
                ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.bam

$SAMTOOLS sort ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.bam \
            -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.bam
$SAMTOOLS index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.bam

$PYTHON $WASP/mapping/rmdup_pe.py ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.bam \
                                  ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.rmdup.bam

$SAMTOOLS sort ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.rmdup.bam \
            -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.rmdup.sorted.bam
$SAMTOOLS index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.rmdup.sorted.bam

# hg38_treatment
$SAMTOOLS merge ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.bam \
                ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.keep.bam \
                ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.bam

$SAMTOOLS sort ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.bam \
            -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.bam
$SAMTOOLS index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.bam

$PYTHON $WASP/mapping/rmdup_pe.py ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.bam \
                                  ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.rmdup.bam

$SAMTOOLS sort ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.rmdup.bam \
            -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.rmdup.sorted.bam
$SAMTOOLS index ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.rmdup.sorted.bam

######################################################################
######################################################################
######################################################################
# OLD ALIGNMENT STEPS TO PERSONALIZED REFERENCE GENOME BELOW for quality control and previous comparison purposes
# alt_input
# Mapping stuff to the base or alternate reference genome
$BWA mem -M -t $THREADS \
               ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.personalized.input.sam.log | \
$SAMTOOLS sort -@ $THREADS -O bam - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.input.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.input.bam

# alt_treatment
# Mapping stuff to the base or alternate reference genome
$BWA mem -M -t $THREADS \
               ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
               ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/${ID}.personalized.treatment.sam.log | \
$SAMTOOLS sort -@ $THREADS -O bam - -o ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.treatment.bam
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.personalized.treatment.bam

######################################################################
# Removing all other temporary files
######################################################################
rm \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.raw_name_sorted.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.raw_name_sorted.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads1.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_reads2.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_supplementary.secondary.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input_singleton.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads1.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_reads2.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_supplementary.secondary.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment_singleton.fastq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.haplotypes.h5 \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_index.h5 \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.snp_tab.h5 \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.to.remap.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.fq1.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.fq2.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.single.fq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.bam.bai \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.remap.keep.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.bam.bai \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.input.keep.merged.sorted.rmdup.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.to.remap.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.fq1.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.fq2.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.single.fq.gz \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.bam.bai \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.remap.keep.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.bam \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.bam.bai \
${OUTPUT_DIR}/${ID}/epi_mark_alignment/${ID}.treatment.keep.merged.sorted.rmdup.bam

rm ${OUTPUT_DIR}/${ID}/epi_mark_alignment/vcf_subsets/*
rmdir ${OUTPUT_DIR}/${ID}/epi_mark_alignment/vcf_subsets
######################################################################
# End message
######################################################################

echo "Done with $(basename ${OUTPUT_DIR})/${ID} Stage 2: Epigenetic Mark Alignment"
