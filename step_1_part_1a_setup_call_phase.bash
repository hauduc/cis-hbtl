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

if [ ! -d ${OUTPUT_DIR}/${ID}/wgs ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}/wgs \
                 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets \
                 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments \
                 ${OUTPUT_DIR}/${ID}/stats
fi

######################################################################
# Check if reference is indexed, and index if not
######################################################################

if [ -e ${REF}.fai ] 
then
    :
else
    $BWA index $REF
fi

######################################################################
# Processing BAM/CRAM file into WGS FASTQ files for alternate reference genome creation
######################################################################
$SAMTOOLS sort -n -@ $THREADS -O bam $BAM -o ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.raw_name_sorted.bam

# Split here and loop over all resulting BAMs from split
$SAMTOOLS split -@ $THREADS -f "${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/%!.%." ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.raw_name_sorted.bam

# Start loop over all splitted BAM files for alignment
for RG_SPLIT_FILE in ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/*.bam; do

# Converting sorted BAM temporary file into FASTQs, including those that aren't paired since I need all of them to map and see if mapping is improved
# removed -t here to not preserve read group names in fastq (should be unnecessary)
$SAMTOOLS fastq -n -@ $THREADS -1 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads1.fastq.gz \
                               -2 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads2.fastq.gz \
                               -0 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).supplementary.secondary.fastq.gz \
                               -s ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).singleton.fastq.gz \
                                  $RG_SPLIT_FILE

$FASTQC ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads1.fastq.gz \
        ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads2.fastq.gz \
        -o ${OUTPUT_DIR}/${ID}/stats

$FASTP --in1 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads1.fastq.gz \
       --in2 ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads2.fastq.gz \
       --html ${OUTPUT_DIR}/${ID}/stats/$(basename -s .bam ${RG_SPLIT_FILE})_fastp.html \
       --json ${OUTPUT_DIR}/${ID}/stats/$(basename -s .bam ${RG_SPLIT_FILE})_fastp.json

# Align each read group to hg38
$BWA mem -M -t $THREADS \
               $REF \
               ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads1.fastq.gz \
               ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/$(basename -s .bam ${RG_SPLIT_FILE}).reads2.fastq.gz \
               2> ${OUTPUT_DIR}/${ID}/stats/$(basename -s .bam ${RG_SPLIT_FILE}).sam.log | \
$SAMTOOLS sort -@ $THREADS -O bam - -o ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments/$(basename ${RG_SPLIT_FILE})
$SAMTOOLS index -@ $THREADS ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments/$(basename ${RG_SPLIT_FILE})
done

# Merge, giving each read group their original reads group name
$SAMTOOLS merge -@ $THREADS \
                -rh <($SAMTOOLS view -H ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.raw_name_sorted.bam | grep "^@RG") \
                ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.bam \
                ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments/*.bam

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Remove name-sorted BAM
rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.raw_name_sorted.bam
# Remove read group subsets
rm ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments/*
rmdir ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/unmerged_alignments
rm ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets/*
rmdir ${OUTPUT_DIR}/${ID}/wgs/read_group_subsets
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Mark dups then sort, then create bai file
$GATK --java-options "-Xmx80G" MarkDuplicatesSpark \
                                --spark-master local[${THREADS}] \
                                -I ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.bam \
                                -M ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.dedup_metrics_hg38.txt \
                                -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Remove unsorted un-rmduped bam
rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.bam
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######################################################################
# Base recalibration round 1
######################################################################
$GATK --java-options "-Xmx80G" BaseRecalibrator \
                               -R $REF \
                               -I ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam \
                               --known-sites $DBSNP_SNPS \
                               --known-sites $ALL_1KG_PHASE3_INDELS \
                               --known-sites $MILLS_AND_1KG_GOLD_INDELS \
                               -O ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.recal_data.table

$GATK --java-options "-Xmx80G" ApplyBQSR \
                               -R $REF \
                               -I ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam \
                               -bqsr ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.recal_data.table \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Remove unrecalibrated BAM and indices
rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam.bai \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.bam.sbi
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######################################################################
# Base recalibration round 2 (for information only)
######################################################################
$GATK --java-options "-Xmx80G" BaseRecalibrator \
                               -R $REF \
                               -I ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
                               --known-sites $DBSNP_SNPS \
                               --known-sites $ALL_1KG_PHASE3_INDELS \
                               --known-sites $MILLS_AND_1KG_GOLD_INDELS \
                               -O ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.post_recal_data.table

# Need R 3.6.3 with certain packages for GATK AnalyzeCovariates installed:
# R> install.packages(c("gsalib", "ggplot2", "reshape", "gplots"))
# Changing your path temporarily to be pointing to R 3.6.3 instead of R 4.1.0 for AnalyzeCovariates
export PATH="/home/ahauduc/.local/bin:/home/ahauduc/software/Python-3.8.6:/home/ahauduc/anaconda3/condabin:/data/bin:/data/bin/blat:/data/scripts:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/home/ahauduc/bin:/gsc/software/linux-x86_64-centos7/samtools-1.10/bin:/gsc/software/linux-x86_64-centos7/bedtools-2.27.1/bin:/gsc/software/linux-x86_64-centos7/bedops-2.4.35/bin:/gsc/software/linux-x86_64-centos7/bcftools-1.13/bin:/home/ahauduc/software/homer/bin:/home/ahauduc/software/MACS2-2.2.7.1/bin:/gsc/software/linux-x86_64-centos7/bwa-0.7.6a:/gsc/software/linux-x86_64-centos7/vcftools-0.1.17/bin:/home/ahauduc/software/gatk-4.1.4.1:/gsc/software/linux-x86_64-centos7/jdk1.8.0_172/bin:/home/ahauduc/software/bart/bart_v2.0/bin:/gsc/software/linux-x86_64-centos7/R-3.6.3/lib64/R/bin:/home/ahauduc/anaconda3/bin:/home/rislam/anaconda2/bin"

$GATK --java-options "-Xmx80G" AnalyzeCovariates \
                               -before ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.recal_data.table \
                               -after  ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.post_recal_data.table \
                               -plots  ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.recalibration_plots.pdf

######################################################################
# Call variants on recalibrated reads
######################################################################
$GATK --java-options "-Xmx80G" HaplotypeCaller \
                               -R $REF \
                               --dbsnp $DBSNP \
                               -I ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal.vcf.gz

# Select SNVs and indels
# SNVs
$GATK --java-options "-Xmx80G" SelectVariants \
                               --select-type-to-include SNP \
                               -R $REF \
                               -V ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal.vcf.gz \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_snvs.vcf.gz

# Indels
$GATK --java-options "-Xmx80G" SelectVariants \
                               --select-type-to-include INDEL \
                               -R $REF \
                               -V ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal.vcf.gz \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_indels.vcf.gz

######################################################################
# Recalibrated variant filtering
######################################################################
# SNVs
$GATK --java-options "-Xmx80G" VariantFiltration \
                               -R $REF \
                               -V ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_snvs.vcf.gz \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_filtered_snvs.vcf.gz \
                               -filter-name "QD_filter" -filter "QD < 2.0" \
                               -filter-name "FS_filter" -filter "FS > 60.0" \
                               -filter-name "MQ_filter" -filter "MQ < 40.0" \
                               -filter-name "SOR_filter" -filter "SOR > 4.0" \
                               -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
                               -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Indels
$GATK --java-options "-Xmx80G" VariantFiltration \
                               -R $REF \
                               -V ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_indels.vcf.gz \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_filtered_indels.vcf.gz \
                               -filter-name "QD_filter" -filter "QD < 2.0" \
                               -filter-name "FS_filter" -filter "FS > 200.0" \
                               -filter-name "SOR_filter" -filter "SOR > 10.0"

######################################################################
# Creating the personalized, haploid, SNV-only FASTA reference
######################################################################
# Temporarily hard-filter everything for GATK FastaAlternateReferenceMaker
$BCFTOOLS view \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_filtered_snvs.vcf.gz \
          --apply-filters PASS \
          --output-type z \
          > ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz

# Index hard filtered output
$BCFTOOLS index --threads $THREADS --tbi ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz

# Create haploid alternate reference
$GATK --java-options "-Xmx80G" FastaAlternateReferenceMaker \
                               -R $REF \
                               -V ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz \
                               -O ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa

# Change chromosome names of personalized reference genome in-place
cat $REF | grep "^>" | awk '{print $1}' > ${OUTPUT_DIR}/${ID}/wgs/hg38_fasta_headers.txt
cat ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa | grep "^>" > ${OUTPUT_DIR}/${ID}/wgs/pers_fasta_headers.txt

awk '
    FILENAME == ARGV[1] { listA[$1] = FNR; next }
    FILENAME == ARGV[2] { listB[FNR] = $1; next }
    {
        for (i = 1; i <= NF; i++) {
            if ($i in listA) {
                $i = listB[listA[$i]]
            }
        }
        print
    }
' ${OUTPUT_DIR}/${ID}/wgs/pers_fasta_headers.txt \
  ${OUTPUT_DIR}/${ID}/wgs/hg38_fasta_headers.txt \
  ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa > ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa.new
mv ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa.new ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa

# Reindexing personalized reference genome
$SAMTOOLS faidx ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa

rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.dict
$JAVA -jar -Xmx80g $PICARD CreateSequenceDictionary \
                           R=${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa \
                           O=${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.dict

$BWA index ${OUTPUT_DIR}/${ID}/wgs/${ID}.hap.snv.fa

######################################################################
# Creating the phased, diploid, SNV-only FASTA paternal & maternal references
######################################################################
# Phase the genome using reads-based phasing, supplying it with hard-filtered passing reads
$WHATSHAP phase \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
          --reference=${REF} \
          -o ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz

# Index output
$BCFTOOLS index --threads $THREADS --tbi ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz

# Get the name of the sample & put that into the strain 
$G2GTOOLS patch \
          --input $REF \
          --vcf ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz \
          --output ${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.fa \
          --strain "$($BCFTOOLS query -l ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_phased_snvs.vcf.gz)" \
          --diploid \
          --pass

# Index new references
$SAMTOOLS faidx ${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.l.fa
$SAMTOOLS faidx ${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.r.fa

$BWA index ${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.l.fa
$BWA index ${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.r.fa

$JAVA -jar -Xmx80g $PICARD CreateSequenceDictionary \
                           R=${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.l.fa \
                           O=${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.l.dict
$JAVA -jar -Xmx80g $PICARD CreateSequenceDictionary \
                           R=${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.r.fa \
                           O=${OUTPUT_DIR}/${ID}/wgs/${ID}.dip_snv.r.dict

######################################################################
# Remove remaining temporary files
######################################################################
# Keep the BAMs
# rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
#    ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bai

# Remove unnecessary VCFs
rm ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal.vcf.gz \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal.vcf.gz.tbi \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_snvs.vcf.gz \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_snvs.vcf.gz.tbi \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_indels.vcf.gz \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_indels.vcf.gz.tbi \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz \
   ${OUTPUT_DIR}/${ID}/wgs/${ID}.recal_hard_filtered_snvs.vcf.gz.tbi

# Remove unnecessary asta headers
rm ${OUTPUT_DIR}/${ID}/wgs/hg38_fasta_headers.txt \
   ${OUTPUT_DIR}/${ID}/wgs/pers_fasta_headers.txt

######################################################################
# Clean up VCF sample names and make R-friendly
######################################################################
for VCF in ${OUTPUT_DIR}/${ID}/wgs/*vcf.gz
do
    $BCFTOOLS query -l $VCF | tr " -/" "___" | awk '{ print "Individual_"$0 }' | bcftools reheader --samples - $VCF --output ${VCF}.new
    mv ${VCF}.new $VCF
    $BCFTOOLS index --tbi --force $VCF
done

######################################################################
# Quick stats on control & treatment
######################################################################
# Control
# Basic stats
$SAMTOOLS flagstat \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
          > ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.flagstat.results.txt

$SAMTOOLS stats \
          ${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
          > ${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.bamstats.results.txt

# Get metrics on alignment
$JAVA -jar $PICARD CollectAlignmentSummaryMetrics \
                   R=$REF \
                   I=${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
                   O=${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.alignment_metrics.txt

$JAVA -jar $PICARD CollectInsertSizeMetrics \
                   INPUT=${OUTPUT_DIR}/${ID}/wgs/${ID}.wgs.dedup.recal.bam \
                   OUTPUT=${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.insert_metrics.txt \
                   HISTOGRAM_FILE=${OUTPUT_DIR}/${ID}/stats/${ID}.wgs.insert_size_histogram.pdf

######################################################################
# End message
######################################################################

echo "Done with $(basename ${OUTPUT_DIR})/${ID} Module 1a: Setup, Calling, & Phasing"
