#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu
library(plyranges)
library(VariantAnnotation)
library(MutationalPatterns)
library(genomation)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(systemPipeR)
library(biomaRt)
library(LDlinkR)
library(AllelicImbalance)
library(broom)
library(gt)
library(gtsummary)
library(UpSetR)
library(tidyverse)

const_canonical_chromosomes <- str_c("chr", c(1:22, "X"))

vector_marks <- c("H3K4me3", "H3K27ac", "H3K4me1","H3K9me3", "H3K27me3", "H3K36me3")
vector_celltypes <- c("BC", "LP", "LC", "SC")
vector_individuals <- c("Individual_14_17_LP_DNA",
                        "Individual_14_18_BC__356K_RNA",
                        "Individual_11_18_BC_RNA_DNA",
                        "Individual_22_18_BC_RNA_DNA",
                        "Individual_24_18_BC_RNA_DNA",
                        "Individual_38_18_LC_RNA_DNA",
                        "Individual_30_18_LC_RNA_DNA",
                        "Individual_15_18_LC_RNA_DNA")
##############################################################################################################################################################################################
##############################################################################################################################################################################################
# H3K36me3
rois_populated_tmp1_H3K36me3 <- 
  make_cohort_roi_table_fraction_overlap(mark = "H3K36me3", rois = roi_breast_normal_enhancers) %>% 
  make_n_top_varying_regions_from_cohort_roi_table_for_fraction_overlap(n_regions = 2000)

