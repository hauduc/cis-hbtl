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

vector_marks <- c("H3K4me3", "H3K4me1", "H3K27ac", "H3K9me3", "H3K27me3", "H3K36me3")
vector_celltypes <- c("BC", "LP", "LC", "SC")
vector_individuals <- c("Individual_14_17_LP_DNA",
                        "Individual_14_18_BC__356K_RNA",
                        "Individual_11_18_BC_RNA_DNA",
                        "Individual_22_18_BC_RNA_DNA",
                        "Individual_24_18_BC_RNA_DNA",
                        "Individual_38_18_LC_RNA_DNA",
                        "Individual_30_18_LC_RNA_DNA",
                        "Individual_15_18_LC_RNA_DNA")

# Initialize list
list_treat_pileup_bw_wasp_files <- list()

for (FILE in list.files("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/filelists/bw", full.names = TRUE)) {
  # Pull out name of the modality
  NAME <- FILE %>% str_remove("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/filelists/bw/breast_normal_") %>% str_remove(".tsv")
  
  # Create vector containing modality files as a correctly-named component of the list
  list_treat_pileup_bw_wasp_files[[NAME]] <- FILE %>% read.table() %>% pull()
  
}

# Import files into R
# Initialize list
list_treat_pileup_bw_wasp <- list()

for (FILE in list.files("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/filelists/bw", full.names = TRUE)) {
  # Pull out name and contents
  NAME <- FILE %>% str_remove("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/filelists/bw/breast_normal_") %>% str_remove(".tsv")
  BW_VECTOR <- FILE %>% read.table() %>% pull()
  
  # Create modality sublist
  list_treat_pileup_bw_wasp[[NAME]] <- GRangesList()
  
  # Read in bigWigs
  for (i in seq_along(BW_VECTOR)) {
    list_treat_pileup_bw_wasp[[NAME]][[names(list_variant_calls)[i]]] <- read_bigwig(BW_VECTOR[i])
  }
}

# rm(list_treat_pileup_bw_wasp_files)

save(list_treat_pileup_bw_wasp, file = "list_treat_pileup_bw_wasp.RData", compress = "bzip2")
