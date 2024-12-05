#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

################################################################################
# Load packages
# CRAN
library(data.table)

# Bioconductor
library(plyranges)
library(VariantAnnotation)
library(MutationalPatterns)
library(genomation)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(MotifDb)
library(systemPipeR)
library(biomaRt)

# CRAN
library(LDlinkR)

# Bioconductor
library(AllelicImbalance)

# CRAN
library(broom)
library(gt)
library(gtsummary)
library(UpSetR)
library(pheatmap)
library(umap)
library(ggdendro)
library(dendextend)
library(pvclust)

# Bioconductor
library(motifbreakR)
library(DiffBind)
library(rGREAT)

# CRAN
library(gprofiler2)

# System
library(parallel)

# GitHub
library(marge)

# Bioconductor
library(monaLisa)
library(universalmotif)

# CRAN
library(tidyverse)
library(dtplyr)

################################################################################
# Set seed
set.seed(0)

################################################################################
# Constants
const_canonical_chromosomes <- stringr::str_c("chr", c(1:22, "X"))
const_canonical_autosomes <- stringr::str_c("chr", 1:22)

################################################################################
# Create metadata dataframe
cohort_metadata <-
  data.frame(sample = 
               c("Individual_1_BC", "Individual_2_BC", "Individual_3_BC", "Individual_4_BC", "Individual_5_BC", "Individual_6_BC", "Individual_7_BC", "Individual_8_BC",
                 "Individual_1_LP", "Individual_2_LP", "Individual_3_LP", "Individual_4_LP", "Individual_5_LP", "Individual_6_LP", "Individual_7_LP", "Individual_8_LP",
                 "Individual_1_LC", "Individual_2_LC", "Individual_3_LC", "Individual_4_LC", "Individual_5_LC", "Individual_6_LC", "Individual_7_LC", "Individual_8_LC",
                 "Individual_1_SC", "Individual_2_SC", "Individual_3_SC", "Individual_4_SC", "Individual_5_SC", "Individual_6_SC", "Individual_7_SC", "Individual_8_SC",
                 "MCF10A", "hTERT_L9", "hTERT_L2"),
             cell_type =
               c(rep("BC", 8), rep("LP", 8), rep("LC", 8), rep("SC", 8), rep("Cell line", 3)),
             individual = 
               c(rep(c("14-17", "14-18", "11-18", "22-18", "24-18", "38-18", "30-18", "15-18"), 4), "MCF10A", "hTERT_L9", "hTERT_L2"))

# Create metadata dataframe
cohort_metadata_v2 <-
  data.frame(sample = 
               c("Individual_1_BC", "Individual_2_BC", "Individual_3_BC", "Individual_4_BC", "Individual_5_BC", "Individual_6_BC", "Individual_7_BC", "Individual_8_BC",
                 "Individual_1_LP", "Individual_2_LP", "Individual_3_LP", "Individual_4_LP", "Individual_5_LP", "Individual_6_LP", "Individual_7_LP", "Individual_8_LP",
                 "Individual_1_LC", "Individual_2_LC", "Individual_3_LC", "Individual_4_LC", "Individual_5_LC", "Individual_6_LC", "Individual_7_LC", "Individual_8_LC",
                 "Individual_1_SC", "Individual_2_SC", "Individual_3_SC", "Individual_4_SC", "Individual_5_SC", "Individual_6_SC", "Individual_7_SC", "Individual_8_SC",
                 "MCF10A", "hTERT_L9", "hTERT_L2"),
             Individual = 
               c("Individual_1", "Individual_2", "Individual_3", "Individual_4", "Individual_5", "Individual_6", "Individual_7", "Individual_8",
                 "Individual_1", "Individual_2", "Individual_3", "Individual_4", "Individual_5", "Individual_6", "Individual_7", "Individual_8",
                 "Individual_1", "Individual_2", "Individual_3", "Individual_4", "Individual_5", "Individual_6", "Individual_7", "Individual_8",
                 "Individual_1", "Individual_2", "Individual_3", "Individual_4", "Individual_5", "Individual_6", "Individual_7", "Individual_8",
                 "MCF10A", "hTERT_L9", "hTERT_L2"),
             cell_type =
               c(rep("BC", 8), rep("LP", 8), rep("LC", 8), rep("SC", 8), rep("Cell line", 3)),
             individual = 
               c(rep(c("14-17", "14-18", "11-18", "22-18", "24-18", "38-18", "30-18", "15-18"), 4), "MCF10A", "hTERT_L9", "hTERT_L2"),
             cemt = 
               c(
                 # BC
                 "CEMT_154",
                 "CEMT_173",
                 "CEMT_177",
                 "CEMT_181",
                 "CEMT_185",
                 "CEMT_190",
                 "CEMT_194",
                 "CEMT_198",
                 
                 # LP
                 "CEMT_155",
                 "CEMT_174",
                 "CEMT_178",
                 "CEMT_182",
                 "CEMT_186",
                 "CEMT_191",
                 "CEMT_195",
                 "CEMT_199",
                 
                 # LC
                 "CEMT_156",
                 "CEMT_175",
                 "CEMT_179",
                 "CEMT_183",
                 "CEMT_187",
                 "CEMT_192",
                 "CEMT_196",
                 "CEMT_200",
                 
                 # SC
                 "CEMT_157",
                 "CEMT_176",
                 "CEMT_180",
                 "CEMT_184",
                 "CEMT_188",
                 "CEMT_193",
                 "CEMT_197",
                 "CEMT_201",
                 
                 # Cell lines
                 "CEMT_7",
                 "CEMT_8",
                 "CEMT_9"))

# Add Misha folders and keep primary samples only
cohort_metadata_v3 <-
  cohort_metadata_v2 %>% 
  mutate(misha_folder_annotation = str_c("pre_", cell_type, "_", str_replace(individual, "-", "."))) %>% 
  filter(cell_type != "Cell line")

################################################################################
# Histone marks (ordered)
vector_marks <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")

# Cell types (ordered)
vector_celltypes <- c("LC", "LP", "BC", "SC")

# Individual names from BAM files
vector_individuals <- c("Individual_14_17_LP_DNA",
                        "Individual_14_18_BC__356K_RNA",
                        "Individual_11_18_BC_RNA_DNA",
                        "Individual_22_18_BC_RNA_DNA",
                        "Individual_24_18_BC_RNA_DNA",
                        "Individual_38_18_LC_RNA_DNA",
                        "Individual_30_18_LC_RNA_DNA",
                        "Individual_15_18_LC_RNA_DNA")

# Numbered individuals
vector_individuals_encoded <- c("Individual_1",
                                "Individual_2",
                                "Individual_3",
                                "Individual_4",
                                "Individual_5",
                                "Individual_6",
                                "Individual_7",
                                "Individual_8")

vector_individuals_encoded_space <- vector_individuals_encoded %>% str_replace("_", " ")
vector_individuals_encoded_space_upset_format <- str_c("  ", vector_individuals_encoded_space)

vector_individuals_encoded_numbers <- 
  c("1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8")

# New constants
# Define modalities
vector_modalities <- 
  c("H3K27ac_LC",  "H3K27ac_LP",  "H3K27ac_BC",  "H3K27ac_SC",  
    "H3K4me1_LC",  "H3K4me1_LP",  "H3K4me1_BC",  "H3K4me1_SC",
    "H3K4me3_LC",  "H3K4me3_LP",  "H3K4me3_BC",  "H3K4me3_SC",
    "H3K36me3_LC", "H3K36me3_LP", "H3K36me3_BC", "H3K36me3_SC",
    "H3K27me3_LC", "H3K27me3_LP", "H3K27me3_BC", "H3K27me3_SC",
    "H3K9me3_LC",  "H3K9me3_LP",  "H3K9me3_BC",  "H3K9me3_SC")

################################################################################
# Set variability categories
ba_variability_categories <- 
  c("blacklisted",
    "unoccupied", 
    "highly variable", 
    "moderately variable", 
    "slighly variable",
    "consistent")

ba_variability_categories_nonabsent <- 
  c("highly variable", 
    "moderately variable", 
    "slighly variable",
    "consistent")

################################################################################
# Important features
vector_important_features <- c("promoter", "enhancer_pellacani", "exon", "intron", "transcript", "intergenic")
vector_important_features_v2 <- c("promoter", "enhancer_pellacani", "transcript")
vector_important_features_simplified <- c("promoter", "enhancer", "exon", "intron", "transcript", "intergenic")

################################################################################
# Colors
vector_celltype_colors <- 
  c("LC" = "#EB8353",
    "LP" = "#FF3739",
    "BC" = "#3A68AE",
    "SC" = "#464546")

vector_celltype_colors_light <-
  c("LC" = "#f0a27e",
    "LP" = "#ff696a",
    "BC" = "#638ccb",
    "SC" = "#b5b4b5")

vector_mark_colors <- 
  c("H3K27ac" = "#3366FF",
    "H3K4me1" = "#FF6633",
    "H3K4me3" = "#FF0000",
    "H3K36me3" = "#990099",
    "H3K27me3" = "#663300",
    "H3K9me3" = "#000066")

# Get primary CEMTs
vector_primary_cemts <- cohort_metadata_v2 %>% filter(cell_type %in% vector_celltypes) %>% pull(cemt) %>% sort()

# Get mark sublists
vector_marks_fig2         <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3")
vector_marks_fig2_v2      <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
vector_marks_remaining    <- c("H3K36me3", "H3K9me3")
vector_marks_remaining_v2 <- c("H3K36me3")
vector_marks_active       <- c("H3K27ac", "H3K4me3")

# Set vectors for features
vector_features_fig_panel_old <- c("promoter", "enhancer", "transcript")
vector_features_fig_panel     <- c("promoters", "enhancers", "transcripts")
vector_features_fig_panel_v2  <- c("promoter", "enhancer_pellacani", "transcript")

