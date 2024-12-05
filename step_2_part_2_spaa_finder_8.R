#!/gsc/software/linux-x86_64-centos7/R-4.1.0/lib64/R/bin/R
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

#################################################
#################################################
##                                             ##
##            SPAA FINDER v1.0                 ##
##             8 INDIVIDUALS                   ##
##                                             ##
#################################################
#################################################

library(tidyverse)
# Setting the right environmental variables to bedr required binaries
PATH_default <- Sys.getenv("PATH")
PATH_to_bedr_tools <- str_c("/gsc/software/linux-x86_64-centos7/bedtools-2.27.1/bin",
                            "/gsc/software/linux-x86_64-centos7/bedops-2.4.35/bin",
                            "/home/ahauduc/anaconda3/bin",
                            sep = ":")
Sys.setenv("PATH" = str_c(PATH_default, 
                          PATH_to_bedr_tools, 
                          sep = ":"))
library(bedr)

# Read arguments
arg <- commandArgs(trailingOnly = TRUE)
# CONSTANTS
const_canonical_chromosomes <- str_c("chr", c(1:22, "X", "Y"))
###############################################
#                                             #
#            DATA LOADING ZONE                #
#                                             #
###############################################
# Directory where output files will be written
arg_OUTPUT_DIR <- arg[1]
# Data source
arg_MARK <- arg[2]
arg_CELLTYPE <- arg[3]
# MACS2 format BED files containing histone peak calls
arg_PEAKS <- arg[4] %>% read.table() %>% pull()

# Handle correct colnames depending on if narrowPeak or broadPeak file
if (arg_MARK == "H3K4me3" | arg_MARK == "H3K27ac" | arg_MARK == "H3K4me1") {
  const_peakfile_colnames <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
} else if (arg_MARK == "H3K9me3" | arg_MARK == "H3K27me3" | arg_MARK == "H3K36me3") {
  const_peakfile_colnames <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
}

# Load in consensus.bed and arg_VCFS_LOADED from setup
load(str_c(arg_OUTPUT_DIR, "setup", "spaa_finder_8_setup.RData", sep = "/"))
###############################################
#                                             #
#           END DATA LOADING ZONE             #
#                                             #
###############################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# PART 2
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
###########################################################################################################################################################################
# Creating lists of BED files loaded into R
###########################################################################################################################################################################
arg_PEAKS_LOADED <- list()
for(FILE in arg_PEAKS) {
  arg_PEAKS_LOADED[[FILE]] <- read.table(FILE, sep = "\t", col.names = const_peakfile_colnames) %>% 
    filter(chr %in% const_canonical_chromosomes) %>% 
    convert2bed() %>% 
    bedr.sort.region(method = "lexicographical", engine = "bedtools")
}

#######################################################
##!!! Fix error with overlapping column name between VCFS_LOADED and consensus.bed that causes error during mean quality calculation
consensus.bed <- consensus.bed %>%
  rename(rs = ID)
##!!!!

#######################################################
# Creating presence matrices for SNV's in each peak BED file within the BED lists
# Create "Peak Present" columns for the individuals
consensus.bed <- consensus.bed %>% 
  mutate(Individual_14_17_LP_DNA_IN_PEAK       = consensus.bed %in.region% arg_PEAKS_LOADED[[1]],
         Individual_14_18_BC__356K_RNA_IN_PEAK = consensus.bed %in.region% arg_PEAKS_LOADED[[2]],
         Individual_11_18_BC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[3]],
         Individual_22_18_BC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[4]],
         Individual_24_18_BC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[5]],
         Individual_38_18_LC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[6]],
         Individual_30_18_LC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[7]],
         Individual_15_18_LC_RNA_DNA_IN_PEAK   = consensus.bed %in.region% arg_PEAKS_LOADED[[8]])

# Testing for presence of SPAAs
# Each SPAA must be a variant where each individual has an SNV and a peak present OR no SNV and no peak present
# Also, SPAAs cannot be areas where all or no individuals have SNV and peak
consensus.bed <- consensus.bed %>% 
  mutate(spaa = Individual_14_17_LP_DNA       +  Individual_14_17_LP_DNA_IN_PEAK       != 1 &
                Individual_14_18_BC__356K_RNA +  Individual_14_18_BC__356K_RNA_IN_PEAK != 1 &
                Individual_11_18_BC_RNA_DNA   +  Individual_11_18_BC_RNA_DNA_IN_PEAK   != 1 &
                Individual_22_18_BC_RNA_DNA   +  Individual_22_18_BC_RNA_DNA_IN_PEAK   != 1 &
                Individual_24_18_BC_RNA_DNA   +  Individual_24_18_BC_RNA_DNA_IN_PEAK   != 1 &
                Individual_38_18_LC_RNA_DNA   +  Individual_38_18_LC_RNA_DNA_IN_PEAK   != 1 &
                Individual_30_18_LC_RNA_DNA   +  Individual_30_18_LC_RNA_DNA_IN_PEAK   != 1 &
                Individual_15_18_LC_RNA_DNA   +  Individual_15_18_LC_RNA_DNA_IN_PEAK   != 1 &
                                            
                Individual_14_17_LP_DNA       +  Individual_14_17_LP_DNA_IN_PEAK       +   
                Individual_14_18_BC__356K_RNA +  Individual_14_18_BC__356K_RNA_IN_PEAK + 
                Individual_11_18_BC_RNA_DNA   +  Individual_11_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_22_18_BC_RNA_DNA   +  Individual_22_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_24_18_BC_RNA_DNA   +  Individual_24_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_38_18_LC_RNA_DNA   +  Individual_38_18_LC_RNA_DNA_IN_PEAK   + 
                Individual_30_18_LC_RNA_DNA   +  Individual_30_18_LC_RNA_DNA_IN_PEAK   + 
                Individual_15_18_LC_RNA_DNA   +  Individual_15_18_LC_RNA_DNA_IN_PEAK   != 0 &
                                            
                Individual_14_17_LP_DNA       +  Individual_14_17_LP_DNA_IN_PEAK       +   
                Individual_14_18_BC__356K_RNA +  Individual_14_18_BC__356K_RNA_IN_PEAK + 
                Individual_11_18_BC_RNA_DNA   +  Individual_11_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_22_18_BC_RNA_DNA   +  Individual_22_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_24_18_BC_RNA_DNA   +  Individual_24_18_BC_RNA_DNA_IN_PEAK   + 
                Individual_38_18_LC_RNA_DNA   +  Individual_38_18_LC_RNA_DNA_IN_PEAK   + 
                Individual_30_18_LC_RNA_DNA   +  Individual_30_18_LC_RNA_DNA_IN_PEAK   + 
                Individual_15_18_LC_RNA_DNA   +  Individual_15_18_LC_RNA_DNA_IN_PEAK   != 16)

# Search for NEGATIVE SPAAs, where lack of variant is associated with a peak, and ignore all-variant or all-peak driven negative SPAAs, which are not SPAAS
consensus.bed <- consensus.bed %>% 
  mutate(spaa_negative = Individual_14_17_LP_DNA  +  Individual_14_17_LP_DNA_IN_PEAK       == 1 &
           Individual_14_18_BC__356K_RNA +  Individual_14_18_BC__356K_RNA_IN_PEAK == 1 &
           Individual_11_18_BC_RNA_DNA   +  Individual_11_18_BC_RNA_DNA_IN_PEAK   == 1 &
           Individual_22_18_BC_RNA_DNA   +  Individual_22_18_BC_RNA_DNA_IN_PEAK   == 1 &
           Individual_24_18_BC_RNA_DNA   +  Individual_24_18_BC_RNA_DNA_IN_PEAK   == 1 &
           Individual_38_18_LC_RNA_DNA   +  Individual_38_18_LC_RNA_DNA_IN_PEAK   == 1 &
           Individual_30_18_LC_RNA_DNA   +  Individual_30_18_LC_RNA_DNA_IN_PEAK   == 1 &
           Individual_15_18_LC_RNA_DNA   +  Individual_15_18_LC_RNA_DNA_IN_PEAK   == 1 &
           
           Individual_14_17_LP_DNA       +
           Individual_14_18_BC__356K_RNA +
           Individual_11_18_BC_RNA_DNA   +  
           Individual_22_18_BC_RNA_DNA   +  
           Individual_24_18_BC_RNA_DNA   + 
           Individual_38_18_LC_RNA_DNA   +  
           Individual_30_18_LC_RNA_DNA   +  
           Individual_15_18_LC_RNA_DNA   != 8 &
           
           Individual_14_17_LP_DNA_IN_PEAK       +   
           Individual_14_18_BC__356K_RNA_IN_PEAK + 
           Individual_11_18_BC_RNA_DNA_IN_PEAK   + 
           Individual_22_18_BC_RNA_DNA_IN_PEAK   + 
           Individual_24_18_BC_RNA_DNA_IN_PEAK   + 
           Individual_38_18_LC_RNA_DNA_IN_PEAK   + 
           Individual_30_18_LC_RNA_DNA_IN_PEAK   + 
           Individual_15_18_LC_RNA_DNA_IN_PEAK   != 8)

###########################################################################################################################################################################
# Subsetting SPAAs
###########################################################################################################################################################################
# Pulling out the number of individuals involved in each SPAA
# First, create more compact spaa-only dataframe with unique identifiers
consensus.spaa.bed <- consensus.bed %>% 
  filter(spaa == TRUE | spaa_negative == TRUE) %>% 
  select(-c(Individual_14_17_LP_DNA_IN_PEAK,
            Individual_14_18_BC__356K_RNA_IN_PEAK, 
            Individual_11_18_BC_RNA_DNA_IN_PEAK, 
            Individual_22_18_BC_RNA_DNA_IN_PEAK, 
            Individual_24_18_BC_RNA_DNA_IN_PEAK, 
            Individual_38_18_LC_RNA_DNA_IN_PEAK, 
            Individual_30_18_LC_RNA_DNA_IN_PEAK, 
            Individual_15_18_LC_RNA_DNA_IN_PEAK))

# Creating a new columns containing the count of number of individuals involved in the SPAA 
consensus.spaa.bed <- consensus.spaa.bed %>% 
  rowwise() %>% 
  mutate(individuals_involved = sum(c(Individual_14_17_LP_DNA,
                                      Individual_14_18_BC__356K_RNA, 
                                      Individual_11_18_BC_RNA_DNA, 
                                      Individual_22_18_BC_RNA_DNA, 
                                      Individual_24_18_BC_RNA_DNA, 
                                      Individual_38_18_LC_RNA_DNA, 
                                      Individual_30_18_LC_RNA_DNA, 
                                      Individual_15_18_LC_RNA_DNA))) %>% 
  as.data.frame()

# Creating a 4th, 5th, and 6th columns delineating unique region identifier for each BED record, for later HOMER compatibility
# Rename ID to prevent conflicts later with the annotation steps
consensus.spaa.bed <- consensus.spaa.bed %>% 
  mutate(region_identifier = str_c(arg_MARK, arg_CELLTYPE, "spaa", rownames(consensus.spaa.bed), sep = "_"), 
         mark = arg_MARK,
         celltype = arg_CELLTYPE,
         .after = "end")

###########################################################################################################################################################################
# Calculating mean quality for each SPAA 
###########################################################################################################################################################################
consensus.spaa.bed <- consensus.spaa.bed %>% 
  mutate(Individual_14_17_LP_DNA_QUAL       = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[1]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_14_18_BC__356K_RNA_QUAL = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[2]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_11_18_BC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[3]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_22_18_BC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[4]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_24_18_BC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[5]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_38_18_LC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[6]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_30_18_LC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[7]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric(),
         Individual_15_18_LC_RNA_DNA_QUAL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[8]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(QUAL) %>% na_if(".") %>% as.numeric()) %>% 
  rowwise() %>% 
  mutate(mean_variant_call_quality = mean(c(Individual_14_17_LP_DNA_QUAL,
                                            Individual_14_18_BC__356K_RNA_QUAL, 
                                            Individual_11_18_BC_RNA_DNA_QUAL, 
                                            Individual_22_18_BC_RNA_DNA_QUAL, 
                                            Individual_24_18_BC_RNA_DNA_QUAL, 
                                            Individual_38_18_LC_RNA_DNA_QUAL, 
                                            Individual_30_18_LC_RNA_DNA_QUAL, 
                                            Individual_15_18_LC_RNA_DNA_QUAL), na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  select(-c(Individual_14_17_LP_DNA_QUAL,
            Individual_14_18_BC__356K_RNA_QUAL, 
            Individual_11_18_BC_RNA_DNA_QUAL, 
            Individual_22_18_BC_RNA_DNA_QUAL, 
            Individual_24_18_BC_RNA_DNA_QUAL, 
            Individual_38_18_LC_RNA_DNA_QUAL, 
            Individual_30_18_LC_RNA_DNA_QUAL, 
            Individual_15_18_LC_RNA_DNA_QUAL))

###########################################################################################################################################################################
# Calculating mean allele frequency for each SPAA 
###########################################################################################################################################################################
consensus.spaa.bed <- consensus.spaa.bed %>% 
  mutate(Individual_14_17_LP_DNA_ACAF       = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[1]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_14_18_BC__356K_RNA_ACAF = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[2]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_11_18_BC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[3]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_22_18_BC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[4]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_24_18_BC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[5]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_38_18_LC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[6]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_30_18_LC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[7]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5),
         Individual_15_18_LC_RNA_DNA_ACAF   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_VCFS_LOADED[[8]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(INFO) %>% na_if(".") %>% str_detect("^AC=2;AF=1") %>% ifelse(1, 0.5)) %>% 
  rowwise() %>% 
  mutate(mean_allele_frequency = mean(c(Individual_14_17_LP_DNA_ACAF,
                                        Individual_14_18_BC__356K_RNA_ACAF, 
                                        Individual_11_18_BC_RNA_DNA_ACAF, 
                                        Individual_22_18_BC_RNA_DNA_ACAF, 
                                        Individual_24_18_BC_RNA_DNA_ACAF, 
                                        Individual_38_18_LC_RNA_DNA_ACAF, 
                                        Individual_30_18_LC_RNA_DNA_ACAF, 
                                        Individual_15_18_LC_RNA_DNA_ACAF), na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  select(-c(Individual_14_17_LP_DNA_ACAF,
            Individual_14_18_BC__356K_RNA_ACAF, 
            Individual_11_18_BC_RNA_DNA_ACAF, 
            Individual_22_18_BC_RNA_DNA_ACAF, 
            Individual_24_18_BC_RNA_DNA_ACAF, 
            Individual_38_18_LC_RNA_DNA_ACAF, 
            Individual_30_18_LC_RNA_DNA_ACAF, 
            Individual_15_18_LC_RNA_DNA_ACAF))

###########################################################################################################################################################################
# Calculating mean peak q value for each SPAA
###########################################################################################################################################################################
consensus.spaa.bed <- consensus.spaa.bed %>% 
  mutate(Individual_14_17_LP_DNA_PQVL       = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[1]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_14_18_BC__356K_RNA_PQVL = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[2]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_11_18_BC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[3]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_22_18_BC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[4]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_24_18_BC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[5]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_38_18_LC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[6]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_30_18_LC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[7]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric(),
         Individual_15_18_LC_RNA_DNA_PQVL   = bedr(input = list(a = consensus.spaa.bed, 
                                                                b = arg_PEAKS_LOADED[[8]]),
                                                   engine = "bedtools",
                                                   method = "intersect",
                                                   params = "-wao", 
                                                   tmpDir = "/tmp") %>% pull(qValue) %>% na_if(".") %>% as.numeric()) %>% 
  rowwise() %>% 
  mutate(mean_peak_q_value = mean(c(Individual_14_17_LP_DNA_PQVL,
                                    Individual_14_18_BC__356K_RNA_PQVL, 
                                    Individual_11_18_BC_RNA_DNA_PQVL, 
                                    Individual_22_18_BC_RNA_DNA_PQVL, 
                                    Individual_24_18_BC_RNA_DNA_PQVL, 
                                    Individual_38_18_LC_RNA_DNA_PQVL, 
                                    Individual_30_18_LC_RNA_DNA_PQVL, 
                                    Individual_15_18_LC_RNA_DNA_PQVL), na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  select(-c(Individual_14_17_LP_DNA_PQVL,
            Individual_14_18_BC__356K_RNA_PQVL, 
            Individual_11_18_BC_RNA_DNA_PQVL, 
            Individual_22_18_BC_RNA_DNA_PQVL, 
            Individual_24_18_BC_RNA_DNA_PQVL, 
            Individual_38_18_LC_RNA_DNA_PQVL, 
            Individual_30_18_LC_RNA_DNA_PQVL, 
            Individual_15_18_LC_RNA_DNA_PQVL))

####################################################################################################################################
# Write SPAA BED to output folder
####################################################################################################################################
write.table(consensus.spaa.bed, 
            file = str_c(arg_OUTPUT_DIR, "/", "spaas", "/", "raw", "/", arg_MARK, "_", arg_CELLTYPE, "_spaa.bed"),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)

save(consensus.spaa.bed, 
     file = str_c(arg_OUTPUT_DIR, "/", "spaas", "/", "raw", "/", arg_MARK, "_", arg_CELLTYPE, "_spaa.bed.RData"),
     compress = "bzip2")

print("Module 2 Step 2: SPAA Finder - Complete")
print(arg_MARK)
print(arg_CELLTYPE)
