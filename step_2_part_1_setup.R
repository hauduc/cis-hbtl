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
# VCFs
arg_VCFS <- arg[2] %>% read.table() %>% pull()

###############################################
#                                             #
#           END DATA LOADING ZONE             #
#                                             #
###############################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# STEP 1: SETUP
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Initialize and create BED list from VCF's
# Then filter out non-canonical chromosomes and sort individual bed files that will be in the list
# in order to optimize for later processing
arg_VCFS_LOADED <- list() # Initiate list
for(FILE in arg_VCFS) {
  arg_VCFS_LOADED[[FILE]] <- read.vcf(FILE) %>% 
    vcf2bed(other = c("ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>% # Include the other rows in your output BED
    filter(chr %in% const_canonical_chromosomes) %>%                               # get rid of the non-canonical chromosomes
    bedr.sort.region(method = "lexicographical", engine = "bedtools")              # Sort lexicographically
}

# Create consensus BED file by initializing and then looping through vcf bodies in the list, then binding each BED file successively to consensus.bed
consensus.bed <- data.frame() # initialize consensus.bed
for(FILE in seq_along(arg_VCFS)) {
  consensus.bed <- arg_VCFS_LOADED %>% 
    pluck(FILE) %>% 
    select(chr, start, end, ID) %>% 
    bind_rows(consensus.bed, .) # Turn your list of converted VCF's into one large bed file by binding rows successively
}

# Taking out locus duplicates in order to have a dataframe of unique SNV's for and non-canonical chromosomes, and sorting
# This is in order to pseudo "pivot_longer" and create a matrix where, from the list of SNV's contained by all 8 individuals, 
consensus.bed <- consensus.bed %>% 
  filter(chr %in% const_canonical_chromosomes) %>%         # filter out non-canonical chromosomes
  distinct(across(c("chr", "start", "end", "ID"))) %>%     # filter out redundant VCF records created by your merger
  bedr.sort.region(method = "lexicographical",             # Sort again just to make sure
                   engine = "bedtools")

# Create "SNV Present" columns for the individuals
# Need to modify so that it automatically created variables named after the sample name
consensus.bed <- consensus.bed %>% 
  mutate(Individual_14_17_LP_DNA       = consensus.bed %in.region% arg_VCFS_LOADED[[1]],
         Individual_14_18_BC__356K_RNA = consensus.bed %in.region% arg_VCFS_LOADED[[2]],
         Individual_11_18_BC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[3]],
         Individual_22_18_BC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[4]],
         Individual_24_18_BC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[5]],
         Individual_38_18_LC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[6]],
         Individual_30_18_LC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[7]],
         Individual_15_18_LC_RNA_DNA   = consensus.bed %in.region% arg_VCFS_LOADED[[8]])

# Export R data for spaa_finder_8.R
save(consensus.bed, arg_VCFS_LOADED, 
     file = str_c(arg_OUTPUT_DIR, "setup", "spaa_finder_8_setup.RData", sep = "/"),
     compress = "bzip2")

print("Module 2 Step 1: Setup - Complete")
