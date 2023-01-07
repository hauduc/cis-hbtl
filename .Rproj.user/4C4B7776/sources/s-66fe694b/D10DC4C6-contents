#!/gsc/software/linux-x86_64-centos7/R-4.1.0/lib64/R/bin/R
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

#################################################
#################################################
##                                             ##
##          SPAA RPKM LINKER v1.0              ##
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
# GSC RPKM files
arg_RPKMS <- arg[4] %>% read.table() %>% pull()
# GTF
arg_GTF <- "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/genome_data/genomic_tracks/genes/Homo_sapiens.GRCh38.102.chr_names.sorted.bed"

# Load in correct consensus.spaa.bed
load(str_c(arg_OUTPUT_DIR, "/", "spaas", "/", "raw", "/", arg_MARK, "_", arg_CELLTYPE, "_spaa.bed.RData"))
###############################################
#                                             #
#           END DATA LOADING ZONE             #
#                                             #
###############################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Intersect SPAAs with closest downstream gene ENSG name in preparation for joining with RPKM files on that column
# Fix colnames and column types
# Remove redundancies
consensus.spaa.rpkm.bed <- bedr(input = list(a = consensus.spaa.bed, 
                                             b = arg_GTF),
     engine = "bedtools",
     method = "closest",
     params = "-D ref -iu -t first", 
     tmpDir = "/tmp") %>% 
  select(1:21, 25) %>% 
  `colnames<-`(c(colnames(consensus.spaa.bed), "ensg_nearest_downstream")) %>% 
  map_at(8:17, as.logical) %>% # fix the individual presence columns, plus spaa and spaa_negative
  map_at(18:21, as.numeric) %>% 
  as.data.frame() %>% 
  distinct(across(c("region_identifier", "ensg_nearest_downstream")), .keep_all = TRUE)

# Left outer join of ENSG RPKM into correct individual column
add_rpkm <- function(rpkm_file_index) {
  read.table(arg_RPKMS[rpkm_file_index], sep = "\t") %>% 
    select(1, 3) %>%
    `colnames<-`(c("RPKM_file_ensgs", "RPKM")) %>% 
    left_join(consensus.spaa.rpkm.bed, ., by = c("ensg_nearest_downstream" = "RPKM_file_ensgs")) %>% 
    select(last_col()) %>% 
    pull()
}

# Running the above function, fill in the correct, corresponding RPKM values for the current modality (e.g. H3K4me3_BC)
consensus.spaa.rpkm.bed <- consensus.spaa.rpkm.bed %>% 
  mutate(Individual_14_17_LP_DNA_RPKM       = add_rpkm(1),
         Individual_14_18_BC__356K_RNA_RPKM = add_rpkm(2),
         Individual_11_18_BC_RNA_DNA_RPKM   = add_rpkm(3),
         Individual_22_18_BC_RNA_DNA_RPKM   = add_rpkm(4),
         Individual_24_18_BC_RNA_DNA_RPKM   = add_rpkm(5),
         Individual_38_18_LC_RNA_DNA_RPKM   = add_rpkm(6),
         Individual_30_18_LC_RNA_DNA_RPKM   = add_rpkm(7),
         Individual_15_18_LC_RNA_DNA_RPKM   = add_rpkm(8))

# In preparation for doing statistical tests along SPAAs, create a proper long dataframe that allows you to do tests across rows
consensus.spaa.rpkm.bed.long <- consensus.spaa.rpkm.bed %>% 
  pivot_longer(cols = c(Individual_14_17_LP_DNA,
                        Individual_14_18_BC__356K_RNA,
                        Individual_11_18_BC_RNA_DNA,
                        Individual_22_18_BC_RNA_DNA,
                        Individual_24_18_BC_RNA_DNA,
                        Individual_38_18_LC_RNA_DNA,
                        Individual_30_18_LC_RNA_DNA,
                        Individual_15_18_LC_RNA_DNA), 
               names_to = "Individual", 
               values_to = "is_present_in_spaa") %>% 
  pivot_longer(cols = c(Individual_14_17_LP_DNA_RPKM,
                        Individual_14_18_BC__356K_RNA_RPKM,
                        Individual_11_18_BC_RNA_DNA_RPKM,
                        Individual_22_18_BC_RNA_DNA_RPKM,
                        Individual_24_18_BC_RNA_DNA_RPKM,
                        Individual_38_18_LC_RNA_DNA_RPKM,
                        Individual_30_18_LC_RNA_DNA_RPKM,
                        Individual_15_18_LC_RNA_DNA_RPKM), 
               names_to = "which_individuals_expression_level", 
               values_to = "expression_level") %>% 
  filter(which_individuals_expression_level == str_c(Individual, "_RPKM"))

# Do Wilcoxon rank-order test across long dataframe by testing the TRUE/FALSE groups within each SPAA
wilcoxon_pvalues <- list()
for (pos in (consensus.spaa.rpkm.bed.long %>% drop_na(expression_level) %>% pull(region_identifier) %>% unique())) {
  
  temp_df <- data.frame(region_identifier = pos,
                        wilcox_test_p_val = wilcox.test((consensus.spaa.rpkm.bed.long %>% filter(region_identifier == pos & is_present_in_spaa == TRUE)  %>% pull(expression_level)),
                                                        (consensus.spaa.rpkm.bed.long %>% filter(region_identifier == pos & is_present_in_spaa == FALSE) %>% pull(expression_level)),
                                                        exact = FALSE)$p.value)
  wilcoxon_pvalues[[pos]] <- temp_df
  rm(temp_df)
  }
# row bind all into a dataframe containing region identifier and p-value for later joining
all_wilcoxon_pvalues <- do.call(rbind, wilcoxon_pvalues)
rm(wilcoxon_pvalues)

error.resistant.t.test.p.val <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

# Do t-test test across long dataframe by testing the TRUE/FALSE groups within each SPAA
t_test_pvalues <- list()
for (pos in (consensus.spaa.rpkm.bed.long %>% drop_na(expression_level) %>% pull(region_identifier) %>% unique())) {
  
  temp_df <- data.frame(region_identifier = pos,
                        t_test_p_val = error.resistant.t.test.p.val((consensus.spaa.rpkm.bed.long %>% filter(region_identifier == pos & is_present_in_spaa == TRUE) %>% pull(expression_level)),
                                                                    (consensus.spaa.rpkm.bed.long %>% filter(region_identifier == pos & is_present_in_spaa == FALSE) %>% pull(expression_level)), 
                                                                    var.equal = TRUE))
  
  t_test_pvalues[[pos]] <- temp_df
  rm(temp_df)
}
# row bind all into a dataframe containing region identifier and p-value for later joining
all_t_test_pvalues <- do.call(rbind, t_test_pvalues)
rm(t_test_pvalues)

# Join p-values with the original (non-long) consensus.spaa.rpkm.bed, and create new columns with adjusted p-values
consensus.spaa.rpkm.bed <- consensus.spaa.rpkm.bed %>% 
  left_join(all_wilcoxon_pvalues, by = "region_identifier") %>% 
  left_join(all_t_test_pvalues, by = "region_identifier") %>% 
  mutate(wilcox_test_p_val_fdr_adj = wilcox_test_p_val %>% p.adjust(method = "fdr"),
         t_test_p_val_fdr_adj      = t_test_p_val      %>% p.adjust(method = "fdr"))
  
# Save output with proper name
write.table(consensus.spaa.rpkm.bed, 
            file = str_c(arg_OUTPUT_DIR, "/", "spaas", "/", "annotated_rpkm", "/", arg_MARK, "_", arg_CELLTYPE, "_spaa.rpkm.bed"),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)

save(consensus.spaa.rpkm.bed, 
     file = str_c(arg_OUTPUT_DIR, "/", "spaas", "/", "annotated_rpkm", "/", arg_MARK, "_", arg_CELLTYPE, "_spaa.rpkm.bed.RData"),
     compress = "bzip2")

print("Module 2 Step 3: RPKM Linker - Complete")
print(arg_MARK)
print(arg_CELLTYPE)
