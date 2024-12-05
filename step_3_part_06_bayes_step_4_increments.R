#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

# Load constants
source("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/source/setup.R")

# Initialize ba-specific list
txdb_grlist_ba <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/txdb_grlist_ba.RDS")

# Load ba object
ba_whole_genome_anno_chr_all <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all.RDS")

# Chrom positions tibble
ba_chrom_positions <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_chrom_positions.RDS")

###########################################################################
###########################################################################
# Set var increments
var_increments1 <- seq(0, 0.500, by = 0.025)
var_increments2 <- seq(0.525, 1, by = 0.025)

# Get number of bins at all non-empty positions
# First half
ba_whole_genome_anno_chr_all_var_increments_all1 <-
  mclapply(vector_modalities, function(current_modality) {
    mclapply(var_increments1, function(current_var_increment) {
      current_mark <-     current_modality %>% str_split_i(pattern = "_", i = 1)
      current_celltype <- current_modality %>% str_split_i(pattern = "_", i = 2)
      current_bins_below_threshold <-
        ba_whole_genome_anno_chr_all %>%
        as_tibble() %>%
        filter(!is.na(get(str_c("bayes_mean_", current_modality))),
               get(str_c("bayes_mean_", current_modality)) != 0) %>%
        pull(str_c("bayes_mean_", current_modality)) %>%
        `<`(current_var_increment) %>%
        sum()
      
      current_valid_bins <-
        ba_whole_genome_anno_chr_all %>%
        as_tibble() %>%
        filter(!is.na(get(str_c("bayes_mean_", current_modality))),
               get(str_c("bayes_mean_", current_modality)) != 0) %>%
        pull(str_c("bayes_mean_", current_modality)) %>%
        length()
      
      tibble(mark = current_mark,
             celltype = current_celltype,
             feature_simplified = "all",
             threshold = current_var_increment,
             fraction_of_bins_below_threshold = current_bins_below_threshold/current_valid_bins) %>%
        mutate(mark = mark %>% factor(vector_marks),
               celltype = celltype %>% factor(vector_celltypes),
               feature_simplified = feature_simplified %>% factor(c("all", vector_important_features_simplified)))
      
    }, mc.cores = detectCores()/length(vector_modalities)) %>% do.call("bind_rows", .)
  }, mc.cores = length(vector_modalities)) %>% do.call("bind_rows", .)

# Pause
message("Done 1/4")
gc()
Sys.sleep(10)

# Get number of bins at all non-empty positions
# Second half
ba_whole_genome_anno_chr_all_var_increments_all2 <-
  mclapply(vector_modalities, function(current_modality) {
    mclapply(var_increments2, function(current_var_increment) {
      
      current_mark <-     current_modality %>% str_split_i(pattern = "_", i = 1)
      current_celltype <- current_modality %>% str_split_i(pattern = "_", i = 2)
      
      current_bins_below_threshold <-
        ba_whole_genome_anno_chr_all %>%
        as_tibble() %>%
        filter(!is.na(get(str_c("bayes_mean_", current_modality))),
               get(str_c("bayes_mean_", current_modality)) != 0) %>%
        pull(str_c("bayes_mean_", current_modality)) %>%
        `<`(current_var_increment) %>%
        sum()
      
      current_valid_bins <-
        ba_whole_genome_anno_chr_all %>%
        as_tibble() %>%
        filter(!is.na(get(str_c("bayes_mean_", current_modality))),
               get(str_c("bayes_mean_", current_modality)) != 0) %>%
        pull(str_c("bayes_mean_", current_modality)) %>%
        length()
      
      tibble(mark = current_mark,
             celltype = current_celltype,
             feature_simplified = "all",
             threshold = current_var_increment,
             fraction_of_bins_below_threshold = current_bins_below_threshold/current_valid_bins) %>%
        mutate(mark = mark %>% factor(vector_marks),
               celltype = celltype %>% factor(vector_celltypes),
               feature_simplified = feature_simplified %>% factor(c("all", vector_important_features_simplified)))
      
    }, mc.cores = detectCores()/length(vector_modalities)) %>% do.call("bind_rows", .)
  }, mc.cores = length(vector_modalities)) %>% do.call("bind_rows", .)

# Combine
ba_whole_genome_anno_chr_all_var_increments_all <- 
  bind_rows(ba_whole_genome_anno_chr_all_var_increments_all1, ba_whole_genome_anno_chr_all_var_increments_all2) %>% 
  arrange(mark, celltype, threshold)

rm(ba_whole_genome_anno_chr_all_var_increments_all1, ba_whole_genome_anno_chr_all_var_increments_all2)

# Pause
message("Done 2/4")
gc()
Sys.sleep(10)

# Add feature-specific posiitons
ba_whole_genome_anno_chr_all_var_increments_features1 <-
  mclapply(vector_modalities, function(current_modality) {
    mclapply(vector_important_features_v2, function(current_feature) {
      mclapply(var_increments1, function(current_var_increment) {
        
        current_mark <-     current_modality %>% str_split_i(pattern = "_", i = 1)
        current_celltype <- current_modality %>% str_split_i(pattern = "_", i = 2)
        
        current_bins_below_threshold <-
          ba_whole_genome_anno_chr_all %>%
          as_tibble() %>%
          filter(!is.na(get(str_c("bayes_mean_", current_modality))),
                 get(str_c("bayes_mean_", current_modality)) != 0) %>%
          filter(get(current_feature) == TRUE) %>%
          pull(str_c("bayes_mean_", current_modality)) %>%
          `<`(current_var_increment) %>%
          sum()
        
        current_valid_bins <-
          ba_whole_genome_anno_chr_all %>%
          as_tibble() %>%
          filter(!is.na(get(str_c("bayes_mean_", current_modality))),
                 get(str_c("bayes_mean_", current_modality)) != 0) %>%
          filter(get(current_feature) == TRUE) %>%
          pull(str_c("bayes_mean_", current_modality)) %>%
          length()
        
        tibble(mark = current_mark,
               celltype = current_celltype,
               feature_simplified = current_feature %>% str_remove("_pellacani$"),
               threshold = current_var_increment,
               fraction_of_bins_below_threshold = current_bins_below_threshold/current_valid_bins) %>%
          mutate(mark = mark %>% factor(vector_marks),
                 celltype = celltype %>% factor(vector_celltypes),
                 feature_simplified = feature_simplified %>% factor(c("all", vector_important_features_simplified)))
        
      }, mc.cores = detectCores()/(length(vector_modalities)*length(vector_important_features_v2))) %>% do.call("bind_rows", .)
    }, mc.cores = length(vector_important_features_v2)) %>% do.call("bind_rows", .)
  }, mc.cores = length(vector_modalities)) %>% do.call("bind_rows", .)

# Pause
message("Done 3/4")
gc()
Sys.sleep(10)

# Add feature-specific posiitons
ba_whole_genome_anno_chr_all_var_increments_features2 <-
  mclapply(vector_modalities, function(current_modality) {
    mclapply(vector_important_features_v2, function(current_feature) {
      mclapply(var_increments2, function(current_var_increment) {
        
        current_mark <-     current_modality %>% str_split_i(pattern = "_", i = 1)
        current_celltype <- current_modality %>% str_split_i(pattern = "_", i = 2)
        
        current_bins_below_threshold <-
          ba_whole_genome_anno_chr_all %>%
          as_tibble() %>%
          filter(!is.na(get(str_c("bayes_mean_", current_modality))),
                 get(str_c("bayes_mean_", current_modality)) != 0) %>%
          filter(get(current_feature) == TRUE) %>%
          pull(str_c("bayes_mean_", current_modality)) %>%
          `<`(current_var_increment) %>%
          sum()
        
        current_valid_bins <-
          ba_whole_genome_anno_chr_all %>%
          as_tibble() %>%
          filter(!is.na(get(str_c("bayes_mean_", current_modality))),
                 get(str_c("bayes_mean_", current_modality)) != 0) %>%
          filter(get(current_feature) == TRUE) %>%
          pull(str_c("bayes_mean_", current_modality)) %>%
          length()
        
        tibble(mark = current_mark,
               celltype = current_celltype,
               feature_simplified = current_feature %>% str_remove("_pellacani$"),
               threshold = current_var_increment,
               fraction_of_bins_below_threshold = current_bins_below_threshold/current_valid_bins) %>%
          mutate(mark = mark %>% factor(vector_marks),
                 celltype = celltype %>% factor(vector_celltypes),
                 feature_simplified = feature_simplified %>% factor(c("all", vector_important_features_simplified)))
        
      }, mc.cores = detectCores()/(length(vector_modalities)*length(vector_important_features_v2))) %>% do.call("bind_rows", .)
    }, mc.cores = length(vector_important_features_v2)) %>% do.call("bind_rows", .)
  }, mc.cores = length(vector_modalities)) %>% do.call("bind_rows", .)

# Combine
ba_whole_genome_anno_chr_all_var_increments_features <- 
  bind_rows(ba_whole_genome_anno_chr_all_var_increments_features1, ba_whole_genome_anno_chr_all_var_increments_features2) %>% 
  arrange(mark, celltype, threshold)

rm(ba_whole_genome_anno_chr_all_var_increments_features1, ba_whole_genome_anno_chr_all_var_increments_features2)

# Pause
message("Done 4/4")
gc()
Sys.sleep(10)

# Save objects
saveRDS(ba_whole_genome_anno_chr_all_var_increments_all, 
        "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_var_increments_all.RDS")

saveRDS(ba_whole_genome_anno_chr_all_var_increments_features, 
        "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_var_increments_features.RDS")

