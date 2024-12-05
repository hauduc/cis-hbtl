#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

# Load constants
source("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/source/setup.R")

# Initialize ba-specific list
txdb_grlist_ba <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/txdb_grlist_ba.RDS")

# Load ba object
ba_whole_genome_anno_chr_all <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all.RDS")

###########################################################################
###########################################################################
# Chrom positions tibble
# Run
# Function to calculate number of positions per chromosome
calculate_ba_whole_genome_anno_chr_all_num_positions_by_chr <-
  function(current_chromosome) {
    tibble(seqnames = current_chromosome,
           num_positions = ba_whole_genome_anno_chr_all %>% 
             plyranges::filter(seqnames == current_chromosome) %>% 
             length())
  }

# Run
ba_chrom_positions <- 
  mclapply(const_canonical_autosomes, calculate_ba_whole_genome_anno_chr_all_num_positions_by_chr, mc.cores = length(const_canonical_autosomes)) %>% 
  do.call("bind_rows", .)

# Calculate cumulative positions
ba_chrom_positions_cumulative <- tibble()

for(i in seq_along(ba_chrom_positions$seqnames)) {
  
  current_start <- ba_chrom_positions$num_positions[1:i-1] %>% sum() %>% `+`(1)
  
  current_end <- ba_chrom_positions$num_positions[1:i] %>% sum()
  
  current_ba_chrom_positions_cumulative <-
    tibble(position_start = current_start,
           position_end = current_end)
  
  ba_chrom_positions_cumulative <- bind_rows(ba_chrom_positions_cumulative, current_ba_chrom_positions_cumulative)
}

ba_chrom_positions <- bind_cols(ba_chrom_positions, ba_chrom_positions_cumulative)
###########################################################################
###########################################################################
# Promoter
# Pivot to longer form
# Chromwise
ba_whole_genome_anno_chr_all_long_promoter <- 
  mclapply(const_canonical_autosomes, function(current_chromosome) {
    
    ba_whole_genome_anno_chr_all %>% 
      as_tibble() %>% 
      filter(seqnames == current_chromosome) %>% 
      select(all_of("promoter"), starts_with("bayes_mean_")) %>%
      mutate(position = filter(ba_chrom_positions, seqnames == current_chromosome)$position_start:filter(ba_chrom_positions, seqnames == current_chromosome)$position_end, .before = 1) %>% 
      lazy_dt() %>% 
      rename("feature_present" = "promoter") %>% 
      pivot_longer(cols = starts_with("bayes_mean_"),
                   names_to = "modality",
                   values_to = "bayes_mean") %>%  
      mutate(mark = str_split_i(modality, pattern = "_", i = 3),
             celltype = str_split_i(modality, pattern = "_", i = 4)) %>%   
      select(-modality) %>% 
      mutate(bin_category = fifelse(is.na(bayes_mean), ba_variability_categories[1], 
                                    fifelse(bayes_mean == 0, ba_variability_categories[2],
                                            fifelse((bayes_mean > 0) & (bayes_mean <= 0.25), ba_variability_categories[3],
                                                    fifelse((bayes_mean > 0.25) & (bayes_mean <= 0.5), ba_variability_categories[4],
                                                            fifelse((bayes_mean > 0.5) & (bayes_mean <= 0.75), ba_variability_categories[5], ba_variability_categories[6])))))) %>% 
      mutate(feature = "promoter",
             feature_present = feature_present %>% factor(c(TRUE, FALSE)),
             mark = mark %>% factor(vector_marks),
             celltype = celltype %>% factor(vector_celltypes),
             bin_category = bin_category %>% factor(ba_variability_categories)) %>% 
      as.data.table()
    
  }, mc.cores = length(const_canonical_autosomes))

ba_whole_genome_anno_chr_all_long_promoter <- rbindlist(ba_whole_genome_anno_chr_all_long_promoter)

ba_whole_genome_anno_chr_all_summarized_promoter <- 
  ba_whole_genome_anno_chr_all_long_promoter %>% 
  group_by(mark, celltype, feature, bin_category, feature_present) %>% 
  summarize(n_bins = n()) %>%
  as_tibble()

saveRDS(ba_whole_genome_anno_chr_all_summarized_promoter, 
        "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_summarized_promoter.RDS")

gc()
###########################################################################
###########################################################################
# Enhancer
# Pivot to longer form
# Chromwise
ba_whole_genome_anno_chr_all_long_enhancer_pellacani <- 
  mclapply(const_canonical_autosomes, function(current_chromosome) {
    
    ba_whole_genome_anno_chr_all %>% 
      as_tibble() %>% 
      filter(seqnames == current_chromosome) %>% 
      select(all_of("enhancer_pellacani"), starts_with("bayes_mean_")) %>%
      mutate(position = filter(ba_chrom_positions, seqnames == current_chromosome)$position_start:filter(ba_chrom_positions, seqnames == current_chromosome)$position_end, .before = 1) %>% 
      lazy_dt() %>% 
      rename("feature_present" = "enhancer_pellacani") %>% 
      pivot_longer(cols = starts_with("bayes_mean_"),
                   names_to = "modality",
                   values_to = "bayes_mean") %>%  
      mutate(mark = str_split_i(modality, pattern = "_", i = 3),
             celltype = str_split_i(modality, pattern = "_", i = 4)) %>%   
      select(-modality) %>% 
      mutate(bin_category = fifelse(is.na(bayes_mean), ba_variability_categories[1], 
                                    fifelse(bayes_mean == 0, ba_variability_categories[2],
                                            fifelse((bayes_mean > 0) & (bayes_mean <= 0.25), ba_variability_categories[3],
                                                    fifelse((bayes_mean > 0.25) & (bayes_mean <= 0.5), ba_variability_categories[4],
                                                            fifelse((bayes_mean > 0.5) & (bayes_mean <= 0.75), ba_variability_categories[5], ba_variability_categories[6])))))) %>% 
      mutate(feature = "enhancer_pellacani",
             feature_present = feature_present %>% factor(c(TRUE, FALSE)),
             mark = mark %>% factor(vector_marks),
             celltype = celltype %>% factor(vector_celltypes),
             bin_category = bin_category %>% factor(ba_variability_categories)) %>% 
      as.data.table()
    
  }, mc.cores = length(const_canonical_autosomes))

ba_whole_genome_anno_chr_all_long_enhancer_pellacani <- rbindlist(ba_whole_genome_anno_chr_all_long_enhancer_pellacani)

ba_whole_genome_anno_chr_all_summarized_enhancer_pellacani <- 
  ba_whole_genome_anno_chr_all_long_enhancer_pellacani %>% 
  group_by(mark, celltype, feature, bin_category, feature_present) %>% 
  summarize(n_bins = n()) %>%
  as_tibble()

saveRDS(ba_whole_genome_anno_chr_all_summarized_enhancer_pellacani, 
        "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_summarized_enhancer_pellacani.RDS")

gc()
###########################################################################
###########################################################################
# Transcript
# Pivot to longer form
# Chromwise
ba_whole_genome_anno_chr_all_long_transcript <- 
  mclapply(const_canonical_autosomes, function(current_chromosome) {
    
    ba_whole_genome_anno_chr_all %>% 
      as_tibble() %>% 
      filter(seqnames == current_chromosome) %>% 
      select(all_of("transcript"), starts_with("bayes_mean_")) %>%
      mutate(position = filter(ba_chrom_positions, seqnames == current_chromosome)$position_start:filter(ba_chrom_positions, seqnames == current_chromosome)$position_end, .before = 1) %>% 
      lazy_dt() %>% 
      rename("feature_present" = "transcript") %>% 
      pivot_longer(cols = starts_with("bayes_mean_"),
                   names_to = "modality",
                   values_to = "bayes_mean") %>%  
      mutate(mark = str_split_i(modality, pattern = "_", i = 3),
             celltype = str_split_i(modality, pattern = "_", i = 4)) %>%   
      select(-modality) %>% 
      mutate(bin_category = fifelse(is.na(bayes_mean), ba_variability_categories[1], 
                                    fifelse(bayes_mean == 0, ba_variability_categories[2],
                                            fifelse((bayes_mean > 0) & (bayes_mean <= 0.25), ba_variability_categories[3],
                                                    fifelse((bayes_mean > 0.25) & (bayes_mean <= 0.5), ba_variability_categories[4],
                                                            fifelse((bayes_mean > 0.5) & (bayes_mean <= 0.75), ba_variability_categories[5], ba_variability_categories[6])))))) %>% 
      mutate(feature = "transcript",
             feature_present = feature_present %>% factor(c(TRUE, FALSE)),
             mark = mark %>% factor(vector_marks),
             celltype = celltype %>% factor(vector_celltypes),
             bin_category = bin_category %>% factor(ba_variability_categories)) %>% 
      as.data.table()
    
  }, mc.cores = length(const_canonical_autosomes))

ba_whole_genome_anno_chr_all_long_transcript <- rbindlist(ba_whole_genome_anno_chr_all_long_transcript)

ba_whole_genome_anno_chr_all_summarized_transcript <- 
  ba_whole_genome_anno_chr_all_long_transcript %>% 
  group_by(mark, celltype, feature, bin_category, feature_present) %>% 
  summarize(n_bins = n()) %>%
  as_tibble()

saveRDS(ba_whole_genome_anno_chr_all_summarized_transcript, 
        "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_summarized_transcript.RDS")

gc()
###########################################################################
###########################################################################
