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
# Non-exclusive (to cell types)
# Initialize
ba_great_enrich_categories_nonexclusive <- tibble()

# Loop
for (current_mark in vector_marks) { 
  for (current_celltype in vector_celltypes) { 
    for (current_ba_var_cat in ba_variability_categories_nonabsent) { 
      
      # Get current relevant granges
      # And add to results
      current_tibble <- 
        ba_whole_genome_anno_chr_all %>% 
        plyranges::filter(get(str_c("bayes_category_", current_mark, "_", current_celltype)) == current_ba_var_cat) %>% 
        reduce_ranges() %>%
        submitGreatJob(gr = ., species = "hg38", request_interval = 1) %>% 
        getEnrichmentTables() %>% 
        map(function(x) x %>% mutate(rank = seq_along(rownames(.)), .before = 1)) %>% 
        bind_rows(.id = "go_domain") %>% 
        as_tibble() %>% 
        mutate(mark = current_mark %>% factor(vector_marks),
               celltype = current_celltype %>% factor(vector_celltypes),
               bayesian_bin_type = current_ba_var_cat %>% factor(ba_variability_categories),
               .before = 1)
      
      # Add rows
      ba_great_enrich_categories_nonexclusive <- bind_rows(ba_great_enrich_categories_nonexclusive, current_tibble)
      
      # Progress
      message(str_c(timestamp(quiet = TRUE), "Done with nonexclusive", current_mark, current_celltype, current_ba_var_cat, sep = " "))
    }
  }
}

###########################################################################
# Exlusive to cell types
# Initialize tibble
ba_great_enrich_categories_celltype_exclusive <- tibble()

# Loop
for (current_mark in vector_marks) { 
  for (current_celltype in vector_celltypes) { 
    for (current_ba_var_cat in ba_variability_categories_nonabsent) { 
      
      # Get current excluded cell types from your current cell type
      current_excluded_celltypes <- 
        vector_celltypes %>% 
        str_subset(current_celltype, negate = TRUE)
      
      # Get current relevant mark, celltype, and variability category ranges
      # And exclude ranges of the same variability category from other cell types
      # In a single filtering expression
      # Then submit to GREAT
      current_tibble <- 
        ba_whole_genome_anno_chr_all %>% 
        plyranges::filter((get(str_c("bayes_category_", current_mark, "_", current_celltype)) == current_ba_var_cat) &
                            ((get(str_c("bayes_category_", current_mark, "_", current_excluded_celltypes[1])) != current_ba_var_cat) &
                               (get(str_c("bayes_category_", current_mark, "_", current_excluded_celltypes[2])) != current_ba_var_cat) &
                               (get(str_c("bayes_category_", current_mark, "_", current_excluded_celltypes[3])) != current_ba_var_cat))) %>% 
        reduce_ranges() %>% 
        submitGreatJob(gr = ., species = "hg38", request_interval = 1) %>% 
        getEnrichmentTables() %>% 
        map(function(x) x %>% mutate(rank = seq_along(rownames(.)), .before = 1)) %>% 
        bind_rows(.id = "go_domain") %>% 
        as_tibble() %>% 
        mutate(mark = current_mark %>% factor(vector_marks),
               celltype = current_celltype %>% factor(vector_celltypes),
               bayesian_bin_type = current_ba_var_cat %>% factor(ba_variability_categories),
               .before = 1)
      
      # Add rows
      ba_great_enrich_categories_celltype_exclusive <- bind_rows(ba_great_enrich_categories_celltype_exclusive, current_tibble)
      
      # Progress
      message(str_c(timestamp(quiet = TRUE), "Done with celltype exclusive", current_mark, current_celltype, current_ba_var_cat, sep = " "))
    }
  }
}

# Save objects
saveRDS(ba_great_enrich_categories_nonexclusive, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_great_enrich_categories_nonexclusive.RDS")
saveRDS(ba_great_enrich_categories_celltype_exclusive, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_great_enrich_categories_celltype_exclusive.RDS")


