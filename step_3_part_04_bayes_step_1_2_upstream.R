#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

# Progress
message(timestamp(quiet = TRUE), " Load constants start")

# Load constants
source("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/source/setup.R")
# Load features
txdb_grlist_ba <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/txdb_grlist_ba.RDS")

# Progress
message(timestamp(quiet = TRUE), " Load constants done, start function set")

###########################################################################
###########################################################################
# Set function
create_bayes_whole_genome_gr_chromwise <- function(current_chromosome) {
  
  # 1 -----------------------------------------------------------------------
  # Set up empty annotated whole genome with positions of genomic features
  ba_whole_genome_anno_chr <- 
    read_tsv("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/refs/chrom_info/hg38_no_alt.chromsizes.txt", 
             col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>% 
    rename("seqnames" = 1, "end" = 2) %>% 
    mutate(start = 1, .after = seqnames) %>% 
    mutate(strand = "*") %>% 
    mutate(end = floor(end/50)*50) %>% 
    filter(seqnames %in% current_chromosome) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = FALSE, 
                             ignore.strand = TRUE, 
                             starts.in.df.are.0based = FALSE,
                             seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
    keepSeqlevels(current_chromosome, pruning.mode = "coarse") %>%
    sort() %>% 
    tile_ranges(50) %>% 
    plyranges::select(-partition)
  
  # Pre-process txdb_grlist_ba
  for (genome_feature in names(txdb_grlist_ba)) {
    txdb_grlist_ba[[genome_feature]] <-
      txdb_grlist_ba[[genome_feature]] %>% 
      keepSeqlevels(current_chromosome, pruning.mode = "coarse")
  }
  
  # Loop through the different sequence feature types
  # Join left and annotate TRUE if present
  for (genome_feature in names(txdb_grlist_ba)) { 
    ba_whole_genome_anno_chr <-
      ba_whole_genome_anno_chr %>% 
      join_overlap_left(txdb_grlist_ba[[genome_feature]])
  }
  # Then change features not present as NA to FALSE
  mcols(ba_whole_genome_anno_chr) <- mcols(ba_whole_genome_anno_chr) %>% map_df(function(x) replace_na(x, FALSE))
  
  # 2 -----------------------------------------------------------------------
  # Loop through marks and create sublist for each mark in original object
  ba_prob_tables_df_chr <- 
    mclapply(vector_modalities, function(current_modality) {
      
      # Set current variables
      current_mark <-     current_modality %>% str_split_i("_", 1)
      current_celltype <- current_modality %>% str_split_i("_", 2)
      
      # Create (empty) tibble to hold final ba probs for each mark + celltype combination
      # Turn raw tibble into GRanges and add into list object 
      # Cleanly read in chromosome-level probability tables as tibble
      # Remove the 50bp using the head command so everything fits within chromosome sizes
      # Make prob tables as a long data frame
      # Create indexed tibble with data from all objects in ba_prob_tables list
      # Parallel implementation of coalescing list into single dataframe, indexed by information columns, for graphing
      str_c("/projects/epigenomics_assembly/IA/mammary_analysis/wasp_prob_tables_failed_removed_prior_updated/bayes_prob_table_", 
            current_celltype, "_Healthy_", current_mark, "_", current_chromosome, ".gz") %>% 
        read_csv(progress = FALSE, show_col_types = FALSE) %>%
        select(1, 2, 3) %>%
        rename(chr = 1, start = 2, bayes_mean = 3) %>% 
        mutate(end = start + 50, .after = start) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE, 
                                 starts.in.df.are.0based = TRUE,
                                 seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
        keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse") %>%
        trim() %>% 
        sort() %>% 
        as_tibble() %>% 
        mutate(position = seq_along(rownames(.)),
               mark = current_mark, 
               celltype = current_celltype, 
               .after = strand) %>% 
        mutate(mark = mark %>% factor(vector_marks),
               celltype = celltype %>% factor(vector_celltypes))
      
    }, mc.cores = floor(detectCores()/length(const_canonical_autosomes))) %>% do.call("bind_rows", .)
  
  # 3 -----------------------------------------------------------------------
  # Add values of different experiment types to your tibble
  # Overlap new column, which is called bayes_mean in all samples, then rename the columna and remove old col
  # Have to do rather than rename due to limits of plyranges
  for (current_mark in vector_marks) {
    for (current_celltype in vector_celltypes) { 
      ba_whole_genome_anno_chr <-
        ba_whole_genome_anno_chr %>% 
        join_overlap_left(ba_prob_tables_df_chr %>% 
                            filter(mark == current_mark, celltype == current_celltype) %>% 
                            select(-width, -strand, -position, -mark, -celltype) %>% 
                            makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
                            keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse")) %>% 
        plyranges::mutate("bayes_mean_{current_mark}_{current_celltype}" := bayes_mean) %>% 
        plyranges::select(-bayes_mean)
    }
  }
  
  # Add category annotations to ba_whole_genome_anno_chr
  ba_whole_genome_anno_chr_categories <- 
    mclapply(vector_modalities, function(current_modality) {
      ba_whole_genome_anno_chr %>% 
        as_tibble() %>% 
        mutate("bayes_category_{current_modality}" := 
                 case_when(is.na(get(str_c("bayes_mean_", current_modality)))                                                           ~ ba_variability_categories[1],
                           (get(str_c("bayes_mean_", current_modality)) == 0)                                                           ~ ba_variability_categories[2],
                           (get(str_c("bayes_mean_", current_modality)) > 0)    & (get(str_c("bayes_mean_", current_modality)) <= 0.25) ~ ba_variability_categories[3],
                           (get(str_c("bayes_mean_", current_modality)) > 0.25) & (get(str_c("bayes_mean_", current_modality)) <= 0.5)  ~ ba_variability_categories[4],
                           (get(str_c("bayes_mean_", current_modality)) > 0.5)  & (get(str_c("bayes_mean_", current_modality)) <= 0.75) ~ ba_variability_categories[5],
                           (get(str_c("bayes_mean_", current_modality)) > 0.75)                                                         ~ ba_variability_categories[6])) %>% 
        select(last_col())
    }, mc.cores = floor(detectCores()/length(const_canonical_autosomes))) %>% do.call("bind_cols", .)
  
  # Assign new mcols
  mcols(ba_whole_genome_anno_chr) <- as_tibble(mcols(ba_whole_genome_anno_chr)) %>% bind_cols(ba_whole_genome_anno_chr_categories)
  
  # Clean up seqinfo
  seqlevels(ba_whole_genome_anno_chr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
  seqinfo(ba_whole_genome_anno_chr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  ba_whole_genome_anno_chr <- ba_whole_genome_anno_chr %>% keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse")
  
  # Remove temp objects
  rm(ba_whole_genome_anno_chr_categories)
  
  # Return chromosome-specific ba_whole_genome_anno
  ba_whole_genome_anno_chr
} 

# Progress
message(timestamp(quiet = TRUE), " Set function done, start mclapply")

###########################################################################
###########################################################################
# Run
ba_whole_genome_anno_chr_all <- 
  mclapply(const_canonical_autosomes, 
           create_bayes_whole_genome_gr_chromwise,
           mc.cores = length(const_canonical_autosomes)) %>% 
  do.call("c", .)

# Save object
saveRDS(ba_whole_genome_anno_chr_all, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all.RDS")
# Progress
message(timestamp(quiet = TRUE), " Save done")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Create v2 where Jonathan 0s are replaced by NAs (update 2024-05-14)
ba_whole_genome_anno_chr_all_zero_rm <- 
  ba_whole_genome_anno_chr_all %>% 
  as_tibble() %>% 
  map_at(.at = vars(starts_with("bayes_mean_")), .f = function(x) ifelse(x == 0, NA, x)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                           seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
  keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse") 

# Progress
message(timestamp(quiet = TRUE), "Start save")
# Save object
saveRDS(ba_whole_genome_anno_chr_all_zero_rm, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all_zero_rm.RDS")
# Progress
message(timestamp(quiet = TRUE), " Save done")
