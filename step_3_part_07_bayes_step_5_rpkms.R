#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

# Load constants
source("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/source/setup.R")

################################################################################################################################################################
################################################################################################################################################################
# Load gene information to test allelic balance across
biomart_gene_info <- readRDS("objects/biomart_gene_info.RDS")

# Make table of gene rpkms for all modalities
# Using your filelists, pulling out each .rpkm file, and then joining them together with bind_cols
ba_rpkms <- 
  mclapply(vector_celltypes, function(current_celltype) {
    mclapply(seq_along(vector_individuals), function(current_individual_i) {
      
      read_tsv(str_c("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/filelists/rpkm/breast_normal_rpkms_", current_celltype, ".tsv"), 
               col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>% 
        pull() %>% # transform into vector
        `[`(current_individual_i) %>% # pick filepath for current individual
        read_tsv(col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>% # read in rpkm file of that file path
        select(1, 3) %>% # pick ensembl gene ids and rpkms
        `colnames<-`(c("ensembl_gene_id", "rpkm")) %>% 
        left_join(as_tibble(mcols(biomart_gene_info)), ., by = join_by(ensembl_gene_id)) %>% # left join by ensembl_gene_id to the mcols of biomart gene info
        select(rpkm) %>% # pull just rpkms from current individual-modality
        rename("{vector_individuals_encoded[current_individual_i]}_{current_celltype}" := rpkm)
      
    }, mc.cores = length(seq_along(vector_individuals))) %>% do.call("bind_cols", .)
  }, mc.cores = length(vector_celltypes)) %>% do.call("bind_cols", .)

# Calculate the relative standard deviations for each ENSG
ba_rpkms_stats <- 
  mclapply(vector_celltypes, function(current_celltype) { 
    
    ba_rpkms %>%
      rowwise %>% 
      mutate("rpkm_mean_{current_celltype}" := mean(c_across(ends_with(str_c("_", current_celltype)))),
             "rpkm_rsd_{current_celltype}"  := 100 * sd(c_across(ends_with(str_c("_", current_celltype)))) / mean(c_across(ends_with(str_c("_", current_celltype))))) %>% 
      ungroup() %>% 
      select(last_col(1), last_col(0))
    
  }, mc.cores = length(vector_celltypes)) %>% do.call("bind_cols", .)

# Use biomart data to get gene positions
# Set new mcols to have gene identity columns + rpkm values from all individuals
# Also filter out NA genes and genes where mean expression across a cell type is below 1 RPKM by creating your temporary mean_ variables for each celltype
ba_rpkms_gr <- 
  biomart_gene_info %>% 
  as_tibble() %>% 
  bind_cols(ba_rpkms_stats, ba_rpkms) %>% 
  filter(seqnames %in% const_canonical_autosomes) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
  keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse")

# Process gene promoters
# Create ENSG-ENST conversion key
# Remove value after dot in ENSG/ENST names
# Drop na and keep only unique genes
ba_ensg_enst_key <- 
  read_gff3("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/data/gencode/gencode.v43.basic.annotation.gff3.gz") %>% 
  as_tibble() %>% 
  select(gene_id, transcript_id) %>% 
  modify(function(x) x %>% str_split("\\.", simplify = TRUE) %>% as.data.frame() %>% pull(1)) %>% 
  drop_na() %>% 
  distinct(gene_id, .keep_all = TRUE)

# Get promoter positions with corresponding ensgs attached
ba_promoter_positions <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene %>% 
  promoters(upstream = 2000, downstream = 200) %>% 
  as_tibble() %>% 
  modify_at(.at = c("tx_name"), function(x) x %>% str_split("\\.", simplify = TRUE) %>% as.data.frame() %>% pull(1)) %>% 
  left_join(ba_ensg_enst_key, by = join_by(tx_name == transcript_id)) %>% 
  filter(seqnames %in% const_canonical_autosomes)

# Replace gene positions with promoter positions in ba_rpkms_gr
ba_promoters_rpkm_anno <-
  ba_rpkms_gr %>% 
  mcols() %>% 
  as_tibble() %>% 
  left_join(x = ba_promoter_positions, y = ., by = join_by(gene_id == ensembl_gene_id)) %>% 
  rename(ensembl_gene_id = gene_id) %>% 
  select(-width, -strand) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                           ignore.strand = TRUE,
                           seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>% 
  keepSeqlevels(const_canonical_autosomes, pruning.mode = "coarse") 

# Repeat steps from above but instead with your promoter data
# Add bayesian values to corresponding gene promoters
# Remove blacklisted regions (for all) and only keep ensembl_gene_id, expression sd, and bayesian values for each modality
# ba_whole_genome_anno_chr_all_rpkms_promoters has gene promoter coordinates and rsd of gene expression for that gene, 
# along with bayesian values for each modality along the length of the promoter
ba_promoters_rpkm_bayes_anno <- 
  ba_whole_genome_anno_chr_all %>% 
  join_overlap_left(ba_promoters_rpkm_anno) %>% 
  as_tibble()

# Standard filter expressions
# Remove zero- and low-expressing genes
ba_promoters_rpkm_bayes_anno_clean <- 
  ba_promoters_rpkm_bayes_anno %>% 
  filter(if_any(starts_with("bayes_category_"), function(x) x != "blacklisted")) %>% # Remove any blacklisted regions
  select(ensembl_gene_id, 
         starts_with("rpkm_mean_"), 
         starts_with("rpkm_rsd_"), 
         starts_with("bayes_mean_")) %>% # Only pick essential colummns
  drop_na() %>% # Remove NAs
  filter(if_all(starts_with("rpkm_mean_"), function(x) x > 1)) # Remove genes where all RPKMs are below 1

# For all genes:
# Summarize rsd and the mean bayes value as it appears across the promoter
# Helper function for processing vector of bayes values across promoter
get_mean_of_nonzero_values <- 
  function(current_bins_vector) {
    # Calculate mean on all nonzero bins
    mean(current_bins_vector[current_bins_vector != 0]) %>% if_else(is.nan(.), 0, .)
  }

get_proportion_of_values_that_are_zero <- 
  function(current_bins_vector) {
    # Get proportion of 0 bins
    length(current_bins_vector[current_bins_vector == 0])/length(current_bins_vector)
  }

# Execute on all
ba_promoters_rpkm_bayes_anno_clean_summarized <- 
  mclapply(vector_modalities, function(current_modality) { 
    
    # Get the strings for modality
    current_mark <- current_modality %>% str_split_i(pattern = "_", i = 1)
    current_celltype <- current_modality %>% str_split_i(pattern = "_", i = 2)
    current_mean_col <- str_c("rpkm_mean_", current_celltype)
    current_rsd_col <- str_c("rpkm_rsd_", current_celltype)
    current_bayes_mean_col <- str_c("bayes_mean_", current_modality)
    
    # Main function
    ba_promoters_rpkm_bayes_anno_clean %>% 
      select(ensembl_gene_id, 
             matches(current_mean_col),
             matches(current_rsd_col), 
             matches(current_bayes_mean_col)) %>%
      group_by(ensembl_gene_id) %>% 
      summarize(gene_rpkm_mean                       = .data[[current_mean_col]] %>% unique(), 
                gene_rpkm_rsd                        = .data[[current_rsd_col]]  %>% unique(), 
                gene_prom_prop_of_zero_bins          = .data[[current_bayes_mean_col]] %>% get_proportion_of_values_that_are_zero(),
                gene_prom_number_of_bins             = .data[[current_bayes_mean_col]] %>% length(),
                gene_prom_mean_theta_raw             = .data[[current_bayes_mean_col]] %>% mean(),
                gene_prom_mean_theta_of_nonzero_bins = .data[[current_bayes_mean_col]] %>% get_mean_of_nonzero_values()) %>% 
      # mutate(gene_prom_mean_theta_raw_quartile =
      #          case_when((gene_prom_mean_theta_raw == 0)                                        ~ ba_variability_categories[2],
      #                    (gene_prom_mean_theta_raw > 0)    & (gene_prom_mean_theta_raw <= 0.25) ~ ba_variability_categories[3],
      #                    (gene_prom_mean_theta_raw > 0.25) & (gene_prom_mean_theta_raw <= 0.5)  ~ ba_variability_categories[4],
      #                    (gene_prom_mean_theta_raw > 0.5)  & (gene_prom_mean_theta_raw <= 0.75) ~ ba_variability_categories[5],
      #                    (gene_prom_mean_theta_raw > 0.75)                                      ~ ba_variability_categories[6]) %>% factor(ba_variability_categories),
      #        gene_prom_mean_theta_of_nonzero_bins_quartile =
      #          case_when((gene_prom_mean_theta_of_nonzero_bins == 0)                                                    ~ ba_variability_categories[2],
      #                    (gene_prom_mean_theta_of_nonzero_bins > 0)    & (gene_prom_mean_theta_of_nonzero_bins <= 0.25) ~ ba_variability_categories[3],
      #                    (gene_prom_mean_theta_of_nonzero_bins > 0.25) & (gene_prom_mean_theta_of_nonzero_bins <= 0.5)  ~ ba_variability_categories[4],
      #                    (gene_prom_mean_theta_of_nonzero_bins > 0.5)  & (gene_prom_mean_theta_of_nonzero_bins <= 0.75) ~ ba_variability_categories[5],
    #                    (gene_prom_mean_theta_of_nonzero_bins > 0.75)                                                  ~ ba_variability_categories[6]) %>% factor(ba_variability_categories)) %>%
    mutate(mark = current_mark %>% factor(vector_marks), 
           celltype = current_celltype %>% factor(vector_celltypes),
           .before = 1)
    
  }, mc.cores = length(vector_modalities)) %>% do.call("bind_rows", .)

# Plot promoters more selectively
# Dotplot
ba_promoters_rpkm_bayes_anno_clean_summarized %>% 
  filter(
    gene_prom_number_of_bins >= 44, # Keep only promoters that haven't been chopped by blacklisting or other filtration (44 or 45 bins = 2,200 bp)
    gene_prom_prop_of_zero_bins < 0.5,       # Keep only promoters that are not mostly empty
  ) %>% 
  ggplot(aes(x = gene_prom_mean_theta_of_nonzero_bins, y = gene_rpkm_rsd, color = celltype)) +
  geom_point(size = 0.5) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", x))) +
  scale_color_manual(values = vector_celltype_colors) +
  facet_wrap(vars(mark), nrow = 1) +
  theme_bw() +
  labs(x = "mean histone 𝜃 estimate at gene promoter", 
       y = "gene RPKM relative standard deviation",
       color = "cell type")
ggsave(str_c("plots/ba_promoters_rpkm_bayes_anno_clean_summarized", "_", "v5", ".png"),
       units = "in",
       width  = 5.5*2.2,
       height = 5.5*0.6)

# ba_whole_genome_anno_chr_all <- readRDS("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_whole_genome_anno_chr_all.RDS")
saveRDS(ba_rpkms, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_rpkms.RDS")
saveRDS(ba_rpkms_stats, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_rpkms_rsds.RDS")
saveRDS(ba_rpkms_gr, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_rpkms_gr.RDS")
saveRDS(ba_ensg_enst_key, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_ensg_enst_key.RDS")
saveRDS(ba_gene_promoter_positions, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_gene_promoter_positions.RDS")
saveRDS(ba_promoters_rpkm_bayes_anno, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_promoters_rpkm_bayes_anno.RDS")
saveRDS(ba_promoters_rpkm_bayes_anno_clean, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_promoters_rpkm_bayes_anno_clean.RDS")
saveRDS(ba_promoters_rpkm_bayes_anno_clean_summarized, "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/objects/ba_promoters_rpkm_bayes_anno_clean_summarized.RDS")









