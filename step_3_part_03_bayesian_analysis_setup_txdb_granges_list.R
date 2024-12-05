#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu

# Load constants
source("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/source/setup.R")

###########################################################################
###########################################################################
# Load in UCSC knownGene features
# Read in standard genomic features
# Add enhancers
txdb_grlist <- 
  genFeatures(TxDb.Hsapiens.UCSC.hg38.knownGene, 
              featuretype = "all", 
              reduce_ranges = TRUE, 
              upstream = 2000, 
              downstream = 200)

txdb_grlist$enhancer_encode <- 
  read_bed("/projects/edcc_new/reference_epigenomes/housekeeping/EDCCProd/Pipes/common/hg38/encode_enhancers_liftover.bed", genome_info = "hg38") %>% 
  keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse") %>% 
  reduce_ranges()

txdb_grlist$enhancer_pellacani <- 
  readxl::read_excel("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/data/1-s2.0-S2211124716314784-mmc6.xlsx") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                           ignore.strand = TRUE,
                           starts.in.df.are.0based = TRUE,
                           seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>%
  keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse") %>%
  reduce_ranges()

names(txdb_grlist) <- c("transcript", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic", "enhancer_encode", "enhancer_pellacani")

# Initialize ba-specific list
txdb_grlist_ba <- list()
# Clean up GRanges into ba-specific list
for (genome_feature in names(txdb_grlist)) { 
  txdb_grlist_ba[[genome_feature]] <- 
    txdb_grlist[[genome_feature]] %>% 
    keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse") %>% 
    sort() %>% 
    plyranges::select(-everything()) %>% 
    reduce_ranges() %>% 
    plyranges::mutate("{genome_feature}" := TRUE)
}

rm(txdb_grlist)