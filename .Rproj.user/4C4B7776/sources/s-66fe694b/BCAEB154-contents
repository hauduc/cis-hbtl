#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu
library(plyranges)
library(VariantAnnotation)
library(MutationalPatterns)
library(genomation)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(systemPipeR)
library(biomaRt)
library(LDlinkR)
library(AllelicImbalance)
library(broom)
library(gt)
library(gtsummary)
library(UpSetR)
library(tidyverse)
library(ggdendro)

const_canonical_chromosomes <- str_c("chr", c(1:22, "X"))

vector_marks <- c("H3K4me3", "H3K4me1", "H3K27ac", "H3K9me3", "H3K27me3", "H3K36me3")
vector_celltypes <- c("BC", "LP", "LC", "SC")
vector_individuals <- c("Individual_14_17_LP_DNA",
                        "Individual_14_18_BC__356K_RNA",
                        "Individual_11_18_BC_RNA_DNA",
                        "Individual_22_18_BC_RNA_DNA",
                        "Individual_24_18_BC_RNA_DNA",
                        "Individual_38_18_LC_RNA_DNA",
                        "Individual_30_18_LC_RNA_DNA",
                        "Individual_15_18_LC_RNA_DNA")

vector_individuals_encoded <- c("Individual_1",
                                "Individual_2",
                                "Individual_3",
                                "Individual_4",
                                "Individual_5",
                                "Individual_6",
                                "Individual_7",
                                "Individual_8")

######################################################################################################################################################################################
######################################################################################################################################################################################
cemt_table <- 
    tibble(cemt = list.files("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/1/output2/H3K4me3/", include.dirs = TRUE) %>% rep(6),
           individual_number = c(rep("1", 4), rep("2", 4), rep("3", 4), rep("4", 4), rep("5", 4), rep("6", 4), rep("7", 4), rep("8", 4)) %>% rep(6),
           mark = c(rep("H3K4me3", 8*4), rep("H3K27ac", 8*4), rep("H3K4me1", 8*4), rep("H3K9me3", 8*4), rep("H3K27me3", 8*4), rep("H3K36me3", 8*4)),
           celltype = rep(c("BC", "LP", "LC", "SC"), 192/4))

# Start with H3K36me3
cemt_table_H3K36me3 <- cemt_table %>% filter(mark == "H3K36me3")
# Modules folder path
modules_path <- "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules"
# Initialize results tibble
ase_depths_H3K36me3 <- tibble()

for (i in seq_along(cemt_table_H3K36me3$cemt)) {
    # Import ChIP-seq treatment vcf file
    ase_treatment_vcf <- 
        readVcfAsVRanges(str_c(modules_path, "/notebook/alignment_analysis/ase/output/", cemt_table_H3K36me3$mark[i], "_", cemt_table_H3K36me3$celltype[i], "_", "Individual", "_", cemt_table_H3K36me3$individual_number[i], "_treatment_cis_hbtls.vcf.gz")) %>%
        as.data.frame() %>%
        select(-width, -strand) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>%
        keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse")
    
    # Scan VCF for heterozygotic variant calls and keep only cis-hBTLs at heterozygotic positions for that individual
    # These are the positions that need to be tested for the particular individual-modality
    ase_het_positions <-
        join_overlap_intersect(consensus.all.spaas.rpkm.bed %>% 
                                   plyranges::filter(mark == cemt_table_H3K36me3$mark[i], celltype == cemt_table_H3K36me3$celltype[i], spaa == TRUE) %>% 
                                   keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse"),
                               readVcfAsVRanges(str_c(modules_path, "/1/output2/", cemt_table_H3K36me3$mark[i], "/", cemt_table_H3K36me3$cemt[i], "/wgs/", cemt_table_H3K36me3$cemt[i], ".recal_hard_filtered_phased_snvs.vcf.gz")) %>%
                                   as.data.frame() %>%
                                   select(-width, -strand) %>%
                                   makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)) %>%
                                   keepSeqlevels(const_canonical_chromosomes, pruning.mode = "coarse") %>% 
                                   plyranges::filter(AF == 0.5)) %>% 
        plyranges::select(-everything())
    
    # Subset the ChIP-seq treatment VCF with proper heterozygotic, modality-relevant positions
    # And select only the positions with some evidence of heterozygotic alignment
    ase_treatment_vcf_subsetted <-
        ase_treatment_vcf %>% 
        join_overlap_intersect(ase_het_positions) %>% 
        plyranges::filter(refDepth != 0 & altDepth != 0)
    
    # Pull out the vectors of depths where ref or alt != 0
    current_tibble <-
        tibble(mark = cemt_table_H3K36me3$mark[i],
               celltype = cemt_table_H3K36me3$celltype[i],
               individual = cemt_table_H3K36me3$individual_number[i],
               ref_depth_mean = ase_treatment_vcf_subsetted$refDepth %>% mean(),
               alt_depth_mean = ase_treatment_vcf_subsetted$altDepth %>% mean())
    
    # Add to tibble
    ase_depths_H3K36me3 <- bind_rows(ase_depths_H3K36me3, current_tibble)
}
