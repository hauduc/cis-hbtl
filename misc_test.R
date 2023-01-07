# test

#!/gsc/software/linux-x86_64-centos7/R-4.1.0/lib64/R/bin/R
# R 4.1.0
# x86_64-centos7-linux-gnu
library(plyranges)
library(VariantAnnotation)
library(genomation)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(tidyverse)

const_canonical_chromosomes <- str_c("chr", c(1:22, "X", "Y"))

##############################################################################################################################################################################################
##############################################################################################################################################################################################

# cell line

seqlevels(cl_variant_calls$mcf10a, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(cl_variant_calls$mcf10a) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(cl_variant_calls$htert_l9, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(cl_variant_calls$htert_l9) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(cl_variant_calls$htert_l2, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(cl_variant_calls$htert_l2) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)


# individuals

seqlevels(list_variant_calls$Individual_14_17_LP_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_14_17_LP_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_14_18_BC__356K_RNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_14_18_BC__356K_RNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_11_18_BC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_11_18_BC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_22_18_BC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_22_18_BC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_24_18_BC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_24_18_BC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_38_18_LC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_38_18_LC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_30_18_LC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_30_18_LC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

seqlevels(list_variant_calls$Individual_15_18_LC_RNA_DNA, pruning.mode = "coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(list_variant_calls$Individual_15_18_LC_RNA_DNA) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)



#######################################################################################################################################################
## SPAAs
# Load in all SPAAs
consensus.rpkm.bed.RData.files <- list.files(path = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/2/output2/spaas/annotated_rpkm", 
                                             pattern = "*RData$", 
                                             full.names = TRUE)

consensus.all.spaas.rpkm.bed <- data.frame()
for (i in seq_along(consensus.rpkm.bed.RData.files)) {
  load(consensus.rpkm.bed.RData.files[i])
  consensus.all.spaas.rpkm.bed <- bind_rows(consensus.all.spaas.rpkm.bed, consensus.spaa.rpkm.bed)
  rm(consensus.spaa.rpkm.bed)
}

# Fix and convert to GRanges
consensus.all.spaas.rpkm.bed <- consensus.all.spaas.rpkm.bed %>% 
  rowwise() %>% 
  mutate(RPKM_total    = sum(c_across(Individual_14_17_LP_DNA_RPKM:Individual_15_18_LC_RNA_DNA_RPKM)),
         RPKM_variance = var(c_across(Individual_14_17_LP_DNA_RPKM:Individual_15_18_LC_RNA_DNA_RPKM))) %>% 
  as.data.frame() %>% 
  mutate(rs = na_if(rs, "NA")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                           ignore.strand = TRUE, 
                           starts.in.df.are.0based = TRUE,
                           seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))





