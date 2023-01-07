# job visualization

#!/gsc/software/linux-x86_64-centos7/R-4.1.3/lib64/R/bin/R
# R 4.1.3
# x86_64-centos7-linux-gnu
# library(plyranges)
# library(VariantAnnotation)
# library(MutationalPatterns)
# library(genomation)
# library(Gviz)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(biomaRt)
# library(LDlinkR)
# library(tidyverse)
##############################################################################################################################################################################################
##############################################################################################################################################################################################

# load("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/genomic.RData")

###################################################################################################
###################################################################################################
###################################################################################################
current_rsID     <- "rs34258884"
current_mark     <- "H3K4me1"
current_celltype <- "LC"

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_example.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)

dev.off()

###################################################################################################
###################################################################################################
###################################################################################################
current_rsID     <- "rs34258884"
current_mark     <- "H3K4me1"
current_celltype <- "LC"

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_cl_example.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)

dev.off()

###################################################################################################
