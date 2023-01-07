# job 1

#!/gsc/software/linux-x86_64-centos7/R-4.1.0/lib64/R/bin/R
# R 4.1.0
# x86_64-centos7-linux-gnu
library(plyranges)
library(VariantAnnotation)
library(MutationalPatterns)
library(genomation)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(LDlinkR)
library(tidyverse)
##############################################################################################################################################################################################
##############################################################################################################################################################################################

# load("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/notebook/alignment_analysis/genomic.RData")

###################################################################################################
###################################################################################################
###################################################################################################
# current_rsID <-           top_hits_ranges$rs[1]
# current_mark <-         top_hits_ranges$mark[1]
# current_celltype <- top_hits_ranges$celltype[1]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[2]
# current_mark <-         top_hits_ranges$mark[2]
# current_celltype <- top_hits_ranges$celltype[2]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[3]
# current_mark <-         top_hits_ranges$mark[3]
# current_celltype <- top_hits_ranges$celltype[3]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[4]
# current_mark <-         top_hits_ranges$mark[4]
# current_celltype <- top_hits_ranges$celltype[4]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[5]
# current_mark <-         top_hits_ranges$mark[5]
# current_celltype <- top_hits_ranges$celltype[5]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ##################################################################################################
# ##################################################################################################
# ##################################################################################################
# current_rsID <-           top_hits_ranges$rs[6]
# current_mark <-         top_hits_ranges$mark[6]
# current_celltype <- top_hits_ranges$celltype[6]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[7]
# current_mark <-         top_hits_ranges$mark[7]
# current_celltype <- top_hits_ranges$celltype[7]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[8]
# current_mark <-         top_hits_ranges$mark[8]
# current_celltype <- top_hits_ranges$celltype[8]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[9]
# current_mark <-         top_hits_ranges$mark[9]
# current_celltype <- top_hits_ranges$celltype[9]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[10]
# current_mark <-         top_hits_ranges$mark[10]
# current_celltype <- top_hits_ranges$celltype[10]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[11]
# current_mark <-         top_hits_ranges$mark[11]
# current_celltype <- top_hits_ranges$celltype[11]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[12]
# current_mark <-         top_hits_ranges$mark[12]
# current_celltype <- top_hits_ranges$celltype[12]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# current_rsID <-           top_hits_ranges$rs[13]
# current_mark <-         top_hits_ranges$mark[13]
# current_celltype <- top_hits_ranges$celltype[13]
# 
# ###################################################################################################
# # Plot cohort
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)
# 
# dev.off()
# ###################################################################################################
# # Plot cohort + cell lines
# png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
#     width = 14,
#     height = 10,
#     units = "in",
#     res = 300)
# 
# plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)
# 
# dev.off()
###################################################################################################
###################################################################################################
###################################################################################################
current_rsID <-           top_hits_ranges$rs[14]
current_mark <-         top_hits_ranges$mark[14]
current_celltype <- top_hits_ranges$celltype[14]

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)

dev.off()
###################################################################################################
# Plot cohort + cell lines
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)

dev.off()
###################################################################################################
###################################################################################################
###################################################################################################
current_rsID <-           top_hits_ranges$rs[15]
current_mark <-         top_hits_ranges$mark[15]
current_celltype <- top_hits_ranges$celltype[15]

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)

dev.off()
###################################################################################################
# Plot cohort + cell lines
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)

dev.off()
###################################################################################################
###################################################################################################
###################################################################################################
current_rsID <-           top_hits_ranges$rs[16]
current_mark <-         top_hits_ranges$mark[16]
current_celltype <- top_hits_ranges$celltype[16]

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)

dev.off()
###################################################################################################
# Plot cohort + cell lines
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)

dev.off()
###################################################################################################
###################################################################################################
###################################################################################################
current_rsID <-           top_hits_ranges$rs[17]
current_mark <-         top_hits_ranges$mark[17]
current_celltype <- top_hits_ranges$celltype[17]

###################################################################################################
# Plot cohort
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort(current_rsID, str_c(current_mark, "_", current_celltype), -2000, +2000)

dev.off()
###################################################################################################
# Plot cohort + cell lines
png(str_c("plots_gviz/", current_rsID, "_", current_mark, "_", current_celltype, "_cohort_incl_cl.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plot_cohort_incl_cl(current_rsID, str_c(current_mark, "_", current_celltype), current_mark, -2000, +2000)

dev.off()
###################################################################################################


