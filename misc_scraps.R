Main visualization of top hits - including cell lines
```{r}
png("plots/rs10502002.png",
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plotTracks(list(IdeogramTrack(genome = gviz_genome,
                              chromosome = top_hits_ranges %>% plyranges::filter(rs == "rs10502002") %>% reduce_ranges() %>% seqnames() %>% as.character(),
                              bands = gviz_cyto_bands),
                GenomeAxisTrack(), 
                
                ##### Cis-hBTL markings
                AnnotationTrack(reduce_ranges(top_hits_ranges), 
                                name = "cis-hBTLs",
                                genome = gviz_genome),
                
                ### Set 1
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[1]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "14-17",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[1]],
                                name = "H3K4me3_LC_Individual_1_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 2
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[2]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "14-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[2]],
                                name = "H3K4me3_LC_Individual_2_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 3
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[3]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "11-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[3]],
                                name = "H3K4me3_LC_Individual_3_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 4
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[4]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "22-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[4]],
                                name = "H3K4me3_LC_Individual_4_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 5
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[5]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "24-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[5]],
                                name = "H3K4me3_LC_Individual_5_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 6
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[6]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "38-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[6]],
                                name = "H3K4me3_LC_Individual_6_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 7
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[7]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "30-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[7]],
                                name = "H3K4me3_LC_Individual_7_Peak_Calls",
                                genome = gviz_genome),
                
                ### Set 8
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[8]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "15-18",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[8]],
                                name = "H3K4me3_LC_Individual_8_Peak_Calls",
                                genome = gviz_genome),
                
                ## Set MCF10A
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/MCF10A/peak_calls/MCF10A.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "H3K4me3_MCF10A_Pileup",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/MCF10A/peak_calls/MCF10A.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "H3K4me3_MCF10A_Peak_Calls",
                                genome = gviz_genome),
                
                ## Set hTERT-L9
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L9/peak_calls/HTERT_L9.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "H3K4me3_HTERT_L9_Pileup",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L9/peak_calls/HTERT_L9.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "H3K4me3_HTERT_L9_Peak_Calls",
                                genome = gviz_genome),
                
                ## Set hTERT-L2
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L2/peak_calls/HTERT_L2.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "H3K4me3_HTERT_L2_Pileup",
                          baseline = 0,
                          ylim = c(0, 20)),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L2/peak_calls/HTERT_L2.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "H3K4me3_HTERT_L2_Peak_Calls",
                                genome = gviz_genome),
                
                ##### Gene regions
                BiomartGeneRegionTrack(biomart = gviz_biomart_hg38_genes,
                                       genome = gviz_genome,
                                       symbol = "MMP7",
                                       name = "ENSEMBL Gene",
                                       filter = list(with_refseq_mrna = TRUE))),
           
           # Global settings
           chromosome = top_hits_ranges %>% plyranges::filter(rs == "rs10502002") %>% reduce_ranges() %>% seqnames() %>% as.character(),
           from =       top_hits_ranges %>% plyranges::filter(rs == "rs10502002") %>% reduce_ranges() %>% start() %>% `+`(-1e3), 
           to =         top_hits_ranges %>% plyranges::filter(rs == "rs10502002") %>% reduce_ranges() %>% start() %>% `+`(+1e3),
           transcriptAnnotation = "symbol", 
           stacking = "squish")

dev.off()
```

Set 2 - rs75071948 - ANXA1
```{r}
png("plots/rs75071948_1.png",
    width = 14,
    height = 10,
    units = "in",
    res = 300)

plotTracks(list(IdeogramTrack(genome = gviz_genome,
                              chromosome = top_hits_ranges %>% plyranges::filter(rs == "rs75071948") %>% reduce_ranges() %>% seqnames() %>% as.character(),
                              bands = gviz_cyto_bands),
                GenomeAxisTrack(), 
                
                ##### Cis-hBTL markings
                AnnotationTrack(reduce_ranges(top_hits_ranges), 
                                name = "cis-hBTLs",
                                genome = gviz_genome),
                
                ### Set 1
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[1]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "14-17",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[1]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 2
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[2]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "14-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[2]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 3
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[3]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "11-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[3]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 4
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[4]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "22-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[4]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 5
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[5]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "24-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[5]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 6
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[6]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "38-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[6]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 7
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[7]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "30-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[7]],
                                name = "",
                                genome = gviz_genome),
                
                ### Set 8
                DataTrack(range = list_treat_pileup_bw_wasp_files$H3K4me3_LC[[8]],
                          genome = gviz_genome,
                          type = "mountain",
                          name = "15-18",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(list_peak_wasp$H3K4me3_LC[[8]],
                                name = "",
                                genome = gviz_genome),
                
                ## Set MCF10A
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/MCF10A/peak_calls/MCF10A.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "MCF10A",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/MCF10A/peak_calls/MCF10A.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "",
                                genome = gviz_genome),
                
                ## Set hTERT-L9
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L9/peak_calls/HTERT_L9.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "184-hTERT-L9",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L9/peak_calls/HTERT_L9.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "",
                                genome = gviz_genome),
                
                ## Set hTERT-L2
                DataTrack(range = "/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L2/peak_calls/HTERT_L2.wasp_treat_pileup.bw",
                          genome = gviz_genome,
                          type = "mountain",
                          name = "184-hTERT-L2",
                          baseline = 0,
                          ylim = c(0, 20),
                          span = 0.5),
                AnnotationTrack(readNarrowPeak("/projects/epigenomics3/epigenomics3_results/users/ahauduc/arp/project/modules/cell_lines/1/output2/H3K4me3/HTERT_L2/peak_calls/HTERT_L2.wasp_peaks.narrowPeak.filtered.bed"),
                                name = "",
                                genome = gviz_genome),
                
                ##### Gene regions
                BiomartGeneRegionTrack(biomart = gviz_biomart_hg38_genes,
                                       genome = gviz_genome,
                                       symbol = "ANXA1",
                                       name = "ENSEMBL Gene",
                                       filter = list(with_refseq_mrna = TRUE))),
           
           # Global settings
           chromosome = top_hits_ranges %>% plyranges::filter(rs == "rs75071948") %>% reduce_ranges() %>% seqnames() %>% as.character(),
           from =       top_hits_ranges %>% plyranges::filter(rs == "rs75071948") %>% reduce_ranges() %>% start() %>% `+`(-1e3), 
           to =         top_hits_ranges %>% plyranges::filter(rs == "rs75071948") %>% reduce_ranges() %>% start() %>% `+`(+1e3),
           transcriptAnnotation = "symbol", 
           stacking = "squish")

dev.off()
```

### test
