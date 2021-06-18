library(GenomicFeatures)
library(Guitar)
txdb = makeTxDbFromGFF('Mus_musculus_MM10_forRNAseq3875_20170608.gtf')

stBedFiles = list("A.merged.bed","B.merged.bed")

pdf("metaplot.tmp.pdf",width=8,height=6)
GuitarPlot(txTxdb = txdb, 
           stBedFiles=stBedFiles,
           headOrtail = FALSE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = c("A","B") 
)+theme(panel.grid.major.y  = element_blank(),
        axis.text.y = element_text(color = 'black',size=10),
        axis.title.y = element_text(size=13),
        legend.position = 'right',
        panel.background = element_blank() 
)+scale_fill_manual(values = c('transparent','transparent'))+
  geom_hline(yintercept =0)+scale_x_continuous(expand = c(0,0))


dev.off()

# R version 3.6.0 (2019-04-26)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

# Matrix products: default

# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936    LC_MONETARY=Chinese (Simplified)_China.936
# [4] LC_NUMERIC=C                               LC_TIME=Chinese (Simplified)_China.936    

# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] Guitar_2.0.0           dplyr_1.0.2            knitr_1.30             ggplot2_3.3.3          magrittr_2.0.1         rtracklayer_1.44.4    
 # [7] GenomicFeatures_1.36.4 AnnotationDbi_1.46.1   Biobase_2.44.0         GenomicRanges_1.36.1   GenomeInfoDb_1.20.0    IRanges_2.18.3        
# [13] S4Vectors_0.22.1       BiocGenerics_0.30.0   

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.6                  lattice_0.20-41             prettyunits_1.1.1           Rsamtools_2.0.3             Biostrings_2.52.0          
 # [6] digest_0.6.27               utf8_1.1.4                  R6_2.5.0                    RSQLite_2.2.2               httr_1.4.2                 
# [11] pillar_1.5.1                zlibbioc_1.30.0             rlang_0.4.10                progress_1.2.2              blob_1.2.1                 
# [16] Matrix_1.3-2                labeling_0.4.2              BiocParallel_1.18.1         stringr_1.4.0               RCurl_1.98-1.2             
# [21] bit_4.0.4                   biomaRt_2.40.5              munsell_0.5.0               DelayedArray_0.10.0         compiler_3.6.0             
# [26] xfun_0.20                   pkgconfig_2.0.3             tidyselect_1.1.0            SummarizedExperiment_1.14.1 tibble_3.1.0               
# [31] GenomeInfoDbData_1.2.1      matrixStats_0.57.0          XML_3.99-0.3                fansi_0.4.2                 crayon_1.4.1               
# [36] withr_2.4.1                 GenomicAlignments_1.20.1    bitops_1.0-6                grid_3.6.0                  gtable_0.3.0               
# [41] lifecycle_1.0.0             DBI_1.1.0                   scales_1.1.1                stringi_1.5.3               farver_2.1.0               
# [46] XVector_0.24.0              ellipsis_0.3.1              vctrs_0.3.6                 generics_0.1.0              tools_3.6.0                
# [51] bit64_4.0.5                 glue_1.4.2                  purrr_0.3.4                 hms_0.5.3                   colorspace_2.0-0           
# [56] memoise_1.1.0