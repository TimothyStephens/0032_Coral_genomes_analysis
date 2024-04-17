---
title: "Complex heatmap of scRNA results"
author: "Timothy Stephens"
date: "20 October, 2023"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 8, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(tibble)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
```

```
## Loading required package: grid
```

```
## ========================================
## ComplexHeatmap version 2.16.0
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite either one:
## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##     genomic data. Bioinformatics 2016.
## 
## 
## The new InteractiveComplexHeatmap package can directly export static 
## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(ComplexHeatmap))
## ========================================
```

```r
options(scipen = 999) #Prevent scientific notation
```




```r
data <- read.table("scRNA.Orthogroups.cell_type_gene_FC1.5_counts.tsv", header=T, row.names=1, sep='\t', check.names=FALSE)
color.celltype <- read.table("../data/Stylophora_pistillata_GAJOv1.celltype_color.tsv", header=T, sep='\t', comment.char='')
color.OG <- read.table("scRNA.Orthogroups.classification.txt", header=T, sep='\t')
```




```r
# Reformat metacell info (reduce to just named cell types, not metacell numbers)
color.celltype <- subset(color.celltype, select=-c(metacell)) %>% unique()
rownames(color.celltype) <- color.celltype$metacell_type

# Get just metadata for metacell types in dataset
color.celltype <- color.celltype[colnames(data),]
```




```r
# Add row names
rownames(color.OG) <- color.OG$orthogroup_id
# Get just the OGs in dataset
color.OG <- color.OG[rownames(data),]
```




```r
plot_ComplexHeatmap <- function(df, 
                                color.columns, color.rows, 
                                row.split.ids, 
                                out.plot.file, 
                                plot.width=12, plot.heigth=12){
  ##
  ## Column (cell type) colors
  ##
  # Get tissue/developmental tissue type colors
  t <- color.columns[, c("tissue", "tissue_color")] %>% unique()
  developmental_stage <- t$tissue_color
  names(developmental_stage) <- t$tissue
  
  # Get metacell type colors
  t <- color.columns[, c("metacell_type", "metacell_type_color")] %>% unique()
  metacell_type <- t$metacell_type_color
  names(metacell_type) <- t$metacell_type
  
  # Get broadcell yupe colors
  t <- color.columns[, c("broadcell_type", "broadcell_type_color")] %>% unique()
  broadcell_type <- t$broadcell_type_color
  names(broadcell_type) <- t$broadcell_type
  
  # Load colors into data.object
  col <- list(
    developmental_stage=developmental_stage,
    metacell_type=metacell_type,
    broadcell_type=broadcell_type
  )
  col.colors <- HeatmapAnnotation(
    developmental_stage=color.columns$tissue,
    metacell_type=color.columns$metacell_type,
    broadcell_type=color.columns$broadcell_type,
    col=col
  )
  
  
  
  ##
  ## Row (OG) colors
  ##
  # Get OG designation colors
  designation <- c(
    "Dark-Restricted"="#984ea3", 
    "Dark-Shared"="#f781bf", 
    "Light-Restricted"="#ff7f00", 
    "Light-Shared"="#ffff33"
  )
  
  # Get no. sequences/OG colors
  no_sequences = circlize::colorRamp2(c(1, 50, 100), c("blue" , "purple", "red"))
  
  # Get best strat taxa
  best_strata_taxa <- c(
    "NA"="#808080",
    "Anthozoa"="#a6cee3",
    "Hexacorallia"="#1f78b4",
    "Cnidaria"="#b2df8a",
    "Astrocoeniina"="#33a02c",
    "Scleractinia"="#fb9a99",
    "Eumetazoa"="#e31a1c",
    "50429"="#6a3d9a",
    "Opisthokonta"="#fdbf6f",
    "Metazoa"="#ff7f00"
  )
  
  # Load colors into data.object
  col <- list(
    designation=designation,
    no_sequences=no_sequences,
    best_strata_taxa=best_strata_taxa
  )
  row.colors <- rowAnnotation(
    designation=color.rows$designation,
    no_sequences=color.rows$no_sequences,
    best_strata_taxa=color.rows$best_strata_taxa,
    col=col
  )
  
  
  
  ##
  ## Get factors to split rows by 
  ##
  row.split <- data.frame(color.rows[, row.split.ids])
  rownames(row.split) <- NULL
  
  
  
  ##
  ## Plot and save heatmap
  ##
  pdf(out.plot.file, width=plot.width, height=plot.heigth)
  p <- Heatmap(df, 
          name = "No. genes FC > 1.5", #title of legend
          column_title = "Metacell types", row_title = "Orthogroups",
          row_names_gp = gpar(fontsize = 1), # Text size for row names
          column_names_gp = gpar(fontsize = 3), # Text size for row names
          cluster_columns=FALSE,
          top_annotation = col.colors,
          left_annotation = row.colors,
          split=row.split,
          use_raster=FALSE,
          col = circlize::colorRamp2(c(0, 1, 10), c("white", "#f781bf", "red")),
          )
  return(p)
  draw(p)
  dev.off()
}
```




```r
color.columns <- color.celltype
color.rows <- color.OG
df <- data

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), "plot_Heatmap.pdf")
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Light-Shared"), ]
df <- data[rownames(data) %in% color.rows$orthogroup_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation", "best_strata_taxa"), "plot_Heatmap_LightShared.pdf")
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Light-Restricted"), ]
df <- data[rownames(data) %in% color.rows$orthogroup_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation", "best_strata_taxa"), "plot_Heatmap_LightRestricted.pdf")
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Dark-Shared"), ]
df <- data[rownames(data) %in% color.rows$orthogroup_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation", "best_strata_taxa"), "plot_Heatmap_DarkShared.pdf")
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Dark-Restricted"), ]
df <- data[rownames(data) %in% color.rows$orthogroup_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation", "best_strata_taxa"), "plot_Heatmap_DarkRestricted.pdf")
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```



# Session Info


```r
sessionInfo()
```

```
## R version 4.3.1 (2023-06-16)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Ventura 13.2.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] ComplexHeatmap_2.16.0 RColorBrewer_1.1-3    reshape2_1.4.4       
## [4] dplyr_1.1.3           tibble_3.2.1         
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.7          utf8_1.2.3          generics_0.1.3     
##  [4] shape_1.4.6         stringi_1.7.12      digest_0.6.33      
##  [7] magrittr_2.0.3      evaluate_0.22       iterators_1.0.14   
## [10] circlize_0.4.15     fastmap_1.1.1       foreach_1.5.2      
## [13] doParallel_1.0.17   plyr_1.8.9          jsonlite_1.8.7     
## [16] GlobalOptions_0.1.2 fansi_1.0.5         codetools_0.2-19   
## [19] jquerylib_0.1.4     cli_3.6.1           rlang_1.1.1        
## [22] crayon_1.5.2        cachem_1.0.8        yaml_2.3.7         
## [25] tools_4.3.1         parallel_4.3.1      colorspace_2.1-0   
## [28] GetoptLong_1.0.5    BiocGenerics_0.46.0 vctrs_0.6.4        
## [31] R6_2.5.1            png_0.1-8           magick_2.8.0       
## [34] matrixStats_1.0.0   stats4_4.3.1        lifecycle_1.0.3    
## [37] stringr_1.5.0       S4Vectors_0.38.2    IRanges_2.34.1     
## [40] clue_0.3-65         cluster_2.1.4       pkgconfig_2.0.3    
## [43] pillar_1.9.0        bslib_0.5.1         glue_1.6.2         
## [46] Rcpp_1.0.11         xfun_0.40           tidyselect_1.2.0   
## [49] rstudioapi_0.15.0   knitr_1.44          rjson_0.2.21       
## [52] htmltools_0.5.6.1   rmarkdown_2.25      compiler_4.3.1
```


