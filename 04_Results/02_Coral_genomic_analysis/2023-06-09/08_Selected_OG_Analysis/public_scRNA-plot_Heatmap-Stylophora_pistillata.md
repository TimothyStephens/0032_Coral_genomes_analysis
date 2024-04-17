---
title: "Complex heatmap of Stylophora pistillata scRNA results"
author: "Timothy Stephens"
date: "21 November, 2023"
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

out.prefix <- "public_scRNA-plot_Heatmap-Stylophora_pistillata"
```



Load data inot `R`

```r
data <- read.table("public_scRNA/data/Stylophora_pistillata_GAJOv1.COMBINED_genes.cell_type_gene_FC.tsv", header=T, row.names=1, sep='\t', check.names=FALSE)
color <- read.table("public_scRNA/data/Stylophora_pistillata_GAJOv1.celltype_color.tsv", header=T, sep='\t', comment.char='')
OG.long <- read.table("../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz", header=T, sep='\t')
OG.class <- read.table("../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz", header=T, sep='\t')
```


Cleanup and melt data.

```r
data.formatted <- subset(data, select=-c(old_name)) %>% replace(is.na(.), 0) %>% replace(. < 2, 0)

#data.formatted[rowSums(is.na(data.formatted)) != ncol(data.formatted), ]
```


Merge OG datasets into a single dataframe, then filter so that we only keep the info of the genes in data

```r
OG <- merge(OG.long, OG.class %>% select(orthogroup_id, designation)) %>% filter(sequence_id %in% rownames(data.formatted))
```


Cleanup cell color metadata dataframe.

```r
# Reformat metacell info (reduce to just named cell types, not metacell numbers)
color.celltype <- subset(color, select=-c(metacell, metacell_color)) %>% unique()
rownames(color.celltype) <- paste(color.celltype$tissue,"_",color.celltype$cell, sep='')

# Get just metadata for metacell types in dataset
color.celltype <- color.celltype[colnames(data.formatted) %>% unique(),]
```




```r
# Add row names
color.OG <- OG %>% select(sequence_id, designation) %>% unique()
rownames(color.OG) <- color.OG$sequence_id

# Get just the OGs in dataset
color.OG <- color.OG[rownames(data.formatted) %>% unique(),]
```




```r
plot_ComplexHeatmap <- function(df, 
                                color.columns, color.rows, 
                                row.split.ids, 
                                cluster_columns=FALSE,
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
  t <- color.columns[, c("cell", "cell_color")] %>% unique()
  cell <- t$cell_color
  names(cell) <- t$cell
  
  # Get broadcell type colors
  t <- color.columns[, c("broadcell", "broadcell_color")] %>% unique()
  broadcell <- t$broadcell_color
  names(broadcell) <- t$broadcell
  
  # Load colors into data.object
  col <- list(
    developmental_stage=developmental_stage,
    cell=cell,
    broadcell=broadcell
  )
  col.colors <- HeatmapAnnotation(
    developmental_stage=color.columns$tissue,
    cell=color.columns$cell,
    broadcell=color.columns$broadcell,
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

  # Load colors into data.object
  col <- list(
    designation=designation
  )
  row.colors <- rowAnnotation(
    designation=color.rows$designation,
    col=col
  )
  
  
  
  ##
  ## Get factors to split rows by 
  ##
  row.split <- data.frame(color.rows[, row.split.ids])
  rownames(row.split) <- NULL
  for (id in row.split.ids){
    row.split[id] <- factor(row.split[[id]], levels=names(get(id)))
  }
  
  
  
  ##
  ## Plot and save heatmap
  ##
  pdf(out.plot.file, width=plot.width, height=plot.heigth)
  p <- Heatmap(df, 
          name = "Gene FC in scRNA-seq data", #title of legend
          column_title = "Metacell types", row_title = "Orthogroups",
          row_names_gp = gpar(fontsize = 1), # Text size for row names
          column_names_gp = gpar(fontsize = 3), # Text size for row names
          cluster_columns=cluster_columns,
          top_annotation = col.colors,
          left_annotation = row.colors,
          row_split=row.split,
          cluster_row_slices=FALSE,
          use_raster=FALSE,
          col = circlize::colorRamp2(c(0, 1, 2), c("white", "#f781bf", "red")),
          )
  return(p)
  draw(p)
  dev.off()
}
```




```r
color.columns <- color.celltype
color.rows <- color.OG
df <- data.formatted

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), FALSE, paste(out.prefix, ".plot_Heatmap.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

```r
plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), TRUE,  paste(out.prefix, ".plot_Heatmap-columnClust.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Light-Shared"), ]
df <- data.formatted[rownames(data.formatted) %in% color.rows$sequence_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), FALSE, paste(out.prefix, ".plot_Heatmap_LightShared.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

```r
plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), TRUE,  paste(out.prefix, ".plot_Heatmap_LightShared-columnClust.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Light-Restricted"), ]
df <- data.formatted[rownames(data.formatted) %in% color.rows$sequence_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), FALSE, paste(out.prefix, ".plot_Heatmap_LightRestricted.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

```r
plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), TRUE,  paste(out.prefix, ".plot_Heatmap_LightRestricted-columnClust.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Dark-Shared"), ]
df <- data.formatted[rownames(data.formatted) %in% color.rows$sequence_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), FALSE, paste(out.prefix, ".plot_Heatmap_DarkShared.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

```r
plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), TRUE,  paste(out.prefix, ".plot_Heatmap_DarkShared-columnClust.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```




```r
color.columns <- color.celltype
color.rows <- color.OG[color.OG$designation %in% c("Dark-Restricted"), ]
df <- data.formatted[rownames(data.formatted) %in% color.rows$sequence_id, ]

plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), FALSE, paste(out.prefix, ".plot_Heatmap_DarkRestricted.pdf", sep=''))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

```r
plot_ComplexHeatmap(df, color.columns, color.rows, c("designation"), TRUE,  paste(out.prefix, ".plot_Heatmap_DarkRestricted-columnClust.pdf", sep=''))
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
##  [1] sass_0.4.7          utf8_1.2.4          generics_0.1.3     
##  [4] shape_1.4.6         stringi_1.8.1       digest_0.6.33      
##  [7] magrittr_2.0.3      evaluate_0.23       iterators_1.0.14   
## [10] circlize_0.4.15     fastmap_1.1.1       foreach_1.5.2      
## [13] doParallel_1.0.17   plyr_1.8.9          jsonlite_1.8.7     
## [16] GlobalOptions_0.1.2 fansi_1.0.5         codetools_0.2-19   
## [19] jquerylib_0.1.4     cli_3.6.1           rlang_1.1.2        
## [22] crayon_1.5.2        withr_2.5.2         cachem_1.0.8       
## [25] yaml_2.3.7          tools_4.3.1         parallel_4.3.1     
## [28] colorspace_2.1-0    GetoptLong_1.0.5    BiocGenerics_0.46.0
## [31] vctrs_0.6.4         R6_2.5.1            png_0.1-8          
## [34] magick_2.8.1        matrixStats_1.1.0   stats4_4.3.1       
## [37] lifecycle_1.0.4     stringr_1.5.1       S4Vectors_0.38.2   
## [40] IRanges_2.34.1      clue_0.3-65         cluster_2.1.4      
## [43] pkgconfig_2.0.3     pillar_1.9.0        bslib_0.5.1        
## [46] glue_1.6.2          Rcpp_1.0.11         xfun_0.41          
## [49] tidyselect_1.2.0    rstudioapi_0.15.0   knitr_1.45         
## [52] rjson_0.2.21        htmltools_0.5.7     rmarkdown_2.25     
## [55] compiler_4.3.1
```


