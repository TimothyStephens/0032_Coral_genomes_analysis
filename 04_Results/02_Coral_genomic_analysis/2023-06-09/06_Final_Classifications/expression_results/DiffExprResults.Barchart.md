---
title: "Differental Expression Data"
author: "Timothy Stephens"
date: "08 February, 2024"
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
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
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
library(ggh4x)
library(gtools)
options(scipen = 999) #Prevent scientific notation
```





# Functions

```r
plot_expression_results <- function(prefix, title, cols){
  data  <- read.table(paste("DiffExprResults.",prefix,".filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv",sep=''),  
                      sep='\t', header=TRUE) %>%
    mutate(Group = if_else(Conditions=="Total", "Total", "Expression_Results")) %>%
    mutate(Group = factor(Group, levels=c("Total", "Expression_Results"))) %>%
    mutate(Conditions = factor(Conditions, levels=mixedsort(unique(Conditions))))
  
  p <- ggplot(data, aes(fill=Designation, x=Conditions, y=Percent)) + 
    geom_bar(position="dodge", stat="identity") +
    facet_wrap2(vars(Group), scales = "free_x") +
    geom_text(aes(label=Count), 
              hjust = -.2, position = position_dodge(.9),
              angle=90, size=2
              ) +
    scale_fill_manual(values=c("Dark" = "#984ea3", "Light" = "#ff7f00")) +
    scale_y_continuous(expand=c(0,0), limits=c(0, max(data$Percent) * 1.05)) +
    theme_bw() +
    ylab("Percent genes with significant differential expression between two conditions (%)") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    force_panelsizes(
        cols = cols
        ) +
    theme(panel.border=element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(strip.text.x.top = element_text(angle = 0, size = 12, face = "bold", hjust = 0.5)) +
    labs(title=title) +
    theme(plot.title=element_text(
      size  = 16,
      hjust = 0.5
      ))
  return(p)
}
```





# Plots


```r
p1 <- plot_expression_results("Mcapitata.3TP", "M. capitata 3TP", c(1, 3))
p2 <- plot_expression_results("Mcapitata.12TP", "M. capitata 12TP", c(1, 33))
p3 <- plot_expression_results("Pacuta.12TP", "P. acuta 12TP", c(1, 30))
```


```r
plot.list <- list(p1, p2, p3)
plots <- marrangeGrob(plot.list, nrow=1, ncol=1)
plots
```

![](DiffExprResults.Barchart_files/figure-html/plot_total_seqs-1.png)<!-- -->![](DiffExprResults.Barchart_files/figure-html/plot_total_seqs-2.png)<!-- -->![](DiffExprResults.Barchart_files/figure-html/plot_total_seqs-3.png)<!-- -->

```r
ggsave(filename="DiffExprResults.Barchart.pdf", plots, width=21, height=29.7, units="cm")
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] gtools_3.9.5       ggh4x_0.2.8        dplyr_1.1.4        gridExtra_2.3     
## [5] cowplot_1.1.3      RColorBrewer_1.1-3 ggplot2_3.4.4      reshape2_1.4.4    
## [9] tibble_3.2.1      
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.4      jsonlite_1.8.8    highr_0.10        compiler_4.3.1   
##  [5] tidyselect_1.2.0  Rcpp_1.0.12       stringr_1.5.1     jquerylib_0.1.4  
##  [9] textshaping_0.3.7 systemfonts_1.0.5 scales_1.3.0      yaml_2.3.8       
## [13] fastmap_1.1.1     R6_2.5.1          plyr_1.8.9        labeling_0.4.3   
## [17] generics_0.1.3    knitr_1.45        munsell_0.5.0     bslib_0.6.1      
## [21] pillar_1.9.0      rlang_1.1.3       utf8_1.2.4        cachem_1.0.8     
## [25] stringi_1.8.3     xfun_0.41         sass_0.4.8        cli_3.6.2        
## [29] withr_3.0.0       magrittr_2.0.3    digest_0.6.34     grid_4.3.1       
## [33] rstudioapi_0.15.0 lifecycle_1.0.4   vctrs_0.6.5       evaluate_0.23    
## [37] glue_1.7.0        farver_2.1.1      ragg_1.2.7        fansi_1.0.6      
## [41] colorspace_2.1-0  rmarkdown_2.25    tools_4.3.1       pkgconfig_2.0.3  
## [45] htmltools_0.5.7
```
