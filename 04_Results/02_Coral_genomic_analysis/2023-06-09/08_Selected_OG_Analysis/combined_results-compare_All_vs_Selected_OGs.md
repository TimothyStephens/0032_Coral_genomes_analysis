---
title: "Plot All vs. Selected Dark OGs"
author: "Timothy Stephens"
date: "01 April, 2024"
output: 
  html_document:
    keep_md: yes
---

# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 12, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(ggridges)
library(ggplot2)
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
library(scales)
library(ggrepel)
library(stringr)
library(glue)
library(ggpubr)
options(scipen = 999) #Prevent scientific notation

text.scale <- 0.5
```

# Load datasets into R

Load the cleaned and processed datasets into R for plotting.


```r
data.selected.OGs <- read.table("combined_results/combined.stats.OGs.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data.selected.OGs <- data.selected.OGs$orthogroup_id

data.OGs  <- read.table("../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>%
  mutate(selected=if_else(orthogroup_id %in% data.selected.OGs, "Yes", "No", NA))
data.seqs <- read.table("combined_results/combined.stats.ALL_seqs.info.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>%
  mutate(selected=if_else(orthogroup_id %in% data.selected.OGs, "Yes", "No", NA))
```

# Generate plot

Plot number of CDS per proteins (for proteins with genome coords available) - ALL Designations.


```r
col <- "num_CDS"
col.name <- "Number of CDS/exons per protein"
x.min <- 1
x.max <- 20

p1 <- merge(data.seqs, 
      data.OGs %>% select(orthogroup_id, designation), 
      all=TRUE,
      by="orthogroup_id"
      ) %>%
  select(!!sym(col), designation) %>%
  rename(value = !!sym(col)) %>%
  rename(group = designation) %>%
  filter(! is.na(value)) %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         ) %>%
  ggplot(aes(x=value, y=group, fill=group)) + 
    geom_density_ridges(scale = 0.95) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, angle = 90, vjust = 0.5, hjust=0, face="bold")
          ) +
    scale_fill_cyclical(values = c("#984ea3", "#bebada", "#ff7f00", "#fdbf6f")) +
    labs(title="Distribution of CDS/exon count per gene",
         x=col.name, 
         y="Orthogorup Designation") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
p1
```

```
## Picking joint bandwidth of 0.204
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_ALL_num_CDS-1.png)<!-- -->

Plot number of CDS per proteins (for proteins with genome coords available) - ALL Designations.


```r
col <- "PEP_length"
col.name <- "Protein length (aa)"
x.min <- 1
x.max <- 1000

p2 <- merge(data.seqs, 
      data.OGs %>% select(orthogroup_id, designation), 
      all=TRUE,
      by="orthogroup_id"
      ) %>%
  select(!!sym(col), designation) %>%
  rename(value = !!sym(col)) %>%
  rename(group = designation) %>%
  filter(! is.na(value)) %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         ) %>%
  ggplot(aes(x=value, y=group, fill=group)) + 
    geom_density_ridges(scale = 0.95) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, angle = 90, vjust = 0.5, hjust=0, face="bold")
          ) +
    scale_fill_cyclical(values = c("#984ea3", "#bebada", "#ff7f00", "#fdbf6f")) +
    labs(title="Distribution of protein lengths",
         x=col.name, 
         y="Orthogorup Designation") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
p2
```

```
## Picking joint bandwidth of 11.5
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_ALL_PEP_length-1.png)<!-- -->

Plot length of proteins.


```r
col <- "num_CDS"
col.name <- "Number of CDS/exons per protein"
x.min <- 1
x.max <- 20
y.labels <- c("Selected Dark Orthogroups", "Orther Orthogroups")
names(y.labels) <- c("Yes", "No")

p3 <- data.seqs %>% 
  select(!!sym(col), selected) %>%
  rename(value = !!sym(col)) %>%
  rename(group = selected) %>%
  filter(! is.na(value)) %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         ) %>%
  ggplot(aes(x=value, y=group, fill=group)) + 
    geom_density_ridges(scale = 0.95) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, angle = 90, vjust = 0.5, hjust=0, face="bold")
          ) +
    labs(title="Distribution of CDS/exon count per gene") +
    scale_fill_cyclical(values = c("#ff7f00", "#984ea3")) +
    labs(x=col.name) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), label=y.labels)
p3
```

```
## Picking joint bandwidth of 0.195
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_num_CDS-1.png)<!-- -->

Plot length of proteins.


```r
col <- "PEP_length"
col.name <- "Protein length (aa)"
x.min <- 1
x.max <- 1000

p4 <- data.seqs %>% 
  select(!!sym(col), selected) %>%
  rename(value = !!sym(col)) %>%
  rename(group = selected) %>%
  filter(! is.na(value)) %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         ) %>%
  ggplot(aes(x=value, y=group, fill=group)) + 
    geom_density_ridges(scale = 0.95) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, angle = 90, vjust = 0.5, hjust=0, face="bold")
          ) +
    labs(title="Distribution of protein lengths") +
    scale_fill_cyclical(values = c("#ff7f00", "#984ea3")) +
    labs(x=col.name) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), label=y.labels)
p4
```

```
## Picking joint bandwidth of 13.8
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_PEP_length-1.png)<!-- -->

Plot abSENSE results.


```r
selected.labs <- c("Selected Dark Orthogroups", "Orther Orthogroups")
names(selected.labs) <- c("Yes", "No")

p5 <- data.OGs %>%
  select(`non-strata_taxa_expectedHDF_false_proportion`, `non-strata_taxa_expectedHDF_true_proportion`, designation, selected) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_false_proportion`)) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_true_proportion`)) %>%
  rename(x = `non-strata_taxa_expectedHDF_true_proportion`) %>%
  rename(y = `non-strata_taxa_expectedHDF_false_proportion`) %>%
  ggplot(aes(x=x, y=y)) + 
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_fixed()+
    labs(title="abSENSE Homology Detection Failure Assessment",
         x="Percent (%) proteins with evidence of HDF", 
         y="Percent (%) proteins without evidence of HDF") +
    theme_minimal() +
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=8*text.scale),
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, face="bold"),
          legend.key.size   = unit(0.5*text.scale, 'cm'), #change legend key size
          legend.key.height = unit(0.5*text.scale, 'cm'), #change legend key height
          legend.key.width  = unit(0.5*text.scale, 'cm'), #change legend key width
          legend.title = element_text(size=12*text.scale, vjust = 1), #change legend title font size
          legend.text  = element_text(size=8*text.scale, angle = 45, hjust = 1),
          legend.position = "bottom", 
          legend.box = "horizontal"
          ) +
    facet_wrap(vars(selected), labeller = labeller(selected = selected.labs))
p5
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_abSENSE_stat_density_2d-1.png)<!-- -->

Plot abSENSE results.


```r
selected.labs <- c("Selected Dark Orthogroups", "Orther Orthogroups")
names(selected.labs) <- c("Yes", "No")

p6 <- data.OGs %>%
  select(`non-strata_taxa_expectedHDF_false_proportion`, `non-strata_taxa_expectedHDF_true_proportion`, designation, selected) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_false_proportion`)) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_true_proportion`)) %>%
  rename(x = `non-strata_taxa_expectedHDF_true_proportion`) %>%
  rename(y = `non-strata_taxa_expectedHDF_false_proportion`) %>%
  ggplot(aes(x=x, y=y)) + 
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    geom_jitter(color="pink", size=0.1, width=0.5, height=0.5) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_fixed()+
    labs(title="abSENSE Homology Detection Failure Assessment",
         x="Percent (%) proteins with evidence of HDF", 
         y="Percent (%) proteins without evidence of HDF",
          ) +
    theme_minimal() +
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=8*text.scale),
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, face="bold"),
          legend.key.size   = unit(0.5*text.scale, 'cm'), #change legend key size
          legend.key.height = unit(0.5*text.scale, 'cm'), #change legend key height
          legend.key.width  = unit(0.5*text.scale, 'cm'), #change legend key width
          legend.title = element_text(size=12*text.scale, vjust = 1), #change legend title font size
          legend.text  = element_text(size=8*text.scale, angle = 45, hjust = 1),
          legend.position = "bottom", 
          legend.box = "horizontal"
          ) +
    facet_wrap(vars(selected), labeller = labeller(selected = selected.labs))
p6
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_abSENSE_stat_density_2d_Points-1.png)<!-- -->


Plot abSENSE results for just selected dark orthogroups.

```r
p7 <- data.OGs %>%
  select(`non-strata_taxa_expectedHDF_false_proportion`, `non-strata_taxa_expectedHDF_true_proportion`, designation, selected) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_false_proportion`)) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_true_proportion`)) %>%
  rename(x = `non-strata_taxa_expectedHDF_true_proportion`) %>%
  rename(y = `non-strata_taxa_expectedHDF_false_proportion`) %>%
  filter(selected == "Yes") %>%
  ggplot(aes(x=x, y=y)) + 
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_fixed()+
    labs(title="abSENSE Homology Detection Failure Assessment",
         x="Percent (%) proteins with evidence of HDF", 
         y="Percent (%) proteins without evidence of HDF") +
    theme_minimal() +
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=8*text.scale),
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, face="bold"),
          legend.key.size   = unit(0.5*text.scale, 'cm'), #change legend key size
          legend.key.height = unit(0.5*text.scale, 'cm'), #change legend key height
          legend.key.width  = unit(0.5*text.scale, 'cm'), #change legend key width
          legend.title = element_text(size=12*text.scale, vjust = 1), #change legend title font size
          legend.text  = element_text(size=8*text.scale, angle = 45, hjust = 1),
          legend.position = "bottom", 
          legend.box = "horizontal"
          )
p7
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_selected_OGs_abSENSE_stat_density_2d-1.png)<!-- -->


```r
p8 <- data.OGs %>%
  select(`non-strata_taxa_expectedHDF_false_proportion`, `non-strata_taxa_expectedHDF_true_proportion`, designation, selected) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_false_proportion`)) %>%
  filter(! is.na(`non-strata_taxa_expectedHDF_true_proportion`)) %>%
  rename(x = `non-strata_taxa_expectedHDF_true_proportion`) %>%
  rename(y = `non-strata_taxa_expectedHDF_false_proportion`) %>%
  filter(selected == "Yes") %>%
  ggplot(aes(x=x, y=y)) + 
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    geom_jitter(color="pink", size=0.1, width=0.5, height=0.5) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_fixed()+
    labs(title="abSENSE Homology Detection Failure Assessment",
         x="Percent (%) proteins with evidence of HDF", 
         y="Percent (%) proteins without evidence of HDF",
          ) +
    theme_minimal() +
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=8*text.scale),
          plot.title   = element_text(size=14*text.scale, hjust = 0.5, face="bold"),
          axis.title.x = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.title.y = element_text(size=10*text.scale, hjust = 0.5, face="bold"),
          axis.text.x  = element_text(size=8*text.scale),
          axis.text.y  = element_text(size=8*text.scale, face="bold"),
          legend.key.size   = unit(0.5*text.scale, 'cm'), #change legend key size
          legend.key.height = unit(0.5*text.scale, 'cm'), #change legend key height
          legend.key.width  = unit(0.5*text.scale, 'cm'), #change legend key width
          legend.title = element_text(size=12*text.scale, vjust = 1), #change legend title font size
          legend.text  = element_text(size=8*text.scale, angle = 45, hjust = 1),
          legend.position = "bottom", 
          legend.box = "horizontal"
          )
p8
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/plot_selected_OGs_abSENSE_stat_density_2d_Points-1.png)<!-- -->



```r
multi_plot <- ggarrange(p1,p2,p3,p4,p7,p8, #plots that are going to be included in this multipanel figure
                        labels = c("A", "B", "C", "D", "E", "F"), #labels given each panel
                        font.label = list(size = 10, color = "black", face = "bold", family = NULL),
                        ncol = 2, nrow = 3, #adjust plot space 
                        common.legend = F) #does the plot have a common legend
```

```
## Picking joint bandwidth of 0.204
```

```
## Picking joint bandwidth of 11.5
```

```
## Picking joint bandwidth of 0.195
```

```
## Picking joint bandwidth of 13.8
```

```r
multi_plot
```

![](combined_results-compare_All_vs_Selected_OGs_files/figure-html/multi_plot-1.png)<!-- -->





```r
col <- "PEP_length"
col.name <- "Protein length (aa)"
x.min <- 1
x.max <- 1000

d <- data.seqs %>% 
  select(!!sym(col), selected) %>%
  rename(value = !!sym(col)) %>%
  rename(group = selected) %>%
  filter(! is.na(value))  %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         )
ggboxplot(data = d, x = "group", y = "value",
          width = 0.5, size = 0.8, 
          xlab = "", ylab = "value") +
  stat_compare_means(comparisons = list(c("Yes", "No")))





shapiro.test(data.seqs$PEP_length)





col <- "num_CDS"
col.name <- "Number of CDS/exons per protein"
x.min <- 1
x.max <- 20

d <- data.seqs %>% 
  select(!!sym(col), selected) %>%
  rename(value = !!sym(col)) %>%
  rename(group = selected) %>%
  filter(! is.na(value))  %>%
  mutate(value = if_else(value < x.min, x.min, 
                      if_else(value > x.max, x.max, value)
                      )
         )

ggboxplot(data = d, x = "group", y = "value",
          width = 0.5, size = 0.8, 
          xlab = "", ylab = "value") +
  stat_compare_means(comparisons = list(c("Yes", "No")))
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
##  [1] ggpubr_0.6.0   glue_1.7.0     stringr_1.5.1  ggrepel_0.9.5  scales_1.3.0  
##  [6] reshape2_1.4.4 dplyr_1.1.4    tibble_3.2.1   ggplot2_3.5.0  ggridges_0.5.6
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.9         utf8_1.2.4         generics_0.1.3     tidyr_1.3.1       
##  [5] rstatix_0.7.2      stringi_1.8.3      digest_0.6.35      magrittr_2.0.3    
##  [9] RColorBrewer_1.1-3 evaluate_0.23      grid_4.3.1         fastmap_1.1.1     
## [13] plyr_1.8.9         jsonlite_1.8.8     backports_1.4.1    purrr_1.0.2       
## [17] fansi_1.0.6        jquerylib_0.1.4    abind_1.4-5        cli_3.6.2         
## [21] rlang_1.1.3        cowplot_1.1.3      munsell_0.5.0      withr_3.0.0       
## [25] cachem_1.0.8       yaml_2.3.8         tools_4.3.1        ggsignif_0.6.4    
## [29] colorspace_2.1-0   broom_1.0.5        vctrs_0.6.5        R6_2.5.1          
## [33] lifecycle_1.0.4    car_3.1-2          MASS_7.3-60.0.1    pkgconfig_2.0.3   
## [37] pillar_1.9.0       bslib_0.6.2        gtable_0.3.4       Rcpp_1.0.12       
## [41] xfun_0.43          tidyselect_1.2.1   highr_0.10         rstudioapi_0.16.0 
## [45] knitr_1.45         farver_2.1.1       htmltools_0.5.8    rmarkdown_2.26    
## [49] carData_3.0-5      labeling_0.4.3     compiler_4.3.1
```
