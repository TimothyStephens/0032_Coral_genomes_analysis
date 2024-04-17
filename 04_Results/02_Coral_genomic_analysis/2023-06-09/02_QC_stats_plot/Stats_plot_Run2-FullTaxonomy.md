---
title: "Plot QC results for each genome/transcriptome used in Orthofinder Run2"
author: "Timothy Stephens"
date: "06/08/2023"
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
library(ggplot2)
library(ggtree)
```

```
## ggtree v3.8.2 For help: https://yulab-smu.top/treedata-book/
## 
## If you use the ggtree package suite in published research, please cite
## the appropriate paper(s):
## 
## Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
## ggtree: an R package for visualization and annotation of phylogenetic
## trees with their covariates and other associated data. Methods in
## Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
## 
## Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
## Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu. Ggtree: A serialized data
## object for visualization of a phylogenetic tree and annotation data.
## iMeta 2022, 1(4):e56. doi:10.1002/imt2.56
## 
## Guangchuang Yu.  Data Integration, Manipulation and Visualization of
## Phylogenetic Trees (1st edition). Chapman and Hall/CRC. 2022,
## doi:10.1201/9781003279242
```

```r
library(aplot)
library(treeio)
```

```
## treeio v1.25.3 For help: https://yulab-smu.top/treedata-book/
## 
## If you use the ggtree package suite in published research, please cite
## the appropriate paper(s):
## 
## LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
## Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
## for phylogenetic tree input and output with richly annotated and
## associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
## doi: 10.1093/molbev/msz240
## 
## Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
## Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu. Ggtree: A serialized data
## object for visualization of a phylogenetic tree and annotation data.
## iMeta 2022, 1(4):e56. doi:10.1002/imt2.56
## 
## S Xu, Z Dai, P Guo, X Fu, S Liu, L Zhou, W Tang, T Feng, M Chen, L
## Zhan, T Wu, E Hu, Y Jiang, X Bo, G Yu. ggtreeExtra: Compact
## visualization of richly annotated phylogenetic data. Molecular Biology
## and Evolution. 2021, 38(9):4039-4042. doi: 10.1093/molbev/msab166
```

```r
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
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:ggtree':
## 
##     rotate
```

```r
options(scipen = 999) #Prevent scientific notation
```





# Functions

Function to plot BUSCO categories - uses the standard BUSCO color scheme and also plots the percent for each category.

```r
format_values_1 = function(values){
  out <- c()
  for(value in values) {
    if(value != 0){
      out <- c(out, sprintf("  %1.1f%%", value))
    } else {
      out <- c(out, "")
    }
  }
  return(out)
}

plot_BUSCO_percentages = function(data, plot_title, 
                                  plot.axis.text.y=FALSE,
                                  plot.legend=FALSE,
                                  plot_size_ratio=1, plot_family="sans",
                                  my_colors=c("#56B4E9", "#3492C7", "#F0E442", "#F04442") # Color pallet for BUSCO gene categories.
                                  )
{
  p <- data %>%
    melt(id.var="SampleID") %>%
    rename(BUSCO_categories = variable, Count = value) %>%
    mutate(BUSCO_categories = factor(BUSCO_categories, c("Complete", "Single-copy", "Duplicated", "Fragmented", "Missing"))) %>%
    filter(BUSCO_categories!="Complete") %>%
    filter(BUSCO_categories!="Total") %>%
    ggplot(aes(y = Count, x = SampleID, fill = BUSCO_categories)) +
      geom_bar(position = position_stack(reverse = TRUE),
               stat="identity") +
      coord_flip() +
      # geom_text_repel: https://ggrepel.slowkow.com/articles/examples.html
      geom_text_repel(aes(label = format_values_1(Count)), # Add two spaces in from of text so it centers off the number not number%
                size = rel(5)*plot_size_ratio,
                fontface = "bold",
                max.overlaps = Inf, #Set max.overlaps = Inf to override this behavior and always show all labels, regardless of too many overlaps.
                position = position_stack(vjust=0.5, reverse = TRUE),
                xlim = c(-Inf, Inf), # Set xlim or ylim to Inf or -Inf to disable repulsion away from the edges of the panel.
                ylim = c(-Inf, Inf),
                point.size = NA, #size of each point for each text label; Set point.size = NA to prevent label repulsion away from data points.
                force=0.5, # force of repulsion between overlapping text labels
                box.padding=0, # padding around the text label
                direction="y") + # move text labels “both” (default), “x”, or “y” directions
      theme_gray(base_size = 8) + 
      scale_y_continuous(labels = c("0","20","40","60","80","100"), 
                         breaks = c(0,20,40,60,80,100)) +
      scale_fill_manual(values = my_colors,
                        labels = c("Complete (C) and single-copy (S)",
                                   "Complete (C) and duplicated (D)",
                                   "Fragmented (F)",
                                   "Missing (M)")) +
      ggtitle(plot_title) +
      xlab("") +
      ylab("%BUSCOs") +
      theme(plot.title = element_text(family = plot_family,
                                      hjust=0.5, 
                                      colour = "black", 
                                      size = rel(2.2)*plot_size_ratio,
                                      face = "bold")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "mm")) +
      theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) +
      theme(panel.grid.minor = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      theme(axis.text.y = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.text.x = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.line = element_line(linewidth = 1*plot_size_ratio, colour = "black")) +
      theme(axis.ticks.length = unit(0.85*plot_size_ratio, "cm")) +
      theme(axis.ticks.y = element_line(colour = "white", linewidth = 0)) +
      theme(axis.ticks.x = element_line(colour = "#222222")) +
      theme(axis.ticks.length = unit(0.4*plot_size_ratio, "cm")) + 
      theme(axis.title.x = element_text(family = plot_family, size = rel(2)*plot_size_ratio)) +
      theme(legend.position = "bottom", legend.title = element_blank()) + 
      theme(legend.text = element_text(family = plot_family, size = rel(1.2)*plot_size_ratio)) +
      theme(legend.key.size = unit(1.5*plot_size_ratio,"line")) +
      guides(fill = guide_legend(override.aes = list(colour = NULL))) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  if (!plot.axis.text.y) {
    p <- p + theme(axis.text.y=element_blank())
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}
```

Function to plot basic single-category bar chart (not stacked). User provided color for bars.

```r
format_values_2 = function(values){
  out <- c()
  for(value in values) {
    if(value != 0){
      out <- c(out, scales::comma(value))
    } else {
      out <- c(out, "")
    }
  }
  return(out)
}

plot_barchart = function(data, column.id, 
                         plot_title, xlab, bar.fill, 
                         divide.by=1, y.log2=FALSE, y.log10=FALSE, 
                         add.geom_text=TRUE, plot.axis.text.y=FALSE,
                         plot_size_ratio=1, plot_family="sans")
{
  p <- data %>%
    melt(id.var="SampleID") %>%
    filter(variable==column.id) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = value/divide.by) %>%
    ggplot(aes(y = value, x = SampleID)) +
      geom_bar(position = position_stack(),
               stat="identity",
               fill=bar.fill) +
      scale_y_continuous(labels = comma) +
      coord_flip() +
      theme_gray(base_size = 8) + 
      ggtitle(plot_title) +
      xlab("") +
      ylab(xlab) +
      theme(plot.title = element_text(family = plot_family,
                                      hjust=0.5, 
                                      colour = "black", 
                                      size = rel(2.2)*plot_size_ratio,
                                      face = "bold")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "mm")) +
      theme(legend.position="top",legend.title = element_blank()) + 
      theme(legend.text = element_text(family = plot_family, size = rel(1.2)*plot_size_ratio)) +
      theme(legend.key.size = unit(1.5*plot_size_ratio,"line")) +
      theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) +
      theme(panel.grid.minor = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      theme(axis.text.y = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.text.x = element_text(family = plot_family, face = "bold", size = rel(1.66)*plot_size_ratio, colour = "black")) +
      theme(axis.line = element_line(colour = "black", linewidth = 1*plot_size_ratio)) +
      theme(axis.ticks.length = unit(0.85*plot_size_ratio, "cm")) +
      theme(axis.ticks.y = element_line(colour = "white", linewidth = 0)) +
      theme(axis.ticks.x = element_line(colour = "#222222")) +
      theme(axis.ticks.length = unit(0.4*plot_size_ratio, "cm")) + 
      theme(axis.title.x = element_text(family = plot_family, size = rel(2)*plot_size_ratio)) +
      guides(fill = guide_legend(override.aes = list(colour = NULL))) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  if (y.log10) {
    p <- p + scale_y_continuous(trans="log10")
  }
  if (y.log2) {
    p <- p + scale_y_continuous(trans="log2")
  }
  if (add.geom_text) {
    p <- p + geom_text(aes(label=format_values_2(value), y=0), size=0.75, hjust='left')
  }
  if (!plot.axis.text.y) {
    p <- p + theme(axis.text.y=element_blank())
  }
  return(p)
}
```

Function to generate a formatted taxa names from samples sheet.

```r
# Format for plotting
# See https://yulab-smu.top/treedata-book/faq.html?q=rename#faq-formatting-label
format_sp_name = function(x){
  n = x[["genera"]]
  
  # Add species if present
  if(x[["species"]] != ""){
    n = paste(n, x[["species"]], sep="~")
  }
  
  # Either print with extra info (not italics) or not
  if(x[["extra"]] != ""){
    e = str_replace_all(x[["extra"]], " ", "~")
    e = str_replace_all(e, ",", "")
    lab = glue(paste("italic(",n,")~",e, sep=''))
  } else {
    lab = glue(paste("italic(",n,")", sep=''))
  }
  return(lab)
}
```





# Load datasets into R

Load the cleaned and processed datasets into R for plotting.

```r
species.tree <- read.iqtree(
  "../01_Orthofinder/Run2.SpeciesTree_rooted_node_labels.tre")
```

```
## Warning in sub("/.*", "", nlabel) %>% as.numeric: NAs introduced by coercion
```

```
## Warning in sub(".*/", "", nlabel) %>% as.numeric: NAs introduced by coercion
```

```r
species.names <- as.phylo(species.tree)$tip.label

samples <- read.csv("../samples.txt", header=TRUE, sep='\t') %>% 
  filter(sample_id %in% species.names)
```
Extract the sample order from dendrogram (using "working" sample names). Also load pretty sample names and replace the labels in the tree with these new names.





# Generate plot

Global setting for plotting

```r
plot_size_ratio <- 0.2
plot.axis.text.y <- FALSE
plot.legend <- FALSE
```

Data.frame with clade positions in tree and colors.

```r
# Dataframe of clades to highlight
clade.colors <- data.frame(
  node  = c(132,        138,         144,         146,       148,        164,            171,            172,          183,            179,                130,        126,        122,          124,                134,       136,          133,             56,          57),
  fill  = c("blue",     "red",       "brown",     "pink",    "yellow",   "purple",       "green",        "orange",     "blue",         "red",             "purple",    "green",    "yellow",     "red",              "blue",    "blue",       "yellow",       "green",        "blue"),
  alpha = c(0.2,        0.4,         0.2,         0.4,       0.2,        0.2,            0.2,            0.2,          0.2,            0.2,                0.5,        0.5,        0.5,          0.5,                0.5,       0.1,          0.2,            0.4,          0.2),
  n     = c("Cnidaria", "Staurozoa", "Scyphozoa", "Cubozoa", "Hydrozoa", "Octocorallia", "Hexacorallia", "Actiniaria", "Scleractinia", "Corallimorpharia", "Placozoa", "Porifera", "Ctenophora", "Choanoflagellata", "Myxozoa", "Medusozoa",  "Endocnidozoa", "Zoantharia", "Scleractinia")
)
```

Format species tree. Color clades and rename taxa into a pretty format (i.e., genera+species italics + extra info not italics).
 - Show tree with full taxonomy to verify that our clades are in the correct position.

```r
samples$name <- apply(samples, 1, format_sp_name )
p <- ggtree(species.tree) %<+% samples + 
  geom_highlight(
    data = clade.colors,
    mapping = aes(
      node=node,
      fill=I(fill),
      alpha=I(alpha),
    ),
    extendto = 1.85,
    to.bottom=TRUE,
  ) +
  geom_cladelab(
    data = clade.colors,
    mapping = aes(
      node=node,
      label=n,
    ),
    align = TRUE,
    offset.text = 0.0,
    hjust = "center",
    fontsize = 2,
    offset = 0.7,
    barsize = 0,
    fontface="bold",
    ) +
  geom_tiplab(aes(label=name),        size=6*plot_size_ratio, linesize=0.2, parse=TRUE, align=TRUE, offset=0.01) +
  geom_tiplab(aes(label=family),      size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.20) +
  geom_tiplab(aes(label=suborder),    size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.27) +
  geom_tiplab(aes(label=order),       size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.35) +
  geom_tiplab(aes(label=subclass),    size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.44) +
  geom_tiplab(aes(label=class),       size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.54) +
  geom_tiplab(aes(label=phylum),      size=6*plot_size_ratio, linesize=0.0, parse=TRUE, align=TRUE, offset=0.62) +
  geom_tippoint(aes(color=data_type), size=6*plot_size_ratio, show.legend=FALSE) +
  hexpand(0.1)

p + geom_text(aes(label=node))
```

![](Stats_plot_Run2-FullTaxonomy_files/figure-html/species_tree_taxonomy-1.png)<!-- -->





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
##  [1] ggpubr_0.6.0   glue_1.6.2     stringr_1.5.0  ggrepel_0.9.4  scales_1.2.1  
##  [6] reshape2_1.4.4 dplyr_1.1.3    tibble_3.2.1   treeio_1.25.3  aplot_0.2.2   
## [11] ggtree_3.8.2   ggplot2_3.4.4 
## 
## loaded via a namespace (and not attached):
##  [1] yulab.utils_0.1.0  sass_0.4.7         utf8_1.2.4         generics_0.1.3    
##  [5] tidyr_1.3.0        rstatix_0.7.2      ggplotify_0.1.2    stringi_1.7.12    
##  [9] lattice_0.22-5     digest_0.6.33      magrittr_2.0.3     evaluate_0.23     
## [13] grid_4.3.1         fastmap_1.1.1      plyr_1.8.9         jsonlite_1.8.7    
## [17] ape_5.7-1          backports_1.4.1    purrr_1.0.2        fansi_1.0.5       
## [21] lazyeval_0.2.2     jquerylib_0.1.4    abind_1.4-5        cli_3.6.1         
## [25] rlang_1.1.1        munsell_0.5.0      tidytree_0.4.5     withr_2.5.2       
## [29] cachem_1.0.8       yaml_2.3.7         tools_4.3.1        parallel_4.3.1    
## [33] ggsignif_0.6.4     memoise_2.0.1      colorspace_2.1-0   broom_1.0.5       
## [37] vctrs_0.6.4        R6_2.5.1           gridGraphics_0.5-1 lifecycle_1.0.3   
## [41] car_3.1-2          fs_1.6.3           ggfun_0.1.3        pkgconfig_2.0.3   
## [45] pillar_1.9.0       bslib_0.5.1        gtable_0.3.4       Rcpp_1.0.11       
## [49] highr_0.10         xfun_0.41          tidyselect_1.2.0   rstudioapi_0.15.0 
## [53] knitr_1.45         farver_2.1.1       htmltools_0.5.7    nlme_3.1-163      
## [57] patchwork_1.1.3    labeling_0.4.3     carData_3.0-5      rmarkdown_2.25    
## [61] compiler_4.3.1
```
