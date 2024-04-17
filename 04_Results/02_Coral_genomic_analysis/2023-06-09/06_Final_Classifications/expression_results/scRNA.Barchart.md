---
title: "scRNA Cell Type Expression Data"
author: "Timothy Stephens"
date: "12 February, 2024"
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





# Function to plot barchart


```r
load_data <- function(sp, ct){
  data <- read.table(paste("scRNA.",sp,".",ct,".filtered.long.no_genes_per_designation_LvD_and_cellType.prop.colors.tsv", sep=''),  
                     sep='\t', header=TRUE)
  
  u <- unique(data$Conditions)
  u <- u[!u %in% c("Total")]

  if("1" %in% u){
    u <- mixedsort(u)
    u <- c("Total", u)
  } else {
    u <- sort(u)
    u <- c("Total", u)
  }

  data <- data %>%
    mutate(Conditions = factor(Conditions, levels=u))

  u <- unique(data$BroadCellType)
  u <- u[!u %in% c("Total")]
  u <- sort(u)
  u <- c("Total", u)

  data <- data %>%
    mutate(BroadCellType = factor(BroadCellType, levels=u))
  
  return(data)
}
```



```r
plot_bar_FC <- function(df, 
                        axis.text.x.size=3,
                        geom_text.size=0.8,
                        plot.title=""
                        ){
  ##
  ## Get colors to use for facet header/strip colors 
  ##
  df.colors <- unique(df[c("BroadCellType", "BroadCellColor")])
  df.colors <- with(df.colors, df.colors[order(BroadCellType, BroadCellColor),])
  t <- df.colors$BroadCellColor
  names(t) <- df.colors$BroadCellType
  strip <- strip_themed(background_x = elem_list_rect(fill = t))
  
  ##
  ## Get number of columns per panel (use either number of columns in Metacell or Cell)
  ##
  num.columns.per.panel <- merge(
    df %>%
        filter(Type=="Metacell") %>%
        group_by(BroadCellType) %>%
        summarize(distinct_points = n()),
    df %>%
        filter(Type=="Cell") %>%
        group_by(BroadCellType) %>%
        summarize(distinct_points = n()),
    by = "BroadCellType", sort=FALSE
  ) %>% 
    mutate(distinct_points = if_else(distinct_points.x > distinct_points.y, 
                                     distinct_points.x, 
                                     distinct_points.y
                                     )
           )
  
  ##
  ## Plot
  ##
  # Setup bar chart
  p <- df %>%
    arrange(-Percent) %>%
    ggplot(aes(fill=Designation, x=Conditions, y=Percent)) + 
    geom_bar(position=position_identity(), stat="identity") +
    facet_grid2(Type~BroadCellType, 
                independent = "x", 
                scales = "free_x", 
                strip=strip,
                switch = "y"
                ) +
    geom_text(aes(label = Count), 
              hjust    = -0.2, 
              position = position_dodge(1.0),
              angle    = 90,
              size     = geom_text.size
              ) +
    scale_fill_manual(values=c("Dark"  = "#984ea3",
                               "Light" = "#ff7f00"
                               )
                      ) +
    scale_y_continuous(expand=c(0,0), limits=c(0, max(df$Percent) * 1.05)) +
    ylab("Percent genes with significant FC in a given cell type (%)") + 
    labs(title=plot.title) +
    theme_bw() +
    theme(axis.text.x=element_text(
      angle = 90, 
      vjust = 0.5, 
      hjust = 1, 
      size  = axis.text.x.size
      ),
    plot.title=element_text(
      size  = 16,
      hjust = 0.5
      )
    ) +
    theme(panel.border=element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.position = "none") + force_panelsizes(
      rows = seq(max(num.columns.per.panel$distinct_points), length(unique(data$Type))),
      cols = num.columns.per.panel$distinct_points
      ) +
    theme(strip.text.x.top = element_text(angle = 90, size = 4, face = "bold"))

  ## Return plot
  return(p)
}
```





# Plots


```r
sp <- "Stylophora_pistillata_GAJOv1"

data.1 <- load_data(sp, "adult_broad_cell_type_gene_FC") %>% mutate(Type="Broad")
data.2 <- load_data(sp, "adult_cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "adult_metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.1, data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Broad","Cell", "Metacell")))

p.SP.adult <- plot_bar_FC(data, 2.5, 0.8, plot.title=gsub("_"," ", paste(sp," Adult",sep='')))
p.SP.adult
```

![](scRNA.Barchart_files/figure-html/plot_Stylophora_pistillata_GAJOv1_adult-1.png)<!-- -->


```r
sp <- "Stylophora_pistillata_GAJOv1"

data.2 <- load_data(sp, "polyp_cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "polyp_metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Cell", "Metacell")))

p.SP.polyp <- plot_bar_FC(data, 5.0, 0.8, plot.title=gsub("_"," ", paste(sp," Polyp",sep='')))
p.SP.polyp
```

![](scRNA.Barchart_files/figure-html/plot_Stylophora_pistillata_GAJOv1_polyp-1.png)<!-- -->


```r
sp <- "Stylophora_pistillata_GAJOv1"

data.2 <- load_data(sp, "larva_cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "larva_metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Cell", "Metacell")))

p.SP.larva <- plot_bar_FC(data, 5.0, 0.8, plot.title=gsub("_"," ", paste(sp," Larva",sep='')))
p.SP.larva
```

![](scRNA.Barchart_files/figure-html/plot_Stylophora_pistillata_GAJOv1_larva-1.png)<!-- -->


```r
sp <- "Nematostella_vectensis_RRUSv1"

data.1 <- load_data(sp, "broad_cell_type_gene_FC") %>% mutate(Type="Broad")
data.2 <- load_data(sp, "cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.1, data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Broad","Cell", "Metacell")))

p.NV <- plot_bar_FC(data, 3.0, 0.8, plot.title=gsub("_"," ", paste(sp,"",sep='')))
p.NV
```

![](scRNA.Barchart_files/figure-html/plot_Nematostella_vectensis_RRUSv1-1.png)<!-- -->


```r
sp <- "Xenia_sp_CTEAv1"

data.1 <- load_data(sp, "broad_cell_type_gene_FC") %>% mutate(Type="Broad")
data.2 <- load_data(sp, "cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.1, data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Broad", "Cell", "Metacell")))

p.Xs <- plot_bar_FC(data, 2.0, 0.8, plot.title=gsub("_"," ", paste(sp,"",sep='')))
p.Xs
```

![](scRNA.Barchart_files/figure-html/plot_Xenia_sp_CTEAv1-1.png)<!-- -->


```r
sp <- "Hydra_vulgaris_MIJPv3"

data.1 <- load_data(sp, "broad_cell_type_gene_FC") %>% mutate(Type="Broad")
data.2 <- load_data(sp, "cell_type_gene_FC")       %>% mutate(Type="Cell")
data.3 <- load_data(sp, "metacell_gene_FC")        %>% mutate(Type="Metacell")

data <- rbind(data.1, data.2, data.3) %>%
  mutate(Type = factor(Type, levels=c("Broad","Cell", "Metacell")))

p.HV <- plot_bar_FC(data, 3.3, 0.8, plot.title=gsub("_"," ", paste(sp,"",sep='')))
p.HV
```

![](scRNA.Barchart_files/figure-html/plot_Hydra_vulgaris_MIJPv3-1.png)<!-- -->



```r
plot.list <- list(p.SP.adult, p.SP.polyp, p.SP.larva, p.NV, p.Xs, p.HV)
plots<-marrangeGrob(plot.list, nrow=1, ncol=1)
#plots
ggsave(filename="scRNA.Barchart.pdf", plots, width=21, height=29.7, units="cm")
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
