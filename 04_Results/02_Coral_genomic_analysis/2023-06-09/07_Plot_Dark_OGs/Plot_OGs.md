---
title: "Plot selected OGs against species tree"
author: "Timothy Stephens"
date: "22 January, 2024"
output: 
  html_document:
    keep_md: yes
---



## Setup

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
## Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods
## for mapping and visualizing associated data on phylogeny using ggtree.
## Molecular Biology and Evolution. 2018, 35(12):3041-3043.
## doi:10.1093/molbev/msy194
## 
## Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
## Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu. Ggtree: A serialized data
## object for visualization of a phylogenetic tree and annotation data.
## iMeta 2022, 1(4):e56. doi:10.1002/imt2.56
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
## Guangchuang Yu.  Data Integration, Manipulation and Visualization of
## Phylogenetic Trees (1st edition). Chapman and Hall/CRC. 2022,
## doi:10.1201/9781003279242
## 
## Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods
## for mapping and visualizing associated data on phylogeny using ggtree.
## Molecular Biology and Evolution. 2018, 35(12):3041-3043.
## doi:10.1093/molbev/msy194
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
library(tidyr)
```

```
## 
## Attaching package: 'tidyr'
```

```
## The following object is masked from 'package:reshape2':
## 
##     smiths
```

```
## The following object is masked from 'package:ggtree':
## 
##     expand
```

```r
library(aplot)
library(knitr)
options(scipen = 999) #Prevent scientific notation
```





## Functions

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



Function to plot tree + heatmap

```r
get_colors <- function(seq.count.max){
  color.list <- c("white", colorRampPalette(c("#ddf1da","#d53e4f"))(5))
  for (i in 1:seq.count.max){
    color.list <- c(color.list, "#d53e4f")
  }
  return(color.list[1:seq.count.max])
}

plot_tree_heatmap <- function(p,
                              OG.class, OG.seqs, OG.annots,
                              select.best_strata_taxa, select.designation, 
                              select.best_strata_species_in_OG_proportion_total,
                              plot.heatmap.rownames=FALSE){
  # Select OGs to plot
  selected.OGs.class <- OG.class %>%
    filter(best_strata_taxa == select.best_strata_taxa) %>%
    filter(designation == select.designation) %>%
    filter(best_strata_species_in_OG_proportion_total >= select.best_strata_species_in_OG_proportion_total) %>%
    arrange(orthogroup_id)
  selected.OGs <- selected.OGs.class %>%
    select(orthogroup_id)
  
  # Check that we have OGs to plot
  if(length(selected.OGs$orthogroup_id) == 0){
    print("No OGs found using the given filtering cutoffs. Returning nothing.")
    return()
  }
  
  # Get sequence membership and count number of seqs per genome/transcriptome
  selected.OGs.seqs <- OG.seqs %>%
    filter(orthogroup_id %in% selected.OGs$orthogroup_id) %>%
    arrange(orthogroup_id, species_id, sequence_id)
  selected.OGs.seqs.counts <- selected.OGs.seqs %>%
    dplyr::group_by(orthogroup_id, species_id) %>%
    dplyr::summarise(n_seq = n(), .groups="keep") %>%
    pivot_wider(names_from = species_id, values_from = n_seq, values_fill=0)
  
  # Get max number of sequences per species (used for setting color range during plotting)
  seq.count.max <- max(selected.OGs.seqs.counts %>% ungroup() %>% select(-orthogroup_id))
  
  # Get sequence annotations
  selected.OGs.annots <- merge(selected.OGs.seqs, OG.annots, by="sequence_id", sort=FALSE)
  
  # Add missing genomes/transcriptomes to matrix - need to do since we are selecting OGs specific to taxonomic subset.
  all.names <- samples$sample_id
  not.in.matrix <- all.names[!all.names %in% colnames(selected.OGs.seqs.counts)]
  for (n in not.in.matrix){
    selected.OGs.seqs.counts[, n] <- 0
  }
  
  # Write results for selected OGs
  write.table(selected.OGs.class, 
            paste("selected_", select.best_strata_taxa, "_", select.designation, "_", select.best_strata_species_in_OG_proportion_total, ".classification.tsv", sep=''), 
            sep='\t',
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE,
            )
  write.table(selected.OGs.seqs, 
          paste("selected_", select.best_strata_taxa, "_", select.designation, "_", select.best_strata_species_in_OG_proportion_total, ".long.tsv", sep=''), 
          sep='\t',
          quote = FALSE, 
          row.names = FALSE, 
          col.names = TRUE,
          )
  write.table(selected.OGs.seqs.counts, 
        paste("selected_", select.best_strata_taxa, "_", select.designation, "_", select.best_strata_species_in_OG_proportion_total, ".seq_counts.tsv", sep=''), 
        sep='\t',
        quote = FALSE, 
        row.names = FALSE, 
        col.names = TRUE,
        )
  write.table(selected.OGs.annots, 
        paste("selected_", select.best_strata_taxa, "_", select.designation, "_", select.best_strata_species_in_OG_proportion_total, ".sequences.nr.top_hits.tsv", sep=''), 
        sep='\t',
        quote = FALSE, 
        row.names = FALSE, 
        col.names = TRUE,
        )
  
  # Order or OGs in heatmap based on clustering of values
  if (nrow(selected.OGs.seqs.counts) >= 2) {
    data <- scale(selected.OGs.seqs.counts[,-1])
    ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
    orthogroup_id.ord <- selected.OGs.seqs.counts$orthogroup_id[ord]
  } else {
    orthogroup_id.ord <- selected.OGs.seqs.counts$orthogroup_id
  }
  
  # Melt matrix into 3-columns
  selected.OGs.seqs.counts <- melt(selected.OGs.seqs.counts, id.vars = "orthogroup_id") %>%
    mutate(value = factor(value)) %>%
    mutate(orthogroup_id = factor(orthogroup_id, levels=orthogroup_id.ord))
  
  # Create additional column where zeros are NA (helpful during plotting since we can specifically set NA as white)
  selected.OGs.seqs.counts$value_nas <- ifelse(selected.OGs.seqs.counts$value==0, NA, selected.OGs.seqs.counts$value)
  
  #Create heatmap
  plot_size_ratio <- 0.2
  
  scaling.factor <- min(c((20 / nrow(selected.OGs)), 1))
  
  plot_family <- "sans"
  plot_title <- paste(length(selected.OGs$orthogroup_id), ' "', select.designation, '" OGs that are specific to "', select.best_strata_taxa, '" and present in >=', select.best_strata_species_in_OG_proportion_total, '% of samples in this group', sep='')
  hm <- ggplot(selected.OGs.seqs.counts, aes(x=orthogroup_id, y=variable)) + 
    geom_tile(aes(fill=cut(value_nas, breaks=0:seq.count.max, labels=0:(seq.count.max-1)))) + 
    scale_fill_manual(drop=FALSE, values=get_colors(seq.count.max), na.value="white", name="No. Genes") +
    theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle(plot_title) +
    theme(axis.text.x = element_text(family = plot_family, 
                                     face = "bold", 
                                     size = rel(scaling.factor * 4)*plot_size_ratio, 
                                     colour = "black",
                                     angle = 90,
                                     )
          ) +
    theme(axis.text.y = element_text(family = plot_family, 
                                     face = "bold", 
                                     size = rel(1.66)*plot_size_ratio, 
                                     colour = "black",
                                     )
          ) +
    theme(plot.title = element_text(family = plot_family, 
                                     face = "bold", 
                                     size = rel(3)*plot_size_ratio, 
                                     colour = "black",
                                     hjust = 0.5,
                                     )
          )
  
  # Remove heatmap row names?
  if(!plot.heatmap.rownames){
    hm <- hm + theme(axis.text.y = element_blank())
  }
  
  # Plot heatmap against tree
  p <- hm %>% insert_left(p, width=1)
  
  # Return plot
  return(p)
}
```





# Load datasets into R

Load tree of samples (to plot).

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



Load OG info.

```r
data.class   <- read.csv("../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz", header=TRUE, sep='\t')
data.seqs    <- read.csv("../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz", header=TRUE, sep='\t')
data.annots  <- read.csv("../06_Final_Classifications/Orthogroups.Run2.sequences.nr.top_hits.tsv.gz", header=TRUE, sep='\t')
```





# Plot dataset tree with colored clades

Global setting for plotting

```r
plot_size_ratio <- 0.2
plot.axis.text.y <- FALSE
plot.legend <- FALSE
```

Data.frame with clade positions in tree and colors.

```r
clade.colors <- data.frame(
  node  = c(132,        138,         144,         146,       148,        164,            171,            172,          183,            179,                130,        126,        122,          124,                134,       136,          133,             56,          57),
  fill  = c("blue",     "red",       "brown",     "pink",    "yellow",   "purple",       "green",        "orange",     "blue",         "red",             "purple",    "green",    "yellow",     "red",              "blue",    "blue",       "yellow",       "green",        "blue"),
  alpha = c(0.2,        0.4,         0.2,         0.4,       0.2,        0.2,            0.2,            0.2,          0.2,            0.2,                0.5,        0.5,        0.5,          0.5,                0.5,       0.1,          0.2,            0.4,          0.2),
  n     = c("Cnidaria", "Staurozoa", "Scyphozoa", "Cubozoa", "Hydrozoa", "Octocorallia", "Hexacorallia", "Actiniaria", "Scleractinia", "Corallimorpharia", "Placozoa", "Porifera", "Ctenophora", "Choanoflagellata", "Myxozoa", "Medusozoa",  "Endocnidozoa", "Zoantharia", "Scleractinia")
)
```

Format species tree. Color clades and rename taxa into a pretty format (i.e., genera+species italics + extra info not italics).

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
    extendto = 1.4,
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
    offset = 0.2,
    barsize = 0,
    fontface="bold",
    ) +
  geom_tiplab(aes(label=name), size=6*plot_size_ratio, linesize=0.2, parse=TRUE, align=TRUE, offset=0.01) +
  geom_tippoint(aes(color=data_type), size=6*plot_size_ratio, show.legend=FALSE) +
  hexpand(.1)

p
```

![](Plot_OGs_files/figure-html/species_tree_plotting-1.png)<!-- -->



Count number of OGs in each combination of groups.

```r
t <- data.class %>%
  filter(best_strata_species_in_OG_proportion_total >= 50) %>%
  filter(no_species > 2) %>%
  group_by(designation, best_strata, best_strata_taxa) %>%
  count() %>%
  arrange(designation, best_strata, desc(n))
kable(t)
```



|designation      |best_strata |best_strata_taxa   |    n|
|:----------------|:-----------|:------------------|----:|
|Dark-Restricted  |clade       |Opisthokonta       |    2|
|Dark-Restricted  |clade1      |Eumetazoa          |    5|
|Dark-Restricted  |class       |Anthozoa           |   23|
|Dark-Restricted  |genus       |Millepora          | 3254|
|Dark-Restricted  |genus       |Porites            |  377|
|Dark-Restricted  |genus       |Montipora          |  178|
|Dark-Restricted  |genus       |Hydra              |  107|
|Dark-Restricted  |genus       |Acropora           |  100|
|Dark-Restricted  |genus       |Pocillopora        |   69|
|Dark-Restricted  |kingdom     |Metazoa            |    1|
|Dark-Restricted  |order       |Stauromedusae      |  733|
|Dark-Restricted  |order       |Alcyonacea         |  182|
|Dark-Restricted  |order       |Anthoathecata      |   63|
|Dark-Restricted  |order       |Scleractinia       |   23|
|Dark-Restricted  |phylum      |Cnidaria           |   32|
|Dark-Restricted  |subclass    |Hydroidolina       |  247|
|Dark-Restricted  |subclass    |Hexacorallia       |  200|
|Dark-Restricted  |subclass    |Octocorallia       |  169|
|Dark-Restricted  |subclass    |Heteroscleromorpha |  112|
|Dark-Restricted  |suborder    |Capitata           |  453|
|Dark-Restricted  |suborder    |Holaxonia          |   71|
|Dark-Restricted  |suborder    |Filifera           |   50|
|Dark-Restricted  |suborder    |Myostaurida        |   41|
|Dark-Restricted  |suborder    |Amyostaurida       |   34|
|Dark-Restricted  |suborder    |Faviina            |    6|
|Dark-Restricted  |suborder    |Fungiina           |    5|
|Dark-Restricted  |suborder    |Astrocoeniina      |    2|
|Dark-Shared      |clade1      |Eumetazoa          |    2|
|Dark-Shared      |class       |Anthozoa           |    3|
|Dark-Shared      |genus       |Millepora          |  405|
|Dark-Shared      |genus       |Hydra              |   14|
|Dark-Shared      |genus       |Porites            |    8|
|Dark-Shared      |genus       |Montipora          |    3|
|Dark-Shared      |genus       |Pocillopora        |    2|
|Dark-Shared      |kingdom     |Metazoa            |    2|
|Dark-Shared      |order       |Stauromedusae      |   22|
|Dark-Shared      |order       |Alcyonacea         |    5|
|Dark-Shared      |order       |Anthoathecata      |    3|
|Dark-Shared      |order       |Scleractinia       |    1|
|Dark-Shared      |phylum      |Cnidaria           |    7|
|Dark-Shared      |subclass    |Heteroscleromorpha |   19|
|Dark-Shared      |subclass    |Hydroidolina       |    9|
|Dark-Shared      |subclass    |Hexacorallia       |    8|
|Dark-Shared      |subclass    |Octocorallia       |    3|
|Dark-Shared      |suborder    |Capitata           |   21|
|Dark-Shared      |suborder    |Amyostaurida       |    4|
|Dark-Shared      |suborder    |Filifera           |    3|
|Dark-Shared      |suborder    |Holaxonia          |    2|
|Dark-Shared      |suborder    |Astrocoeniina      |    1|
|Dark-Shared      |suborder    |Faviina            |    1|
|Light-Restricted |clade       |Opisthokonta       |    1|
|Light-Restricted |clade1      |Eumetazoa          |    1|
|Light-Restricted |class       |Anthozoa           |   36|
|Light-Restricted |genus       |Millepora          |   70|
|Light-Restricted |genus       |Porites            |   67|
|Light-Restricted |genus       |Hydra              |   46|
|Light-Restricted |genus       |Acropora           |   33|
|Light-Restricted |genus       |Pocillopora        |   32|
|Light-Restricted |genus       |Montipora          |   28|
|Light-Restricted |kingdom     |Metazoa            |    3|
|Light-Restricted |order       |Alcyonacea         |  121|
|Light-Restricted |order       |Stauromedusae      |   55|
|Light-Restricted |order       |Anthoathecata      |   20|
|Light-Restricted |order       |Scleractinia       |   11|
|Light-Restricted |phylum      |Cnidaria           |   38|
|Light-Restricted |subclass    |Hexacorallia       |  163|
|Light-Restricted |subclass    |Octocorallia       |   82|
|Light-Restricted |subclass    |Hydroidolina       |   52|
|Light-Restricted |subclass    |Heteroscleromorpha |   42|
|Light-Restricted |suborder    |Holaxonia          |   88|
|Light-Restricted |suborder    |Capitata           |   57|
|Light-Restricted |suborder    |Filifera           |   12|
|Light-Restricted |suborder    |Faviina            |    8|
|Light-Restricted |suborder    |Myostaurida        |    7|
|Light-Restricted |suborder    |Amyostaurida       |    4|
|Light-Shared     |clade       |Opisthokonta       | 4995|
|Light-Shared     |clade1      |Eumetazoa          |  526|
|Light-Shared     |class       |Anthozoa           |  224|
|Light-Shared     |genus       |Millepora          | 2226|
|Light-Shared     |genus       |Porites            |  212|
|Light-Shared     |genus       |Hydra              |  209|
|Light-Shared     |genus       |Montipora          |   83|
|Light-Shared     |genus       |Acropora           |   27|
|Light-Shared     |genus       |Pocillopora        |   27|
|Light-Shared     |kingdom     |Metazoa            | 1860|
|Light-Shared     |order       |Stauromedusae      |  485|
|Light-Shared     |order       |Alcyonacea         |  234|
|Light-Shared     |order       |Anthoathecata      |   79|
|Light-Shared     |order       |Scleractinia       |   23|
|Light-Shared     |phylum      |Cnidaria           |  740|
|Light-Shared     |subclass    |Hexacorallia       |  465|
|Light-Shared     |subclass    |Heteroscleromorpha |  297|
|Light-Shared     |subclass    |Hydroidolina       |  187|
|Light-Shared     |subclass    |Octocorallia       |  161|
|Light-Shared     |suborder    |Capitata           |  394|
|Light-Shared     |suborder    |Holaxonia          |  105|
|Light-Shared     |suborder    |Filifera           |   36|
|Light-Shared     |suborder    |Amyostaurida       |   22|
|Light-Shared     |suborder    |Myostaurida        |   22|
|Light-Shared     |suborder    |Faviina            |   11|
|Light-Shared     |suborder    |Fungiina           |    9|


Plot heatmap of species presence in each Scleractinia Dark-Restricted OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Scleractinia", "Dark-Restricted", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Scleractinia_Dark-Restricted-1.png)<!-- -->



Plot heatmap of species presence in each Hexacorallia Dark-Restricted OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Hexacorallia", "Dark-Restricted", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Hexacorallia_Dark-Restricted-1.png)<!-- -->



Plot heatmap of species presence in each Cnidaria Dark-Restricted OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Cnidaria", "Dark-Restricted", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Cnidaria_Dark-Restricted-1.png)<!-- -->



Plot heatmap of species presence in each Opisthokonta Dark-Restricted OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Opisthokonta", "Dark-Restricted", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Opisthokonta_Dark-Restricted-1.png)<!-- -->





Plot heatmap of species presence in each Scleractinia Dark-Shared OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Scleractinia", "Dark-Shared", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Scleractinia_Dark-Shared-1.png)<!-- -->



Plot heatmap of species presence in each Hexacorallia Dark-Shared OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Hexacorallia", "Dark-Shared", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Hexacorallia_Dark-Shared-1.png)<!-- -->



Plot heatmap of species presence in each Cnidaria Dark-Shared OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Cnidaria", "Dark-Shared", 50)
```

![](Plot_OGs_files/figure-html/data_plots-Cnidaria_Dark-Shared-1.png)<!-- -->



Plot heatmap of species presence in each Opisthokonta Dark-Shared OG.

```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Opisthokonta", "Dark-Shared", 50)
```

```
## [1] "No OGs found using the given filtering cutoffs. Returning nothing."
```

```
## NULL
```










Plot heatmap of species presence in each Scleractinia Dark-Restricted OG (>=5%).

```r
t <- data.class %>%
  filter(best_strata_species_in_OG_proportion_total >= 5) %>%
  filter(no_species > 2) %>%
  group_by(designation, best_strata, best_strata_taxa) %>%
  count() %>%
  arrange(designation, best_strata, desc(n)) %>%
  filter(best_strata_taxa == "Scleractinia")
kable(t)
```



|designation      |best_strata |best_strata_taxa |    n|
|:----------------|:-----------|:----------------|----:|
|Dark-Restricted  |order       |Scleractinia     |  874|
|Dark-Shared      |order       |Scleractinia     |   49|
|Light-Restricted |order       |Scleractinia     |  644|
|Light-Shared     |order       |Scleractinia     | 1166|


```r
plot_tree_heatmap(p, data.class, data.seqs, data.annots, "Scleractinia", "Dark-Restricted", 5)
```

![](Plot_OGs_files/figure-html/data_plots-Scleractinia_Dark-Restricted_5-1.png)<!-- -->






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
##  [1] knitr_1.45     tidyr_1.3.0    ggpubr_0.6.0   glue_1.6.2     stringr_1.5.1 
##  [6] ggrepel_0.9.4  scales_1.2.1   reshape2_1.4.4 dplyr_1.1.3    tibble_3.2.1  
## [11] treeio_1.25.3  aplot_0.2.2    ggtree_3.8.2   ggplot2_3.4.4 
## 
## loaded via a namespace (and not attached):
##  [1] yulab.utils_0.1.0  sass_0.4.7         utf8_1.2.4         generics_0.1.3    
##  [5] rstatix_0.7.2      ggplotify_0.1.2    stringi_1.8.1      lattice_0.22-5    
##  [9] digest_0.6.33      magrittr_2.0.3     evaluate_0.23      grid_4.3.1        
## [13] fastmap_1.1.1      plyr_1.8.9         jsonlite_1.8.7     ape_5.7-1         
## [17] backports_1.4.1    purrr_1.0.2        fansi_1.0.5        lazyeval_0.2.2    
## [21] jquerylib_0.1.4    abind_1.4-5        cli_3.6.1          rlang_1.1.2       
## [25] munsell_0.5.0      tidytree_0.4.5     withr_2.5.2        cachem_1.0.8      
## [29] yaml_2.3.7         tools_4.3.1        parallel_4.3.1     ggsignif_0.6.4    
## [33] memoise_2.0.1      colorspace_2.1-0   broom_1.0.5        vctrs_0.6.4       
## [37] R6_2.5.1           gridGraphics_0.5-1 lifecycle_1.0.4    car_3.1-2         
## [41] fs_1.6.3           ggfun_0.1.3        pkgconfig_2.0.3    pillar_1.9.0      
## [45] bslib_0.5.1        gtable_0.3.4       Rcpp_1.0.11        highr_0.10        
## [49] xfun_0.41          tidyselect_1.2.0   rstudioapi_0.15.0  farver_2.1.1      
## [53] htmltools_0.5.7    nlme_3.1-163       patchwork_1.1.3    labeling_0.4.3    
## [57] carData_3.0-5      rmarkdown_2.25     compiler_4.3.1
```


