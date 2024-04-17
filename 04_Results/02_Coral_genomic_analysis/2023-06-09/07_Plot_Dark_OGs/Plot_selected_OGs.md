---
title: "Plot selected OGs against species tree"
author: "Timothy Stephens"
date: "22 February, 2024"
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
## S Xu, Z Dai, P Guo, X Fu, S Liu, L Zhou, W Tang, T Feng, M Chen, L
## Zhan, T Wu, E Hu, Y Jiang, X Bo, G Yu. ggtreeExtra: Compact
## visualization of richly annotated phylogenetic data. Molecular Biology
## and Evolution. 2021, 38(9):4039-4042. doi: 10.1093/molbev/msab166
## 
## LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
## Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
## for phylogenetic tree input and output with richly annotated and
## associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
## doi: 10.1093/molbev/msz240
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
## G Yu. Data Integration, Manipulation and Visualization of Phylogenetic
## Trees (1st ed.). Chapman and Hall/CRC. 2022. ISBN: 9781032233574
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
      # geom_text_repel(aes(label = format_values_1(Count)), # Add two spaces in from of text so it centers off the number not number%
      #          size = rel(5)*plot_size_ratio,
      #          fontface = "bold",
      #          max.overlaps = Inf, #Set max.overlaps = Inf to override this behavior and always show all labels, regardless of too many overlaps.
      #          position = position_stack(vjust=0.5, reverse = TRUE),
      #          xlim = c(-Inf, Inf), # Set xlim or ylim to Inf or -Inf to disable repulsion away from the edges of the panel.
      #          ylim = c(-Inf, Inf),
      #          point.size = NA, #size of each point for each text label; Set point.size = NA to prevent label repulsion away from data points.
      #          force=0.5, # force of repulsion between overlapping text labels
      #          box.padding=0, # padding around the text label
      #          direction="y") + # move text labels “both” (default), “x”, or “y” directions
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
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
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
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
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
                              plot.heatmap.rownames=FALSE,
                              plot_size_ratio=0.2){
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
  
  # Melt matrix into 3-columns
  selected.OGs.seqs.counts <- melt(selected.OGs.seqs.counts, id.vars = "orthogroup_id") %>%
    mutate(value = factor(value))
  
  # Create additional column where zeros are NA (helpful during plotting since we can specifically set NA as white)
  selected.OGs.seqs.counts$value_nas <- ifelse(selected.OGs.seqs.counts$value==0,NA, selected.OGs.seqs.counts$value)
  
  scaling.factor <- min(c((20 / nrow(selected.OGs)), 1))
  
  plot_family <- "sans"
  plot_title <- paste(select.best_strata_taxa, select.designation, sep=' ')
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
          ) +
    theme(legend.text = element_text(family = plot_family, size = rel(1.2)*plot_size_ratio)) +
    theme(legend.key.size = unit(1.5*plot_size_ratio,"line")) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # Remove heatmap row names?
  if(!plot.heatmap.rownames){
    hm <- hm + theme(axis.text.y = element_blank())
  }
  
  # Plot heatmap against tree
  #p <- hm %>% insert_left(p, width=1)
  
  # Return plot
  return(hm)
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

genome.stats <- read.table(
  "../02_QC_stats_plot/all_genomes-01_stats-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>% 
  filter(SampleID %in% species.names)
transcriptome.stats <- read.table(
  "../02_QC_stats_plot/all_transcriptomes-01_stats-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>% 
  filter(SampleID %in% species.names) %>%
  rename(`Number of genes` = `Number of CDS`)

genomeData.busco.genome.eukaryota <- read.table(
  "../02_QC_stats_plot/all_genomes-02_busco-genome.fa.busco_eukaryota_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomeData.busco.genome.metazoa <- read.table(
  "../02_QC_stats_plot/all_genomes-02_busco-genome.fa.busco_metazoa_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomeData.busco.protein.eukaryota <- read.table(
  "../02_QC_stats_plot/all_genomes-02_busco-pep.faa.busco_eukaryota_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomeData.busco.protein.metazoa <- read.table(
  "../02_QC_stats_plot/all_genomes-02_busco-pep.faa.busco_metazoa_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)

transcriptomeData.busco.protein.eukaryota <- read.table(
  "../02_QC_stats_plot/all_transcriptomes-02_busco-pep.faa.busco_eukaryota_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
transcriptomeData.busco.protein.metazoa <- read.table(
  "../02_QC_stats_plot/all_transcriptomes-02_busco-pep.faa.busco_metazoa_odb10-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
```

Merge transcriptome and genome datasets together so that all species from the tree are represented (species without this data get zero values).

```r
# Genome/Transcriptome stats
seq.stats <- merge(genome.stats, transcriptome.stats, all=TRUE) %>% replace(is.na(.), 0)

# BUSCO results
busco.protein.eukaryota <- rbind(genomeData.busco.protein.eukaryota, transcriptomeData.busco.protein.eukaryota)
busco.protein.metazoa <- rbind(genomeData.busco.protein.metazoa, transcriptomeData.busco.protein.metazoa)

busco.genome.eukaryota <- data.frame(SampleID=species.names)
busco.genome.eukaryota <- merge(x=busco.genome.eukaryota, y=genomeData.busco.genome.eukaryota, all.x=TRUE) %>% replace(is.na(.), 0)

busco.genome.metazoa <- data.frame(SampleID=species.names)
busco.genome.metazoa <- merge(x=busco.genome.metazoa, y=genomeData.busco.genome.metazoa, all.x=TRUE) %>% replace(is.na(.), 0)
```

Load OG info.

```r
data.class   <- read.csv("../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz", header=TRUE, sep='\t')
data.seqs    <- read.csv("../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz", header=TRUE, sep='\t')
data.annots  <- read.csv("../06_Final_Classifications/Orthogroups.Run2.sequences.nr.top_hits.tsv.gz", header=TRUE, sep='\t')
```





# Plot data phylogeny

Global setting for plotting

```r
plot_size_ratio <- 0.2
plot.axis.text.y <- FALSE
plot.legend <- FALSE
plot_family <- "sans"
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

#p
```



Count number of sequences in each designation per species.

```r
data.seqs.class <- merge(data.seqs, data.class)
data.species.OGtype.count <- data.seqs.class %>% 
  select(orthogroup_id, species_id, sequence_id, designation) %>% 
  unique() %>% 
  group_by(species_id, designation) %>%
  count()
```

Count number of sequences in each designation per species.

```r
data.seqs.class %>% 
  select(orthogroup_id, species_id, sequence_id, designation) %>% 
  unique() %>% 
  group_by(designation) %>%
  count()
```

```
## # A tibble: 4 × 2
## # Groups:   designation [4]
##   designation            n
##   <chr>              <int>
## 1 Dark-Restricted   384502
## 2 Dark-Shared        51367
## 3 Light-Restricted   65872
## 4 Light-Shared     3276409
```


```r
p.OGtype.absolute <- data.species.OGtype.count %>%
  ggplot(aes(y = n, x = species_id)) +
    geom_bar(aes(fill=designation),
             position = position_stack(),
             stat="identity") +
    scale_fill_manual(values = c("#984ea3", "#bebada", "#ff7f00", "#fdbf6f")) +
    scale_y_continuous(labels = comma) +
    theme_gray(base_size = 8) + 
    ggtitle("No of genes") +
    xlab("") +
    ylab("No. of genes") +
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
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    #geom_text(aes(label = n), size=1.0, position = position_stack(vjust = 0.5, reverse = TRUE)) +
    theme(axis.text.y=element_blank()) +
    coord_flip() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p.OGtype.relative <- data.species.OGtype.count %>%
  ggplot(aes(y = n, x = species_id)) +
    geom_bar(aes(fill=designation),
             stat="identity",
             position = 'fill') +
    scale_fill_manual(values = c("#984ea3", "#bebada", "#ff7f00", "#fdbf6f")) +
    scale_y_continuous(labels = comma) +
    theme_gray(base_size = 8) + 
    ggtitle("Relative No. of genes") +
    xlab("") +
    ylab("Relative No. of genes") +
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
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    #geom_text(aes(label = n), size=1.0, position = position_stack(vjust = 0.5, reverse = TRUE)) +
    theme(axis.text.y=element_blank()) +
    coord_flip() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p.gene.stats.count <- plot_barchart(data=seq.stats,
                                    column.id="Number of genes", 
                                    plot_title="Protein-coding genes", 
                                    xlab="Number of genes", 
                                    bar.fill="#d9d9d9", 
                                    plot.axis.text.y=plot.axis.text.y,
                                    plot_size_ratio=plot_size_ratio,
                                    add.geom_text=FALSE)

p.busco.protein.eukaryota <- plot_BUSCO_percentages(
                                    data=busco.protein.eukaryota, 
                                    plot_title="BUSCO Eukaryota - Protein", 
                                    plot.axis.text.y=plot.axis.text.y,
                                    plot.legend=plot.legend,
                                    plot_size_ratio=plot_size_ratio)


# Plot heatmap of species presence in each Scleractinia Dark-Shared OG.
p.ScleractiniaDarkshared     <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Scleractinia", "Dark-Shared",     50, plot_size_ratio=0.022)

# Plot heatmap of species presence in each Scleractinia Dark-Restricted OG.
p.ScleractiniaDarkRestricted <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Scleractinia", "Dark-Restricted", 50, plot_size_ratio=0.025)

# Plot heatmap of species presence in each Hexacorallia Dark-Shared OG.
p.HexacoralliaDarkShared     <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Hexacorallia", "Dark-Shared",     50, plot_size_ratio=0.020)

# Plot heatmap of species presence in each Hexacorallia Dark-Restricted OG.
p.HexacoralliaDarkRestricted <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Hexacorallia", "Dark-Restricted", 50, plot_size_ratio=0.200)

# Plot heatmap of species presence in each Cnidaria Dark-Shared OG.
p.CnidariaDarkShared         <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Cnidaria", "Dark-Shared",         50, plot_size_ratio=0.020)

# Plot heatmap of species presence in each Cnidaria Dark-Restricted OG.
p.CnidariaDarkRestricted     <- plot_tree_heatmap(combined.plot, data.class, data.seqs, data.annots, "Cnidaria", "Dark-Restricted",     50, plot_size_ratio=0.032)
```



```r
p.busco.protein.eukaryota %>% 
  insert_left(p, width=6) %>%
  insert_right(p.gene.stats.count) %>%
  insert_right(p.OGtype.absolute) %>%
  insert_right(p.OGtype.relative) %>%
  insert_right(p.CnidariaDarkRestricted,     width = 1.00) %>%
  insert_right(p.CnidariaDarkShared,         width = 0.25) %>%
  insert_right(p.HexacoralliaDarkRestricted, width = 6.00) %>%
  insert_right(p.HexacoralliaDarkShared,     width = 0.25) %>%
  insert_right(p.ScleractiniaDarkRestricted, width = 0.70) %>%
  insert_right(p.ScleractiniaDarkshared,     width = 0.04)
```

![](Plot_selected_OGs_files/figure-html/plot_distribution_Dark-Restricted_OGs-1.png)<!-- -->





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
##  [1] knitr_1.45     tidyr_1.3.1    ggpubr_0.6.0   glue_1.7.0     stringr_1.5.1 
##  [6] ggrepel_0.9.5  scales_1.3.0   reshape2_1.4.4 dplyr_1.1.4    tibble_3.2.1  
## [11] treeio_1.25.3  aplot_0.2.2    ggtree_3.8.2   ggplot2_3.4.4 
## 
## loaded via a namespace (and not attached):
##  [1] yulab.utils_0.1.4  sass_0.4.8         utf8_1.2.4         generics_0.1.3    
##  [5] rstatix_0.7.2      ggplotify_0.1.2    stringi_1.8.3      lattice_0.22-5    
##  [9] digest_0.6.34      magrittr_2.0.3     evaluate_0.23      grid_4.3.1        
## [13] fastmap_1.1.1      plyr_1.8.9         jsonlite_1.8.8     ape_5.7-1         
## [17] backports_1.4.1    purrr_1.0.2        fansi_1.0.6        lazyeval_0.2.2    
## [21] jquerylib_0.1.4    abind_1.4-5        cli_3.6.2          rlang_1.1.3       
## [25] munsell_0.5.0      tidytree_0.4.6     withr_3.0.0        cachem_1.0.8      
## [29] yaml_2.3.8         tools_4.3.1        parallel_4.3.1     ggsignif_0.6.4    
## [33] memoise_2.0.1      colorspace_2.1-0   broom_1.0.5        vctrs_0.6.5       
## [37] R6_2.5.1           gridGraphics_0.5-1 lifecycle_1.0.4    car_3.1-2         
## [41] fs_1.6.3           ggfun_0.1.4        pkgconfig_2.0.3    pillar_1.9.0      
## [45] bslib_0.6.1        gtable_0.3.4       Rcpp_1.0.12        highr_0.10        
## [49] xfun_0.41          tidyselect_1.2.0   rstudioapi_0.15.0  farver_2.1.1      
## [53] htmltools_0.5.7    nlme_3.1-164       patchwork_1.2.0    labeling_0.4.3    
## [57] carData_3.0-5      rmarkdown_2.25     compiler_4.3.1
```


