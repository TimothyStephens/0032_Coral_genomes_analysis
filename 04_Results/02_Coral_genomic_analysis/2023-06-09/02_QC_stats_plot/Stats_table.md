---
title: "Build table of QC results for each genome/transcriptome"
author: "Timothy Stephens"
date: "24/07/2023"
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
library(scales)
library(stringr)
library(writexl)
options(scipen = 999) #Prevent scientific notation
```





# Functions

Function to generate a formatted taxa names from samples sheet.

```r
# Format for plotting
# See https://yulab-smu.top/treedata-book/faq.html?q=rename#faq-formatting-label
format_sp_name = function(x){
  n = x[["genera"]]
  
  # Add species if present
  if(x[["species"]] != ""){
    n = paste(n, x[["species"]], sep=" ")
  }
  
  # Either print with extra info (not italics) or not
  if(x[["extra"]] != ""){
    e = x[["extra"]]
    lab = paste(n, e, sep=' ')
  } else {
    lab = paste(n, sep=' ')
  }
  return(lab)
}
```





# Load datasets into R

Load the cleaned and processed datasets into R for plotting.

```r
samples <- read.csv("../samples.txt", header=TRUE, sep='\t')
species.names <- samples$sample_id

genome.stats <- read.table(
  "all_genomes-01_stats-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>% 
  filter(SampleID %in% species.names)
transcriptome.stats <- read.table(
  "all_transcriptomes-01_stats-results.tsv", sep='\t', 
  header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>% 
  filter(SampleID %in% species.names)

genomes.busco4table.genome.eukaryota <- read.table("all_genomes-02_busco-genome.fa.busco_eukaryota_odb10-results4table.tsv",
  sep='\t', header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomes.busco4table.genome.metazoa <- read.table("all_genomes-02_busco-genome.fa.busco_metazoa_odb10-results4table.tsv", 
  sep='\t', header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomes.busco4table.protein.eukaryota <- read.table("all_genomes-02_busco-pep.faa.busco_eukaryota_odb10-results4table.tsv", 
  sep='\t', header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
genomes.busco4table.protein.metazoa <- read.table("all_genomes-02_busco-pep.faa.busco_metazoa_odb10-results4table.tsv", 
  sep='\t', header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)

transcriptomes.busco4table.protein.eukaryota <- read.table("all_transcriptomes-02_busco-pep.faa.busco_eukaryota_odb10-results4table.tsv", 
  sep='\t', header=TRUE, check.names=FALSE) %>% 
  filter(SampleID %in% species.names)
transcriptomes.busco4table.protein.metazoa <- read.table("all_transcriptomes-02_busco-pep.faa.busco_metazoa_odb10-results4table.tsv", 
  sep='\t', header=TRUE, check.names=FALSE)%>% 
  filter(SampleID %in% species.names)
```

Extract just the columns (in a logical order) from each dataset to output into an excel sheet for publication. Add extra empty columns which will separate the data from different sections (saves having to add this formatting later by hand). Also add "XX--" prefix to BUSCO columns names to make then unique when joining the tables. Will remove later.


## Genomes stats


```r
g1 <- genome.stats %>%
  select("SampleID", 
         "Total scaffold length (bp)", 
         "Total contig length (bp)",
         "Number of scaffolds",
         "Number of contigs",
         "N50 of scaffolds (bp)",
         "N50 of contigs (bp)",
         "Percent gaps",
         "Percent GC") %>%
  mutate("Genome assembly stats"="") %>%
  relocate(last_col(), .before = "Total scaffold length (bp)") %>%
  mutate(across(c("Total scaffold length (bp)", 
                  "Total contig length (bp)",
                  "Number of scaffolds",
                  "Number of contigs",
                  "N50 of scaffolds (bp)",
                  "N50 of contigs (bp)"), comma))

g2 <- genomes.busco4table.genome.metazoa %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"genome\" completeness (v5.0; metazoa_odb10; 954 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("GM--Complete (no.)" = "Complete") %>%
  rename("GM--  Single-copy (no.)" = "Single-copy") %>%
  rename("GM--  Duplicated (no.)" = "Duplicated") %>%
  rename("GM--Fragmented (no.)" = "Fragmented") %>%
  rename("GM--Missing (no.)" = "Missing")
  
g3 <- genomes.busco4table.genome.eukaryota %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"genome\" completeness (v5.0; eukaryota_odb10; 255 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("GE--Complete (no.)" = "Complete") %>%
  rename("GE--  Single-copy (no.)" = "Single-copy") %>%
  rename("GE--  Duplicated (no.)" = "Duplicated") %>%
  rename("GE--Fragmented (no.)" = "Fragmented") %>%
  rename("GE--Missing (no.)" = "Missing")
 
g4 <- genome.stats %>%
  select("SampleID", 
         "Number of genes", 
         "Average gene length (bp)",
         "Gene percent GC",
         "Average transcript length (bp)",
         "Average number of CDS per gene/transcript",
         "Average CDS length (bp)",
         "Number of single-CDS transcripts",
         "Percent single-CDS transcripts",
         "CDS percent GC",
         "Number of introns",
         "Average intron length (bp)",
         "Intron percent GC",
         "Number of intergenic regions",
         "Average intergenic region length (bp)",
         "Intergenic region percent GC") %>%
  mutate("Protein-coding gene stats (combined CDS+introns)"="") %>%
  relocate(last_col(), .before = "Number of genes") %>%
  mutate("Protein-coding transcript stats (based on CDS features)"="") %>%
  relocate(last_col(), .before = "Average transcript length (bp)") %>%
  mutate("CDS stats"="") %>%
  relocate(last_col(), .before = "Average CDS length (bp)") %>%
  mutate("Intron stats (predicted between CDS features)"="") %>%
  relocate(last_col(), .before = "Number of introns") %>%
  mutate("Intergenic stats (predicted between genes built from CDS features)"="") %>%
  relocate(last_col(), .before = "Number of intergenic regions") %>%
  mutate(across(c("Number of genes", 
                  "Average gene length (bp)",
                  "Average transcript length (bp)",
                  "Average number of CDS per gene/transcript",
                  "Average CDS length (bp)",
                  "Number of single-CDS transcripts",
                  "Number of introns",
                  "Average intron length (bp)",
                  "Number of intergenic regions",
                  "Average intergenic region length (bp)"), comma))

g5 <- genomes.busco4table.protein.metazoa %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; metazoa_odb10; 954 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PM--Complete (no.)" = "Complete") %>%
  rename("PM--  Single-copy (no.)" = "Single-copy") %>%
  rename("PM--  Duplicated (no.)" = "Duplicated") %>%
  rename("PM--Fragmented (no.)" = "Fragmented") %>%
  rename("PM--Missing (no.)" = "Missing")
  
g6 <- genomes.busco4table.protein.eukaryota %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; eukaryota_odb10; 255 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PE--Complete (no.)" = "Complete") %>%
  rename("PE--  Single-copy (no.)" = "Single-copy") %>%
  rename("PE--  Duplicated (no.)" = "Duplicated") %>%
  rename("PE--Fragmented (no.)" = "Fragmented") %>%
  rename("PE--Missing (no.)" = "Missing")

g <- data.frame(SampleID=samples$sample_id)
g$Name <- apply(samples, 1, format_sp_name)

g <- merge(g, g1, by="SampleID", sort=FALSE, no.dups=FALSE)
g <- merge(g, g2, by="SampleID", sort=FALSE, no.dups=FALSE)
g <- merge(g, g3, by="SampleID", sort=FALSE, no.dups=FALSE)
g <- merge(g, g4, by="SampleID", sort=FALSE, no.dups=FALSE)
g <- merge(g, g5, by="SampleID", sort=FALSE, no.dups=FALSE)
g <- merge(g, g6, by="SampleID", sort=FALSE, no.dups=FALSE)

g <- g %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(rowname = str_remove(rowname, "GM--")) %>%
  mutate(rowname = str_remove(rowname, "GE--")) %>%
  mutate(rowname = str_remove(rowname, "PM--")) %>%
  mutate(rowname = str_remove(rowname, "PE--"))
```


```r
write_xlsx(g, "all_genomes-Combined-results.xlsx")
```


## Transcriptomes stats


```r
t1 <- transcriptome.stats %>%
  select("SampleID", 
         "Number of CDS", 
         "Total CDS length (bp)",
         "Average CDS length (bp)",
         "N50 of CDS (bp)",
         "Percent GC") %>%
  mutate("Transcriptome assembly stats"="") %>%
  relocate(last_col(), .before = "Number of CDS") %>%
  mutate(across(c("Number of CDS",
                  "Total CDS length (bp)",
                  "Average CDS length (bp)",
                  "N50 of CDS (bp)",
                  "Percent GC"), comma))

t2 <- transcriptomes.busco4table.protein.metazoa %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; metazoa_odb10; 954 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PM--Complete (no.)" = "Complete") %>%
  rename("PM--  Single-copy (no.)" = "Single-copy") %>%
  rename("PM--  Duplicated (no.)" = "Duplicated") %>%
  rename("PM--Fragmented (no.)" = "Fragmented") %>%
  rename("PM--Missing (no.)" = "Missing")

t3 <- transcriptomes.busco4table.protein.eukaryota %>%
  select(-c("Total")) %>%
  mutate("BUSCO \"protein\" completeness (v5.0; eukaryota_odb10; 255 total genes)"="") %>%
  relocate(last_col(), .before = "Complete") %>%
  rename("PE--Complete (no.)" = "Complete") %>%
  rename("PE--  Single-copy (no.)" = "Single-copy") %>%
  rename("PE--  Duplicated (no.)" = "Duplicated") %>%
  rename("PE--Fragmented (no.)" = "Fragmented") %>%
  rename("PE--Missing (no.)" = "Missing")

t <- data.frame(SampleID=samples$sample_id)
t$Name <- apply(samples, 1, format_sp_name)

t <- merge(t, t1, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t, t2, by="SampleID", sort=FALSE, no.dups=FALSE)
t <- merge(t, t3, by="SampleID", sort=FALSE, no.dups=FALSE)

t <- t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(rowname = str_remove(rowname, "PM--")) %>%
  mutate(rowname = str_remove(rowname, "PE--"))
```


```r
write_xlsx(t, "all_transcriptomes-Combined-results.xlsx")
```





# Session Info


```r
sessionInfo()
```

```
## R version 4.2.3 (2023-03-15)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Ventura 13.2.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] writexl_1.4.2 stringr_1.5.0 scales_1.2.1  dplyr_1.1.1   tibble_3.2.1 
## 
## loaded via a namespace (and not attached):
##  [1] rstudioapi_0.14  knitr_1.42       magrittr_2.0.3   munsell_0.5.0   
##  [5] tidyselect_1.2.0 colorspace_2.1-0 R6_2.5.1         rlang_1.1.0     
##  [9] fastmap_1.1.1    fansi_1.0.4      tools_4.2.3      xfun_0.38       
## [13] utf8_1.2.3       cli_3.6.1        withr_2.5.0      jquerylib_0.1.4 
## [17] htmltools_0.5.5  yaml_2.3.7       digest_0.6.31    lifecycle_1.0.3 
## [21] sass_0.4.5       vctrs_0.6.1      glue_1.6.2       cachem_1.0.7    
## [25] evaluate_0.20    rmarkdown_2.21   stringi_1.7.12   compiler_4.2.3  
## [29] bslib_0.4.2      pillar_1.9.0     generics_0.1.3   jsonlite_1.8.4  
## [33] pkgconfig_2.0.3
```
