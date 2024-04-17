# Explore largest components

Following: https://elifesciences.org/articles/67667#s4

## Setup analysis directory

Link data files.

```bash
ln -s ../31_Network_Analysis/hhdb_hhr_reformatted_filtered.tab.90ID.Component.OGcount.tsv Component.OGcount.tsv
ln -s ../06_Final_Classifications/expression_results/Orthogroups.Run2.classification.RNA_results.tsv.gz
```

Setup bash environment.

```bash
conda activate py27
```



## Analysis

Copy files for components with >100 OGs.

```bash
awk -F'\t' '$2>100{print $1".SFDP.SBM"}' Component.OGcount.tsv > Component.selected.txt
cat Component.selected.txt \
  | while read P;
    do 
      ln -s ../31_Network_Analysis/${P}.gt.xz
      ln -s ../31_Network_Analysis/${P}.vertices.tsv
    done
```







