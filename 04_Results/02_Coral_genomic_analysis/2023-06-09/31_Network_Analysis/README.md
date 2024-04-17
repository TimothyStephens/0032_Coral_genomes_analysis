# Network of OG Deep Homology search

Following: https://elifesciences.org/articles/67667#s4

## Setup analysis directory

Link data files.

```bash
ln -s ../30_OG_Deep_Homology_Search/hhdb_hhr_reformatted_filtered.tab
ln -s ../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz
```

Setup bash environment.

```bash
conda activate graph-tool
export PYTHONPATH=/home/timothy/miniconda3/envs/graph-tool/lib/python3.11/site-packages

SCRIPTS=/home/timothy/GitHub/Network_Analysis_Scripts
```



## Run Network Analysis

```bash
./run_01_load_edges.sh
./run_02_extract_connected_components.sh
```

Get list of componenets to run (including full graph).

```bash
ls -1 hhdb_hhr_reformatted_filtered.tab.90ID*.gt.xz | sed -e 's/.gt.xz//' | sort | uniq > component_prefixes.txt
```

Run layout and structure calculations and plot each connected component.

```bash
./run_03_compute_layout_SFDP_springBlock.sh
./run_04_compute_structure_SBM.sh
./run_05_print_vertices.sh
```

Get number of OGs per connected component.

```bash
while read F;
do
  N=$(cat $F.SFDP.SBM.vertices.tsv | wc -l)
  N=$((N-1))
  echo -e "$F\t$N"
done < component_prefixes.txt \
  | sort -k 2,2nr \
  > hhdb_hhr_reformatted_filtered.tab.90ID.Component.OGcount.tsv
```

