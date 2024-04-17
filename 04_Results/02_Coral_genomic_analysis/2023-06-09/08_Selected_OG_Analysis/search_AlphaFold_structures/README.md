# Search AlphaFold2 structures

Search the AlphaFold2 strctures that we generated against various databases of 3D structures using `fold-seek`.



## Setup analysis directory

Setup bash environment.

```bash
conda activate main
```



## Combine Coral Structurs

Combine AlphaFold2 structures together before we run them through `fold-seek`.

```bash
cat ../Orthogroups_analysis/*.representative.fa.alphafold.pdb.gz > Orthogroups_analysis.alphafold.pdb.gz
```



## Analysis 2

Analysis

```bash
analysis
```

Analysis

## Results

### Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.py" --exclude="*" \
 ${WD}/./dir1 \
 ${WD}/./dir2 \
 . --dry-run
```

### Visulize results

etc. 

