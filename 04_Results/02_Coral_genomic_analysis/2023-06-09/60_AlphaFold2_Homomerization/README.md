# Assess quaternary structure of representative proteins

## Setup analysis directory

Setup bash environment.

```bash
conda activate main
```

## Select proteins to analyze

Select which representative proteins have surficient structural quality to bother rerunning in this analysis.

```bash
cat ../08_Selected_OG_Analysis/combined_results/combined.stats.OGs.tsv \
  | awk -F'\t' '{print $1"\t"$13"\t"$15}' \
  | awk -F'\t' 'NR==1 || ($2+$3)>50' \
  | awk -F'\t' 'NR!=1{print $1}' \
  > selected_OGs.txt
```

Link to represenative sequences.

```bash
while read F;
do
  ln -s ../08_Selected_OG_Analysis/Orthogroups_analysis/$F.representative.fa
done < selected_OGs.txt
```



## Analysis 2

Analysis

```bash
/opt/conda/envs/ColabFold_v1.5.5-rev1/bin/colabfold_search --mmseqs /opt/conda/envs/ColabFold_v1.5.5-rev1/bin/mmseqs OG0007517.representative.fa /scratch/timothy/20240220 msas
/opt/conda/envs/ColabFold_v1.5.5-rev1/bin/colabfold_batch msas predictions
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

