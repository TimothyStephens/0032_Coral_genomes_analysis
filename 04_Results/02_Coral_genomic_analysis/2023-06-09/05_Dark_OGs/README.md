# Identify dark OGs

Perform analysis designed to identify "dark" and "light" genes and orthogroups.



## Classify `nr` sequences

Working directory: `05_Dark_OGs/`

Classify each sequence in the `nr` database that will be used for dark gene analysis.

Export info about each sequence in the `nr` database and classify them as in/out taxon target taxon group AND annotated/unannotated. 

```bash
./run_01_classify_nr_seqs.sh
```

Results file: `nr_seqs_classification.tsv.gz`



---

Create SQLite database with nr classification information. 

```bash
./run_02_create_database.sh
```

Results file: `nr_seqs_classification.tsv.gz.sqlite3`



## Annotate target proteins

Annotate each sequence in protein set with best information from `nr`.

Link to `fasta` protein seqs and `diamond blastp` against `nr` results for each genome/transcriptome dataset used.

```bash
mkdir diamond_results; cd diamond_results/
ln_loop ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/*/10_functional_annotation/*.pep.faa.diamond_blastp_nr.outfmt6_short.gz
ln_loop ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/*/10_functional_annotation/*.pep.faa
```

Run annotation workflow.

- Will assign top hits:
    - That aren't functionally ambiguous (i.e., has known function; not "known" or "hypothetical" etc.)
    - That are from outside our target group (i.e., not a coral). 
- We do this since "function" of protein is a separate question to its "taxonomic distribution"

Each dataset will have its own batch submission script (to make rerunning, or running now datasets easier).

```bash
for F in diamond_results/*.faa;
do
  P=$(basename $F)
  sed -e "s/<<<QUERY>>>/${P}/" run_03_add_classification.sh > run_03_add_classification.sh.${P}
done
chmod +x run_03_add_classification.sh*
```

**Run each job separately**

Combine results files together.

```bash
cat diamond_results/*.gz.top > diamond_results.nr.top_hits
```



---

Add `nr` information to each rebuilt OG.

```bash
ln -s ../04_Process_OGs/OG_analysis.combined.results.tsv
ln -s ../04_Process_OGs/Run2.Orthogroups_ALL.long.tsv
```

Run script to add `nr` annotations to each seq_id in each cluster.

```bash
./run_04_add_classification_to_clusters.sh
```

Check that there are no sequences with "Missing" info - this would indicate that something went wrong with the `nr` classification step.

```bash
grep 'Missing' Run2.Orthogroups_ALL.long.nr.top_hits.tsv
```

Classify OGs.

```bash
cat Run2.Orthogroups_ALL.long.nr.top_hits.tsv \
  | awk -F'\t' '{ 
      if (!seen[$1]++){
        C[$1]["Annotated"]="No";
        C[$1]["OutsideGroup"]="No";
      }
      if($4!="NA") {C[$1]["Annotated"]="Yes";}
      if($9!="NA") {C[$1]["OutsideGroup"]="Yes";}
    } END {
      for(i in C){
        if(C[i]["Annotated"]=="Yes" && C[i]["OutsideGroup"]=="Yes"){print i"\tLight-Shared"};
        if(C[i]["Annotated"]=="Yes" && C[i]["OutsideGroup"]=="No" ){print i"\tLight-Restricted"};
        if(C[i]["Annotated"]=="No"  && C[i]["OutsideGroup"]=="Yes"){print i"\tDark-Shared"};
        if(C[i]["Annotated"]=="No"  && C[i]["OutsideGroup"]=="No" ){print i"\tDark-Restricted"};
      }
    }' \
  > Run2.Orthogroups_ALL.long.nr.top_hits.type.tsv
```

Count the number of each type of cluster.

```bash
cut -f2 Run2.Orthogroups_ALL.long.nr.top_hits.type.tsv | sort | uniq -c
```

> 213214 Dark-Restricted
>      11422 Dark-Shared
>        9708 Light-Restricted
>     88889 Light-Shared



Link to other strata info for each OG.

```bash
python scripts/add_value_to_table_SQLite3.py \
  -a <(awk 'BEGIN{print "cluster_id\tDesignation"} {print}' \
         Run2.Orthogroups_ALL.long.nr.top_hits.type.tsv) \
  -d "MissingDesignation" \
  -i OG_analysis.combined.results.tsv \
  -o OG_analysis.combined.results.type.tsv
```





