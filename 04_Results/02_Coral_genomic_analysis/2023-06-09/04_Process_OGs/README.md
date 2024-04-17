# Refine and analyze orthogroups produced by `orthofinder`

Take the orthogroups produced by `orthofinder` and:

- Orthogroup taxonomic distribution - Classify at which taxonomic level the orthogroup is specific to - i.e., just species X or shared by all.

- Rate of evolution - Check if the rate of evolution (using `abSENSE`) of the genes in a group supports the absence/lack of detection in other species.

## Prepare orthogroup information

Working directory: `04_Process_OGs/`

Link to orthogroup results to process.

```bash
ln -s ../01_Orthofinder/Run2.Orthogroups_ALL.long.tsv
ln -s ../01_Orthofinder/Run2.Orthogroups_ALL.long.stats.seqs_per_OG.txt
ln -s ../01_Orthofinder/Run2.SpeciesTree_rooted.tre
```

Link to taxi info for each of the species in our analysis.

```bash
python ../../../../02_Scripts/grepf_column.py \
  --keep_header \
  -f ../01_Orthofinder/Run2.SpeciesTree_rooted_node_labels.tip_labels.txt \
  -i ../03_Lineage_Info/species_taxids.lineage.tsv \
  -o species_taxids.lineage.tsv
```

Link to seq files used for `orthofinder` in the `seq_files/` sub directory. Unlink files from datasets that we didn't use in the final clustering.

```bash
mkdir seq_files; cd seq_files/

ln_loop ../../01_Orthofinder/files2cluster/*.pep.faa

## Low quality
unlink Alatina_alata_BNNLv1.genes.pep.faa
unlink Montastraea_faveolata_REEFv1.transcripts.pep.faa
unlink Myxobolus_pendula_ONCAv1.transcripts.pep.faa
unlink Pachycerianthus_borealis_NANAv1.transcripts.pep.faa
unlink Porites_astreoides_REEFv1.transcripts.pep.faa
unlink Renilla_reniformis_FLUSv1.genes.pep.faa
unlink Thelohanellus_kitauei_TJCNv1.genes.pep.faa

## Duplicates
unlink Calvadosia_cruxmelitensis_CRGBv3.2-TRANS.transcripts.pep.faa
unlink Cassiopea_xamachana_T1-Av2-TRANS.transcripts.pep.faa

cd ../
```



## Refine orthogroups

Split list of OGs into 100 separate parts so we can process everything in a highly parallel fashion. Also get list of parts.

```bash
mkdir -p IDs_split

cat Run2.Orthogroups_ALL.long.stats.seqs_per_OG.txt \
  | sort -k2,2nr \
  | awk '{print $1}' \
  | split -a 3 --numeric-suffixes=1 -n r/100 - IDs_split/

ls -1 --color=none IDs_split/* > files2run.txt
```

Create distance file so we know how far each of the species are from each other in the tree (used for de novo gene detection).

```bash
conda activate r4_env-ape
export R_LIBS="/home/timothy/miniconda3/envs/r4_env-ape/lib/R/library"
R
```

Use `R` to get distance matrix.

```R
library(ape)
tree<-read.tree(file = "Run2.SpeciesTree_rooted.tre")
distances<-cophenetic.phylo(x = tree)
write.table(distances,"pairwise_distances.matrix.tsv", row.names = TRUE, sep = "\t")
```

Manually fix header by adding "Name" to first column.

```bash
sed -e '1 s/^/Name\t/' \
    -e 's/"//g' \
    -i "pairwise_distances.matrix.tsv"

awk -F'\t' '{
  if(NR==1){
    for(i=2; i<=NF; i++){
      H[i]=$i
    }
  } else {
    for(i=2; i<=NF; i++){
      print $1"\t"H[i]"\t"$i
    }
  }
}' "pairwise_distances.matrix.tsv" \
 > "pairwise_distances.flat.tsv"
```



Install required packages to run custom workflow.

```bash
P="/home/timothy/miniconda3/envs/OG_analysis"
mamba create  -p $P python=2.7
mamba install -p $P -c anaconda scipy
mamba install -p $P -c conda-forge matplotlib
mamba install -p $P -c conda-forge dill

# copy install info and setup on other server
conda list --export > package-list.txt
#mamba create --name OG_analysis --file package-list.txt
```



## Process each HOG

Copy files to Amarel server to run.

```bash
scp -r timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/04_Process_OGs .
```

Run each split set of OG IDs.

```bash
sbatch run_01_OG_analysis.sh
```

Check that all parts have finished by checking for the output files.

```bash
while read ID;
do
  if [[ ! -s "OG_analysis/${ID}.results.tsv" ]];
  then
    echo $ID
  fi
done < <(cat IDs_split/*) > not_finished_1.txt
```

Looks like 7 parts failed becuase they were single amino acid sequences that caused `diamond` to think that they were nucleotide and not proteins.

```bash
grep 'Parallel finished with exit status:' OG_analysis.slurm_out.* | grep -v 'status: 0'
```

> OG_analysis.slurm_out.30475287-00.12:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475306-00.31:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475310-00.35:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475314-00.39:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475325-00.50:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475337-00.62:Parallel finished with exit status: 1
> OG_analysis.slurm_out.30475359-00.84:Parallel finished with exit status: 1

Need to rerun these parts using a version of the script that uses a different version of `diamond` for database construction. This version allows us to surpress the warning. 

```bash
while read ID;
do
  ./OG_analysis_multiVersionDIAMOND.sh --id "$ID" --long Run2.Orthogroups_ALL.long.tsv --seqs seq_files --dist pairwise_distances.flat.tsv --lineage species_taxids.lineage.tsv > logs/"$ID".log 2>&1
done < not_finished_1.txt
```

Check that the parts we reran have finished correctly.

```bash
while read ID;
do
  if [[ ! -s "OG_analysis/${ID}.results.tsv" ]];
  then
    echo $ID
  fi
done < not_finished_1.txt > not_finished_2.txt
```



## Check logs

Check that all parts have finished by filtering each of the part logs.

```bash
./run_02_check_logs.sh

# Check for logs with a specific phrase
while read ID;
do
  grep -l 'RuntimeWarning:' "logs/${ID}.log"
done < <(cat IDs_split/*)
```







## Combine results

Combine and parse the results of the reclustering analysis. Also, compress and zip results files once we are dont combining results.

```bash
./run_03_merge_parts.sh
```

Check that every cluster the we had in the input are in the output and that we have the correct number of columns.

```bash
## Check that every row has the same number of columns
awk -F'\t' '{print NF}' OG_analysis.combined.results.tsv | uniq -c
#323234 27

## Check all inout IDs are represented in the output
awk -F'\t' 'NR>1{print $1}' OG_analysis.combined.results.tsv | awk -F'.' '{ if(NF==3){print $1"."$2}else{print $1} }' | sort | uniq | md5sum
#496cb9f27db6f698aacb7040e805a192  -
cat IDs_split/* | sort | md5sum
#496cb9f27db6f698aacb7040e805a192  -
```



## Compress results

Compress results so that we can copy it back to the Coral server.

```bash
./run_04_compress_results.sh
```



