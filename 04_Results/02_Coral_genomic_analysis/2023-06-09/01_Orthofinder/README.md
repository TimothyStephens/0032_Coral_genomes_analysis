# Generate protein orthogroups

Run `Orthofinder` using the predicted proteins from the coral/outgroup genomes/transcriptomes that we have generated/downloaded.

## Setup analysis directory

Link protein sequence files for the coral genomes. Link in a separate directory; tell `orthofinder` to look in here for datasets to run with.

```bash
mkdir files2cluster; cd files2cluster/
ln_loop ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*_genomes/*/00_databases/*.genes.pep.faa
ln_loop ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*_transcriptomes/*/00_databases/*.transcripts.pep.faa
```

Unlink samples that we don't want to analyze.

```bash
# Not part of this analysis
unlink Montipora_sp1_aff_capitata_ULFMv1.genes.pep.faa
```

Setup bash environment.

```bash
conda activate py27
```

## Setup analysis directory

Run `orthofinder` from a lower directory just incase we run into absolute path overflow issues. Can happen although probably not with this run.

Run from: `/scratch/timothy/Orthofinder/`

```bash
mv /scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-02-20/01_Orthofinder/run_Orthofinder.sh .

mkdir files2cluster; cd files2cluster/
ln_loop /scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-02-20/01_Orthofinder/files2cluster/*
```

# Run 1 (Full dataset)

## Run `orthofinder`

Run `orthofinder` using standard (default) parameters, except for `-S diamond_ultra_sens` which runs `diamond` using ultra-sensitive parameters. Had to set `ulimit -n 10000` so `orthofinder` could have more file handles open at the same time.

```bash
./run_Orthofinder-Run1.sh
```

## Copy results

```bash
cp files2cluster/OrthoFinder/Results_Jul14/Orthogroups/Orthogroups.GeneCount.tsv Run1.Orthogroups.GeneCount.tsv
cp files2cluster/OrthoFinder/Results_Jul14/Orthogroups/Orthogroups.tsv Run1.Orthogroups.tsv
cp files2cluster/OrthoFinder/Results_Jul14/Orthogroups/Orthogroups_UnassignedGenes.tsv Run1.Orthogroups_UnassignedGenes.tsv
cp files2cluster/OrthoFinder/Results_Jul14/Phylogenetic_Hierarchical_Orthogroups/N0.tsv Run1.Phylogenetic_Hierarchical_Orthogroups.N0.tsv
cp files2cluster/OrthoFinder/Results_Jul14/Species_Tree/SpeciesTree_rooted_node_labels.txt Run1.SpeciesTree_rooted_node_labels.tre
cp files2cluster/OrthoFinder/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt Run1.SpeciesTree_rooted.tre
dos2unix Run1*
```

## Cleanup species tree

Cleanup names in species tree produced by Orthofinder.

```bash
sed -i -e 's/.genes.pep//g' -e 's/.transcripts.pep//g' \
    Run1.SpeciesTree_rooted_node_labels.tre
sed -i -e 's/.genes.pep//g' -e 's/.transcripts.pep//g' \
    Run1.SpeciesTree_rooted.tre
```

## Cleanup species names in orthogroups `tsv` files

```bash
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run1.Orthogroups.GeneCount.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run1.Orthogroups.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run1.Orthogroups_UnassignedGenes.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.tsv
```

## Get lists of tips in trees

```bash
conda activate r4_env-ape
export PATH="/home/timothy/miniconda3/envs/r4_env-ape/bin:$PATH"
export R_LIBS="/home/timothy/miniconda3/envs/r4_env-ape/lib/R/library"

Rscript ~/GitHub/Utils/Phylo_Tree_Tools/print_tip_labels.R \
    Run1.SpeciesTree_rooted_node_labels.tre | awk '{print $2}' \
  > Run1.SpeciesTree_rooted_node_labels.tip_labels.txt
```



## Reformat OGs

Convert to "long" format (i.e., one gene name per row, each OG can/will occupy multiple rows).

Output format: OG_id [tab] genome_name [tab] seq_id

```bash
# Orthogroups
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=2; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run1.Orthogroups.tsv \
 > Run1.Orthogroups.long.tsv

# Singletons
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=2; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run1.Orthogroups_UnassignedGenes.tsv \
 > Run1.Orthogroups_UnassignedGenes.long.tsv
```

Combine the multi-sequence OGs and singleton OGs long-formatted files together.

```bash
cat Run1.Orthogroups.long.tsv Run1.Orthogroups_UnassignedGenes.long.tsv > Run1.Orthogroups_ALL.long.tsv

# Check no. seqs in long formatted files.
wc -l Run1.Orthogroups_ALL.long.tsv
#4138912 Run1.Orthogroups_ALL.long.tsv
```



## Reformat HOGs

Convert to "long" format (i.e., one gene name per row, each HOG can/will occupy multiple rows).

Output format: HOG_id [tab] genome_name [tab] seq_id

```bash
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=4; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.tsv \
 > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv
```

Infer singleton HOGs.

```bash
N=0
for F in $(cat files2cluster/OrthoFinder/Results_Jul14/WorkingDirectory/SpeciesIDs.txt | awk '{print "files2cluster/"$2}');
do
  P=$(basename $F | sed -e 's/.genes.pep.faa//' -e 's/.transcripts.pep.faa//')
  cat "${F}" \
    | grep '>' \
    | sed -e 's/>//' \
    | ~/scripts/grepf_column.py -v -f <(cut -f3 "Run1.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv") \
    | while read ID; do echo -e "$P\t$ID"; done
done \
  | while read LINE; do printf "SOG%07.f\t%s\n" "$N" "$LINE"; N=$((N+1)); done \
  > "Run1.Phylogenetic_Hierarchical_Orthogroups.N0.singletons.long.tsv"
```

Combine multi-seq + singleton HOGs.

```bash
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv \
    Run1.Phylogenetic_Hierarchical_Orthogroups.N0.singletons.long.tsv \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv

# Check no. seqs in long formatted files.
wc -l Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv
#4138912 Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv
```



## Stats OGs

TOTAL PROTEINS: 4,138,912

- (1) no. genes per group

Output format: OG_name [tab] count_seqs_in_OG

```bash
cat Run1.Orthogroups_ALL.long.tsv \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_OG.txt
```

- (2) no. species per OG

Output format: OG_name [tab] count_species_in_OG

```bash
cat Run1.Orthogroups_ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.species_per_OG.txt
```

- (3) no. genes per species in each OG

Output format: OG_name [tab] species_name [tab] count_seqs_per_species_in_OG

```bash
cat Run1.Orthogroups_ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species_per_OG.txt
```

- (4) no. genes per species (sanity check; should match input pep files)

Output format: species_name [tab] count_seqs_per_species

```bash
cat Run1.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species.txt
```

- (5) no. genes in multi-sequence OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run1.Orthogroups_ALL.long.stats.seqs_per_OG.txt) \
    -i Run1.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species_multiSeqOGs.txt
```

- (6) no. genes in single-sequence OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run1.Orthogroups_ALL.long.stats.seqs_per_OG.txt) \
    -i Run1.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species_singleSeqOGs.txt
```

- (7) no. genes in multi-species OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run1.Orthogroups_ALL.long.stats.species_per_OG.txt) \
    -i Run1.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt
```

- (8) no. genes in single-species OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run1.Orthogroups_ALL.long.stats.species_per_OG.txt) \
    -i Run1.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt
```

- (x) combine per species results files

```bash
echo -e "species_name\tcount_seqs_per_species\tno_seqs_multiSeqOGs\tno_seqs_singleSeqOGs\tno_seqs_multiSpOGs\tno_seqs_singleSpOGs" \
  > Run1.Orthogroups_ALL.long.stats.seqs_per_species.combined.txt
cat Run1.Orthogroups_ALL.long.stats.seqs_per_species.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Orthogroups_ALL.long.stats.seqs_per_species_multiSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Orthogroups_ALL.long.stats.seqs_per_species_singleSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Orthogroups_ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Orthogroups_ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt \
  >> Run1.Orthogroups_ALL.long.stats.seqs_per_species.combined.txt
```



## Stats HOGs

TOTAL PROTEINS: 4,138,912

- (1) no. genes per group

Output format: HOG_name [tab] count_seqs_in_HOG

```bash
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt
```

- (2) no. species per HOG

Output format: HOG_name [tab] count_species_in_HOG

```bash
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt
```

- (3) no. genes per species in each HOG

Output format: HOG_name [tab] species_name [tab] count_seqs_per_species_in_HOG

```bash
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_per_OG.txt
```

- (4) no. genes per species (sanity check; should match input pep files)

Output format: species_name [tab] count_seqs_per_species

```bash
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.txt
```

- (5) no. genes in multi-sequence HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt) \
    -i Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSeqOGs.txt
```

- (6) no. genes in single-sequence HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt) \
    -i Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSeqOGs.txt
```

- (7) no. genes in multi-species HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt) \
    -i Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt
```

- (8) no. genes in single-species HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt) \
    -i Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt
```

- (x) combine per species results files

```bash
echo -e "species_name\tcount_seqs_per_species\tno_seqs_multiSeqOGs\tno_seqs_singleSeqOGs\tno_seqs_multiSpOGs\tno_seqs_singleSpOGs" \
  > Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.combined.txt
cat Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt \
  >> Run1.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.combined.txt
```



# Run 2 (Removed low quality datasets)

Re-run `orthofinder` using standard (default) parameters **BUT** without the really incomplete datasets. 

We define incomplete as any dataset with > 50% missing BUSCO proteins (in the "protein" dataset) using both the metazoa and eukaryota datasets. NOTE: we also removed "Pachycerianthus_borealis_NANAv1" as it almost had >50% missing BUSCOs in both datasets (50.1% and 49.8%) and was initially misplaced in the tree (likely as a result of it being generally incomplete).

Also remove the two datasets which we have duplicates of. Pick the "best", preferencing the genome data were possible.

```bash
#20: Alatina_alata_BNNLv1.genes.pep.faa
#69: Montastraea_faveolata_REEFv1.transcripts.pep.faa
#79: Myxobolus_pendula_ONCAv1.transcripts.pep.faa
#82: Pachycerianthus_borealis_NANAv1.transcripts.pep.faa
#102: Porites_astreoides_REEFv1.transcripts.pep.faa
#113: Renilla_reniformis_FLUSv1.genes.pep.faa
#121: Thelohanellus_kitauei_TJCNv1.genes.pep.faa

## Duplicates
#28: Calvadosia_cruxmelitensis_CRGBv3.2-TRANS.transcripts.pep.faa
#31: Cassiopea_xamachana_T1-Av2-TRANS.transcripts.pep.faa
```

Comment (`#`) out these entries in `files2cluster/OrthoFinder/Results_Jul14/WorkingDirectory/SpeciesIDs.txt` and rerun Orthofinder to recompute orthogroups. NOTE: Had to run orthofinder from `/scratch/timothy/Orthofinder/` as otherwise the absolute paths used internally were too long and caused errors.

```bash
./run_Orthofinder-Run2.sh
```

Link to new results file in `files2cluster/OrthoFinder/`

```bash
ln -s Results_Jul14/WorkingDirectory/OrthoFinder/Results_Jul24/ Results_Jul24
```

## Results

The orthogroups are listed in `files2cluster/OrthoFinder/Results_Jul24/Orthogroups/Orthogroups.tsv`. Orthofinder reccomend usiong the "Phylogenetic Hierarchical Orthogroups" but this approach seems to split almost all (>95%) of the *M. capitata* KBHIv3 proteins into singletons.

## Copy results

```bash
cp files2cluster/OrthoFinder/Results_Jul24/Orthogroups/Orthogroups.GeneCount.tsv Run2.Orthogroups.GeneCount.tsv
cp files2cluster/OrthoFinder/Results_Jul24/Orthogroups/Orthogroups.tsv Run2.Orthogroups.tsv
cp files2cluster/OrthoFinder/Results_Jul24/Orthogroups/Orthogroups_UnassignedGenes.tsv Run2.Orthogroups_UnassignedGenes.tsv
cp files2cluster/OrthoFinder/Results_Jul24/Phylogenetic_Hierarchical_Orthogroups/N0.tsv Run2.Phylogenetic_Hierarchical_Orthogroups.N0.tsv
cp files2cluster/OrthoFinder/Results_Jul24/Species_Tree/SpeciesTree_rooted_node_labels.txt Run2.SpeciesTree_rooted_node_labels.tre
cp files2cluster/OrthoFinder/Results_Jul24/Species_Tree/SpeciesTree_rooted.txt Run2.SpeciesTree_rooted.tre
dos2unix Run2*
```

## Cleanup species tree

Cleanup names in species tree produced by Orthofinder.

```bash
sed -i -e 's/.genes.pep//g' -e 's/.transcripts.pep//g' \
    Run2.SpeciesTree_rooted_node_labels.tre
sed -i -e 's/.genes.pep//g' -e 's/.transcripts.pep//g' \
    Run2.SpeciesTree_rooted.tre
```

## Cleanup species names in orthogroups `tsv` files

```bash
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run2.Orthogroups.GeneCount.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run2.Orthogroups.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run2.Orthogroups_UnassignedGenes.tsv
sed -i -e '1 s/.genes.pep//g' -e '1 s/.transcripts.pep//g' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.tsv
```

## Get lists of tips in trees

```bash
conda activate r4_env-ape
export PATH="/home/timothy/miniconda3/envs/r4_env-ape/bin:$PATH"
export R_LIBS="/home/timothy/miniconda3/envs/r4_env-ape/lib/R/library"

Rscript ~/GitHub/Utils/Phylo_Tree_Tools/print_tip_labels.R \
    Run2.SpeciesTree_rooted_node_labels.tre | awk '{print $2}' \
  > Run2.SpeciesTree_rooted_node_labels.tip_labels.txt
```



## Reformat OGs

Convert to "long" format (i.e., one gene name per row, each OG can/will occupy multiple rows).

Output format: OG_id [tab] genome_name [tab] seq_id

```bash
# Orthogroups
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=2; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run2.Orthogroups.tsv \
 > Run2.Orthogroups.long.tsv

# Singletons
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=2; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run2.Orthogroups_UnassignedGenes.tsv \
 > Run2.Orthogroups_UnassignedGenes.long.tsv
```

Combine the multi-sequence OGs and singleton OGs long-formatted files together.

```bash
cat Run2.Orthogroups.long.tsv Run2.Orthogroups_UnassignedGenes.long.tsv > Run2.Orthogroups_ALL.long.tsv

# Check no. seqs in long formatted files.
wc -l Run2.Orthogroups_ALL.long.tsv
#3778150 Run2.Orthogroups_ALL.long.tsv
```



## Reformat HOGs

Convert to "long" format (i.e., one gene name per row, each HOG can/will occupy multiple rows).

Output format: HOG_id [tab] genome_name [tab] seq_id

```bash
awk -F'\t' '{
    if(NR==1){
      split($0,HEADER,"\t")
  } else { 
    for(i=4; i<=NF; i++){ 
      if($i!="") {
        split($i,a,", "); 
        for(j=1; j<=length(a); j++){
          print $1"\t"HEADER[i]"\t"a[j]
        }
      }
    }
  }
}' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.tsv \
 > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv
```

Infer singleton HOGs.

```bash
N=0
for F in $(grep -v '#' files2cluster/OrthoFinder/Results_Jul14/WorkingDirectory/SpeciesIDs.txt | awk '{print "files2cluster/"$2}');
do
  P=$(basename $F | sed -e 's/.genes.pep.faa//' -e 's/.transcripts.pep.faa//')
  cat "${F}" \
    | grep '>' \
    | sed -e 's/>//' \
    | ~/scripts/grepf_column.py -v -f <(cut -f3 "Run2.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv") \
    | while read ID; do echo -e "$P\t$ID"; done
done \
  | while read LINE; do printf "SOG%07.f\t%s\n" "$N" "$LINE"; N=$((N+1)); done \
  > "Run2.Phylogenetic_Hierarchical_Orthogroups.N0.singletons.long.tsv"
```

Combine multi-seq + singleton HOGs.

```bash
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.long.tsv \
    Run2.Phylogenetic_Hierarchical_Orthogroups.N0.singletons.long.tsv \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv

# Check no. seqs in long formatted files.
wc -l Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv
#3778150 Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv
```



## Stats OGs

TOTAL PROTEINS: 3,778,150

- (1) no. genes per group

Output format: OG_name [tab] count_seqs_in_OG

```bash
cat Run2.Orthogroups_ALL.long.tsv \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_OG.txt
```

- (2) no. species per OG

Output format: OG_name [tab] count_species_in_OG

```bash
cat Run2.Orthogroups_ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.species_per_OG.txt
```

- (3) no. genes per species in each OG

Output format: OG_name [tab] species_name [tab] count_seqs_per_species_in_OG

```bash
cat Run2.Orthogroups_ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species_per_OG.txt
```

- (4) no. genes per species (sanity check; should match input pep files)

Output format: species_name [tab] count_seqs_per_species

```bash
cat Run2.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species.txt
```

- (5) no. genes in multi-sequence OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run2.Orthogroups_ALL.long.stats.seqs_per_OG.txt) \
    -i Run2.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species_multiSeqOGs.txt
```

- (6) no. genes in single-sequence OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run2.Orthogroups_ALL.long.stats.seqs_per_OG.txt) \
    -i Run2.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species_singleSeqOGs.txt
```

- (7) no. genes in multi-species OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run2.Orthogroups_ALL.long.stats.species_per_OG.txt) \
    -i Run2.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt
```

- (8) no. genes in single-species OGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run2.Orthogroups_ALL.long.stats.species_per_OG.txt) \
    -i Run2.Orthogroups_ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt
```

- (x) combine per species results files

```bash
echo -e "species_name\tcount_seqs_per_species\tno_seqs_multiSeqOGs\tno_seqs_singleSeqOGs\tno_seqs_multiSpOGs\tno_seqs_singleSpOGs" \
  > Run2.Orthogroups_ALL.long.stats.seqs_per_species.combined.txt
cat Run2.Orthogroups_ALL.long.stats.seqs_per_species.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Orthogroups_ALL.long.stats.seqs_per_species_multiSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Orthogroups_ALL.long.stats.seqs_per_species_singleSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Orthogroups_ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Orthogroups_ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt \
  >> Run2.Orthogroups_ALL.long.stats.seqs_per_species.combined.txt
```



## Stats HOGs

TOTAL PROTEINS: 4,138,912

- (1) no. genes per group

Output format: HOG_name [tab] count_seqs_in_HOG

```bash
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt
```

- (2) no. species per HOG

Output format: HOG_name [tab] count_species_in_HOG

```bash
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt
```

- (3) no. genes per species in each HOG

Output format: HOG_name [tab] species_name [tab] count_seqs_per_species_in_HOG

```bash
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f1,2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_per_OG.txt
```

- (4) no. genes per species (sanity check; should match input pep files)

Output format: species_name [tab] count_seqs_per_species

```bash
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.txt
```

- (5) no. genes in multi-sequence HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt) \
    -i Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSeqOGs.txt
```

- (6) no. genes in single-sequence HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_OG.txt) \
    -i Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSeqOGs.txt
```

- (7) no. genes in multi-species HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2>1{print $1}' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt) \
    -i Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt
```

- (8) no. genes in single-species HOGs per species

```bash
~/scripts/grepf_column.py \
    -f <(awk -F'\t' '$2==1{print $1}' Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.species_per_OG.txt) \
    -i Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.tsv \
  | cut -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$1}' \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt
```

- (x) combine per species results files

```bash
echo -e "species_name\tcount_seqs_per_species\tno_seqs_multiSeqOGs\tno_seqs_singleSeqOGs\tno_seqs_multiSpOGs\tno_seqs_singleSpOGs" \
  > Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.combined.txt
cat Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSeqOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_multiSpeciesOGs.txt \
  | ~/scripts/add_value_to_table.py -a Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species_singleSpeciesOGs.txt \
  >> Run2.Phylogenetic_Hierarchical_Orthogroups.N0.ALL.long.stats.seqs_per_species.combined.txt
```



# Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/01_Orthofinder"

rsync -zarv --prune-empty-dirs --relative \
  --exclude="*/files2cluster/*" \
  --include="*/" \
  --include="*.sh*" --include="Run*" \
  --include="md5sum_list*" \
  --exclude="*" \
  ${WD}/./ \
  . --dry-run
```

