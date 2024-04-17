# Run 3D Structure Analysis

Predict 3D structure using AlphaFold.







## Setup `R` environment

```bash
mamba create -n plotCoralGenes_0.1.4 -c bioconda -c conda-forge r-base r-essentials r-devtools r-biocmanager

mamba activate plotCoralGenes_0.1.4
mamba install -c bioconda -c conda-forge bioconductor-Biostrings r-RColorBrewer r-ape r-cowplot r-dplyr bioconductor-ggmsa r-ggplot2 r-ggpubr bioconductor-ggtree r-gridExtra r-janitor r-knitr r-patchwork r-phytools r-png r-reshape2 r-shiny r-tibble r-tidyverse

~/miniforge3/envs/plotCoralGenes_0.1.4/bin/R
devtools::install_github("swsoyee/r3dmol")
install.packages("plotCoralGenes_0.1.4.tar.gz", repos = NULL, type="source")
```







## Format "Dark" OGs

Extract the "Dark" OGs from each set of groups that we have filtered and identified.

```bash
cat ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Shared_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Shared_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Shared_50.long.tsv \
  | awk -F'\t' '$1!="orthogroup_id" {print $1}' \
  | sort | uniq \
  > Orthogroups.txt

cat ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Shared_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Shared_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Restricted_50.long.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Shared_50.long.tsv \
  | awk -F'\t' 'NR==1 || $1!="orthogroup_id"' \
  > Orthogroups.long.tsv

cat ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Restricted_50.classification.tsv \
    ../07_Plot_Dark_OGs/selected_Cnidaria_Dark-Shared_50.classification.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Restricted_50.classification.tsv \
    ../07_Plot_Dark_OGs/selected_Hexacorallia_Dark-Shared_50.classification.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Restricted_50.classification.tsv \
    ../07_Plot_Dark_OGs/selected_Scleractinia_Dark-Shared_50.classification.tsv \
  | awk -F'\t' 'NR==1 || $1!="orthogroup_id"' \
  > Orthogroups.classification.tsv
```



Make analysis directory.

```bash
mkdir Orthogroups.analysis
```



Get all sequences so that it is quicker to extract them for each orthogroup.

```bash
mkdir all_seqs_data; cd all_seqs_data/
cat ../../05_Dark_OGs/diamond_results/*.faa | gzip -c > all_combined.pep.faa.gz
```

Get all eggNOG + InterProScan annotations so that it is quicker to extract them for each orthogroup.

```bash
cat ../../01_Orthofinder/Run2.SpeciesTree_rooted_node_labels.tip_labels.txt \
  | while read PREFIX;
    do
      cat ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/${PREFIX}/10_functional_annotation/${PREFIX}*.pep.faa.InterProScan.gff3
    done \
  | awk -F'\t' '$1!~"#" && NF==9' \
  | gzip -c \
  > all_combined.pep.faa.InterProScan.gff3.gz

cat ../../01_Orthofinder/Run2.SpeciesTree_rooted_node_labels.tip_labels.txt \
  | while read PREFIX;
    do
      cat ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/${PREFIX}/10_functional_annotation/${PREFIX}*.pep.faa.emapper.annotations
    done \
  | awk 'NR==1 || $1!~"#"' \
  | sed -e '1 s/#//' \
  | gzip -c \
  > all_combined.pep.faa.emapper.tsv.gz
```







## 01 - Extract Protein Sequences

Extract protein sequences for genes in each orthogroup.

```bash
./run_01_extract_sequences.sh
```







## 02 - Run sequence alignment

Use `mafft-linsi` to align the sequences within an orthogroup.

```bash
./run_02_align_sequences.sh
```

Download results.

```bash
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Gene_Structure/./*/*.aln . --dry-run
```







## 03 - Run AlphaFold

Run AlphaFold on each representative sequence.

```bash
split -a 1 --numeric-suffixes=1 -n r/4 Orthogroups.txt Orthogroups.txt.
```

Download results.

```bash
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Gene_Structure/./*/*.alphafold.pdb.gz . --dry-run
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Gene_Structure/./*/*.alphafold.pLDDT_summary.tsv . --dry-run
```







## 04 - Build trees

Biuld trees from `mafft` alignments using `iqtree`

```bash
./run_04_build_tree.sh
```

Download results.

```bash
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Gene_Structure/./*/*.contree . --dry-run
```







## 05 - Extract additonal functional annotations

Extract additional functional annotations (eggNOG-Mapper or InterProScan).

```bash
./run_05_extract_annotations.sh
```







## 06 - BLASTP nr

BLASTP dark orthogroup protein sequences against NCBI's nr database.

```bash
./run_06_blastp_nr.sh.log
```







## 07 - Visualize gene models

Visualize *M. capitata* gene models to help validate dark gene structure using IGV.

### Extract *M. capitata* gene models

```bash
./Orthogroups_analysis_07_extract_Mcap_gene_models.sh
```



### Download files for IGV analysis.

RNA-seq mapping.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/transcriptomic/PRJNA694677-Illumina-HISAT2_StringTie2/MC-289_All_combined.HISAT2_RNAseq_mapping.coordSorted.bam* .
```

All gene modes.

```bash
rsync -avzPL timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/00_databases/Montipora_capitata_KBHIv3.genes.gff3 .
```

Reference genome.

```bash
rsync -avzPL timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/00_databases/Montipora_capitata_KBHIv3.assembly.fasta* .
```



Extrated gene modesl per set.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Selected_OG_Analysis/Orthogroups_analysis/*.genes.gff3 . --dry-run
```



```bash
for GFF3 in *.genes.gff3;
do
  LOCUS=$(awk -F'\t' '$3=="transcript"{E=int(($5-$4)*0.20); print $1":"$4-E"-"$5+E}' $GFF3)
  STRAND=$(awk -F'\t' '$3=="transcript"{ if($7=="+"){print "POSITIVE"}else{print "NEGATIVE"} }' $GFF3)
  sed -e "s/<<<GFF3>>>/$GFF3/g" -e "s/<<<LOCUS>>>/${LOCUS}/g" -e "s/<<<STRAND>>>/${STRAND}/g" \
      ../igv_data/template.igv_session.xml \
    > ${GFF3%*.genes.gff3}.igv_session.xml
done
```



Setup an IGV batch file (https://software.broadinstitute.org/software/igv/batch) to auto snapshot each gene. Create on the system where you will run IGV (so that we can set the absolute path to where we want the snapshots copied to).

```bash
echo "snapshotDirectory '${PWD}'" > batch_igv_sessions.txt
for IGV in *.igv_session.xml;
do
  echo "load ${IGV}"
  echo "snapshot ${IGV%*.xml}.svg"
  echo "snapshot ${IGV%*.xml}.png"
done >> batch_igv_sessions.txt
```



Upload to server.

```bash
rsync -avzP *.igv_session.* timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Selected_OG_Analysis/Orthogroups_analysis/ --dry-run
```







## 08 - Plot orthogroup phylogeny

Plot phylogeny of orthogroup sequences using tree build in `04 - Build trees`.

```bash
./run_08_Plot_OG_Phylogeny.sh
```







## 10 - Plot Expression Patterns *M. capitata*

Plot accumulation values of *M. capitata* genes.

```bash
./run_10_Plot_Expression_Patterns_Mcapitata.sh
```







## 11 - Plot Expression Patterns *P. acuta*

Plot accumulation values of *P. acuta* genes.

```bash
./run_11_Plot_Expression_Patterns_Pacuta.sh
```







## 20 - Plot *Stylophora pistillata* adult scRNA-seq

Plot *Stylophora pistillata* adult scRNA-seq gene expression data.

```bash
./run_20_Plot_scRNA_Stylophora_pistillata_adult.sh
```









## 21 - Plot *Stylophora pistillata* **larva** scRNA-seq

Plot *Stylophora pistillata* adult scRNA-seq gene expression data.

```bash
./run_21_Plot_scRNA_Stylophora_pistillata_larva.sh
```









## 22 - Plot *Stylophora pistillata* polyp scRNA-seq

Plot *Stylophora pistillata* adult scRNA-seq gene expression data.

```bash
./run_22_Plot_scRNA_Stylophora_pistillata_polyp.sh
```









## 23 - Plot *Hydra vulgaris scRNA-seq

Plot *Hydra vulgaris* scRNA-seq gene expression data.

```bash
./run_23_Plot_scRNA_Hydra_vulgaris.sh
```







## 24 - Plot *Nematostella vectensis scRNA-seq

Plot *Nematostella vectensis* scRNA-seq gene expression data.

```bash
./run_24_Plot_scRNA_Nematostella_vectensis.sh
```







## 25 - Plot *Xenia* sp. scRNA-seq

Plot *Xenia* sp. scRNA-seq gene expression data.

```bash
./run_25_Plot_scRNA_Xenia_sp.sh
```







## 99 - Combine results

Combine all results from previous analysis together, both as an interactive HTML file, and as a table.



### Extract stats for seqs in OG

Extract no. CDS and protein lengths for each sequence in selected OGs.

```bash
mkdir -p combined_results
mkdir -p all_seqs_data; cd all_seqs_data/
```

Get No. CDS.

```bash
cat ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/*/00_databases/*.gff3 \
  | awk -F'\t' '$3=="CDS"{print $9}' \
  | sed -e 's/.*Parent=//' -e 's//;.*/' \
  | sort | uniq -c \
  | awk '{print $2"\t"$1}' \
  > all_combined_no_CDS.tsv
```

Get protein lengths.

```bash
cat ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/*/*/00*/*.pep.faa \
  | ~/programs/seqkit_v2.3.1/seqkit fx2tab \
  | awk -F'\t' '{print $1"\t"length($2)}' \
  > all_combined_PEP_length.tsv 
```

Combine CDS and PEP results for all genes in selected OGs.

```bash
zcat ../../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz \
  | python ../../../../../02_Scripts/add_value_to_table_SQLite3.py \
      -c 3 -d "NA" \
      -a <(awk -F'\t' 'BEGIN{print "sequence_id\tnum_CDS"}{print}' all_combined_no_CDS.tsv) \
  | python ../../../../../02_Scripts/add_value_to_table_SQLite3.py \
      -c 3 \
      -a <(awk -F'\t' 'BEGIN{print "sequence_id\tPEP_length"}{print}' all_combined_PEP_length.tsv) \
  > ../combined_results/combined.stats.ALL_seqs.info.tsv
```



Extract stats for each OG.

```bash
cd Orthogroups_analysis/

while read OG;
do
  cat "${OG}.long.tsv" \
    | awk -F'\t' 'NR>1{print $3}' \
    | python "../../../../../02_Scripts/add_value_to_table_SQLite3.py" -d 0 -a "../all_seqs_data/all_combined_no_CDS.tsv" \
    > "${OG}.no_CDS.tsv"
done < Orthogroups.txt

while read OG;
do
  cat "${OG}.long.tsv" \
    | awk -F'\t' 'NR>1{print $3}' \
    | python "../../../../../02_Scripts/add_value_to_table_SQLite3.py" -d 0 -a "../all_seqs_data/all_combined_PEP_length.tsv" \
    > "${OG}.PEP_length.tsv"
done < Orthogroups.txt
```



Extract "classification" results.

```bash
while read OG;
do
  cat "../Orthogroups.classification.tsv" \
    | python "../../../../../02_Scripts/grepf_column.py" --keep_header -f <(echo $OG) \
    > "${OG}.classification.tsv"
done < Orthogroups.txt
```







Run scripts to combine results files together into a single seq of tables and HTML document. 

```bash
./Orthogroups_analysis_99_Combine_Results.sh
```







### Combine "selected" datasets

Aggrigate combind `*.html` files together. 

```bash
cd combined_results/; mkdir -p Orthogroups_analysis
rsync -avzP --exclude="*representative.fa*" --exclude="run*" \
  ../Orthogroups_analysis/*.html \
  Orthogroups_analysis/
```



Copy the `stats.OGs.tsv` and `stats.seqs.tsv` files over to  `combined_results/`.

```bash
cp ../Orthogroups_analysis/stats.OGs.tsv  combined.stats.OGs.tsv
cp ../Orthogroups_analysis/stats.seqs.tsv combined.stats.seqs.tsv
```



Download results

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Selected_OG_Analysis/combined_results .
```













## 99 - `GREP` extra results for manuscript



### Combine additional annotations.

**Egg-NOG Mapper**

```bash
for F in Orthogroups_analysis/*.emapper.annotations;
do
  P=$(basename ${F%*.emapper.annotations}); 
  awk -vP=$P '{ if(NR==1){print "orthogroup_id\t"$0}else{print P"\t"$0} }' $F; 
done \
  | awk 'NR==1 || $1!="orthogroup_id"' \
  > combined_results/combined.annotations.emapper.tsv
```

> 1	orthogroup_id
> 2	query
> 3	seed_ortholog
> 4	evalue
> 5	score
> 6	eggNOG_OGs
> 7	max_annot_lvl
> 8	COG_category
> 9	Description
> 10	Preferred_name
> 11	GOs
> 12	EC
> 13	KEGG_ko
> 14	KEGG_Pathway
> 15	KEGG_Module
> 16	KEGG_Reaction
> 17	KEGG_rclass
> 18	BRITE
> 19	KEGG_TC
> 20	CAZy
> 21	BiGG_Reaction
> 22	PFAMs

Get just the rows without blank "Description", "GOs", or "KEGG_ko" columns (i.e., has some described function).

```bash
cat combined_results/combined.annotations.emapper.tsv \
  | awk -F'\t' '$9!="-" || $11!="-" || $13!="-"' \
  > combined_results/combined.annotations.emapper.withDescFunc.tsv
```

Group and count the nunber of unique descriptions.

```bash
cat combined_results/combined.annotations.emapper.withDescFunc.tsv \
  | awk -F'\t' 'NR>1{print $1"\t"$10"\t"$9}' \
  | sort | uniq \
  | awk 'BEGIN{print "OG\tGene_ID\tDescription"}{print}' \
  > combined_results/combined.annotations.emapper.withDescFunc.summary.tsv
```



**InterProScan**

```bash
for F in Orthogroups_analysis/*.InterProScan.gff3;
do
  P=$(basename ${F%*.InterProScan.gff3});
  awk -vD=$D -vP=$P '{print P"\t"$0}' $F;
done \
  | awk -F'\t' '$4!="polypeptide"' \
  | awk 'BEGIN{print "orthogroup_id\tseqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes"}{print}' \
  > combined_results/combined.annotations.InterProScan.tsv
```

Get just "functional" annotation.

```bash
cat combined_results/combined.annotations.InterProScan.tsv \
  | awk -F'\t' 'NR==1 || ($3!="MobiDBLite" && $3!="Coils" && $3!="Phobius" && $3!="TMHMM")' \
  > combined_results/combined.annotations.InterProScan.functional.tsv
```

Get just "structural" annotations ["MobiDBLite", "Coils, "Phobius", and "TMHMM"].

```bash
cat combined_results/combined.annotations.InterProScan.tsv \
  | awk -F'\t' 'NR==1 || ($3=="MobiDBLite" || $3=="Coils" || $3=="Phobius" || $3=="TMHMM")' \
  > combined_results/combined.annotations.InterProScan.structural.tsv
```

Get just the transmembrane regions as predicted by TMHMM.

```bash
awk -F'\t' 'NR==1 || $3=="TMHMM"' \
    combined_results/combined.annotations.InterProScan.structural.tsv \
  > combined_results/combined.annotations.InterProScan.structural.TMHMM.tsv
```

Get InterProScan annotation coordinates and overlap the "functional" annotations with "transmambrane" regions. Get bed regions.

```bash
export PATH="/home/timothy/programs/bedtools-2.29.2/bin:$PATH"

cat combined_results/combined.annotations.InterProScan.functional.tsv \
  | awk -F'\t' 'NR>1{print $1"-"$2"\t"$5-1"\t"$6"\t"$10}' \
  | sed -e 's/[^\t]*signature_desc=\([^;]*\);.*Name=\([^;]*\);*.*/\2:\1/' -e 's/[^\t]*Name=\([^;]*\);*.*/\1/' \
  | bedtools sort \
  | bedtools merge -c 4 -o distinct \
  > combined_results/combined.annotations.InterProScan.functional.bed

cat combined_results/combined.annotations.InterProScan.structural.TMHMM.tsv \
  | awk -F'\t' 'NR>1{print $1"-"$2"\t"$5-1"\t"$6}' \
  | bedtools sort \
  | bedtools merge\
  > combined_results/combined.annotations.InterProScan.structural.TMHMM.bed
```

Get descriptions and IDs of features that are >30bp after removing regions which overlap with TMHMM predictions (i.e., functional regins which are not probable transmembrane regions).

```bash
bedtools subtract \
    -a combined_results/combined.annotations.InterProScan.functional.bed \
    -b combined_results/combined.annotations.InterProScan.structural.TMHMM.bed \
  | awk -F'\t' '($3-$2)>30' \
  | awk -F'\t' '{split($1,a,"-"); split($4,b,","); for(i=1; i<=length(b); i++){print a[1]"\t"b[i]} }' \
  | sort | uniq \
  | awk -F'\t' 'BEGIN{print "OG\tDescription"}{print}' \
  > combined_results/combined.annotations.InterProScan.functional.NoTHMHHoverlap.tsv
```



**BLASTP**

Extract BLASTP results, excluding taxids of the species in our study set.

```bash
for F in Orthogroups_analysis/*.fa.blastp_nr.outfmt6.gz;
do
  P=$(basename ${F%*.fa.blastp_nr.outfmt6.gz}); 
  python ../../../../02_Scripts/grepf_column.py -v -c 16 \
    -i <(zcat $F | awk 'BEGIN{OFS=FS="\t"}{split($16,a,";"); $16=a[1]; print}') \
    -f ../03_Lineage_Info/target_taxids.txt \
  | awk -F'\t' '$11<1e-5' \
  | awk -vP=$P '{print P"\t"$0}'
done \
  | awk 'BEGIN{print "orthogroup_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tsalltitles"}{print}' \
  > combined_results/combined.annotations.blastp_nr.tsv
```



See which "types" of dark orthogroups have blastp hits.

```bash
cat combined_results/combined.annotations.blastp_nr.tsv | cut -f1 | uniq | python ~/scripts/add_value_to_table_SQLite3.py -a <(zcat ../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}')
```







### AlphaFold2 search against structure databases

Get hits to structures from non-target species that have <0.01 e-values (Cutoff chosen based on https://github.com/steineggerlab/foldseek/issues/167).

```bash
cat search_AlphaFold_structures/foldseek.PDB.outfmt6 \
  | python ../../../../02_Scripts/grepf_column.py -v -c 16 \
      -f ../03_Lineage_Info/target_taxids.txt \
  | awk -F'\t' '$11<0.01' \
  | sed -e 's/.representative.fa.alphafold.pdb.gz//' \
  >  combined_results/foldseek.filtered.outfmt6.tmp

cat search_AlphaFold_structures/foldseek.Alphafold-UniProt50.outfmt6 \
  | python ../../../../02_Scripts/grepf_column.py -v -c 16 \
      -f ../03_Lineage_Info/target_taxids.txt \
  | awk -F'\t' '$11<0.01' \
  | sed -e 's/.representative.fa.alphafold.pdb.gz//' \
  >> combined_results/foldseek.filtered.outfmt6.tmp

cat combined_results/foldseek.filtered.outfmt6.tmp \
  | sort -k1,1 -k12,12nr \
  | awk 'BEGIN{print "orthogroup_id\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\ttheader\ttaxid\ttaxname\ttaxlineage"}{print}' \
  > combined_results/foldseek.filtered.outfmt6

rm combined_results/foldseek.filtered.outfmt6.tmp
```





### Overlap functional and strutural annotations



```bash
awk 'NR>1{print $1}' combined.annotations.emapper.withDescFunc.summary.tsv \
  | sort | uniq \
  | awk 'BEGIN{print "OG\tEggNOG_Descriptions"}{print $1"\tyes"}' \
  > t.egg

cat combined.annotations.InterProScan.functional.bed \
  | sed -e 's/-.*//' \
  | sort | uniq \
  | awk 'BEGIN{print "OG\tInterProScan_functional_Annotations"}{print $1"\tyes"}' \
  > t.func

cat combined.annotations.InterProScan.structural.TMHMM.bed \
  | sed -e 's/-.*//' \
  | sort | uniq \
  | awk 'BEGIN{print "OG\tInterProScan_TMHMM_Annotations"}{print $1"\tyes"}' \
  > t.tmhmm

awk -F'\t' 'NR>1{print $1}' foldseek.filtered.outfmt6 \
  | sort | uniq \
  | awk -F'\t' 'BEGIN{print "OG\tFoldseek_hits"}{print $1"\tyes"}' \
  > t.alpha
```



```bash
cat t.* \
  | awk '$1!="OG"{print $1}' \
  | sort | uniq \
  | awk 'BEGIN{print "OG"}{print}' \
  | python ~/scripts/add_value_to_table.py -a t.tmhmm \
  | python ~/scripts/add_value_to_table.py -a t.func  \
  | python ~/scripts/add_value_to_table.py -a t.egg   \
  | python ~/scripts/add_value_to_table.py -a t.alpha \
  > combined.annotations.summary.tsv
```



Get number of OGs with specific yes/no combinations (modify as needed to answer question at hand from big table).

```bash
awk -F'\t' 'NR>1 && ($3=="yes" || $4=="yes"){print $3"\t"$4}' combined.annotations.summary.tsv | sort | uniq -c
```

>    5 	       yes
>    51 yes	
>    17 yes	yes







### *M. capitata* and *P. acuta* diff. accum. results

Grep the number of OGs with Mcap or Pacuta diff expression results (i.e., how many show stress responsiveness).

```bash
# M. capitata HI2018_12TP
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $44"\n"$45}' \
  | sed -e '/^$/d' -e 's/;/\n/g' -e 's/x.*//' \
  | sort | uniq -c \
  | awk 'BEGIN{print "Comparison\tCount"}{print $2"\t"$1}' \
  > "summary.Montipora_capitata-HI2018_12TP.tsv"

# P. acuta HI2018
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $49"\n"$50}' \
  | sed -e '/^$/d' -e 's/;/\n/g' -e 's/x.*//' \
  | sort | uniq -c \
  | awk 'BEGIN{print "Comparison\tCount"}{print $2"\t"$1}' \
  > "summary.Pocillopora_acuta-HI2018.tsv"

# M. capitata HI2019_3TP_Prot
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $46}' \
  | sed -e '/^$/d' -e 's/;/\n/g' -e 's/x.*//' \
  | sort | uniq -c \
  | awk 'BEGIN{print "Comparison\tCount"}{print $2"\t"$1}' \
  > "summary.Montipora_capitata-HI2019_3TP_Prot.tsv"

# M. capitata HI2019_3TP_Trans
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $47}' \
  | sed -e '/^$/d' -e 's/;/\n/g' -e 's/x.*//' \
  | sort | uniq -c \
  | awk 'BEGIN{print "Comparison\tCount"}{print $2"\t"$1}' \
  > "summary.Montipora_capitata-HI2019_3TP_Trans.tsv"
```



Count the number of OGs with diff expre/accum results across any time points or treatments.

```bash
# TOTAL: 271

# M. capitata HI2018_12TP
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $44"\t"$45}' | awk -F'\t' '$1!="" || $2!=""' \
  | wc -l
# 46

# P. acuta HI2018
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $49"\t"$50}' | awk -F'\t' '$1!="" || $2!=""' \
  | wc -l
# 128

# M. capitata HI2019_3TP_Prot
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $46}' | awk -F'\t' '$1!=""' \
  | wc -l
# 0

# M. capitata HI2019_3TP_Trans
cat "combined.stats.OGs.tsv" \
  | awk -F'\t' 'NR>1{print $47}' | awk -F'\t' '$1!=""' \
  | wc -l
# 14
```



Grep the total number of diff expressed genes (so we can check against with OG results)

```bash
# M. capitata HI2018_12TP
cat "/scratch/timothy/projects/0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/transcriptomic/PRJNA731596-Illumina-salmon_Metaproteome/Montipora_capitata_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.tsv" \
  | grep 'Montipora_capitata_KBHIv3___' \
  | awk -F'\t' '$7<0.05 && $3!="NA" && ($3>0.5 || $3<-0.5) {print $8"\t"$9}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > "summary.Montipora_capitata-HI2018_12TP.full.tsv"

# M. capitata HI2019_3TP_Prot
cat "/scratch/timothy/projects/0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/proteomics/MSV000088443-LC_MSMS-Proteome_Discoverer_Metaproteome/Montipora_capitata_MetaProteome.protein_abundance.normalized.DiffAccumResults.tsv" \
  | grep 'Montipora_capitata_KBHIv3___' \
  | awk -F'\t' '$4<0.05 && $2!="NA" && ($2>0.5 || $2<-0.5) {print $6"\t"$7}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > "summary.Montipora_capitata-HI2019_3TP_Prot.full.tsv"

# M. capitata HI2019_3TP_Trans
cat "/scratch/timothy/projects/0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/transcriptomic/PRJNA694677-Illumina-salmon_Metaproteome/Montipora_capitata_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.tsv" \
  | grep 'Montipora_capitata_KBHIv3___' \
  | awk -F'\t' '$7<0.05 && $3!="NA" && ($3>0.5 || $3<-0.5) {print $8"\t"$9}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > "summary.Montipora_capitata-HI2019_3TP_Trans.full.tsv"

# P. acuta HI2018
cat "/scratch/timothy/projects/0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Pocillopora_acuta_KBHIv2/11_omics_data/transcriptomic/PRJNA731596-Illumina-salmon_Metaproteome/Pocillopora_acuta_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.tsv" \
  | grep 'Pocillopora_acuta_KBHIv2___' \
  | awk -F'\t' '$7<0.05 && $3!="NA" && ($3>0.5 || $3<-0.5) {print $8"\t"$9}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' \
  > "summary.Pocillopora_acuta-HI2018.full.tsv"
```



```bash
cat summary.Montipora_capitata-HI2018_12TP.full.tsv \
  | awk -F'\t' '{split($1,a,"_"); split($2,b,"_"); if(a[3]==b[3]){print} }'
cat summary.Pocillopora_acuta-HI2018.full.tsv \
  | awk -F'\t' '{split($1,a,"_"); split($2,b,"_"); if(a[3]==b[3]){print} }'
```



Extract just the OGs with scRNAseq "broadcell" or "meta cell" annotations.

```bash
cat combined.stats.OGs.tsv \
  | awk -F'\t' '{print $1"\t"$55"\t"$58"\t"$56"\t"$59"\t"$62"\t"$60"\t"$63"\t"$65}' \
  > summary.scRNA.tsv
```



Number of OGs with significant FC in "meta cells" in Nematostella_vectensis

```bash
cat summary.scRNA.tsv | awk -F'\t' 'NR!=1 && $4!="" {print $1"\t"$4}' | wc -l
# 28
```

Number of OGs with significant FC in "meta cells" in Stylophora_pistillata ADULT

```bash
cat summary.scRNA.tsv | awk -F'\t' 'NR!=1 && $7!="" {print $1"\t"$7}' | wc -l
# 35
```

Overlap between two species.

```bash
cat summary.scRNA.tsv | awk -F'\t' 'NR!=1 && $4!="" && $7!="" {print $4"\t"$7}' | wc -l
# 5
```









