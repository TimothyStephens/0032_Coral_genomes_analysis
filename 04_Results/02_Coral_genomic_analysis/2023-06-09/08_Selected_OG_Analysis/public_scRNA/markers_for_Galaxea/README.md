# Reformat publically avaliable ssRNA-seq data

Reformat data from the *Stylophora pistillata* ssRNA-Seq manuscript (https://doi.org/10.1016/j.cell.2021.04.005).

> Spis -> Stylophora pistillata (this study)
> Xesp -> Xenia sp. (10.1038/s41586-020-2385-7)
> Nvec -> Nematostella vectensis (10.1016/j.cell.2018.05.019)
> Hvul -> Hydra vulgaris (10.1126/science.aav9314)



Link to "new" datasets.

```bash
ln -s ../../../01_Orthofinder/files2cluster/Stylophora_pistillata_GAJOv1.genes.pep.faa
ln -s ../../../01_Orthofinder/files2cluster/Xenia_sp_CTEAv1.genes.pep.faa
ln -s ../../../01_Orthofinder/files2cluster/Nematostella_vectensis_RRUSv1.genes.pep.faa
ln -s ../../../01_Orthofinder/files2cluster/Hydra_vulgaris_MIJPv3.genes.pep.faa
```



Setup env.

```bash
conda activate /home/timothy/miniconda3/envs/py27
SCRIPTS="../../../../../../02_Scripts"
```



## Orthofinder

Get all orthogroups which contain:

- Hydra_vulgaris_MIJPv3
- Nematostella_vectensis_RRUSv1
- Stylophora_pistillata_GAJOv1
- Xenia_sp_CTEAv1

Which are our target taxa since they are the ones which have scRNA.



Get Orthogroup IDs, Orthogroups in "long" format, and the classtification details of each Orthogroup.

```bash
python "$SCRIPTS/grepf_column.py" \
    -i <(zcat "../../../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz") \
    -f <(echo -e "Hydra_vulgaris_MIJPv3\nNematostella_vectensis_RRUSv1\nStylophora_pistillata_GAJOv1\nXenia_sp_CTEAv1") \
    -c 2 \
  | cut -f1 \
  | sort \
  | uniq \
  > "scRNA.Orthogroups.txt"

python "$SCRIPTS/grepf_column.py" \
    -i <(zcat "../../../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz") \
    -f "scRNA.Orthogroups.txt" \
    -c 1 --keep_header \
  > "scRNA.Orthogroups.long.tsv"

python "$SCRIPTS/grepf_column.py" \
    -i <(zcat "../../../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz") \
    -f "scRNA.Orthogroups.txt" \
    -c 1 --keep_header \
  > "scRNA.Orthogroups.classification.txt"
```



Add scRNA-seq FC values to Orthogroup "long" information.

```bash
cat "scRNA.Orthogroups.long.tsv" \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 3 -a <(sed -e 's/new_name/sequence_id/' ../data/Stylophora_pistillata_GAJOv1.COMBINED_genes.cell_type_gene_FC.tsv | cut -f1,3-) \
  | awk -F'\t' '$4!=""' \
  > "scRNA.Orthogroups.cell_type_gene_FC.tsv"
```



Count the number of genes per Orthogroup which have FC > 1.5

```bash
cat "scRNA.Orthogroups.cell_type_gene_FC.tsv" \
  | awk -F'\t' '{
      if(NR==1){
        for(i=4; i<=NF; i++){
          HEADER[i]=$i
        }
      } else {
        for(i=4; i<=NF; i++){
          if($i>1.5){
            RESULTS[$1][HEADER[i]]++
          }
        }
      }
    } END{
      S="orthogroup_id"
      for(i in HEADER){
        S=S"\t"HEADER[i]
      }
      print S
      for(i in RESULTS){
        S=i
        for(j in HEADER){
          C=RESULTS[i][HEADER[j]]
          if(C==""){C=0}
          S=S"\t"C
        }
        print S
      }
    }' > "scRNA.Orthogroups.cell_type_gene_FC1.5_counts.tsv"
```







## Kevin top genes

Cleanup results from Kevin

```bash
sed -e 's/[^,]*,//' -e 's/"//g' -e 's/,/\t/g' cell.markers.spec.csv > cell.markers.spec.cleaned.csv
sed -e 's/[^,]*,//' -e 's/"//g' -e 's/,/\t/g' cell.markers.broad.csv > cell.markers.broad.cleaned.csv
```



Link old names that Kevin sent with new names that I have. 

```bash
cat ../data/*n2o.tsv \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"}{print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 2 \
      -a <(sed -e '1 s/gene/old_name/' cell.markers.broad.cleaned.csv) \
  | awk -F'\t' '$3!=""' \
  > "cell.markers.broad.cleaned.renamed.tsv"

cat ../data/*n2o.tsv \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"}{print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 2 \
      -a <(sed -e '1 s/gene/old_name/' cell.markers.spec.cleaned.csv) \
  | awk -F'\t' '$3!=""' \
  > "cell.markers.spec.cleaned.renamed.tsv"
```



Add Kevin's data to the "long" formatted Orthogroups - Just the OG info for the 4 taxa in which we have scRNA-seq data

```bash
cat "scRNA.Orthogroups.long.tsv" \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 3 \
      -a <(sed -e '1 s/new_name/sequence_id/' "cell.markers.spec.cleaned.renamed.tsv") \
  | awk -F'\t' '$4!=""' \
  > "cell.markers.spec.cleaned.renamed.long.tsv"

cat "scRNA.Orthogroups.long.tsv" \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 3 \
      -a <(sed -e '1 s/new_name/sequence_id/' "cell.markers.broad.cleaned.renamed.tsv") \
  | awk -F'\t' '$4!=""' \
  > "cell.markers.broad.cleaned.renamed.long.tsv"
```



Add Kevin's data to the "long" formatted Orthogroups - Just the OG info for the 4 taxa in which we have scRNA-seq data + *Galaxea* - Showing all genes, not just the ones which have scRNA-seq data.

```bash
cat "scRNA.Orthogroups.long.tsv" \
  | python "$SCRIPTS/grepf_column.py" \
    -f <(awk 'NR>1{print $1}' "cell.markers.broad.cleaned.renamed.long.tsv" | uniq) \
    -c 1 \
    --keep_header \
  | python "$SCRIPTS/grepf_column.py" \
    -f <(echo -e "Hydra_vulgaris_MIJPv3\nNematostella_vectensis_RRUSv1\nStylophora_pistillata_GAJOv1\nXenia_sp_CTEAv1\nGalaxea_fascicularis_REEFv1") \
    -c 2 \
    --keep_header \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -a <(sed -e '1 s/new_name/sequence_id/' "cell.markers.broad.cleaned.renamed.tsv") \
      -c 3 \
  > "cell.markers.broad.cleaned.renamed.long.all.tsv"

cat "scRNA.Orthogroups.long.tsv" \
  | python "$SCRIPTS/grepf_column.py" \
    -f <(awk 'NR>1{print $1}' "cell.markers.spec.cleaned.renamed.long.tsv" | uniq) \
    -c 1 \
    --keep_header \
  | python "$SCRIPTS/grepf_column.py" \
    -f <(echo -e "Hydra_vulgaris_MIJPv3\nNematostella_vectensis_RRUSv1\nStylophora_pistillata_GAJOv1\nXenia_sp_CTEAv1\nGalaxea_fascicularis_REEFv1") \
    -c 2 \
    --keep_header \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -a <(sed -e '1 s/new_name/sequence_id/' "cell.markers.spec.cleaned.renamed.tsv") \
      -c 3 \
  > "cell.markers.spec.cleaned.renamed.long.all.tsv"
```



Test: Extract *Galaxea* and *Stylophora* genes for `OG0000014` and get 1-to-1 orthologs. This will let us see how easy/hard it will be to connect the two datasets and transfer over the cell type annotations.

```bash
awk '$1=="OG0000014" && $2=="Galaxea_fascicularis_REEFv1"{print $3}' cell.markers.spec.cleaned.renamed.long.all.tsv | xargs ~/programs/samtools-1.11/bin/samtools faidx ../data/Galaxea_fascicularis_REEFv1.genes.pep.faa > OG0000014.Galaxea_fascicularis_REEFv1
awk '$1=="OG0000014" && $2=="Stylophora_pistillata_GAJOv1"{print $3}' cell.markers.spec.cleaned.renamed.long.all.tsv | xargs ~/programs/samtools-1.11/bin/samtools faidx ../data/Stylophora_pistillata_GAJOv1.genes.pep.faa > OG0000014.Stylophora_pistillata_GAJOv1

perl "../data/Find_Orthologs/find_orthologs.py" \
  -i1 "OG0000014.Stylophora_pistillata_GAJOv1" \
  -i2 "OG0000014.Galaxea_fascicularis_REEFv1" \
  -o "OG0000014.RBH" \
  -t p \
  -b "/home/timothy/programs/ncbi-blast-2.13.0+/bin" \
  -c 48
```



Calculate 1-to-1 orthologs between all *Galaxea* and *Stylophora* genes.

```bash
perl "../data/Find_Orthologs/find_orthologs.py" \
  -i1 "all.Stylophora_pistillata_GAJOv1.genes.pep.faa" \
  -i2 "all.Galaxea_fascicularis_REEFv1.genes.pep.faa" \
  -o "all.RBH" \
  -t p \
  -b "/home/timothy/programs/ncbi-blast-2.13.0+/bin" \
  -c 48

# Filter by >60% query pr >60% subject coverage
awk -F'\t' 'NR==1 || $9>60 || $13>60' "all.RBH.results.tsv" > "all.RBH.filtered.results.tsv"
awk -F'\t' 'NR>1{print $1"\t"$2}' "all.RBH.filtered.results.tsv" > "all.RBH.filtered.pairs.tsv"
```



Add OG information to filtered pairs.

```bash
cat "all.RBH.filtered.pairs.tsv" \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 1 \
      -a <(awk '{print $3"\t"$1}' "scRNA.Orthogroups.long.tsv") \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 2 \
      -a <(awk '{print $3"\t"$1}' "scRNA.Orthogroups.long.tsv") \
  | awk -F'\t' '{print $1"\t"$3"\t"$2"\t"$4}' \
  | awk -F'\t' 'BEGIN{print "galaxea_geneid\tgalaxea_geneid_OG\tstylophora_geneid\tstylophora_geneid_OG"}{print}' \
  > "all.RBH.filtered.pairs.OGsAdded.tsv"
```



Add scRNA-seq FC info to filtered 1-to-1 orthologs.

```bash
cat "all.RBH.filtered.pairs.OGsAdded.tsv" \
  | awk -F'\t' 'NR==1 || $2==$4' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" \
      -c 3 \
      -a <(sed -e 's/new_name/stylophora_geneid/' -e 's/old_name/stylophora_old_geneid/' "../data/Stylophora_pistillata_GAJOv1.COMBINED_genes.cell_type_gene_FC.tsv") \
  | awk -F'\t' '$5!=""' \
  > "all.RBH.pairs.COMBINED_genes.cell_type_gene_FC.tsv"
```









Get list of genes, and their respective cell types, with FC > 2.5 with no other cell types between 2.5-1.5 FC.

```bash
TOP_FC_CUTOFF="2.5"
MID_FC_CUTOFF="1.5"
```



```bash
cat "all.RBH.pairs.COMBINED_genes.cell_type_gene_FC.tsv" \
  | cut -f1-42 \
  | awk -F'\t' -vTOP_FC_CUTOFF="${TOP_FC_CUTOFF}" -vMID_FC_CUTOFF="${MID_FC_CUTOFF}" '{
      if(NR==1){
        for(i=6; i<=NF; i++){
          HEADER[i]=$i
        }
      } else {
        for(i=6; i<=NF; i++){
          if($i>TOP_FC_CUTOFF && $i!="NA"){
            TOP_FC[$1][HEADER[i]]=$i
            split(HEADER[i],a,"_")
            TOP_TYPE[$1][a[1]"_"a[2]]++
          } else if($i<=TOP_FC_CUTOFF && $i>MID_FC_CUTOFF && $i!="NA"){
            MID_FC[$1][HEADER[i]]=$i
            MID_TYPE[$1][a[1]"_"a[2]]++
          }
        }
      }
    } END {
      print "new_name\tcell_type_count\tbroad_type_count\tcell_type\tbroad_type\tFC\tmid_range_cell_type_count\tmid_range_broad_type_count"
      for(i in TOP_FC){
        for(j in TOP_FC[i]){
          split(j,a,"_")
          print i"\t"length(TOP_FC[i])"\t"length(TOP_TYPE[i])"\t"j"\t"a[1]"_"a[2]"\t"TOP_FC[i][j]"\t"length(MID_FC[i])"\t"length(MID_TYPE[i])
        }
      }
    }' > "all.RBH.pairs.COMBINED_genes.cell_type_gene_FC${TOP_FC_CUTOFF}_noMid${MID_FC_CUTOFF}_counts.tsv"
```



Filter for genes that are "specific" and "broad" cell type markers.

```bash
F="all.RBH.pairs.COMBINED_genes.cell_type_gene_FC${TOP_FC_CUTOFF}_noMid${MID_FC_CUTOFF}_counts"
cat "$F.tsv" \
  | awk -F'\t' 'NR==1 || ($2==1 && $7==0)' \
  > "$F.specific_celltype_markers.tsv"
cat "$F.tsv" \
  | awk -F'\t' 'NR==1 || ($3==1 && $8<=1)' \
  > "$F.broad_celltype_markers.tsv"
```



Get the number of genes that are "markers" for each specific cell type and each broad cell type. 

```bash
cat "../data/Stylophora_pistillata_GAJOv1.celltype_color.tsv" \
  | grep 'adult' \
  | cut -f4 \
  | uniq \
  | python "$SCRIPTS/add_value_to_table.py" \
      -a <(cat "$F.specific_celltype_markers.tsv" \
             | cut -f4 \
             | sort \
             | uniq -c \
             | awk '{print $2"\t"$1}'
         ) \
      -d 0 \
  | awk -F'\t' 'BEGIN{print "specific_cell_type\tmarker_gene_count"}{print}' \
  > "$F.specific_celltype_markers.counts.tsv"

cat "../data/Stylophora_pistillata_GAJOv1.celltype_color.tsv" \
  | grep 'adult' \
  | cut -f6 \
  | uniq \
  | python "$SCRIPTS/add_value_to_table.py" \
      -a <(cat "$F.broad_celltype_markers.tsv" \
             | cut -f5 \
             | sort \
             | uniq -c \
             | awk '{print $2"\t"$1}'
         ) \
      -d 0 \
  | awk -F'\t' 'BEGIN{print "broad_cell_type\tmarker_gene_count"}{print}' \
  > "$F.broad_celltype_markers.counts.tsv"
```



Download results.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/08_Gene_Structure/public_ssRNA/orthogroups/all.RBH.pairs.COMBINED_genes.cell_type_gene_FC2* . --dry-run
```









