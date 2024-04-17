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





## Format *Xenia* sp.

Manuscript that our data are from the same publicaton. Convert "old" IDs in manuscript dataset to the "new" IDs that I used for the manuscript.

```bash
zcat Levy_2021_ManuscriptData/Xesp_gene_annotation.tsv.gz \
  | cut -f1 \
  | awk -F'\t' '{print "Xenia_sp_CTEAv1___Xe_"$1"-T1\t"$1}' \
  | sed -e 's/Xesp_//' \
  | python "$SCRIPTS/grepf_column.py" -f <(grep '>' Xenia_sp_CTEAv1.genes.pep.faa | sed -e 's/>//') \
  > Xenia_sp_CTEAv1.genes.n2o.tsv
```

Extract and rename data from Levy 2021 using the "new" names that we have now linked to the "old" names.  

```bash
P1="Xenia_sp_CTEAv1"
P2="Xesp"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_broad_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.broad_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.metacell_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_transcription_factors.tsv.gz" | awk 'BEGIN{print "old_name\ttranscription_factor_1\tTtranscription_factor_2"} {print}') \
  > "$P1.genes.transcription_factors.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_gene_annotation.tsv.gz" | awk 'BEGIN{print "old_name\tannotation_1\tannotation_2"} {print}') \
  > "$P1.genes.gene_annotation.tsv"

```



Metacell color schemes used by Levy 2021.

```bash
echo -e "tissue\ttissue_color\t\tmetacell_color\tcell\tcell_color\tbroadcell\tbroadcell_color" \
  > Xenia_sp_CTEAv1.celltype_color.tsv

awk -F'\t' 'NR>1{print "adult\t#fb9a99\t"$0}' \
    Levy_2021_GithubData/clustering_xenia/scdb/Xesp_metacell_annotation \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$5}' \
  | ../../../../../../02_Scripts/add_value_to_table.py \
      -c 5 \
      -a Levy_2021_GithubData/clustering_xenia/scdb/Xesp_broad_cell_type_annotation \
    >> Xenia_sp_CTEAv1.celltype_color.tsv
```











## Format *Nematostella vectensis*

Manuscript that our data are from (https://doi.org/10.1101/2020.10.30.359448) is the "new" chomosome level assembly. The "old" assembly if from https://doi.org/10.1016/j.cell.2018.05.019

Convert "old" IDs in manuscript dataset to the "new" IDs that I used for the manuscript.

NOTE: Needed to download the "old" dataset from https://tanaylab.github.io/old_resources/pages/724.html as the old lab website is no longer avaliable.

```bash
wget https://arnau-sebe-pedros.s3.eu-west-1.amazonaws.com/Sebe-Pedros2018/reproduce_results.tar.gz
wget https://arnau-sebe-pedros.s3.eu-west-1.amazonaws.com/Single_cell_datasets.tar.gz
```

Reformat "old" names so that they match what was used by Levy 2021.

```bash
cat Sebe-Pedros_2018_WebsiteData/Single_cell_datasets/Nematostella/Nematostella_proteins.fasta \
  | sed -e 's/>/>Nvec_/' \
  > Nvec.pep.faa
```

Check that these gene names overlap perfectly with the ones in the Levy 2021 dataset.

```bash
cat \
  <(grep '>' Nvec.pep.faa | sed -e 's/>//') \
  <(zcat Levy_2021_ManuscriptData/Nvec_gene_annotation.tsv.gz | cut -f1) \
  | sort | uniq -u | wc -l
# 0
```

> All gene names match perfectly.

Use reciprocal-best-hit blastp to find the 1-to-1 orthologs between the"old" and "new" datasets. 

```bash
python Find_Orthologs/find_orthologs.py \
  -i1 Nematostella_vectensis_RRUSv1.genes.pep.faa \
  -i2 Nvec.pep.faa \
  -o  Nvec.pep.faa.RBH \
  -t p \
  -b /home/timothy/programs/ncbi-blast-2.13.0+/bin \
  -c 48

# Filter by >60% query pr >60% subject coverage
awk -F'\t' 'NR==1 || $9>60 || $13>60' \
    Nvec.pep.faa.RBH.results.tsv \
  > Nvec.pep.faa.RBH.filtered.results.tsv
awk 'NR>1{print $1"\t"$2}' \
    Nvec.pep.faa.RBH.filtered.results.tsv \
  > Nvec.pep.faa.RBH.filtered.pairs.tsv
```

> Total 1-to-1 orthologs: 13260
> Filtered 1-to-1 orthologs: 12924

NOTE: There is a significant difference in the number of "old" vs "new" genes, which explains the low number of 1-to-1 orthologs.

> Nvec.pep.faa:32280
> Nematostella_vectensis_RRUSv1.genes.pep.faa:15355

Reformat into "new" vs "old"

```bash
awk -F'\t' '{print $2"\t"$1}' \
    Nvec.pep.faa.RBH.filtered.pairs.tsv \
  > Nematostella_vectensis_RRUSv1.genes.n2o.tsv
```

Extract and rename data from Levy 2021 using the "new" names that we have now linked to the "old" names.  

```bash
P1="Nematostella_vectensis_RRUSv1"
P2="Nvec"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_adult_broad_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.broad_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_adult_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_adult_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.metacell_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_transcription_factors.tsv.gz" | awk 'BEGIN{print "old_name\ttranscription_factor_1\tTtranscription_factor_2"} {print}') \
  > "$P1.genes.transcription_factors.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_gene_annotation.tsv.gz" | awk 'BEGIN{print "old_name\tannotation_1\tannotation_2"} {print}') \
  > "$P1.genes.gene_annotation.tsv"

```



Metacell color schemes used by Levy 2021.

```bash
echo -e "tissue\ttissue_color\tmetacell\tmetacell_color\tcell\tcell_color\tbroadcell\tbroadcell_color" \
  > Nematostella_vectensis_RRUSv1.celltype_color.tsv

awk -F'\t' 'NR>1{print "adult\t#fb9a99\t"$0}' \
    Levy_2021_GithubData/clustering_nematostella/scdb/Nvec_metacell_annotation \
  | sed -e 's/gastrodermis_parietal_circular_prog/gastrodermis_muscle_parietal_circular_prog/' \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$5}' \
  | ../../../../../../02_Scripts/add_value_to_table.py \
      -c 5 \
      -a Levy_2021_GithubData/clustering_nematostella/scdb/Nvec_broad_cell_type_annotation \
    >> Nematostella_vectensis_RRUSv1.celltype_color.tsv
```











## Format *Hydra vulgaris*

Download Hydra genes (v2) from https://research.nhgri.nih.gov/hydra

```bash
wget https://research.nhgri.nih.gov/hydra/download/genemodels_proteins/hydra2.0_genemodels.aa.gz
```

Use reciprocal-best-hit blastp to find the 1-to-1 orthologs between the"old" and "new" datasets. 

```bash
python Find_Orthologs/find_orthologs.py \
  -i1 Hydra_vulgaris_MIJPv3.genes.pep.faa \
  -i2 hydra2.0_genemodels.aa \
  -o  hydra2.0_genemodels.aa.RBH \
  -t p \
  -b /home/timothy/programs/ncbi-blast-2.13.0+/bin \
  -c 48

# Filter by >60% query pr >60% subject coverage
awk -F'\t' 'NR==1 || $9>60 || $13>60' hydra2.0_genemodels.aa.RBH.results.tsv > hydra2.0_genemodels.aa.RBH.filtered.results.tsv
awk 'NR>1{print $1"\t"$2}' hydra2.0_genemodels.aa.RBH.filtered.results.tsv > hydra2.0_genemodels.aa.RBH.filtered.pairs.tsv
```

> Total 1-to-1 orthologs: 12477
> Filtered 1-to-1 orthologs: 12057

NOTE: There is a significant difference in the number of "old" vs "new" genes, which explains the low number of 1-to-1 orthologs.

> hydra2.0_genemodels.aa:36059
> Hydra_vulgaris_MIJPv3.genes.pep.faa:21385

Reformat into "new" vs "old"

```bash
awk -F'\t' '{print $2"\t"$1}' \
    hydra2.0_genemodels.aa.RBH.filtered.pairs.tsv \
  | sed -e 's/Sc4wPfr_.*\.g/Hvul_g/' -e 's/\.t/_/' \
  > Hydra_vulgaris_MIJPv3.genes.n2o.tsv
```

Extract and rename data from Levy 2021 using the "new" names that we have now linked to the "old" names.  

```bash
P1="Hydra_vulgaris_MIJPv3"
P2="Hvul"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_broad_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.broad_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.metacell_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_transcription_factors.tsv.gz" | awk 'BEGIN{print "old_name\ttranscription_factor_1\tTtranscription_factor_2"} {print}') \
  > "$P1.genes.transcription_factors.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_gene_annotation.tsv.gz" | awk 'BEGIN{print "old_name\tannotation_1\tannotation_2"} {print}') \
  > "$P1.genes.gene_annotation.tsv"

```

WARNING: Some of the genes in `Levy_2021_ManuscriptData/${P2}_gene_annotation.tsv.gz` have multiple (non-unique) annotations. Take just the first in the file.



Metacell color schemes used by Levy 2021.

```bash
echo -e "tissue\ttissue_color\tmetacell\tmetacell_color\tcell\tcell_color\tbroadcell\tbroadcell_color" \
  > Hydra_vulgaris_MIJPv3.celltype_color.tsv

awk -F'\t' 'NR>1{print "adult\t#fb9a99\t"$0}' \
    Levy_2021_GithubData/clustering_hydra/scdb/Hvul_metacell_annotation \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$5}' \
  | ../../../../../../02_Scripts/add_value_to_table.py \
      -c 5 \
      -a <(sed -e 's/ecEp_head_hypostome/cnidocyte_ecEp_head_hypostome/' Levy_2021_GithubData/clustering_hydra/scdb/Hvul_broad_cell_type_annotation) \
    >> Hydra_vulgaris_MIJPv3.celltype_color.tsv
```











## Format *Stylophora pistillata*

Download proteins used by Levy 2021 from their R Shiny app (https://sebe-lab.shinyapps.io/Stylophora_cell_atlas/) using the "Download FASTA" button in the "Genes heatmap" subtabs under each of the "Adult coral", "Polyp", and "Larva" top level tabs. Data in `Levy_2021_RShiny/`

Combine datasets and get unique list of proteins.

```bash
cat Levy_2021_RShiny/selected_genes_* \
  | ~/programs/seqkit_v2.3.1/seqkit fx2tab \
  | awk -F'\t' '{print $1"\t"$2}' \
  | sort \
  | uniq \
  | ~/programs/seqkit_v2.3.1/seqkit tab2fx \
  > "Levy_2021_RShiny.combined.faa"
```

Check that we have captured all of the proteins in the Levy 2021 dataset.

```bash
zcat "Levy_2021_ManuscriptData/Spis_adult_broad_cell_type_gene_FC.tsv.gz" \
  | grep -v -f <(grep '>' "Levy_2021_RShiny.combined.faa" | sed -e 's/>//') \
  | grep -v 'orphan_peak_\|Spis_LOC' \
  | cut -f1
```

> Excluding de novo transcripts (ignore for now) the only protein not in our sequence file is "Spis10006_1" - Not sure why it is missing from the RShiny dataset, but we can ignore for now since it unlikely to significantly affect our results.

Use reciprocal-best-hit blastp to find the 1-to-1 orthologs between the"old" and "new" datasets. 

```bash
python Find_Orthologs/find_orthologs.py \
  -i1 Stylophora_pistillata_GAJOv1.genes.pep.faa \
  -i2 Levy_2021_RShiny.combined.faa \
  -o  Levy_2021_RShiny.combined.faa.RBH \
  -t p \
  -b /home/timothy/programs/ncbi-blast-2.13.0+/bin \
  -c 48

# Filter by >60% query pr >60% subject coverage
awk -F'\t' 'NR==1 || $9>60 || $13>60' Levy_2021_RShiny.combined.faa.RBH.results.tsv > Levy_2021_RShiny.combined.faa.RBH.filtered.results.tsv
awk 'NR>1{print $1"\t"$2}' Levy_2021_RShiny.combined.faa.RBH.filtered.results.tsv > Levy_2021_RShiny.combined.faa.RBH.filtered.pairs.tsv
```

> Total 1-to-1 orthologs: 18901
> Filtered 1-to-1 orthologs: 18557

NOTE: The the two protein sets are roghly similar in terms of protein counts, although its possible that some of the NCBI-derived genes do not exist in the GeefGenomics dataset, explaining the slightly lower number of 1-to-1 orthologs.

> Levy_2021_RShiny.combined.faa:24213
> Stylophora_pistillata_GAJOv1.genes.pep.faa:25769

Reformat into "new" vs "old"

```bash
awk -F'\t' '{print $2"\t"$1}' \
    Levy_2021_RShiny.combined.faa.RBH.filtered.pairs.tsv \
  > Stylophora_pistillata_GAJOv1.genes.n2o.tsv
```

Extract and rename data from Levy 2021 using the "new" names that we have now linked to the "old" names.  

```bash
P1="Stylophora_pistillata_GAJOv1"
P2="Spis"

##
## Annotations
##
grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_transcription_factors.tsv.gz" | awk 'BEGIN{print "old_name\ttranscription_factor_1\tTtranscription_factor_2"} {print}') \
  > "$P1.genes.transcription_factors.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_gene_annotation.tsv.gz" | awk 'BEGIN{print "old_name\tannotation_1\tannotation_2\tannotation_3\tannotation_4"} {print}') \
  > "$P1.genes.gene_annotation.tsv"


##
## Adult
##
P3="adult"
grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_broad_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_broad_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_metacell_gene_FC.tsv"


##
## Larva
##
P3="larva"
grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_metacell_gene_FC.tsv"


##
## Polyp
##
P3="polyp"
grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_cell_type_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_cell_type_gene_FC.tsv"

grep '>' "$P1.genes.pep.faa" \
  | sed -e 's/>//' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 1 -d "NA" \
      -a "$P1.genes.n2o.tsv" \
  | awk -F'\t' 'BEGIN{print "new_name\told_name"} {print}' \
  | python "$SCRIPTS/add_value_to_table_SQLite3.py" -c 2 \
      -d $'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA' \
      -a <(zcat "Levy_2021_ManuscriptData/${P2}_${P3}_metacell_gene_FC.tsv.gz" | awk '{ if(NR==1){print "old_name\t"$0}else{print} }') \
  > "$P1.genes.${P3}_metacell_gene_FC.tsv"

```



Combine results from all celltypes from the three tissues together.

```bash
cat Stylophora_pistillata_GAJOv1.genes.n2o.tsv \
  | awk 'BEGIN{print "new_name\told_name"}{print}' \
  | python $SCRIPTS/add_value_to_table_SQLite3.py \
      -a <(cat Stylophora_pistillata_GAJOv1.genes.adult_cell_type_gene_FC.tsv \
             | awk 'BEGIN{OFS=FS="\t"} { if(NR==1){ for(i=3; i<=NF; i++){$i="adult_"$i}; print }else{print} }' \
             | cut -f1,3-) \
  | python $SCRIPTS/add_value_to_table_SQLite3.py \
      -a <(cat Stylophora_pistillata_GAJOv1.genes.larva_cell_type_gene_FC.tsv \
             | awk 'BEGIN{OFS=FS="\t"} { if(NR==1){ for(i=3; i<=NF; i++){$i="larva_"$i}; print }else{print} }' \
             | cut -f1,3-) \
  | python $SCRIPTS/add_value_to_table_SQLite3.py \
      -a <(cat Stylophora_pistillata_GAJOv1.genes.polyp_cell_type_gene_FC.tsv \
             | awk 'BEGIN{OFS=FS="\t"} { if(NR==1){ for(i=3; i<=NF; i++){$i="polyp_"$i}; print }else{print} }' \
             | cut -f1,3-) \
  > Stylophora_pistillata_GAJOv1.COMBINED_genes.cell_type_gene_FC.tsv
```



Combine the color schemes (used by Levy 2021) for each cell type from the three tissues together.

```bash
echo -e "tissue\tmetacell\tmetacell_color\tcell\tcell_color\tbroadcell" \
  > Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp

awk -F'\t' 'NR>1{P="adult"; split($2,a,"_"); print P"\t"$1"\t"$3"\t"$2"\t"$3"\t"a[1]}' \
  Levy_2021_RShiny/metacell_annotation_coral.tsv \
    >> Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp
awk -F'\t' 'NR>1{P="larva"; split($2,a,"_"); print P"\t"$1"\t"$3"\t"$2"\t"$3"\t"a[1]}' \
  Levy_2021_RShiny/metacell_annotation_larva.tsv \
    >> Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp
awk -F'\t' 'NR>1{P="polyp"; split($2,a,"_"); print P"\t"$1"\t"$3"\t"$2"\t"$3"\t"a[1]}' \
  Levy_2021_RShiny/metacell_annotation_polyp.tsv \
    >> Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp

cat Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp \
  | python $SCRIPTS/add_value_to_table_SQLite3.py -a colors.txt \
  | awk -F'\t' '{ if(NR==1){print $0"\t"$6}else{print $0"\t"$1"_"$6} }' \
  | python $SCRIPTS/add_value_to_table_SQLite3.py -c 8 -a colors.txt  \
  | awk -F'\t' '{print $1"\t"$7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9}' \
  > Stylophora_pistillata_GAJOv1.celltype_color.tsv

rm -f Stylophora_pistillata_GAJOv1.celltype_color.tsv.tmp
```











