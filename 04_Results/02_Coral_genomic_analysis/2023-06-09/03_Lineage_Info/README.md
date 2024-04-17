# Lineage information

Sort and pre-process the NCBIs taxon ID information/dataset so that we know which IDs blog to our target high level groups.

## Lineage information

Working directory: `03_Lineage_Info/`

```bash
mkdir 03_Lineage_Info; cd 03_Lineage_Info/
ln -s ../../../../02_Scripts/add_value_to_table_SQLite3.py
```

Get lineage information for each taxon ID.

Link to tax dump files.

```bash
ln -s /scratch/databases/taxid_dump/nodes.dmp
ln -s /scratch/databases/taxid_dump/names.dmp 
```

Run `ncbitax2lin` to get lineage information.

```bash
./run_ncbitax2lin.sh
```

Extract the following columns from comma-separated file.

>  tax_id,species,genus,subfamily,family,suborder,order,subclass,class,phylum,clade1,kingdom,clade,superkingdom,no rank

```bash
zcat lineages.txt.gz \
  | awk -F',' '{print $1"\t"$8"\t"$7"\t"$57"\t"$6"\t"$60"\t"$5"\t"$55"\t"$4"\t"$3"\t"$11"\t"$38"\t"$10"\t"$2"\t"$40}' \
  | awk -F'\t' 'BEGIN{OFS=FS="\t"}{ for(i=2; i<=NF; i++){ if($i==""){$i="NA"} }; print }' \
  | gzip -c \
  > lineages.slim.tsv.gz
```

Results file: `lineages.slim.tsv.gz`



Get lineage information for each species in `orthofinder` analysis. 

Create file linking: specie_id <-> taxid (Add isolate column for cases where a orthogroup is specific to a dataset)

```bash
awk -F'\t' '{print $1"\t"$8}' ../samples.txt > species_taxids.txt

./add_value_to_table_SQLite3.py \
    -c 2 -i species_taxids.txt \
    -a <(zcat lineages.slim.tsv.gz) \
  | awk 'BEGIN{OFS=FS="\t"}{ if(NR==1){$2=$2""FS"isolate"}else{$2=$2""FS""$1}; print }' \
  | grep 
  > species_taxids.lineage.tsv
```

Results file: `species_taxids.lineage.tsv`



Extract `tax_id`s from `phylum` of `Cnidaria,Ctenophora,Placozoa,Porifera` AND `class` of `Choanoflagellata` so that we can ignore them downstream since they are the focus of our orthofinder analysis.

```bash
rm -f target_taxa.txt target_taxids.txt

zcat lineages.slim.tsv.gz | awk -F'\t' '$10=="Cnidaria"'        >> target_taxa.txt
zcat lineages.slim.tsv.gz | awk -F'\t' '$10=="Ctenophora"'      >> target_taxa.txt
zcat lineages.slim.tsv.gz | awk -F'\t' '$10=="Placozoa"'        >> target_taxa.txt
zcat lineages.slim.tsv.gz | awk -F'\t' '$10=="Porifera"'        >> target_taxa.txt
zcat lineages.slim.tsv.gz | awk -F'\t' '$9=="Choanoflagellata"' >> target_taxa.txt

cut -f1 target_taxa.txt > target_taxids.txt

# Check numbers
grep -ic 'Cnidaria' target_taxa.txt
#13649
grep -ic 'Ctenophora' target_taxa.txt
#313
grep -ic 'Placozoa' target_taxa.txt
#115
grep -ic 'Porifera' target_taxa.txt
#5863
grep -ic 'Choanoflagellata' target_taxa.txt
#155
```

Results file: `target_taxids.txt`

