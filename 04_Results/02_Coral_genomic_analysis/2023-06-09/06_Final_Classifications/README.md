# Compile and cleanup final results files

Gather together the results files to make sure we use the right files to use downstream and that all of the column names are formatted correctly.



OG groupings.

```bash
cat ../05_Dark_OGs/Run2.Orthogroups_ALL.long.tsv \
  | awk -F'\t' 'BEGIN{print "orthogroup_id\tspecies_id\tsequence_id"}{print}' \
  | gzip -c \
  > Orthogroups.Run2.long.tsv.gz
```



OG classifications.

```bash
cat ../05_Dark_OGs/OG_analysis.combined.results.type.tsv \
  | sed -e '1 s/cluster_id/orthogroup_id/' -e '1 s/Designation/designation/' \
  | gzip -c \
  > Orthogroups.Run2.classification.tsv.gz
```



Sequence `nr` top hits.

```bash
cat ../05_Dark_OGs/diamond_results.nr.top_hits \
  | awk -F'\t' 'BEGIN{print "sequence_id\ttopAnnotation_sseqid\ttopAnnotation_evalue\ttopAnnotation_bitscore\ttopAnnotation_staxids\ttopAnnotation_title\ttopTaxon_sseqid\ttopTaxon_evalue\ttopTaxon_bitscore\ttopTaxon_staxids\ttopTaxon_title"}{print}' \
  | gzip -c \
  > Orthogroups.Run2.sequences.nr.top_hits.tsv.gz
```



Number of orthogroups in each functional and taxonomic group.

```bash
gunzip -c Orthogroups.Run2.classification.tsv.gz \
  | awk -F'\t' 'NR>1{print $28"\t"$9"\t"$10}' \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$4"\t"$1}' \
  | python ~/scripts/add_value_to_table_SQLite3.py -c 3 \
      -a <(cat ../03_Lineage_Info/species_taxids.lineage.tsv | awk -F'\t' 'NR>1{if($4=="NA"){print $2"\t"$2" ("$3")"}else{print $2"\t"$2" ("$4")"}}' \
             | sort | uniq) \
  | awk -F'\t' 'BEGIN{print "Designation\tTaxonomic level\tTaxonomic group\tOrthogroup count"} {if($2=="tax_id"){print $1"\t"$2"\t"$5"\t"$4}else{print $1"\t"$2"\t"$3"\t"$4}}' \
  > Orthogroups.Run2.classification.count.tsv
```





## Download results

Download results from Coral server.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/06_Final_Classifications/* . --dry-run
```









## Get expresion results for each final OG

```bash
mkdir -p expression_results; cd expression_results/
```

Produce a table with the number of genes with significant expression under a given condition 

```bash
PREFIX="DiffExprResults.Mcapitata.3TP"
FILE="../../../../../../0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/transcriptomic/PRJNA694677-Illumina-salmon_Metaproteome/Montipora_capitata_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.tsv"

# Get all genes with DiffExpr results > cutoff values
cat "${FILE}" \
  | awk -F'\t' 'NR==1 || $1~"^Montipora_capitata_KBHIv3"' \
  | awk -F'\t' '$3!="NA" && ($3<-0.5 || $3>0.5) && $7!="NA" && $7<0.05 {
      print $1"\t"$8"vs"$9"\t"$3"\t"$7
    }' \
  | sed -e 's/MC-289_TP1-Amb/Ambient1/g' \
  | sed -e 's/MC-289_TP1-HiT/HighTemp1/g' \
  | sed -e 's/MC-289_TP3-Amb/Ambient3/g' \
  | sed -e 's/MC-289_TP3-HiT/HighTemp3/g' \
  | sed -e 's/MC-289_TP5-Amb/Ambient5/g' \
  | sed -e 's/MC-289_TP5-HiT/HighTemp5/g' \
  | grep 'Ambient1vsHighTemp1\|Ambient3vsHighTemp3\|Ambient5vsHighTemp5' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  | sort \
  > "${PREFIX}.filtered.long.tsv"





# Get total number of genes per designation (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Montipora_capitata_KBHIv3" {print $3"\t"$4}' \
  | awk -F'\t' '{seen[$2]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation.tsv"

# Get total number of genes per high level designation (Light vs. Dark) (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Montipora_capitata_KBHIv3" {print $3"\t"$4}' \
  | awk -F'\t' '{split($2,a,"-"); seen[a[1]]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation_LvD.tsv"





# Get the total number of (unique) genes per designation
cat "${PREFIX}.filtered.long.tsv" \
	| awk -F'\t' '{print $5"\t"$1}' \
	| sort | uniq \
	| awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation.tsv"

# Get the total number of (unique) genes per designation AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Montipora_capitata_KBHIv3"{print $8"vs"$9}' \
  | sed -e 's/MC-289_TP1-Amb/Ambient1/g' \
  | sed -e 's/MC-289_TP1-HiT/HighTemp1/g' \
  | sed -e 's/MC-289_TP3-Amb/Ambient3/g' \
  | sed -e 's/MC-289_TP3-HiT/HighTemp3/g' \
  | sed -e 's/MC-289_TP5-Amb/Ambient5/g' \
  | sed -e 's/MC-289_TP5-HiT/HighTemp5/g' \
  | grep 'Ambient1vsHighTemp1\|Ambient3vsHighTemp3\|Ambient5vsHighTemp5' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark-Shared,"$1"\nLight-Restricted,"$1"\nDark-Restricted,"$1"\nLight-Shared,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{print $5"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"





# Get the total number of (unique) genes per high level designation (Light vs. Dark)
cat "${PREFIX}.filtered.long.tsv" \
  | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$1}' \
  | sort | uniq \
  | awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv"

# Get the total number of (unique) genes per high level designation (Light vs. Dark) AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Montipora_capitata_KBHIv3"{print $8"vs"$9}' \
  | sed -e 's/MC-289_TP1-Amb/Ambient1/g' \
  | sed -e 's/MC-289_TP1-HiT/HighTemp1/g' \
  | sed -e 's/MC-289_TP3-Amb/Ambient3/g' \
  | sed -e 's/MC-289_TP3-HiT/HighTemp3/g' \
  | sed -e 's/MC-289_TP5-Amb/Ambient5/g' \
  | sed -e 's/MC-289_TP5-HiT/HighTemp5/g' \
  | grep 'Ambient1vsHighTemp1\|Ambient3vsHighTemp3\|Ambient5vsHighTemp5' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark,"$1"\nLight,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"
```

```bash
PREFIX="DiffExprResults.Mcapitata.12TP"
FILE="../../../../../../0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/11_omics_data/transcriptomic/PRJNA731596-Illumina-salmon_Metaproteome/Montipora_capitata_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.Ctl_vs_Treat.tsv"

# Get all genes with DiffExpr results > cutoff values
cat "${FILE}" \
  | awk -F'\t' 'NR==1 || $1~"^Montipora_capitata_KBHIv3"' \
  | awk -F'\t' '$3!="NA" && ($3<-0.5 || $3>0.5) && $7!="NA" && $7<0.05 {
      sub("Mcapitata_","",$8); sub("_TP","",$8); 
      sub("Mcapitata_","",$9); sub("_TP","",$9); 
      print $1"\t"$8"vs"$9"\t"$3"\t"$7
    }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  | sort \
  > "${PREFIX}.filtered.long.tsv"





# Get total number of genes per designation (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Montipora_capitata_KBHIv3" {print $3"\t"$4}' \
  | awk -F'\t' '{seen[$2]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation.tsv"

# Get total number of genes per high level designation (Light vs. Dark) (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Montipora_capitata_KBHIv3" {print $3"\t"$4}' \
  | awk -F'\t' '{split($2,a,"-"); seen[a[1]]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation_LvD.tsv"





# Get the total number of (unique) genes per designation
cat "${PREFIX}.filtered.long.tsv" \
	| awk -F'\t' '{print $5"\t"$1}' \
	| sort | uniq \
	| awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation.tsv"

# Get the total number of (unique) genes per designation AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Montipora_capitata_KBHIv3"{
        sub("Mcapitata_","",$8); sub("_TP","",$8); 
        sub("Mcapitata_","",$9); sub("_TP","",$9); 
        print $8"vs"$9
      }' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark-Shared,"$1"\nLight-Restricted,"$1"\nDark-Restricted,"$1"\nLight-Shared,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{print $5"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"





# Get the total number of (unique) genes per high level designation (Light vs. Dark)
cat "${PREFIX}.filtered.long.tsv" \
  | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$1}' \
  | sort | uniq \
  | awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv"

# Get the total number of (unique) genes per high level designation (Light vs. Dark) AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Montipora_capitata_KBHIv3"{
      sub("Mcapitata_","",$8); sub("_TP","",$8); 
      sub("Mcapitata_","",$9); sub("_TP","",$9); 
      print $8"vs"$9
    }' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark,"$1"\nLight,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"
```

```bash
PREFIX="DiffExprResults.Pacuta.12TP"
FILE="../../../../../../0031_Coral_genomes_NEW/03_Analysis/2022-02-03/coral_genomes/Pocillopora_acuta_KBHIv2/11_omics_data/transcriptomic/PRJNA731596-Illumina-salmon_Metaproteome/Pocillopora_acuta_MetaTranscriptome.cds.fna.numreads.matrix.DiffExprResults.Ctl_vs_Treat.tsv"

# Get all genes with DiffExpr results > cutoff values
cat "${FILE}" \
  | awk -F'\t' 'NR==1 || $1~"^Pocillopora_acuta_KBHIv2"' \
  | awk -F'\t' '$3!="NA" && ($3<-0.5 || $3>0.5) && $7!="NA" && $7<0.05 {
      sub("Pacuta_","",$8); sub("_TP","",$8); 
      sub("Pacuta_","",$9); sub("_TP","",$9); 
      print $1"\t"$8"vs"$9"\t"$3"\t"$7
    }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  | sort \
  > "${PREFIX}.filtered.long.tsv"





# Get total number of genes per designation (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Pocillopora_acuta_KBHIv2" {print $3"\t"$4}' \
  | awk -F'\t' '{seen[$2]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation.tsv"

# Get total number of genes per high level designation (Light vs. Dark) (regardness of expression results)
zcat ../Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '$2=="Pocillopora_acuta_KBHIv2" {print $3"\t"$4}' \
  | awk -F'\t' '{split($2,a,"-"); seen[a[1]]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.no_genes_per_designation_LvD.tsv"





# Get the total number of (unique) genes per designation
cat "${PREFIX}.filtered.long.tsv" \
	| awk -F'\t' '{print $5"\t"$1}' \
	| sort | uniq \
	| awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation.tsv"

# Get the total number of (unique) genes per designation AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Pocillopora_acuta_KBHIv2"{
        sub("Pacuta_","",$8); sub("_TP","",$8); 
        sub("Pacuta_","",$9); sub("_TP","",$9); 
        print $8"vs"$9
      }' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark-Shared,"$1"\nLight-Restricted,"$1"\nDark-Restricted,"$1"\nLight-Shared,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{print $5"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_and_exprCond.prop.tsv"





# Get the total number of (unique) genes per high level designation (Light vs. Dark)
cat "${PREFIX}.filtered.long.tsv" \
  | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$1}' \
  | sort | uniq \
  | awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv"

# Get the total number of (unique) genes per high level designation (Light vs. Dark) AND expression conditions combination
cat "${FILE}" \
  | awk -F'\t' '$1~"^Pocillopora_acuta_KBHIv2"{
      sub("Pacuta_","",$8); sub("_TP","",$8); 
      sub("Pacuta_","",$9); sub("_TP","",$9); 
      print $8"vs"$9
    }' \
  | sort | uniq \
  | awk -F'\t' '{ print "Dark,"$1"\nLight,"$1 }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "${PREFIX}.filtered.long.tsv" \
             | awk -F'\t' '{split($5,a,"-"); print a[1]"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "${PREFIX}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "${PREFIX}.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "${PREFIX}.filtered.long.no_genes_per_designation_LvD_and_exprCond.prop.tsv"
```



```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/06_Final_Classifications/expression_results/DiffExprResults.*_exprCond.prop.tsv .
```





```bash
zcat ../Orthogroups.Run2.long.tsv.gz \
  | grep 'Montipora_capitata_KBHIv3' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '{print $4}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$1}'
```

>Dark-Restricted	3893
>Dark-Shared	     150
>Light-Restricted	1087
>Light-Shared	     49254

```bash
zcat ../Orthogroups.Run2.long.tsv.gz \
  | grep 'Pocillopora_acuta_KBHIv2' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
  | awk -F'\t' '{print $4}' \
  | sort | uniq -c \
  | awk '{print $2"\t"$1}'
```

>Dark-Restricted	 2551
>Dark-Shared	      63
>Light-Restricted	1113
>Light-Shared	     30003











```bash
#for FILE in ../../08_Selected_OG_Analysis/public_scRNA/data/*.genes.*_gene_FC.tsv;
for FILE in ../../08_Selected_OG_Analysis/public_scRNA/data/Hydra_vulgaris_MIJPv3.genes.metacell_gene_FC.tsv;
do
  echo "## $FILE"
  F=$(basename ${FILE%*.tsv}); 
  SP=${F%*.gene*}; 
  TYPE=${F#*.genes.*};
  echo "## $SP --- $TYPE"
  CUTOFF="1.5"

# Get all genes with at least one cell type with > cutoff FC values
cat "${FILE}" \
  | awk -F'\t' -vCUTOFF="${CUTOFF}" '{
      if(NR==1){
        print
      } else {
        ABOVE=0
        for(i=3; i<=NF; i++){
          if($i!="NA" && $i>CUTOFF){
            ABOVE=1
          }
        }
        if(ABOVE == 1){print}
      }
    }' \
  > "scRNA.${SP}.${TYPE}.filtered.tsv"

# Convert to long (R melted) format - Add designation to each line, remove cell types without > cutoff FC values
cat "scRNA.${SP}.${TYPE}.filtered.tsv" \
  | awk -F'\t' -vCUTOFF="${CUTOFF}" '{
      if(NR==1){
        for(i=3; i<=NF; i++){COLNAMES[i]=$i}
      } else {
        for(i=3; i<=NF; i++){
          if($i!="NA" && $i>CUTOFF){
            print $1"\t"COLNAMES[i]"\t"$i
          }
        }
      }
    }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  > "scRNA.${SP}.${TYPE}.filtered.long.tsv"





# Get total number of genes per designation (regardness of cell type expression results)
cat "${FILE}" \
  | awk -F'\t' '{
      if(NR>1){
        print $1
      }
    }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  | awk -F'\t' '{seen[$2]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "scRNA.${SP}.${TYPE}.no_genes_per_designation.tsv"

# Get total number of genes per high level designation (Light vs. Dark) (regardness of cell type expression results)
cat "${FILE}" \
  | awk -F'\t' '{
      if(NR>1){
        print $1
      }
    }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a <(zcat ../Orthogroups.Run2.long.tsv.gz \
             | python ~/scripts/add_value_to_table_SQLite3.py \
                 -a <(zcat ../Orthogroups.Run2.classification.tsv.gz | awk -F'\t' '{print $1"\t"$28}') \
             | awk -F'\t' '{print $3"\t"$4}' \
         ) \
  | awk -F'\t' '{split($2,a,"-"); seen[a[1]]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "scRNA.${SP}.${TYPE}.no_genes_per_designation_LvD.tsv"





# Get the total number of (unique) genes per designation
cat "scRNA.${SP}.${TYPE}.filtered.long.tsv" \
	| awk -F'\t' '{print $4"\t"$1}' \
	| sort | uniq \
	| awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation.tsv"

# Get the total number of (unique) genes per designation AND cell type combo
cat "${FILE}" \
  | awk -F'\t' 'NR==1{ for(i=3; i<=NF; i++){print "Dark-Shared,"$i"\nLight-Restricted,"$i"\nDark-Restricted,"$i"\nLight-Shared,"$i} }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "scRNA.${SP}.${TYPE}.filtered.long.tsv" \
             | awk -F'\t' '{print $4"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "scRNA.${SP}.${TYPE}.no_genes_per_designation.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.prop.tsv"

# Add colors to values
COLOR_TYPE=$(awk -vTYPE="$TYPE" 'BEGIN{if(TYPE~"^larva"){print "larva"}else if(TYPE~"^polyp"){print "polyp"}else{print "adult"} }')
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.prop.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -c 2 -d $'Total\twhite' \
      -a <(cat "../../08_Selected_OG_Analysis/public_scRNA/data/${SP}.celltype_color.tsv" \
             | awk -F'\t' -vCOLOR_TYPE="${COLOR_TYPE}" 'NR>1 && $1==COLOR_TYPE{print $7"\t"$7"\t"$8; print $5"\t"$7"\t"$8; print $3"\t"$7"\t"$8}' \
             | sort | uniq \
             | awk 'BEGIN{print "Conditions\tBroadCellType\tBroadCellColor"}{print}'\
          ) \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_and_cellType.prop.colors.tsv"





# Get the total number of (unique) genes per high level designation (Light vs. Dark)
cat "scRNA.${SP}.${TYPE}.filtered.long.tsv" \
  | awk -F'\t' '{split($4,a,"-"); print a[1]"\t"$1}' \
  | sort | uniq \
  | awk -F'\t' '{seen[$1]++}END{ for(d in seen){print d"\t"seen[d]} }' \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD.tsv"

# Get the total number of (unique) genes per high level designation (Light vs. Dark) AND cell type combo
cat "${FILE}" \
  | awk -F'\t' 'NR==1{ for(i=3; i<=NF; i++){print "Dark,"$i"\nLight,"$i} }' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -d 0 \
      -a <(cat "scRNA.${SP}.${TYPE}.filtered.long.tsv" \
             | awk -F'\t' '{split($4,a,"-"); print a[1]"\t"$2"\t"$1}' \
             | sort | uniq \
             | awk -F'\t' '{seen[$1][$2]++; cont[$2]++}END{ for(d in seen){ for(c in cont){n=0; if(seen[d][c]!=""){n=seen[d][c]} print d","c"\t"n} } }' \
          ) \
  | sed -e 's/,/\t/' \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.tsv"

# Add header to (4 column) file
echo -e "Designation\tConditions\tCount\tPercent" \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.prop.tsv"

# Get proportion of genes in each designation-cellType combo, out of all genes in
# a designation with FC > cutoff (so not out of total, just total selected)
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.prop.tsv"

# Get proportion of genes in each designation, out of all genes in a designation
# across the whole proteome (so out of the total genes, not just total selected)
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\tTotal\t"$2}' \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -a "scRNA.${SP}.${TYPE}.no_genes_per_designation_LvD.tsv" \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($3/$4)*100}' \
  >> "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.prop.tsv"

# Add colors to values
COLOR_TYPE=$(awk -vTYPE="$TYPE" 'BEGIN{if(TYPE~"^larva"){print "larva"}else if(TYPE~"^polyp"){print "polyp"}else{print "adult"} }')
cat "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.prop.tsv" \
  | python ~/scripts/add_value_to_table_SQLite3.py \
      -c 2 -d $'Total\twhite' \
      -a <(cat "../../08_Selected_OG_Analysis/public_scRNA/data/${SP}.celltype_color.tsv" \
  | awk -F'\t' -vCOLOR_TYPE="${COLOR_TYPE}" 'NR>1 && $1==COLOR_TYPE{print $7"\t"$7"\t"$8; print $5"\t"$7"\t"$8; print $3"\t"$7"\t"$8}' \
  | sort | uniq \
  | awk 'BEGIN{print "Conditions\tBroadCellType\tBroadCellColor"}{print}') \
  > "scRNA.${SP}.${TYPE}.filtered.long.no_genes_per_designation_LvD_and_cellType.prop.colors.tsv"

done
```





```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/06_Final_Classifications/expression_results/scRNA.*_LvD_and_cellType.prop.colors.tsv .
```





Get PER SEQUENCE rna-seq results

```bash
zcat ../Orthogroups.Run2.long.tsv.gz > t
for F in *.filtered.long.tsv;
do
  P=${F%*.filtered.long.tsv}
  cat t \
    | python ../../../../../02_Scripts/add_value_to_table_SQLite3.py -c 3 \
        -a <(python ../../../../../02_Scripts/groupby.py \
               -i <(awk -F'\t' '{print $1"\t"$2}' $F) --delim_groups ';' \
                      | awk -vP=$P 'BEGIN{print "sequence_id\t"P}{print}') \
        -o tt
  mv tt t
done
mv t Orthogroups.Run2.long.RNA_results.tsv
pigz Orthogroups.Run2.long.RNA_results.tsv
```

Get PER ORTHOGROUP rna-seq results

```bash
zcat ../Orthogroups.Run2.classification.tsv.gz > t
for F in *.filtered.long.tsv;
do
  P=${F%*.filtered.long.tsv}
  zcat ../Orthogroups.Run2.long.tsv.gz \
    | python ../../../../../02_Scripts/add_value_to_table_SQLite3.py -c 3 \
        -a <(python ../../../../../02_Scripts/groupby.py -i <(awk -F'\t' '{print $1"\t"$2}' $F) --delim_groups ';' \
               | awk -vP=$P 'BEGIN{print "sequence_id\t"P}{print}' \
            ) \
    | awk -F'\t' '$4!="" {print $1"\t"$4}' \
    | python ../../../../../02_Scripts/groupby.py --delim_groups ';' \
    | awk -F'\t' '{ if(NR==1){print}else{ split($2,a,";"); for(i in a){T[$1][a[i]]++}; L=""; for(i in T[$1]){L=L""i"x"T[$1][i]";"}; print $1"\t"L } }' \
    > ttt
  
  cat t \
    | python ../../../../../02_Scripts/add_value_to_table_SQLite3.py -c 1 \
        -a ttt \
        -o tt
  mv tt t
done
mv t Orthogroups.Run2.classification.RNA_results.tsv
pigz Orthogroups.Run2.classification.RNA_results.tsv
```





