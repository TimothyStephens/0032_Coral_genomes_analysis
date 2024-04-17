# Plot QC stats against orthogroup tree

Plot `busco` stats and genome assembly stats against the orthogroup tree from `orthofinder`.

## Check stats files

Check stats files md5sum's and that the order of the rows are the same across files. This will help check that when we combine everything it will be in the correct order.

```bash
# GeneStats
while read LINE;
do 
  ID=$( echo -e "$LINE" | awk -F'\t' '{print $1}')
  DIR=$(echo -e "$LINE" | awk -F'\t' '{print $4}')
  md5sum "../${DIR}/01_stats/${ID}.GeneStats.tsv"
done < <(awk -F'\t' '$3=="genome" || $3=="transcriptome"' ../samples.txt) > md5sum_list.GeneStats.txt

# Genome BUSCO
for LINEAGE in "genome.fa.busco_eukaryota_odb10" "genome.fa.busco_metazoa_odb10" "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  while read F;
  do
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    md5sum "../${DIR}/02_busco/${LINEAGE}.results.txt"
  done < <(awk -F'\t' '$3=="genome"' ../samples.txt)
done > md5sum_list.BUSCO_Genomes.txt

# Transcriptome BUSCO
for LINEAGE in "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    md5sum "../${DIR}/02_busco/${LINEAGE}.results.txt"
  done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt)
done > md5sum_list.BUSCO_Transcriptomes.txt
```





```bash
# GeneStats
while read LINE;
do 
  ID=$( echo -e "$LINE" | awk -F'\t' '{print $1}')
  DIR=$(echo -e "$LINE" | awk -F'\t' '{print $4}')
  FILE="../${DIR}/01_stats/${ID}.GeneStats.tsv"
  echo -e $(awk -F'\t' '{print $1}' "${FILE}" | md5sum | cut -d' ' -f1)"\t$FILE"
done < <(awk -F'\t' '$3=="genome"' ../samples.txt) > md5sum_header.GeneStats_Genomes.txt
while read LINE;
do 
  ID=$( echo -e "$LINE" | awk -F'\t' '{print $1}')
  DIR=$(echo -e "$LINE" | awk -F'\t' '{print $4}')
  FILE="../${DIR}/01_stats/${ID}.GeneStats.tsv"
  echo -e $(awk -F'\t' '{print $1}' "${FILE}" | md5sum | cut -d' ' -f1)"\t$FILE"
done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt) > md5sum_header.GeneStats_Transcriptomes.txt

# BUSCO
for LINEAGE in "genome.fa.busco_eukaryota_odb10" "genome.fa.busco_metazoa_odb10" "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    FILE="../${DIR}/02_busco/${LINEAGE}.results.txt"
    echo -e $(awk -F'\t' '{print $1}' "${FILE}" | md5sum | cut -d' ' -f1)"\t$FILE"
  done < <(awk -F'\t' '$3=="genome"' ../samples.txt)
done > md5sum_header.BUSCO_Genomes.txt
for LINEAGE in "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    FILE="../${DIR}/02_busco/${LINEAGE}.results.txt"
    echo -e $(awk -F'\t' '{print $1}' "${FILE}" | md5sum | cut -d' ' -f1)"\t$FILE"
  done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt)
done > md5sum_header.BUSCO_Transcriptomes.txt
```



## Combine stats files

Combine genome stats files.

```bash
# Use first entry in list to create headers for combined file
ID=$( awk -F'\t' '$3=="genome"{print $1}' ../samples.txt | head -n1)
DIR=$(awk -F'\t' '$3=="genome"{print $4}' ../samples.txt | head -n1)
awk -F'\t' 'BEGIN {L="SampleID"} $1!~"^#" && $0!=""{L=L"\t"$1} END {print L} ' \
    "../${DIR}/01_stats/${ID}.GeneStats.tsv" > "all_genomes-01_stats-results.tsv"

# For each "genome" dataset print stats
while read LINE;
do 
  ID=$( echo -e "$LINE" | awk -F'\t' '{print $1}')
  DIR=$(echo -e "$LINE" | awk -F'\t' '{print $4}')
  awk -F'\t' 'BEGIN {L=""} $1!~"^#" && $0!=""{
      if(L=="") {
        L=$2
      } else {
        L=L"\t"$2
      }
      if($1=="Genome file") {
        gsub(".assembly.fasta", "", $2); 
        L=$2"\t"L
      }
    } END {print L}' "../${DIR}/01_stats/${ID}.GeneStats.tsv" \
      >> "all_genomes-01_stats-results.tsv"
done < <(awk -F'\t' '$3=="genome"' ../samples.txt)
```

Combine transcriptome stats files.

```bash
# Use first entry in list to create headers for combined file
ID=$( awk -F'\t' '$3=="transcriptome"{print $1}' ../samples.txt | head -n1)
DIR=$(awk -F'\t' '$3=="transcriptome"{print $4}' ../samples.txt | head -n1)
awk -F'\t' 'BEGIN {L="SampleID"} $1!~"^#" && $0!=""{L=L"\t"$1} END {print L} ' \
    "../${DIR}/01_stats/${ID}.GeneStats.tsv" > "all_transcriptomes-01_stats-results.tsv"

# For each "transcriptome" dataset print stats
while read LINE;
do 
  ID=$( echo -e "$LINE" | awk -F'\t' '{print $1}')
  DIR=$(echo -e "$LINE" | awk -F'\t' '{print $4}')
  awk -F'\t' 'BEGIN {L=""} $1!~"^#" && $0!=""{
      if(L=="") {
        L=$2
      } else {
        L=L"\t"$2
      }
      if($1=="CDS file") {
        gsub(".transcripts.cds.fna", "", $2); 
        L=$2"\t"L
      }
      if($1=="PEP file") {
        gsub(".transcripts.pep.faa", "", $2); 
        L=$2"\t"L
        E="\tNA"
      }
    } END {print L""E}' "../${DIR}/01_stats/${ID}.GeneStats.tsv" \
      >> "all_transcriptomes-01_stats-results.tsv"
done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt)
```

## Combine `busco` results

Combine genome `busco` results files.

```bash
for LINEAGE in "genome.fa.busco_eukaryota_odb10" "genome.fa.busco_metazoa_odb10" "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  # Use first entry in list to create headers for combined file
  ID=$( awk -F'\t' '$3=="genome"{print $1}' ../samples.txt | head -n1)
  DIR=$(awk -F'\t' '$3=="genome"{print $4}' ../samples.txt | head -n1)
  awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' \
    "../${DIR}/02_busco/${LINEAGE}.results.txt" \
    > "all_genomes-02_busco-${LINEAGE}-results.tsv"
  
  # For each "genome" dataset print stats
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    awk -F'\t' -v ID="$ID" 'BEGIN {L=ID} {gsub("%", "", $3); L=L"\t"$3} END {print L}' \
      "../${DIR}/02_busco/${LINEAGE}.results.txt" \
      >> "all_genomes-02_busco-${LINEAGE}-results.tsv"
  done < <(awk -F'\t' '$3=="genome"' ../samples.txt)
done
```

Combine transcriptome `busco` results files.

```bash
for LINEAGE in "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  # Use first entry in list to create headers for combined file
  ID=$( awk -F'\t' '$3=="transcriptome"{print $1}' ../samples.txt | head -n1)
  DIR=$(awk -F'\t' '$3=="transcriptome"{print $4}' ../samples.txt | head -n1)
  awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' \
    "../${DIR}/02_busco/${LINEAGE}.results.txt" \
    > "all_transcriptomes-02_busco-${LINEAGE}-results.tsv"
  
  # For each "transcriptome" dataset print stats
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    awk -F'\t' -v ID="$ID" 'BEGIN {L=ID} {gsub("%", "", $3); L=L"\t"$3} END {print L}' \
      "../${DIR}/02_busco/${LINEAGE}.results.txt" \
      >> "all_transcriptomes-02_busco-${LINEAGE}-results.tsv"
  done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt)
done
```

## Combine stats tables.

Combine `busco` results for the different datasets/modes together into a single matrix. Take only the % and no. of genes recovered in each category from results file.

Combine genomes `busco` results.

```bash
for LINEAGE in "genome.fa.busco_eukaryota_odb10" "genome.fa.busco_metazoa_odb10" "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  # Use first entry in list to create headers for combined file
  ID=$( awk -F'\t' '$3=="genome"{print $1}' ../samples.txt | head -n1)
  DIR=$(awk -F'\t' '$3=="genome"{print $4}' ../samples.txt | head -n1)
  awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' \
    "../${DIR}/02_busco/${LINEAGE}.results.txt" \
    > "all_genomes-02_busco-${LINEAGE}-results4table.tsv"
  
  # For each "genome" dataset print stats
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    awk -F'\t' -v ID="$ID" 'BEGIN {L=ID} {L=L"\t"$3" ("$2")"} END {print L}' \
      "../${DIR}/02_busco/${LINEAGE}.results.txt" \
      >> "all_genomes-02_busco-${LINEAGE}-results4table.tsv"
  done < <(awk -F'\t' '$3=="genome"' ../samples.txt)
done
```

Combine transcriptomes `busco` results.

```bash
for LINEAGE in "pep.faa.busco_eukaryota_odb10" "pep.faa.busco_metazoa_odb10";
do
  # Use first entry in list to create headers for combined file
  ID=$( awk -F'\t' '$3=="transcriptome"{print $1}' ../samples.txt | head -n1)
  DIR=$(awk -F'\t' '$3=="transcriptome"{print $4}' ../samples.txt | head -n1)
  awk -F'\t' 'BEGIN {L="SampleID"} {L=L"\t"$1} END {print L} ' \
    "../${DIR}/02_busco/${LINEAGE}.results.txt" \
    > "all_transcriptomes-02_busco-${LINEAGE}-results4table.tsv"
  
  # For each "transcriptome" dataset print stats
  while read F;
  do 
    ID=$( echo -e "$F" | awk -F'\t' '{print $1}')
    DIR=$(echo -e "$F" | awk -F'\t' '{print $4}')
    awk -F'\t' -v ID="$ID" 'BEGIN {L=ID} {L=L"\t"$3" ("$2")"} END {print L}' \
      "../${DIR}/02_busco/${LINEAGE}.results.txt" \
      >> "all_transcriptomes-02_busco-${LINEAGE}-results4table.tsv"
  done < <(awk -F'\t' '$3=="transcriptome"' ../samples.txt)
done
```

