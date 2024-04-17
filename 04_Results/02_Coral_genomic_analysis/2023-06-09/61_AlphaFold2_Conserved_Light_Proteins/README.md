# Run AlphaFold2 on M. capitata Conserved Light Proteins


Link to *M. capitata* proteins and BUSCO protein-mode results (used to pick single copy conserved proteins)
```bash
ln -s /scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/02_busco/pep.faa.busco_eukaryota_odb10.tar.gz
ln -s /scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Montipora_capitata_KBHIv3/02_busco/Montipora_capitata_KBHIv3.genes.pep.faa
```

Unzip BUSCO results.
```bash
tar -zxvf pep.faa.busco_eukaryota_odb10.tar.gz
```

Extract "Complete" proteins as dictated by BUSCO.
```bash
~/miniforge3/envs/py27/bin/python ./../../../../02_Scripts/grepf_fasta.py \
  -f <(cat pep.faa.busco_eukaryota_odb10/run_eukaryota_odb10/full_table.tsv | awk -F'\t' '$2=="Complete" {print $3}') \
  -i Montipora_capitata_KBHIv3.genes.pep.faa \
  -o proteins.faa
```

Get each protein in its own file to make downstream processing easier.
```bash
mkdir proteins

cat proteins.faa \
  | grep '>' \
  | sed -e 's/>//' \
  | while read NAME;
      do
      ~/miniforge3/envs/py27/bin/python ./../../../../02_Scripts/grepf_fasta.py \
        -f <(echo "$NAME") \
        -i proteins.faa \
        -o "proteins/${NAME}.faa"
done

ls -1 proteins/*.faa > proteins.parts
```

Run DIAMOND BLASTP on protein sequences against UniProt SwissProt to assign gene names + IDs.
```bash
./run_diamond_UniProt_SwissProt.sh
```

Get top SwissProt hit as gene annotation.
```bash
zcat proteins.faa.diamond_blastp_UniProt_SwissProt.outfmt6.gz \
  | awk -F'\t' 'BEGIN{print "sequence_id\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle"}!seen[$1]++{print}' \
  | awk -F'\t' '{print $15"\t"$0}' \
  | sed -e 's/sp|[^ ]* //' \
  | awk 'BEGIN{OFS=FS="\t"}{print $0"\t"$1}' \
  | cut -f2-15,17 \
  > proteins.parts.diamond_blastp_UniProt_SwissProt.tophit.outfmt6
```

Extract OrthoGroup info from previous analysis.
```bash
zcat ../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz \
  | python ~/scripts/grepf_column.py --keep_header -c 3 -f <(cat proteins.parts | sed -e 's@.*/@@' -e 's/.faa//') \
  > proteins.parts.long.tsv
```

Extract OrthoGroup classification.
```bash
zcat ../06_Final_Classifications/Orthogroups.Run2.classification.tsv.gz \
  | python ~/scripts/grepf_column.py --keep_header -c 1 -f <(cut -f1 proteins.parts.long.tsv) \
  > proteins.parts.classification.tsv
```

Combine Orthogroup results together.
```bash
cat proteins.parts.long.tsv \
  | python ~/scripts/add_value_to_table.py -c 1 -a proteins.parts.classification.tsv \
  | python ~/scripts/add_value_to_table.py -c 3 -a proteins.parts.diamond_blastp_UniProt_SwissProt.tophit.outfmt6 \
  > proteins.parts.combined.tsv
```

Run ColabFold on Amarel Server.
```bash
./run_ColabFold_Search.sh

# Check it finished correctly
grep -L 'ExitStatus: 0' ColabFold_Search.slurm_out.* | sed -e 's/.*00.//' | xargs printf "%s,"
grep -L 'ExitStatus: 0' ColabFold_Search.slurm_out.* | xargs rm

grep -l 'FAILED_PRECONDITION\|from jax._src import dispatch' ColabFold_Batch.slurm_out.* > t
cat t | sed -e 's/.*00.//' | xargs printf "%s,"
cat t | sed -e 's/.*00.//' | while read N; do awk -vN=$N 'NR==N {print $1".colabfold_batch.done"}' proteins.parts; done | xargs ls
cat t | xargs rm

mv proteins.parts2run proteins.parts2run.old
for F in `cat proteins.parts2run.old`; do if [[ ! -e "$F.colabfold_batch.done" ]]; then echo $F; fi; done > proteins.parts2run
```


Search best structure aginst PDB datbase.
```bash
./run_ColabFold_PDB_search.sh
```

Combine results
```bash
for F in proteins/*.faa.colabfold/plots/0.pdb;
do
  cat ${F%*.pdb}.PDB_search.fmt4.tsv \
    | sed -e "s@0.pdb@$F@"
done \
  | awk -F'\t' '!seen[$1]++' \
  | awk -F'\t' '{ if(NR==1){print $0"\tqcovpct\ttcovpct"}else{print $0"\t"( ($8-($7-1)) / $13 )*100"\t"( ($10-($9-1)) / $14 )*100} }' \
  > proteins.faa.PDB_search.fmt4.tsv

cat proteins/*.faa.colabfold/plots/0.confidence.tsv | awk -F'\t' '!seen[$1]++' > proteins.faa.confidence.tsv
```

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/61_AlphaFold2_Conserved_Light_Proteins/proteins.parts.diamond_blastp_UniProt_SwissProt.tophit.outfmt6 .
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/61_AlphaFold2_Conserved_Light_Proteins/proteins.faa.PDB_search.fmt4.tsv .
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/61_AlphaFold2_Conserved_Light_Proteins/proteins.faa.confidence.tsv .
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/61_AlphaFold2_Conserved_Light_Proteins/./proteins/*/plots .
```










