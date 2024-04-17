# Explore deep homology of refined orthofinder clusters

Following: https://elifesciences.org/articles/67667#s4

## Setup analysis directory

Link data files.

```bash
ln -s ../../../../02_Scripts/grepf_fasta.py
ln -s ../../../../02_Scripts/grepf_column.py

ln -s ../06_Final_Classifications/Orthogroups.Run2.long.tsv.gz
```

Setup bash environment.

```bash
conda activate py27
```



## Install hh-suite

```bash
git clone https://github.com/soedinglab/hh-suite.git
cd hh-suite && mkdir build && cd build
cmake -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx ..
make

INSTALL_BASE_DIR=$(dirname $PWD)
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_BASE_DIR" ..
make install
```

## Run hh-suite analysis

Extract cluster IDs for later processing.

```bash
./run_01_extract_ids.sh
```

Extract cluster sequences.

```bash
./run_02_extract_seqs.sh
```

Align seqs in cluster.

```bash
./run_03_align_seqs.sh
```

NOTE: There looks to be 93 sequences/OGs which `mafft` ignores becuase there are only 1 sequence in the input. We need to cat the sequences for these clusters into the alignment file (jast as `mafft` does for the others) to prevent errors in downstream analysis.

Get list of problem alignment files and `cat` the cluster fasta files to the ends of these alignment files.

```bash
./run_04_correct_single_seq_OGs.sh
```

Check that all proteins in the total protein file are represented in the alignment file.

```bash
./run_05_check_OG_seqs.sh
```

Run `ffindex_build` to combine alignments into a single 

```bash
./run_06_ffindex_build.sh
```

Compress cluster files.

```bash
./run_07_compress_clusters.sh
```

Run `hhconsensus` to build a consensus sequence per alignment

```bash
./run_08_hhconsensus.sh
```

Run `hhmake` to build a hmm from the alignment + consensus sequence info.

```bash
./run_09_hhmake.sh
```

Run `cstranslate` to translate a sequence/alignment into an abstract state alphabet (whatever that means)

```bash
./run_10_cstranslate.sh
```

Run `ffindex_order` to make sure the database entries are ordered correctly

```bash
./run_11_ffindex_order.sh
```



Split into 5000 parst so that we can run it on the Amarel server.

```bash
mkdir -p hhdb_a3m.split
cut -f1 hhdb_a3m.ffindex | split -a 4 --numeric-suffixes=1 -n r/5000 - hhdb_a3m.split/
for N in {0001..5000};
do
  echo $N
  ffindex_order hhdb_a3m.split/${N} hhdb_a3m.ffdata hhdb_a3m.ffindex hhdb_a3m.split/${N}_a3m.ffdata hhdb_a3m.split/${N}_a3m.ffindex
done
```



Run a self-vs-self comparison of the database using `hhblits` on the Coral server.

```bash
./run_09_hhblits_self-vs-self.sh
```



Copy files to Amarel and run `hhblits` on there.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/30_OG_Deep_Homology_Search/hhdb_* .
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/30_OG_Deep_Homology_Search/run_12_hhblits_self-vs-self_Amarel.sh .
```

Run each of the 5000 parts as separate jobs on amarel.

```bash
run_12_hhblits_self-vs-self_Amarel.sh
```

Push from Amarel to Coral.

```bash
rsync -avzP * timothy@coral.rutgers.edu:/scratch/timothy/projects/0032_Coral_genomes_analysis/03_Analysis/02_Coral_genomic_analysis/2023-06-09/30_OG_Deep_Homology_Search/
```



Convert hhr to BLAST-like format.

```bash
./run_13_hhblits_hhr_reformat.sh
```

Filter BLAST-like results.

```bash
./run_14_hhblits_hhr_filter.sh
```



Combine filtered BLAST-like results.

```bash
./run_15_combine_results.sh
```

Compress split results.

```bash
./run_16_compress_results.sh
```

