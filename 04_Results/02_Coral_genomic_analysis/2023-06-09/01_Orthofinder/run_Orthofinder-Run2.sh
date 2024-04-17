#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

ulimit -n 100000

export PATH="$PATH:/home/timothy/programs/OrthoFinder_v2.5.4"
export PATH="$PATH:/home/timothy/programs/OrthoFinder_v2.5.4/bin"

T_NCPUS=12
A_NCPUS=12

#### Start Script
run_cmd "md5sum files2cluster/*.pep.faa | tee md5sum_list-Run2.txt"
run_cmd "orthofinder -t $T_NCPUS -a $A_NCPUS -b files2cluster/OrthoFinder/Results_Jul14/WorkingDirectory"


