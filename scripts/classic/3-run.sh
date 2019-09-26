#!/bin/bash

# Submit one job per library, in which the job is to map the fastq files to a 
# reference with CellRanger and call the feature-barcode matrices.
# $1 = project name (e.j BCLLATLAS), $2 = subproject name (e.j. BCLLATAS_10)
jobs_path="/scratch/devel/rmassoni/$1/$2/jobs"
echo $jobs_path
cd "$jobs_path"
for dir in ./*; do
  cd "$dir"
  mnsubmit "${dir}.cmd"
  cd "$jobs_path"
done


