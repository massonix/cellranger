# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files (2 links, as it is pair-end). If it one
# library is sequence in multiple flowcells/lanes, it concatenates the fastqs into a single file.
# 2. output: dir which contains the files with the (1) standard output and (2) standard error of cellranger
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger

# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re

# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--reference",
		    dest = "reference",
		    action = "store",
		    default = None,
		    help = "Reference genome to use (human or mouse")
parser.add_argument("--subproject",
		    dest = "subproject",
		    action = "store",
		    default = None,
		    help = "Subproject we are working on (i.e. BCLLATLAS_10") 
options = parser.parse_args()

# Read the lims output table and the hashtag references
lims = pd.read_csv("info.txt", sep = "\t", header = 0)
hash_xls = pd.ExcelFile("{}_hashing_summary.xlsx".format(options.subproject))
print("Successfully read!")
print(lims.head())

# Define important paths and directories
fastq_path = "/project/production/fastq"
reference = options.reference
if reference == "human":
    ref_path = "/scratch/devel/rmassoni/reference/human/refdata-cellranger-GRCh38-3.0.0/"
elif reference == "mouse":
    ref_path = "/scratch/devel/rmassoni/reference/mouse/refdata-cellranger-mm10-3.0.0"
if not os.path.exists("jobs"):
    os.mkdir("jobs")
print("Reference path is {}".format(ref_path))

# Create cellranger-friendly symmlinks to fastq files
sample_ids = np.unique(lims["SampleName"])
for iden in sample_ids:
    print("Current library is {}".format(iden))

    # Define and create directories
    cwd = os.getcwd()
    regex = re.compile("(_cDNA$|_HTO$)")
    iden_dir = regex.sub("", iden)
    jobs_dir = "{}/jobs/{}".format(cwd, iden_dir)
    fastq_dir = "{}/fastq".format(jobs_dir)
    output_dir = "{}/output".format(jobs_dir)
    dirs_to_create = [jobs_dir, fastq_dir, output_dir]
    print("The directories to create are {}".format(dirs_to_create))
    for new_dir in dirs_to_create:
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
    bool_mask = (lims["SampleName"] == iden) & (lims["libraryPassFail"] != "fail") & (lims["LanePassFail"] != "fail")
    lims_sub = lims.loc[bool_mask]
    print(lims_sub)
    library = lims_sub.loc[lims_sub.index[0], "library"]
    index = lims_sub.loc[lims_sub.index[0], "index"]
    if "IlluminaTruSeq" in index:
        library_type = "HTO"
        library_type_feat = "Antibody Capture"
    else:
        library_type = "cDNA"    
        library_type_feat = "Gene Expression"
    fastq_subdir = "{}/{}".format(fastq_dir, library_type)
    if not os.path.exists(fastq_subdir):
        os.mkdir(fastq_subdir)

    # Fastq files are under /project/production/fastq/FC/Lane/fastq/FC_Lane_Index_1/2.fastq.gz
    # The CellRanger convention establishes that the symlinks should be named: lib_S1_L00(Lane)_R1/2_001.fastq.gz
    count_id = np.sum(lims_sub["SampleName"] == iden)
    if count_id > 1:
        # Concatenates fastq files from same sample but in different flow cells and/or lanes
        fc_lane_list = [(lims_sub["flowcell"][x], lims_sub["lane"][x]) for x in lims_sub.index]
        print(fc_lane_list)
        fastq_path_list_r1 = ["{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index) for fc, lane in fc_lane_list]
        print(fastq_path_list_r1)
        fastq_path_list_r2 = ["{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index) for fc, lane in fc_lane_list]
        print(fastq_path_list_r2)
        fastq_path_list_r1.insert(0, "cat")
        fastq_path_list_r2.insert(0, "cat")
        subprocess.run(fastq_path_list_r1, stdout = open("{}/{}_S1_L001_R1_001.fastq.gz".format(fastq_subdir, iden), "w"))
        subprocess.run(fastq_path_list_r2, stdout = open("{}/{}_S1_L001_R2_001.fastq.gz".format(fastq_subdir, iden), "w"))
    else:
        # As this sample is only present in one flowcell and lane, there is no need to concatenate fastq
        # Create symlinks
        fc = lims_sub["flowcell"].values[0]
        lane = lims_sub["lane"].values[0]
        fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
        fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
        subprocess.run(["ln", "-s", fastq_path_r1, "{}/{}_S1_L001_R1_001.fastq.gz".format(fastq_subdir, iden_dir)])
        subprocess.run(["ln", "-s", fastq_path_r2, "{}/{}_S1_L001_R2_001.fastq.gz".format(fastq_subdir, iden_dir)]) 

    # Create the libraries.csv file, which will specifies cellranger the fastq directory and the type of library (HTO or cDNA)
    if not os.path.exists("jobs/{}/libraries.csv".format(iden_dir)):
        lib_csv = open("{}/libraries.csv".format(jobs_dir), "w")
        lib_csv_str = "fastqs,sample,library_type\n{},{},{}\n".format(fastq_subdir, iden, library_type_feat)
        lib_csv.write(lib_csv_str)
        lib_csv.close()
    else:
        lib_csv = open("{}/libraries.csv".format(jobs_dir), "a")
        lib_csv_str = "{},{},{}".format(fastq_subdir, iden, library_type_feat)
        lib_csv.write(lib_csv_str)
        lib_csv.close()
    
    # Create the feature reference csv file to identify each hashtag with each experimental condition
    hash_xls_iden = hash_xls.parse(iden_dir)
    hash_xls_iden.to_csv("{}/feature_reference.csv".format(jobs_dir), sep = ",", index_label = False, index = False)
    
    # Create job script with the call to cellranger
    job_script_file = open("{}/{}.cmd".format(jobs_dir, iden_dir), "w")
    job_script = """#!/bin/bash 

# @ initialdir = . 
# @ error = ./output/{}.err 
# @ output = ./output/{}.out 
# @ cpus_per_task = 12 
# @ wall_clock_limit = 16:00:00 

module load PYTHON/2.7.5 
module load lims/1.2 

/scratch/devel/rmassoni/cellranger-3.0.2/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id {} --chemistry SC3Pv3 --expect-cells 5000 --localcores 12 --localmem 64 --transcriptome {};
    """.format(iden_dir, iden_dir, iden_dir, ref_path)
    job_script_file.write(job_script)
    job_script_file.close()


