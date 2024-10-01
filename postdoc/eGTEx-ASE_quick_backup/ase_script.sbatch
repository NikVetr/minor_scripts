#!/bin/bash

#SBATCH --job-name=rstan_parallel  # Job name
#SBATCH --cpus-per-task=128         # Number of total CPU cores
#SBATCH --partition=batch          # Partition (billed)
#SBATCH --account=smontgom         # SLURM account (billing)
#SBATCH --time=99:00:00            # Time limit (adjust as needed)
#SBATCH --mem-per-cpu=2G           # Memory per CPU core
#SBATCH --output=output_%j.log     # Standard output log
#SBATCH --error=error_%j.log       # Error log

# Load R and required modules
module load R/4.2.2                # Adjust to match the available R version on SCG

# Define output directory
OUTPUT_DIR="/oak/stanford/groups/smontgom/nikgvetr/output"
mkdir -p $OUTPUT_DIR               # Create output directory if it doesn't exist

# Run the R script
Rscript ase_script_rstan_GTEx.R
