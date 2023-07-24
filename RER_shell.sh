#!/bin/bash
#SBATCH --job-name=rerconverge    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --output=rerconverge.log   # Standard output and error log
pwd; hostname; date

conda activate phangorn

Rscript RER_revision_analysis.R
