#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=1-00:00:00
#SBATCH --output=log/main_%j
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load mamba
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run analysis.nf \
    -params-file LmrCD_analysis.json \
    -resume \
    -with-report \
    -with-timeline \
    -ansi-log false
