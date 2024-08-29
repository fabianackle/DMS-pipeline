#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=1-00:00:00
#SBATCH --output=log/main_%j

module load mamba
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run main.nf