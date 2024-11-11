# DMS-pipeline
DMS-pipeline is a Nextflow pipeline for processing deep mutational scanning (DMS) data. This pipeline implements the sample processing and analysis by [Gianmarco Meier et al.](https://www.nature.com/articles/s41589-022-01205-1) For the original code, please see Gianmarco's repo [DMS_ABC](https://github.com/giameier/DMS_ABC).

## Running the pipeline
1. Clone the repository.
2. Create or edit a new parameter file e.g. `LmrCD.json`.
3. Edit the `run_DMS-pipeline.slurm` script and referer to the respective parameter file.
4. Run the pipeline on the cluster:
`sbatch run_DMS-pipeline.slurm`

> [!NOTE]
> You can also perform a dry-run of the pipeline locally with `./run_DMS-pipeline.stub`.
