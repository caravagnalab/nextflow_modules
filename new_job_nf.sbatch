#!/bin/sh

#SBATCH --job-name=NF_SEQUENZA
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00

module load singularity
module load java/19.0.1


export TMP_DIR=/orfeo/LTS/LADE/LT_storage/lvaleriani/tmp

config=/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow_modules/modules/config/nextflow_config_slurm

/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow/nextflow run /orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow_modules/main.nf -c $config \
--samples small_input.csv

