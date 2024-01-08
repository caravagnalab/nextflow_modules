#!/bin/sh

#SBATCH --job-name=NF_CNAqc
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00

module load singularity
module load java/19.0.1

/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow/nextflow run /orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow_modules/main.nf \
 -c /orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow_modules/modules/config/nextflow_config_slurm \
 --samples small_input.csv

