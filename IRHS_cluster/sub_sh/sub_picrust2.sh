#!/bin/bash
#SBATCH --job-name=picrust2_mapp
#SBATCH --output=picrust2_mapp-%N-%j.out
#SBATCH --error=picrust2_mapp-%N-%j.err
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=tristan.lafontrapnouil@inrae.fr
#SBATCH --mail-type=ALL

# To be launched wtih sbatch scripts/sub_traitor.sh from a directory containing /images and /scripts with the current file in /scripts 
# Load conda
echo "load conda"
source /usr/local/applis/conda/etc/profile.d/conda.sh

# Environment
echo "activate environment"
conda activate /mnt/home2/tlafontrapn/tests/conda_envs/envs/venv_picrust2
conda activate venv_picrust2

# Test env
echo "test picrust2 installation"
pytest

# Run full pipeline
echo "run full pipe"
picrust2_pipeline.py -s data/fasta4picrust2.fasta -i data/biom4picrust2.biom -o outs/picrust2_out_pipeline_split -p 1

echo "end script"
