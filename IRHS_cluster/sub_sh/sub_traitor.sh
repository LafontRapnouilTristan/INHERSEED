#!/bin/bash
#SBATCH --job-name=traitor_inherseed
#SBATCH --output=traitor_inherseed-%N-%j.out
#SBATCH --error=traitor_inherseed-%N-%j.err
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=05:00:00
#SBATCH --mail-user=tristan.lafontrapnouil@inrae.fr
#SBATCH --mail-type=ALL

# To be launched wtih sbatch scripts/sub_traitor.sh from a directory containing /images and /scripts with the current file in /scripts 
# Load conda
echo "load conda"
source /usr/local/applis/conda/etc/profile.d/conda.sh

# Environment
echo "activate environment"
conda activate /mnt/home2/tlafontrapn/tests/conda_envs/envs/venv_traitor

# Test env
echo "test traitor installation"
traitor -h

# Run extraction
echo "run extraction"
traitor extract -i "./images" -o "./images_extracted" -u -b --rm_bg -p 10

# Run alignement
echo "run alignement"
traitor align -i "./images" -m "./images_extracted" -o "./images_aligned"

# Run measurements
echo "run measurements"
traitor measure -i "./images_aligned" -o "./objects_measurements.csv" -n 8 -t 8

echo "end script"
