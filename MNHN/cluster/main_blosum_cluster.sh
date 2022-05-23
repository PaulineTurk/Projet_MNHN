#!/bin/bash

# Job name
#SBATCH --job-name=Blosum_no_context

# Number pf nodes
#SBATCH --nodes=1

# Number of processor per node
#SBATCH --ntasks-per-node=16

# Nomber of RAM per node
#SBATCH --mem=200Go

# Type of machines requested (type1 or 2)
#SBATCH --partition=type_2

# Name of output file
#SBATCH --output=Blosum_no_context.out

# Calculation times
#SBATCH --time=1-00:00:00

# Add email
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pauline.turk@agroparistech.fr


# loading modules
module load userspace/tr17.10
module load python/3.6.3




########### script ##########
date
echo -e "\nStart Pfam treatment\n"

path_folder_fasta=()                # chemin seed d'apprentissage
path_folder_pid=()      # chemin des pid de seed
path_non_contextual_result=() # chemin ou mettre les résultats

cd /trinity/home/pturk/Projet_MNHN/MNHN/
python ./main_blosum_cluster.py $path_folder_fasta $path_folder_pid $path_non_contextual_result
date

echo -e "\nDone\n"