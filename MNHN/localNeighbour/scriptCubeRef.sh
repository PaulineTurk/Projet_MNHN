#!/bin/bash

# Job name
#SBATCH --job-name=cube

# Number pf nodes
#SBATCH --nodes=1

# Number of processor per node
#SBATCH --ntasks-per-node=16

# Nomber of RAM per node
#SBATCH --mem=200Go

# Type of machines requested (type1 or 2)
#SBATCH --partition=type_2

# Name of output file
#SBATCH --output=cube.out

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
echo -e "\nStart calcul cubes\n"

path_folder_pID=(/mnt/beegfs/pturk/PID_couple)           # à changer qd ils seront dans le cluster
path_folder_data_split=(/mnt/beegfs/pturk/PfamSplit_50)  # à changer qd ils seront dans le cluster
path_new_folder=(/mnt/beegfs/pturk/Cubes_10_Pfam_50_A)   # à changer qd ils seront dans le cluster

cd /trinity/home/pturk
python ./script_cubes.py $path_folder_pID $path_folder_data_split $path_new_folder
date

echo -e "\nDone\n"
