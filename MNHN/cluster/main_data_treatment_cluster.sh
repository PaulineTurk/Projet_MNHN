#!/bin/bash

# Job name
#SBATCH --job-name=Pfam_treatment

# Number pf nodes
#SBATCH --nodes=1

# Number of processor per node
#SBATCH --ntasks-per-node=16

# Nomber of RAM per node
#SBATCH --mem=200Go

# Type of machines requested (type1 or 2)
#SBATCH --partition=type_2

# Name of output file
#SBATCH --output=treatment.out

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

path_file=(/mnt/beegfs/pturk/Pfam_Sample)                               # chemin fichier Stockholm téléchargé depuis Pfam
path_folder_save=(/mnt/beegfs/pturk/Pfam_Stockholm)                     # chemin a fixer ou mettre les fichier stockholms (1/seed)
path_folder_fasta=(/mnt/beegfs/pturk/Pfam_Fasta)                        # chemin ou mettre les fichiers fasta (1/seed)
path_data_corrected=(/mnt/beegfs/pturk/Pfam_Upper)                      # chemin ou mettre les fichiers fasta avec tous les caractères en majuscule
path_folder_pid=(/mnt/beegfs/pturk/PID)                                 # chemin ou mettre les numpy de pid de chaque seed
path_folder_fasta_nonRedondant=(/mnt/beegfs/pturk/Pfam_nonRedondant)    # chemin ou mettre les seed après clustering et retrait de redondance
result_folder=(/mnt/beegfs/pturk/data_Result)                           # chemin ou mettre Pfam splittée à 50/50 par rapport aux seeds entre train/test


cd /trinity/home/pturk/Projet_MNHN/MNHN/
python ./main_data_treatment_cluster.py $path_file $path_folder_save $path_folder_fasta $path_data_corrected $path_folder_pid $path_folder_fasta_nonRedondant $result_folder
date

echo -e "\nDone\n"
