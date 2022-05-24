import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[2]
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.folder as folder
import MNHN.treatment.description as description
import MNHN.treatment.stockholm as stockholm
import MNHN.treatment.capitalizer as capitalizer
import MNHN.treatment.split as split
import MNHN.treatment.pid as pid
import MNHN.treatment.redundancy as redundancy


# Séparation multi-Stockholm en un fichier Fasta par seed
list_standard_aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

path_file = sys.argv[1] # fichier Stockholm téléchargé depuis Pfam

path_folder_save = sys.argv[2] # chemin à choisir
stockholm.stockholm_separator(path_file, path_folder_save)
path_folder_stockholm = path_folder_save

path_folder_fasta = sys.argv[3]  # chemin à choisir
stockholm.multi_stockholm_to_fasta(path_folder_stockholm, path_folder_fasta)
path_data = path_folder_fasta
print(f"\nVisualisation of {path_data}")
list_residu = list_standard_aa
residu_count, total_residu, character_count, total_character = description.data_count(path_data, list_residu)
description.bar_plot_data_count(path_data, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_data, character_count, total_character , "Character")


# Capitalisation des caractères de Pfam
path_data_corrected = sys.argv[4]  # chemin à choisir
capitalizer.multi_capitalization(path_data, path_data_corrected)
path_data = path_data_corrected
print(f"\nVisualisation of {path_data}")
list_residu = list_standard_aa
residu_count, total_residu, character_count, total_character = description.data_count(path_data, list_residu)
description.bar_plot_data_count(path_data, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_data, character_count, total_character , "Character")


# Calcul des PID
path_folder_fasta = path_data_corrected # upper correction only
path_folder_pid = sys.argv[5]    # chemin à choisir
list_inclusion = list_standard_aa
pid.save_pid(path_folder_fasta, path_folder_pid, list_inclusion)


# Clustering
path_folder_fasta = path_data_corrected # still upper correction only
path_folder_fasta_nonRedondant = sys.argv[6]   # chemin à choisir
list_residu = list_standard_aa
redundancy.multi_non_redundancy_correction(path_folder_fasta, path_folder_fasta_nonRedondant, list_residu, pid_sup = 99)


# Split des seeds en seed train et seed test (50%/50% de Pfam) 
path_folder_data = path_folder_fasta_nonRedondant      # redondant corrected
result_folder = sys.argv[7]    # chemin à choisir
folder.creat_folder(result_folder)
path_folder_data_split = f"{result_folder}/Pfam_split" 
percentage_A = 50
name_data_A, name_data_B = "Pfam_train", "Pfam_test"
split.data_split(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B)