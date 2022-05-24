import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]
sys.path.append(str(package_root_directory_MNHN))  


import MNHN.blosum.blosumfonction as blosumfonction
import MNHN.utils.folder as folder
import os


path_folder_fasta = sys.argv[1]   # chemin des seeds d'entrainement
name_folder_fasta =  os.path.basename(path_folder_fasta)
path_folder_pid = sys.argv[2]  # chemin des pid par seeds
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

count_AA, nb_AA, count_coupleAA, nb_coupleAA = blosumfonction.multi_count_for_blosum(path_folder_fasta, path_folder_pid, 
                                                                                     list_residu, pid_inf = 62)

path_non_contextual_result = sys.argv[3]  # chemin à choisir
path_folder_Result = folder.creat_folder(path_non_contextual_result)
freq_AA, freq_coupleAA = blosumfonction.freq_for_blosum(count_AA, nb_AA, count_coupleAA, nb_coupleAA, path_folder_Result)


# blosum
blosum = blosumfonction.blosum_score(freq_AA, freq_coupleAA, path_folder_Result, scale_factor = 2)
blosumfonction.blosum_visualisation(blosum)

matrix_diff, pid_inf_ref, average_euclidean_d, average_diff = blosumfonction.blosum_difference(blosum, pid_inf_ref = 62)
title_heatmap = f"Heatmap of Blosum({name_folder_fasta}) - Blosum{pid_inf_ref}Ref:\nmean difference \
                 = {average_diff}, mean euclidean distance = {average_euclidean_d}"
blosumfonction.blosum_heatmap(blosum, path_folder_Result, title_heatmap, size_annot = 5)


# conditional proba
cond_proba = blosumfonction.blosum_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result)
blosumfonction.sum_line(cond_proba)
title_heatmap = f"Heatmap of the conditional probability matrix \n computed on {name_folder_fasta}"
blosumfonction.blosum_heatmap(cond_proba, path_folder_Result, title_heatmap)