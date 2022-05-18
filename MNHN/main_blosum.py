import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))  


import MNHN.blosum.blosumfonction as blosumfonction
import MNHN.utils.folder as folder
import os


path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_10"   # seed train
name_folder_fasta =  os.path.basename(path_folder_fasta)
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

count_AA, nb_AA, count_coupleAA, nb_coupleAA = blosumfonction.multi_count_for_blosum(path_folder_fasta, path_folder_pid, 
                                                                                     list_residu, pid_inf = 62)

path_folder_Result = folder.creat_folder("/Users/pauline/Desktop/data_Result/Pfam_split/NonContextual_Result")
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