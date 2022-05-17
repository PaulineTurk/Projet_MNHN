import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))  


import MNHN.blosum.matrix as blosumCode
from MNHN.utils.folder import creatFolder
import os


path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_10"   # seed train
name_folder_fasta =  os.path.basename(path_folder_fasta)
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

count_AA, nb_AA, count_coupleAA, nb_coupleAA = blosumCode.multicountBlosum(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62)

path_folder_Result = creatFolder("/Users/pauline/Desktop/data_Result/Pfam_split/NonContextual_Result")
freq_AA, freq_coupleAA = blosumCode.freqBlosum(count_AA, nb_AA, count_coupleAA, nb_coupleAA, path_folder_Result)


# blosum
blosum = blosumCode.blosumScore(freq_AA, freq_coupleAA, path_folder_Result, scale_factor = 2)
blosumCode.blosumVisualisation(blosum)

matrix_diff, pid_inf_ref, average_euclidean_d, average_diff = blosumCode.diffBlosum(blosum, pid_inf_ref = 62)
title_heatmap = f"Heatmap of Blosum({name_folder_fasta}) - Blosum{pid_inf_ref}Ref:\nmean difference = {average_diff}, mean euclidean distance = {average_euclidean_d}"
blosumCode.heatmapBlosum(blosum, path_folder_Result, title_heatmap, size_annot = 5)


# conditional proba
cond_proba = blosumCode.blosumConditionalProba(freq_AA, freq_coupleAA, path_folder_Result)
blosumCode.sumLine(cond_proba)
title_heatmap = f"Heatmap of the conditional probability matrix \n computed on {name_folder_fasta}"
blosumCode.heatmapBlosum(cond_proba, path_folder_Result, title_heatmap)