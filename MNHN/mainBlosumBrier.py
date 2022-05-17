import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.brier.brierScore as brier
from MNHN.brier.predictor import predictorBlosum

path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_10"
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


brier_score = brier.multiBrier01(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62)
print("\nworst predictor")
print(brier_score)

brier_score = brier.multiBrierPerfect(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62)
print("\nbest predictor")
print(brier_score)


path_cond_proba = "/Users/pauline/Desktop/data_Result/Pfam_split/NonContextual_Result/Blosum_proba_cond.npy"
unit_Brier = predictorBlosum(path_cond_proba, list_residu)
print("\npredictor Blosum")
brier_score = brier.multiBrierMatrix(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62)
print(brier_score)