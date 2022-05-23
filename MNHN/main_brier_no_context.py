import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.brier.brierScore as brierScore
import MNHN.brier.predictor as predictor

path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_train"
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


brier_score = brierScore.brier_score_predictor_01(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62)
print("\nworst predictor")
print(brier_score)

brier_score = brierScore.brier_score_predictor_perfect(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62)
print("\nbest predictor")
print(brier_score)


path_cond_proba = "/Users/pauline/Desktop/data_Result/Pfam_split/NonContextual_Result/Blosum_proba_cond.npy"
unit_Brier = predictor.predictor_blosum(path_cond_proba, list_residu)
print("\npredictor Blosum")
brier_score = brierScore.brier_score_matrix(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62)
print(brier_score)


# Version of Brier Score per amino-acid:
print("\npredictor Blosum per amino acid")
brier_score,_,_ = brierScore.brier_score_matrix_v2(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62)
#print(brier_score)