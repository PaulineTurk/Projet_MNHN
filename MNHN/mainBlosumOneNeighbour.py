import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.blosumNeighbour.matrix as blosumNCode
from MNHN.utils.folder import creatFolder

 
path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_10"  # data train
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
delay_num = -1
kp_SeqChoice = "k"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
triplet_count, name_folder_fasta = blosumNCode.multiTripletCount(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf = 62)


path_NeighborRes = creatFolder("/Users/pauline/Desktop/data_Result/Pfam_split/OneNeighbour_Result")
cond_proba, path_proba_cond = blosumNCode.cubeOneNeighbour(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)
blosumNCode.blosumVisualisation(cond_proba)