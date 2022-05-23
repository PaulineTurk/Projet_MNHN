import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.localNeighbour.localNeighbourfonction as localNeighbourfonction
import MNHN.utils.folder as folder

 
path_folder_fasta = "/Users/pauline/Desktop/data_Result/Pfam_split/Pfam_train"  # data train
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
delay_num = -1
kp_SeqChoice = "k"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
triplet_count, name_folder_fasta = localNeighbourfonction.multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf = 62)


path_NeighborRes = folder.creat_folder("/Users/pauline/Desktop/data_Result/Pfam_split/OneNeighbour_Result")
cond_proba, path_proba_cond = localNeighbourfonction.cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)
localNeighbourfonction.dico_visualisation(cond_proba)