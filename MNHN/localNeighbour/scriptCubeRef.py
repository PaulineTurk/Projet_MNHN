import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [3]
sys.path.append(str(package_root_directory_MNHN))
import MNHN.localNeighbour.localNeighbourfonction as localNeighbourfonction

path_folder_pid = sys.argv[1]    # chemin des pid
path_folder_fasta = sys.argv[2]  # chemin des seeds d'entrainement
path_NeighborRes = sys.argv[3]   # folder a créer avant de lancer le calcul des cubes
delay_num = -1
kp_SeqChoice = "k"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

triplet_count, name_folder_fasta = localNeighbourfonction.multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf = 62)
cond_proba, path_proba_cond = localNeighbourfonction.cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)