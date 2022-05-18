import mainNeighbor
import sys



path_folder_pId = sys.argv[1]
path_folder_data_split = sys.argv[2]
path_new_folder = sys.argv[3]




list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
name_NeighborRes = "NeighborRes"
name_data_train = "Pfam_A"
path_folder_fasta_train = f"{path_folder_data_split}/{name_data_train}"
path_NeighborRes = f"{path_new_folder}/{name_NeighborRes}"
percentage_train = 50


path_proba_cond =  mainNeighbor.simpleContextualBlosum(path_folder_fasta_train, percentage_train, path_folder_pId, path_NeighborRes, delay_num, kp_SeqChoice, list_residu, pid_inf = 62)
print(path_folder_pId)