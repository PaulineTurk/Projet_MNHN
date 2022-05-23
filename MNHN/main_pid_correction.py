import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[1]
sys.path.append(str(package_root_directory_MNHN))



import MNHN.treatment.pid as pid
import MNHN.utils.dicoVisualizer as dicoVisualizer


list_standard_aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
path_data_corrected = "/Users/pauline/Desktop/Structuration/Pfam_verif_cube"
list_inclusion = list_standard_aa
path_folder_fasta = path_data_corrected # upper correction only
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID_v2"
pid.save_pid_v2(path_folder_fasta, path_folder_pid, list_inclusion)

dicoVisualizer.dico_visualizer("/Users/pauline/Desktop/data_Test/PID_v2/PF00002.27.pid.npy")