import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[1]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.folder as folder
import MNHN.treatment.description as description
import MNHN.treatment.stockholm as stockholm
import MNHN.treatment.capitalizer as capitalizer
import MNHN.treatment.split as split
import MNHN.treatment.pid as pid
import MNHN.treatment.redundancy as redundancy


path_file = "/Users/pauline/Desktop/data_Test/Pfam_Sample"
path_folder_save = "/Users/pauline/Desktop/data_Test/Pfam_Stockholm"

stockholm.stockholm_separator(path_file, path_folder_save)
path_folder_stockholm = path_folder_save

path_folder_fasta = "/Users/pauline/Desktop/data_Test/Pfam_Fasta"
stockholm.multi_stockholm_to_fasta(path_folder_stockholm, path_folder_fasta)
path_data = path_folder_fasta
print(f"\nVisualisation of {path_data}")
residu_count, total_residu = description.data_count(path_data)
description.bar_plot_data_count(path_data, residu_count, total_residu)

path_data_corrected = "/Users/pauline/Desktop/data_Test/Pfam_Upper"
capitalizer.multi_capitalization(path_data, path_data_corrected)
path_data = path_data_corrected
print(f"\nVisualisation of {path_data}")
residu_count, total_residu = description.data_count(path_data)
description.bar_plot_data_count(path_data, residu_count, total_residu)


path_folder_fasta = path_data_corrected # upper correction only
path_folder_pid = "/Users/pauline/Desktop/data_Test/PID"
pid.save_pid(path_folder_fasta, path_folder_pid)



path_folder_fasta = path_data_corrected # still upper correction only
path_folder_fasta_nonRedondant = "/Users/pauline/Desktop/data_Test/Pfam_nonRedondant"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
redundancy.multi_non_redundancy_correction(path_folder_fasta, path_folder_fasta_nonRedondant, 
                          list_residu, pid_sup = 99)

path_folder_data = path_folder_fasta_nonRedondant  # redondant corrected
result_folder = "/Users/pauline/Desktop/data_Result"
folder.creat_folder("/Users/pauline/Desktop/data_Result")
path_folder_data_split = f"{result_folder}/Pfam_split" 
percentage_A = 90
name_data_A, name_data_B = "Pfam_90", "Pfam_10"
split.data_split(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B)