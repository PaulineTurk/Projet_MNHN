import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[1]
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.folder as folder
import MNHN.treatment.description as description
import MNHN.treatment.stockholm as stockholm
import MNHN.treatment.capitalizer as capitalizer
import MNHN.treatment.split as split
import MNHN.treatment.pid as pid
import MNHN.treatment.redundancy as redundancy


list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

path_data = "/Users/pauline/Desktop/pfam_description_soutenance_mi_stage/Pfam_fasta_99"

print(f"\nVisualisation of {path_data}")
residu_count, total_residu, character_count, total_character = description.data_count(path_data, list_residu)
description.bar_plot_data_count(path_data, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_data, character_count, total_character , "Character")
