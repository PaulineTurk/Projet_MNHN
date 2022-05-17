from pickle import FALSE
from pathlib import Path


import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

from MNHN.utils.timer import Timer
import MNHN.utils.folder as folder


def capitalization(path_file, path_file_corrected):
    """
    Convert all the lowercase residu into uppercase.

    path_file: path of the fasta file to correct
    path_file_corrected; path of the fasta file corrected
    """
    with open(path_file, "r") as file:
        with open(path_file_corrected, "w") as file_corrected:
            for line in file:
                if line[0] != ">":   
                    line = line.upper()
                file_corrected.write(line)



def multi_capitalization(path_data, path_data_corrected):
    """
    Convert all the lowercase residu into uppercase.

    path_data: path of the folder of fasta file to correct
    path_data_corrected: folder in which the fasta file corrected are saved
    """
    t = Timer()
    t.start()

    folder.creat_folder(path_data_corrected)

    files = Path(path_data).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        path_file_corrected = f"{path_data_corrected}/{accession_num}.fasta.upper"
        capitalization(file, path_file_corrected)
    t.stop("Correction upper files")

