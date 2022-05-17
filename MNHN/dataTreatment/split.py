import shutil

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

from MNHN.utils.folder import getAccessionNb



def addFileFromFolder(file_name, folder_name_source, folder_path_target, extension_file_name_target):
    path_file_source = f"{folder_name_source}/{file_name}"
    accession_num = getAccessionNb(path_file_source)

    path_file_target = f"{folder_path_target}/{accession_num}.{extension_file_name_target}"
    shutil.copy2(path_file_source, path_file_target)