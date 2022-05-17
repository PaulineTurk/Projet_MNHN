from pathlib import Path
import matplotlib.pyplot as plt
import os

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

from MNHN.utils.fastaReader import readFastaMul





def dataCount(path_data):
    """
    path_folder: path of the folder containing fasta files to describe
    """

    # initialisation
    nbre_seed = 0
    nbre_seq = 0
    total_position = 0
    total_residu = 0
    residu_count = {}   # consider all residus (dico construction along the way)

    # count
    files = Path(path_data).iterdir()
    for file in files:
        nbre_seed += 1

        data_Pfam = readFastaMul(file)
        len_seq = len(data_Pfam[0][1])
        total_position += len_seq
   
        for _, seq in data_Pfam:
            nbre_seq += 1
            total_residu += len_seq 
            for aa in seq:
                if aa in residu_count:
                    residu_count[aa] += 1
                else:
                    residu_count[aa] = 1

    print("nbre_seed:", '{:_.2f}'.format(nbre_seed))
    print("nbre_seq:", '{:_.2f}'.format(nbre_seq))
    print("nbre_position:", '{:_.2f}'.format(total_position))
    print("total_residu:", '{:_.2f}'.format(total_residu))

    # mean len sequence
    if nbre_seq != 0:
        mean_len_seq = round(total_residu/nbre_seq, 2)
        print("mean_len_seq:", '{:_.2f}'.format(mean_len_seq))
    else: 
        print("no sequence")

    # mean nbre sequence per seed
    if nbre_seed != 0:
        mean_nbre_seq = round(nbre_seq/nbre_seed, 2)
        print("mean_nbre_seq:", '{:_.2f}'.format(mean_nbre_seq))
    else:
        print('no seed')

    return residu_count, total_residu







def barPlot(path_folder_to_describe: str, descriptor: dict, feature: str):
    """
    path_folder_to_describe: to name the graph and the figure according to the folder described
    descriptor: dictionary that contains the feature information for each residu
    feature: can be "count" or "percentage"
    """
    plt.bar(list(descriptor.keys()), descriptor.values(), color='g')
    plt.xlabel('Residus')
    plt.ylabel('Percentage')
    dir_image = os.path.dirname(path_folder_to_describe)
    name_dir = os.path.basename(path_folder_to_describe)
    title_graph = f"Residu {feature} in {name_dir}"
    title_graph_object = f"{dir_image}/{title_graph}"
    plt.title(title_graph)
    plt.savefig(title_graph_object)
    plt.close()