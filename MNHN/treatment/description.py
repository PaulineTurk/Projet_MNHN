from pathlib import Path
import matplotlib.pyplot as plt
import os

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

import MNHN.utils.fastaReader as fastaReader
from MNHN.utils.timer import Timer





def data_count(path_data: str, list_residu):
    """
    path_folder: path of the folder of fasta files to describe
    """
    # initialisation
    nbre_seed = 0
    nbre_seq = 0
    total_position = 0
    total_residu = 0            # en fait il s'agit des caractères inclus
    total_character = 0
    character_count = {}        # consider all residus (dico constructed along the way)

    residu_count = {}  
    for aa in list_residu:
        residu_count[aa] = 0 

    # counting
    files = Path(path_data).iterdir()
    for file in files:
        nbre_seed += 1
        data_Pfam = fastaReader.read_multi_fasta(file)
        len_seq = len(data_Pfam[0][1])
        total_position += len_seq
        for _, seq in data_Pfam:
            nbre_seq += 1
            total_character += len_seq
            for aa in seq:
                if aa in character_count:
                    character_count[aa] += 1
                else:
                    character_count[aa] = 1
                
                if aa in list_residu:
                    total_residu += 1 
                    residu_count[aa] += 1
                    

    print("nbre_seed:               ", '{:_.2f}'.format(nbre_seed))
    print("nbre_seq:                ", '{:_.2f}'.format(nbre_seq))
    print("nbre_position:           ", '{:_.2f}'.format(total_position))
    print("total_character:         ", '{:_.2f}'.format(total_character))
    print("total_character_included:", '{:_.2f}'.format(total_residu))

    # mean len sequence
    if nbre_seq != 0:
        mean_len_seq = round(total_residu/nbre_seq, 2)
        print("mean_len_seq:            ", '{:_.2f}'.format(mean_len_seq))
    else: 
        print("no sequence")

    # mean nbre sequence per seed
    if nbre_seed != 0:
        mean_nbre_seq = round(nbre_seq/nbre_seed, 2)
        print("mean_nbre_seq:           ", '{:_.2f}'.format(mean_nbre_seq))
    else:
        print('no seed')

    return residu_count, total_residu, character_count, total_character


def bar_plot_data_count(path_folder_to_describe: str,  residu_count: dict, total_residu: float, entity_name):
    """
    path_folder_to_describe: to name the graph and the figure according to the folder described
    descriptor: dictionary that contains the feature information for each residu
    feature: can be "count" or "percentage"
    """
    residu_percentage = {k: round(100*v / total_residu, 2) for k, v in residu_count.items()}

    # dico sorted by descending order
    residu_percentage_sorted = dict(sorted(residu_percentage.items(), key=lambda item: item[1], reverse=True))

    plt.bar(list(residu_percentage_sorted.keys()), residu_percentage_sorted.values(), color='g')
    plt.xlabel(entity_name)
    plt.ylabel('Percentage')

    dir_image = os.path.dirname(path_folder_to_describe)
    name_dir = os.path.basename(path_folder_to_describe)

    title_graph = f"{entity_name} percentage in {name_dir}"
    title_graph_object = f"{dir_image}/{title_graph}"

    plt.title(title_graph)
    plt.savefig(title_graph_object)
    plt.close()