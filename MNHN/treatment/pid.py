import numpy as np


import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

import MNHN.utils.fastaReader as fastaReader
from MNHN.utils.timer import Timer
import MNHN.utils.folder as folder


# calcul correct du pid (ne compter que les matchs d'acides aminés standards 
# et diviser par la longeur min d'acide aminé standard entre les 2 séquences)
def pid(seq_1, seq_2, list_inclusion):  
    """
    Return the percentage of identity between the two sequences: seq_1 and seq_2
    list_inclusion: liste des caractères inclus
    """
    pid = 0
    nb_included_character_seq_1 = 0
    nb_included_character_seq_2 = 0

    len_seq = len(seq_1)
    for indice_aa in range(len_seq):
        if seq_1[indice_aa] in list_inclusion:
            nb_included_character_seq_1 += 1
        if seq_2[indice_aa] in list_inclusion:
            nb_included_character_seq_2 += 1
            
        if seq_1[indice_aa] in list_inclusion and seq_2[indice_aa] in list_inclusion and seq_1[indice_aa] == seq_2[indice_aa]:
                pid += 1

    return 100*pid/min(nb_included_character_seq_1, nb_included_character_seq_2)


def pid_two_seq(path_fasta_file, path_file_pId, list_inclusion):
    liste_seq = fastaReader.read_multi_fasta(path_fasta_file)
    pid_couple = {}
    for name_1, seq_1 in liste_seq:
        pid_couple[name_1] = {}
        for name_2, seq_2 in liste_seq:
            pid_couple[name_1][name_2] = pid(seq_1, seq_2, list_inclusion)
    np.save(path_file_pId, pid_couple) 


def save_pid(path_folder_fasta, path_folder_pid, list_inclusion):
    """
    For each fasta file in path_folder_fasta, compute the dictionary of pid 
    for each couple of sequences and save it in path_folder_pid
    """
    t = Timer()
    t.start()

    folder.creat_folder(path_folder_pid)

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        path_file_pid = f"{path_folder_pid}/{accession_num}.pid"
        pid_two_seq(file, path_file_pid, list_inclusion)
        
    t.stop("Compute and save the pid files")