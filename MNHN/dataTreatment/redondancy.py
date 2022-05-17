import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

from MNHN.utils.fastaReader import readFastaMul
from MNHN.dataTreatment.pid import pid


def lenSeqCorrected(seq, list_residu):
    """
    Return the number of residus in seq that are included in list_residu
    """
    len_seq_corrected = 0
    for aa in seq:
        if aa in list_residu:
            len_seq_corrected += 1
    return len_seq_corrected



def clusterAntiRedundancy(liste_seq, file_seq_non_redondant, included_residue, pid_sup):    
    """
    Return a partition of liste_seq of sequences with a percentage of identity greater or equal than pid_sup
    """
    cluster = {}
    if liste_seq:   # if the list is not empty
        name_0, seq_0 = liste_seq[0] 
        len_seq_real_0 = lenSeqCorrected(seq_0, included_residue)
        cluster[0] = [(name_0, seq_0, len_seq_real_0)]

        for name_1, seq_1 in liste_seq:
            len_seq_real_1 = lenSeqCorrected(seq_1, included_residue)
            group = 0
            indice = 0

            while group <= len(cluster) - 1 and indice <= len(cluster[group]) - 1:
                seq_2 = cluster[group][indice][1]
                pourcentage_id = pid(seq_1, seq_2) 
                if pourcentage_id < pid_sup:
                    group += 1
                    indice = 0
                else:
                    if indice == len(cluster[group]) - 1:
                        cluster[group].append((name_1, seq_1, len_seq_real_1))
                        indice += 2 # avoid infinite loop
                    else:
                        indice += 1
            if group == len(cluster):
                cluster[group] = [(name_1, seq_1, len_seq_real_1)]
    else:
        print(file_seq_non_redondant)
    return cluster



def representativeNonRedundant(cluster):
    """
    Select the first sequence with the longest length in the cluster as the cluster representative
    """
    seq_non_redundant = []
    for group in cluster:
        current_group = cluster[group]
        representative = current_group[0]   
        for elem in current_group:
            if elem[2] > representative[2]: # 2 stands for the coorected lenght of a sequence
                representative = elem
        seq_non_redundant.append(representative[0])   # 0 stands for the name of the sequence
    return seq_non_redundant




def nonRedundant(path_file_fasta, path_file_seq_non_redundant, list_residu, pid_sup):
    """
    Rewrite the fasta file by correcting the issue of redundancy according to pid_sup.
    """
    seed = readFastaMul(path_file_fasta)
    cluster = clusterAntiRedundancy(seed, path_file_seq_non_redundant, list_residu, pid_sup)
    seq_non_redundant = representativeNonRedundant(cluster)

    with open(path_file_fasta, "r") as file:
        with open(path_file_seq_non_redundant, "w") as file_corrected:
            flag_write = False
            for line in file:
                if line[0] == ">":   
                    if line[1:-1].split(" ")[0] in seq_non_redundant:    # keep the name only 
                        flag_write = True
                    else:
                        flag_write = False
                if flag_write == True:
                    file_corrected.write(line)