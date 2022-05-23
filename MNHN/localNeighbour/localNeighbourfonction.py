import numpy as np
from pathlib import Path
import pandas as pd
import os, shutil


import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


from MNHN.utils.timer import Timer
import MNHN.utils.fastaReader as fastaReader
import MNHN.utils.folder as folder



def triplet_count_for_cube_initialisation(list_residu):
    """
    Initialisation of the count of triplets at 1 for each valid triplet
    to avoid issues when the triplet is not in data_train.

    Notice: a triplet is the succession of 
    the amino acid of the known sequence considered (k), the aligned amino acid to predict (p), the neighbor amino acid considered (c).
    """
    triplet_count = {}
    for aa_k in list_residu:
        triplet_count[aa_k] = {}
        for aa_p in list_residu:
            triplet_count[aa_k][aa_p] = {}
            for aa_c in list_residu:
                triplet_count[aa_k][aa_p][aa_c]  = 1   
    return triplet_count



    

def triplet_count_for_cube(file_fasta, path_folder_pid, triplet_count, delay_num, kp_SeqChoice, list_residu, pid_inf):    
    """
    Count the valid triplets of amino acid in a seed.

    Notice: a triplet is the succession of 
    the amino acid of the known sequence considered (k), the aligned amino acid to predict (p), the neighbor amino acid considered (c).

    The triplet is counted when its 3 amino acids belongs to 2 sequences with pid > pid_sup 
    and they are in list_residu.
    """
    # k: known
    # p: predict
    # c: context
    # kp_SeqChoice: choice the reference sequence to look at its neighbors
    accession_num = folder.get_accession_number(file_fasta)
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pId.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])

    if seed:
        for name_k, seq_k in seed:
            for name_p, seq_p in seed:
                if name_k != name_p:
                    if pid_couple[name_k][name_p] >= pid_inf:
                        if delay_num < 0:   # before
                            index_range = - delay_num, len_seq
                        elif delay_num > 0: # after
                            index_range = 0, len_seq - delay_num
                        
                        if kp_SeqChoice == "k":
                            seq_c = seq_k
                        elif kp_SeqChoice == "p":
                            seq_c = seq_p

                        for aa_index in range(index_range[0], index_range[1]):       
                            aa_k = seq_k[aa_index] 
                            aa_p = seq_p[aa_index] 
                            if all(x in list_residu for x in [aa_k, aa_p]):     
                                index_neighbor = aa_index + delay_num
                                aa_c = seq_c[index_neighbor] 
                                if aa_c in list_residu: 
                                    triplet_count[aa_k][aa_p][aa_c] += 1
    else:
        print(accession_num)

    return triplet_count


def triplet_conditional_proba(list_residu, triplet_count):
    # pseudo_count idea removed because we have enough data
    intra_couple_count = {}
    for aa_k in list_residu:
        intra_couple_count[aa_k] = {}
        for aa_c in list_residu: 
            intra_couple_count[aa_k][aa_c] = 0
            for aa_p in list_residu:
                intra_couple_count[aa_k][aa_c] += triplet_count[aa_k][aa_p][aa_c]

    cond_proba = {}
    for aa_k in list_residu:
        cond_proba[aa_k] = {}
        for aa_p in list_residu:
            cond_proba[aa_k][aa_p] = {}
            for aa_c in list_residu: 
                if intra_couple_count[aa_k][aa_c] != 0:
                    cond_proba[aa_k][aa_p][aa_c] = (triplet_count[aa_k][aa_p][aa_c]) / (intra_couple_count[aa_k][aa_c])  
                else:
                    cond_proba[aa_k][aa_p][aa_c] = 0

    return cond_proba                 



def sum_line(cond_proba, list_residu, aa_k, aa_c):
    """
    The sum of the conditional probabilities on each line (for aa_k and aa_c fixed) 
    must be equal to 1 to respect the total probability formula.
    """
    sum_line = 0
    for aa_p in list_residu:
        sum_line += cond_proba[aa_k][aa_p][aa_c]
    return sum_line



def sum_plate(cond_proba):
    """
    cond_proba: cube of the conditional probabilities of each valid triplet.

    The sum of the conditional probabilities on each horizontal level of the cube
    must be equal to the length of an edge of the cube to respect the total probability formula.
    """
    for aa_1 in cond_proba:
        sum_plateau = 0
        for aa_2 in cond_proba[aa_1]:
            for aa_3 in cond_proba[aa_1][aa_2]:
                sum_plateau += cond_proba[aa_1][aa_2][aa_3]
        print(f"{aa_1}, {sum_plateau}")


# path_NeighborRes  # Create target Directory if don't exist dans le main
def multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf = 62):
    t = Timer()
    t.start()

    triplet_count = triplet_count_for_cube_initialisation(list_residu)

    files = Path(path_folder_fasta).iterdir()
    name_folder_fasta = os.path.basename(path_folder_fasta)

    for file in files:
            triplet_count = triplet_count_for_cube(file, path_folder_pid, triplet_count, delay_num, kp_SeqChoice, list_residu, pid_inf)

    t.stop("Compute the valid triplets count")

    return triplet_count, name_folder_fasta
    

def cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu):
    t = Timer()
    t.start()
    print(f"Conditional probability cube: {delay_num},{kp_SeqChoice}")
    cond_proba = triplet_conditional_proba(list_residu, triplet_count)
    path_proba_cond = f"{path_NeighborRes}/proba_cond_({str(delay_num)},{kp_SeqChoice})_{name_folder_fasta}"
    np.save(path_proba_cond, cond_proba) 
    path_proba_cond = f"{path_proba_cond}.npy"
    t.stop("Compute a cube i.e the conditional probability matrix with 1 neighbour")

    return cond_proba, path_proba_cond


def dico_visualisation(dico): # to factorise
    """
    Visualisation of the matrix
    """
    df_dico = np.transpose(pd.DataFrame.from_dict(dico))  
    print(df_dico)

    return df_dico