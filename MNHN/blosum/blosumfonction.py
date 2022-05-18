import pandas as pd
from math import log2
import numpy as np
import os 
from math import sqrt
import blosum as bl
import seaborn as sb
import matplotlib.pyplot as plt

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  
sys.path.append(str(package_root_directory_MNHN))  
from MNHN.utils.timer import Timer
import MNHN.utils.fastaReader as fastaReader
import MNHN.utils.folder as folder



def count_for_blosum(num_accession, path_folder_pid, seed, pid_inf,  
                count_AA, nb_AA, count_coupleAA, nb_coupleAA, list_residu):   
    """
    In the seed with id num_accession, count the number of each valid amino acid 
    and the number of valid couple
    
    num_accession: pid of the seed
    name_folder_pid: the pid file of this seed
    seed: the (name, seq) tuples of this seed
    pid_inf: smaller pid to validate a couple of sequence
    count_AA: dictionary of the count of each valid amino acid to cumulate this count on many seeds
    nb_AA: count of all valid amino acids
    count_coupleAA: dictionary of the count of each valid couple of amino acids to cumulate this count on many seeds
    nb_coupleAA: count of all valid couple of amino acids
    list_residu: list of valid amino acids
    """

    pid_couple = np.load(f"{path_folder_pid}/{num_accession}.pId.npy", allow_pickle='TRUE').item()
    nb_seq = len(seed)
    for i in range(nb_seq):
        name_1, seq_1 = seed[i]
        for j in range(i + 1, nb_seq):
            name_2 ,seq_2 = seed[j] 
            if pid_couple[name_1][name_2] >= pid_inf:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in list_residu and aa_2 in list_residu:

                        # count AA
                        count_AA[aa_1] += 1
                        count_AA[aa_2] += 1
                        nb_AA += 2

                        # count couple AA
                        if aa_1 == aa_2:
                            count_coupleAA[aa_1][aa_2] += 2
                        else:
                            count_coupleAA[aa_1][aa_2] += 1
                            count_coupleAA[aa_2][aa_1] += 1
                        nb_coupleAA += 2

    return count_AA, nb_AA, count_coupleAA, nb_coupleAA


def multi_count_for_blosum(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62):
    """
    Iterate count_for_blosum on the files included in path_folder_fasta.
    """
    t = Timer()
    t.start()

    # intialisation of the count of each valid residu
    count_AA = {}  
    for aa in list_residu:
        count_AA[aa] = 0
    nb_AA = 0

    # intialisation of the count of each valid couple of residus
    count_coupleAA = {}
    for aa_1 in list_residu:
        count_coupleAA[aa_1] = {}
        for aa_2 in list_residu:
                count_coupleAA[aa_1][aa_2] = 0
    nb_coupleAA = 0

    # files to use in order to evaluate Blosum
    files = Path(path_folder_fasta).iterdir()

    for file in files:
            accession_num = folder.get_accession_number(file)
            seed_train = fastaReader.read_multi_fasta(file)
            count_AA, nb_AA, count_coupleAA, nb_coupleAA = count_for_blosum(accession_num, path_folder_pid, seed_train, pid_inf,
                                                                            count_AA, nb_AA, count_coupleAA, nb_coupleAA, list_residu)

    t.stop("Compute the count of amino acid and couple of amino acids")

    return count_AA, nb_AA, count_coupleAA, nb_coupleAA




def freq_for_blosum(count_AA, nb_AA, count_coupleAA, nb_coupleAA, path_folder_Result):
    """
    Compute and save the frequences of each valid amino acid 
    and each valid couple of amino acids.
    """
    # get the list of valid residus
    list_residu = count_AA.keys()

    t = Timer()
    t.start()

    # frequence of each valid amino acid
    freq_AA = {}
    for aa in list_residu:
        if nb_AA != 0:
            freq_AA[aa] = count_AA[aa]/nb_AA
        else: 
            freq_AA[aa] = 0
    
    # frequence of each couple of amino acid
    freq_coupleAA = {}
    for aa_1 in list_residu:
        freq_coupleAA[aa_1] = {}
        for aa_2 in list_residu:
            if nb_coupleAA != 0:
                freq_coupleAA[aa_1][aa_2] = count_coupleAA[aa_1][aa_2]/nb_coupleAA
            else:
                freq_coupleAA[aa_1][aa_2] = 0
    t.stop("Compute the frequence of amino acid and couple of amino acids")

    path_freqAA = f"{path_folder_Result}/Blosum_freq_AA"
    np.save(path_freqAA, freq_AA) 

    return freq_AA, freq_coupleAA





def blosum_score(freq_AA, freq_coupleAA, path_folder_Result, scale_factor = 2):
    """
    Compute and save the Blosum matrix 
    """
    list_residu = freq_AA.keys()

    t = Timer()
    t.start()
    blosum = {}
    for aa_1 in list_residu:
        blosum[aa_1] = {}
        for aa_2 in list_residu:
            if freq_coupleAA[aa_1][aa_2] != 0:
                blosum[aa_1][aa_2] = round(scale_factor * log2(freq_coupleAA[aa_1][aa_2]    
                                                            / (freq_AA[aa_1] * freq_AA[aa_2])))
            else:
                blosum[aa_1][aa_2] = 0 
    t.stop("Compute the blosum matrix")

    path_matrix = f"{path_folder_Result}/Blosum_score"
    np.save(path_matrix, blosum) 
    return blosum



def blosum_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result):
    """
    Compute and save the matrix of conditional probabilities
    """
    list_residu = freq_AA.keys()

    t = Timer()
    t.start()
    cond_proba = {}
    for aa_1 in list_residu:
        cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            if freq_coupleAA[aa_1][aa_2] != 0:
                cond_proba[aa_1][aa_2] = freq_coupleAA[aa_1][aa_2]/freq_AA[aa_1]
            else:
                cond_proba[aa_1][aa_2] = 0
    t.stop("Compute the blosum conditional probability matrix")

    path_cond_proba = f"{path_folder_Result}/Blosum_proba_cond"
    np.save(path_cond_proba, cond_proba)
    return cond_proba



def blosum_heatmap(matrix, path_folder_Result, title, size_annot = 3):
    """
    Save the heatmap of the matrix in path_folder_Result
    """
    heatmap_matrix = pd.DataFrame(matrix).T.fillna(0)  # à transposer?
    heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": size_annot}, fmt = '.2g')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(title)
    plt.close()
    path_save_fig = f"{path_folder_Result}/{title}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)


def blosum_visualisation(blosum):
    """
    Visualisation of the matrix
    """
    df_blosum = np.transpose(pd.DataFrame.from_dict(blosum))  
    print(df_blosum)
    return df_blosum



def sum_line(blosum):
    """
    To check that the somme of the line is equal to one 
    for the conditional probability matrix
    """
    df_blosum = np.transpose(pd.DataFrame.from_dict(blosum)) 
    sum_line = df_blosum.sum(axis=1)
    print("Sum of the line:\n", sum_line)




def blosum_difference(blosum, pid_inf_ref = 62):
    """
    Quantify the distance between the blosum computed and a blosum of reference
    """

    list_residu = blosum.keys()

    # blosum ref importation
    blosum_ref = bl.BLOSUM(pid_inf_ref) 

    # initialisation
    matrix_diff = {}
    average_diff = 0
    count = 0
    euclidean_d = 0

    # evaluation of the differences
    for aa1 in list_residu:
        matrix_diff[aa1] = {}
        for aa2 in list_residu:
            matrix_diff[aa1][aa2] = int(blosum[aa1][aa2] - blosum_ref[aa1 + aa2])
            average_diff += matrix_diff[aa1][aa2]
            euclidean_d += (matrix_diff[aa1][aa2])**2
            count += 1

    average_euclidean_d = round(sqrt(euclidean_d/count), 2)
    average_diff = round(average_diff/count, 2)

    return matrix_diff, pid_inf_ref, average_euclidean_d, average_diff 