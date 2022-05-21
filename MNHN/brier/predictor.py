from random import choice
import numpy as np
import copy


def predictor_01(seed, num_accession, path_folder_pid, brier_score, count, list_residu, pid_inf):
    """
    Worst predictor that always predicts with a probability of 1 a random uncorrect amino acid.

    brier_score: somme of the score for each prediction
    count: number of predictions to normalise the brier_score
    """
    pid_couple = np.load(f"{path_folder_pid}/{num_accession}.pId.npy", allow_pickle='TRUE').item()
  
    for name_1, seq_1 in seed:
        for name_2 ,seq_2 in seed:
            if name_1 != name_2:     
                if pid_couple[name_1][name_2] >= pid_inf:
                    for (aa_1, aa_2) in zip(seq_1, seq_2):
                        if aa_1 in list_residu and aa_2 in list_residu:
                            count += 1

                            # random selection of the wrong prediction
                            new_list_residu = list_residu.copy()
                            new_list_residu.remove(aa_2)
                            aa_predicted = choice(new_list_residu)

                            # brier score computing for one prediction
                            for j in list_residu:
                                if aa_2 == j:
                                    # probability of predicting the correct amino acid = 0
                                    proba = 0
                                else:
                                    if j == aa_predicted:
                                        # probability of predicting the wrong amino acid predicted = 1
                                        proba = 1
                                    else:
                                        # probability of predicting any other amino acid  = 0
                                        proba = 0
                                brier_score += (proba - int(aa_2 == j))**2

    return brier_score, count


def predictor_perfect(seed, num_accession, path_folder_pid, brier_score, count, list_residu, pid_inf):
    """
    Perfect predictor that always predict with a probability of 1 the correct amino acid.
    """
    pid_couple = np.load(f"{path_folder_pid}/{num_accession}.pId.npy", allow_pickle='TRUE').item()
    for name_1, seq_1 in seed:
        for name_2 ,seq_2 in seed:
            if name_1 != name_2:     
                if pid_couple[name_1][name_2] >= pid_inf:
                    for (aa_1, aa_2) in zip(seq_1, seq_2):
                        if aa_1 in list_residu and aa_2 in list_residu:
                            count += 1
                            for j in list_residu:
                                if aa_2 == j:
                                    # probability of predicting the correct amino acid = 1
                                    proba = 1
                                else:
                                    # probability of predicting the wrong amino acid = 0
                                    proba = 0
                                brier_score += (proba - int(aa_2 == j))**2

    return brier_score, count


def predictor_blosum(path_cond_proba, list_residu):
    """
    Predictor that uses the conditional probability blosum matrix.
    """
    cond_proba = np.load(path_cond_proba ,allow_pickle='TRUE').item()
    unit_Brier = brier_unit(cond_proba, list_residu)

    return unit_Brier


def predicteur_equiprobable(list_residu):
    """
    Predictor that predict each amino acid with the same probability.
    """
    nb_aa = len(list_residu)
    equiproba = 1/nb_aa

    cond_proba = {}
    for aa_1 in list_residu:
        cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            cond_proba[aa_1][aa_2] = equiproba

    unit_Brier = brier_unit(cond_proba, list_residu)

    return unit_Brier


def predictor_stationary(path_freq_aa, list_residu):
    """
    Predictor that predict an amino acid with the frequency of that amino acid. 
    """
    freq_aa = np.load(path_freq_aa, allow_pickle='TRUE').item()

    cond_proba = {}
    for aa_1 in list_residu:
        cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            cond_proba[aa_1][aa_2] = freq_aa[aa_2]

    unit_Brier = brier_unit(cond_proba)

    return unit_Brier


def predictor_identity(list_residu):
    """
    Predictor that predict with a probabilty of 1 that the amino acid will not be substituted.
    """
    cond_proba = {}
    for aa_1 in list_residu:
        cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            if aa_1 == aa_2:
                cond_proba[aa_1][aa_2] = 1
            else:
                cond_proba[aa_1][aa_2] = 0
    
    unit_brier = brier_unit(cond_proba, list_residu)

    return unit_brier


def brier_unit(cond_proba, list_residu):  
    unit_Brier = {}
    for aa_1 in list_residu:
        unit_Brier[aa_1] = {}
        for aa_2 in list_residu:
            unit = 0
            for j in list_residu: 
                unit += (cond_proba[aa_1][j] - int(aa_2 == j))**2 
            unit_Brier[aa_1][aa_2] = unit

    return unit_Brier


def brier_matrix(unit_Brier, seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf):
    """
    Compute brier_score and count for matrix predictors
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
    nb_seq = len(seed)

    for i in range(nb_seq):
        name_1, seq_1 = seed[i]
        for j in range(i + 1, nb_seq):
            name_2 ,seq_2 = seed[j] 
            if pid_couple[name_1][name_2] >= pid_inf:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in list_residu and aa_2 in list_residu:
                        count += 2
                        brier_score += unit_Brier[aa_1][aa_2]
                        brier_score += unit_Brier[aa_2][aa_1]

    return brier_score, count


def brier_matrix_v2(unit_Brier, seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf):
    """
    Compute brier_score and count for matrix predictors
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
    nb_seq = len(seed)

    for i in range(nb_seq):
        name_1, seq_1 = seed[i]
        for j in range(i + 1, nb_seq):
            name_2 ,seq_2 = seed[j] 
            if pid_couple[name_1][name_2] >= pid_inf:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in list_residu and aa_2 in list_residu:
                        count[aa_1] += 1
                        count[aa_2] += 1
                        brier_score[aa_1] += unit_Brier[aa_1][aa_2]
                        brier_score[aa_2] += unit_Brier[aa_2][aa_1]

    return brier_score, count