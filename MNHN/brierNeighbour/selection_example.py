from pathlib import Path
import numpy as np
import os
import pandas as pd
import random
import math
from sklearn.utils import compute_sample_weight


import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]
sys.path.append(str(package_root_directory_MNHN))

import MNHN.utils.fastaReader as fastaReader 
import MNHN.utils.folder as folder
from MNHN.utils.timer import Timer


def one_seed_selection(accession_num, seed, pid, pid_inf, nb_non_info_seed, 
                     dico_seed, nb_valid_aa_couple_global, list_residu):
    """
    accession_num: pour retrouver le seed dans le dico
    """

    len_align = len(seed[0][1])   # à calculer qu'une seule fois car une longeur d'alignement dans un seed
    nb_seq = len(seed)
    print(f"seed accession num, nb seq, len alignment: {accession_num}, {len_align}, {nb_seq}")

    dico_seq = {}
    nb_valid_aa_couple_seed = 0

    for i in range(nb_seq):
        name_A, seq_A = seed[i]   # à tirer au sort après 1/2 que seq_1 = seq_A
        for j in range(i+1, nb_seq):
            name_B, seq_B = seed[j]

            if pid[name_A][name_B] >= pid_inf:
                # si le couple vérifie la condition sur pid_inf
                # initialisation de la liste des position valides hors context
                list_valid_position = []
                for index, aa_A in enumerate(seq_A):
                    if aa_A in list_residu:
                        if seq_B[index] in list_residu:
                            list_valid_position.append(index)
                
                nb_valid_position = len(list_valid_position)
                        

                if nb_valid_position != 0:    # ajouter le couple name_A, name_B au dico_seq en cours de remplissage 
                                              # ssi il contient au moins une séquence non contextuelle valide

                    dico_seq[(name_A, name_B)] = {}
                    dico_seq[(name_A, name_B)]["nb valid aa couple (seq)"] = nb_valid_position
                    dico_seq[(name_A, name_B)]["sequences A, B"] = (seq_A, seq_B)
                    dico_seq[(name_A, name_B)]["valid position"] = set(list_valid_position)   # ensemble des positions valides
                                          
                nb_valid_aa_couple_seed += nb_valid_position  # compté sur les paires de séquences possibles et non les couples


    if nb_valid_aa_couple_seed == 0: # seed inutil
        nb_non_info_seed += 1       
    else:
        dico_seed[accession_num] = {}
        dico_seed[accession_num]["nb valid aa couple (seed)"] = nb_valid_aa_couple_seed
        dico_seed[accession_num]["len alignment"] = len_align 

    nb_valid_aa_couple_global += nb_valid_aa_couple_seed

    return dico_seed, nb_non_info_seed, nb_valid_aa_couple_global, dico_seq


def multi_seeds_selection(path_folder_seed, path_folder_pid, pid_inf, list_residu, path_folder_dico_seq, path_file_dico_seed):
    """
    path_folder_dico_seq: dossier créé s'il n'existe pas, effacé puis recréé sinon.
                          Les dictionnaires de seed y sont enregistrés.
    path_file_dico_seed: dossier créé s'il n'existe pas, effacé puis recréé sinon.
                         Le dictionnaire sur les seeds y est enregistré 
                         ainsi que sa version normalisée par la somme de tous les comptes.
    """
    t = Timer()
    t.start()

    # initialisation des "compteurs" cumulés sur tous les seeds
    dico_seed = {}   # len(dico_seed) = nb_seed_valid
    nb_non_info_seed = 0  # pour vérifier ce compte + len(dico_seed) = len(path_folder_seed)
    nb_valid_aa_couple_global = 0

    # creation of path_folder_dico_seq
    path_folder_dico_seq = folder.creat_folder(path_folder_dico_seq)

    # creation of path_folder_dico_seed
    path_file_dico_seed = folder.creat_folder(path_file_dico_seed)

    files = Path(path_folder_seed).iterdir()
    for file in files:
        accession_num =  folder.get_accession_number(file)
        pid = np.load(f"{path_folder_pid}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
        seed = fastaReader.read_multi_fasta(file)

        dico_seed, nb_non_info_seed, nb_valid_aa_couple_global, dico_seq = one_seed_selection(accession_num, seed, pid, pid_inf, nb_non_info_seed, 
                                                                                            dico_seed, nb_valid_aa_couple_global, list_residu)
        
        #print(f"\n{accession_num}")
        #print(np.transpose(pd.DataFrame.from_dict(dico_seq)))

        # sauvegarde des dico_seq
        path_file_dico_seq = f"{path_folder_dico_seq}/{accession_num}.seq"
        np.save(path_file_dico_seq, dico_seq) 

    # sauvegarde du dico_seed
    path_save_file_dico_seed = f"{path_file_dico_seed}/seed"
    np.save(path_save_file_dico_seed, dico_seed) 

    # normalisation de dico_seed
    dico_seed_normalised = {}
    for key in dico_seed:
        dico_seed_normalised[key] = {}
        dico_seed_normalised[key]["proportion valid aa couple (seed)"] = dico_seed[key]["nb valid aa couple (seed)"]/nb_valid_aa_couple_global
        dico_seed_normalised[key]["len alignment"] = dico_seed[key]["len alignment"] 
    
    #print("dico_seed_normalised")
    #print(np.transpose(pd.DataFrame.from_dict(dico_seed_normalised)))
        
    # sauvegarde du dico_seed_normalised
    path_save_file_dico_seed_normalised = f"{path_file_dico_seed}/seed_normalised"
    np.save(path_save_file_dico_seed_normalised, dico_seed_normalised) 

    print(f"\nnb de seed valides, invalides: {len(dico_seed)}, {nb_non_info_seed}")

    t.stop("Traitement de la fraction Pfam entrainement")
    

def example_number_per_seed(path_dico_seed_normalised, nb_exemple_test, path_dico_exemple):
    """
    déterminer pseudo-aléatoirement le nombre d'exemple valide à piocher de chaque seed. 

    path_dico_seed: pré-calculé (à loader une fois)
    nb_tes_exemple: ordre de grandeur du nombre de couples valides à piocher au total de Pfam test.

    path_dico_exemple: chemin au fichier ou le file des exemples à prendre de chaque seed est enregistré
    """
    dico_seed_normalised = np.load(path_dico_seed_normalised, allow_pickle='TRUE').item()

    # conversion des proportions associées à chaque seed en nombre d'exemples à prendre de chaque seed

    ## multiplication de chaque proportion par le nombre d'exemples souhaité
    dico_exemple = {}
    for key in dico_seed_normalised:
        dico_exemple[key] = {}
        nb_ex_estimation = dico_seed_normalised[key]["proportion valid aa couple (seed)"] * nb_exemple_test

        # partie entière inf + int(random.uniform(0, 1) < partie décimale)      
        decimal_part, integer_part = math.modf(nb_ex_estimation)
        proba = random.uniform(0, 1)
        nb_ex_exact = integer_part + int(proba < decimal_part)   # généralisation en prévision de l'envoie d'exemple par batch 
                                                                     # au réseau de neurones (respecter l'ordre de grandeur du batch demandé ...) 
        dico_exemple[key]["nb exemples tests"] = nb_ex_exact
        dico_exemple[key]["len alignment"] = dico_seed_normalised[key]["len alignment"]

    path_dico_exemple = f"{path_dico_exemple}/exemple_seed"
    np.save(path_dico_exemple, dico_exemple) 

    print(np.transpose(pd.DataFrame.from_dict(dico_exemple)))


def get_bound_position(len_seq, context_kl, context_kr, context_pl, context_pr):  # à bien couvrir toutes les situations car il y a aussi des courts peptides
                                                                                  # et on envisage d'aller voir au moins jusqu'à 10 voisins
    """
    Bound of position that are valid according to the contextual window defined

    return:
        valid_interval: ensemble définissant les positions valides selon la fenetre de contexte local choisi
    """
    valid_interval = []
    max_left_window = max(context_kl, context_pl)
    max_right_window = max(context_kr, context_pr)
    for index in range(0, len_seq):
        if 0 <= index - max_left_window <= len_seq -1 and 0 <= index + max_right_window <= len_seq -1:
            valid_interval.append(index)
    valid_interval_ensemble = set(valid_interval) # conversion en ensemble

    return valid_interval_ensemble


def random_example_selection(list_example, dico_seq, valid_interval, context_kl, context_kr, context_pl, context_pr, example_selected_count):  
    """
    A ce stade on suppose qu on a selectionné un seed pour y piocher ce qu'on veut et on met cette fonction
    dans une boucle autant de fois qu'on doit prendre d'exemple de ce seed.

    context_kl, context_kr, context_pl, context_pr: pour récupérer le bon voisinage
    example_selected_count: pour le suivi du nombre d'exemple selectionner (pour en selectionner exactement ce qui a été déterminé pseudo- aléatoirement)

    return:
        list_example: liste de tuple. Chaque tuple est associé à un exemple qui sera utilisé pour la calcul du score de Brier.
    """

    list_pair_seq_name = tuple(dico_seq.keys())
    # récupération du poid de chaque paire de séquence de dico_seq
    weights_list = []  # random.choices(weights = list?/tuple?) --> choix de list pour le moment
    for key in dico_seq:
        weights_list.append(dico_seq[key]["nb valid aa couple (seq)"])  # à trouver un meilleur nom de variable

    pair_name = random.choices(list_pair_seq_name, weights = weights_list, k = 1)[0] 
    list_valid_position = dico_seq[pair_name]["valid position"]  # ensemble des positions valides


    final_list_valid_position = list(list_valid_position.intersection(valid_interval))   # je ne sais pas s'il y a un sens plus avantageux qu'un autre ...
                                                                                         # sachant que valid_interval est continue et pas forcément list_valid_position

    
    if final_list_valid_position:   # s'il existe au moins une position valide, piocher jusqu'à la trouver (au pire cas, ca devrait tout de meme etre rapide)
        # récupérer la paire dont on sait qu'on va en piocher un exemple
        pair_seq = dico_seq[pair_name]["sequences A, B"]

        # orienter la paire en 1 couple aléatoire
        proba = random.uniform(0, 1)
        if proba <= 0.5:
            seq_1, seq_2 = pair_seq
        else:
            seq_2, seq_1 = pair_seq      
            
        # random selection d'une position valide selon le context
        position_selected = random.choice(final_list_valid_position)

        
        # récupérer le couple et son voisinage et le mettre dans un liste de tuples
        example_selected = (seq_1[position_selected], seq_2[position_selected], seq_1[position_selected-context_kl: position_selected], 
                    seq_1[position_selected: position_selected+context_kr], seq_2[position_selected-context_pl: position_selected], seq_2[position_selected: position_selected+context_pr]) 

        list_example.append(example_selected)
        example_selected_count += 1
        #print("example_selected_count:", example_selected_count)

    return list_example, example_selected_count


def multi_random_example_selection(path_folder_seed, path_dico_exemple_complet, path_dico_seq, 
                                   context_kl, context_kr, context_pl, context_pr):
    """
    path_folder_seed: à relire tous les seed max une fois (à voir comment articuler)
    path_dico_exemple_complet: loader une fois le dico de l'info sur les seed et nb d'exemples à prendre de chacun
    path_dico_seq: path du folder contenant les dico_seq de chaque seed

    context_kl: considérer les acides aminés à gauche de l'acide aminé connu jusqu'à context_lk positions à gauche (kl = known left)
    context_kr: considérer les acides aminés à gauche de l'acide aminé connu jusqu'à context_rk positions à droite (kr = known right)

    context_pl: considérer les acides aminés à gauche de l'acide aminé à prédire jusqu'à context_lp positions à gauche (pl = to predict left)
    context_pr: considérer les acides aminés à gauche de l'acide aminé à prédire jusqu'à context_rp positions à droite (pr = to predict right)
    """
    t = Timer()
    t.start()
    
    # initialisation de la liste d'exemples
    list_example = []

    # load de dico_exemple une seule fois
    dico_exemple = np.load(path_dico_exemple_complet, allow_pickle='TRUE').item()

    # "load" des files fasta de seeds
    files = Path(path_folder_seed).iterdir()

    for file in files:
        accession_num =  folder.get_accession_number(file)
        if accession_num in dico_exemple: 
            nb_ex_test = dico_exemple[accession_num]["nb exemples tests"]
            if nb_ex_test != 0:
                # fixer l'intervalle d'indices compatibles selon len_align du seed et la fenetre définie
                len_align = dico_exemple[accession_num]["len alignment"]
                valid_interval_ensemble = get_bound_position(len_align, context_kl, context_kr, context_pl, context_pr)

                if valid_interval_ensemble != []: 
                    # en loader le dico_seq associé grace à l'accession_num
                    dico_seq = np.load(f"{path_dico_seq}/{accession_num}.seq.npy", allow_pickle='TRUE').item()
                    # je crois qu'il y a tjs des cas ultra atypiques qui peuvent encore se produire ... 

                    # selection des nb_ex_test exemples dans le seed
                    example_selected_count = 0
                    while example_selected_count < nb_ex_test:
                        list_example, example_selected_count = random_example_selection(list_example, dico_seq, valid_interval_ensemble,
                                                     context_kl, context_kr, context_pl, context_pr, example_selected_count)  
                    print("example_selected_count :",example_selected_count)

    t.stop("Selection des exemples tests")
    return list_example
                       
                    



if __name__ == '__main__':

    path_folder_seed = "/Users/pauline/Desktop/coupleSeqTest"                # chemin vers les seeds tests
    path_folder_pid = "/Users/pauline/Desktop/Overfitting_test/PID_couple"   # chemin vers les pid de seeds
    pid_inf = 62
    list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    
    # /Users/pauline/Desktop/data_Result/selection_test_brier_voisin  # dossier à crééer avant de faire n'importe quel test
    path_folder_dico_seq = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seq" # chemin ou on veut enregistrer les dico_seq
    path_file_dico_seed = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed" # chemin ou on veut enregistrer les dico_seed

    # # calcul du dico_seed et des dico_seq (1 par seed informatif) A CALCULER QU UNE SEULE FOIS POUR TOUS LES TESTS FUTUR, SUR LA FRACTION
    # DE PFAM D'INTERET
    multi_seeds_selection(path_folder_seed, path_folder_pid, pid_inf, list_residu, path_folder_dico_seq, path_file_dico_seed)



    # A CALCULER UNE SEULE FOIS PAR CONTEXTE SOUHAITÉ
    path_dico_seed_normalised = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed/seed_normalised.npy"
    nb_exemple_test = 10
    path_dico_exemple = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed"
    example_number_per_seed(path_dico_seed_normalised, nb_exemple_test, path_dico_exemple)

    path_dico_seq = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seq"
    path_dico_exemple_complet = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed/exemple_seed.npy"
    context_kl, context_kr, context_pl, context_pr = 2, 3, 1, 5

    list_example = multi_random_example_selection(path_folder_seed, path_dico_exemple_complet, path_dico_seq, 
                                   context_kl, context_kr, context_pl, context_pr)
    print(list_example)

