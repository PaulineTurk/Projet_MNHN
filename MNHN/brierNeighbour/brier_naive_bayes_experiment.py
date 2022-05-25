import numpy as np
import matplotlib.pyplot as plt
import os

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]
sys.path.append(str(package_root_directory_MNHN))

import MNHN.brierNeighbour.selection_example as selection_example
from MNHN.utils.timer import Timer

def cube_loader(max_relative_distance, k_or_p, l_or_r, path_cube_folder):
    """
    max_relative_distance: indice le plus lointain dans quart de fenetre
    k_or_p: séquence d'origine (k) ou de destination (p)
    l_or_r: voisinage à gauche ou à droite
    """
    # initialisation de la liste des cubes pour un quart de fenetre contextuelle
    list_cube_quarter_window = []

    for i in range(1, max_relative_distance + 1): 
            if l_or_r == "l":
                path_cube = f"{path_cube_folder}/proba_cond_(-{i},{k_or_p}).npy" 
                list_cube_quarter_window.append(np.load(path_cube, allow_pickle='TRUE').item())

            if l_or_r == "r":  # je sais qu'un else marche mais if pour le moment
                path_cube = f"{path_cube_folder}/proba_cond_({i},{k_or_p}).npy"  # à vérifier le nom, sinon le simplifier ...
                list_cube_quarter_window.append(np.load(path_cube, allow_pickle='TRUE').item())

    return list_cube_quarter_window
    

def vecteur_from_cube(aa_1, contexte: list, list_cube: list, list_inclusion: list):
    """
    aa_1: acide amnié de départ
    contexte: liste des caractères d'une fraction de voisinage local
    list_cube: list des cubes pré)calculés et pré-chargés nécessaires pour ce contexte
    """
    len_list_inclusion = len(list_inclusion)
    list_vect = [np.array(len_list_inclusion*[1])]
    for index, aa_c in enumerate(contexte):
        if aa_c in list_inclusion:
            vect = []
            for aa in list_inclusion:
                vect.append(list_cube[index][aa_1][aa][aa_c])
            list_vect.append(np.array(vect))

    return list_vect


def unit_brier_naive_bayes(vect, aa_2, list_inclusion):
    """
    vect: vecteur de distribution de probabilité pré-calculé pour chaque exemple selon son contexte
    aa_2: acide aminé de destination
    """
    unit_brier = 0
    for index, aa in enumerate(list_inclusion):
        unit_brier += (vect[index] - int(aa_2 == aa))**2
    return unit_brier


def naive_bayes_brier(list_example, context_kl, context_kr, context_pl, context_pr, path_cube_folder, path_blosum_proba, list_inclusion):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_cube_folder: chemin du dossier contenant les cubes
    path_blosum_proba: chemin de la blosum probabiliste
    """
    t = Timer()
    t.start()

    # initialisation du score de Brier (Bayésien naif)
    score_brier_naive_bayes = 0

    # chargement des cubes et de la matrice de probabilité conditionnelle
    blosum_cond_proba = np.load(path_blosum_proba, allow_pickle='TRUE').item()
    list_cube_quarter_window_kl = cube_loader(context_kl, "k", "l", path_cube_folder)
    list_cube_quarter_window_kr = cube_loader(context_kr, "k", "r", path_cube_folder)
    list_cube_quarter_window_pl = cube_loader(context_pl, "p", "l", path_cube_folder)
    list_cube_quarter_window_pr = cube_loader(context_pr, "p", "r", path_cube_folder)

    for example in list_example:
        total_list_vect = []

        aa_1 = example[0]
        aa_2 = example[1]
        aa_c_kl = example[2]
        aa_c_kr = example[3]
        aa_c_pl = example[4]
        aa_c_pr = example[5]

        # initialisation du vecteur de distribution de probabilité de mutation 
        # de aa_1 en chacun des 20 aa standards
        vect_distribution = []
        for aa in list_inclusion:
            vect_distribution.append(blosum_cond_proba[aa_1][aa])
        
        total_list_vect.append(np.array(vect_distribution))
        #print("init total list with blosum proba:\n", total_list_vect)

        list_vect_kl = vecteur_from_cube(aa_1, aa_c_kl, list_cube_quarter_window_kl, list_inclusion)
        #print("list_vect_kl:\n", list_vect_kl)
        list_vect_kr = vecteur_from_cube(aa_1, aa_c_kr, list_cube_quarter_window_kr, list_inclusion)
        #print("list_vect_kr:\n", list_vect_kr)
        list_vect_pl = vecteur_from_cube(aa_1, aa_c_pl, list_cube_quarter_window_pl, list_inclusion)
        #print("list_vect_pl:\n", list_vect_pl)
        list_vect_pr = vecteur_from_cube(aa_1, aa_c_pr, list_cube_quarter_window_pr, list_inclusion)        
        #print("list_vect_pr:\n", list_vect_pr)
        
        total_list_vect = total_list_vect + list_vect_kl + list_vect_kr + list_vect_pl + list_vect_pr   # vect of arrays
        #print("total_list_vect:\n", total_list_vect)

        final_vector = np.prod(np.vstack(total_list_vect), axis=0) # pas encore normalisé
        final_vector_normalized = final_vector/np.sum(final_vector)

        #print("final_vector_normalized\n", final_vector_normalized)
        #print("sum vector normalized:", np.sum(final_vector_normalized))

        score_brier_naive_bayes += unit_brier_naive_bayes(final_vector_normalized, aa_2, list_inclusion)

    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example
    print("nombre d'exemples tirés:", nb_example)
    print("score_brier_naive_bayes   :", score_brier_naive_bayes)
    t.stop("Calcul score de Brier avec Naive Bayes")

    return score_brier_naive_bayes

  


if __name__ == '__main__':
    # selection des examples
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
    #selection_example.multi_seeds_selection(path_folder_seed, path_folder_pid, pid_inf, list_residu, path_folder_dico_seq, path_file_dico_seed)



    # A CALCULER UNE SEULE FOIS PAR CONTEXTE SOUHAITÉ
    #list_context = [(1, 0, 0, 0), (2, 0, 0, 0)]  # à définir en fonction des cubes disponibles
    list_context = [(1, 0, 0, 0)]  # à définir en fonction des cubes disponibles
    list_nb_example = [i for i in np.arange(10_000, 100_000, 1000)]

    path_dico_seed_normalised = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed/seed_normalised.npy"
    path_dico_exemple = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed"

    path_dico_seq = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seq"
    path_dico_exemple_complet = "/Users/pauline/Desktop/data_Result/selection_test_brier_voisin/seed/exemple_seed.npy"

    path_cube_folder = "/Users/pauline/Desktop/cube_test"
    path_blosum_proba = "/Users/pauline/Desktop/Overfitting_test/test_1/BlosumRes/BlosumRes_50_A/Blosum_proba_cond.npy"
    list_inclusion = list_residu

    path_image = "/Users/pauline/Desktop/data_Result"  # ou sauvegarder les graphes en sortie

    for context in list_context:
        context_kl, context_kr, context_pl, context_pr = context

        list_score_brier_bayes_naif = []

        for nb_exemple_test in list_nb_example:
            selection_example.example_number_per_seed(path_dico_seed_normalised, nb_exemple_test, path_dico_exemple)
            list_example = selection_example.multi_random_example_selection(path_folder_seed, path_dico_exemple_complet, path_dico_seq, 
                                    context_kl, context_kr, context_pl, context_pr)
            #print(list_example) 
            print("context :", context)
            print("nombre d'exemples demandés:", nb_exemple_test )
            score_brier_naive_bayes = naive_bayes_brier(list_example, context_kl, context_kr, context_pl, context_pr, path_cube_folder, path_blosum_proba, list_inclusion)
            list_score_brier_bayes_naif.append(score_brier_naive_bayes)
        print(list_score_brier_bayes_naif)
        plt.scatter(list_nb_example, list_score_brier_bayes_naif, label = f"{context_kl}, {context_kr}, {context_pl}, {context_pr}", alpha = 0.75, s = 5)
        plt.xlabel("Nombre d'exemples test")
        plt.ylabel('Brier Score')
        title = f"Score de Brier avec naive Bayes calculé sur {os.path.basename(path_folder_seed)}"
        plt.title(title)
        plt.legend()
        plt.savefig(f"{path_image}/{title}.png")