import numpy as np


def cube_loader(max_relative_distance, k_or_p, l_or_r, path_cube_folder):
    """
    max_relative_distance: indice le plus lointain dans quart de fenetre
    k_or_p: séquence d'origine (k) ou de destination (p)
    l_or_r: voisinage à gauche ou à droite
    """
    if max_relative_distance != 0: 
        for i in range(1, max_relative_distance): # prob si context_kl == 1 ? ajouter +1?
            if l_or_r == "l":
                path_cube = f"{path_cube_folder}/proba_cond_(-{i},{k_or_p}).npy"  # à vérifier le nom, sinon le simplifier ... à changer nom dans mon calculateur de cube
            if l_or_r == "r":  # je sais qu'un else marche mais if pour le moment
                path_cube = f"{path_cube_folder}/proba_cond_({i},{k_or_p}).npy"  # à vérifier le nom, sinon le simplifier ...

            np.load(path_cube, allow_pickle='TRUE').item() 

        #à donner des noms !!!!! pour pouvoir les appeler
    




def naive_bayes_brier(list_example, context_kl, context_kr, context_pl, context_pr, path_cube_folder, path_blosum_proba):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_cube_folder: chemin du dossier contenant les cubes
    path_blosum_proba: chemin de la blosum probabiliste
    """
    # charger les cubes nécessaires 1 seule fois, à donner des noms !!!!! pour pouvoir les appeler
    cube_loader(context_kl, "k", "l", path_cube_folder)
    cube_loader(context_kr, "k", "r", path_cube_folder)
    cube_loader(context_pl, "p", "l", path_cube_folder)
    cube_loader(context_pr, "p", "r", path_cube_folder)

    # et charger la blosum probabiliste (non contextuelle) dans tous les cas 
    # car possibilité de n'avoir pour voisin que des gaps
    np.load(path_blosum_proba, allow_pickle='TRUE').item() 


    for example in list_example:
        for index in range(2, 5): # chaque 4 quarts de fenetre
            pass
  


