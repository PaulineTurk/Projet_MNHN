import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

import MNHN.utils.folder as folder
import MNHN.utils.fastaReader as fastaReader
from MNHN.utils.timer import Timer

import MNHN.brier.predictor as predictor 


def brier_score_predictor_01(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62):   
    """
    Compute Brier Score on files with the predictor 0/1.
    """
    t = Timer()
    t.start()

    # intialisation
    brier_score = 0
    count = 0

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        seed = fastaReader.read_multi_fasta(file)
        brier_score, count = predictor.predictor_01(seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf)
    if count != 0:
        brier_score = brier_score/count
    else:
        brier_score = None
        print("Brier Score not computable because 0 prediction was done")

    t.stop("Brier Score with predictor 0/1 (worst case scenario)")

    return brier_score


def brier_score_predictor_perfect(path_folder_fasta, path_folder_pid, list_residu, pid_inf = 62):   
    """
    Compute Brier Score on files with the perfect predictor.
    """
    t = Timer()
    t.start()

    # intialisation
    brier_score = 0
    count = 0

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        seed = fastaReader.read_multi_fasta(file)
        brier_score, count = predictor.predictor_perfect(seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf)
    
    if count != 0:
        brier_score = brier_score/count
    else:
        brier_score = None
        print("Brier Score not computable because 0 prediction was done")

    t.stop("Brier Score with perfect predictor")

    return brier_score





def brier_score_matrix(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62):  
    """
    Compute Brier Score on files with a predictor from the list: ["Blosum Predictor", "Equiprobable Predictor",  
                                                                  "Stationary Predictor", "Identity Predictor"]
    """
    t = Timer()
    t.start()

    # intialisation
    brier_score = 0
    count = 0

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        seed = fastaReader.read_multi_fasta(file)
        brier_score, count = predictor.brier_matrix(unit_Brier, seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf)
    if count != 0:
        brier_score = brier_score/count
    else:
        brier_score = None
    t.stop("Brier Score")   # intéret de donner un nom à chaque prédicteur 
    
    return brier_score


def brier_score_matrix_v2(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62):  
    """
    Compute Brier Score on files with a predictor from the list: ["Blosum Predictor", "Equiprobable Predictor",  
                                                                  "Stationary Predictor", "Identity Predictor"]
    For each amino acid and the total brier score
    """
    t = Timer()
    t.start()

    # intialisation
    brier_score = {}
    count = {}
    for aa in list_residu:
        brier_score[aa] = 0
        count[aa] = 0

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = folder.get_accession_number(file)
        seed = fastaReader.read_multi_fasta(file)
        brier_score, count = predictor.brier_matrix_v2(unit_Brier, seed, accession_num, path_folder_pid, brier_score, count, list_residu, pid_inf)
    for aa in count:
        if count[aa] != 0:
            brier_score[aa] /= count[aa] 
        else:
            brier_score[aa]  = 0
    t.stop("Brier Score")   # intéret de donner un nom à chaque prédicteur 

    total_count = sum(count.values())

    fraction_count = {}
    for aa in count:
        fraction_count[aa] = count[aa]/total_count

    total_brier_score = 0
    for aa in brier_score:
        total_brier_score += brier_score[aa] * fraction_count[aa]
    print("Brier Score for each amino-acid:\n", brier_score)
    print("Count of predictions for each amino-acid:\n", count)
    print("Total Brier Score:\n", total_brier_score)
    
    return brier_score, count, total_brier_score


# def overfittingTest(path_data, percentage_A, path_pid, path_BlosumRes, name_data_train, name_data_test, list_residu):
#     folder_fasta_train = f"{path_data}/PfamSplit_{str(percentage_A)}/{name_data_train}" 
#     folder_fasta_test = f"{path_data}/PfamSplit_{str(percentage_A)}/{name_data_test}" 
#     print("folder_fasta_train:", folder_fasta_train)
#     print("folder_fasta_test:", folder_fasta_test)

#     path_matrix_cond_proba = f"{path_BlosumRes}/BlosumRes_{str(percentage_A)}_{name_data_train}/Blosum_proba_cond.npy" 


#     predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba, list_residu)
#     Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta_test, path_pid, unit_Brier_Blosum, list_residu, pid_inf = 62) 
#     print("Blosum Predictor Brier Score:", Brier_Score_global)
#     print("")