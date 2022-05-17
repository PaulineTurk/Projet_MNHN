from pickle import FALSE
from pathlib import Path
from sklearn.model_selection import train_test_split

import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))

from MNHN.utils.timer import Timer
from MNHN.utils.folder import creatFolder, getAccessionNb

from MNHN.dataTreatment.descriptionTool import dataCount, barPlot
from MNHN.dataTreatment.lowerToUpperTool import LowerToUpper
from MNHN.dataTreatment.stockholm import stockholmToFasta
from MNHN.dataTreatment.pid import pidCoupleSeq
from MNHN.dataTreatment.redondancy import nonRedundant
from MNHN.dataTreatment.split import addFileFromFolder

def stockholmSeparator(path_file, path_folder_save):
    """
    Separate a multiStockholm file into monoStockholm files.

    path_file_name: path of the multiStockholm file
    path_folder_save: path of th folder where the monoStockholm files generated are saved 
    """
    t = Timer()
    t.start()

    input_handle = open(path_file)
    creatFolder(path_folder_save)

    # collect the accession number of each alignment
    list_accession_num = []
    for l in input_handle:
        if l[0:7] == "#=GF AC":
            init_accession_num = l.index('PF')
            accession_num = l[init_accession_num: -1]
            list_accession_num.append(accession_num)
    input_handle.close()
    nbre_file = len(list_accession_num)

    # generate the monoStockholm files named after the accession number's alignment
    input_handle = open(path_file)
    file_out_nbre = 0
    path_file_out = f"{path_folder_save}/{list_accession_num[file_out_nbre]}.stockholm"
    output_handle = open(path_file_out, "w")

    for l in input_handle:
        output_handle.write(l)
        if l[0:2] == "//" and file_out_nbre <= nbre_file - 2: # avoid generating an empty file at the end
            output_handle.close()
            file_out_nbre += 1
            path_file_out = f"{path_folder_save}/{list_accession_num[file_out_nbre]}.stockholm"
            output_handle = open(path_file_out, "w")

    output_handle.close()
    input_handle.close()
    t.stop("Separation of the multiStockholm file into monoStockholm files")





def multiStockholmToFasta(path_folder_stockholm, path_folder_fasta):
    """
    Convert Stockholm files into Fasta files
    """
    t = Timer()
    t.start()

    creatFolder(path_folder_fasta)

    files_stockholm = Path(path_folder_stockholm).iterdir()
    for file_stockholm in files_stockholm:
        accession_num = getAccessionNb(file_stockholm)
        path_file_fasta = f"{path_folder_fasta}/{accession_num}.fasta"
        stockholmToFasta(file_stockholm, path_file_fasta)
    t.stop("Conversion of Stockholm files into Fasta files")







def dataVisualisation(path_data):
    """
    path_folder: path of the folder containing fasta files to describe
    """
    t = Timer()
    t.start()

    residu_count, total_residu = dataCount(path_data)
    #barPlot(path_data, residu_count, "count")

    residu_percentage_distribution = {k: round(100*v / total_residu, 2) for k, v in residu_count.items()}
    barPlot(path_data, residu_percentage_distribution, "percentage")

    t.stop("Time for data description")






def multiLowerToUpper(path_data, path_data_corrected):
    """
    Convert all the lowercase residu into uppercase.

    path_data: path of the folder of fasta file to correct
    path_data_corrected: folder in which the fasta file corrected are saved
    """
    t = Timer()
    t.start()

    creatFolder(path_data_corrected)

    files = Path(path_data).iterdir()
    for file in files:
        accession_num = getAccessionNb(file)
        path_file_corrected = f"{path_data_corrected}/{accession_num}.fasta.upper"
        LowerToUpper(file, path_file_corrected)
    t.stop("Correction upper files")







def savePid(path_folder_fasta, path_folder_pid):
    """
    For each fasta file in path_folder_fasta, compute the dictionary of pid 
    for each couple of sequences and save it in path_folder_pid
    """
    t = Timer()
    t.start()

    creatFolder(path_folder_pid)

    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = getAccessionNb(file)
        path_file_pid = f"{path_folder_pid}/{accession_num}.pid"
        pidCoupleSeq(file, path_file_pid)
        
    t.stop("Compute and save the pId files")
    



def savePIdNonRedondant(path_folder_fasta, path_folder_fasta_non_redondant, list_residu, pid_sup = 99):
    """
    Create a folder of fasta files with the redondant issue corrected

    path_folder_fasta: original fasta folder
    path_folder_fasta_non_redondant: fasta folder corrected
    list_residu: list of valid residu to evaluate the corrected len of each sequence
    pid_sup: percentage of identity for the clustering
    """
    t = Timer()
    t.start()
    creatFolder(path_folder_fasta_non_redondant)
    
    files = Path(path_folder_fasta).iterdir()
    for file in files:
        accession_num = getAccessionNb(file)
        path_fasta_non_redondant = f"{path_folder_fasta_non_redondant}/{accession_num}.fasta.nonRedundant"
        nonRedundant(file, path_fasta_non_redondant, list_residu, pid_sup)
    t.stop("Compute and save non-redundant files")



def dataSplit(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B):
    """
    Random split files in categorie A and B according to percentage_A

    path_folder_data: folder to split
    path_folder_data_split: folder created where the data splitted are saved
    percentage_A: percentage of files in category A
    name_data_A: name of the folder in path_folder_data_split for the files in category A
    name_data_B: name of the folder in path_folder_data_split for the files in category B
    """
    t = Timer()
    t.start()

    creatFolder(path_folder_data_split)
    path_folder_data_A = f"{path_folder_data_split}/{name_data_A}"
    path_folder_data_B = f"{path_folder_data_split}/{name_data_B}"
    creatFolder(path_folder_data_A)
    creatFolder(path_folder_data_B)

    files = Path(path_folder_data).iterdir()
    data_name = []

    for file_path in files:
        file_name = str(file_path).split("/")[-1]
        data_name.append(file_name)

    fraction_A = percentage_A/100
    x_A ,x_B = train_test_split(data_name, train_size = fraction_A)  
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html

    for file_name in x_A:
        addFileFromFolder(file_name, path_folder_data, path_folder_data_A, "A")
    
    for file_name in x_B:
        addFileFromFolder(file_name, path_folder_data, path_folder_data_B, "B")

    t.stop("Split data_total in data_A and data_B")