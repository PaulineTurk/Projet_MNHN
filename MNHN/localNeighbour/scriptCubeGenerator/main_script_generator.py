import os
import shutil


def copyFileRenameComplete(path_folder, path_original_py, path_original_sh, list_delay_number):
    """
    path_folder: path to save the generated py and sh scripts
    path_original_py: path of the py script to complete
    path_original_sh: path of the sh script to complete
    list_delay_number: list of distant neighbor to generate cubes
    """

    if os.path.isdir(path_folder):
        shutil.rmtree(path_folder) 
    os.mkdir(path_folder)

    for delay_num in list_delay_number:
        for kp_SeqChoice in ["k", "p"]:

            # py file
            name_py_file = f"cube{delay_num}_{kp_SeqChoice}.py"
            py_to_execute = f"{path_folder}/{name_py_file}"
            shutil.copy(path_original_py, py_to_execute)

            file = open(py_to_execute, "r")
            list_line = file.readlines()
            list_line[9] = f"delay_num = {delay_num}\n"
            list_line[10] = f"kp_SeqChoice = '{kp_SeqChoice}'\n"

            file = open(py_to_execute, "w")
            file.writelines(list_line)
            file.close()


            # sh file
            name_sh_file = f"cube{delay_num}_{kp_SeqChoice}.sh"
            path_target = f"{path_folder}/{name_sh_file}"
            shutil.copy(path_original_sh, path_target)

            file = open(path_target, "r")
            list_line = file.readlines()
            list_line[3] = f"#SBATCH --job-name=cube{delay_num}_{kp_SeqChoice}\n"
            list_line[18] = f"#SBATCH --output=cube{delay_num}_{kp_SeqChoice}_10.out\n"
            list_line[43] = "cd /trinity/home/pturk/scriptCube\n"
            list_line[44] = f"python ./{name_py_file} $path_folder_pID $path_folder_data_split $path_new_folder\n"

            file = open(path_target, "w")
            file.writelines(list_line)
            file.close()






if __name__ == '__main__': 
    path_folder = '/Users/pauline/Desktop/Structuration/MNHN/blosumNeighbour/scriptCubeGenerator/scriptCube'
    path_original_py = '/Users/pauline/Desktop/Structuration/MNHN/blosumNeighbour/scriptCubeGenerator/scriptCubeRef.py'
    path_original_sh = '/Users/pauline/Desktop/Structuration/MNHN/blosumNeighbour/scriptCubeGenerator/scriptCubeRef.sh'
    list_delay_number = [k for k in range(-10, 11) if k!=0]

    copyFileRenameComplete(path_folder, path_original_py, path_original_sh, list_delay_number)
