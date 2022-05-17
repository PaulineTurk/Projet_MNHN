def LowerToUpper(path_file, path_file_corrected):
    """
    Convert all the lowercase residu into uppercase.

    path_file: path of the fasta file to correct
    path_file_corrected; path of the fasta file corrected
    """
    with open(path_file, "r") as file:
        with open(path_file_corrected, "w") as file_corrected:
            for line in file:
                if line[0] != ">":   
                    line = line.upper()
                file_corrected.write(line)