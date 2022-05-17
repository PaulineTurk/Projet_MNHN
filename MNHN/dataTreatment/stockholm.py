from Bio import AlignIO


def stockholmToFasta(path_stockholm, path_fasta) :
    """Convert a Stockholm file into a fasta file"""
    with open(path_stockholm, "r") as file_stockholm:
        with open(path_fasta, "w") as file_fasta:
            alignments = AlignIO.parse(file_stockholm, "stockholm")
            for alignment in alignments:
                AlignIO.write([alignment], file_fasta, "fasta")