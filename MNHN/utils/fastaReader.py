def read_multi_fasta(nomFi):
    seq=""
    nom=""
    lesSeq=[]
    with open(nomFi,"r") as f:  
        for l in f:
            if l[0] == '>':
                if seq != "":
                    tmp=(nom,seq)
                    lesSeq.append(tmp)
                nom=l[1:-1]
                seq=""
            else:
                seq = seq + l.strip()
    if seq != "":
        tmp=(nom,seq)
        lesSeq.append(tmp)
    return lesSeq