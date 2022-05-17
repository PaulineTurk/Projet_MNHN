import numpy as np
import math


list_1 = np.array([1, 3, 4])
print(list_1)
list_2 = [4, 2, 5]

vector_proba = np.multiply(list_1, list_2) # multiplication of 2 arrays or lists
print(vector_proba)
type = type(vector_proba)
print(type)
print(vector_proba[0])






print("intersection")


list_residu1 = set([1, 3, 2, 4])  # set to define ensemble (union, intersection, ... get possible)
list_residu2 = set([1, 3,  4, 5, 7])
print(list_residu1.intersection(list_residu2))
intersection = list(list_residu1.intersection(list_residu2))
print(intersection)
index = intersection[2]
print(index)




div = 3
dico = {"A":1, "B":2, "C": 3}
dico_norm = {k: v/div for k, v in dico.items()}
print("dico: ", dico)
print("dico_norm: ", dico_norm)



test = {}
test[(3, 2)] = {}
test[(3, 2)]["oko oko"] = 4
test[(3, 2)]["OUI OUI"] =['J', 'OK'] 
for key in test:
    type_key = type(key)
    print(f"type key {type_key}")


import random
list_mot = ['cnfiogcuhdnigundfhv', 'gundfhv']
for i in range(10):
    print(random.choice(list_mot))


numberList = [111, 222, 333, 444, 555]
print(random.choices(numberList, weights=(10, 20, 30, 40, 50), k=5))  # use relative weights
# Output [555, 222, 555, 222, 555]




for i in range(10):
    print(random.randint(1, 10))


print(math.modf(0.24))
print(math.modf(3.24))


print("traitement Pfam pour selection des exemples d'apprentissage")
print((19600*155/51)/(2*3600))



import blosum as bl
blosum_ref = bl.BLOSUM(62) 
print(blosum_ref)