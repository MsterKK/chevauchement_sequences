# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:50:57 2020

@author: kevin
"""

import FM_index as fm

def calcul_D(seq, l, count_ref, Rev_occ_ref):
	"""Calcule l'array D(.) -- 
	D[i] contient le nombre minimal de mismatch entre la sequence de reference et seq[0,i]
	---
	input: 
		seq: string -- chaine de caractères que l'on veut comparer à une séquence de référence
		l: int -- longueur de la séquence de référence
		count_ref: table C des caractères lexicographiquement plus petit (générée avec count_occ) de la sequence de reference
		Rev_occ_ref: table des occurences (générée avec count_table) 
			de la BWT de l'inverse de la sequence de reference
	"""
	k = 1
	z = 0
	D = [0] * len(seq)

	for i in range(len(seq)):
		k = count_ref[seq[i]] + Rev_occ_ref[seq[i]][k-1]+1
		l = count_ref[seq[i]] + Rev_occ_ref[seq[i]][l]
		if k > l :
			k = 1
			l = len(ref)
			z = z + 1
		D[i]=z
	return D



def alignement_inexact(query, seq_ref):
    pass

def Inex_rec(seq, i, z, k, l, D, count_ref, occ_ref):
	print('longueur',i)
    
	#Dans le cas où le nombre minimal de mismatch entre seq[0,i] est supérieur au nombre
	#maximal de mismatch z permis par l'utilisateur, l'algorithme renvoie une liste vide
	if z < D[i]:
		return []
	#Dans le cas où i, la longueur du préfixe de seq que l'on regarde, est inférieure à 0,
	#on renvoie l'intervalle SA
	if i < 0:
		print('k,l',(k,l))
		return [k,l]
	I = []
	#I = I + [Inex_rec(seq, i - 1 , z - 1, k, l, D, count_ref, occ_ref)]
	print('1',I)
	#on parcourt les lettres de l'alphabet, en recherche d'un match
	for letter in ['A','C','G','T']:
#		print('l',l)
		k = count_ref[letter] + occ_ref[letter][k-1]
		l = count_ref[letter] + occ_ref[letter][l-1]
#		print((k,l))
#		print('lettre1',letter)
		#On vérifie que l'intervalle où l'on regarde n'est pas vide
		if k <= l:
			print('lettre2',letter)
			#si on a un match, on continue l'algorithme en décrémentant i
			if letter == seq[i]:
				I = I + Inex_rec(seq, i - 1, z, k, l, D, count_ref, occ_ref)
			#si on a pas de match, on continue l'algorithme en décrémentant i et z (nombre de mismatchs)
			else: 
				I = I + Inex_rec(seq, i - 1, z - 1, k, l, D, count_ref, occ_ref)
	print(I)
	return I








class Alignement() :
    
    def __init__(self, seq_ref, query) :
        
        # Initialisation des constants
        self.seq_ref = seq_ref
        self.query = query
        
        # Prétraitement de la séquence de référence
        self.pretraitement()
    
    def pretraitement(self) :
        
        # FMindex : BW, table de C, table d'occurence de BW
        self.ref_BWT, self.ref_C, self.ref_occ = fm.FMindex(self.seq_ref)

        # Reverse de seq_ref
        self.seq_ref_rev = self.seq_ref[::-1]
        
        # Table de 
        self.ref_rev_BWT = fm.BWT(self.ref_rev, fm.suffix_array(self.ref_rev)[1])
        self.ref_rev_occ = fm.count_table(self.ref_rev_BWT)
        
    
    def calcul_D(self) :
        """Calcule l'array D(.) -- 
    	D[i] contient le nombre minimal de mismatch entre la sequence de reference et seq[0,i]
    	---
    	input: 
    		seq: string -- chaine de caractères que l'on veut comparer à une séquence de référence
    		l: int -- longueur de la séquence de référence
    		count_ref: table C des caractères lexicographiquement plus petit (générée avec count_occ) de la sequence de reference
    		Rev_occ_ref: table des occurences (générée avec count_table) 
    			de la BWT de l'inverse de la sequence de reference
    	"""
        k = 1
        l = len(self.seq_ref)
        z = 0
        self.D = [0] * len(self.query)
        
        for i in range(len(self.seq_ref)):
            k = self.ref_C[self.query[i]] + self.ref_rev_occ[self.query[i]][k - 1] + 1
            l = self.ref_C[self.query[i]] + self.ref_rev_occ[self.query[i]][l]
            if k > l :
                k = 1
                l = len(self.seq_ref)
                z = z + 1
            self.D[i] = z
    
    def calcul_D2(self) :
        z = 0
        j = 0
        self.D = [0]*len(self.query)
        for i in range(len(self.query)) :
            if self.query[j:i+1] not in self.seq_ref :
                z += 1
                j = i + 1
            self.D[i] = z

    def Inex_rec(self, i, z, k, l) :
        
        #Dans le cas où le nombre minimal de mismatch entre seq[0,i] est supérieur au nombre
        #maximal de mismatch z permis par l'utilisateur, l'algorithme renvoie une liste vide
        if z < self.D[i] :
            return {}
        #Dans le cas où i, la longueur du préfixe de seq que l'on regarde, est inférieure à 0,
    	#on renvoie l'intervalle SA
        if i < 0 :
            return [k,l]
        
        
    def alignement_inexacte(self, z) :
        
        self.calcul_D()
        self.Inex_rec(len(self.query) - 1, z, 1, len(self.seq_ref) - 1) 





#test 
seq = "ATGAGA"
ref = "ATGCTA"
l = len(ref)


ref_BWT, count_ref, occ_ref = fm.FMindex(ref)

print('count_ref',count_ref)

Rev_ref = ref[::-1]
print('Rev_ref',Rev_ref)
Rev_ref_BWT = fm.BWT(Rev_ref, fm.suffix_array(Rev_ref)[1])
Rev_occ_ref = fm.count_table(Rev_ref_BWT)
print('Rev_occ_ref',Rev_occ_ref)

D = calcul_D(seq, l, count_ref, Rev_occ_ref)

#print('ref_BWT',ref_BWT)
#print('Rev_ref_BWT',Rev_ref_BWT)
#D=[0,0,0,0]
#I = Inex_rec(seq, len(seq) - 1,1,1,l,D,count_ref,occ_ref)
#
#print(count_ref)
#print(occ_ref)

