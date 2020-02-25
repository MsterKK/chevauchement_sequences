# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:50:57 2020

@author: kevin
"""

from FM_index import *

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
		print('----------')
		print('k',k,'l',l)
		print('count_ref',count_ref[seq[i]] )
		print('Rev_occ_ref',Rev_occ_ref[seq[i]][k-1])
		k = count_ref[seq[i]] + Rev_occ_ref[seq[i]][k-1]
		l = count_ref[seq[i]] + Rev_occ_ref[seq[i]][l]
		print('(k,l)',(k,l))
		if k > l :
			k = 1
			l = len(ref)
			z = z + 1
			print(z)
		D[i]=z
	return D


def alignement_inexact():
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


#test 
seq = "ATG"
ref = "ATGAGA"
l = len(ref)


ref_BWT, count_ref, occ_ref = FMindex(ref)

print('count_ref',count_ref)

Rev_ref = ref[::-1]
print('Rev_ref',Rev_ref)
Rev_ref_BWT = BWT(Rev_ref)
Rev_occ_ref = count_table(Rev_ref_BWT)
print('Rev_occ_ref',Rev_occ_ref)

D = calcul_D(seq, l, count_ref, Rev_occ_ref)

#print('ref_BWT',ref_BWT)
#print('Rev_ref_BWT',Rev_ref_BWT)
#D=[0,0,0,0]
#I = Inex_rec(seq, len(seq) - 1,1,1,l,D,count_ref,occ_ref)
#
#print(count_ref)
#print(occ_ref)

