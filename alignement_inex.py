# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:50:57 2020

@author: kevin
"""

import FM_index as fm

def calcul_D(seq, l, C_ref, Rev_occ_ref):
	"""Calcule l'array D(.) -- 
	D[i] contient le nombre minimal de mismatch entre la sequence de reference et seq[0,i]
	---
	input: 
		seq: string -- chaine de caractères que l'on veut comparer à une séquence de référence
		l: int -- longueur de la séquence de référence
		C_ref: table C des caractères lexicographiquement plus petit (générée avec count_occ) de la sequence de reference
		Rev_occ_ref: table des occurences (générée avec count_table) 
			de la BWT de l'inverse de la sequence de reference
	"""
	k = 1
	z = 0
	taille_ref = l
	taille_seq = len(seq)
	
	#on crée une case de trop à D[i] afin que D[0] = 0 ce qui permet l'arrêt de la fonction récursive Inex_rec
	D = [0] * (taille_seq+1)

	for i in range(taille_seq):
		k = C_ref[seq[i]] + Rev_occ_ref[seq[i]][k-1]+1
		l = C_ref[seq[i]] + Rev_occ_ref[seq[i]][l]
		if k > l :
			k = 1
			l = taille_ref
			z = z + 1
		D[i+1]=z
	return D

def alignement_inexact(query, seq_ref, C_ref, occ_ref, Rev_occ_ref, z):
	"""Fonction qui réalise l'alignemement de la sequence query sur la sequence de reférence
	avec un nombre maximal de z mismatchs
	--- Inputs ---
	query: string -- séquence à aligner sur la séquence de référence
	C_ref: table C de la séquence de référence
	occ_ref: tableau des occurences de la BWT de séquence de référence
	Rev_occ_ref: tableau des occurences de l'inverse de la BWT de la séquence de référence
	z: int -- nombre maximal de mismatchs permis
	--- Output ---
	I: ensemble contenant les intervalles SA où query match avec les préfixes 
	du suffixe array de la séquence de référence
	"""
	
	#initialisation des variables
	l = len(seq_ref)
	
	#Calcul du D(.) array
	D = calcul_D(query, l, C_ref, Rev_occ_ref)
	I = Inex_rec(query, len(query) - 1, z, 1, l, D,C_ref,occ_ref)
	
	return I

def Inex_rec(seq, i, z, k, l, D, C_ref, occ_ref):
	"""
	seq = query : séquence à analyser
			str
	i : index du dernier caractère considéré dans la sous-seq en partant de la fin
	z : nombre de mismatch
	k : borne inf de l'intervalle SA
	l : borne sup de l'intervalle SA
	D : matrice LB nb mismatch
	C_ref : table C
	occ_ref : table des occurences pour chaque caractère dans la BWT de seq_ref
	"""
	#Dans le cas où le nombre minimal de mismatch entre seq[0,i] est supérieur au nombre
	#maximal de mismatch z permis par l'utilisateur, l'algorithme renvoie une liste vide
	if z < D[i+1]:
		return set()  #retourne un ensemble
	
	#Dans le cas où i, la longueur du préfixe de seq que l'on regarde, est inférieure à 0,
	#on renvoie l'intervalle SA => Condition d'arrêt lorsque seq a été complètement parcourue
	if i < 0:
		new_set = set()
		new_set.add((k,l))
		return new_set
	
	I = set()
	
	#Appel récursif de la fonction permettant de gérer une délétion dans seq par rapport la sequence de réf
	I = I.union(Inex_rec(seq, i - 1, z - 1, k, l, D, C_ref, occ_ref))
	
	#on parcourt les lettres de l'alphabet, en recherche d'un match
	for letter in ['A','C','G','N','T']:
		k_bis = C_ref[letter] + occ_ref[letter][k-1]+1 
		l_bis = C_ref[letter] + occ_ref[letter][l]
		
		#On vérifie que l'intervalle SA où l'on regarde n'est pas vide
		if k_bis <= l_bis:

			#appel récursif de la fonction permettant de gérer une déletion dans la sequence de ref par rapport à seq
			I = I.union(Inex_rec(seq, i , z - 1, k_bis, l_bis, D, C_ref, occ_ref))

			#si on a un match, on continue l'algorithme en décrémentant i
			if letter == seq[i] and letter != 'N':
				I = I.union(Inex_rec(seq, i - 1, z, k_bis, l_bis, D, C_ref, occ_ref))
				
			#si on a pas de match, on continue l'algorithme en décrémentant i et z (nombre de mismatchs)
			else: 
				I = I.union(Inex_rec(seq, i - 1, z - 1, k_bis, l_bis, D, C_ref, occ_ref))
	return I


class Alignement_inex() :
	
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
		self.D = [0] * len(self.query+1)
		
		for i in range(len(self.seq_ref)):
			k = self.ref_C[self.query[i]] + self.ref_rev_occ[self.query[i]][k - 1] + 1
			l = self.ref_C[self.query[i]] + self.ref_rev_occ[self.query[i]][l]
			if k > l :
				k = 1
				l = len(self.seq_ref)
				z = z + 1
			self.D[i+1] = z
	
	def calcul_D2(self) :
		"""Alternative au calcul de la matrice D
		Plus lente, mais plus compréhensible, n'utilise pas BW
		"""
		z = 0 #nb de mismatchs
		j = 0 #index de départ de la 
		self.D = [0]*len(self.query)
		for i in range(len(self.query)) :
			if self.query[j:i+1] not in self.seq_ref :
				z += 1
				j = i + 1
			self.D[i] = z
	
	def alignement_inexact(self, z):
		"""Fonction qui réalise l'alignemement de la sequence query sur la sequence de reférence
		avec un nombre maximal de z mismatchs
		--- Inputs ---
		query: string -- séquence à aligner sur la séquence de référence
		C_ref: table C de la séquence de référence
		occ_ref: tableau des occurences de la BWT de séquence de référence
		Rev_occ_ref: tableau des occurences de l'inverse de la BWT de la séquence de référence
		z: int -- nombre maximal de mismatchs permis
		--- Output ---
		I: ensemble contenant les intervalles SA où query match avec les préfixes 
		du suffixe array de la séquence de référence
		"""
		
		#Calcul du D(.) array
		self.calcul_D()
		I = Inex_rec(self.query, len(self.query) - 1, z, 1, len(self.seq_ref))
		
		return I
	
	def Inex_rec(seq, i, z, k, l):
		"""
		seq = query : séquence à analyser
				str
		i : index du dernier caractère considéré dans la sous-seq en partant de la fin
		z : nombre de mismatch
		k : borne inf de l'intervalle SA
		l : borne sup de l'intervalle SA
		D : matrice LB nb mismatch
		C_ref : table C
		occ_ref : table des occurences pour chaque caractère dans la BWT de seq_ref
		"""
		#Dans le cas où le nombre minimal de mismatch entre seq[0,i] est supérieur au nombre
		#maximal de mismatch z permis par l'utilisateur, l'algorithme renvoie une liste vide
		if z < self.D[i+1]:
			return set()  #retourne un ensemble
		
		#Dans le cas où i, la longueur du préfixe de seq que l'on regarde, est inférieure à 0,
		#on renvoie l'intervalle SA => Condition d'arrêt lorsque seq a été complètement parcourue
		if i < 0:
			new_set = set()
			new_set.add((k,l))
			return new_set
		
		I = set()
		
		#Appel récursif de la fonction permettant de gérer une délétion dans seq par rapport la sequence de réf
		I = I.union(self.Inex_rec(seq, i - 1, z - 1, k, l))
		
		#on parcourt les lettres de l'alphabet, en recherche d'un match
		for letter in ['A','C','G','N','T']:
			k_bis = self.ref_C[letter] + self.ref_occ[letter][k-1]+1 
			l_bis = self.ref_C[letter] + self.ref_occ[letter][l]
			
			#On vérifie que l'intervalle SA où l'on regarde n'est pas vide
			if k_bis <= l_bis:
	
				#appel récursif de la fonction permettant de gérer une déletion dans la sequence de ref par rapport à seq
				I = I.union(self.Inex_rec(seq, i , z - 1, k_bis, l_bis))
	
				#si on a un match, on continue l'algorithme en décrémentant i
				if letter == seq[i] and letter != 'N':
					I = I.union(self.Inex_rec(seq, i - 1, z, k_bis, l_bis))
					
				#si on a pas de match, on continue l'algorithme en décrémentant i et z (nombre de mismatchs)
				else: 
					I = I.union(self.Inex_rec(seq, i - 1, z - 1, k_bis, l_bis))
		return I 




#test 
seq = "GNC"
ref = "ATGTAGCNAGGTC"
l = len(ref) + 1
z = 1

ref_BWT, C_ref, occ_ref = fm.FMindex(ref)

print('C_ref',C_ref)
print('occ_ref', occ_ref)

Rev_ref = ref[::-1]
print('Rev_ref',Rev_ref)
Rev_ref_BWT, Rev_C, Rev_occ_ref = fm.FMindex(Rev_ref)
print('Rev_occ_ref',Rev_occ_ref)

D = calcul_D(seq, l, C_ref, Rev_occ_ref)

print('D',D)

I = Inex_rec(seq, len(seq) - 1,1,1,l,D,C_ref,occ_ref)
print(' -- end of algorithm --')
print('I',I)

I = alignement_inexact(seq, ref, C_ref, occ_ref, Rev_occ_ref, z)
print ('I',I)
