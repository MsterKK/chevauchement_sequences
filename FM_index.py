# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:34:16 2020

Certaines parties du code sont inspirées d'algorithmes trouvés sur Github:
- notamment https://gist.github.com/Puriney/6324227
"""


alphabet = ["A","C","G","T"]

def suffix_array(seq):
	"""Fonction qui construit la table des suffixes associée à la chaîne de caractères seq
	--- 
	output: 
		la liste contenant les rangs des suffixes triés
	"""
	#on ajoute le caractère $ à la fin de la chaîne de caractères
	seq = seq +'$'
	#tri de la liste composée des suffixes ainsi que du rang du suffixe
	liste_suffixes = sorted([(seq[i:],i) for i in range(len(seq))])
	#Construction de la liste contenant les rangs de la liste des suffixes
	S = [liste_suffixes[i][1] for i in range(len(seq))]
	return S

def BWT(seq):
	'''Fonction qui réalise la transformation de Burrows-Wheeler à partir de la table des suffixes
	---
	output:
		bw: string
	'''
	bw = ''
	#ajout du caractère $ à la chaîne de caractère originale
	seq2 = seq +'$'
	#récuperation des rangs de la table des suffixes pour construire la BWT
	for c in suffix_array(seq):
		bw += seq2[c]
	return bw

def count_occ(seq):
	"""Fonction qui, pour un caractère c, compte les occurences de tous les caractères lexicographiquement plus 
	petits que c
	--- 
	output
		dic: dictionnaire
	"""
	
	dic = {}
	seq = seq + "$"
	alphabet = ["$","A","C","G","T"]
	#initialisation du dictionnaire
	for letter in alphabet:
		dic[letter]=0
	#comptage du nombre de lettre de chaque type
	for letter in seq:
		dic[letter] += 1
	#inverse la liste dans l'ordre lexicographique
	alphabetRev = alphabet[::-1]
	nb_carac = len(seq)
	#ajout du nb de caractères précédant chaque lettre en partant de la derniere lettre dans l'ordre lexicographique
	for letter in alphabetRev:
			dic[letter] = nb_carac - dic[letter]
			nb_carac = dic[letter]
	return dic

def count_table(seq):
	"""Fonction qui renvoie la table des occurences d'une sequence, çad un dictionnaire qui, pour chaque lettre de l'alphabet,
	contient la ligne du tableau des occurences (voir rapport pour le détail)
	-
	output
		dic_table: dictionnaire contenant des listes pour chaque clef
	"""
	dic_table = {}
	alphabet =["$","A","C","G","T"]
	taille_seq = len(seq)
	#initialisation du dictionnaire
	for letter in alphabet:
		dic_table[letter] = [0]*taille_seq
	#parcours de la chaine de caractères pour trouver les occurences de chaque caractère
	for k in range(taille_seq):
		c = seq[k]
		dic_table[c][k] = 1
	#concatenation des occurences trouvées afin de former le tableau
	for letter in alphabet:
		for k in range(1,taille_seq):
			dic_table[letter][k] += dic_table[letter][k-1]
	return dic_table


def FMindex(seq):
	"""Fonction qui renvoie la transformation de Burrows-Wheeler ainsi que la table C[c] (ici un dictionnaire)
	contenant les occurences des caractère lexicographiquement plus petit que c, et la table des occurences (voir rapport
	pour détail)
	output:
		(bwt,dic)
		où - bw: string
		- dic: dictionnaire
		- tables_occurences
	"""
	bw = BWT(seq)
	dic = count_occ(seq)
	table_occurences = count_table(seq)
	return bw, dic, table_occurences
	
	
	
	
	