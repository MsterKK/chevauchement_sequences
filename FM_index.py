# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:34:16 2020

@author: kevin
"""


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
	"""Fonction qui compte les occurences de chaque caractères de seq
	--- 
	output
		dic: dictionnaire
	"""
	
	dic = {}
	seq = seq + "$"
	#alphabet trié utilisé pour former la chaîne seq
	alphabet = sorted(set(seq))
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


def FM_index(seq):
	"""Fonction qui renvoie la transformation de Burrows-Wheeler ainsi que la table C (ici un dictionnaire)
	contenant les occurences de chaque caractère 
	output:
		(bwt,dic)
		où - bw: string
		- dic: dictionnaire
	"""
	bw = BWT(seq)
	dic = count_occ(seq)
	return (bw,dic)
	
	
	
	
	