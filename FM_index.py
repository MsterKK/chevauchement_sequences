# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:34:16 2020

Certaines parties du code sont inspirées d'algorithmes trouvés sur Github:
- notamment https://gist.github.com/Puriney/6324227
"""


def suffix_array(seq) :
    """Fonction qui construit la table des suffixes associée à la chaîne de caractères seq
	--- 
	output: 
		la liste contenant les rangs des suffixes triés
	"""
    
	# on ajoute le caractère $ à la fin de la chaîne de caractères
    seq = seq + '$'
    
    # tri de la liste composée des suffixes ainsi que du rang du suffixe
    liste_suffixes = sorted([(seq[i:],i) for i in range(len(seq))])
    
    #Construction de la liste contenant les rangs de la liste des suffixes
    S = [] ; L = ''
    for suffixe in liste_suffixes :
        S.append(suffixe[1])
        L += suffixe[0][0]
    
    return L, S

def BWT(seq, liste_rank):
	'''Fonction qui réalise la transformation de Burrows-Wheeler à partir de la table des suffixes
	---
	output:
		bw: string
	'''
	bw = ''
	#ajout du caractère $ à la chaîne de caractère originale
	seq2 = seq +'$'
	#récuperation des rangs de la table des suffixes pour construire la BWT
	for c in liste_rank:
		bw += seq2[c-1]
	return bw

#def count_occ(seq):
#	"""Fonction qui, pour un caractère c, compte les occurences de tous les caractères lexicographiquement plus 
#	petits que c
#	--- 
#	output
#		dic: dictionnaire
#	"""
#	
#	dic = {}
#	seq = seq + "$"
#	alphabet = ["$","A","C","G","T"]
#	#initialisation du dictionnaire
#	for letter in alphabet:
#		dic[letter]=0
#	#comptage du nombre de lettre de chaque type
#	for letter in seq:
#		dic[letter] += 1
#	#inverse la liste dans l'ordre lexicographique
#	alphabetRev = alphabet[::-1]
#	nb_carac = len(seq)
#	#ajout du nb de caractères précédant chaque lettre en partant de la derniere lettre dans l'ordre lexicographique
#	for letter in alphabetRev:
#			dic[letter] = nb_carac - dic[letter]
#			nb_carac = dic[letter]
#	return dic


def count_occ(seq_sorted) :
    """Fonction qui, pour une chaine de caractère seq, après avoir été transformée par BWT / Suffix Array
    compte les occurences de tous les caractères lexicographiquement plus 
	petits que c
	--- 
	output
		dic: dictionnaire
	"""
    # Initialisation du dictionnaire et du compteur d'occurence
    table_C = {}
    compteur = 0
    
    # Parcours de la chaine de caractère
    for letter in seq_sorted :
        if letter not in table_C :
            table_C[letter] = compteur
        compteur += 1
    
    # Ajout du dernier élément
    table_C['#'] = compteur
    
    return table_C


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


def FMindex(seq) :
    """Fonction qui renvoie la transformation de Burrows-Wheeler ainsi que la table C[c] (ici un dictionnaire)
	contenant les occurences des caractère lexicographiquement plus petit que c, et la table des occurences (voir rapport
	pour détail)
	output:
		(bwt,dic)
		où - bw: string
		- liste_L: dictionnaire
		- tables_occurences
	"""
    liste_L, liste_rank = suffix_array(seq)
    bw = BWT(seq, liste_rank)
    table_C = count_occ(liste_L)
    table_occurences = count_table(bw)
    return bw, table_C, table_occurences


def Backward_count(seq, query) :
    
    # Initialisation de FM-index de la séquence de référence
    bw, table_C, table_occurrences = FMindex(seq)
    keys_C = list(table_C.keys())
    nxt_key = {keys_C[i] : keys_C[i + 1] for i in range(1, len(keys_C) - 1)}
    
    # Initialisation de l'algo
    query = query[::-1]
    sub_query = query[0]
    if sub_query not in keys_C :
        dict_bw = {}
    else :
        sp = table_C[sub_query]
        ep = table_C[nxt_key[sub_query]] - 1
        dict_bw = {sub_query : (sp, ep)}
        
        # Parcours du query
        index = 1
        while index <= len(query) - 1 and sp <= ep :
            
            sub_query = query[index]
            sp = table_C[sub_query] + table_occurrences[sub_query][sp - 1]
            ep = table_C[sub_query] + table_occurrences[sub_query][ep] - 1
            dict_bw[query[:index]] = (sp, ep)
            index += 1
            
    return dict_bw
        
import generation_sequences as gs
liste_adn = gs.gen_seq()
seq = liste_adn[0]
Backward_count('ACTTTAC', 'TT')
