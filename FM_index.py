# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:34:16 2020
"""

global alphabet
alphabet = ["$","A","C","G", "N", "T"] 

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
    S = [] ; seq_sorted = ''
    for suffixe in liste_suffixes :
        S.append(suffixe[1])
        seq_sorted += suffixe[0][0]
    
    return seq_sorted, S, liste_suffixes

def BWT(seq, liste_rank):
    '''Fonction qui réalise la transformation de Burrows-Wheeler à partir de la table des suffixes
    ---
    output:
        bw: string
    '''
    bw = ''
    #ajout du caractère $ à la chaîne de caractère originale
    seq2 = seq + '$'
    #récuperation des rangs de la table des suffixes pour construire la BWT
    for c in liste_rank:
        bw += seq2[c-1]
    return bw


def calcul_C(seq_sorted) :
    """Fonction qui, pour une chaine de caractère seq, après avoir été transformée par BWT / Suffix Array
    compte les occurences de tous les caractères lexicographiquement plus 
    petits que c
    --- 
    output
        dic: dictionnaire
    """
    alphabet = ["$","A","C","G", "N", "T","#"]                             
    
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

    #Rajout des caractères non présents dans seq_sorted (s'il en existe)
    for k in range(len(alphabet)-2,0,-1):
        if alphabet[k] not in table_C:
            table_C[alphabet[k]] = table_C[alphabet[k+1]]                                                                    
    
    return table_C


def count_table(seq):
    """Fonction qui renvoie la table des occurences d'une sequence, çad un dictionnaire qui, pour chaque lettre de l'alphabet,
    contient la ligne du tableau des occurences (voir rapport pour le détail)
    -
    output
        dic_table: dictionnaire contenant des listes pour chaque clef
    """
    dic_table = {}
    global alphabet
    taille_seq = len(seq)
    #initialisation du dictionnaire
    for letter in alphabet:
        dic_table[letter] = [0]*(taille_seq+1)
    #parcours de la chaine de caractères pour trouver les occurences de chaque caractère
    for k in range(1,taille_seq+1):
        c = seq[k-1]
        dic_table[c][k] = 1
    #concatenation des occurences trouvées afin de former le tableau
    for letter in alphabet:
        for k in range(1,taille_seq+1):
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
    liste_L, liste_rank, liste_suffixes = suffix_array(seq)
    bw = BWT(seq, liste_rank)
    table_C = calcul_C(liste_L)
    table_occurences = count_table(bw)
    return bw, table_C, table_occurences


    


