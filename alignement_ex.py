# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 21:36:21 2020

@author: ADMIN
"""

import generation_sequences as gs
import FM_index as fm

liste_adn = gs.gen_seq()
seq = liste_adn[0]

seq="ATCATCG"
query = "ATC"


class Alignement_ex() :

    def __init__(self, seq_ref, query) :
        self.seq_ref = seq_ref
        self.query = query
    
    
def Backward_count(seq, query) :
    """
    La fonction qui sert à faire l'alignement exact, donne une dictionnaire des intervalles SA des éléments de query (par ex ATC : C, CT, CTA) 
    """
    
    # Initialisation de FM-index de la séquence de référence
    liste_L, liste_rank, liste_suffixes = fm.suffix_array(seq)
    bw, table_C, table_occurrences = fm.FMindex(seq)
    
    # Construction de la table des occurrences au final
    for key in table_occurrences :
        table_occurrences[key].pop(0)
    
    # Dictionnaire de la clé suivante pour la table C
    keys_C = list(table_C.keys())
    nxt_key = {keys_C[i] : keys_C[i + 1] for i in range(1, len(keys_C) - 1)}
    
    # Reverser le query
    query = query[::-1]
    sub_query = query[0]
    
    # Initialisation du dictionnaire
    if sub_query not in keys_C :
        dict_bw = {}
    else :
        sp = table_C[sub_query]
        ep = table_C[nxt_key[sub_query]] - 1

        dict_bw = {sub_query : (sp, ep)}
    
        # Parcours du query
        index = 1
        while index <= len(query) - 1 and sp <= ep :
            
            # Calcul de l'intervalle du sub_query
            sub_query = query[index]
            sp = table_C[sub_query] + table_occurrences[sub_query][sp - 1]
            ep = table_C[sub_query] + table_occurrences[sub_query][ep] - 1
            
            # Enregistrement dans le dictionnaire
            dict_bw[query[:index + 1]] = (sp, ep)
            
            # Passage à l'élément suivant
            index += 1
            
    return dict_bw



Backward_count(seq, query)
#Backward_count(seq, 'ACTTTAC')