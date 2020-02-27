# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 21:36:21 2020

@author: ADMIN
"""

import generation_sequences as gs
import FM_index as fm

liste_adn = gs.gen_seq()
seq = liste_adn[0]

seq="ATCGATCG"
query = "ATCGTCG"


class Alignement_ex() :

    def __init__(self, seq_ref, query, index_ini = 0) :
    
        self.seq_ref = seq_ref
        self.query = query
        self.index_ini = index_ini
        self.index = index_ini
    
        # Prétraitement des séquences
        self.pretraitement()
        
    def pretraitement(self) :
        
        # Préparation des éléments par BWT / Suffix Array
        self.liste_L, self.liste_rank, self.liste_suffixes = fm.suffix_array(self.seq_ref)
        self.bw = fm.BWT(self.seq_ref, self.liste_rank)
        self.table_C = fm.calcul_C(self.liste_L)
        self.table_occurrences = fm.count_table(self.bw)
        
        # Construction de la table finale des occurrences (rank)
        for key in self.table_occurrences :
            self.table_occurrences[key].pop(0)
    
    def Backward_count(self) :
        """
        La fonction qui sert à faire l'alignement exact, donne une dictionnaire des intervalles SA des éléments de query (par ex ATC : C, CT, CTA) 
        """
        
        # Dictionnaire de la clé suivante pour la table C
        keys_C = list(self.table_C.keys())
        self.nxt_key = {keys_C[i] : keys_C[i + 1] for i in range(1, len(keys_C) - 1)}
        
        # Reverser le query
        query_rev = self.query[::-1]
        sub_query = query_rev[self.index]
        
        # Initialisation du dictionnaire
        if sub_query not in keys_C :
            self.dict_SA = {}
        else :
            sp = self.table_C[sub_query]
            ep = self.table_C[self.nxt_key[sub_query]] - 1
    
            self.dict_SA = {(sub_query, self.index_ini, self.index) : (sp, ep)}
        
            # Parcours du query
            self.index = self.index + 1
            while self.index <= len(query_rev) - 1 and sp <= ep :
                
                # Calcul de l'intervalle du sub_query
                sub_query = query_rev[self.index]
                sp = self.table_C[sub_query] + self.table_occurrences[sub_query][sp - 1]
                ep = self.table_C[sub_query] + self.table_occurrences[sub_query][ep] - 1

                # Enregistrement du sub_query précédent dans le dictionnaire
                self.dict_SA[(query_rev[self.index_ini:self.index + 1], self.index_ini, self.index)] = (sp, ep)
                
                # Passage à l'élément suivant
                self.index += 1
            
            # Récursivité
            if self.index <= len(query_rev) - 1 :
                alig_nxt = Alignement_ex(self.seq_ref, self.query, self.index - 1)
                add_dict = alig_nxt.Backward_count()
                self.dict_SA.update(add_dict)
        
        return self.dict_SA

    def position_alignement(self) :
        
        # Retourner le dictionnaire des intervalles SA dict_SA
        self.Backward_count()
        self.dict_SA = {(k[0][::-1], len(self.query) - k[2] - 1, len(self.query) - k[1] - 1) : v for k, v in self.dict_SA.items() if v[0] <= v[1]}
        
        # Dictionnaire du positionnement de l'alignement
        self.pos_align = {}
        for key in self.dict_SA :
            pos_SA = self.dict_SA[key]
            self.pos_align[key] = (self.liste_suffixes[pos_SA[0]][1], self.liste_suffixes[pos_SA[1]][1])
            
    def find_subalig(self, pos_align) : # essaie à maximiser l'alignement
        
        # Clémaximale
        max_align = 0
        cle = ()
        for k, v in pos_align.items() :
            if len(k[0]) >= max_align :
                max_align = len(k[0])
                cle, value = k, v
                
        return cle, value
     
    def alignement(self) :
        # Réalisation BWA exacte
        self.index_ini = self.index = 0
        self.position_alignement()
        
        seq_ref_ch = ''
        query_ch = ''
        
        pos_align = self.pos_align
        
        while len(pos_align) > 0 :
            key, val = self.find_subalig(pos_align)
            
            for pos_ref in val :
                dictx = 
                y = 
            
            
            # Nouveau dictionnaire du positionnement d'alignement
            # Enlever les sub_query de même union
            pos_align = {k : v for k, v in pos_align.items() if k[2] != key[2]}
            
            for pos_ref in val :
                seq_ref_ch.join(self.seq_ref[:pos_ref] + key[1]*'-' + key[0])
                query_ch.join(pos*'-' + )

            
        return seq_ref_ch, query_ch
        #self.find_subalig(pos_align)
        
        
    
al_in = Alignement_ex(seq, query)
al_in.alignement()
al_in.pos_align

al_in.position_alignement()
al_in.pos_align
al_in.find_subalig(al_in.pos_align)

al_in.dict_SA
liste_suffixes = al_in.liste_suffixes
test_occ = al_in.table_occurrences
liste_suffixes = al_in.liste_suffixes

test = Alignement_ex('ATCGATCG', 'GTGAT')
test.position_alignement()
test.pos_align

#Backward_count(seq, 'ACTTTAC')