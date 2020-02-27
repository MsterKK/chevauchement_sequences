# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:33:30 2020

@author: kevin
"""

from generation_sequences import gen_seq
import alignement_inex as ai

#génération des chaînes à traiter
liste_chaines = gen_seq(nb_seq = 1000, l_min = 5, l_max= 10)
tailles_chaines =[len(chaine) for chaine in liste_chaines]

#affichage des chaînes
for chaine in liste_chaines:
	print(chaine)
	
	#Paramètres pour l'alignement
#Nombre de mismatch maximal entre 2 séquences
z = 2


# Pré-traitement des chaînes
dic_pre_traitement = {}
nb_chaines = len(liste_chaines)
for num_chaine in range(1,nb_chaines+1):
	nom_chaine = 'S'+str(num_chaine)
	dic_pre_traitement[nom_chaine] = ai.pretraitement #liste_chaines[num_chaine-1]
dic_chevauchement = {}
#Parcours de toutes les combinaisons de chaînes distinctes possibles
for n1 in range(1,nb_chaines):
	for n2 in range(n1 + 1, nb_chaines+1):
		
		id_ch1 = 'S'+str(n1)
		id_ch2 = 'S' + str(n2)
		#La chaine la plus grande est désignée comme étant la reférence
		if tailles_chaines[n1-1] < tailles_chaines[n2-1]:
			key_ref = id_ch2
			seq_ref = liste_chaines[n2-1]
			key_query = id_ch1
			query = liste_chaines[n1-1]
			
		else: 
			key_ref = id_ch1
			seq_ref = liste_chaines[n1-1]
			key_query = id_ch2
			query = liste_chaines[n2-1]
		
		#récupération des données sur la séquence de ref
		C_ref, occ_ref, Rev_occ_ref = dic_pre_traitement[key_ref]
		
		#alignement entre les deux séquences
		I = ai.alignement_inexact(query, seq_ref, C_ref, occ_ref, Rev_occ_ref, z)
		clef = id_ch1+'-'+id_ch2
		if I != set():
			dic_chevauchement[clef] = I

print(dic_chevauchement)
		
		

#alignement_inexact(query, seq_ref, C_ref, occ_ref, Rev_occ_ref, z)
