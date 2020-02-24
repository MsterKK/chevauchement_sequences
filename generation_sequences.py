# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:25:37 2020

@author: kevin
"""

import random as rd


def gen_seq(nb_seq = 10, l_min = 30, l_max= 50):
	"""Génère un nombre nb_seq de chaines de caractères aléatoires contenant les caractères {A,C,G,T} 
	d'une longueur allant de l_min à l_max. Ces chaînes de caractères sont retournées dans une liste.
	"""
	
	liste_chaines = []
	
	for k in range(nb_seq):
		letters = 'ACGT'
		length = rd.randint(l_min, l_max)
		rd_seq = ''.join((rd.choice(letters) for i in range(length)))
		liste_chaines.append(rd_seq)
	
	return liste_chaines