# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:06:36 2020

@author: ADMIN
"""


def base_compare(seq1,seq2,i,j):
    if seq1[i] == seq2[j]:
        return 2
    else:
        return -1
    
# Needleman-Wunsch Alignment 
def NWalignment(seq1, seq2):
    
    # Construction de la matrice de similarité
    m = len(seq1)
    n = len(seq2)
    g = -3
    matrix = []
    for i in range(0, m):
        matrix.append([0 for j in range(0, n)])
    for sii in range(0, m):
        matrix[sii][0] = sii*g
    for sjj in range(0, n):
        matrix[0][sjj] = sjj*g
    for siii in range(1, m):
        for sjjj in range(1, n):
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + base_compare(seq1,seq2,siii, sjjj), matrix[siii][sjjj-1] + g)
    
    # Construction de l'alignement entre deux séquences
    sequ1 = [seq1[m-1]]
    sequ2 = [seq2[n-1]]
    
    while m > 1 and n > 1:
        if max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-2][n-2]:
            m -= 1
            n -= 1
            sequ1.append(seq1[m-1])
            sequ2.append(seq2[n-1])
        elif max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-1][n-2]:
            n -= 1
            sequ1.append('-')
            sequ2.append(seq2[n-1])
        else:
            m -= 1
            sequ1.append(seq1[m-1])
            sequ2.append('-')
    sequ1.reverse()
    sequ2.reverse()
    align_seq1 = ''.join(sequ1)
    align_seq2 = ''.join(sequ2)
    
    # Calcul du score de l'alignement
    align_score = 0
    for k in range(0, len(align_seq1)):
        if align_seq1[k] == align_seq2[k]:
            align_score += 1
    align_score = align_score/len(align_seq1)
    return align_seq1, align_seq2, align_score


sequence1 = 'GTGCCCCGGCGCCACGANGG'
sequence2 = 'GTGCCCGGCTGCAAGCANGG'
sequence3 = 'CGCCCCCTCGTGGCGCCNGG'

NWalignment(sequence1, sequence2)


