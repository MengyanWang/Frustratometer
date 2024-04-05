import random
from math import log, factorial
from collections import Counter
import numpy as np
import math

def sequence_mutation(sequence):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    sequence = sequence.copy()
    res = random.randint(0, len(sequence) - 1)
    aa_old = sequence[res]
    aa_new = random.choice(amino_acids)
    sequence[res] = aa_new
    # count for all aa in new seq
    amino_acid_counts = Counter(sequence)
    
    n_old = amino_acid_counts.get(aa_old, 0) + 1  
    n_new = amino_acid_counts[aa_new] - 1 
    
    het_diff = log(factorial(n_old) * factorial(n_new) / (factorial(n_old - 1) * factorial(n_new + 1)))
    
    return sequence, het_diff

def heterogeneity(sequence):
    N = len(sequence)
    _, counts = np.unique(sequence, return_counts=True)
    denominator = np.prod(np.array([np.math.factorial(count) for count in counts]))
    # denominator = np.prod(np.array([math.factorial(count) for count in counts]))
    het = np.math.factorial(N) / denominator
    # het = math.factorial(N) / denominator
    return np.log(het)
    # return het

# def np(seq):
#     N = len(seq)
#     amino_acid_counts = Counter(seq)
#     np_numerator = factorial(N)
#     np_denominator = 1
#     for count in amino_acid_counts.values():
#         np_denominator *= factorial(count)
#     return np_numerator / np_denominator

def het_diff_simplify(mutated_sequence):
    aa_old = "I"
    aa_new = "A"
    amino_acid_counts = Counter(mutated_sequence)
    n_old = amino_acid_counts.get(aa_old, 0) + 1  
    n_new = amino_acid_counts[aa_new] - 1 
    het_diff = log(factorial(n_old) * factorial(n_new) / (factorial(n_old - 1) * factorial(n_new + 1)))

    return het_diff

def het_diff_simplify2(mutated_sequence):
    aa_old = "I"
    aa_new = "A"
    amino_acid_counts = Counter(mutated_sequence)
    n_old = amino_acid_counts.get(aa_old, 0) + 1  
    n_new = amino_acid_counts[aa_new] - 1 
    # het_diff = log(factorial(n_old) * factorial(n_new) / (factorial(n_old - 1) * factorial(n_new + 1)))
    het_diff = log(n_old /(n_new+1) )

    return het_diff

sequence = list("SISSRVKSKRIQLGLNQAELAQKVGTTQQ")
mutated_sequence = list("SASSRVKSKRIQLGLNQAELAQKVGTTQQ")#mut the first aa from S to A
# het_diff =  log(heterogeneity(mutated_sequence)) -  log(heterogeneity(sequence))
het_diff =  heterogeneity(mutated_sequence) - heterogeneity(sequence)
# het_diff = log(np(mutated_sequence)/np(sequence))
print('Before Hetergeneity Simplification', het_diff)
# mutated_sequence, het_diff = sequence_mutation(sequence)
het_diff_simplify = het_diff_simplify(mutated_sequence=mutated_sequence)
print('After Hetergeneity Simplification', het_diff_simplify)
assert np.allclose(het_diff,het_diff_simplify)


    
