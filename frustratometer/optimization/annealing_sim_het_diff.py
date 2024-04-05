import random
import numpy as np
import frustratometer
import pandas as pd  # Import pandas for data manipulation
import time
from pathlib import Path
import argparse
import random
from math import log, factorial
from collections import Counter


def sequence_permutation(sequence):
    sequence = sequence.copy()
    res1, res2 = random.sample(range(len(sequence)), 2)
    sequence[res1], sequence[res2] = sequence[res2], sequence[res1]
    het_diff = 1
    return sequence, het_diff

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

# def heterogeneity(sequence):
#     N = len(sequence)
#     _, counts = np.unique(sequence, return_counts=True)
#     denominator = np.prod(np.array([np.math.factorial(count) for count in counts]))
#     het = np.math.factorial(N) / denominator
#     return np.log(het)

def heterogeneity_approximation(sequence):
    """
    Uses Stirling's approximation to calculate the heterogeneity of a sequence
    """
    N = len(sequence)
    _, counts = np.unique(sequence, return_counts=True)
    def stirling_log(n):
        return 0.5 * np.log(2 * np.pi * n) + n * np.log(n / np.e)
    
    log_n_factorial = stirling_log(N)
    log_denominator = sum(stirling_log(count) for count in counts)
    het = log_n_factorial - log_denominator
    return het

def montecarlo_steps(temperature, model, sequence, Ep=100, n_steps = 100):
    kb = 0.001987
    energy = model.native_energy(sequence) # Get energy and het every 1000 steps
    het = heterogeneity_approximation(sequence)
    for _ in range(n_steps):
        # new_sequence, energy_difference, het_difference = sequence_permutation(sequence) if random.random() > 0.5 else sequence_mutation(sequence)
        new_sequence, het_difference = sequence_permutation(sequence) if random.random() > 0.5 else sequence_mutation(sequence)
        new_energy = model.native_energy(new_sequence)
        # new_het = heterogeneity_approximation(new_sequence)
        energy_difference = new_energy - energy
        # het_difference = new_het - het
        exponent=(-energy_difference + Ep * het_difference) / (kb * temperature)
        acceptance_probability = np.exp(min(0, exponent))
        if random.random() < acceptance_probability:
            sequence = new_sequence
            energy = new_energy
            # het = new_het
    return sequence, energy, het

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--to", default="./", help="location of output file")
    args = parser.parse_args()

    # native_pdb = "tests/data/1r69.pdb"
    native_pdb = "./1r69.pdb"
    structure = frustratometer.Structure.full_pdb(native_pdb, "A")
    model = frustratometer.AWSEM(structure, distance_cutoff_contact=10, min_sequence_separation_contact=2)
    sequence = list("SISSRVKSKRIQLGLNQAELAQKVGTTQQSIEQLENGKTKRPRFLPELASALGVSVDWLLNGT")
    
    toPath = Path(args.to)
    toPath.mkdir(exist_ok=True)
    start_time = time.time()

    simulation_data = []
    for temp in range(800, 199, -1):
        sequence, energy, het = montecarlo_steps(temp, model, sequence, Ep=10, n_steps=1000)
        simulation_data.append({'Temperature': temp, 'Sequence': ''.join(sequence), 'Energy': energy, 'Heterogeneity': het})
        print(temp, ''.join(sequence), energy, het)
    simulation_df = pd.DataFrame(simulation_data)
    simulation_df.to_csv(toPath/"mcso_simulation_results.csv", index=False)

    time_taken = time.time() - start_time  # time_taken is in seconds
    timeFile = toPath / "time.dat"
    with open(timeFile, "w") as out:
        out.write(str(time_taken)+"\n")