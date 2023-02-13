# Adam Langowski 147896

import random
import numpy as np


def fitness(used_oligonucleotides, n, k):
    """
    Calculates the fitness score of a solution based on the number of used oligonucleotides in a DNA spectrum.

    :param used_oligonucleotides: The number of oligonucleotides that have been used in the DNA spectrum
    :param n: The length of the DNA sequence
    :param k: The length of each oligonucleotide in the DNA spectrum
    :return: The fitness score, represented as the proportion of used oligonucleotides to the expected number of oligonucleotides in the DNA spectrum
    """
    oligonucleotides_expected_count = n - k + 1
    score = used_oligonucleotides / oligonucleotides_expected_count
    return score


def single_crossover(parent_1, parent_2, idx=0) -> object:
    """
    Performs single-point crossover on two parent solutions to produce two offspring solutions.

    :param parent_1: One parent solution
    :param parent_2: Another parent solution
    :param idx: The index of the single point where the crossover should occur. If idx is not specified, it is set
     to a default value of 0.25 * len(parent_1).
    :return: Two offspring solutions created by combining the first idx elements of parent_1 and the last
     len(parent_1) - idx elements of parent_2, and vice versa.
    """
    if idx == 0:
        idx = int(0.25 * len(parent_1))
    child_1 = np.append(parent_1[:idx], parent_2[idx:])
    child_2 = np.append(parent_2[:idx], parent_1[idx:])
    return child_1, child_2


def crossover(parent_1, parent_2):
    """
    Performs multiple single-point crossovers on two parent solutions to produce two offspring solutions.

    :param parent_1: One parent solution
    :param parent_2: Another parent solution
    :return: Two offspring solutions created by performing multiple single-point crossovers on the parent solutions.
     The specific points at which the crossovers occur are determined by the values in the x array.
    """
    child_1 = child_2 = ''
    x = np.array([2, (0.9 * len(parent_1))])
    for i in x:
        child_1, child_2 = single_crossover(parent_1, parent_2, int(i))
    return child_1, child_2


def tournament(candidates, scores):
    """
    Select the best candidate from a group of candidates based on their fitness scores.

    :param candidates: List of candidate solutions
    :param scores: List of fitness scores for the candidates
    :return: Tuple of the winning candidate and its corresponding fitness score
    """
    max_fitness = max(scores)
    max_index = scores.index(max_fitness)
    winner = candidates[max_index]
    return winner, max_fitness


def mutation(spectrum, mutations_percentage=10):
    """
    Mutates the input `spectrum` by a given percentage of mutations.

    :param spectrum: A list of sequences that represent the input spectrum.
    :param mutations_percentage: The percentage of mutations to perform on the input `spectrum` (default is 10).
    :return: None
    The input `spectrum` is modified in place.
    """
    mutations_count = int((len(spectrum) * mutations_percentage) // 100)
    sequences_mutations_count = int((len(spectrum) * 5) // 100)
    options = ['insertion', 'deletion', 'substitution']

    for _ in range(mutations_count):
        random_seq = random.randrange(len(spectrum))
        chosen_seq = spectrum[random_seq]
        for _ in range(sequences_mutations_count):
            random_oligonucleotide = random.randrange(len(chosen_seq) - 1)
            type_of_mutation = random.choice(options)

            chosen_oligonucleotide = list(spectrum[random_seq][random_oligonucleotide].strip(" "))
            '''checks that result in not an empty string'''
            if len(chosen_oligonucleotide) <= 2:
                continue

            random_nuc = random.randrange(0, len(chosen_oligonucleotide) - 1)
            if type_of_mutation == "insertion":
                add = random.choice(["A", "C", "G", "T"])
                chosen_oligonucleotide.append(add)
            elif type_of_mutation == "deletion":
                del chosen_oligonucleotide[random_nuc]
            elif type_of_mutation == "substitution":
                sub = random.choice(["A", "C", "G", "T"])
                chosen_oligonucleotide[random_nuc] = sub

            spectrum[random_seq][random_oligonucleotide] = ''.join(chosen_oligonucleotide)
