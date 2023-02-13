# Adam Langowski 147896

import instance_generator
import genetic_algorithm_functions
import time
import itertools
import random

'''Parameters with possible modifications'''
given_length_n = 100
given_k = 8
permutations_count = 100000
population_count = 300
max_generations_count = 50
errors_percentage = 5
positive_errors_percentage = 0.5
negative_errors_percentage = 0.5

'''List storing fitnesses of sequences'''
mean_fitness = []

'''Lists for further applications'''
permutation_fitness = list()
permutation_list = list()
permutation_count = list()
next_gen = list()
next_gen_fitness = list()


def average(given_list):
    """
    This function calculates the average of the elements in the given list.

    :param given_list: A list of numbers whose average is to be calculated.
    :return: The average of the elements in the given list.
    """
    return sum(given_list) / len(given_list)


def group(iterable, n):
    """
    Divides the input iterable into groups of `n` elements.

    :param iterable: The iterable to be divided into groups.
    :param n: The number of elements in each group.
    :return: A list of tuples, where each tuple represents a group of `n` elements from the input iterable.
    """
    return [iterable[i:i + n] for _ in range(0, len(iterable), n)]


def max_overlapping_suffix_length(a, b):
    """
    This function takes in strings `a` and `b`, and returns the length of the maximum common suffix between `a` and `b`.

    :param a: A string representing the first input.
    :param b: A string representing the second input.
    :return: An integer representing the length of the maximum common suffix between `a` and `b`.
    """
    return max(j for j in range(len(b) + 1) if a.endswith(b[:j]))


def overlap(permutations):
    """
    Given a list of strings, `permutations`, the function finds the concatenation of strings in the list
    such that the overlapping suffixes between two consecutive strings are maximized. The function returns
    the concatenated string and the count of strings used in the concatenation.

    :param permutations: list of strings
    :return: tuple of (concatenated string, count of strings used)
    """
    _count = 1
    _result = permutations[0]
    for pm in permutations[1:]:
        overlap_length = max_overlapping_suffix_length(_result, pm)
        if overlap_length != 0:
            _count += 1
            _result += pm[overlap_length:]
    return _result, _count


print("This is program based on genetic algorithm reconstructing DNA from oligonucleotides.")
print("If you want to modify any parameters you can change them in lines 9-16 in main.py file")

'''Starting calculating execution time'''
start = time.time()

print("Starting algorithm...")
print("Creating first random generation...")

spectrum = instance_generator.generate(given_length_n, given_k, errors_percentage, positive_errors_percentage,
                                       negative_errors_percentage)

print('Shuffling spectrum to make it non-trivial problem.')
random.shuffle(spectrum)

print('Generating a sequence of permutations_count permutations of the input spectrum list.')
for permutation in itertools.islice(itertools.permutations(spectrum), 0, permutations_count):
    permutation_list.append(permutation)
    result, count = overlap(permutation)
    permutation_count.append(count)

print('Calculating fitnesses.')
for org in permutation_count:
    fitness = genetic_algorithm_functions.fitness(org, given_length_n, given_k)
    permutation_fitness.append(fitness)

mean_fitness.append(average(permutation_fitness))  
print(f'First population fitness: {mean_fitness[0]}.')

print(f'Creating new sequences - solutions...')
for generations in range(max_generations_count):
    print(f"Generation number: [{generations + 1}] ")
    population = 0

    '''Tournament selection'''
    while len(next_gen) < population_count:
        selected_list = list()
        selected_fitness = list()
        for i in range(8):
            x = random.randrange(len(permutation_list))
            selected_list.append(permutation_list[x])
            selected_fitness.append(permutation_fitness[x])
        winner, winner_score = genetic_algorithm_functions.tournament(selected_list, selected_fitness)
        next_gen.append(list(winner))
        next_gen_fitness.append(winner_score)
        population += 1

    '''Crossing-over'''
    for one, two in group(next_gen, 2):
        options = ['single', 'multi', 'none']
        choice = random.choices(options, weights=[0.6, 0.3, 0.1])[0]
        if choice == 'single':
            child1, child2 = genetic_algorithm_functions.single_crossover(one, two)
            child1 = child1.tolist()
            child2 = child2.tolist()
            next_gen.append(child1)
            result, score = overlap(child1)
            next_gen_fitness.append(genetic_algorithm_functions.fitness(score, given_length_n, given_k))
            next_gen.append(child2)
            _result, score = overlap(child2)
            next_gen_fitness.append(genetic_algorithm_functions.fitness(score, given_length_n, given_k))

        if choice == 'multi':
            child1, child2 = genetic_algorithm_functions.crossover(one, two)
            child1 = child1.tolist()
            child2 = child2.tolist()
            next_gen.append(child1)
            result, score = overlap(child1)
            next_gen_fitness.append(genetic_algorithm_functions.fitness(score, given_length_n, given_k))
            next_gen.append(child2)
            _result, score = overlap(child2)
            next_gen_fitness.append(genetic_algorithm_functions.fitness(score, given_length_n, given_k))

        if choice == 'none':
            pass
    mean_fitness.append(average(next_gen_fitness))
    print(average(next_gen_fitness))

    '''Mutating'''
    genetic_algorithm_functions.mutation(next_gen)
    next_gen_fitness.clear()

    '''Evaluating the fitness of a set of sequences'''
    for org in next_gen:
        res, score = overlap(org)
        next_gen_fitness.append(genetic_algorithm_functions.fitness(score, given_length_n, given_k))
    mean_fitness.append(average(next_gen_fitness))

    '''Ordering lists of new generations'''
    tmp_lst = []
    for seq, score in zip(next_gen, next_gen_fitness):
        tmp_lst.append([seq, score])
    tmp_lst.sort(key=lambda row: row[1], reverse=True)
    next_gen.clear()
    next_gen_fitness.clear()
    permutation_list.clear()
    permutation_fitness.clear()
    while len(permutation_list) < population_count:
        queue = 0
        permutation_list.append(tmp_lst[0][0])
        permutation_fitness.append(tmp_lst[0][1])
        queue += 1

end = time.time()

'''Printing results'''
print(f'\nFirst population fitness: {mean_fitness[0]}')
print(f'Best fitness score after {max_generations_count} generations: {max(mean_fitness)}')
print(f'Fitness upgraded by {max(mean_fitness) - mean_fitness[0]}')
print("\nAlgorithm time: %s seconds" % (end - start))
