# Adam Langowski 147896

import random


def fitness(dna, spectrum, error=0.01):
    # tworzenie listy oligonukleotydów dla sekwencji DNA
    dna_oligos = []
    for i in range(len(dna) - 3):
        oligo = dna[i:i + 4]
        dna_oligos.append(oligo)

    # liczenie liczby zgodnych oligonukleotydów
    matching_oligos = 0
    for oligo in spectrum:
        if oligo in dna_oligos:
            matching_oligos += 1

    # uwzględnienie błędu (pozytywnego lub negatywnego)
    error_value = len(spectrum) * error
    matching_oligos += error_value

    # obliczenie "fitness" jako procent zgodnych oligonukleotydów
    fitness = matching_oligos / len(spectrum)
    return fitness


def crossover(dna1, dna2):
    # losowanie punktu cięcia
    cut_point = random.randint(0, len(dna1))

    # sklejanie części sekwencji rodziców w nową sekwencję
    child = dna1[:cut_point] + dna2[cut_point:]
    return child


def mutate(dna):
    # mutacja sekwencji DNA
    pass


def reconstruct_dna(spectrum, dna_length):
    # tworzenie początkowej populacji sekwencji DNA
    population = [random.randint(0, 4 ** dna_length - 1) for _ in range(100)]

    # pętla główna algorytmu
    while True:
        # obliczenie "fitness" dla każdej sekwencji w populacji
        fitness_values = [fitness(dna, spectrum) for dna in population]

        # wybieranie najlepszych sekwencji (rodziców)
        parents = [population[i] for i in range(len(population)) if fitness_values[i] > 0.99]

        # jeśli mamy już odpowiednią liczbę rodziców, kończymy działanie algorytmu
        if len(parents) >= dna_length:
            break

        # tworzenie potomków przez krzyżowanie i mutację rodziców
        children = []
        while len(children) < len(population):
            dna1 = random.choice(parents)
            dna2 = random.choice(parents)
            child = crossover(dna1, dna2)
            child = mutate(child)
            children.append(child)

        # nowa populacja zastępuje poprzednią
        population = children

    # zwracanie najlepszych sekwencji (odpowiadających spektrowi)
    return parents[:dna_length]
