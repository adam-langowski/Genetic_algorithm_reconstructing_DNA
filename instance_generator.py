import random


def generate(length, k, p, pos_errors, neg_errors):
    """
    Generates a DNA spectrum from a random DNA sequence with the given parameters.

    :param length: int, the length of the DNA sequence
    :param k: int, the length of the substrings
    :param p: int, the percentage of errors
    :param pos_errors: float, the percentage of positive errors
    :param neg_errors: float, the percentage of negative errors
    :return: list of str, the generated DNA spectrum
    """
    spectrum = list()
    sequence = ''.join(random.choice('ACGT') for _ in range(length))

    for substring in splice(sequence, k):
        spectrum.append(substring)
    mutate(spectrum, k, p, pos_errors, neg_errors)
    return spectrum


def generate_oligonucleotide(k):
    """
    Generates a random oligonucleotide

    :param k: int, oligonucleotide length
    :return: str, oligonucleotide
    """
    oligonucleotide = ''.join(random.choice('ACGT') for _ in range(k))
    return oligonucleotide


def splice(sequence, k):
    """
    Splices a given DNA sequence into overlapping substrings of length k.

    :param sequence: str, The DNA sequence to splice
    :param k: int, The length of the substrings
    :return: list of str, The spliced substrings
    """
    for i in range((len(sequence) - (k-1))):
        yield sequence[i:k+i]


def mutate(spectrum, k, p, pos_errors, neg_errors):
    """
    Introduces mutations into the given DNA spectrum.

    :param spectrum: A list of DNA substrings (oligonucleotides)
    :param k: Length of the oligonucleotides in the spectrum
    :param p: The overall percentage of mutations to introduce
    :param pos_errors: The percentage of positive mutations (additions)
    :param neg_errors: The percentage of negative mutations (deletions)
    :return: None
    The input `spectrum` is modified in place.
    """
    errors = ['positive', 'negative']
    mutation_count = (len(spectrum) * p) // 100
    for mutation in range(mutation_count):
        mutation_type = random.choices(errors, weights=(pos_errors, neg_errors))
        if ''.join(mutation_type) == 'positive':
            added = False
            while added is not True:
                oligonucleotide = generate_oligonucleotide(k)
                if generate_oligonucleotide(k) not in spectrum:
                    spectrum.append(oligonucleotide)
                    added = True
        elif ''.join(mutation_type) == 'negative':
            x = random.choice(spectrum)
            spectrum.remove(x)


def oligonucleotide_count(spectrum):
    """
    Counts the frequency of each oligonucleotide in a DNA spectrum.

    :param spectrum: A list of DNA substrings (oligonucleotides)
    :return: A dictionary where each key is an oligonucleotide in the spectrum and each value is the number of times
     that oligonucleotide occurs in the spectrum
    """
    counts_for_oligonucleotides = {}
    for oligonucleotide in spectrum:
        if oligonucleotide not in counts_for_oligonucleotides:
            counts_for_oligonucleotides[oligonucleotide] = 0
        counts_for_oligonucleotides[oligonucleotide] += 1
    return counts_for_oligonucleotides
