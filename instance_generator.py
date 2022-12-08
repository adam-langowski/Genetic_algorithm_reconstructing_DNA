import random


def generate_instance(dna_length, oligo_length, error=0.01):
    # tworzenie sekwencji DNA
    alphabet = "ATCG"
    dna = ''.join([random.choice(alphabet) for _ in range(dna_length)])

    # tworzenie listy oligonukleotydów dla sekwencji DNA
    dna_oligos = []
    for i in range(len(dna) - oligo_length + 1):
        oligo = dna[i:i + oligo_length]
        dna_oligos.append(oligo)

    # losowe dodanie błędów (pozytywnych i negatywnych) do spektrum
    error_count = round(len(dna_oligos) * error)
    error_oligos = random.sample(dna_oligos, error_count)
    for i in range(error_count):
        if random.random() < 0.5:
            # dodanie błędu pozytywnego
            error_oligos[i] = ''.join([random.choice(alphabet) for _ in range(oligo_length)])
        else:
            # dodanie błędu negatywnego
            error_oligos[i] = random.choice(dna_oligos)

    # tworzenie spektrum z błędami
    spectrum = dna_oligos + error_oligos

    # losowe wybranie początkowego oligonukleotydu
    start_oligo = random.choice(spectrum)

    # zwracanie d
