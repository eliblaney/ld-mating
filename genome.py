import random

class Genome:

    name = "Unknown Genome"
    seq1 = ""
    seq2 = ""
    exons = []

    def __init__(self, name, seq, exons):
        self.name = name
        self.seq1 = seq
        self.seq2 = rc(seq)
        self.exons = exons

    def get_name(self):
        return self.name

    def get_sequences(self):
        return self.seq1, self.seq2

    def mutate(self):
        raise NotImplementedError("Genome.mutate() is not implemented yet.")

    def crossover(self, other):
        raise NotImplementedError("Genome.crossover() is not implemented yet.")

    def get_fitness(self):
        raise NotImplementedError("Genome.get_fitness() is not implemented yet.")

    def get_allele(self):
        raise NotImplementedError("Genome.get_allele() is not implemented yet.")

    def get_biallelic(self, allele):
        # Returns something like [0, 1]
        raise NotImplementedError("Genome.get_biallelic() is not implemented yet.")


default_probabilities = [
    ([0.27, 0.20, 0.23, 0.30], [0.99, 0.1]),     # Intron
    ([0.0005, 0.0001, 0.9993, 0.0001], [0, 1]), # Splice1
    ([0.0001, 0.0069, 0.0001, 0.9929], [0, 1]), # Splice2
    ([0.2, 0.3, 0.3, 0.2], [0.99, 0.1]),         # Exon
    ([0.0005, 0.0001, 0.9993, 0.0001], [0, 1]), # Splice1
    ([0.0001, 0.0069, 0.0001, 0.9929], [0, 1]), # Splice2
    ([0.27, 0.20, 0.23, 0.30], [0.99, 0.1])      # Intron
]

def markov_model(probabilities=default_probabilities):
    states = "ACGT*>"
    return [dict(zip(states, [*p[0], *p[1]])) for p in default_probabilities]

def rc(seq):
    """Find the reverse complement of a sequence.

    Keyword arguments:
    seq -- A DNA sequence
    """
    return seq[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

def random_genome(min_length=100000, name="Random Genome", markov=markov_model):
    seq = ""
    length = 0
    exons = []

    mm = markov()
    num_stages = len(mm)
    stage = 0

    while length < min_length:
        l = 0
        threshold = mm[stage]['*']
        while random.random() < threshold:
            l += 1
        if stage == 3: # Mark exon bounds
            exons.append((length, length+l))
        seq += ''.join(random.choices("ACGT", weights=list(mm[stage].values())[:4], k=l))
        length += l
        stage = (stage + 1) % num_stages
    return Genome(name, seq, exons)
