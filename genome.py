import random

class Genome:

    name = "Unknown Genome"
    seq1 = ""
    seq2 = ""

    def __init__(self, name, seq):
        self.name = name
        self.seq1 = seq
        self.seq2 = seq

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


def random_genome():
    raise NotImplementedError("Genome.random_genome() is not implemented yet.")
