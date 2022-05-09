import math
import random
import genome
import allel

class Population:

    name = "Unknown Population"
    ms = None
    N = 100
    generation_num = 0
    population = []
    num_worst = 20
    total = 1

    def __init__(self, name, mating_structure, population_size=100, num_worst=5):
        self.name = name
        self.ms = mating_structure
        self.N = population_size
        self.num_worst = num_worst 

        # Create initial population at random
        for i in range(self.N):
            self.population.append(genome.random_genome())

    def run(self, generations=50):
        for i in range(generations):
            self.new_generation()

    def get_population(self):
        # sort models in population by fitness ascending
        self.population.sort(key=lambda g: g.get_fitness())
        return self.population

    def new_generation(self):
        self.generation_num = self.generation_num + 1

        self.get_population()
        
        # replace the worst models with children based on fitness
        self.population = self.population[self.num_worst:]

        self.ms(self)

        # do mutations on remaining models
        for g in self.population:
            g.mutate()

    def linkage_disequilibrium_r2(self):
        """Estimates the linkage disequilibrium parameter r^2 for each individual in the population."""
        g = allel.GenotypeArray([g.get_biallelic() for g in self.population], dtype='i1')
        gn = g.to_n_alt(fill=-1)
        r = allel.rogers_huff_r(gn)
        return r ** 2

    def linkage_disequilibrium_avg_r2(self):
        """Estimates the average linkage disequilibrium parameter r^2 for the entire population."""
        r2 = self.linkage_disequilibrium_r2().tolist()
        return sum(r2) / len(r2)

    def allele_frequencies(self, pop=None):
        """Finds the allelic frequencies for each variant in the population."""
        if pop == None:
            pop = self.population
        haplotypes = genome.haplotypes
        counts = []
        for i in range(len(haplotypes)):
            variants = haplotypes[i]
            for j in range(len(variants)):
                c = 0
                for g in [g.get_biallelic() for g in pop]:
                    for k in g[i]:
                        if k == j:
                            c += 1
                counts.append(c)
        total = sum(counts)
        freqs = [c / total for c in counts]
        return freqs

    def linkage_disequilibrium_D(self):
        """Estimates the linkage disequilibrium parameter D for each individual in the population."""
        r2 = self.linkage_disequilibrium_r2()
        D = []
        for i in range(self.N):
            D.append(r2[i] * math.prod(self.allele_frequencies(pop=self.population[i])))
        return D

    def linkage_disequilibrium_avg_D(self):
        """Estimates the average linkage disequilibrium parameter D for the entire population."""
        return self.linkage_disequilibrium_avg_r2() * math.prod(self.allele_frequencies())

    def get_generation(self):
        return self.generation_num
