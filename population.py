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
        """Estimates the linkage disequilibrium parameter r^2 for each individual in the population"""
        g = allel.GenotypeArray([g.get_biallelic() for g in self.population], dtype='i1')
        gn = g.to_n_alt(fill=-1)
        r = allel.rogers_huff_r(gn)
        return r ** 2

    def get_generation(self):
        return self.generation_num
