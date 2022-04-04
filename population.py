import random
import genome
import allel

class Population:

    name = "Unknown Population"
    ms = None
    N = 20
    max_generations = 50
    generation_num = 1
    population = []
    num_elite = 5
    num_worst = 5
    total = 1

    def __init__(self, name, mating_structure, population_size=20, generations=50, num_elite=5, num_worst=5):
        self.name = name
        self.ms = mating_structure
        self.N = population_size
        self.max_generations = generations
        self.num_elite = num_elite 
        self.num_worst = num_worst 
        self.total = self.N * self.max_generations

        # Create initial population at random
        for i in range(self.N):
            self.population.append(genome.random_genome())

    def new_generation(self):
        self.generation_num = self.generation_num + 1
        
        # sort models in population by fitness ascending
        self.population.sort(key=self.get_fitness, reverse=True)

        # replace the worst models with children based on fitness
        survivors = self.population[self.num_worst:]
        fitnesses = [s.get_fitness() for s in survivors]

        parents = random.choices(
                survivors,
                weights=fitnesses,
                k=2*self.num_worst
                )

        for i in range(self.num_worst):
            p1 = parents[i*2]
            p2 = parents[i*2+1]
            # Perform crossover mutations
            # TODO: Decide if crossovers should be weighted?
            # w1 = self.get_fitness(p1)
            # w2 = self.get_fitness(p2)
            p1.crossover(p2)

        # do mutations on remaining non-elite models
        for i in range(self.N - self.num_elite - self.num_worst):
            j = i + self.num_worst
            self.population[j].mutate()

    def linkage_disequilibrium_r2(self, allele):
        """Estimates the linkage disequilibrium parameter r^2 for each individual in the population for a given allele."""
        g = allel.GenotypeArray([genome.get_biallelic(allele) for genome in self.population], dtype='i1')
        gn = g.to_n_alt(fill=-1)
        r = allel.rogers_huff_r(gn)
        return r ** 2

    def get_generation(self):
        return self.generation_num
