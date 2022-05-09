from genome import Genome
from population import Population
from mating import random_mating, monogamy, polygamy, polyandry

# Set up mating structures
mating_structures = [random_mating, monogamy, polygamy, polyandry]

# Run simulations for each mating structure, recording regular time points
for ms in mating_structures:
    pop = Population(ms.__name__, ms, population_size=100)

    generations = 100
    with open(pop.name + '.txt', 'w') as f:
        print("Simulating", pop.name + "...")
        f.write("{};{};{}\n".format(pop.name, len(pop.get_population()), generations))
        for i in range(100):
            pop.run(1)

            f.write("{};{}\n".format(pop.get_generation(), pop.linkage_disequilibrium_r2().tolist()))
            
            # Analyze LD between different structures, run ANOVA and possibly other tests
            # f.write("{};{};{}\n".format(pop.name, len(p), pop.get_generation()))
            # f.write("{}\n".format(pop.linkage_disequilibrium_r2().tolist()))
            
            # [f.write("{:.2f};{};{}\n".format(g.get_fitness(), g.exons1, g.exons2)) for g in p]