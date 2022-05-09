from genome import Genome
from population import Population
from mating import random_mating, monogamy, polygamy, polyandry

# Set up mating structures
mating_structures = [random_mating, monogamy, polygamy, polyandry]

# Run simulations for each mating structure, recording regular time points
for ms in mating_structures:
    pop = Population(ms.__name__, ms, population_size=100)

    generations = 100
    with open(pop.name + '.csv', 'w') as f:
        print("Simulating {} with {} generations and a population size of {}...".format(pop.name, generations, len(pop.get_population())))
        f.write("generation,Dmin,Dmax\n")
        for i in range(generations):
            if i % (generations / 5) == 0:
                print("{:.0f}% complete...".format(i / generations * 100))
            pop.run(1)

            Dmin, Dmax = pop.linkage_disequilibrium_avg_D()
            f.write("{},{},{}\n".format(pop.get_generation(), Dmin, Dmax))
            
            # Analyze LD between different structures, run ANOVA and possibly other tests
            # f.write("{};{};{}\n".format(pop.name, len(p), pop.get_generation()))
            # f.write("{}\n".format(pop.linkage_disequilibrium_r2().tolist()))
            
            # [f.write("{:.2f};{};{}\n".format(g.get_fitness(), g.exons1, g.exons2)) for g in p]