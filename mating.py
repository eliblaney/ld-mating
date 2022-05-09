import math
import random
import genome

def random_mating(population):
    fitnesses = lambda gs: [g.get_fitness() for g in gs]
    pop = population.get_population()
    num = population.N
    fitnesses = fitnesses(pop)
    parents = random.choices(
            pop,
            weights=fitnesses,
            k=2*num
            )

    pop.clear()
    for i in range(num):
        p1 = parents[i*2]
        p2 = parents[i*2+1]
        child = genome.from_parents(p1, p2)
        pop.append(child)

def monogamy(population):
    fitnesses = lambda gs: [g.get_fitness() for g in gs]
    pop = population.get_population()
    num = population.N
    male = [g for g in pop if g.male]
    male_fitnesses = fitnesses(male)
    female = [g for g in pop if not g.male]
    female_fitnesses = fitnesses(female)
    get_relative_fitness = lambda g, fs: len([0 for f in fs if f <= g.get_fitness()]) / len(fs)
    
    pop.clear()
    while num > 0:
        p1 = male.pop()
        p2 = female.pop()
        num_children = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
        for i in range(num_children):
            child = genome.from_parents(p1, p2)
            pop.append(child)
        num -= num_children

def polygamy(population):
    fitnesses = lambda gs: [g.get_fitness() for g in gs]
    pop = population.get_population()
    num = population.N
    male = [g for g in pop if g.male]
    male_fitnesses = fitnesses(male)
    female = [g for g in pop if not g.male]
    female_fitnesses = fitnesses(female)
    get_relative_fitness = lambda g, fs: len([0 for f in fs if f <= g.get_fitness()]) / len(fs)
    pop.clear()
    
    p1 = male.pop()
    i = math.ceil(get_relative_fitness(p1, male_fitnesses) * population.N)
    while num > 0:
        if female:
            p2 = female.pop()
        num_children = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
        for i in range(num_children):
            child = genome.from_parents(p1, p2)
            pop.append(child)
        num -= num_children
        i -= 1
        if i == 0:
            if male:
                p1 = male.pop()
            i = math.ceil(get_relative_fitness(p1, male_fitnesses) * population.N)


def polyandry(population):
    fitnesses = lambda gs: [g.get_fitness() for g in gs]
    pop = population.get_population()
    num = population.N
    male = [g for g in pop if g.male]
    male_fitnesses = fitnesses(male)
    female = [g for g in pop if not g.male]
    female_fitnesses = fitnesses(female)
    get_relative_fitness = lambda g, fs: len([0 for f in fs if f <= g.get_fitness()]) / len(fs)
    pop.clear()
    
    p2 = female.pop()
    i = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
    num_children = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
    while num > 0:
        if male:
            p1 = male.pop()
        for i in range(num_children):
            child = genome.from_parents(p1, p2)
            pop.append(child)
        num -= num_children
        i -= 1
        if i == 0:
            if female:
                p2 = female.pop()
            i = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)

# def polyandry(population):
#     fitnesses = lambda gs: [g.get_fitness() for g in gs]
#     pop = population.get_population()
#     num = population.N
#     male = [g for g in pop if g.male]
#     male_fitnesses = fitnesses(male)
#     female = [g for g in pop if not g.male]
#     female_fitnesses = fitnesses(female)
#     get_relative_fitness = lambda g, fs: len([0 for f in fs if f <= g.get_fitness()]) / len(fs)
#     pop.clear()
#     
#     p2 = female.pop()
#     best = male.pop()
#     i = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
#     while num > 0:
#         if male:
#             male.pop()
#         i = i - 1
#         if i == 0:
#             child = genome.from_parents(best, p2)
#             pop.append(child)
#             num -= 1
#             if female:
#                 p2 = female.pop()
#             if male:
#                 best = male.pop()
#             i = math.ceil(get_relative_fitness(p2, female_fitnesses) * population.N)
 
