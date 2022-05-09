import random
import re
import math

haplotypes = [["ATGCACTGTGGATCCAGGTGCCGATATTTATACCAGGGGAACAGCACAACATGA", "ATGCACTGTGGATCAAAAAAAAAAAAAAAAAACCAGGGGAACAGCACAACATGA"], ["ATGGTGATGGACCATGGATTCCAGGATTGGACTAGACCATGGCCCACGTGCTAG", "ATGGTGATGGACCATGGATTCCATTTTTTTTTTTTTTTTTTTTTCACGTGCTAG"]]

class Genome:

    seq1 = ""
    seq2 = ""
    exons1 = []
    exons2 = []
    male = True

    def __init__(self, seq, exons):
        self.seq1 = seq
        self.seq2 = seq
        self.exons1 = exons
        self.exons2 = exons.copy()
        self.male = random.random() < 0.5
        self.insert_haplotypes()

    def get_sequences(self):
        return self.seq1, self.seq2

    def get_exons1(self):
        return [self.seq1[i[0]:i[1]] for i in self.exons1]

    def get_exons2(self):
        return [self.seq2[i[0]:i[1]] for i in self.exons2]

    def set_exons1(self, exons):
        self.exons1 = exons

    def set_exons2(self, exons):
        self.exons2 = exons

    def get_exon_positions(self):
        return [self.exons1, self.exons2]

    def mutate(self, mutation_rate=0.001):
        self.seq1 = ''.join([random.choice('ACGT') if random.random() < mutation_rate else base for base in self.seq1])
        self.seq2 = ''.join([random.choice('ACGT') if random.random() < mutation_rate else base for base in self.seq2])

    def crossover(self, other):
        p = random.randrange(min(len(self.seq1), len(other.seq1)))
        temp = self.seq1[:p] + other.seq1[p:]
        other.seq1 = other.seq1[:p] + self.seq1[p:]
        self.seq1 = temp

        self_exons_copy = self.exons1.copy()
        other_exons_copy = other.exons2.copy()
        self.update_exons(other, p, self_exons_copy, other_exons_copy, self.set_exons1)
        other.update_exons(self, p, other_exons_copy, self_exons_copy, other.set_exons2)

    def update_exons(self, other, p, self_exons_copy, other_exons_copy, set_exons):
        self_exon_index = -1
        self_in_exon = False
        for i in range(len(self_exons_copy)):
            e = self_exons_copy[i]
            if e[0] <= p and e[1] > p:
                self_exon_index = i
                self_in_exon = True
                break
            elif e[0] > p:
                self_exon_index = i
                break
        other_exon_index = -1
        other_in_exon = False
        for i in range(len(other_exons_copy)):
            e = other_exons_copy[i]
            if e[0] <= p and e[1] > p:
                other_exon_index = i
                other_in_exon = True
                break
            elif e[0] > p:
                other_exon_index = i
                break

        if self_exon_index == -1 and other_exon_index == -1:
            return

        if self_exon_index == -1:
            new_exons = self_exons_copy
            if other_in_exon:
                new_exons += other_exons_copy[other_exon_index + 1:]
                if new_exons and new_exons[0][1] < p:
                    new_exons[0][1] = p
            else:
                new_exons += other_exons_copy[other_exon_index:]
            set_exons(new_exons)
            return

        if other_exon_index == -1:
            new_exons = self_exons_copy[:self_exon_index]
            if self_in_exon:
                new_exons.append([self_exons_copy[self_exon_index][0], len(self.seq1)])
            set_exons(new_exons)
            return

        new_exons = []
        if self_in_exon:
            exon_start = self_exons_copy[self_exon_index][0]
            exon_end = p
            try:
                exon_end = other_exons_copy[other_exon_index][1]
            except: pass
            new_exons = self_exons_copy[:self_exon_index] + [[exon_start, exon_end]] + other_exons_copy[other_exon_index + 1:]
        else:
            new_exons = self_exons_copy[:self_exon_index] + other_exons_copy[other_exon_index + 1:]

        if new_exons and type(new_exons[0]) != list:
            new_exons = [new_exons]
        new_exons = []
        for e in new_exons:
            if e not in new_exons:
                new_exons.append(e)
        set_exons(new_exons)

    # def get_fitness(self):
    #     fitness = 1
    #     for exons in [self.get_exons1(), self.get_exons2()]:
    #         for e in exons:
    #             good_streaks = re.findall(r"(?:A{5,}|C{5,})", e)
    #             bad_streaks = re.findall(r"(?:G{10,}|T{10,})", e)
    #             fitness += math.sqrt(sum(map(len, good_streaks))) - sum(map(len, bad_streaks))
    #     return fitness

    # def get_biallelic(self):
    #     alleles = []
    #     for allele_vars in ["AC", "GT"]:
    #         allele = []
    #         for exons in [self.get_exons1(), self.get_exons2()]:
    #             streaks1 = 0
    #             streaks2 = 0
    #             for e in exons:
    #                 streaks1 += len(re.findall(allele_vars[0] + "{5,}", e))
    #                 streaks2 += len(re.findall(allele_vars[1] + "{5,}", e))
    #             allele.append(1 if streaks2 > streaks1 else 0)
    #         alleles.append(allele)
    #     return alleles
            
    def get_fitness(self):
        fitness += sum(len(re.findall(r"[GC]")))
        for alleles in self.get_biallelic():
            for allele in alleles:
                for x in allele:
                    if x >= 0:
                        fitness *= 1.2
                    else:
                        fitness /= 1.2
        return fitness

    def get_biallelic(self, haplotypes=haplotypes):
        alleles = []
        for haplotype in len(haplotypes):
            allele = []
            for exon in genome.get_exons1():
                matched = False
                for i in range(len(haplotype)):
                    variant = haplotype[i]
                    if matches(exon, h):
                        allele.append(i)
                        matched = True
                        break
                if not matched:
                    allele.append(-1)
            for exon in genome.get_exons2():
                matched = False
                for i in range(len(haplotype)):
                    variant = haplotype[i]
                    if matches(exon, h):
                        allele.append(i)
                        matched = True
                        break
                if not matched:
                    allele.append(-1)
            alleles.append(allele)
        return haplotype

    def insert_haplotypes(self):
        for haplotype in haplotypes:
            first_strand = True
            for chromatid in self.get_exon_positions():
                for start, end in chromatid:
                    if start + end > max(map(len, haplotype)):
                        h = random.choice(haplotype)
                        if first_strand:
                            first_strand = False
                            self.seq1 = self.seq1[:start] + h + self.seq1[start + len(h):]
                        else:
                            self.seq2 = self.seq2[:start] + h + self.seq2[start + len(h):]
                        break

def from_parents(p1, p2, crossover_rate=0.1):
    child = Genome(p1.seq1, p1.exons1)
    child.seq2 = p2.seq2
    child.exons2 = p2.exons2
    if random.random() < crossover_rate:
        child.crossover(child)
    return child

default_probabilities = [
    ([0.27, 0.20, 0.23, 0.30], [0.99, 0.1]),    # Intron
    ([0.0005, 0.0001, 0.9993, 0.0001], [0, 1]), # Splice1
    ([0.0001, 0.0069, 0.0001, 0.9929], [0, 1]), # Splice2
    ([0.2, 0.3, 0.3, 0.2], [0.99, 0.1]),        # Exon
    ([0.0005, 0.0001, 0.9993, 0.0001], [0, 1]), # Splice1
    ([0.0001, 0.0069, 0.0001, 0.9929], [0, 1]), # Splice2
    ([0.27, 0.20, 0.23, 0.30], [0.99, 0.1])     # Intron
]

def markov_model(probabilities=default_probabilities):
    states = "ACGT*>"
    return [dict(zip(states, [*p[0], *p[1]])) for p in default_probabilities]

def random_haplotype_gene():
    return random.choice(random.choice(haplotypes))

def matches(a, b, threshold=0.9):
    la = len(a)
    lb = len(b)
    l = min(la, lb)
    diff = abs(la - lb)
    count = 0
    for i in l:
        j = a[i]
        k = b[i]
        if j == k:
            count += 1
    return count / l

def random_genome(min_length=100000, markov=markov_model):
    seq = ""
    length = 0
    exons = []

    mm = markov()
    num_stages = len(mm)
    stage = 0

    while length < min_length:
        l = 1
        threshold = mm[stage]['*']
        while random.random() < threshold:
            l += 1
        if stage == 3: # Mark exon bounds
            exons.append([length, length+l])
        seq += ''.join(random.choices("ACGT", weights=list(mm[stage].values())[:4], k=l))
        length += l
        stage = (stage + 1) % num_stages
    return Genome(seq, exons)
