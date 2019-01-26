from time import sleep
from copy import deepcopy as copy
from secrets import choice
from random import randint
from itertools import combinations
import logging
from logging.config import dictConfig

import numpy as np

from tools import shuffle, chronometer
from chromosome import Schedule

__author__ = '@arthurj'


logger = logging.getLogger()

def merge(genes_parent_a, genes_parent_b):
    crosspoint = randint(1, len(genes_parent_a) - 1)

    gene_set = set(genes_parent_a + genes_parent_b)
    genes_a = list(shuffle(genes_parent_a))
    genes_b = list(shuffle(genes_parent_b))

    genes_f = [genes_a[:crosspoint], genes_b[:crosspoint]]
    adjuncts = [genes_b[crosspoint:], genes_a[crosspoint:]]

    for genes in genes_f:
        adjunct = adjuncts.pop(0)
        for gene in adjunct:
            if gene not in genes:
                genes.append(gene)
            else:
                genes.append(gene_set.difference(genes + adjunct).pop())
    return tuple(genes_f[0]), tuple(genes_f[1])


def fit(x):
    return x

def penality(x):
    return 5*x

def sex(table, spots, candidates, couples):
    offspring = set()
    for couple in couples:
        combined = merge(couple[0].genes, couple[1].genes)
        try:
            some_a = Schedule(table, spots, candidates, fit, penality, combined[0])
            offspring.add(some_a)
        except Exception as e:
            if e.args[0] != 'CRASH':
                raise e
        try:
            some_b = Schedule(table, spots, candidates, fit, penality, combined[1])
            offspring.add(some_b)
        except Exception as e:
            if e.args[0] != 'CRASH':
                raise e
    return offspring


def assex(chromosomes, table, spots, candidates):
    clones = set()
    for some in chromosomes:
        try:
            clones.add(Schedule(table,
                                spots,
                                candidates,
                                fit, 
                                penality,
                                shuffle(list(some.genes))))
        except Exception as e:
            if e.args[0] != 'CRASH':
                raise e
    return clones


class Population:

    def __init__(self, table, spots, candidates, chromosomes, pool,
                    workers=4, max_cs_size=85):

        self.workers = workers
        self.pool = pool

        if chromosomes is None or len(chromosomes) < self.workers:
            raise Exception('Not enough chromosomes to set a population.')

        self.table = np.array(table)
        self.spots = spots
        self.candidates = candidates

        self.limits = dict([('max', max_cs_size),
                            ('min', int(max_cs_size * .3)),
                            ('top', int(max_cs_size * .05))])

        self.select(chromosomes)

        if len(self.chromosomes) is 0 or self.chromosomes is None:
            raise Exception(
                'The population could not overcome the retrictions.')

    def select(self, chromosomes):
        cs = [copy(s) for s in chromosomes]
        cs.sort(reverse=True)

        set_cs = set(cs[:self.limits['min']])

        for i in range(len(cs)):
            chosen = choice(cs)
            cs.remove(chosen)
            set_cs.add(chosen)
            if len(set_cs) >= self.limits['max']:
                break
        self.chromosomes = sorted(list(set_cs), reverse=True)

    @property
    def statistics(self):
        fitness = [somo.fitness for somo in self.chromosomes]
        statistics = [0, 0, 0, 0]
        statistics[0] = self.chromosomes[0].fitness
        statistics[1] = self.chromosomes[-1].fitness
        statistics[2] = np.mean(fitness)
        statistics[3] = np.std(fitness)
        return tuple(statistics)

    @chronometer
    def next_generation(self):
        new_chromosomes = set()
        [new_chromosomes.add(somo)
            for somo in self.chromosomes[:self.limits['top']]]  # keep the best of the last generation

        couples = list(combinations(self.chromosomes, 2))

        procs = []

        step = int(len(self.chromosomes) / self.workers) or 1
        slices = [[i, i + step] for i in range(0, len(self.chromosomes), step)
                    if i < len(self.chromosomes)]
        slices[-1][-1] = -1
        [procs.append(self.pool.apply_async(assex,
                                            (
                                                self.chromosomes[i:j], 
                                                self.table,
                                                self.spots, 
                                                self.candidates)))
            for i, j in slices]

        step = int(len(couples) / self.workers) or 1
        slices = [[i, i + step] for i in range(0, len(couples), step)
                    if i < len(couples)]
        slices[-1][-1] = -1
        [procs.append(self.pool.apply_async(sex, (
                                                self.table,
                                                self.spots,
                                                self.candidates,
                                                couples[i:j])))
            for i, j in slices]
        
        while True:
            if np.all([proc.ready() for proc in procs]):
                [new_chromosomes.update(proc.get()) for proc in procs]
                break
            sleep(.027)
        logger.info(f'{len(new_chromosomes)} new subjects created successfully.')

        self.chromosomes = list()
        self.select(new_chromosomes)

        return self

    def __str__(self):
        s = f'\tPopulation size: {len(self.chromosomes)}.\n'
        s += '\tBest:{0:.2f}, Worst:{1:.2f}, Mean:{2:.2f}, Standard Deviantion:{3:.3f};\n'\
                .format(*self.statistics)
        s += f'\tBest subject at current generation:\n{self.chromosomes[0]}'
        return s

    def __repr__(self):
        string = self.__str__() + '\n'
        for item in self.chromosomes:
            string += str(item) + '\n'
        return string
