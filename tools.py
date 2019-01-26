import sys
from secrets import choice
from random import random
from datetime import datetime

import logging
from logging.config import dictConfig

import numpy as np

logger = logging.getLogger()

def shuffle(original_genes, prob=0.038):
    genes = list(original_genes)
    for i, j in enumerate(genes):
        if prob > random():
            k = genes.index(choice(genes))
            genes[i], genes[k] = genes[k], genes[i]
    return tuple(genes)


def chronometer(func):
    def wrap(*args):
        t = datetime.now()
        res = func(*args)
        tempo = datetime.now() - t
        logger.info(f'Running time: {tempo}')
        return res

    return wrap


def read(file_name):
    """
        Lê o arquivo, desconsidera linhas de comentário (iniciadas por #), e retira whitespaces. 
    """
    table = list()

    with open(file_name, 'r') as content:
        table = [[cell.strip() for cell in line.split(';')]
                    for line in content.readlines()
                    if not line.startswith('#')]
    return table


def process_input(table):
    """
        Cria a matriz e as relações entre horarios-linhas
    e entre professores-colunas.
        Retorna a matriz e os 2 dicionarios de relacionamento.

    :param nome_arq: Nome do arquivo
    """

    spots_t = [[spot[0].split('!')[0], indice - 1]
                    for indice, spot in enumerate(table)
                    if spot[0] and not spot[0].startswith('!')]
    spots = dict()
    for i, spot_t in enumerate(spots_t):
        if i + 1 != len(spots_t):
            final = spots_t[i + 1][1]
        else:
            final = len(table) - 1
        spots[spot_t[0]] = tuple(list(range(spot_t[1], final)))

    candidates_t = [[name, index - 1]
                        for index, name in enumerate(table[0])
                        if name]

    candidates = dict()
    for i, name_t in enumerate(candidates_t):
        if i + 1 != len(candidates_t):
            final = candidates_t[i + 1][1]
        else:
            final = len(table) - 1
        candidates[name_t[0]] = tuple(list(range(name_t[1], final)))

    table = np.array([[float(v) for v in value[1:] if v != '']
                        for value in table[1:]])
    return table, spots, candidates


def report(chromosome, table, spots, candidates, spot_descr=None, f=sys.stdout):
    print(f'Mean Satisfaction: {chromosome.fitness * 10:.2f}%', file=f)
    print(chromosome.genes, file=f)
    if chromosome.zero_counter > 0:
        print('Zeros given:', chromosome.zero_counter, file=f)
    inv_spot = {spots[key]: key for key in spots}
    inv_candidates = {candidates[key]: key for key in candidates}
    for cols in inv_candidates:
        print(f'\n{inv_candidates[cols]} ({len(cols)} spots),', file=f)
        mean = 0
        for col in cols:
            line = chromosome.genes.index(col)
            for lines in inv_spot:
                if line in lines:
                    mean += table[line][col]
                    detail = f'(Line: {line});'
                    if spot_descr and spot_descr[line]:
                        detail = f'({spot_descr[line]});'
                    print(f'{inv_spot[lines]} {detail} satisfaction:'.ljust(49, '.'),
                            f'{table[line][col] * 10}%', 
                            file=f)
        mean /= len(cols)
        print(f"Subjects' mean satisfaction:\t{mean * 10:.1f}%\n", end='', file=f)
        print('.' * 60, file=f)
    print('\n', file=f)
