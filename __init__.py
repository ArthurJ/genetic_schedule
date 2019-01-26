
import os
from datetime import datetime
from multiprocessing import cpu_count, Pool, freeze_support

import tools
from population import Population, fit
from chromosome import Schedule


__author__ = '@arthurj'

raw_table = tools.read('ex_realistic_case.csv')
table, spots, candidates = tools.process_input(raw_table)
chromosomes = set()
maxi = 100
count = 0
retries_on_stability = 32
max_fitness = 10

qtd_process = cpu_count()*2

while len(chromosomes) < 200:
    count += 1
    if count > 1000:
        print("[WARN!] Attempts' limit to create subjects reached.")
        break
    try:
        cs = Schedule(table, spots, candidates, fitness_func=fit, fitness_penality_zero=fit)
        chromosomes.add(cs)
    except Exception as e:
        if e.args[0] not in ('CRASH',):
            raise e

print(f'Initial population size: {len(chromosomes)}', '.' * 60, sep='\n')


@tools.chronometer
def iterate(g: Population, 
            max_wait_4_new_fitness=retries_on_stability, 
            max_fitness=max_fitness):
    good_old_days = list()
    for i in range(1, 1000):
        print(f'-- {str(i)}ª Generation --')
        g.next_generation()
        print(g)
        good_old_days.append(g.chromosomes[0].fitness)
        if good_old_days[-1] == max_fitness:
            max_wait_4_new_fitness = int(max_wait_4_new_fitness)/2
        print(f'[{str(good_old_days.count(g.chromosomes[0].fitness))}' +
              'ª occurency of this fitness]', '.' * 60, sep='\n', end='\n')
        if good_old_days.count(g.chromosomes[0].fitness) >= max_wait_4_new_fitness:
            break
    return g


if __name__ == '__main__':
    freeze_support()

    pool = Pool(processes=qtd_process)

    population = Population(table, spots, candidates, list(chromosomes),
                            pool=pool, max_cs_size=maxi, workers=qtd_process)
    
    population = iterate(population)

    result_folder = f"Result date: {datetime.now().strftime('%Y-%m-%d__%H_%M_%S')}"
    
    os.mkdir(result_folder)
    os.chdir(result_folder)
    
    for i, cs in enumerate(population.chromosomes):
        spot_descr=[]

        for x in raw_table[1:]:
            if '!' in x[0]:
                descr = x[0].split('!')[1]
            spot_descr.append(descr)

        with open(f'Schedule {i + 1} (Fitness:{cs.fitness:.2f}).txt', 'w') as output:
            tools.report(cs, table, spots, candidates, f=output,
                            spot_descr=spot_descr)
