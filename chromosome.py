from tools import shuffle

__author__ = '@arthurj'

class Schedule:
    '''

    '''
    def __init__(self, table, spots, candidates, genes=None):
        if not genes:
            genes = list(range(len(table)))
            genes = tuple(shuffle(genes, prob=.8))
        if Schedule._crash_detection(spots, candidates, genes):
            raise Exception('CRASH')
        
        self.genes = genes
        self.fitness = 0

        self.zero_counter = 0
        self._evaluate_fitness(table)
    
    @classmethod
    def _crash_detection(cls, spots, candidates, genes):
        for hour in spots:
            for name in candidates:
                count = 0
                for i in candidates[name]:
                    for j in spots[hour]:
                        if i is genes[j]:
                            count += 1
                if count >= 2:
                    return True
        return False

    def _evaluate_fitness(self, table, penality_factor=5):
        fitness = 0
        penality_index = len(table)
        for i, j in enumerate(self.genes):
            if not table[i][j]:
                self.zero_counter += 1
                fitness -= penality_factor * penality_index
            else:
                fitness += table[i][j]
        self.fitness = fitness/len(table)
    
    def __le__(self, other):
        if self.fitness == other.fitness and \
                        self.zero_counter > other.zero_counter:
            return True
        if self.fitness <= other.fitness:
            return True
        else:
            return False
    
    def __lt__(self, other):
        if self.fitness == other.fitness and \
                        self.zero_counter > other.zero_counter:
            return True
        if self.fitness < other.fitness:
            return True
        else:
            return False
    
    def __ge__(self, other):
        if self.fitness == other.fitness and \
                        self.zero_counter < other.zero_counter:
            return True
        if self.fitness >= other.fitness:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.fitness == other.fitness and \
                        self.zero_counter < other.zero_counter:
            return True
        if self.fitness > other.fitness:
            return True
        else:
            return False
    
    def __eq__(self, other):
        if self.fitness == other.fitness \
            and self.genes == other.genes:
            return True
        else:
            return False
    
    def __hash__(self):
        return hash(self.genes)

    def __str__(self):
        return   f'  {self.genes} Fitness:{self.fitness:.2f}' + \
                (f'\t({self.zero_counter} zeros)' if self.zero_counter else '')