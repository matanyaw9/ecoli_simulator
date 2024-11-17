import random               
import config
from EColi import EColi

class ExperimentManager:
    def __init__(self, initial_pop_size:int = 1, max_generations: int = 1, fitness:float=config.DEFAULT_FITNESS) -> None:
        """Initialize the experiment with a population of E. coli and set constraints.

        Args:
            initial_pop_size (int): starting number of E. coli. Defaults to 1.
            max_generations (int): max number of generations. Defaults to 1.
        """
        self.current_genration_number = 0
        self.max_generations = max_generations
        self.population = set()
        for _ in range(initial_pop_size):
            self.population.add(EColi(fitness=fitness))

    def get_population(self):
        return self.population
    
    def get_population_size(self)-> int:
        return len(self.population)
    
    def print_population(self) -> None:
        for e in self.population:
            print(e)
    
    def get_generation_number(self)-> int:
        return self.current_genration_number
    

    def run_generation(self):
        """This function runs a single generation of ecoli. for each ecoli it will allow it to survive & reproduce by the probability 
        of it's fitness. 
        """
        new_generation = set()
        while self.population: 
            cur_ecoli = self.population.pop()
            dice = random.random()
            if dice < cur_ecoli.calc_fitness(): # if bacteria survives
                child1, child2 = cur_ecoli.multiply()
                new_generation.update({child1, child2})
        self.population = new_generation
        self.current_genration_number += 1



