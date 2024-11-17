import random
import config

DNA_LENGTH: 900 # number of amino acids in our model
DUPLICATION_ERROR_RATE: 1e-10 


class EColi:

    ecoli_counter = 0
    
    def __init__(self, fitness:float=config.DEFAULT_FITNESS) -> None:
        """Creates a new EColi instance

        Args:
            fitness (float): The fitness of the ecoli, will determin it's probability to survive & reproduce or die. 
        """
        self.id = EColi.ecoli_counter
        self.fitness = fitness
        EColi.ecoli_counter += 1

    
    def calc_fitness(self) -> int: 
        return self.fitness
    
    def multiply(self) -> None:
        parent_fitness = self.fitness
        child_1 = EColi(fitness=parent_fitness)
        child_2 = EColi(fitness=parent_fitness)
        return child_1, child_2
    
    def __str__(self) -> str:
        return f"E. Coli \tid: {self.id}\tfitness: {self.fitness}"
    
    def __repr__(self) -> str:
        return f"EColi(id={self.id}, fitness={self.fitness})"
