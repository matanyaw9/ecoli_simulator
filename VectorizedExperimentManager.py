import numpy as np
import config

class VectorizedExperimentManager:
    def __init__(self) -> None:
        self.current_generation_number = 0
        self.max_generatoins = config.MAX_NUMBER_OF_GENERATIONS

        self.fitness_vec = None

    
    def generate_initial_population(self, size_of_initial_population=config.DEFAULT_SIZE_OF_INITIAL_POPULATION):
        self.fitness_vec = np.full(size_of_initial_population, config.DEFAULT_FITNESS)
        # TODO might be redundant to return the fitness vector
        return self.fitness_vec
    
    def run_generation(self):
        survived = self.fitness_vec > np.random.rand(self.fitness_vec.size)
        
        self.fitness_vec = self.fitness_vec[survived]
        if (self.fitness_vec.size == 0):
            print("! The Bacteria has died out !")
            return
        self.fitness_vec = np.hstack((self.fitness_vec, self.fitness_vec))
        print(f"population size: {self.fitness_vec.size}")

