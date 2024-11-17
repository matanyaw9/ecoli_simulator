import random               

class ExperimentManager:
    def __init__(self, initial_pop_size:int = 1, max_generations: int = 1) -> None:
        """Initialize the experiment with a population of E. coli and set constraints.

        Args:
            initial_pop_size (int): starting number of E. coli. Defaults to 1.
            max_generations (int): max number of generations. Defaults to 1.
        """
        self.current_genration = 0
        self.max_generations = max_generations
        self.population = set()

    
    def run_generation():
        NotImplemented


experiment = ExperimentManager(1, 5)