import random

class EColi:
    
    all_ecoli = []
    def __init__(self, id, dna) -> None:
        self.id = id
        self.dna = self.generate_random_dna()
    
    def calc_fitness() -> int: 
        return 1
    
    def multiply() -> None:
        NotImplemented

    def generate_random_dna(self) -> str:
        # Generate a 900-base DNA string using the bases A, T, C, and G
        return ''.join(random.choice('ACTG') for _ in range(900))
    
    def get_dna(self) -> str:
        return self.dna
    
    
    
    
    @staticmethod
    def get_all_ecoli():
        return EColi.all_ecoli