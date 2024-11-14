import random

DNA_LENGTH: 900 # number of amino acids in our model
DUPLICATION_ERROR_RATE: 1e-10 


class EColi:
    
    all_ecoli = []
    def __init__(self, id:int, dna:str=None) -> None:
        """Creates a new EColi instance

        Args:
            id (int): The identification of this specific EColi instance
            dna (str, optional): The Gene that will determin the EColi's fitness.
        """
        self.id = id
        if dna: 
            self.dna = dna
        else:
            self.dna = self.generate_random_dna()
    
    def calc_fitness() -> int: 
        return 1
    
    def multiply() -> None:
        NotImplemented

    def generate_random_dna(self) -> str:
        # TODO this function will be deleted later
        # Generate a 900-base DNA string using the bases A, T, C, and G
        return ''.join(random.choice('ACTG') for _ in range(DNA_LENGTH))
    
    def get_dna(self) -> str:
        return self.dna
    

    @staticmethod
    def gene_copy(ecoli:EColi):
        old_dna = ecoli.get_dna()

    
    
    @staticmethod
    def get_all_ecoli():
        return EColi.all_ecoli