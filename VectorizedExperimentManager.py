import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, EngFormatter


class VectorizedExperimentManager:
    def __init__(self) -> None:
        print("Created a new experiment!")

        self.current_generation_number = 0
        self.max_generatoins = config.MAX_NUMBER_OF_GENERATIONS

        self.fitness_vec = None


    def get_population_size(self):
        return self.fitness_vec.size

    def generate_initial_population(self,
                                    size_of_initial_population=config.DEFAULT_SIZE_OF_INITIAL_POPULATION,
                                    beta_production=config.BETA_LACTAMAS_PRODUCTION,
                                    initial_beta_stored=config.INITIAL_BETA_STORED,
                                    multiplication_rate=config.DEFAULT_MULTIPLICATION_RATE):

        self.beta_production = beta_production
        self.multiplication_rate = multiplication_rate

        # Population Feature Vectors
        self.fitness_vec = np.full(shape= size_of_initial_population,
                                   fill_value= min(initial_beta_stored, 1),
                                   dtype=np.float64)
        self.beta_stored_vec = np.full(shape= size_of_initial_population,
                                       fill_value= initial_beta_stored,
                                       dtype= np.float64)

        # Followup data:
        self.pop_size_per_generation = [size_of_initial_population]
        self.max_stored_beta_per_generation = [initial_beta_stored]
        self.avg_stored_beta_per_generation = [initial_beta_stored]

        return

    def run_generation(self, print_info=False):
        self.current_generation_number += 1

        # Kill stochastically
        survived = self.fitness_vec > np.random.rand(self.fitness_vec.size)
        self.fitness_vec = self.fitness_vec[survived]
        self.beta_stored_vec = self.beta_stored_vec[survived]

        # Produce Beta
        self.beta_stored_vec += self.beta_production

        # Split stochastically
        multiplied = self.multiplication_rate > np.random.rand(
            self.fitness_vec.size)
        self.beta_stored_vec[multiplied] /= float(2)

        self.fitness_vec = np.minimum(self.beta_stored_vec, 1)

        self.fitness_vec = np.hstack(
            (self.fitness_vec, self.fitness_vec[multiplied]))
        self.beta_stored_vec = np.hstack(
            (self.beta_stored_vec, self.beta_stored_vec[multiplied]))

        if (self.fitness_vec.size == 0):
            print("! The Bacteria has died out !")
            return

        if print_info:
            print(f"Generation: {self.current_generation_number}\t population size: {self.fitness_vec.size}")
            print(f"Max amount of protein: {np.max(self.beta_stored_vec):.2f}")
            print(f"Min amount of protein: {np.min(self.beta_stored_vec):.2f}")
            # ecoli_data = pd.DataFrame({
            #     "Fitness": self.fitness_vec,
            #     "beta_stored": self.beta_stored_vec
            # })
            # print(ecoli_data)

        # Followup data
        self.pop_size_per_generation.append(self.fitness_vec.size)
        self.max_stored_beta_per_generation.append(np.max(self.beta_stored_vec))
        self.avg_stored_beta_per_generation.append(np.mean(self.beta_stored_vec))


# Plotting

    def plot_max_stored_protein(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Stored Protein Over Generations"
        ylabel = "Max Protein Stored"
        self.general_plot(self.max_stored_beta_per_generation, title, ylabel )
        
    
    def plot_pop_sizes(self):
        """This is a shortcut function for plottin the populaiton size over generations.
        """
        title = "E. coli Population Over Generations"
        ylabel = "Population Size"
        self.general_plot(self.pop_size_per_generation, title, ylabel)

    def plot_avg_stored_protein(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Average Protein Over Generations"
        ylabel = "Average Protein Stored"
        self.general_plot(self.avg_stored_beta_per_generation, title, ylabel )

    def general_plot(self, data: list, title, y_label):
        plt.plot(data, marker='o')
        plt.xlabel("Generation")
        plt.ylabel(y_label)
        plt.ylim(0, max(data) * 1.1)
        plt.gca().yaxis.set_major_formatter(EngFormatter())
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.show()