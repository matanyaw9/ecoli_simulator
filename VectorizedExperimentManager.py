import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter


class VectorizedExperimentManager:
    def __init__(self, max_generations=config.MAX_NUMBER_OF_GENERATIONS,
                 antibiotics_amount=config.DEFAULT_ANTIBIOTICS) -> None:
        print("Created a New Experiment!")

        # TODO add the antibiotics of the experiment
        self.antibiotics_amount = antibiotics_amount
        self.total_beta_in_medium = 0
        self.current_generation_number = 0
        # TODO is it redundant?
        self.max_generatoins = max_generations

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
        self.total_beta_in_medium = size_of_initial_population * initial_beta_stored

        # Population Feature Vectors
        self.fitness_vec = np.full(shape=size_of_initial_population,
                                   fill_value=min(initial_beta_stored, 1),
                                   dtype=np.float64)
        self.beta_stored_vec = np.full(shape=size_of_initial_population,
                                       fill_value=initial_beta_stored,
                                       dtype=np.float64)

        # Followup data:
        self.population_size_history = [size_of_initial_population]
        self.max_beta_stored_history = [initial_beta_stored]
        self.avg_beta_storage_history = [initial_beta_stored]
        self.antibiotics_in_medium_history = [self.antibiotics_amount]
        self.total_beta_in_medium_history = [self.total_beta_in_medium]

        return

    def run_generation(self, print_info=False):
        self.current_generation_number += 1

        # Kill stochastically
        survived = self.fitness_vec > np.random.rand(self.fitness_vec.size)
        self.fitness_vec = self.fitness_vec[survived]
        self.beta_stored_vec = self.beta_stored_vec[survived]

        if (self.fitness_vec.size == 0):
            print(
                f"! The Bacteria has Died Out After {self.current_generation_number} Generations!")
            return

        # Produce Beta
        self.beta_stored_vec += self.beta_production
        self.total_beta_in_medium += self.beta_production * self.fitness_vec.size

        # Antibiotics removal
        self.antibiotics_amount -= self.total_beta_in_medium
        self.antibiotics_amount = max(0, self.antibiotics_amount)
        if self.antibiotics_amount < 0:
            print(
                f"! The Bacteria Wins - No more AntiBiotics After {self.current_generation_number} Generations !")

        # Split stochastically
        multiplied = self.multiplication_rate > np.random.rand(
            self.fitness_vec.size)
        self.beta_stored_vec[multiplied] /= float(2)

        self.fitness_vec = np.minimum(
            self.beta_stored_vec, 1)  # fitness is probability

        self.fitness_vec = np.hstack(
            (self.fitness_vec, self.fitness_vec[multiplied]))
        self.beta_stored_vec = np.hstack(
            (self.beta_stored_vec, self.beta_stored_vec[multiplied]))

        if print_info:
            print(
                f"Generation: {self.current_generation_number}\t population size: {self.fitness_vec.size}")
            print(f"Max amount of protein: {np.max(self.beta_stored_vec):.2f}")
            print(f"Min amount of protein: {np.min(self.beta_stored_vec):.2f}")
            # ecoli_data = pd.DataFrame({
            #     "Fitness": self.fitness_vec,
            #     "beta_stored": self.beta_stored_vec
            # })
            # print(ecoli_data)

        # Followup data
        self.population_size_history.append(self.fitness_vec.size)
        self.max_beta_stored_history.append(np.max(self.beta_stored_vec))
        self.avg_beta_storage_history.append(np.mean(self.beta_stored_vec))
        self.antibiotics_in_medium_history.append(self.antibiotics_amount)
        self.total_beta_in_medium_history.append(self.total_beta_in_medium)


# Plotting History

    def plot_max_stored_protein_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Stored Protein Over Generations"
        ylabel = "Max Protein Stored"
        self.general_plot(self.max_beta_stored_history, title, ylabel)

    def plot_pop_sizes_history(self):
        """This is a shortcut function for plottin the populaiton size over generations.
        """
        title = "E. coli Population Over Generations"
        ylabel = "Population Size"
        self.general_plot(self.population_size_history, title, ylabel)

    def plot_avg_stored_protein_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Average Protein Over Generations"
        ylabel = "Average Protein Stored"
        self.general_plot(self.avg_beta_storage_history, title, ylabel)

    def plot_antibiotics_amount_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "Antibiotics Amount Over Generations"
        ylabel = "Antibiotics Amount"
        self.general_plot(self.antibiotics_in_medium_history, title, ylabel)

    def plot_total_beta_stored_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Stored Protein Over Generations"
        ylabel = "Protein Stored"
        self.general_plot(self.total_beta_in_medium_history, title, ylabel)

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
