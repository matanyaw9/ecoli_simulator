import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from collections.abc import Iterable


class VectorizedExperimentManager:
    def __init__(self, max_generations=config.MAX_NUMBER_OF_GENERATIONS,
                 initial_antibiotic_amount=config.DEFAULT_ANTIBIOTICS) -> None:
        print("Created a New Experiment!")

        self.antibiotics_amount_in_medium = initial_antibiotic_amount
        self.total_beta_in_medium = 0
        self.current_generation_number = 0
        # TODO is it redundant?
        self.max_generatoins = max_generations

        # self.fitness_vec = None
        self.bacterial_death_likelihood_vec = None    # easier way to approach fitness? it's (1 - fitness)

    def get_population_size(self):
        return self.bacterial_death_likelihood_vec.size

    def generate_initial_population(self,
                                    size_of_initial_population=config.DEFAULT_SIZE_OF_INITIAL_POPULATION,
                                    beta_production=config.BETA_LACTAMAS_PRODUCTION,
                                    initial_beta_stored=config.INITIAL_BETA_STORED,
                                    multiplication_rate=config.DEFAULT_MULTIPLICATION_RATE):

        self.generation_when_mic_passed = None
        self.beta_production = beta_production
        self.multiplication_rate = multiplication_rate
        self.total_beta_in_medium = size_of_initial_population * initial_beta_stored

        # Population Feature Vectors
        # self.fitness_vec = np.full(shape=size_of_initial_population,fill_value=min(initial_beta_stored, 1),dtype=np.float64)
        self.bacterial_death_likelihood_vec = np.full(shape=size_of_initial_population,fill_value=1-min(initial_beta_stored, 1),dtype=np.float64)
        self.bacterial_beta_stored_vec = np.full(shape=size_of_initial_population,fill_value=initial_beta_stored,dtype=np.float64)

        # Followup data:
        self.population_size_history = [size_of_initial_population]
        self.max_beta_stored_history = [initial_beta_stored]
        self.avg_beta_storage_history = [initial_beta_stored]
        self.antibiotics_in_medium_history = [self.antibiotics_amount_in_medium]
        self.total_beta_in_medium_history = [self.total_beta_in_medium]
        self.survival_rate_history = []

        return

    def run_generation(self, print_info=False):
        self.current_generation_number += 1

        #  ---------- Kill Stochastically ----------
        survived = self.bacterial_death_likelihood_vec < np.random.rand(self.bacterial_death_likelihood_vec.size)
        self.survival_rate_history.append(survived.sum()/survived.size)
        self.bacterial_death_likelihood_vec = self.bacterial_death_likelihood_vec[survived]
        self.bacterial_beta_stored_vec = self.bacterial_beta_stored_vec[survived]

        # Followup data
        self.population_size_history.append(self.bacterial_death_likelihood_vec.size)

        if (self.bacterial_death_likelihood_vec.size == 0):
            print(
                f"! The Bacteria has Died Out After {self.current_generation_number} Generations !")
            return

        # ---------- Produce Beta ----------
        self.bacterial_beta_stored_vec += self.beta_production
        self.total_beta_in_medium += self.beta_production * self.bacterial_death_likelihood_vec.size

        # ---------- Antibiotics removal ----------
        self.antibiotics_amount_in_medium -= self.total_beta_in_medium
        self.antibiotics_amount_in_medium = max(0, self.antibiotics_amount_in_medium)
        
        # Followup data
        self.antibiotics_in_medium_history.append(self.antibiotics_amount_in_medium) 
        self.max_beta_stored_history.append(np.max(self.bacterial_beta_stored_vec))
        self.avg_beta_storage_history.append(np.mean(self.bacterial_beta_stored_vec))
        self.total_beta_in_medium_history.append(self.total_beta_in_medium)

        if self.antibiotics_amount_in_medium == 0:
            print(
                f"! The Bacteria Wins - No more AntiBiotics After {self.current_generation_number} Generations !")
            return

        # ---------- Split Stochastically ----------
        multiplied = self.multiplication_rate > np.random.rand(self.bacterial_death_likelihood_vec.size)
        self.bacterial_beta_stored_vec[multiplied] /= float(2)

        self.bacterial_death_likelihood_vec = np.maximum(1 - self.bacterial_beta_stored_vec, 0)  # death likelihood is 1 - fitness
        # If there is not much antibiotics on the plate, the bacteria will be more resistant
        if self.antibiotics_amount_in_medium < config.MINIMUM_INHIBITORY_CONCENTRATION:
            self.generation_when_mic_passed = self.current_generation_number if self.generation_when_mic_passed is None else self.generation_when_mic_passed
            remaining_of_MIC = self.antibiotics_amount_in_medium / config.MINIMUM_INHIBITORY_CONCENTRATION
            # *************************************************************************************************************
            self.bacterial_death_likelihood_vec *= remaining_of_MIC
            # *************************************************************************************************************


        # self.fitness_vec = np.hstack(
            # (self.fitness_vec, self.fitness_vec[multiplied]))
        self.bacterial_death_likelihood_vec = np.hstack(
            (self.bacterial_death_likelihood_vec, self.bacterial_death_likelihood_vec[multiplied]))
        self.bacterial_beta_stored_vec = np.hstack(
            (self.bacterial_beta_stored_vec, self.bacterial_beta_stored_vec[multiplied]))

        if print_info:
            print(
                f"Generation: {self.current_generation_number}\t population size: {self.bacterial_death_likelihood_vec.size}")
            print(f"Max amount of protein: {np.max(self.bacterial_beta_stored_vec):.2f}")
            print(f"Min amount of protein: {np.min(self.bacterial_beta_stored_vec):.2f}")
            # ecoli_data = pd.DataFrame({
            #     "Fitness": self.fitness_vec,
            #     "beta_stored": self.beta_stored_vec
            # })
            # print(ecoli_data)




# Plotting History

    def plot_survival_rate_history(self):
        """This is a shortcut function for plotting the survival rate over generations.
        """
        title = "E. coli Survival Rate Over Generations"
        ylabel = "Survival Rate"
        self.general_plot(self.survival_rate_history, title, ylabel, show_mic_gen_pass=True, prettify_numbers=False)

    def plot_max_stored_protein_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Stored Protein Over Generations"
        ylabel = "Max Protein Stored"
        self.general_plot(self.max_beta_stored_history, title, ylabel, show_mic_gen_pass=True)

    def plot_pop_sizes_history(self):
        """This is a shortcut function for plottin the populaiton size over generations.
        """
        title = "E. coli Population Over Generations"
        ylabel = "Population Size"
        self.general_plot(self.population_size_history, title, ylabel, show_mic_gen_pass=True)

    def plot_avg_stored_protein_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Max Average Protein Over Generations"
        ylabel = "Average Protein Stored"
        self.general_plot(self.avg_beta_storage_history, title, ylabel, show_mic_gen_pass=True)

    def plot_antibiotics_amount_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "Antibiotics Amount Over Generations"
        ylabel = "Antibiotics Amount"
        self.general_plot(self.antibiotics_in_medium_history, title, ylabel, show_mic_value=True)

    def plot_total_beta_stored_history(self):
        """This is a shortcut function for plotting the max stored protein over generations.
        """
        title = "E. coli Stored Protein Over Generations"
        ylabel = "Protein Stored"
        self.general_plot(self.total_beta_in_medium_history, title, ylabel, show_mic_gen_pass=True)

    def plot_bacteria_amount_history_right_after_mic_passed(self, show_n_gens_before=3):
        title = "E. coli Population Right After MIC Passed"
        ylabel = "Population Size"
        start_index = max(0, self.generation_when_mic_passed - show_n_gens_before)
        x_data = range(start_index, len(self.population_size_history))
        self.general_plot(self.population_size_history[start_index:], title, ylabel, x_data=x_data,
                          show_mic_gen_pass=True)
        
    def plot_antibiotics_amount_history_right_after_mic_passed(self, show_n_gens_before=3):
        title = "Antibiotics Amount Right After MIC Passed"
        ylabel = "Antibiotics Amount"
        start_index = max(0, self.generation_when_mic_passed - show_n_gens_before)
        x_data = range(start_index, len(self.antibiotics_in_medium_history))
        self.general_plot(self.antibiotics_in_medium_history[start_index:], title, ylabel, x_data=x_data,
                          show_mic_gen_pass=True)

    def general_plot(self, y_data: list, title, y_label, x_data=None, show_mic_value:bool=False, show_mic_gen_pass:bool=False, 
                     prettify_numbers:bool=True):
        if x_data is None:
            x_data = range(len(y_data))        
        plt.plot(x_data, y_data, marker='o')
        plt.xlabel("Generation")
        plt.ylabel(y_label)
        plt.ylim(0, max(y_data) * 1.1)
        if show_mic_value:
            plt.axhline(y=config.MINIMUM_INHIBITORY_CONCENTRATION, color='r', linestyle='--', linewidth=1, label=f"MIC = {show_mic_value}")
        if show_mic_gen_pass:
            plt.axvline(x=self.generation_when_mic_passed, color='r', linestyle='--', linewidth=1,
                         label=f"MIC Passed at {self.generation_when_mic_passed}")
        if prettify_numbers:
            plt.gca().yaxis.set_major_formatter(EngFormatter())
        plt.title(title)
        plt.grid(True)
        plt.legend()

        plt.show()
