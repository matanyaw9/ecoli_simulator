import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from collections.abc import Iterable


# def severity_of_antibiotics(amount_of_antibiotics: float) -> float:
#     # This is effectively a sigmoid function
#     return 1 - 1 / (1 + np.exp(amount_of_antibiotics - config.MINIMUM_INHIBITORY_CONCENTRATION))


def pseudo_sigma(x: float, shift:float=0, steepness:float=1) -> float:
    processed_x = steepness * (x - shift)
    if processed_x < -10:
        return 0
    if processed_x > 10:
        return 1
    return 1 - 1 / (1 + np.exp(processed_x))



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
                                    initial_available_beta=config.INITIAL_AVAILABLE_BETA,
                                    multiplication_rate=config.DEFAULT_MULTIPLICATION_RATE):

        self.INITIAL_SIZE_OF_POPULATION = size_of_initial_population
        self.BETA_PTODUCTION_RATE = beta_production
        self.INITIAL_AVAILABLE_BETA = initial_available_beta
        self.MULTIPLICATION_RATE = multiplication_rate
        self.INITIAL_ANTIBIOITCS = self.antibiotics_amount_in_medium

        self.generation_when_mic_passed = None
        self.total_beta_in_medium = size_of_initial_population * initial_available_beta

        # Population Feature Vectors
        self.bacterial_death_likelihood_vec = np.full(shape=size_of_initial_population,fill_value=1-min(initial_available_beta, 1),dtype=np.float64)
        self.bacterial_beta_available_vec = np.full(shape=size_of_initial_population,fill_value=initial_available_beta,dtype=np.float64)

        # Followup data:
        self.population_size_history = [size_of_initial_population]
        self.max_beta_available_history = [initial_available_beta]
        self.avg_beta_storage_history = [initial_available_beta]
        self.antibiotics_in_medium_history = [self.antibiotics_amount_in_medium]
        self.total_beta_in_medium_history = [self.total_beta_in_medium]
        self.survival_rate_history = []

        return

    def run_generation(self, print_info=False):
        self.current_generation_number += 1

        #  ---------- Kill Stochastically ----------
        antibiotics_severeness = pseudo_sigma(self.antibiotics_amount_in_medium, config.MINIMUM_INHIBITORY_CONCENTRATION)
        self.bacterial_death_likelihood_vec *= antibiotics_severeness

        survived = self.bacterial_death_likelihood_vec < np.random.rand(self.bacterial_death_likelihood_vec.size)
        self.survival_rate_history.append(survived.sum()/survived.size)
        self.bacterial_death_likelihood_vec = self.bacterial_death_likelihood_vec[survived]
        self.bacterial_beta_available_vec = self.bacterial_beta_available_vec[survived]

        # Followup data
        self.population_size_history.append(self.bacterial_death_likelihood_vec.size)

        if (self.bacterial_death_likelihood_vec.size == 0):
            print(
                f"! The Bacteria has Died Out After {self.current_generation_number} Generations !")
            return

        # ---------- Produce Beta ----------
        self.bacterial_beta_available_vec += self.BETA_PTODUCTION_RATE
        self.total_beta_in_medium += self.BETA_PTODUCTION_RATE * self.bacterial_death_likelihood_vec.size

        # ---------- Antibiotics removal ----------
        self.antibiotics_amount_in_medium -= self.total_beta_in_medium
        self.antibiotics_amount_in_medium = max(0, self.antibiotics_amount_in_medium)
        
        # Followup data
        self.antibiotics_in_medium_history.append(self.antibiotics_amount_in_medium) 
        self.max_beta_available_history.append(np.max(self.bacterial_beta_available_vec))
        self.avg_beta_storage_history.append(np.mean(self.bacterial_beta_available_vec))
        self.total_beta_in_medium_history.append(self.total_beta_in_medium)

        if self.antibiotics_amount_in_medium == 0:
            print(
                f"! The Bacteria Wins - No more AntiBiotics After {self.current_generation_number} Generations !")
            return

        # ---------- Split Stochastically ----------
        multiplied = self.MULTIPLICATION_RATE > np.random.rand(self.bacterial_death_likelihood_vec.size)
        self.bacterial_beta_available_vec[multiplied] /= float(2)

        self.bacterial_death_likelihood_vec = np.maximum(1 - self.bacterial_beta_available_vec, 0)  # death likelihood is 1 - fitness
        
    

        self.bacterial_death_likelihood_vec = np.hstack(
            (self.bacterial_death_likelihood_vec, self.bacterial_death_likelihood_vec[multiplied]))
        self.bacterial_beta_available_vec = np.hstack(
            (self.bacterial_beta_available_vec, self.bacterial_beta_available_vec[multiplied]))

        if print_info:
            print(
                f"Generation: {self.current_generation_number}\t population size: {self.bacterial_death_likelihood_vec.size}")
            print(f"Max amount of protein: {np.max(self.bacterial_beta_available_vec):.2f}")
            print(f"Min amount of protein: {np.min(self.bacterial_beta_available_vec):.2f}")

   
    def run_simulation(self, print_info=False):
        while self.antibiotics_amount_in_medium > 0:
            if self.get_population_size() == 0:
                break
            self.run_generation(print_info=print_info)

    def plot_all(self):
        self.plot_pop_sizes_history()
        self.plot_antibiotics_amount_history()
        self.plot_survival_rate_history()
        self.plot_total_beta_in_medium_history()

# Plotting History

    def plot_survival_rate_history(self):
        """This is a shortcut function for plotting the survival rate over generations.
        """
        title = "E. coli Survival Rate Over Generations"
        ylabel = "Survival Rate"
        self.general_plot(self.survival_rate_history, title, ylabel, show_mic_gen_pass=True, prettify_numbers=False)

    def plot_max_available_protein_history(self):
        """This is a shortcut function for plotting the max available protein over generations.
        """
        title = "E. coli Max Available Protein Over Generations"
        ylabel = "Max Protein Available"
        self.general_plot(self.max_beta_available_history, title, ylabel, show_mic_gen_pass=True)

    def plot_pop_sizes_history(self, plot_right_after_MIC=True, show_n_gens_before=0):
        """This is a shortcut function for plottin the populaiton size over generations.
        """
        title = "E. coli Population Over Generations"
        ylabel = "Population Size"
        self.general_plot(self.population_size_history, title, ylabel, show_mic_gen_pass=True)

        if plot_right_after_MIC:
            if not self.generation_when_mic_passed:
                print("MIC has not been passed yet!")
                return
            title = "E. coli Population Right After MIC Passed"
            ylabel = "Population Size"
            start_index = max(0, self.generation_when_mic_passed - show_n_gens_before)
            x_data = range(start_index, len(self.population_size_history))
            self.general_plot(self.population_size_history[start_index:], title, ylabel, x_data=x_data,
                            show_mic_gen_pass=True)

    def plot_avg_available_protein_history(self):
        """This is a shortcut function for plotting the max available protein over generations.
        """
        title = "E. coli Max Average Protein Over Generations"
        ylabel = "Average Protein Available"
        self.general_plot(self.avg_beta_storage_history, title, ylabel, show_mic_gen_pass=True)

    def plot_antibiotics_amount_history(self, plot_right_after_MIC=True, show_n_gens_before=0):
        """This is a shortcut function for plotting the max available protein over generations.
        """
        title = "Antibiotics Amount Over Generations"
        ylabel = "Antibiotics Amount"
        self.general_plot(self.antibiotics_in_medium_history, title, ylabel, show_mic_value=True)
        
        if plot_right_after_MIC:
            if not self.generation_when_mic_passed:
                print("MIC has not been passed yet!")
                return
            # New plot for right after MIC passed
            title = "Antibiotics Amount Right After MIC Passed"
            ylabel = "Antibiotics Amount"
            start_index = max(0, self.generation_when_mic_passed - show_n_gens_before)
            x_data = range(start_index, len(self.antibiotics_in_medium_history))
            self.general_plot(self.antibiotics_in_medium_history[start_index:], title, ylabel, x_data=x_data,
                            show_mic_gen_pass=True)

    def plot_total_beta_in_medium_history(self):
        """This is a shortcut function for plotting the max available protein over generations.
        """
        title = "E. coli Available Protein Over Generations"
        ylabel = "Protein Available"
        self.general_plot(self.total_beta_in_medium_history, title, ylabel, show_mic_gen_pass=True)

    def general_plot(self, y_data: list, title: str, y_label: str, x_data=None, show_mic_value:bool=False, show_mic_gen_pass:bool=False, 
                     prettify_numbers:bool=True):
        if x_data is None:
            x_data = range(len(y_data))    
        if self.generation_when_mic_passed is None:
                self.generation_when_mic_passed = self.find_mic_passing_generation()
        plt.plot(x_data, y_data, marker='o')
        plt.xlabel("Generation")
        plt.ylabel(y_label)
        plt.ylim(0, max(y_data) * 1.1)
        if show_mic_value:
            plt.axhline(y=config.MINIMUM_INHIBITORY_CONCENTRATION, color='r', linestyle='--', linewidth=1,
                         label=f"MIC = {config.MINIMUM_INHIBITORY_CONCENTRATION}")
        if show_mic_gen_pass and self.generation_when_mic_passed:
            plt.axvline(x=self.generation_when_mic_passed, color='r', linestyle='--', linewidth=1,
                         label=f"MIC Passed at {self.generation_when_mic_passed}")
        if prettify_numbers:
            plt.gca().yaxis.set_major_formatter(EngFormatter())
        plt.title(title)
        plt.grid(True)
        plt.legend()
        NUMBER_FORMATTER = EngFormatter()

        # Add model details as a text box
        model_details = (
            f"Initial Population: {NUMBER_FORMATTER(self.INITIAL_SIZE_OF_POPULATION)}\n"
            f"Multiplication Rate: {self.MULTIPLICATION_RATE}\n"
            f"Initial BL: {NUMBER_FORMATTER(self.INITIAL_AVAILABLE_BETA)}\n"
            f"Initial Antibiotics: {NUMBER_FORMATTER(self.INITIAL_ANTIBIOITCS)}\n"
            f"Beta Production Rate: {self.BETA_PTODUCTION_RATE}"
        )
        # Position the text box in the plot (adjust x and y as needed)
        plt.gca().text(0.05, 0.95, model_details, transform=plt.gca().transAxes,
                    fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))


        plt.show()

    def find_mic_passing_generation(self):
        pass_gen = np.argmax(np.array(self.antibiotics_in_medium_history) < config.MINIMUM_INHIBITORY_CONCENTRATION)
        return pass_gen if self.antibiotics_in_medium_history[pass_gen] < config.MINIMUM_INHIBITORY_CONCENTRATION else None