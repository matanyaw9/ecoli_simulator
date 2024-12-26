import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from collections.abc import Iterable

TYPES_OF_PLASMIDS = 3
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

    production_per_plasmid = np.array([1, 0.1, 0])

    # The values are to be changed:
    # # Mutation Rates            to:  wt     m1     m2     from:
    # mutation_rate_matrix = np.array([[0.99, 0.005, 0.005], # wt
    #                                  [0.005, 0.99, 0.005], # m1
    #                                  [0.005, 0.005, 0.99]]) # m2
    
    # Mutation Rates            to:  wt     m1     m2     from:
    # MUTATION_RATE_MATRIX = np.array([[1, 0, 0], # wt
    #                                  [0.5, 0.5, 0], # m1
    #                                  [0, 0.5, 0.5]]) # m2


    def __init__(self, 
                 mutation_rate_matrix,
                 max_generations=config.MAX_NUMBER_OF_GENERATIONS,
                 initial_antibiotic_amount=config.DEFAULT_ANTIBIOTICS,
                 MIC=config.MINIMUM_INHIBITORY_CONCENTRATION
                 ) -> None:
        print("Created a New Experiment!")
        
        self.MUTATION_RATE_MATRIX = mutation_rate_matrix
        self.antibiotics_amount_in_medium = initial_antibiotic_amount
        self.MIC = MIC
        self.total_beta_in_medium = 0
        self.current_generation_number = 0
        # TODO is it redundant?
        self.max_generatoins = max_generations

        # self.fitness_vec = None
        self.bacterial_death_likelihood_vec = None    # easier way to approach fitness? it's (1 - fitness)
        

    def get_population_size(self):
        return self.bacterial_death_likelihood_vec.size

    def generate_initial_population(self, wt_num, m1_num, m2_num,
                                    size_of_initial_population=config.DEFAULT_SIZE_OF_INITIAL_POPULATION,
                                    initial_available_beta=config.INITIAL_AVAILABLE_BETA,
                                    multiplication_rate=config.DEFAULT_MULTIPLICATION_RATE):
        
        size_of_initial_population = int(size_of_initial_population)
        self.INITIAL_SIZE_OF_POPULATION = size_of_initial_population
        # self.BETA_PTODUCTION_RATE = beta_production
        self.INITIAL_AVAILABLE_BETA = initial_available_beta
        self.MULTIPLICATION_RATE = multiplication_rate
        self.INITIAL_ANTIBIOITCS = self.antibiotics_amount_in_medium

        self.generation_when_mic_passed = None
        self.total_beta_in_medium = (size_of_initial_population) * initial_available_beta

        # Population Feature Vectors
        self.bacterial_death_likelihood_vec = np.full(shape=size_of_initial_population,fill_value=1-min(initial_available_beta, 1),dtype=np.float64)
        self.bacterial_beta_available_vec = np.full(shape=size_of_initial_population,fill_value=initial_available_beta,dtype=np.float64)
        wt_plasmids = np.ones(size_of_initial_population) * wt_num
        m1_plasmids = np.ones(size_of_initial_population) * m1_num
        m2_plasmids = np.ones(size_of_initial_population) * m2_num
        self.plasmids = np.column_stack((wt_plasmids, m1_plasmids, m2_plasmids))




        # Followup data:
        self.population_size_history = [size_of_initial_population]
        self.max_beta_available_history = [initial_available_beta]
        self.avg_beta_storage_history = [initial_available_beta]
        self.antibiotics_in_medium_history = [self.antibiotics_amount_in_medium]
        self.total_beta_in_medium_history = [self.total_beta_in_medium]
        self.avg_wt_per_cell = [wt_num]
        self.avg_m1_per_cell = [m1_num]
        self.avg_m2_per_cell = [m2_num]
        self.survival_rate_history = []
        self.print_info()
        return

    def mutation_matrix(self, single_cell_plasmids,):
        mut_matrix = np.vstack([np.random.multinomial(single_cell_plasmids[i], self.MUTATION_RATE_MATRIX[i], 1) for i in range(TYPES_OF_PLASMIDS)])
        # print(pass_matrix.shape)
        return mut_matrix

    
    def run_generation(self, print_info=False):
        self.current_generation_number += 1

        #  ---------- Calc Fitness ----------
        antibiotics_severeness = pseudo_sigma(self.antibiotics_amount_in_medium, self.MIC)
        self.bacterial_death_likelihood_vec *= antibiotics_severeness
        

        #  ---------- Kill Stochastically ----------
        survived = self.bacterial_death_likelihood_vec < np.random.rand(self.bacterial_death_likelihood_vec.size)
        
        self.bacterial_death_likelihood_vec = self.bacterial_death_likelihood_vec[survived]
        self.bacterial_beta_available_vec = self.bacterial_beta_available_vec[survived]
        self.plasmids = self.plasmids[survived]

        # Followup data
        self.survival_rate_history.append(survived.sum()/survived.size)
        self.population_size_history.append(self.bacterial_death_likelihood_vec.size)

        if self.bacterial_death_likelihood_vec.size == 0:
            print(f"! The Bacteria has Died Out After {self.current_generation_number} Generations !")
            return

        # ---------- Produce Beta ----------
        produced_beta_vec = np.dot(self.plasmids, self.production_per_plasmid)
        self.bacterial_beta_available_vec += produced_beta_vec
        self.total_beta_in_medium += produced_beta_vec.sum()

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

        # ---------- Stochastic Cell Division  ----------
        multiplied = self.MULTIPLICATION_RATE > np.random.rand(self.bacterial_death_likelihood_vec.size)
        if multiplied.any():
            self.bacterial_beta_available_vec[multiplied] /= float(2)

            #  Update Fitness
            self.bacterial_death_likelihood_vec = np.maximum(1 - self.bacterial_beta_available_vec, 0)  # death likelihood is 1 - fitness

            self.bacterial_death_likelihood_vec = np.hstack((self.bacterial_death_likelihood_vec, self.bacterial_death_likelihood_vec[multiplied]))
            self.bacterial_beta_available_vec = np.hstack((self.bacterial_beta_available_vec, self.bacterial_beta_available_vec[multiplied]))
            self.plasmids = np.row_stack((self.plasmids, self.plasmids[multiplied]))

            new_cells = np.hstack((multiplied, multiplied[multiplied])) # a way to follow what cells are new and what are old

            # Mutate Plasmids:
            new_plasmids = self.mutate_plasmids(new_cells)
            self.plasmids[new_cells] = new_plasmids 

        # Followup plasmids
        self.avg_wt_per_cell.append(np.mean(self.plasmids[:, 0]))
        self.avg_m1_per_cell.append(np.mean(self.plasmids[:, 1]))
        self.avg_m2_per_cell.append(np.mean(self.plasmids[:, 2]))
        
        if print_info:
            self.print_info()
    
    def print_info(self):
        print(f"Generation: {self.current_generation_number}\t pop size: {self.bacterial_death_likelihood_vec.size}")
        # print(self.plasmids)
        # print(f"Max amount of protein: {np.max(self.bacterial_beta_available_vec):.2f}")
        # print(f"Min amount of protein: {np.min(self.bacterial_beta_available_vec):.2f}")

    def mutate_plasmids(self, new_cells):
         # Mutate Plasmids:
        mut_matrixes = []
        for cell in self.plasmids[new_cells]:
            mut_matrixes.append(self.mutation_matrix(cell))
        mut_matrixes = np.stack(mut_matrixes)
        new_plasmids = np.sum(mut_matrixes, axis=1)
        return new_plasmids
   
    def run_simulation(self, print_info=False):
        while self.antibiotics_amount_in_medium > 0:
            if self.get_population_size() == 0:
                break
            if self.current_generation_number >= self.max_generatoins:
                print(f"! The Bacteria Survived {self.max_generatoins} Generations - Experiment Over !")
                break
            self.run_generation(print_info=print_info)

    def plot_all(self):
        self.plot_pop_sizes_history()
        self.plot_antibiotics_amount_history()
        self.plot_survival_rate_history()
        self.plot_total_beta_in_medium_history()
        self.plot_plasmids_history()

# Plotting History

    def plot_plasmids_history(self):
        """This is a shortcut function for plotting the plasmids over generations.
        """
        title = "E. coli Plasmids Over Generations"
        ylabel = "Plasmids"
        plasmids_history = np.array([self.avg_wt_per_cell, self.avg_m1_per_cell, self.avg_m2_per_cell])
        self.general_plot(plasmids_history, title, ylabel, show_mic_gen_pass=True, prettify_numbers=False, labels=["WT", "M1", "M2"])

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

    # # def general_plot_old(self, y_data: list, title: str, y_label: str,
    #                      x_data=None, show_mic_value:bool=False,
    #                      show_mic_gen_pass:bool=False, prettify_numbers:bool=True):
    #     if x_data is None:
    #         x_data = range(len(y_data))    
    #     if self.generation_when_mic_passed is None:
    #             self.generation_when_mic_passed = self.find_mic_passing_generation()
    #     plt.plot(x_data, y_data, marker='o')
    #     plt.xlabel("Generation")
    #     plt.ylabel(y_label)
    #     plt.ylim(0, max(y_data) * 1.1)
    #     if show_mic_value:
    #         plt.axhline(y=self.MIC, color='r', linestyle='--', linewidth=1,
    #                      label=f"MIC = {self.MIC}")
    #     if show_mic_gen_pass and self.generation_when_mic_passed:
    #         plt.axvline(x=self.generation_when_mic_passed, color='r', linestyle='--', linewidth=1,
    #                      label=f"MIC Passed at {self.generation_when_mic_passed}")
    #     if prettify_numbers:
    #         plt.gca().yaxis.set_major_formatter(EngFormatter())
    #     plt.title(title)
    #     plt.grid(True)
    #     plt.legend()
    #     NUMBER_FORMATTER = EngFormatter()

    #     # Add model details as a text box
    #     model_details = (
    #         f"Initial Population: {NUMBER_FORMATTER(self.INITIAL_SIZE_OF_POPULATION)}\n"
    #         f"Multiplication Rate: {self.MULTIPLICATION_RATE}\n"
    #         f"Initial BL: {NUMBER_FORMATTER(self.INITIAL_AVAILABLE_BETA)}\n"
    #         f"Initial Antibiotics: {NUMBER_FORMATTER(self.INITIAL_ANTIBIOITCS)}\n"
    #         # f"Beta Production Rate: {self.BETA_PTODUCTION_RATE}"
    #     )
    #     # Position the text box in the plot (adjust x and y as needed)
    #     plt.gca().text(0.05, 0.95, model_details, transform=plt.gca().transAxes,
    #                 fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))


    #     plt.show()


    def general_plot(self, y_data: list, title: str, y_label: str, x_data=None, show_mic_value: bool = False, 
                 show_mic_gen_pass: bool = False, prettify_numbers: bool = True,
                 labels=None):
        """
        A general plot function to handle either single or multiple plots on the same figure.
        
        :param y_data: A list of y-values (single list or list of lists for multiple plots).
        :param title: Title of the plot.
        :param y_label: Label for the y-axis.
        :param x_data: x-values (can be None, in which case indices will be used).
        :param show_mic_value: Whether to show the MIC value line.
        :param show_mic_gen_pass: Whether to show the MIC passing generation line.
        :param prettify_numbers: Whether to format the numbers for better readability.
        """
        if x_data is None:
            # Handle x_data for multiple plots
            if isinstance(y_data[0], Iterable):
                x_data = [range(len(y)) for y in y_data]
            else:
                x_data = range(len(y_data))

        # Ensure y_data is treated as a list of lists for consistency
        if not isinstance(y_data[0], Iterable):
            y_data = [y_data]
            x_data = [x_data]

        if self.generation_when_mic_passed is None:
            self.generation_when_mic_passed = self.find_mic_passing_generation()

        # Plot each dataset
        for i, (x, y) in enumerate(zip(x_data, y_data)):
            plt.plot(x, y, marker='o', label=f"{labels[i]}" if labels else None)

        NUMBER_FORMATTER = EngFormatter()
        plt.xlabel("Generation")
        plt.ylabel(y_label)
        plt.ylim(0, np.max(y_data) * 1.1)

        if show_mic_value:
            plt.axhline(y=self.MIC, color='r', linestyle='--', linewidth=1,
                        label=f"MIC = {self.MIC}")
        if show_mic_gen_pass and self.generation_when_mic_passed:
            plt.axvline(x=self.generation_when_mic_passed, color='r', linestyle='--', linewidth=1,
                        label=f"MIC Passed at {NUMBER_FORMATTER(self.generation_when_mic_passed)}")

        if prettify_numbers:
            plt.gca().yaxis.set_major_formatter(EngFormatter())

        plt.title(title)
        plt.grid(True)
        plt.legend()

        # Add model details as a text box
        model_details = (
            f"Initial Population: {NUMBER_FORMATTER(self.INITIAL_SIZE_OF_POPULATION)}\n"
            f"Multiplication Rate: {self.MULTIPLICATION_RATE}\n"
            f"Initial BL: {NUMBER_FORMATTER(self.INITIAL_AVAILABLE_BETA)}\n"
            f"Initial Antibiotics: {NUMBER_FORMATTER(self.INITIAL_ANTIBIOITCS)}\n"
        )
        plt.gca().text(0.05, 0.95, model_details, transform=plt.gca().transAxes,
                    fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))

        plt.show()

    def find_mic_passing_generation(self):
        pass_gen = np.argmax(np.array(self.antibiotics_in_medium_history) < self.MIC)
        return pass_gen if self.antibiotics_in_medium_history[pass_gen] < self.MIC else None