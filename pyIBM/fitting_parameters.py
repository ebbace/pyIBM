import random
import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable
from concurrent.futures import ProcessPoolExecutor
import time
import datetime as dt
from functions.plot import plot_levels
from functions.E import give_predictions, calculate_rmse
from functions.read import read_experimental_data, read_variables, read_cg_coefficients


"""
Genetic Algorithm for Nuclear Parameter Fitting
-----------------------------------------------
This script applies a genetic algorithm to fit IBM parameters to experimental
nuclear energy levels. It evolves parameter sets across generations to minimize
the RMSE between theoretical predictions and experimental data.


Workflow:
1. Define fixed and free parameters.
2. Generate an initial population of parameter sets.
3. Evaluate fitness using RMSE against experimental data.
4. Apply selection, crossover, and mutation to evolve populations.
5. Track best/worst fitness per generation.
6. Plot fitness evolution and predicted vs experimental energy levels.
"""

# --- Settings--------------------------------------------------------------
nucleus = "Te"  # Target nucleus
A = 108  # Mass number
N = 4  # Boson number

# Free parameters: specify bounds [low, high]. Fixed parameters: specify constant value.
fixed = {"chi": 0.4,
         "a1": [-1000, 0],
         "a2": [-20, 0],
         "a3": 0,
         "a4": [-20, 0],
         "a5": 0,
         "a6": 0,
         "aF": 0,
         }

# --- Algorithm parameters--------------------------------------------------
population_size = 150 # Number of individuals per generation (>= 3)
generations = 20 # Maximum number of generations
target_fitness = -10 # Stopping criterion: fitness target (closer to 0 is better)
mutation_rate = 1 # Mutation probability (0 = none, 1 = always)

# --- Initialization--------------------------------------------------------
print("Starting fitting of {}-{}".format(nucleus, A))
cg = read_cg_coefficients(N)

# Check if nucleus has even mass number
even = (A % 2 == 0)


def fitness_function(params):
    """
    Evaluate fitness of a parameter set.
    Fitness = RMSE (negative) between predicted and experimental energy levels.

    Args:
    params (tuple): Parameter values in order according to `fixed` dict.

    Returns:
    float: Fitness value (RMSE, negative).
    """
    variables = {"N": N,
                 "chi": params[0],
                 "a1": params[1],
                 "a2": params[2],
                 "a3": params[3],
                 "a4": params[4],
                 "a5": params[5],
                 "a6": params[6],
                 "aF": params[7]}
    targets, spins, bandhead = read_experimental_data(A, nucleus)
    predictions, _ = give_predictions(A, nucleus, variables, cg, bandhead, spins, even)
    rmse = calculate_rmse(predictions, targets)
    return rmse  # Negated, find maximum


# --- Population Initialization---------------------------------------------
def create_initial_population(size, parameters):
    """
    Create an initial population of individuals.

    Args:
    size (int): Population size.
    parameters (dict): Parameter bounds or fixed values.

    Returns:
    list[tuple]: List of parameter sets.
    """
    population = []
    for _ in range(size):
        individual = np.array([])
        for p in parameters.keys():
            if isinstance(fixed[p], list):
                individual = np.append(individual, random.uniform(fixed[p][0], fixed[p][1]))
            else:
                individual = np.append(individual, fixed[p])
        population.append(tuple(individual))
    return population


# --- Selection-------------------------------------------------------------
def selection(population, fitnesses, tournament_size=3):
    """
    Tournament selection: choose best among random subsets.

    Args:
    population (list[tuple]): Current population.
    fitnesses (list[float]): Fitness values.
    tournament_size (int): Number of individuals per tournament.

    Returns:
    list[tuple]: Selected individuals.
    """
    selected = []
    for _ in range(len(population)):
        tournament = random.sample(list(zip(population, fitnesses)), tournament_size)
        winner = max(tournament, key=lambda x: x[1])[0]
        selected.append(winner)
    return selected


# --- Crossover-------------------------------------------------------------
def crossover(parent1, parent2, fixed):
    """
    Blend crossover between two parents.

    Args:
    parent1 (tuple): First parent.
    parent2 (tuple): Second parent.
    fixed (dict): Parameter constraints.

    Returns:
    tuple: Two offspring individuals.
    """
    alpha = random.random()
    c1 = np.array([])
    c2 = np.array([])
    for j, p in enumerate(fixed.keys()):
        if isinstance(fixed[p], list):
            c1 = np.append(c1, alpha * parent1[j] + (1 - alpha) * parent2[j])
            c2 = np.append(c2, alpha * parent2[j] + (1 - alpha) * parent1[j])
        else:
            c1 = np.append(c1, fixed[p])
            c2 = np.append(c2, fixed[p])
    child1 = tuple(c1)
    child2 = tuple(c2)
    return child1, child2


# --- Mutation--------------------------------------------------------------
def mutation(individual, mutation_rate, fixed):
    """
    Mutate individual parameters with given probability.

    Args:
    individual (tuple): Parameter set.
    mutation_rate (float): Probability of mutation.
    fixed (dict): Parameter constraints.

    Returns:
    tuple: Mutated individual.
    """
    individual = list(individual)
    for i in range(len(individual)):
        p = list(fixed.keys())[i]
        if isinstance(fixed[p], list):
            if random.random() < mutation_rate:
                m = (fixed[p][1] - fixed[p][0])/10
                mutation_amount = random.uniform(-m, m)
                individual[i] += mutation_amount
                # Ensure the individual stays within bounds
                individual[i] = max(min(individual[i], fixed[p][1]), fixed[p][0])
        else:
            pass
    return tuple(individual)


# --- Genetic Algorithm Main Loop-------------------------------------------
def genetic_algorithm(population_size, generations, mutation_rate, fixed, A, nucleus):
    """
    Run the genetic algorithm optimization.

    Args:
    population_size (int): Number of individuals per generation.
    generations (int): Maximum number of generations.
    mutation_rate (float): Probability of mutation.
    fixed (dict): Parameter bounds and constants.
    A (int): Mass number.
    nucleus (str): Nucleus symbol.

    Returns:
    tuple: Best solution found.
    """
    population = create_initial_population(population_size, fixed)

    # Tracking best/worst fitness per generation
    best_performers = []
    worst_performers = []

    # Pretty table to summarize progress
    table = PrettyTable()
    table.field_names = ["Generation", "chi", "a1", "a2", "a3", "a4", "a5", "a6", "aF", "Fitness"]
    best_fitness = -10000
    generation = 0
    with ProcessPoolExecutor() as executor:
        while best_fitness < target_fitness and generation < generations:
            start_gen = time.time()

            # Evaluate fitness of all individuals in parallel
            fitnesses = list(executor.map(fitness_function, population))

            # Track best and worst performers
            max_arg = np.argmax(fitnesses)
            best_individual = population[max_arg]
            best_fitness = fitnesses[max_arg]
            min_arg = np.argmin(fitnesses)
            worst_individual = population[min_arg]
            worst_fitness = fitnesses[min_arg]

            best_performers.append((best_individual, best_fitness))
            worst_performers.append((worst_individual, worst_fitness))
            table.add_row([generation + 1, *best_individual, best_fitness])

            # Selection
            population = selection(population, fitnesses)

            # Generate next population
            next_population = []
            for i in range(0, len(population)-1, 2):
                parent1 = population[i]
                parent2 = population[i + 1]

                child1, child2 = crossover(parent1, parent2, fixed)

                next_population.append(mutation(child1, mutation_rate,  fixed))
                next_population.append(mutation(child2, mutation_rate,  fixed))

            # Elitism: preserve best individual
            next_population.append(best_individual)
            population = next_population
            generation +=1
            end_gen = time.time()
            print("[{}] Estimated time left: {:.2f} min".format(dt.datetime.now().time().strftime('%H:%M'),
                                                                (end_gen - start_gen) / 60 * (
                                                                            generations - generation - 1)))
            print("Current best solution found: {}".format(best_individual))
            print("Current best fitness: {}".format(best_fitness))

    # Print summary table
    print(table)

    # Plot fitness evolution
    generations_list = range(1, len(best_performers) + 1)
    best_fitness_values = [fit[1] for fit in best_performers]
    worst_fitness_values = [fit[1] for fit in worst_performers]

    fig, ax = plt.subplots()
    ax.plot(generations_list, best_fitness_values, label='Best Fitness', color='black')
    ax.fill_between(generations_list, worst_fitness_values, best_fitness_values, color = "black", alpha = 0.5)
    ax.set_xlabel('Generation')
    ax.set_ylabel('Fitness')
    ax.set_title('Fitness Over Generations')
    ax.legend()

    # Plot calculated vs experimental levels using best solution
    experimental, spins, bandhead = read_experimental_data(A, nucleus)
    variables = {
        "N": N,
        "chi": best_individual[0],
        "a1": best_individual[1],
        "a2": best_individual[2],
        "a3": best_individual[3],
        "a4": best_individual[4],
        "a5": best_individual[5],
        "a6": best_individual[6],
        "aF": best_individual[7]
    }
    calculated, _ = give_predictions(A, nucleus, variables, cg, bandhead, spins, even)
    plot_levels(calculated, experimental)
    return best_individual

# --- Run the Algorithm-----------------------------------------------------
start = time.time()
print("[{}] Starting generation 1".format(dt.datetime.now().time().strftime('%H:%M')))
best_solution = genetic_algorithm(population_size, generations, mutation_rate, fixed, A, nucleus)
end = time.time()

print("[{}] Finished! Total time elapsed: {:.1f} min".format(dt.datetime.now().time().strftime('%H:%M'),
                                                              (end - start)/60))
print(f"Best solution found: {best_solution}")
# Save results to file
try:
    with open("Fit_results_{}.txt".format(A), "w") as text_file:
        text_file.write("Best solution:\n{}".format(best_solution))
except:
    pass

# Show plots
plt.show()