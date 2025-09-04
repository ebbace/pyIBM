import matplotlib.pyplot as plt
import numpy as np
from functions.E import give_predictions, calculate_rmse
from functions.read import read_experimental_data, read_variables, read_cg_coefficients
from functions.BE2 import calculate_B

"""
IBM Contour Plotting Script

This script generates a contour plot of a chosen nuclear property (Z-function) 
for a given nucleus in the Interacting Boson Model (IBM). The user can vary two 
Hamiltonian parameters and visualize their effect on energy levels, level ratios, 
RMSE, or B(E2) ratios.
"""


# --- Settings ----------------------------------------
nucleus = "Te"
A = 109

#Set the two parameters you want to be free to a list [min, max], set the rest to fixed values
parameters = {"chi": 0.40,
              "a1": -554.35189,
              "a2": 0,
              "a3": [0, 250],
              "a4": [0, -20],
              "a5": 0,
              "a6": 0,
              "aF": 437.24}

# Resolution of contour grid. Higher = finer but slower.
dim = 20

# --- Code begins here -------------------------------
# --- Variables --------------------------------------
variables = read_variables(A, nucleus) # Load nuclear parameters from data file
N = variables["N"] # Boson number

# ---------------------------------------------------
# check if even
if A % 2 == 0:
    iseven = True
else:
    iseven = False

# ---------------------------------------------------
# Load Clebsch-Gordan coefficients
cg = read_cg_coefficients(N)

# --- Helper functions for Z-axis calculations -----------------------------

def ask_for_input():
    """Prompt user to select the Z-function to plot in the contour graph."""
    z_function = input("Select what to plot:"
                       "\n(1) Energy level energy "
                       "\n(2) Ratio of energy levels "
                       "\n(3) RMSE between calculated and experimental energy levels"
                       "\n(4) B(E2) ratio\n")
    return z_function


def energy_level(A, nucleus, spin, order, variables, cg, bh, even):
    """
    Compute the energy of a given nuclear level.

    Args:
    A (int): Mass number
    nucleus (str): Element
    spin (float): Spin of the level.
    order (int): Order of the level (1 = yrast).
    variables (dict): parameter values
    cg: Clebsch-Gordan coefficients.
    bh (float): band head spin
    even (bool): Whether A is even.

    Returns:
    level (float): Energy of the level.
    title (str): Label for plotting.
    """
    spin = int(float(spin) - bh)
    eigenvalues, eigenvectors = give_predictions(A, nucleus, variables, cg, bh, [spin], even)
    if even:
        label = "{}_{}".format(spin, order)
    else:
        label = "{}/2_{}".format(int((spin + bh) * 2), order)
    level = eigenvalues[label]
    title = "E({})".format(label)
    return level, title


def energy_level_ratio(A, nucleus, spins, orders, variables, cg, bh, even):
    """
    Compute the ratio of two energy levels.

    Args:
    spins [float]: List of spins of the levels.
    orders [int]: List of orders of the levels (1 = yrast).

    Returns:
    ratio (float): Ratio between the energy levels.
    title (str): Label for plotting.
    """
    spins = [int(float(s) - bh) for s in spins]
    eigenvalues, eigenvectors = give_predictions(A, nucleus, variables, cg, bh, spins, even)
    if even:
        title = "E({}_{})/E({}_{})".format(spins[0], orders[0], spins[1], orders[1])
        ratio = eigenvalues["{}_{}".format(spins[0], orders[0])] / eigenvalues["{}_{}".format(spins[1], orders[1])]
    else:
        title = "E({}/2_{})/E({}/2_{})".format(int(spins[0]*2 + bh*2), orders[0], int(spins[1]*2 + bh*2), orders[1])
        ratio = eigenvalues["{}/2_{}".format(int(spins[0]*2 + bh*2), orders[0])] / \
                eigenvalues["{}/2_{}".format(int(spins[1]*2 + bh*2), orders[1])]
    return ratio, title


def BE2_ratio(A, nucleus, spins, orders, variables, cg, bandhead, even):
    """
    Compute the ratio of B(E2) transition probabilities.

    Args:
    A (int): Mass number.
    nucleus (str): Nucleus symbol (e.g., "Te").
    spins [float]: List of spins of the levels.
    orders [int]: List of orders of the levels (1 = yrast).
    variables (dict): parameter values
    cg: Clebsch-Gordan coefficients.
    bandhead (float): band head spin
    even (bool): Whether A is even.

    Returns:
    B_ratio (float): ratio between B values
    title (str): Title of Z-function
    """
    spins = [int(float(s) - bandhead) for s in spins]
    B_values = np.zeros(2)
    labels = ["", ""]
    eigenvalues, eigenvectors = give_predictions(A, nucleus, variables, cg, bandhead, spins, even)

    for i, L in enumerate(spins):
        if even:
            initial = "{}_{}".format(int(L), orders[i])
            final = "{}_{}".format(int(L-2), orders[i])
            labels[i] = "{} -> {}".format(initial, final)
        else:
            initial = "{}/2_{}".format(int((L + bandhead) * 2), orders[i])
            final = "{}/2_{}".format(int((L - 2 + bandhead) * 2), orders[i])
            labels[i] = "{} -> {}".format(initial, final)
        B = calculate_B(eigenvectors[initial], eigenvectors[final], int(L), int(L - 2),
                        variables["N"], cg, variables["chi"])
        B_values[i] = B

    title = "B({})/B({})".format(labels[1], labels[0])
    B_ratio = B_values[1]/B_values[0]
    return B_ratio, title


def rmse(A, nucleus, variables, cg, even):
    """
    Compute the negative RMSE between experimental and predicted levels.
    Used for fitting (lower RMSE = better).

    Args:
    A (int): Mass number.
    nucleus (str): Nucleus symbol (e.g., "Te").
    variables (dict): parameter values
    cg: Clebsch-Gordan coefficients.
    even (bool): Whether A is even.

    Returns:
    rmse (float): Fitness value, positive
    title (str): Title of Z-function
    """
    title = "RMSE, {}-{}".format(nucleus, A)
    targets, spins, bandhead = read_experimental_data(A, nucleus)
    predictions, _ = give_predictions(A, nucleus, variables, cg, bandhead, spins, even)
    rmse = -calculate_rmse(predictions, targets)
    return rmse, title


# --- Main function ----------------------------


def contour(nucleus, A, N, cg, even, parameters, dim):
    """
    Generate a contour plot for two free parameters in the IBM Hamiltonian.

    Args:
    nucleus (str): Nucleus symbol (e.g., "Te").
    A (int): Mass number.
    N (int): Boson number.
    cg: Clebsch-Gordan coefficients.
    even (bool): Whether A is even.
    parameters (dict): Fixed and free parameters.
    dim (int): Resolution of contour grid.
    """
    Z = np.zeros([dim,dim])  # Z-values for contour

    free_keys = np.array([])  # Names of free parameters
    fixed_params = {}  # Store fixed parameters
    ranges = []  # Parameter ranges for mesh

    # ---- Choose Z function-------
    selected = False
    while selected == False:
        z_function = ask_for_input()
        if z_function == "1":
            spin = input("Spin of level (for non-integer spin, write in decimals: 11/2 = 5.5): ")
            order = int(input("Order of level (1 = yrast):"))
            selected = True
        elif z_function == "2":
            spin1 = input("Spin of numerator (for non-integer spin, write in decimals: 11/2 = 5.5): ")
            order1 = input("Order of numerator (1 = yrast):")
            spin2 = input("Spin of denominator (for non-integer spin, write in decimals: 11/2 = 5.5): ")
            order2 = input("Order of denominator (1 = yrast):")
            spin = [spin1, spin2]
            order = [int(order1), int(order2)]
            selected = True
        elif z_function == "3":
            selected = True
        elif z_function == "4":
            print("Calculating B(E2; K_i -> K_i - 2) / B(E2; J_j -> J_j - 2)")
            spin1 = input("Spin K (for non-integer spin, write in decimals: 11/2 = 5.5): ")
            order1 = input("Order i (1 = yrast):")
            spin2 = input("Spin J (for non-integer spin, write in decimals: 11/2 = 5.5): ")
            order2 = input("Order j (1 = yrast):")
            spins = [spin2, spin1]
            orders = [int(order2), int(order1)]
            selected = True
        else:
            print("Not a valid input")

    bandhead = float(input("Spin of bandhead (for non-integer spin, write in decimals: 11/2 = 5.5): "))

    print("\nConstructing contour plot, this could take some time...")

    # Identify free vs fixed parameters
    for p in parameters.keys():
        if isinstance(parameters[p], list):
            free_keys = np.append(free_keys, p)
            ranges.append(np.linspace(parameters[p][0], parameters[p][1], dim))
        else:
            fixed_params[p] = parameters[p]
    curr_params = parameters
    for i, a1 in enumerate(ranges[0]):
        curr_params[free_keys[0]] = a1
        for j, a2 in enumerate(ranges[1]):
            curr_params[free_keys[1]] = a2
            variables = {"N": N,
                         "chi": curr_params["chi"],
                         "a1": curr_params["a1"],
                         "a2": curr_params["a2"],
                         "a3": curr_params["a3"],
                         "a4": curr_params["a4"],
                         "a5": curr_params["a5"],
                         "a6": curr_params["a6"],
                         "aF": curr_params["aF"]}
            if z_function == "1":
                z, title = energy_level(A, nucleus, spin, order, variables, cg, bandhead, even)
                lvls = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
            elif z_function=="2":
                z, title = energy_level_ratio(A, nucleus, spin, order, variables,
                                               cg, bandhead, even)
                lvls = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
            elif z_function =="3":
                z, title = rmse(A, nucleus, variables, cg, even)
                lvls = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
            elif z_function =="4":
                z, title = BE2_ratio(A, nucleus, spins, orders, variables, cg, bandhead, even)
                lvls = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]

            Z[j][i] = z

    # --- Plotting ---
    X, Y = np.meshgrid(ranges[0], ranges[1])
    plt.contourf(X, Y, Z, cmap='rainbow', levels = lvls, extend="both")
    plt.xlabel(free_keys[0])
    plt.ylabel(free_keys[1])
    plt.colorbar()
    plt.title("{}-{}\n{}\n{}".format(nucleus, A, title, fixed_params), fontsize=12)
    plt.tight_layout()
    fig = plt.gcf()
    plt.show()
    fig.savefig("Contour_plot_{}{}_{}_{}_{}.png".format(A, nucleus, z_function, free_keys[0], free_keys[1]))

# --- Run script -------------------------------------------------
contour(nucleus, A, N, cg, iseven, parameters, dim)
