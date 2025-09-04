import matplotlib.pyplot as plt
from functions.BE2 import calculate_B
from functions.plot import plot_levels
from functions.E import give_predictions, calculate_rmse
from functions.read import read_experimental_data, read_variables, read_cg_coefficients


"""
Nuclear Energy Level and B(E2) Transition Analysis
--------------------------------------------------
This script computes theoretical nuclear energy levels using the Interacting Boson Model (IBM),
compares them with experimental data, and calculates B(E2) transition probabilities.


Workflow:
1. Read input parameters and model variables.
2. Calculate energy eigenvalues/eigenvectors.
3. Compare with experimental levels and compute RMSE.
4. Print and plot energy levels.
5. Calculate B(E2) transition strengths and ratios.
"""


# ---Settings----------------------------------------
# User input for the nucleus under study
nucleus = input("Nucleus: ") # Example: "Te"
A = int(input("Mass number: ")) # Example: 120


# List of spins above the bandhead spin (i.e. the spin of the states are given by (spins + bandhead spin))
# that will be used in the calculations (spin 0 is always included)
spins = [2, 3, 4, 6, 8, 10]
maxspin = max(spins)

verbose = True  # Set to False to minimize console output

# ---Variables---------------------------------------
variables = read_variables(A, nucleus) # Read IBM model variables
N = variables["N"] # Number of bosons
chi = variables["chi"] # chi

# ---Check if A is even-------------------------------
# Determines if the nucleus has even or odd mass number
if A % 2 == 0:
    even = True
else:
    even = False

# ---Read CG coefficients----------------------------
cg = read_cg_coefficients(N, verbose=verbose)

# ---Calculate Energy levels-------------------------
experimental, _, bandhead = read_experimental_data(A, nucleus) # Experimental levels

# Theoretical predictions: eigenvalues (energies) and eigenvectors (wavefunctions)
eigenvalues, eigenvectors = give_predictions(A, nucleus, variables, cg, bandhead, spins, even, verbose)

# Root mean square error between experimental and predicted levels
rmse = calculate_rmse(eigenvalues, experimental)
print("\nRMSE of energy levels fit: " + str(-rmse))

# Plot calculated vs experimental energy levels
plot_levels(eigenvalues, experimental)

# ---Print Calculated Energy Levels-------------------
# These print functions can be changed to print the desired output. Use the format:
# print(eigenvalues['L_o'])), where L is the spin and o is the order.
# For odd-A nuclei, spins are labeled as fractions, ex "7/2"
print("\nEigenvalues (keV):")

# Print yrast energy levels (order = 1)
for spin in spins:
    if not even:
        # For odd-A nuclei, spins are half-integers, labeled e.g. "7/2"
        state = str(int((spin + bandhead)*2)) + "/2"
    else:
        state = str(spin)
    try:
        print("{}_1: {:.2f}".format(state, eigenvalues['{}_1'.format(state)]))
    except:
        pass

# Print yrare energy levels for (order = 2)
for spin in [0, 2, 4]:
    if not even:
        state = str(int((spin + bandhead)*2)) + "/2"
    else:
        state = str(spin)
    try:
        print("{}_2: {:.2f}".format(state, eigenvalues['{}_2'.format(state)]))
    except:
        pass

# Print yrare energy levels (order = 3)
for spin in [0]:
    if not even:
        state = str(int((spin + bandhead)*2)) + "/2"
    else:
        state = str(spin)
    try:
        print("{}_3: {:.2f}".format(state, eigenvalues['{}_3'.format(state)]))
    except:
        pass

# ---Calculate and print B(E2) values ------------------------------------
# These print functions can be changed to print the desired output. Use the format:
# B = calculate_B(eigenvectors["Ji_oi"], eigenvectors["Jf_of"], Li, Lf, N, cg, chi)
# print("B({} -> {}): {:.2f}".format("Ji_oi", "Jf_of", B))
# Where Ji and Jf are the initial and final spin of the states (Ex. Ji = 19/2, Jf = 15/2), oi and of are the orders
# of the inital and final states (Ex. 1 for yrast state), and Li and Lf are the angular momentum above the band head
# (Ex. 4 and 2 if the band head is 11/2)

print("\nBE2 values:")
maxorder = 1 # Maximum order of excited bands to consider
B_values = {} # Store calculated B(E2) values
for L in range(0, int(maxspin), 2):
    for order in range(1, maxorder+1):
        if even:
            initial = "{}_{}".format(int(L + 2),  order)
            final = "{}_{}".format(int(L),  order)
        else:
            initial = "{}/2_{}".format(int((L + 2 + bandhead)*2),  order)
            final = "{}/2_{}".format(int((L + bandhead) * 2), order)
        try:
            B = calculate_B(eigenvectors[initial], eigenvectors[final], L+2, L, N, cg, chi)
            B_values[("{} -> {}").format(initial, final)] = B   # Save results
            print("B({} -> {}): {:.2f}".format(initial, final, B))  # Print results
        except:
            pass

# --- calculate and print B(E2) ratios------------------------------------
print("\nB-ratios:")
i = 0
order = 1
for L in range(0, int(maxspin)-2, 2):
    if even:
        initial_top = "{}_{}".format(int(L + 4), order)
        final_top = "{}_{}".format(int(L+2), order)

        initial_bottom = "{}_{}".format(int(L + 2), order)
        final_bottom = "{}_{}".format(int(L), order)
    else:
        initial_top = "{}/2_{}".format(int((L + 4 + bandhead) * 2), order)
        final_top = "{}/2_{}".format(int((L + 2 + bandhead) * 2), order)

        initial_bottom = "{}/2_{}".format(int((L + 2 + bandhead) * 2), order)
        final_bottom = "{}/2_{}".format(int((L + bandhead) * 2), order)
    try:
        B_lower = B_values["{} -> {}".format(initial_bottom, final_bottom)]
        B_higher = B_values["{} -> {}".format(initial_top, final_top)]
        print("B({} -> {})/B({} -> {}): {:.2f}".format(initial_top, final_top, initial_bottom,
                                                           final_bottom, B_higher/B_lower))
    except:
        pass
    i += 1

# Show the level diagram plot
plt.show()
