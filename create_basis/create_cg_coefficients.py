import numpy as np
import pandas as pd
import subprocess
import shlex
import os
from io import StringIO

basis_dir = "./basis/"  # Path to the CG coefficients library (list1, list2, list3 directories)
cg_code_dir = "./Su3cgvcs/"  # path to the cgu3o3drv code, used for calculating CG coeffs if not already in library.
verbose = True

def read_and_create_cg_coefficients(N, basis_dir, cg_code_dir, verbose=False):
    files = ["{}list1/N={}.txt".format(basis_dir, N),
             "{}list2/N={}.txt".format(basis_dir, N),
             "{}list3/N={}.txt".format(basis_dir, N),
             "{}list4/N={}.txt".format(basis_dir, N),
             "{}list5/N={}.txt".format(basis_dir, N),
             "{}list6/N={}.txt".format(basis_dir, N)
             ]

    #Read cg coefficients, create the lists if they don't exist in the database

    try:
        if verbose: print("Reading cg file: " + files[0])
        cg1 = read_and_format_cg(files[0])

    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=1)
        cg1 = read_and_format_cg(files[0])

    try:
        if verbose: print("Reading cg file: " + files[1])
        cg2 = read_and_format_cg(files[1])
    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=2)
        cg2 = read_and_format_cg(files[1])
    try:
        if verbose: print("Reading cg file: " + files[2])
        cg3 = read_and_format_cg(files[2])
    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=3)
        cg3 = read_and_format_cg(files[2])
    try:
        if verbose: print("Reading cg file: " + files[3])
        cg4 = read_and_format_cg(files[3])
    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=4)
        cg4 = read_and_format_cg(files[3])
    try:
        if verbose: print("Reading cg file: " + files[4])
        cg5 = read_and_format_cg(files[4])
    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=5)
        cg5 = read_and_format_cg(files[4])
    try:
        if verbose: print("Reading cg file: " + files[5])
        cg6 = read_and_format_cg(files[5])
    except:
        if verbose: print("Cg file could not be read, creating missing cg file...")
        create_CG_database(basis_dir, cg_code_dir, N, list=6)
        cg6 = read_and_format_cg(files[5])

    return [cg1, cg2, cg3, pd.DataFrame(cg4), pd.DataFrame(cg5), pd.DataFrame(cg6)]


def generate_basis(N):
    basis = pd.DataFrame(columns=["lmb", "mu", "L"])
    lmb_array = np.array([])
    mu_array = np.array([])
    i = 0
    while 2 * N - 6 * i >= 0:
        j = 0
        while 2 * N - 6 * i - 4 * j >= 0:
            lmb = 2 * N - 6 * i - 4 * j
            mu = 2 * j
            lmb_array = np.append(lmb_array, lmb)
            mu_array = np.append(mu_array, mu)
            j += 1
        i+=1


    # Degenerate lambda and mu based on K and L
    lambda_expanded = np.array([])
    mu_expanded = np.array([])
    L_expanded = np.array([])

    for i, lmb in enumerate(lmb_array):
        max_val = int(lmb + mu_array[i])
        for L in range(0, max_val + 1):
            lambda_expanded = np.append(lambda_expanded, lmb)
            mu_expanded = np.append(mu_expanded, mu_array[i])
            L_expanded = np.append(L_expanded, L)
    basis["lmb"] = lambda_expanded
    basis["mu"] = mu_expanded
    basis["L"] = L_expanded
    return basis


def construct_basis(N, list):
    basis = pd.DataFrame(columns=["lm1", "mu1", "L1", "lm2", "mu2", "L2", "lm3", "mu3", "L3"])
    if list == 1:
        basis_N3 = generate_basis(N)
        basis_N2 = generate_basis(1)
        basis_N1 = generate_basis(N-1)
    if list == 2:
        basis_N3 = generate_basis(N-1)
        basis_N2 = generate_basis(2)
        basis_N1 = generate_basis(N)
    if list == 3:
        basis_N3 = generate_basis(N)
        basis_N2 = pd.DataFrame({"lmb": [1], "mu": [1], "L": [2]})
        basis_N1 = basis_N3
    if list == 4:
        basis_N3 = generate_basis(N)
        basis_N2 = generate_basis(1)
        basis_N1 = generate_basis(N - 1)
    if list == 5:
        basis_N3 = generate_basis(N - 1)
        basis_N2 = generate_basis(1)
        basis_N1 = generate_basis(N - 2)
    if list == 6:
        basis_N3 = generate_basis(N - 2)
        basis_N2 = generate_basis(1)
        basis_N1 = generate_basis(N - 3)
    row = 0
    for i in basis_N1.index:
        for j in [basis_N2.index[-1]]:
            for k in basis_N3.index:
                basis.loc[row] = np.array([basis_N1["lmb"][i], basis_N1["mu"][i], basis_N1["L"][i],
                                         basis_N2["lmb"][j], basis_N2["mu"][j], basis_N2["L"][j],
                                         basis_N3["lmb"][k], basis_N3["mu"][k], basis_N3["L"][k]])
                row +=1
    basis = basis.sort_values(["lm1", "lm3", "L2"], kind="mergesort", ascending = [False, False, True])
    return basis


def create_CG_database(basis_dir, cg_code_dir, N, list):
    coupled_basis = construct_basis(N, list = list)
    first_time = True
    for i in coupled_basis.index:
        lm1 = int(coupled_basis["lm1"][i])
        mu1 = int(coupled_basis["mu1"][i])
        L1 = int(coupled_basis["L1"][i])
        lm2 = int(coupled_basis["lm2"][i])
        mu2 = int(coupled_basis["mu2"][i])
        L2 = int(coupled_basis["L2"][i])
        lm3 = int(coupled_basis["lm3"][i])
        mu3 = int(coupled_basis["mu3"][i])
        L3 = int(coupled_basis["L3"][i])
        s = "{} {} {} {} {} {}\n{} {} {}".format(lm1, mu1, lm2, mu2, lm3, mu3, L1, L2, L3)

        with open(cg_code_dir + "inso3.txt", "w") as f:
            f.write(s)

        proc = subprocess.Popen("{}cgu3o3drv < {}inso3.txt".format(cg_code_dir, cg_code_dir), stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        output = out.decode("utf-8").split("::\n\n")[1]
        output = output.split(" WARNING")[0]
        if output.split("\n")[1].split()[-1] == "exist":
            pass
        else:
            if first_time: # Open file and overwrite first time
                filename="{}list{}/N={}.txt".format(basis_dir, list, N)
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                with open(filename, "w") as f:
                    f.write(output)
                first_time = False
            else: # Open file and append
                with open("{}list{}/N={}.txt".format(basis_dir, list, N), "a") as f:
                    f.write(output.split('\n', 1)[1])


def read_and_format_cg(path):
    # The order of the parameters in the file should be lm1 mu1 k1 L1 lm2 mu2 k2 L2 rh lm3 mu3 k3 L3
    # Read the file line by line
    with open(path, "r") as file:
        lines = file.readlines()

    # Split each line by whitespace and convert it into a list of lists
    data = [line.split() for line in lines]

    # Convert to DataFrame, filling missing values with NaN to handle rows of different lengths
    t1 = pd.DataFrame(data)

    # Convert all elements to numeric where possible, keeping others as NaN
    t1 = t1.apply(pd.to_numeric, errors='coerce')

    t2 = t1[t1.apply(lambda x: 8 <= x.count() <= 14, axis=1)]
    data = []
    for row in t2.values:
        row_list = row.tolist()

        # Check the number of non-NaN elements
        non_nan_count = sum(pd.notna(row_list))

        if non_nan_count == 14:
            new_row = [row_list[0], row_list[1], row_list[2], row_list[3], row_list[4], row_list[5], row_list[6],
                       row_list[7], row_list[9], row_list[10], row_list[11], row_list[12], row_list[8],
                       row_list[13]]
            data.append(new_row)
        else:
            # Pad with 'a' in the positions mentioned.
            new_row = ['a', 'a', row_list[0], row_list[1], 'a', 'a', row_list[2], row_list[3],
                       'a', 'a', row_list[5], row_list[6], row_list[4], row_list[7]]
            data.append(new_row)


    # Convert to DataFrame
    data = pd.DataFrame(data)
    # This is to remove the padding and replace with the mu and lambda belonging to that group
    g1 = np.where(t2.apply(lambda x: len(x.dropna()) == 14, axis=1))[0] + 1
    g2 = [i - 1 for i in g1 if i > 1]
    g3 = [0] + g2 + [len(t2)]
    g4 = [data.iloc[g3[k]:g3[k + 1]] for k in range(len(g3) - 1)]

    cg = []
    for group in g4:
        first_row = group.iloc[0]
        for _, row in group.iterrows():
            new_row = [
                first_row[0], first_row[1], row[2], row[3], first_row[4], first_row[5],
                row[6], row[7], first_row[8], first_row[9], row[10], row[11], row[12], row[13]
            ]
            cg.append(new_row)

    return cg

N = int(input("Boson number: "))
read_and_create_cg_coefficients(N, basis_dir, cg_code_dir, verbose)
