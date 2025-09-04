import pandas as pd
import numpy as np


def read_paths(A=None, nucleus=None):
    """
    Generate file paths for basis, matrix elements, variables, and experimental data.

    Args:
    A (int | None): Mass number of the nucleus (optional, only needed for matrix/variable/experimental paths).
    nucleus (str | None): Name of the nucleus (optional, only needed for matrix/variable/experimental paths).

    Returns:
    paths (dict): Dictionary containing directory paths for required data.
    """
    # ---Edit these if necessary------------------------------
    basis_dir = "basis/"  # Path to the CG coefficients library (list1, list2, list3 directories)
    matrix_dir = "matrix_elements_library/{}/{}/".format(nucleus, A)  # Path to the matrix elements (Ha, Hb, Hc files)
    variables_dir = "variables_library/{}/".format(nucleus)  # Path to the variables files
    experimental_data_dir = "experimental_data_library/{}/".format(nucleus)  # Path to the experimental data

    # ---------------------------------------------------------

    paths = {"basis": basis_dir,
             "matrix elements": matrix_dir,
             "variables": variables_dir,
             "experimental data": experimental_data_dir}

    return paths


def read_cg_coefficients(N, verbose=False):
    """
    Read Clebsch-Gordan coefficients from precomputed text files and format them.

    Args:
    N (int): Boson number.
    verbose (bool): If True, print which files are being read.

    Returns:
    List[Any]: List containing formatted CG coefficients as lists and DataFrames.
    """
    paths = read_paths()
    basis_dir = paths["basis"]

    files = ["{}list1/N={}.txt".format(basis_dir, N),
             "{}list2/N={}.txt".format(basis_dir, N),
             "{}list3/N={}.txt".format(basis_dir, N),
             "{}list4/N={}.txt".format(basis_dir, N),
             "{}list5/N={}.txt".format(basis_dir, N),
             "{}list6/N={}.txt".format(basis_dir, N)
             ]

    def read_and_format_cg(path):
        """
        Read and format a CG coefficient file.

        Args:
        path (str): Path to the CG coefficient file.

        Returns:
        List[List[Any]]: Formatted CG coefficients.
        """
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

    # Read cg coefficients, create the lists if they don't exist in the database
    if verbose: print("Reading cg file: " + files[0])
    cg1 = read_and_format_cg(files[0])

    if verbose: print("Reading cg file: " + files[1])
    cg2 = read_and_format_cg(files[1])

    if verbose: print("Reading cg file: " + files[2])
    cg3 = read_and_format_cg(files[2])

    if verbose: print("Reading cg file: " + files[3])
    cg4 = read_and_format_cg(files[3])

    if verbose: print("Reading cg file: " + files[4])
    cg5 = read_and_format_cg(files[4])

    if verbose: print("Reading cg file: " + files[5])
    cg6 = read_and_format_cg(files[5])

    return [cg1, cg2, cg3, pd.DataFrame(cg4), pd.DataFrame(cg5), pd.DataFrame(cg6)]


def read_variables(A, nucleus):
    """
    Read nuclear variables from a CSV file.

    Args:
    A (int): Mass number.
    nucleus (str): Nucleus name.

    Returns:
    Dict[str, Any]: Dictionary of variables.
    """
    paths = read_paths(A, nucleus)
    variables_dir = paths["variables"]
    try:
        variables = pd.read_csv(variables_dir + "{}.csv".format(A)).to_dict("records")[0]
    except:
        print("Variables could not be read, check path:\n{}")
    return variables



def read_experimental_data(A, nucleus):
    """
    Read experimental nuclear level data from CSV file.

    Args:
    A (int): Mass number.
    nucleus (str): Nucleus name.

    Returns:
    Tuple[Dict[str, float], List[int], float]:
    - Dict of energy levels relative to the bandhead.
    - List of normalised spin values.
    - Band head spin value.
    """
    paths = read_paths(A, nucleus)
    exp_data_path = paths["experimental data"]
    levels_dict = pd.read_csv(exp_data_path + "{}.csv".format(A)).to_dict("records")[0]
    N_dict = {}
    spins = []
    spins_normalised = []
    for k in levels_dict.keys():
        spin = int(k.split("_")[0].split("/")[0])
        try:
            denominator = int(k.split("_")[0].split("/")[1])
        except:
            denominator = 1
        spins.append(spin/denominator)
    bandhead = min(spins)
    for spin in spins:
        if spin != bandhead:
            spin_norm = int(spin - bandhead)
            spins_normalised.append(spin_norm)
    spins = list(set(spins_normalised)) #remove duplicate spins


    if denominator == 1:
        bandhead_label = "{}_1".format(int(bandhead))
    else:
        bandhead_label = "{}/{}_1".format(int(bandhead*denominator), denominator)
    for key in levels_dict.keys():
        N_dict[key] = levels_dict[key] - levels_dict[bandhead_label]
    return N_dict, spins, bandhead

