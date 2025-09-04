import pandas as pd
import numpy as np
import sympy as sp
from functions.Ha import construct_Ha_even, construct_Ha_odd
from functions.Hb import construct_Hb_even, construct_Hb_odd
from functions.Hc import construct_Hc

def construct_and_solve_H_even(matrix_dir, variables, L, cg, min_eigenvalue = False, verbose = False):
    """
    Construct and solve the even Hamiltonian matrix H = Ha + Hb.

    Input:
    matrix_dir (str): Directory where matrix CSV files are stored or should be created.
    variables (dict): Dictionary of symbolic variables and their numerical values.
    L (int): Angular momentum quantum number.
    cg (list): Clebsch-Gordan coefficients
    min_eigenvalue (bool or float): If False, compute ground state energy. Otherwise, use provided value.
    verbose (bool): Print diagnostic messages if True.

    Returns:
    s [tuples]: List of (eigenvalue, eigenvector) pairs with ground state energy subtracted.
    min_eigenvalue (float): Minimum eigenvalue (ground state energy).
    """
    L = int(L)
    N = variables["N"]
    # import matrix elements H
    def read_matrix_elements(matrix_dir, N, L, cg):
        """Read or construct Ha and Hb matrices from CSV files."""
        try:
            if verbose: print("Reading matrix elements: " + matrix_dir + "h{}a.csv".format(L))
            Ha = pd.read_csv(matrix_dir + "h{}a.csv".format(L), header=None)
        except:
            if verbose: print("Couldn't read matrix element, creating element...")
            construct_Ha_even(matrix_dir, N, L, cg)
            Ha = pd.read_csv(matrix_dir + "h{}a.csv".format(L), header=None)
        try:
            if verbose: print("Reading matrix elements: " + matrix_dir + "h{}b.csv".format(L))
            Hb = pd.read_csv(matrix_dir + "h{}b.csv".format(L), header=None)
        except:
            if verbose: print("Couldn't read matrix element, creating element...")
            construct_Hb_even(matrix_dir, N, L, cg)
            Hb = pd.read_csv(matrix_dir + "h{}b.csv".format(L), header=None)
        return Ha, Hb
    Ha, Hb = read_matrix_elements(matrix_dir, N, L, cg)
    def solve_H(Ha, Hb, variables):
        """Solve eigenvalue problem for H = Ha + Hb."""
        # Define variables
        m, chi, a1, a2, a3, a4, a5, a6 = sp.symbols("m chi a1 a2 a3 a4 a5 a6")
        # apply symbolic evaluation and ensure variables are correctly substituted.
        def evaluate_expression(expr):
            """Safely evaluate symbolic expressions with substitution."""
            try:
                if isinstance(expr, str):
                    expr = expr.strip()

                    symbolic_expr = sp.sympify(expr)
                    substituted_expr = symbolic_expr.subs(variables)

                    return substituted_expr.evalf(16)
                else:
                    return expr
            except Exception as e:
                print(f"Exception occurred for expression: {expr}, Error: {e}")
                return expr


        # Attempt to convert to float
        def to_float(x):
            """Convert to float, fallback to NaN."""
            try:
                return float(x)
            except (ValueError, TypeError):
                return float("NaN")


        # Compute the eigenvalues
        def compute_eigensystem(matrix_a, matrix_b):
            """Compute eigenvalues and eigenvectors of H = A + B."""
            matrix_a = matrix_a.applymap(evaluate_expression).applymap(to_float).fillna(0)
            matrix_b = matrix_b.applymap(evaluate_expression).applymap(to_float).fillna(0)
            matrix_sum = matrix_a + matrix_b
            eigenvalues, eigenvectors = np.linalg.eig(matrix_sum)
            return eigenvalues, np.transpose(eigenvectors)
        sev = compute_eigensystem(Ha, Hb)
        return sev

    # Compute the minimum of the eigenvalues
    if verbose: print("Solving Hamiltonian... for L = {}".format(L))
    sev = solve_H(Ha,Hb, variables)
    if min_eigenvalue is False:
        min_eigenvalue = min(sev[0])

    # Subtracts the ground state energy, pairs eigenvalues with eigenvectors
    # and sorts based on the real part of the eigenvalues
    def adjust_and_sort(eigenvalues, eigenvectors, min_value):
        """Shift eigenvalues by ground state energy and sort by real part."""
        adjusted = [(ev - min_value, ev_vec) for ev, ev_vec in zip(eigenvalues, eigenvectors)]

        sorted_adjusted = sorted(adjusted, key=lambda x: np.real(x[0]))

        return sorted_adjusted
    s = adjust_and_sort(sev[0], sev[1], min_eigenvalue)
    return s, min_eigenvalue

def construct_and_solve_H_odd(matrix_dir, variables, J, J_bh, cg, min_eigenvalue = False, verbose = False):
    """
    Construct and solve the odd Hamiltonian matrix H = Ha + Hb + Hc.

    Input:
    matrix_dir (str): Directory where matrix CSV files are stored or should be created.
    variables (dict): Dictionary of symbolic variables and their numerical values.
    J (float): Total angular momentum (can be half-integer).
    J_bh (float): angular momentum of band head (used in odd systems).
    cg (list): Clebsch-Gordan coefficients
    min_eigenvalue (bool or float): If False, compute ground state energy. Otherwise, use provided value.
    verbose (bool): Print diagnostic messages if True.

    Returns:
    s [tuples]: List of (eigenvalue, eigenvector) pairs with ground state energy subtracted.
    min_eigenvalue (float): Minimum eigenvalue (ground state energy).
    """
    N = variables["N"]
    # import matrix elements H
    def read_matrix_elements_odd(matrix_dir, N, J, J_bh, cg):
        """Read or construct Ha, Hb, Hc matrices from CSV files."""
        try:
            if verbose: print("Reading matrix elements: " + matrix_dir + "h{}a.csv".format(int(J*2)))
            Ha = pd.read_csv(matrix_dir + "h{}a.csv".format(int(J*2)), header=None)
        except:
            if verbose: print("Couldn't read matrix element, creating element...")
            construct_Ha_odd(matrix_dir, N, J, cg, J_bh)
            Ha = pd.read_csv(matrix_dir + "h{}a.csv".format(int(J*2)), header=None)
        try:
            if verbose: print("Reading matrix elements: " + matrix_dir + "h{}b.csv".format(int(J*2)))
            Hb = pd.read_csv(matrix_dir + "h{}b.csv".format(int(J*2)), header=None)
        except:
            if verbose: print("Couldn't read matrix element, creating element...")
            construct_Hb_odd(matrix_dir, N, J, cg, J_bh)
            Hb = pd.read_csv(matrix_dir + "h{}b.csv".format(int(J*2)), header=None)
        try:
            if verbose: print("Reading matrix elements: " + matrix_dir + "h{}c.csv".format(int(J*2)))
            Hc = pd.read_csv(matrix_dir + "h{}c.csv".format(int(J*2)), header=None)
        except:
            if verbose: print("Couldn't read matrix element, creating element...")
            construct_Hc(matrix_dir, N, J, cg, J_bh)
            Hc = pd.read_csv(matrix_dir + "h{}c.csv".format(int(J*2)), header=None)
        return Ha, Hb, Hc

    Ha, Hb, Hc = read_matrix_elements_odd(matrix_dir, N, J, J_bh, cg)

    def solve_H(Ha, Hb, Hc, variables):
        """Solve eigenvalue problem for H = Ha + Hb + Hc."""
        # Define variables
        #m, chi, a1, a2, a3, a4, a5, a6, a7 = sp.symbols("m chi a1 a2 a3 a4 a5 a6 a7")
        # apply symbolic evaluation and ensure variables are correctly substituted.
        def evaluate_expression(expr):
            try:
                if isinstance(expr, str):
                    expr = expr.strip()

                    symbolic_expr = sp.sympify(expr)
                    substituted_expr = symbolic_expr.subs(variables)

                    return substituted_expr.evalf(16)
                else:
                    return expr
            except Exception as e:
                print(f"Exception occurred for expression: {expr}, Error: {e}")
                return expr


        # Attempt to convert to float
        def to_float(x):
            try:
                return float(x)
            except (ValueError, TypeError):
                return float("NaN")


        # Compute the eigenvalues
        def compute_eigensystem(matrix_a, matrix_b, matrix_c):
            matrix_a = matrix_a.applymap(evaluate_expression).applymap(to_float).fillna(0)
            matrix_b = matrix_b.applymap(evaluate_expression).applymap(to_float).fillna(0)
            matrix_c = matrix_c.applymap(evaluate_expression).applymap(to_float).fillna(0)
            matrix_sum = matrix_a + matrix_b + matrix_c
            matrix_sum = matrix_sum.fillna(0)
            eigenvalues, eigenvectors = np.linalg.eig(matrix_sum)
            return eigenvalues, np.transpose(eigenvectors)

        sev = compute_eigensystem(Ha, Hb, Hc)
        return sev

    # Compute the minimum of the eigenvalues
    if verbose: print("Solving Hamiltonian... for J = {}/2".format(int(J*2)))
    sev = solve_H(Ha,Hb,Hc, variables)
    if min_eigenvalue is False:
        min_eigenvalue = min(sev[0])

    # Subtracts the ground state energy, pairs eigenvalues with eigenvectors
    # and sorts based on the real part of the eigenvalues
    def adjust_and_sort(eigenvalues, eigenvectors, min_value):
        adjusted = [(ev - min_value, ev_vec) for ev, ev_vec in zip(eigenvalues, eigenvectors)]

        sorted_adjusted = sorted(adjusted, key=lambda x: np.real(x[0]))

        return sorted_adjusted
    s = adjust_and_sort(sev[0], sev[1], min_eigenvalue)
    return s, min_eigenvalue