from functions.H import construct_and_solve_H_even, construct_and_solve_H_odd
import numpy as np
from functions.read import read_paths


def give_predictions(A, nucleus, variables, cg, bh, spins, even, verbose=False):
    """
    Compute eigenvalues and eigenvectors of the Hamiltonian for a given nucleus.

    Args:
    A (int): Mass number of the nucleus.
    nucleus (str): Symbol of the nucleus (e.g., 'Te').
    variables (Dict[str, Any]): Model parameters for the Hamiltonian.
    cg (List[Any]): Clebsch-Gordan coefficients.
    bh (int): Spin of the bandhead.
    spins (List[int]): List of spins to compute beyond the bandhead.
    even (bool): Whether the nucleus has even mass number.
    verbose (bool, optional): Print information. Defaults to False.

    Returns:
    Tuple[Dict[str, float], Dict[str, np.ndarray]]:
    - eigenvalues: Dictionary of state labels to energy values.
    - eigenvectors: Dictionary of state labels to eigenvector arrays.
    """
    paths = read_paths(A, nucleus)
    m_dir = paths["matrix elements"]

    eigenvalues = {}
    eigenvectors = {}
    if even:
        s0, min = construct_and_solve_H_even(m_dir, variables, bh, cg, verbose=verbose)
        for order, p in enumerate(s0):
            lbl = "{}_{}".format(int(bh),  order+1)
            eigenvalues[lbl] = p[0]
            eigenvectors[lbl] = p[1]
    else:
        s0, min = construct_and_solve_H_odd(m_dir, variables, bh, bh, cg, verbose=verbose)
        for order, p in enumerate(s0):
            lbl = "{}/2_{}".format(int((bh)*2), order+1)
            eigenvalues[lbl] = p[0]
            eigenvectors[lbl] = p[1]

    for L in spins:
        if even:
            try:
                s, _ = construct_and_solve_H_even(m_dir, variables, L, cg, min, verbose=verbose)
                for order, p in enumerate(s):
                    lbl = "{}_{}".format(L, order+1)
                    eigenvalues[lbl] = p[0]
                    eigenvectors[lbl] = p[1]
            except:
                if verbose: print("Could not construct matrix element, skipping")
        else:
            try:
                s, _ = construct_and_solve_H_odd(m_dir, variables, L + bh, bh, cg, min,
                                                 verbose=verbose)
                for order, p in enumerate(s):
                    lbl = "{}/2_{}".format(int((L+bh)*2), order+1)
                    eigenvalues[lbl] = p[0]
                    eigenvectors[lbl] = p[1]
            except:
                if verbose: print("Could not construct matrix element, skipping")
    return eigenvalues, eigenvectors


def calculate_rmse(predictions, targets):
    """
    Calculate the root mean square error (RMSE) between predicted and experimental values.

    Args:
    predictions (Dict[str, float]): Predicted energy levels keyed by state label.
    targets (Dict[str, float]): Experimental energy levels keyed by state label.

    Returns:
    float: Negative RMSE (negated so that optimization routines can maximize it).
    """
    diff = np.array([])
    for state in targets.keys():
        try:
            diff = np.append(diff, targets[state] - predictions[state])
        except:
            pass
    rmse = np.sqrt((np.abs(diff) ** 2).mean())
    return -rmse  # Negated to find the maximum


