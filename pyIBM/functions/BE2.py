import numpy as np
from math import sqrt
from sympy.physics.wigner import wigner_6j
from functions.definitions import integer_part, kronecker_delta, sign


def calculate_B(svi, svf, Li, Lf, m, cg, x):
    """
    Compute the B transition probability between two states using Clebsch–Gordan coefficients,
    basis construction, and Wigner 6-j symbols.

    Args:
    svi (List[float]): State vector for initial state basis.
    svf (List[float]): State vector for final state basis.
    Li (int): Angular momentum of the initial state.
    Lf (int): Angular momentum of the final state.
    m (int): Boson number
    cg (List[Any]): List of Clebsch–Gordan coefficient tables.
    x (float): chi.


    Returns:
    float: The computed B transition probability.
    """
    cg1, cg2, cg3, = cg[0], cg[1], cg[2]
    ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1) for u in range(0, 2 * m - 6 * b + 1, 2) if
          2 * m - 6 * b - 2 * u >= 0]

    def basis(L):
        """
        Construct basis states for a given angular momentum L.

        Args:
        L (int): Angular momentum.

        Returns:
        List[tuple]: Basis states as (r, u, s, L).
        """
        basis = [(r[0], r[1], s, L)
              for r in ru
              for s in range(1, integer_part((r[0] + r[1] + 2 - L) / 2) *
                             kronecker_delta(sign((r[0] + r[1] + 2 - L) / 2), 1) -
                             integer_part((r[0] + 1 - L) / 2) *
                             kronecker_delta(sign((r[0] + 1 - L) / 2), 1) -
                             integer_part((r[1] + 1 - L) / 2) *
                             kronecker_delta(sign((r[1] + 1 - L) / 2), 1) + 1)]
        return basis

    ibasis = basis(Li)
    fbasis = basis(Lf)

    cg1L = [item for item in cg1 if item[11] == Lf]  # Using zero-based index
    cg4L = [item for item in cg2 if item[3] == Li]
    cg5L = [item for item in cg3 if (item[3] == Li and item[11] == Lf)]

    aa = []

    # Functions for basis1
    def ra(h):
        return fbasis[h][0]  # Using zero-based indexing


    def ua(h):
        return fbasis[h][1]


    def ka(h):
        return fbasis[h][2]


    def La(h):
        return fbasis[h][3]

    # Functions for basis2
    def rb(j):
        return ibasis[j][0]


    def ub(j):
        return ibasis[j][1]


    def kb(j):
        return ibasis[j][2]


    def Lb(j):
        return ibasis[j][3]


    def rbb(jb):
        return ibasis[jb][0]


    def ubb(jb):
        return ibasis[jb][1]


    def kbb(jb):
        return ibasis[jb][2]


    def Lbb(jb):
        return ibasis[jb][3]



    for h in range(len(fbasis)):
        for j in range(len(ibasis)):
            # Find all matching entries for the current (h, j)
            entries = [entry for entry in cg5L
                       if (entry[0] == rb(j) and entry[1] == ub(j) and entry[2] == kb(j) and entry[3] == Lb(j) and
                           entry[8] == ra(h) and entry[9] == ua(h) and entry[10] == ka(h) and entry[11] == La(h))]

            if entries:
                aa.extend(entries)  # Add the found entries
            else:
                aa.append(0)
    aa1 = []
    for s in range(len(aa)):
        if isinstance(aa[s], list) and len(aa[s]) > 0:  # Check if aa[s] is a list and not empty
            value = aa[s][13] * (sqrt(2 * aa[s][11] + 1) * 3 / (2 * sqrt(6)) *
                                 (-1) ** 1 *
                                 # (-1) ** (1 + (1 if aa[s][1] == 0 else 0)) *  # alternatively: KroneckerDelta
                                 sqrt(4 / 3 * (aa[s][0] ** 2 + aa[s][1] ** 2 +
                                               3 * aa[s][0] + 3 * aa[s][1] +
                                               aa[s][0] * aa[s][1])))
            aa1.append(value)  # Append the calculated value
        else:
            aa1.append(0)

    ff = [
        [
            [entry for entry in cg4L
             if entry[0] == rb(j) and entry[1] == ub(j) and entry[2] == kb(j) and entry[3] == Lb(j)],
            [entry for entry in cg1L
             if entry[8] == ra(h) and entry[9] == ua(h) and entry[10] == ka(h) and entry[11] == La(h)]
        ]
        for h in range(len(fbasis)) for j in range(len(ibasis))
    ]

    ff1 = []

    for s in range(len(ff)):
        sum_value = 0  # Accumulator for the sum for each s

        # Loop over elements in ff[s][1]
        for k in range(len(ff[s][0])):

            for kp in range(len(ff[s][1])):

                # Check the conditions
                if (ff[s][0][k][8] == ff[s][1][kp][0] and
                        ff[s][0][k][9] == ff[s][1][kp][1] and
                        ff[s][0][k][10] == ff[s][1][kp][2] and
                        ff[s][0][k][11] == ff[s][1][kp][3]):

                    term1 = sqrt(5) * (-1) ** (ff[s][1][kp][11] + ff[s][0][k][3])
                    term2 = (wigner_6j(2, 2, 2, ff[s][0][k][3], ff[s][1][kp][11], ff[s][1][kp][3]))
                    term3 = sqrt((2 * ff[s][1][kp][11] + 1) * (2 * ff[s][1][kp][3] + 1))
                    term4 = ff[s][1][kp][13] * ff[s][0][k][13]

                    # Calculate the inner product of all terms
                    combined_term = term1 * term2 * term3 * term4

                    # Nested conditions for the second term involving different comparisons
                    condition1 = (ff[s][1][kp][8] == (ff[s][1][kp][0] + 2) and
                                  ff[s][1][kp][9] == ff[s][1][kp][1])
                    term5 = (sqrt((2 * (m - 1) + 2 * ff[s][1][kp][0] + ff[s][1][kp][1] + 12) *
                                  (ff[s][1][kp][0] + 2) * (ff[s][1][kp][0] + ff[s][1][kp][1] + 3) /
                                  (6 * (ff[s][1][kp][8] + 1) * (ff[s][1][kp][8] + ff[s][1][kp][9] + 2)))
                             if condition1 else 0)

                    condition2 = (ff[s][1][kp][8] == (ff[s][1][kp][0] - 2) and
                                  ff[s][1][kp][9] == (ff[s][1][kp][1] + 2))
                    term6 = (sqrt((2 * (m - 1) - ff[s][1][kp][0] + ff[s][1][kp][1] + 9) *
                                  ff[s][1][kp][0] * (ff[s][1][kp][1] + 2) /
                                  (6 * (ff[s][1][kp][8] + 1) * (ff[s][1][kp][9] + 1)))
                             if condition2 else 0)

                    condition3 = (ff[s][1][kp][8] == ff[s][1][kp][0] and
                                  ff[s][1][kp][9] == (ff[s][1][kp][1] - 2))
                    term7 = (sqrt((2 * (m - 1) - ff[s][1][kp][0] - 2 * ff[s][1][kp][1] + 6) *
                                  ff[s][1][kp][1] * (ff[s][1][kp][0] + ff[s][1][kp][1] + 1) /
                                  (6 * (ff[s][1][kp][8] + ff[s][1][kp][9] + 2) * (ff[s][1][kp][9] + 1)))
                             if condition3 else 0)

                    final_term = term5 + term6 + term7
                    if (ff[s][0][k][0] == ff[s][0][k][8] + 2 and
                            ff[s][0][k][1] == ff[s][0][k][9]):
                        term8 = sqrt(
                            (2 * (m - 1) + 2 * ff[s][0][k][8] + ff[s][0][k][9] + 12) *
                            (ff[s][0][k][8] + 2) *
                            (ff[s][0][k][8] + ff[s][0][k][9] + 3) /
                            (6 * (ff[s][0][k][0] + 1) * (ff[s][0][k][0] + ff[s][0][k][1] + 2))
                        )
                    else:
                        term8 = 0
                    if (ff[s][0][k][0] == ff[s][0][k][8] - 2 and
                            ff[s][0][k][1] == ff[s][0][k][9] + 2):
                        term9 = sqrt(
                            (2 * (m - 1) - ff[s][0][k][8] + ff[s][0][k][9] + 9) *
                            ff[s][0][k][8] *
                            (ff[s][0][k][9] + 2) /
                            (6 * (ff[s][0][k][0] + 1) * (ff[s][0][k][1] + 1))
                        )
                    else:
                        term9 = 0
                    if (ff[s][0][k][0] == ff[s][0][k][8] and
                            ff[s][0][k][1] == ff[s][0][k][9] - 2):
                        term10 = sqrt(
                            (2 * (m - 1) - ff[s][0][k][8] - 2 * ff[s][0][k][9] + 6) *
                            ff[s][0][k][9] *
                            (ff[s][0][k][8] + ff[s][0][k][9] + 1) /
                            (6 * (ff[s][0][k][0] + ff[s][0][k][1] + 2) * (ff[s][0][k][1] + 1))
                        )
                    else:
                        term10 = 0
                    add_term = (term8 + term9 + term10) * sqrt(
                        (ff[s][0][k][0] + 1) *
                        (ff[s][0][k][1] + 1) *
                        (ff[s][0][k][0] + ff[s][0][k][1] + 2) /
                        ((ff[s][0][k][8] + 1) * (ff[s][0][k][9] + 1) * (ff[s][0][k][8] + ff[s][0][k][9] + 2))
                    )

                    # Add the product of combined_term and final_term to the accumulator
                    sum_value += combined_term * final_term * add_term

        # Append the accumulated sum to the ff1 list
        ff1.append(float(sum_value))

    length_ibasis = len(ibasis)
    length_fbasis = len(fbasis)

    # Initialize B
    B = 0

    # Loop over h and j for summation
    for h in range(length_fbasis):
        for j in range(length_ibasis):
            index_aa1 = j + (h) * length_ibasis
            index_ff1 = j + (h) * length_ibasis
            B += svi[j] * svf[h] * (aa1[index_aa1] + (x + np.sqrt(7) / 2.0) * ff1[index_ff1])

    # Final calculation of B
    B = (B) ** 2 / (2 * Li + 1)
    return B