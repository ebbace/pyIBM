import numpy as np
import pandas as pd
from sympy.physics.wigner import wigner_6j
from functions.definitions import integer_part, kronecker_delta, sign
from math import sqrt
import sympy as sp
import os

# This script calculates the coupling to the odd fermion

def construct_Hc(matrix_dir, m, Jt, cg, jp):
    cg1, cg2, cg3, = cg[0], cg[1], cg[2]

    ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1) for u in range(0, 2 * m - 6 * b + 1, 2) if
          2 * m - 6 * b - 2 * u >= 0]

    basis = [[r[0], r[1], s, a]
              for a in range(int(abs(Jt - jp)), int(Jt + jp + 1)) # plus one due to Python excluding the last number in the range
              for r in ru
              for s in range(1, integer_part((r[0] + r[1] + 2 - a) / 2) *
                             kronecker_delta(sign((r[0] + r[1] + 2 - a) / 2), 1) -
                             integer_part((r[0] + 1 - a) / 2) *
                             kronecker_delta(sign((r[0] + 1 - a) / 2), 1) -
                             integer_part((r[1] + 1 - a) / 2) *
                             kronecker_delta(sign((r[1] + 1 - a) / 2), 1) + 1)]


    cg1L = sum([
        [item for item in cg1 if item[8] == b[0] and item[9] == b[1] and item[10] == b[2] and item[11] == b[3]]
        for b in basis
    ], [])


    cg2L = sum([
        [item for item in cg2 if item[0] == b[0] and item[1] == b[1] and item[2] == b[2] and item[3] == b[3]]
        for b in basis
    ], [])


    cg5L = sum([
        [item for item in cg3 if item[8] == b[0] and item[9] == b[1] and item[10] == b[2] and item[11] == b[3]]
        for b in basis
    ], [])


    # Functions for basis1
    def ra(h):
        return basis[h][0]  # Using zero-based indexing

    def ua(h):
        return basis[h][1]

    def ka(h):
        return basis[h][2]

    def La(h):
        return basis[h][3]

    def rb(j):
        return basis[j][0]

    def ub(j):
        return basis[j][1]

    def kb(j):
        return basis[j][2]

    def Lb(j):
        return basis[j][3]

    def threebar1(s1):
        result = 0

        # Condition 1
        if cg1L[s1][8] == cg1L[s1][0] + 2 and cg1L[s1][9] == cg1L[s1][1]:
            result += np.sqrt(
                (2 * (m - 1) + 2 * cg1L[s1][0] + cg1L[s1][1] + 12) *
                (cg1L[s1][0] + 2) *
                (cg1L[s1][0] + cg1L[s1][1] + 3) /
                (6 * (cg1L[s1][8] + 1) * (cg1L[s1][8] + cg1L[s1][9] + 2))
            )

        # Condition 2
        if cg1L[s1][8] == cg1L[s1][0] - 2 and cg1L[s1][9] == cg1L[s1][1] + 2:
            result += np.sqrt(
                (2 * (m - 1) - cg1L[s1][0] + cg1L[s1][1] + 9) *
                cg1L[s1][0] *
                (cg1L[s1][1] + 2) /
                (6 * (cg1L[s1][8] + 1) * (cg1L[s1][9] + 1))
            )

        # Condition 3
        if cg1L[s1][8] == cg1L[s1][0] and cg1L[s1][9] == cg1L[s1][1] - 2:
            result += np.sqrt(
                (2 * (m - 1) - cg1L[s1][0] - 2 * cg1L[s1][1] + 6) *
                cg1L[s1][1] *
                (cg1L[s1][0] + cg1L[s1][1] + 1) /
                (6 * (cg1L[s1][8] + cg1L[s1][9] + 2) * (cg1L[s1][9] + 1))
            )

        return result

    def threebar2(s2):
        result = 0

        # Condition 1
        if cg2L[s2][0] == cg2L[s2][8] + 2 and cg2L[s2][1] == cg2L[s2][9]:
            result += np.sqrt(
                (2 * (m - 1) + 2 * cg2L[s2][8] + cg2L[s2][9] + 12) *
                (cg2L[s2][8] + 2) *
                (cg2L[s2][8] + cg2L[s2][9] + 3) /
                (6 * (cg2L[s2][0] + 1) * (cg2L[s2][0] + cg2L[s2][1] + 2))
            )

        # Condition 2
        if cg2L[s2][0] == cg2L[s2][8] - 2 and cg2L[s2][1] == cg2L[s2][9] + 2:
            result += np.sqrt(
                (2 * (m - 1) - cg2L[s2][8] + cg2L[s2][9] + 9) *
                cg2L[s2][8] *
                (cg2L[s2][9] + 2) /
                (6 * (cg2L[s2][0] + 1) * (cg2L[s2][1] + 1))
            )

        # Condition 3
        if cg2L[s2][0] == cg2L[s2][8] and cg2L[s2][1] == cg2L[s2][9] - 2:
            result += np.sqrt(
                (2 * (m - 1) - cg2L[s2][8] - 2 * cg2L[s2][9] + 6) *
                cg2L[s2][9] *
                (cg2L[s2][8] + cg2L[s2][9] + 1) /
                (6 * (cg2L[s2][0] + cg2L[s2][1] + 2) * (cg2L[s2][1] + 1))
            )

        return result

    y1 = []
    for h in range(len(basis)):
        row = []
        for j in range(len(basis)):
            common_part = basis[h][:4] + basis[j][:4]
            prefactor = (-sqrt(5)) * (-1) ** (Lb(j) + jp + Jt) * wigner_6j(La(h), jp, Jt, jp, Lb(j), 2)
            sum_value = 0
            for s5 in range(len(cg5L)):
                c = cg5L[s5]
                if (c[0] == rb(j) and c[1] == ub(j) and c[2] == kb(j) and c[3] == Lb(j) and
                        c[8] == ra(h) and c[9] == ua(h) and c[10] == ka(h) and c[11] == La(h) and
                        c[0] == c[8] and c[1] == c[9]):
                    delta = kronecker_delta(c[1], 0)
                    factor = (sqrt(2 * c[11] + 1) * 3 / (2 * sqrt(6)) *
                              (-1) ** (1 + delta) *
                              sqrt(4 / 3 * (c[0] ** 2 + c[1] ** 2 + 3 * c[0] + 3 * c[1] + c[0] * c[1])))
                    sum_value += factor * c[13]
            row.append(common_part + [prefactor * sum_value])
        y1.append(row)

    # y3 translation
    y3 = []
    for h in range(len(basis)):
        row = []
        for j in range(len(basis)):
            common_part = basis[h][:4] + basis[j][:4]
            condition = (abs(basis[h][0] - basis[j][0]) <= 4 and
                         abs(basis[h][1] - basis[j][1]) <= 4 and
                         abs(basis[h][3] - basis[j][3]) <= 2)
            if condition:
                sum_value = 0
                for s1 in range(len(cg1L)):
                    for s2 in range(len(cg2L)):
                        c1 = cg1L[s1]
                        c2 = cg2L[s2]
                        if (c2[0] == rb(j) and c2[1] == ub(j) and c2[2] == kb(j) and c2[3] == Lb(j) and
                                c1[8] == ra(h) and c1[9] == ua(h) and c1[10] == ka(h) and c1[11] == La(h) and
                                c1[0] == c2[8] and c1[1] == c2[9] and c1[2] == c2[10] and c1[3] == c2[11]):
                            prefactor = (-sqrt(5)) * (-1) ** (Lb(j) + jp + Jt) * wigner_6j(La(h), jp, Jt, jp, Lb(j), 2)
                            more_factors = (sqrt(5) * (-1) ** (La(h) + Lb(j) + 2) *
                                            wigner_6j(2, 2, 2,Lb(j), La(h), c1[3]) *
                                            sqrt((2 * c1[11] + 1) * (2 * c1[3] + 1)))
                            sqrt_frac = sqrt((c2[0] + 1) * (c2[1] + 1) * (c2[0] + c2[1] + 2) /
                                             ((c2[8] + 1) * (c2[9] + 1) * (c2[8] + c2[9] + 2)))
                            sum_value += (prefactor * more_factors * c1[13] * c2[13] *
                                          threebar1(s1) * threebar2(s2) * sqrt_frac)
                row.append(common_part + [sum_value])
            else:
                row.append(common_part + [0])
        y3.append(row)

    chi = sp.symbols('chi')
    aF = sp.symbols('aF')
    Qq = [[(y1[i][j][8] + (chi + sqrt(7) / 2.0) * y3[i][j][8]) for j in range(len(basis))]
          for i in range(len(basis))]

    H = Qq
    H_df = aF*pd.DataFrame(H)

    filename = matrix_dir + "h{}c.csv".format(int(Jt*2))
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    H_df.to_csv(filename, index=False, header=False)
