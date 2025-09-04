import pandas as pd
import sympy as sp
import os
from sympy.physics.wigner import wigner_6j
from functions.definitions import integer_part, kronecker_delta, sign


def construct_Ha_even(matrix_dir, m, L, cg):
    ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1) for u in range(0, 2 * m - 6 * b + 1, 2) if
          2 * m - 6 * b - 2 * u >= 0]

    H, _ = construct_H(m, L, cg, ru)

    # Convert the matrix H to a DataFrame for easier export
    H_df = pd.DataFrame(H)

    filename = matrix_dir + "h{}a.csv".format(L)
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    H_df.to_csv(filename, index=False,header=False)


def construct_Ha_odd(matrix_dir, m, Jt, cg, jp):
    tt = []

    ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1) for u in range(0, 2 * m - 6 * b + 1, 2) if
          2 * m - 6 * b - 2 * u >= 0]

    for L in range(int(abs(Jt - jp)), int(Jt + jp + 1)):
        try:
            H, basis1 = construct_H(m, L, cg, ru)
            Ht = []
            for h in range(len(basis1)):
                row = []
                for k in range(len(basis1)):
                    row.append(basis1[h][:4] + basis1[k][:4] + [H[h][k]])
                Ht.append(row)
            tt.append([L, basis1.copy(), Ht])
        except:
            tt.append([L, [], []])

    # Flatten tt to get all base vectors and matrix entries
    base = [item for t in tt for item in t[1]]
    ggy3 = [item for t in tt for item in t[2]]

    hhy3 = []
    for j in range(len(ggy3)):
        row = []
        for k in range(len(base)):
            try:
                match = next(
                    x for x in ggy3[j]
                    if base[k][0] == x[4] and base[k][1] == x[5] and
                    base[k][2] == x[6] and base[k][3] == x[7]
                )
                row.append(match[8])  # 9th element (index 8)
            except (StopIteration):
                row.append(0)
        hhy3.append(row)
    H_all = pd.DataFrame(hhy3)

    filename = matrix_dir + "h{}a.csv".format(int(Jt*2))
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    H_all.to_csv(filename, index=False, header=False)


def construct_H(m, L, cg, ru):
    import numpy as np
    cg1, cg2, cg3, = cg[0], cg[1], cg[2]

    basis1 = [[r[0], r[1], s, L]
              for r in ru
              for s in range(1, integer_part((r[0] + r[1] + 2 - L) / 2) *
                             kronecker_delta(sign((r[0] + r[1] + 2 - L) / 2), 1) -
                             integer_part((r[0] + 1 - L) / 2) *
                             kronecker_delta(sign((r[0] + 1 - L) / 2), 1) -
                             integer_part((r[1] + 1 - L) / 2) *
                             kronecker_delta(sign((r[1] + 1 - L) / 2), 1) + 1)]

    basis2 = [[r[0], r[1], s, a]
              for a in range(abs(L - 2), L + 3)
              for r in ru
              for s in range(1, integer_part((r[0] + r[1] + 2 - a) / 2) *
                             kronecker_delta(sign((r[0] + r[1] + 2 - a) / 2), 1) -
                             integer_part((r[0] + 1 - a) / 2) *
                             kronecker_delta(sign((r[0] + 1 - a) / 2), 1) -
                             integer_part((r[1] + 1 - a) / 2) *
                             kronecker_delta(sign((r[1] + 1 - a) / 2), 1) + 1)]

    basis3 = [[r[0], r[1], s, a]
              for a in range(abs(L - 3), L + 4)
              for r in ru
              for s in range(1, integer_part((r[0] + r[1] + 2 - a) / 2) *
                             kronecker_delta(sign((r[0] + r[1] + 2 - a) / 2), 1) -
                             integer_part((r[0] + 1 - a) / 2) *
                             kronecker_delta(sign((r[0] + 1 - a) / 2), 1) -
                             integer_part((r[1] + 1 - a) / 2) *
                             kronecker_delta(sign((r[1] + 1 - a) / 2), 1) + 1)]

    cg1L = [item for item in cg1 if item[11] == L]  # Using zero-based index
    cg2L = [item for item in cg2 if item[3] == L]
    cg3L = [item for Lm in range(abs(L - 2), L + 3) for item in cg1 if item[11] == Lm]
    cg4L = [item for Lm in range(abs(L - 2), L + 3) for item in cg2 if item[11] == Lm]
    cg44L = [item for Lm in range(abs(L - 2), L + 3) for item in cg2 if item[3] == Lm]

    cg5L = [item for item in cg3 if item[11] == L]
    cg55L = [item for Lm in range(abs(L - 2), L + 3) for item in cg3 if item[11] == Lm]
    cg555L = [item for item in cg3 if item[3] == L]

    # Assuming basis1 and basis2 are already defined lists

    # Functions for basis1
    def ra(h):
        return basis1[h][0]  # Using zero-based indexing

    def ua(h):
        return basis1[h][1]

    def ka(h):
        return basis1[h][2]

    def La(h):
        return basis1[h][3]

    # Functions for basis2
    def rb(j):
        return basis2[j][0]

    def ub(j):
        return basis2[j][1]

    def kb(j):
        return basis2[j][2]

    def Lb(j):
        return basis2[j][3]

    def rbb(jb):
        return basis2[jb][0]

    def ubb(jb):
        return basis2[jb][1]

    def kbb(jb):
        return basis2[jb][2]

    def Lbb(jb):
        return basis2[jb][3]

    # Functions for another basis1
    def rc(k):
        return basis1[k][0]

    def uc(k):
        return basis1[k][1]

    def kc(k):
        return basis1[k][2]

    def Lc(k):
        return basis1[k][3]

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

    LL = np.array([[L * (L + 1) if (ra(h) == rc(k) and ua(h) == uc(k) and ka(h) == kc(k) and La(h) == Lc(k))
                    else 0
                    for k in range(len(basis1))]
                   for h in range(len(basis1))])

    # 2. delta21 Calculation
    delta21 = np.array([[1 if (basis2[h][0] == basis1[j][0] and
                               basis2[h][1] == basis1[j][1] and
                               basis2[h][2] == basis1[j][2] and
                               basis2[h][3] == basis1[j][3])
                         else 0
                         for j in range(len(basis1))]
                        for h in range(len(basis2))])

    # 3. aa Calculation
    aa = []
    for h in range(len(basis1)):
        for j in range(len(basis2)):
            # Find all matching entries for the current (h, j)
            entries = [entry for entry in cg5L
                       if (entry[0] == rb(j) and entry[1] == ub(j) and entry[2] == kb(j) and entry[3] == Lb(j) and
                           entry[8] == ra(h) and entry[9] == ua(h) and entry[10] == ka(h) and entry[11] == La(h))]

            if entries:
                aa.extend(entries)  # Add the found entries
            else:
                aa.append(0)
            # 4. bb Calculation
    bb = []
    for j in range(len(basis2)):
        for jb in range(len(basis2)):
            # Find all matching entries for the current (j, jb)
            entries = [entry for entry in cg55L
                       if
                       (entry[0] == rbb(jb) and entry[1] == ubb(jb) and entry[2] == kbb(jb) and entry[3] == Lbb(jb) and
                        entry[8] == rb(j) and entry[9] == ub(j) and entry[10] == kb(j) and entry[11] == Lb(j))]

            if entries:
                bb.extend(entries)  # Add the found entries
            else:
                bb.append(0)

    # 5. cc Calculation
    cc = []
    for jb in range(len(basis2)):
        for k in range(len(basis1)):
            # Find all matching entries for the current (jb, k)
            entries = [entry for entry in cg55L
                       if (entry[0] == rc(k) and entry[1] == uc(k) and entry[2] == kc(k) and entry[3] == Lc(k) and
                           entry[8] == rbb(jb) and entry[9] == ubb(jb) and entry[10] == kbb(jb) and entry[11] == Lbb(
                            jb))]

            if entries:
                cc.extend(entries)  # Add the found entries
            else:
                cc.append(0)

            # 6. aa1 Calculation
    aa1 = []
    for s in range(len(aa)):
        if isinstance(aa[s], list) and len(aa[s]) > 0:  # Check if aa[s] is a list and not empty
            value = aa[s][13] * (np.sqrt(2 * aa[s][11] + 1) * 3 / (2 * np.sqrt(6)) *
                                 (-1) ** (1 + (1 if aa[s][1] == 0 else 0)) *  # KroneckerDelta equivalent
                                 np.sqrt(4 / 3 * (aa[s][0] ** 2 + aa[s][1] ** 2 +
                                                  3 * aa[s][0] + 3 * aa[s][1] +
                                                  aa[s][0] * aa[s][1])))
            aa1.append(value)  # Append the calculated value
        else:
            aa1.append(0)

    # 7. y1 Calculation
    y1 = np.array([[aa1[j + (h - 1) * len(basis2)] if (j + (h - 1) * len(basis2)) < len(aa1) else 0
                    for j in range(len(basis2))]
                   for h in range(1, len(basis1) + 1)])

    # 1. ff Calculation
    ff = [
        [
            [entry for entry in cg4L
             if entry[0] == rb(j) and entry[1] == ub(j) and entry[2] == kb(j) and entry[3] == Lb(j)],
            [entry for entry in cg1L
             if entry[8] == ra(h) and entry[9] == ua(h) and entry[10] == ka(h) and entry[11] == La(h)]
        ]
        for h in range(len(basis1)) for j in range(len(basis2))
    ]

    # This preserves the structure (120, 2) for ff

    # 2. gg Calculation
    gg = [
        [
            [entry for entry in cg44L
             if entry[0] == rbb(jb) and entry[1] == ubb(jb) and entry[2] == kbb(jb) and entry[3] == Lbb(jb)],
            [entry for entry in cg3L
             if entry[8] == rb(j) and entry[9] == ub(j) and entry[10] == kb(j) and entry[11] == Lb(j)]
        ]
        for j in range(len(basis2)) for jb in range(len(basis2))
    ]

    # This preserves the structure (120, 2) for gg

    # 3. hh Calculation
    hh = [
        [
            [entry for entry in cg2L
             if entry[0] == rc(k) and entry[1] == uc(k) and entry[2] == kc(k) and entry[3] == Lc(k)],
            [entry for entry in cg3L
             if entry[8] == rbb(jb) and entry[9] == ubb(jb) and entry[10] == kbb(jb) and entry[11] == Lbb(jb)]
        ]
        for jb in range(len(basis2)) for k in range(len(basis1))
    ]
    # This preserves the structure (120, 2) for hh

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

                    term1 = np.sqrt(5) * (-1) ** (ff[s][1][kp][11] + ff[s][0][k][3])
                    term2 = (wigner_6j(2, 2, 2, ff[s][0][k][3], ff[s][1][kp][11], ff[s][1][kp][3]))
                    term3 = np.sqrt((2 * ff[s][1][kp][11] + 1) * (2 * ff[s][1][kp][3] + 1))
                    term4 = ff[s][1][kp][13] * ff[s][0][k][13]

                    # Calculate the inner product of all terms
                    combined_term = term1 * term2 * term3 * term4

                    # Nested conditions for the second term involving different comparisons
                    condition1 = (ff[s][1][kp][8] == (ff[s][1][kp][0] + 2) and
                                  ff[s][1][kp][9] == ff[s][1][kp][1])
                    term5 = (np.sqrt((2 * (m - 1) + 2 * ff[s][1][kp][0] + ff[s][1][kp][1] + 12) *
                                     (ff[s][1][kp][0] + 2) * (ff[s][1][kp][0] + ff[s][1][kp][1] + 3) /
                                     (6 * (ff[s][1][kp][8] + 1) * (ff[s][1][kp][8] + ff[s][1][kp][9] + 2)))
                             if condition1 else 0)

                    condition2 = (ff[s][1][kp][8] == (ff[s][1][kp][0] - 2) and
                                  ff[s][1][kp][9] == (ff[s][1][kp][1] + 2))
                    term6 = (np.sqrt((2 * (m - 1) - ff[s][1][kp][0] + ff[s][1][kp][1] + 9) *
                                     ff[s][1][kp][0] * (ff[s][1][kp][1] + 2) /
                                     (6 * (ff[s][1][kp][8] + 1) * (ff[s][1][kp][9] + 1)))
                             if condition2 else 0)

                    condition3 = (ff[s][1][kp][8] == ff[s][1][kp][0] and
                                  ff[s][1][kp][9] == (ff[s][1][kp][1] - 2))
                    term7 = (np.sqrt((2 * (m - 1) - ff[s][1][kp][0] - 2 * ff[s][1][kp][1] + 6) *
                                     ff[s][1][kp][1] * (ff[s][1][kp][0] + ff[s][1][kp][1] + 1) /
                                     (6 * (ff[s][1][kp][8] + ff[s][1][kp][9] + 2) * (ff[s][1][kp][9] + 1)))
                             if condition3 else 0)

                    final_term = term5 + term6 + term7
                    if (ff[s][0][k][0] == ff[s][0][k][8] + 2 and
                            ff[s][0][k][1] == ff[s][0][k][9]):
                        term8 = np.sqrt(
                            (2 * (m - 1) + 2 * ff[s][0][k][8] + ff[s][0][k][9] + 12) *
                            (ff[s][0][k][8] + 2) *
                            (ff[s][0][k][8] + ff[s][0][k][9] + 3) /
                            (6 * (ff[s][0][k][0] + 1) * (ff[s][0][k][0] + ff[s][0][k][1] + 2))
                        )
                    else:
                        term8 = 0
                    if (ff[s][0][k][0] == ff[s][0][k][8] - 2 and
                            ff[s][0][k][1] == ff[s][0][k][9] + 2):
                        term9 = np.sqrt(
                            (2 * (m - 1) - ff[s][0][k][8] + ff[s][0][k][9] + 9) *
                            ff[s][0][k][8] *
                            (ff[s][0][k][9] + 2) /
                            (6 * (ff[s][0][k][0] + 1) * (ff[s][0][k][1] + 1))
                        )
                    else:
                        term9 = 0
                    if (ff[s][0][k][0] == ff[s][0][k][8] and
                            ff[s][0][k][1] == ff[s][0][k][9] - 2):
                        term10 = np.sqrt(
                            (2 * (m - 1) - ff[s][0][k][8] - 2 * ff[s][0][k][9] + 6) *
                            ff[s][0][k][9] *
                            (ff[s][0][k][8] + ff[s][0][k][9] + 1) /
                            (6 * (ff[s][0][k][0] + ff[s][0][k][1] + 2) * (ff[s][0][k][1] + 1))
                        )
                    else:
                        term10 = 0
                    add_term = (term8 + term9 + term10) * np.sqrt(
                        (ff[s][0][k][0] + 1) *
                        (ff[s][0][k][1] + 1) *
                        (ff[s][0][k][0] + ff[s][0][k][1] + 2) /
                        ((ff[s][0][k][8] + 1) * (ff[s][0][k][9] + 1) * (ff[s][0][k][8] + ff[s][0][k][9] + 2))
                    )

                    # Add the product of combined_term and final_term to the accumulator
                    sum_value += combined_term * final_term * add_term

        # Append the accumulated sum to the ff1 list
        ff1.append(float(sum_value))
    # Assuming ff1, basis1, and basis2 have been defined as per the Mathematica structure
    y3 = [[ff1[j + (h - 1) * len(basis2)] for j in range(len(basis2))] for h in range(1, len(basis1) + 1)]
    y3 = np.array(y3)

    # y3333
    def compute_y3333(basis1, basis2, ff1, L):
        y3333 = np.zeros((len(basis1), len(basis2)))  # Initialize array

        for h in range(1, len(basis1) + 1):
            for j in range(1, len(basis2) + 1):
                J = basis2[j - 1][3]  # basis2[[j]][[4]] in Mathematica (1-based index)
                sqrt_term = np.sqrt(J * (J + 1) * (2 * J + 1))

                six_j1 = float(wigner_6j(1, 2, 1, J, L, L).doit())
                six_j2 = float(wigner_6j(1, 2, 1, L, J, J).doit())

                y3333[h - 1, j - 1] = np.sqrt(sqrt_term * six_j1 * six_j2) * ff1[j - 1 + (h - 1) * len(basis2)]

        return y3333

    y3333 = compute_y3333(basis1, basis2, ff1, L)

    # y1111
    import numpy as np
    from sympy.physics.quantum.cg import Wigner6j

    def compute_y1111(basis1, basis2, aa1, L):
        y1111 = np.zeros((len(basis1), len(basis2)))  # Initialize the array

        for h in range(1, len(basis1) + 1):
            for j in range(1, len(basis2) + 1):
                J = basis2[j - 1][3]  # basis2[[j]][[4]] in Mathematica (1-based index)
                sqrt_term = np.sqrt(J * (J + 1) * (2 * J + 1))

                six_j1 = float(wigner_6j(1, 2, 1, J, L, L).doit())
                six_j2 = float(wigner_6j(1, 2, 1, L, J, J).doit())

                y1111[h - 1, j - 1] = np.sqrt(sqrt_term * six_j1 * six_j2) * aa1[j - 1 + (h - 1) * len(basis2)]

        return y1111

    y1111 = compute_y1111(basis1, basis2, aa1, L)
    nd = np.zeros((len(basis1), len(basis1)))
    for h in range(len(basis1)):
        for k in range(len(basis1)):
            summation = 0
            for s1 in range(len(cg1L)):
                for s2 in range(len(cg2L)):
                    if (cg1L[s1][8] == ra(h) and
                            cg1L[s1][9] == ua(h) and
                            cg1L[s1][10] == ka(h) and
                            cg2L[s2][0] == rc(k) and
                            cg2L[s2][1] == uc(k) and
                            cg2L[s2][2] == kc(k) and
                            cg1L[s1][0] == cg2L[s2][8] and
                            cg1L[s1][1] == cg2L[s2][9] and
                            cg1L[s1][2] == cg2L[s2][10] and
                            cg1L[s1][3] == cg2L[s2][11]):
                        term = ((-1) ** (L - cg1L[s1][3]) *
                                np.sqrt((2 * cg2L[s2][11] + 1) / (2 * L + 1)) *
                                cg1L[s1][13] * cg2L[s2][13] *
                                threebar1(s1) * threebar2(s2) *
                                np.sqrt((cg2L[s2][0] + 1) * (cg2L[s2][1] + 1) *
                                        (cg2L[s2][0] + cg2L[s2][1] + 2) /
                                        ((cg1L[s1][0] + 1) * (cg1L[s1][1] + 1) *
                                         (cg1L[s1][0] + cg1L[s1][1] + 2))))
                        summation += term
            nd[h][k] = summation

    # Define the shape of the matrices
    num_rows, num_cols = 4, 30

    # Create symbolic matrix representations for y1 and y3

    # Define the symbolic variable chi
    chi = sp.symbols('chi')
    a1, a2, a4, a5, a6 = sp.symbols('a1 a2 a4 a5 a6')

    # Define QQ symbolically
    QQ = (1 / (2 * L + 1)) * (y1 + (chi + np.sqrt(7) / 2) * y3).dot((y1 + (chi + np.sqrt(7) / 2) * y3).T)

    LddL = L * (L + 1) * wigner_6j(1, 2, 1, L, L, L) * np.dot(y3, delta21)

    LQL = L * (L + 1) * wigner_6j(1, 2, 1, L, L, L) * np.dot((y1 + (chi + np.sqrt(7) / 2.0) * y3), delta21)

    term = y1111 + (chi + np.sqrt(7) / 2.0) * y3333
    LQQL = 3 * np.sqrt(L * (L + 1.) / (2 * L + 1)) * (term @ term.T)

    # Define the symbolic expression for H
    # + v51 * LddL # extra term that could be added
    H = a1 * nd + a2 * QQ + a5 * LL + a4 * LQL + a6 * LQQL
    return H, basis1



