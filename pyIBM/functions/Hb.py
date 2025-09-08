import pandas as pd
import numpy as np
from math import floor,fabs
from sympy.physics.wigner import wigner_6j
import sympy as sp
import os
from functions.definitions import integer_part, kronecker_delta, sign


def construct_Hb_even(matrix_dir, m, L, cg):
    def generate_ru(m):
        ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1)
              for u in range(0, 2 * m - 6 * b + 1, 2)
              if 2 * m - 6 * b - 2 * u >= 0]
        return ru

    ru1 = generate_ru(m)
    ru2 = generate_ru(m - 1)
    ru3 = generate_ru(m - 2)
    ru4 = generate_ru(m - 3)

    sixd, _ = construct_sixd(m, L, cg, ru1, ru2, ru3, ru4)

    a3 = sp.symbols('a3')
    H3d = a3 * sixd
    H3d_df = pd.DataFrame(H3d)
    filename = matrix_dir + "h{}b.csv".format(L)
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    H3d_df.to_csv(filename, index=False, header=False)
    return H3d_df


def construct_Hb_odd(matrix_dir, m, Jt, cg, jp):
    tt = []

    def generate_ru(m):
        ru = [(2 * m - 6 * b - 2 * u, u) for b in range(0, m // 3 + 1)
              for u in range(0, 2 * m - 6 * b + 1, 2)
              if 2 * m - 6 * b - 2 * u >= 0]
        return ru

    ru1 = generate_ru(m)
    ru2 = generate_ru(m - 1)
    ru3 = generate_ru(m - 2)
    ru4 = generate_ru(m - 3)

    for L in range(int(abs(Jt - jp)), int(Jt + jp + 1)):
        try:
            sixdd, basis1 = construct_sixd(m, L, cg, ru1, ru2, ru3, ru4)
            sixd = []
            for h in range(len(basis1)):
                row = []
                for k in range(len(basis1)):
                    row.append(basis1[h][:4] + basis1[k][:4] + [sixdd[h][k]])
                sixd.append(row)
            tt.append([L, basis1.copy(), sixd])
        except:
            tt.append([L, [], []])


    # Flatten tt to get all base vectors and matrix entries
    base = [item for t in tt for item in t[1]]
    ggy3 = [item for t in tt for item in t[2]]

    # Build matrix hhy3
    hhy3 = []
    for j in range(len(ggy3)):
        row = []
        for k in range(len(base)):
            try:
                match = next(
                    entry for entry in ggy3[j]
                    if base[k][0] == entry[4] and base[k][1] == entry[5] and
                    base[k][2] == entry[6] and base[k][3] == entry[7]
                )
                row.append(match[8])
            except StopIteration:
                row.append(0)
        hhy3.append(row)

    # H matrix
    HM= hhy3

    # Select non-zero matrix entries
    sixdeven = [
        base[j][:4] + base[k][:4] + [HM[j][k]]
        for j in range(len(base))
        for k in range(len(base))
        if abs(HM[j][k]) > 0
    ]

    # Build ru table
    ru = [
        [2 * m - 6 * b - 2 * u, u]
        for b in range(int(m / 3) + 1)
        for u in range(0, 2 * m - 6 * b + 1, 2)
        if 2 * m - 6 * b - 2 * u >= 0
    ]

    # Construct basisodd
    basisodd = []
    for L in range(int(abs(Jt - jp)), int(Jt + jp + 1)):
        for j in range(len(ru)):
            m1 = ru[j][0]
            m2 = ru[j][1]

            term1 = integer_part((m1 + m2 + 2 - L) / 2) * kronecker_delta(int((m1 + m2 + 2 - L) / 2 > 0), 1)
            term2 = integer_part((m1 + 1 - L) / 2) * kronecker_delta(int((m1 + 1 - L) / 2 > 0), 1)
            term3 = integer_part((m2 + 1 - L) / 2.0) * kronecker_delta(int((m2 + 1 - L) / 2.0 > 0), 1)

            upper_s = term1 - term2 - term3

            for s in range(1, upper_s + 1):
                basisodd.append([m1, m2, s, L])

    # Compute sixdodd
    sixdodd = []
    for j in range(len(basisodd)):
        row = []
        for k in range(len(basisodd)):
            val = sum(
                entry[8]
                for entry in sixdeven
                if entry[:4] == basisodd[j] and entry[4:8] == basisodd[k]
            )
            row.append(val)
        sixdodd.append(row)

    # Export to CSV
    a3 = sp.symbols('a3')
    H_df = a3*pd.DataFrame(sixdodd)

    filename = matrix_dir + "h{}b.csv".format(int(Jt*2))
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    H_df.to_csv(filename, index=False, header=False)


def construct_sixd(m, L, cg, ru1, ru2, ru3, ru4):
    def generate_basis1(ru, L):
        basis = []
        for j in range(len(ru)):
            ru1_j1 = ru[j][0]
            ru1_j2 = ru[j][1]
            for s in range(1,
                           floor((ru1_j1 + ru1_j2 + 2 - L) / 2) *
                           kronecker_delta(sign((ru1_j1 + ru1_j2 + 2 - L) / 2), 1)
                           - floor((ru1_j1 + 1 - L) / 2) *
                           kronecker_delta(sign((ru1_j1 + 1 - L) / 2), 1)
                           - floor((ru1_j2 + 1 - L) / 2) *
                           kronecker_delta(sign((ru1_j2 + 1 - L) / 2), 1) + 1):
                basis.append([ru1_j1, ru1_j2, s, L])
        return basis

    def generate_basis2(ru, L):
        basis = []
        for a in range(int(fabs(L - 2)), L + 3):  # a runs from Abs(L - 2) to L + 2
            for j in range(len(ru)):
                ru_j1 = ru[j][0]
                ru_j2 = ru[j][1]
                for s in range(1,
                               floor((ru_j1 + ru_j2 + 2 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + ru_j2 + 2 - a) / 2), 1)
                               - floor((ru_j1 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + 1 - a) / 2), 1)
                               - floor((ru_j2 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j2 + 1 - a) / 2), 1) + 1):
                    basis.append([ru_j1, ru_j2, s, a])
        return basis

    def generate_basis3(ru, L):
        basis = []
        for a in range(int(fabs(L - 2)), L + 3):  # a runs from Abs(L - 2) to L + 2
            for j in range(len(ru)):
                ru_j1 = ru[j][0]
                ru_j2 = ru[j][1]
                for s in range(1,
                               floor((ru_j1 + ru_j2 + 2 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + ru_j2 + 2 - a) / 2), 1)
                               - floor((ru_j1 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + 1 - a) / 2), 1)
                               - floor((ru_j2 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j2 + 1 - a) / 2), 1) + 1):
                    basis.append([ru_j1, ru_j2, s, a])
        return basis

    def generate_basis4(ru, L):
        basis = []
        for a in range(int(fabs(L - 3)), L + 4):  # a runs from Abs(L - 3) to L + 3
            for j in range(len(ru)):
                ru_j1 = ru[j][0]
                ru_j2 = ru[j][1]
                for s in range(1,
                               floor((ru_j1 + ru_j2 + 2 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + ru_j2 + 2 - a) / 2), 1)
                               - floor((ru_j1 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j1 + 1 - a) / 2), 1)
                               - floor((ru_j2 + 1 - a) / 2) *
                               kronecker_delta(sign((ru_j2 + 1 - a) / 2), 1) + 1):
                    basis.append([ru_j1, ru_j2, s, a])
        return basis

    basis1 = generate_basis1(ru1, L)
    basis2 = generate_basis2(ru2, L)
    basis3 = generate_basis3(ru3, L)
    basis4 = generate_basis4(ru4, L)

    # Function to select rows where the 12th column (index 11) equals L
    def select_by_L_df(df, L, column_index=11):
        return df.loc[df.iloc[:, column_index] == L]

    # Function to select rows within a range of L values for the 12th column (index 11)
    def select_by_range_L_df(df, L, column_index=11):
        # Initialize an empty DataFrame to store the ordered results
        ordered_df = pd.DataFrame()

        # Iterate over the range in the required order (L-2 to L+2)
        for i in range(L - 2, L + 3):
            # Select the rows where the value in the column matches i
            temp_df = df.loc[df.iloc[:, column_index] == i]
            # Concatenate to the ordered DataFrame
            ordered_df = pd.concat([ordered_df, temp_df], ignore_index=True)

        return ordered_df

    cg1_df, cg2_df, cg3_df, = cg[3], cg[4], cg[5]

    # Apply the function for cg2L_df and cg3L_df
    cg1L_df = select_by_L_df(cg1_df, L)
    cg2L_df = select_by_range_L_df(cg2_df, L)
    cg3L_df = select_by_range_L_df(cg3_df, L)

    # Assume cg1L_df, basis1, basis2, and m are defined beforehand

    def calculate_term1(cg1L_row, m):
        if cg1L_row[8] == cg1L_row[0] + 2 and cg1L_row[9] == cg1L_row[1]:
            numerator = (2 * (m - 1) + 2 * cg1L_row[0] + cg1L_row[1] + 12) * (cg1L_row[0] + 2) * (
                        cg1L_row[0] + cg1L_row[1] + 3)
            denominator = 6 * (cg1L_row[8] + 1) * (cg1L_row[8] + cg1L_row[9] + 2)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_term2(cg1L_row, m):
        if cg1L_row[8] == cg1L_row[0] - 2 and cg1L_row[9] == cg1L_row[1] + 2:
            numerator = (2 * (m - 1) - cg1L_row[0] + cg1L_row[1] + 9) * cg1L_row[0] * (cg1L_row[1] + 2)
            denominator = 6 * (cg1L_row[8] + 1) * (cg1L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_term3(cg1L_row, m):
        if cg1L_row[8] == cg1L_row[0] and cg1L_row[9] == cg1L_row[1] - 2:
            numerator = (2 * (m - 1) - cg1L_row[0] - 2 * cg1L_row[1] + 6) * cg1L_row[1] * (
                        cg1L_row[0] + cg1L_row[1] + 1)
            denominator = 6 * (cg1L_row[8] + cg1L_row[9] + 2) * (cg1L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_matrix_element(cg1L_row, basis1_h, basis2_j, m):
        if np.all(cg1L_row[8:12] == basis1_h) and np.all(cg1L_row[0:4] == basis2_j):
            term1 = calculate_term1(cg1L_row, m)
            term2 = calculate_term2(cg1L_row, m)
            term3 = calculate_term3(cg1L_row, m)
            return np.sqrt(2 * cg1L_row[11] + 1) * cg1L_row[13] * (term1 + term2 + term3)
        return 0

    # Initialize matrix d12
    d12 = np.zeros((len(basis1), len(basis2)))

    # Loop over basis1, basis2, and cg1L_df rows to populate d12 matrix
    for h in range(len(basis1)):
        for j in range(len(basis2)):
            d12[h, j] = sum(
                calculate_matrix_element(cg1L_df.iloc[s1], basis1[h], basis2[j], m) for s1 in range(len(cg1L_df)))

    def calculate_d22_term1(cg2L_row, m):
        if cg2L_row[8] == cg2L_row[0] + 2 and cg2L_row[9] == cg2L_row[1]:
            numerator = (2 * ((m - 1) - 1) + 2 * cg2L_row[0] + cg2L_row[1] + 12) * (cg2L_row[0] + 2) * (
                        cg2L_row[0] + cg2L_row[1] + 3)
            denominator = 6 * (cg2L_row[8] + 1) * (cg2L_row[8] + cg2L_row[9] + 2)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d22_term2(cg2L_row, m):
        if cg2L_row[8] == cg2L_row[0] - 2 and cg2L_row[9] == cg2L_row[1] + 2:
            numerator = (2 * ((m - 1) - 1) - cg2L_row[0] + cg2L_row[1] + 9) * cg2L_row[0] * (cg2L_row[1] + 2)
            denominator = 6 * (cg2L_row[8] + 1) * (cg2L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d22_term3(cg2L_row, m):
        if cg2L_row[8] == cg2L_row[0] and cg2L_row[9] == cg2L_row[1] - 2:
            numerator = (2 * ((m - 1) - 1) - cg2L_row[0] - 2 * cg2L_row[1] + 6) * cg2L_row[1] * (
                        cg2L_row[0] + cg2L_row[1] + 1)
            denominator = 6 * (cg2L_row[8] + cg2L_row[9] + 2) * (cg2L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d22_matrix_element(cg2L_row, basis2_j, basis3_jb, m, L):
        if np.all(cg2L_row[8:12] == basis2_j) and np.all(cg2L_row[0:4] == basis3_jb):
            six_j = float(wigner_6j(2, 2, 2, int(cg2L_row[3]), L, int(cg2L_row[11])))
            factor = np.sqrt(2 * cg2L_row[11] + 1) * cg2L_row[13]
            term1 = calculate_d22_term1(cg2L_row, m)
            term2 = calculate_d22_term2(cg2L_row, m)
            term3 = calculate_d22_term3(cg2L_row, m)
            return six_j * factor * (term1 + term2 + term3)
        return 0

    # Initialize d22 matrix
    d22 = np.zeros((len(basis2), len(basis3)))

    # Loop to populate d22
    for j in range(len(basis2)):
        for jb in range(len(basis3)):
            d22[j, jb] = sum(calculate_d22_matrix_element(cg2L_df.iloc[s2], basis2[j], basis3[jb], m, L) for s2 in
                             range(len(cg2L_df)))

    def calculate_d23_term1(cg3L_row, m):
        if cg3L_row[8] == cg3L_row[0] + 2 and cg3L_row[9] == cg3L_row[1]:
            numerator = (2 * ((m - 2) - 1) + 2 * cg3L_row[0] + cg3L_row[1] + 12) * (cg3L_row[0] + 2) * (
                        cg3L_row[0] + cg3L_row[1] + 3)
            denominator = 6 * (cg3L_row[8] + 1) * (cg3L_row[8] + cg3L_row[9] + 2)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d23_term2(cg3L_row, m):
        if cg3L_row[8] == cg3L_row[0] - 2 and cg3L_row[9] == cg3L_row[1] + 2:
            numerator = (2 * ((m - 2) - 1) - cg3L_row[0] + cg3L_row[1] + 9) * cg3L_row[0] * (cg3L_row[1] + 2)
            denominator = 6 * (cg3L_row[8] + 1) * (cg3L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d23_term3(cg3L_row, m):
        if cg3L_row[8] == cg3L_row[0] and cg3L_row[9] == cg3L_row[1] - 2:
            numerator = (2 * ((m - 2) - 1) - cg3L_row[0] - 2 * cg3L_row[1] + 6) * cg3L_row[1] * (
                        cg3L_row[0] + cg3L_row[1] + 1)
            denominator = 6 * (cg3L_row[8] + cg3L_row[9] + 2) * (cg3L_row[9] + 1)
            return np.sqrt(numerator / denominator)
        return 0

    def calculate_d23_matrix_element(cg3L_row, basis3_jb, basis4_kb, m, L):
        if np.all(cg3L_row[8:12] == basis3_jb) and np.all(cg3L_row[0:4] == basis4_kb):
            six_j = float(wigner_6j(2, 2, 3, int(cg3L_row[3]), L, int(cg3L_row[11])))
            factor = -np.sqrt(35) * (-1) ** (cg3L_row[3] + cg3L_row[11]) * np.sqrt(2 * cg3L_row[11] + 1) * cg3L_row[13]
            term1 = calculate_d23_term1(cg3L_row, m)
            term2 = calculate_d23_term2(cg3L_row, m)
            term3 = calculate_d23_term3(cg3L_row, m)
            return six_j * factor * (term1 + term2 + term3)
        return 0

    # Initialize d23 matrix
    d23 = np.zeros((len(basis3), len(basis4)))

    # Loop to populate d23
    for jb in range(len(basis3)):
        for kb in range(len(basis4)):
            d23[jb, kb] = sum(calculate_d23_matrix_element(cg3L_df.iloc[s3], basis3[jb], basis4[kb], m, L) for s3 in
                              range(len(cg3L_df)))

    ddddag = np.dot(np.dot(d12, d22), d23)
    sixd = -7 / 30 * np.dot(ddddag, ddddag.T) / (2 * L + 1)

    return sixd, basis1

