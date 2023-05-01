import numpy as np
import os
from copy import deepcopy
import pandas as pd

def find_shared_indexes(arr):
    shared_indexes = []
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            val = arr[i][j]
            if val not in [x[0] for x in shared_indexes]:
                for k in range(i + 1, len(arr)):
                    if val in arr[k]:
                        shared_indexes.append((val, i, k))
                        break
    return shared_indexes

def flip_repeated_values(A, B):
    # Create a dictionary to keep track of repeated values and how many times they've appeared
    repeated_values = {}
    # Loop through each row of the matrix
    for i in range(len(A)):
        # Loop through each element in the row
        for j in range(len(A[i])):
            # If the element has already appeared once before
            if A[i][j] in repeated_values and repeated_values[A[i][j]] == 1:
                # Flip the sign of the element to negative
                A[i][j] = -A[i][j]
            # Otherwise, add the element to the dictionary of repeated values
            else:
                repeated_values[A[i][j]] = 1
    ones_matrix = [[1 if val >= 0 else -1 for val in row] for row in A]
    flipped = np.multiply(ones_matrix, B)
    return flipped

def deltaQ(Q, K):
    sum_numerator = np.zeros_like(Q)
    sum_denominator = np.zeros_like(Q)
    for i in range(len(sum_denominator)):
        sum_numerator[i] = K[i] * Q[i] * abs(Q[i]) ** (n - 1)
        sum_denominator[i] = K[i] * abs(Q[i]) ** (n - 1)
    return -np.sum(sum_numerator) / (n * np.sum(sum_denominator))


def apply_correction(arr, Q_loops, dQ):
    n_loops = len(arr)
    shared_indexes = find_shared_indexes(arr)
    updated_Q_loops = np.zeros((n_loops, len(arr[0])))
    for idx, row in enumerate(arr):
        updated_Q_loops[idx] = Q_loops[idx] + dQ[idx]
        for val, i, j in shared_indexes:
            if idx == i:
                updated_Q_loops[idx][arr[i].index(val)] += - dQ[j]
            elif idx == j:
                updated_Q_loops[idx][arr[j].index(val)] += -dQ[i]
    return updated_Q_loops


def solve_hybrid_pipes(Q_loops, K_loops, loop_idexes,tol=1e-7, max_iter=100):
    Q_loops_old = Q_loops.copy()
    dQ = np.ones(n_loops)
    count = 0
    residuals = np.sum(np.absolute(dQ))
    while residuals > tol:
        print(f"Iteration {count}: Residual {residuals}")
        for i in range(n_loops):
            dQ[i] = deltaQ(Q_loops_old[i], K_loops[i])
        Q_loops_new = apply_correction(loop_idexes, Q_loops_old, dQ)
        Q_loops_old = Q_loops_new
        residuals = np.sum(np.absolute(dQ))
        count += 1
    return Q_loops_old, dQ, residuals

print("""
#######################################################################################################################
Hybrid Piping Line Solver (Hardy Cross Method)\n
By: Alexander Perez and Joshua Nenadich\n
This program can solve a hybrid piping system with n symetrical loops. Each loop must have the same number of lines. The
User provides the type of pipe, the units of the volumetric flow rate Q, the number of loops in the system and the
number of lines per loop.

The User will also provide a file name of a csv (comma separated value) file that contains the following:
Column 1 --> Volumetric Flow Rates for each line
Column 2 --> Lengths of the Pipes in feet or meters
Column 3 --> Diameters of the pipes in feet or meters

There should be additional columns i.e. Column 4, 5, 6... that contain the lines that make up each loop
#######################################################################################################################
""")
run = True



while run:
    try:
        n = 1.852
        while True:
            # n = float(input("Enter a value for n: "))
            try:
                C = float(input("Enter a value for constant C for the type of Pipe: "))

            except ValueError:
                print("Please enter a valid decimal or integer value for C")
                # better try again... Return to the start of the loop
                continue

            else:
                break

        while True:
            try:
                k1 = str(input(
                    "Specify the units of the Volumetric Flow Rate Q, type cfs for ft^3/s, mgd for million gals/day or cms for m^3/s: "))

                if k1 == "cfs":
                    k1 = 4.727
                elif k1 == "mgd":
                    k1 = 10.63
                elif k1 == "cms":
                    k1 == 10.466
                else:
                    raise ValueError

            except ValueError:
                print("Please enter cfs, mgd or cms")

            else:
                break

        while True:
            try:
                n_loops = int(input("Enter the number of loops in the system: "))

            except ValueError:
                print("Please enter a valid integer for number of loops")
                # better try again... Return to the start of the loop
                continue

            else:
                break

        while True:
            try:
                n_lines = int(input("Enter the number of lines per loop: "))

            except ValueError:
                print("Please enter a valid integer for number of lines")
                # better try again... Return to the start of the loop
                continue

            else:
                break

        while True:
            try:
                file_name = str(input("Enter file name for pipe data. Ensure its in the same directory as this program: "))
                file_exists = os.path.exists(file_name)
                if not file_exists:
                    raise FileNotFoundError
                # print(file_exists)

            except FileNotFoundError:
                print("Please enter a valid file name and ensure it is in the same directory as this program")
                # better try again... Return to the start of the loop
                continue

            else:
                break

        # C = 100
        # k1 = 4.727
        # n_loops = 5
        # n_lines = 3
        #
        # file_name = "Piping_Network_Data.csv"
        L = np.loadtxt(file_name, delimiter=",", usecols=1, skiprows=1)
        D = np.loadtxt(file_name, delimiter=",", usecols=2, skiprows=1)
        Q = np.loadtxt(file_name, delimiter=",", usecols=0, skiprows=1)
        K = np.zeros_like(L)
        loop_idexes = [0 for i in range(n_loops)]
        for i in range(n_loops):
            col = int(4 + i)
            a = np.loadtxt(file_name, delimiter=',', usecols=col, skiprows=1, max_rows=n_lines)
            loop_idexes[i] = a

        loop_idexes = np.vstack(loop_idexes)

        for i in range(len(K)):
            Length = L[i]
            Diameter = D[i]
            a = (k1 * Length) / (C ** n * Diameter ** (4.8704))
            K[i] = float(a)

        K_loops = np.zeros_like(loop_idexes)
        Q_loops = np.zeros_like(loop_idexes)
        for j in range(len(loop_idexes)):
            for i in range(len(loop_idexes[0])):
                a = loop_idexes[j][i]
                K_loops[j][i] = K[int(a)]
                Q_loops[j][i] = Q[int(a)]

        Q_loops = Q_loops.tolist()
        K_loops = K_loops.tolist()
        loop_idexes = loop_idexes.tolist()
        copy_loop_idexes = deepcopy(loop_idexes)
        Q_loops = flip_repeated_values(copy_loop_idexes, Q_loops)
        Q_solution, dQ, Residuals = solve_hybrid_pipes(Q_loops, K_loops, loop_idexes)

        count = 0
        df_q_loops = []
        dat_dq = []
        dat_dq_names = []
        df_q_loops_indexes = []
        print("""
########################################################################################################################
Volumetric Flowrates:
        """)
        for i in range(len(loop_idexes)):
            dat_dq.append(dQ[i])
            dat_dq_names.append(f"dQ{i + 1}")
            for j in range(len(loop_idexes[0])):
                df_q_loops.append(Q_solution[i][j])
                df_q_loops_indexes.append(loop_idexes[i][j] + 1)
                print(f"Q{int(loop_idexes[i][j] + 1)} = {round(Q_solution[i][j], 3)}")

        dat = pd.DataFrame(df_q_loops, index=df_q_loops_indexes, columns=["Q[cfs]"])
        file_name_res = file_name.split(".", 1)[0]
        name_1 = file_name_res + "_results.csv"
        dat.to_csv(name_1)
        print(f"""
#######################################################################################################################
dQ values for each loop:
                """)

        for i in range(len(dat_dq)):
            print(f"Q{dat_dq_names[i]} = {dat_dq[i]}")

        dat2 = pd.DataFrame(dat_dq,index=dat_dq_names)
        name_2 = file_name_res + "_residuals.csv"
        dat2.to_csv(name_2)
        print(f"""
The residuals have been saved to a save file named {file_name_res}_residuals.csv
#######################################################################################################################
        """)
        a = input("Run again? (y/n): ")
        if a == "y":
            continue
        else:
            run = False

    except:
        print("It seems an error ocurred, please verify that the number of loops and lines entered mathces the "
              "data in the csv file and that the csv file is made as detailed in the instructions.")