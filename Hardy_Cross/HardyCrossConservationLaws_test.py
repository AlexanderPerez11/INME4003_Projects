import numpy as np
from copy import deepcopy

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

n = 1.852
C = 100
k1 = 4.727

L = np.loadtxt("Piping_Network_Data_test2.csv", delimiter=",", usecols=1, skiprows=1)
D = np.loadtxt("Piping_Network_Data_test2.csv", delimiter=",", usecols=2, skiprows=1)
Q = np.loadtxt("Piping_Network_Data_test2.csv", delimiter=",", usecols=0, skiprows=1)
# print(L)
K = np.zeros_like(L)

n_loops = 2
n_lines = 3
loop_idexes = [0 for i in range(n_loops)]
# print(loop_idexes)
for i in range(n_loops):
    col = int(4+i)
    a = np.loadtxt("Piping_Network_Data_test2.csv", delimiter=',', usecols=col, skiprows=1, max_rows=n_lines)
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
        print(a)
        K_loops[j][i] = K[int(a)]
        Q_loops[j][i] = Q[int(a)]

Q_loops = Q_loops.tolist()
K_loops =  K_loops.tolist()
loop_idexes = loop_idexes.tolist()
copy_loop_idexes = deepcopy(loop_idexes)
Q_loops = flip_repeated_values(copy_loop_idexes,  Q_loops)
Q_solution, dQ, Residuals = solve_hybrid_pipes(Q_loops, K_loops, loop_idexes)
print(Q_solution)
