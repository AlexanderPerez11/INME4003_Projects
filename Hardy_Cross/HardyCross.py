import numpy as np

L1 = 2000.0
L2 = 2000.0
L3 = 2000.0
L4 = 2000.0
L5 = 2000.0
L6 = 2000.0
L7 = 2828.0
L8 = 2000.0
L9 = 2000.0
L10 = 2000.0
L11 = 2000.0
L12 = 2000.0

D1 = 6
D2 = 6
D3 = 8
D4 = 8
D5 = 8
D6 = 8
D7 = 12
D8 = 6
D9 = 6
D10 = 6
D11 = 6
D12 = 6

L = np.array([L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12])
D = np.array([D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12])

K = np.zeros_like(L)

n = 1.852
C = 100
k1 = 4.727

Q1 = 5.0
Q4 = 2.0
Q6 = 2.0
Q8 = 6.0
Q12 = 5.0
Q2 = 6 + Q1
Q5 = 4 - Q4 - Q8
Q10 = Q12 - 4
Q11 = Q12 + 4
Q9 = 2 +Q8 + Q11 - Q10
Q7 = 3+Q6 -Q5-Q9
Q3 = Q1 - Q4-Q7 -3
Q = np.array([Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12])
print(Q)
n_loops = 5
loop1_index = [0, 1, 2]
loop2_index = [2, 5, 6]
loop3_index = [3, 4, 6]
loop4_index = [4,7, 8]
loop5_index = [9, 10, 11]
loop_idexes = [loop1_index, loop2_index, loop3_index, loop4_index, loop5_index]

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

    # Return the modified matrix A
    ones_matrix = [[1 if val >= 0 else -1 for val in row] for row in A]
    flipped = np.multiply(ones_matrix,B)
    return flipped


def find_shared_rows(matrix, n):
    # create a set to keep track of the shared values
    shared_values = set(matrix[n])

    # find the index of the rows that share a value with the fixed row index n
    shared_rows = []
    for i, row in enumerate(matrix):
        if i == n:
            continue
        if any(value in shared_values for value in row):
            shared_rows.append(i)

    # return the result
    if shared_rows:
        return shared_rows
    else:
        return None


Q_loop1 = [Q[i] for i in loop1_index]
Q_loop2 = [Q[i] for i in loop2_index]
Q_loop3 = [Q[i] for i in loop3_index]
Q_loop4 = [Q[i] for i in loop4_index]
Q_loop5 = [Q[i] for i in loop5_index]
Q_loops = [Q_loop1, Q_loop2, Q_loop3, Q_loop4, Q_loop5]
print(Q_loops)
Q_loops = flip_repeated_values(loop_idexes, Q_loops)
print(Q_loops)
for i in range(len(K)):
    a = (k1 * L[i]) / ((C ** (n)) * ((D[i]/12) ** (4.8704)))
    print(a)
    K[i] = float(a)

K_loop1 = [float(K[i]) for i in loop1_index]
K_loop2 = [float(K[i]) for i in loop2_index]
K_loop3 = [float(K[i]) for i in loop3_index]
K_loop4 = [float(K[i]) for i in loop4_index]
K_loop5 = [float(K[i]) for i in loop5_index]

K_loops = np.array([K_loop1, K_loop2, K_loop3, K_loop4, K_loop5])
print(K_loops)

def deltaQ(Q, K):
    sum_numerator = np.zeros_like(Q)
    sum_denominator = np.zeros_like(Q)
    for i in range(len(sum_denominator)):
        sum_numerator[i] = -K[i] * Q[i] * abs(Q[i]) ** (n - 1)
        sum_denominator[i] = K[i] * abs(Q[i]) ** (n - 1)

    return np.sum(sum_numerator) / (n * np.sum(sum_denominator))


Q_loops = [Q_loop1, Q_loop2, Q_loop3, Q_loop4, Q_loop5]
K_loops = [K_loop1, K_loop2, K_loop3, K_loop4, K_loop5]
print(K_loops)

def solve_hybrid_pipes(Q_loops, tol=1e-7, max_iter=100):
    Q_loops_old = Q_loops.copy()
    Q_loops_new = np.zeros_like(Q_loops)
    dQ = [1, 1, 1, 1, 1]

    count = 0
    residuals = np.sum(np.absolute(dQ))
    #abs(dQ[0]) > tol and abs(dQ[1]) > tol and abs(dQ[2]) > tol and abs(dQ[3]) > tol and abs(dQ[4]) > tol
    while residuals> tol:
        for i in range(n_loops):
            dQ[i] = deltaQ(Q_loops_old[i], K_loops[i])
        for k in range(n_loops):
            correction_indexes = find_shared_rows(loop_idexes, k)
            if correction_indexes == None:
                shared_correction = 0
            else:
                shared_corrections = []
                for i in correction_indexes:
                    shared_corrections.append(dQ[i])
                shared_correction = sum(shared_corrections)

            for j in range(len(Q_loops[k])):

                Q_loops_new[k][j] = Q_loops_old[k][j] + dQ[k] - shared_correction
        residuals = np.sum(np.absolute(dQ))
        Q_loops_old = Q_loops_new

        count+=1


    return Q_loops_old, dQ, residuals


a, b, c = solve_hybrid_pipes(Q_loops)
print(Q_loops)
print(a, b, c)
