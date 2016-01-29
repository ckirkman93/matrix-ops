import numpy as np
import random

from util import *

A = np.array([[1., 1./2, 1./3], [1./2, 1., 1./4], [1./3, 1./4, 1.]])
b = np.array([[.1], [.1], [.1]])
x_exact = np.array([[9./190], [28./475], [33./475]])
eps = 0.00005
M = 100

def iterativeMethod():
    x0_list = []

    jacobi_error_list = []
    gs_error_list = []
    jacobi_N_list = []
    gs_N_list = []

    for i in range(100):
        x0 = np.array([[random.uniform(-1, 1)], [random.uniform(-1, 1)], [random.uniform(-1, 1)]])
        x0_list.append(x0)

    # part 2b
    jacobi_results = []
    gs_results = []
    for x0 in x0_list:
        x_approx1 = jacobi_iter(x0, eps, M)
        jacobi_results.append(x_approx1)
        x_approx2 = gs_iter(x0, eps, M)
        gs_results.append(x_approx2)

    # part 2c
    jacobi_average = np.array([[0], [0], [0]])
    jacobi_N = 0

    for x0, approx, N in jacobi_results:
        jacobi_average = np.add(approx, jacobi_average)
        jacobi_N += N
        jacobi_error_list.append(compute_error(x0, x_exact))
        jacobi_N_list.append(N)
    for i in range(len(jacobi_average)):
        jacobi_average[i] = jacobi_average[i] / 100.
    jacobi_N /= 100
    jacobi_error = compute_error(jacobi_average, x_exact)

    print 'Jacobi Results:'
    print
    print 'Average Result'
    print jacobi_average
    print
    print 'Average Numbers of Iterations'
    print jacobi_N
    print
    print 'Error'
    print jacobi_error

    gs_average = np.array([[0], [0], [0]])
    gs_N = 0

    for x0, approx, N in gs_results:
        gs_average = np.add(approx, gs_average)
        gs_N += N
        gs_error_list.append(compute_error(x0, x_exact))
        gs_N_list.append(N)
    for i in range(len(gs_average)):
        gs_average[i] = gs_average[i] / 100.
    gs_N /= 100

    gs_error = compute_error(gs_average, x_exact)

    print
    print "======="
    print
    print 'GS Results:'
    print
    print 'Average Result'
    print gs_average
    print
    print 'Average Numbers of Iterations'
    print gs_N
    print
    print 'Error'
    print gs_error

    print
    print '================='
    print

    stepsRatio = (1.0)*jacobi_N / gs_N

    print
    print 'Steps Ratio (Jacobi : GS)'
    print stepsRatio

def jacobi_iter(x0, eps, M):
    x_n = x0
    # Assumes A is an n x n matrix
    n = A.shape[0]

    S_inv = np.zeros([n, n])
    for i in range(n):
        if A[i][i] == 0:
            S_inv[i][i] = 0
        else:
            S_inv[i][i] = 1. / A[i][i]

    LU = np.copy(A)
    for i in range(n):
        for j in range(n):
            if i == j:
                LU[i][j] = 0

    S_inv_LU = multMxN(S_inv, LU) * -1
    S_inv_b = multMxN(S_inv, b)

    N = 0
    hitMax = True
    while N < M:
        N += 1
        x_n_minus1 = x_n
        x_n = multMxN(S_inv_LU, x_n) + S_inv_b
        error = compute_error(x_n, x_n_minus1)
        if error <= eps:
            hitMax = False
            break

    if hitMax:
        return None

    return [x0, x_n, N]

def gs_iter(x0, eps, M):
    x_n = x0
    # Assumes A is an n x n matrix
    n = A.shape[0]

    S = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            if i >= j:
                S[i][j] = A[i][j]

    # Hella hard-coded
    S_inv = np.zeros([n, n])
    S_inv[0][0] = 1/S[0][0]
    S_inv[1][0] = -S[1][0]/(S[0][0]*S[1][1])
    S_inv[1][1] = 1/S[1][1]
    S_inv[2][0] = ((-S[1][1]*S[2][0])+(S[1][0]*S[2][1]))/(S[0][0]*S[1][1]*S[2][2])
    S_inv[2][1] = -S[2][1]/(S[1][1]*S[2][2])
    S_inv[2][2] = 1/S[2][2]

    U = np.copy(A)
    for i in range(n):
        for j in range(n):
            if i >= j:
                U[i][j] = 0

    S_inv_U = multMxN(S_inv, U) * -1
    S_inv_b = multMxN(S_inv, b)

    N = 0
    hitMax = True
    while N < M:
        N += 1
        x_n_minus1 = x_n
        x_n = multMxN(S_inv_U, x_n) + S_inv_b
        error = compute_error(x_n, x_n_minus1)
        if error <= eps:
            hitMax = False
            break

    if hitMax:
        return None

    return [x0, x_n, N]

iterativeMethod()
