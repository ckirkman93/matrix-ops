import numpy as np
import random
import sys
from util import *

def powerMethod():
    ## change these for testing! ##
    v = np.array([[1], [1]])
    eps = 0.00005
    N = 100

    onehundredclub = 0
    trace_list = []
    trace_inv_list = []
    det_list = []
    det_inv_list = []
    n_list = []
    n_inv_list = []

    if len(sys.argv) > 1:
        file = sys.argv[1]
        with open(file) as f:
            lines = f.readlines()
            n = len(lines)
            A = np.zeros([n, n])
            for i in range(n):
                counter = 0
                k = 0
                row = []
                entry = ''
                while k < len(lines[i]) and lines[i][k] != '\n':
                    if lines[i][k] != ' ':
                        entry += lines[i][k]
                    else:
                        row.append(int(entry))
                        entry = ''
                    k += 1
                row.append(float(entry))
                for j in range(len(row)):
                    if j < n:
                        A[i][j] = row[j]


    if len(sys.argv) > 1:
        it = 1
    else:
        it = 1000
    for i in range(it):
        A = np.zeros([2, 2])
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                A[i][j] = random.uniform(-2, 2)
        det = (A[0][0] * A[1][1]) - (A[0][1] * A[1][0])
        if det == 0:
            i -= 1
            continue

        # Building A inverse
        A_inv = np.zeros([2, 2])
        a = A[0][0]
        A_inv[0][0] = A[1][1]
        A_inv[1][1] = a
        A_inv[0][1] = -A[0][1]
        A_inv[1][0] = -A[1][0]
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                A_inv[i][j] = A_inv[i][j] * (1.0/det)
        det_inv = (A_inv[0][0] * A_inv[1][1]) - (A_inv[0][1] * A_inv[1][0])

        print 'A =\n', A
        print
        print 'A inverse =\n', A_inv
        print

        result = power_method(A, v, eps, N)
        if result != None:
            value, vector, numIterations = result
            n_list.append(numIterations)

            print 'largest eigenvalue:'
            print value
            print 'accompanying eigenvector:'
            print vector
            print 'number of iterations:'
            print numIterations

        else:
            #print None
            n_list.append(100)
            onehundredclub += 1

        #print

        result = power_method(A_inv, v, eps, N)
        if result != None:
            value_inv, vector_inv, numIterations_inv = result
            n_inv_list.append(numIterations_inv)

            print 'smallest eigenvalue:'
            print value_inv
            print 'accompanying eigenvector:'
            print vector_inv
            print 'number of iterations:'
            print numIterations_inv

        else:
            #print None
            n_inv_list.append(100)

        tr = trace(A)
        tr_inv = trace(A_inv)
        trace_list.append(tr)
        trace_inv_list.append(tr_inv)

        print
        print 'trace of A / trace of A inverse'
        print tr, "/", tr_inv
        print

        det_list.append(det)
        det_inv_list.append(det_inv)

        print 'determinant of A / determinant of A inverse'
        print det, "/", det_inv
        print
        print '==============='
        print

    """
    for d in det_list:
        print d
    print '====='
    for t in trace_list:
        print t
    print '====='
    for n in n_list:
        print n
    print '====='
    print '====='
    for d_i in det_inv_list:
        print d_i
    print '====='
    for t_i in trace_inv_list:
        print t_i
    print '====='
    for n_i in n_inv_list:
        print n_i
    """

def power_method(A, v, eps, N):
    wt = np.array([[1, 0]])
    v_list = [v]
    eigenvalues = []
    eigenvectors = []
    numIterations = 0
    while numIterations < N:
        numIterations += 1
        eigenvectors.append(v / magnitude(v))

        v = multMxN(A, v)
        v_list.append(v)

        # Eigenvalues
        prev = v_list[-2] # Gets the second to last entry in v_list
        current = v_list[-1]
        eigenvalue = np.dot(wt, current) / np.dot(wt, prev)
        eigenvalues.append(eigenvalue)

        if len(eigenvalues) > 1:
            if abs(eigenvalues[-1] - eigenvalues[-2]) < eps:
                break

    if numIterations >= N:
        return None
    return [eigenvalues[-1][0][0], eigenvectors[-1], numIterations]

powerMethod()
