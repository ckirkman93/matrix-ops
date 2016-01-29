import numpy as np
import math
import sys

def LUdecomp():
    # Assumes an n x n or an n x n augmented matrix formatted
    # exactly as the examples provided to us, eg:
    # 1 2 4
    # 2 3 2
    # 5 6 7
    # or else it'll die and cry in that order.

    if len(sys.argv) > 1:
        file = sys.argv[1]
        with open(file) as f:
            lines = f.readlines()
            n = len(lines)
            A = np.zeros([n, n])
            b = np.zeros([n, 1])
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
                    else:
                        b[i][0] = row[j]

        # Assumes the second argument is a file containing 'b'
        if len(sys.argv) > 2:
            file = sys.argv[2]
            with open(file) as f:
                lines = f.readlines()
                n = len(lines)
                b = np.zeros([n, 1])
                for i in range(n):
                    entry = ''
                    j = 0
                    while j < len(lines[i]) and lines[i][j] != '\n':
                        entry += lines[i][j]
                        j += 1
                    b[i][0] = float(entry)

        print "A =\n", A
        print "b =\n", b
        solve_lu_b(A, b)
        solve_qr_b(A, b)

    else:
        for n in range(2, 13):
            A = get_pascal_matrix(n)
            b = get_b(n)

            print "**** n =", n
            print "P =\n", A
            print "b =\n", b
            solve_lu_b(get_pascal_matrix(n), b)
            solve_qr_b(get_pascal_matrix(n), b)
            print
            print
            print '===================='
            print
            print

def solve_lu_b(A, b):
    originalA = np.copy(A)
    L, U, error = lu_fact(A)
    numRows = A.shape[0]
    numColumns = A.shape[1]

    # First solve Ly = b for y.
    y = []
    for i in range(numRows):
        lhs = []
        # Makes perfect sense, promise
        for j in range(numColumns - (numColumns - i) + 1):
            lhs.append(L[i][j])
        for j in range(len(y)):
            lhs[j] = lhs[j] * y[j]
        constant = 0
        for j in range(len(lhs) - 1):
            constant += lhs[j]
        rhs = b[i]
        rhs = rhs - constant
        yi = rhs / lhs[i]
        y.append(yi[0])

    # Next, solve Ux = y for x.
    x = []
    for i in range(numRows):
        lhs = []
        for j in range(numColumns - (numColumns - i) + 1):
            lhs.insert(0, U[numRows - i - 1][numColumns - j - 1])
        for j in range(1, len(x) + 1):
            lhs[j] = lhs[j] * x[j - 1]
        constant = 0
        for j in range(1, len(lhs)):
            constant += lhs[j]
        rhs = y[numRows - i - 1]
        rhs = rhs - constant
        xi = rhs / lhs[0]
        x.insert(0, xi)

    x_array = np.array([])
    for el in x:
        x_array = np.append(x_array, [el])
    x_array = x_array[np.newaxis].T

    Ax = multMxN(originalA, x_array)
    solutionError = compute_error(Ax, b)

    print
    print "+++++ LU decomposition"
    print "*** L =\n", L
    print "*** U =\n", U
    print "***x solution =\n", x_array
    print "*** ||LU - P|| (error) =", error
    print "*** ||Px - b|| (error) =", solutionError

    return x_array

def solve_qr_b(A, b):
    Q = []
    R = []
    error = []

    Q1, R1, error1 = qr_fact_househ(A)
    Q2, R2, error2 = qr_fact_givens(A)

    Q.append(Q1)
    Q.append(Q2)
    R.append(R1)
    R.append(R2)
    error.append(error1)
    error.append(error2)

    numRows = A.shape[0]
    numColumns = A.shape[1]

    for z in range(2):
        y = multMxN(Q[z].transpose(), b)
        # Solve Rx = y for y
        x = []
        for i in range(numRows):
            lhs = []
            for j in range(numColumns - (numColumns - i) + 1):
                lhs.insert(0, R[z][numRows - i - 1][numColumns - j - 1])
            for j in range(1, len(x) + 1):
                lhs[j] = lhs[j] * x[j - 1]
            constant = 0
            for j in range(1, len(lhs)):
                constant += lhs[j]
            rhs = y[numRows - i - 1]
            rhs = rhs - constant
            xi = rhs / lhs[0]
            x.insert(0, xi)

        x_array = np.array([])
        for el in x:
            x_array = np.append(x_array, [el])
        x_array = x_array[np.newaxis].T

        Ax = multMxN(A, x_array)
        solutionError = compute_error(Ax, b)

        print
        if z == 0:
            print "+++++ QR decomposition (Householder Reflections)"
        else:
            print "+++++ QR decomposition (Givens Rotations)"
        print "*** Q =\n", Q[z]
        print "*** R =\n", R[z]
        print "*** x solution =\n", x_array
        print "*** ||QR - P|| (error) =", error[z]
        print "*** ||Px - b|| (error) =", solutionError

def lu_fact(A):
    # First check that A is n x n.
    # A.shape return the size of the matrix (n x m)
    if len(A.shape) != 2:
        return # Not a matrix
    elif A.shape[0] != A.shape[1]:
        return # Not an n x n
    n = A.shape[0]

    originalA = np.copy(A)

    # A[row][col]

    pivotPos = 0
    opMatrices = []
    L = np.eye(n)
    for i in range(n): # For each column
        # Puts the row with the largest entry of the ith column at the top of the matrix.
        # Todo: What if the entire ith column is zero?
        for j in range(i+1, n): # For each row/entry below the ith pivot
            if A[i][i] == 0:
                continue
            if A[j][i] != 0: # Then we want to make it zero
                opMatrix = np.eye(n) # Makes an n x n identity matrix (operation matrix)
                # xRi + Rj = 0; x = -Rj / Ri
                x = float(-A[j][i]) / A[i][i]
                opMatrix[j][i] = x
                L[j][i] = -x
                opMatrices.append(opMatrix)
                for k in range(n): # For each entry in the jth + 1 row
                    # Rj <= xRi + Rj
                    A[j][k] = x * A[i][k] + A[j][k]

    U = A

    #compute error
    LU = multNxN(L, U)
    error = compute_error(LU, originalA)

    return [L, U, error]

def qr_fact_househ(A):
    H_matrices = []
    numColumns = A.shape[0]
    Ai = np.copy(A)
    for n in range(numColumns-1):
        # This get the entries of the nth column, pivot and below
        v = Ai[n:,n]

        Hai = compute_Hai(v)

        # Build Hi
        if n > 0:
            Hi = np.eye(numColumns)
            for row in range(numColumns - (numColumns-len(v))):
                Hai_row = Hai[row]
                numLeadingZeros = numColumns - len(v)
                leadingZeroes = []
                for i in range(numLeadingZeros):
                    leadingZeroes.append(0)
                newRow = np.insert(Hai_row, 0, leadingZeroes) # Places zero(es) in front of the row
                Hi[row+(numColumns - len(v))] = newRow
        else:
            Hi = Hai

        H_matrices.append(Hi)
        Ai = multNxN(Hi, Ai)

    # Compute Q = H1 * H2 * ... * Hn
    Q = np.eye(numColumns)
    for Hi in H_matrices:
        Q = multNxN(Q, Hi)

    R = Ai

    QR = multNxN(Q, R)
    error = compute_error(QR, A)

    return[Q, R, error]

def qr_fact_givens(A):
    Gi_matrices = []
    numColumns = A.shape[0]
    numRows = numColumns # Assumes n x n matrix
    Ai = np.copy(A)
    for n in range(numColumns):
        for m in range(numRows-n):
            v = Ai[n:,n]
            x = v[0]
            xHeight = n

            allBelowPivotZero = True
            for i in range(1, len(v)):
                if v[i] != 0:
                    y = v[i]
                    yHeight = n + i
                    allBelowPivotZero = False
                    break
            if allBelowPivotZero:
                continue

            cos =  compute_cos(x, y)
            sin = compute_sin(x, y)
            width = yHeight - xHeight

            Gi = np.eye(numColumns)
            Gi[xHeight][n] = cos
            Gi[yHeight][n] = sin
            Gi[xHeight][n+width] = -sin
            Gi[yHeight][n+width] = cos
            Gi_matrices.append(Gi)

            Ai = multNxN(Gi, Ai)

    R = Ai

    #compute Q
    Q = np.eye(numColumns)
    for Gi in Gi_matrices:
        Gi_inv = Gi.transpose()
        Q = multNxN(Q, Gi_inv)

    #compute error
    QR = multNxN(Q, R)
    error = compute_error(QR, A)

    return[Q, R, error]

def get_pascal_matrix(n):
    pascal_matrix = np.zeros([n, n])
    for i in range(n): # for each row
        for j in range(n): # for each column
            pascal_matrix[i][j] = factorial(((i+1) - 1) + ((j+1) - 1)) / (factorial((i+1) - 1) * factorial((j+1) - 1))
    return pascal_matrix

def get_b(n):
    b = np.zeros([n, 1])
    for i in range(n):
        b[i][0] = 1./(i+1)
    return b

### For Givens' ###
def compute_cos(x, y):
    return 1.0*x / ((x**2 + y**2)**(.5))

def compute_sin(x, y):
    return -1.0*y / ((x**2 + y**2)**(.5))

### For Householder's ###
def compute_Hai(v):
    n = len(v)
    I = np.eye(n)
    u = compute_u(v)
    u_times_u = multMxN(u, u.T)
    Hi = I - (2.0/(magnitude(u)**2) * u_times_u)
    return Hi

# Expects v to be an 1D np.array
# Returns u used in Householder Reflections
def compute_u(v):
    n = len(v)
    magV = magnitude(v)

    e1 = np.array([1])
    for i in range(n-1):
        e1 = np.append(e1, [0])

    # Treating near zero value as zero
    if abs(v[0]) < 1.e-6:
        sign_v1 = 1
    else:
        sign_v1 = v[0] / abs(v[0])
    u = (v + (sign_v1 * magV * e1))
    return u[np.newaxis].T

### General use ###
def magnitude(v):
    mag = 0
    for i in range(len(v)):
        mag += v[i]**2
    return mag**(.5)

def normalize(v):
    return v / magnitude(v)

def compute_error(A, B):
    difference = abs(A - B)
    maxDifference = -float('inf') # Negative infinity
    for entry in np.nditer(difference):
        if entry > maxDifference:
            maxDifference = entry
    return maxDifference

def multNxN(A, B):
    n = A.shape[0]
    product = np.empty([n, n]) # Creating a new n x n matrix to store the product

    for i in range(n): # For each row i of A
        for j in range(n): # For each column j of B
            productEntry = 0
            for k in range(n): # For each entry in that row of A / column of B... if that makes sense
                productEntry += A[i][k] * B[k][j]
            product[i][j] = productEntry
    return product

def multMxN(A, B):
    if B.shape[0] != A.shape[1]:
        return None # Get outta here

    # (mA x n) , (n x mB)
    n = A.shape[1] # A's width / B's height
    mA = A.shape[0] # A's height
    mB = B.shape[1] # B's width

    # The resulting matrix will be mA x mB
    product = np.zeros((mA, mB))

    for i in range(mA): # For each row of A
        for j in range(mB): # For each column of B
            productEntry = 0
            for k in range(n): # For entry in that row/column pair
                product[i][j] += A[i][k] * B[k][j]
    return product

def invert(A):
    # Assumes A is n x n
    n = A.shape[0]

    Ainv = np.eye(n)
    for i in range(n): # For each column of A
        for j in range(i+1, n): # For each row/entry below the ith pivot
            if A[j][i] != 0:
                x = float(-A[j][i]) / A[i][i]
                for k in range(n):
                    A[j][k] = x * A[i][k] + A[j][k]
                    Ainv[j][k] = x * Ainv[i][k] + Ainv[j][k]
    # Now A is in Row Echelon Form
    for i in range(1, n): # For each column of A, except the first
        for j in range(0, i): # For each entry ABOVE the ith pivot
            if A[j][i] != 0:
                x = float(-A[j][i]) / A[i][i]
                for k in range(n):
                    A[j][k] = x * A[i][k] + A[j][k]
                    Ainv[j][k] = x * Ainv[i][k] + Ainv[j][k]
    # Now to divide the rows by the pivot entry...
    for i in range(n): # For each pivot entry
        if A[i][i] != 0 and A[i][i] != 1:
            d = A[i][i]
            for j in range(n):
                A[i][j] = A[i][j] / d
                Ainv[i][j] = Ainv[i][j] / d

    return Ainv

def factorial(n):
    result = 1
    for x in range(1, n + 1):
        result = result * x
    return result

def partialPivot(A, pivotPos):
    n = A.shape[0]
    maxEntryInCol = -float('inf')
    swapRow = None
    columnIsZero = True
    # For each entry in column 'i'
    for i in range(n - pivotPos):
        # Only check the entries below the pivot
        if abs(A[pivotPos+i][pivotPos]) > maxEntryInCol:
            maxEntryInCol = abs(A[pivotPos+i][pivotPos])
            swapRow = pivotPos+i
            if A[pivotPos+i][pivotPos] != 0:
                columnIsZero = False

    if columnIsZero:
        return None

    # If this is false, the largest entry in the column is
    # already in the pivot position -> no need to swap rows
    if swapRow is not None and swapRow != pivotPos:
        # copies A[swapRow] into temp
        temp = list(A[swapRow])
        A[swapRow] = A[pivotPos]
        A[pivotPos] = temp
    return A

LUdecomp()
