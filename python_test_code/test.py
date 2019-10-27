import scipy.io
import scipy.linalg
import numpy as np

import sys

def printf(format, *args):
    sys.stdout.write(format % args)

def dump_non_zeros(array):
    for row_n, row in enumerate(array):
        for col_n, element in enumerate(row):
            if element == 0:
                continue
            
            printf("%d %d %.6f\n", row_n, col_n, element)

mat = scipy.io.mmread("../src/data/mcfe.mtx").toarray()
ones = np.ones(mat.shape[0])
b = mat.dot(ones)

#print(b)

lu, piv = scipy.linalg.lu_factor(mat);
P, L, U = scipy.linalg.lu(mat);

np.set_printoptions(threshold=sys.maxsize)
#print(P[0:20,0:20])

##test algo
b_p = P.dot(b)

p = np.arange(0, len(piv))
for i,v in enumerate(piv):
    v = piv[i]
    temp = p[i]
    p[i] = p[v]
    p[v] = temp

b2 = np.copy(b)
for i,c in enumerate(p):
    b2[c] = b[i]

error = False
for ba, bb in zip(b_p,b2):
    if ba != bb:
        error = True

print(error)
print(p[:20])
print(b2[:20])
print(b_p[:20])


#dump_non_zeros(P)
#dump_non_zeros(lu)

#print(P@b)
#print(piv);
