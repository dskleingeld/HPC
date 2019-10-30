import scipy.io
import scipy.linalg
import numpy as np
from numpy.linalg import inv

import sys

def printf(format, *args):
    sys.stdout.write(format % args)

def dump_non_zeros(array):
    for row_n, row in enumerate(array):
        for col_n, element in enumerate(row):
            if element == 0:
                continue
            
            printf("%d %d %.6f\n", row_n, col_n, element)

def print_C_array(array):
    for element in array:
        printf(",%.6f", element)

mat = scipy.io.mmread("../src/data/mcfe.mtx").toarray()
#mat = scipy.io.mmread("../src/tests/sparse").toarray()

pattern = np.array([5.0, -5.0])
length = mat.shape[0]

pat = np.tile(pattern, (length+1)//2)
pat = pat[0:length]
#ones = np.ones(length)

b = mat.dot(pat)

print(b)
#print(pat)

lu, piv = scipy.linalg.lu_factor(mat);
P, L, U = scipy.linalg.lu(mat);

np.set_printoptions(threshold=sys.maxsize)
#print(P[0:20,0:20])

##test algo
b_p = inv(P).dot(b)
#print(inv(P)[0:20,0:20])
#print(piv[0:20])
#the transpose of a permutation matrix P is its inverse.

#print_C_array(b_p)

p = np.arange(0, len(piv))
b2 = np.copy(b)
for i,v in enumerate(piv):
    v = piv[i]
    temp = p[i]
    p[i] = p[v]
    p[v] = temp

#print(p[0:20])

#for i,c in enumerate(p):
#    b2[c] = b[i]

for i in range(len(b2)):
    b2[i] = b[p[i]]

error = False
for ba, bb in zip(b_p,b2):
    if ba != bb:
        error = True

#print(error)
#print(p[:20])
#print(b2)
#print(b_p[:20])


#dump_non_zeros(P)
#dump_non_zeros(lu)

#print(P@b)
#print(piv);


#x= scipy.linalg.lu_solve((lu,piv),b)
#print(x)