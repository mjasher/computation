import time
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh

# correlation kernel/function parameters
sigma = 0.8	# variance or sill
L = 5000 # correlation length

# grid parameters
nrow = 109
ncol = 63
Lx = 1000.
Ly = 900.

# initialise swig wrapped cpp class
import correlation
C_matrix = correlation.C_matrix(sigma, L)

# equivalent python functions, including timing
def py_correlation(x1, x2):
    N = len(x1)
    L_array = np.array([ L for i in range(N)]) # correlation lengths (ND)
    return sigma**2 * np.exp( -N*np.sum((x1-x2)**2 / L_array**2) )

x = np.linspace( Lx/2, Lx*(nrow-0.5), nrow )
y = np.linspace( Ly/2, Ly*(ncol-0.5), ncol )
X,Y = np.meshgrid( x, y )
coordinates = np.dstack((np.rollaxis(X,1).flatten(), np.rollaxis(Y,1).flatten()))[0]


t0= time.time()

Cs = np.zeros((nrow*ncol, nrow*ncol))
for i in range(nrow*ncol):
    for j in range(i, nrow*ncol):
        Cs[i,j] = C_matrix.correlation(coordinates[i], coordinates[j])
        Cs[j,i] = Cs[i,j]

# C = np.zeros((nrow*ncol, nrow*ncol))
# for i in range(nrow*ncol):
#     for j in range(i, nrow*ncol):
#         C[i,j] = py_correlation(coordinates[i], coordinates[j])
#         C[j,i] = C[i,j]
# assert np.allclose(C,Cs)

t1= time.time()

num_eig = 20
eigh_range = ((Cs.shape[0]-num_eig+0), Cs.shape[0]-1)
eig_vals, eig_vecs = eigh( Cs, turbo = True, eigvals = eigh_range )
# e,v = eigh( C )

t2 = time.time()


print "=================================="
print "python (with c++ kernel function) took ", t2-t1, "seconds for SVD and "
print t1-t0, "seconds to build C"
print "=================================="


t0 = time.time()

C_matrix.set_dims(nrow, ncol, Lx, Ly)

out_vec = np.zeros((nrow*ncol), 'd')
def Av(v):
	C_matrix.av_no_C(nrow*ncol, v, out_vec)
	return out_vec

A = LinearOperator((nrow*ncol,nrow*ncol), matvec=Av, dtype='d')
# A.matvect(scipy.ones(nrow*ncol))

t1 = time.time()

s_eig_vals, s_eig_vecs = eigsh( A, k=num_eig)

t2 = time.time()

assert np.allclose(eig_vals, s_eig_vals)
# sometimes one gets -v rather than v?
for i in range(num_eig):
	assert np.allclose(eig_vecs[:,i], s_eig_vecs[:,i]) or np.allclose(eig_vecs[:,i], -1*s_eig_vecs[:,i])


print "=================================="
print "C 'not in memory' SVD took ", t2-t1, "seconds and "
print t1-t0, "seconds to set up"
print "no memory approach seems to match pure python implementation"
print "=================================="



# test results 
# =====================
# correlation kernel/function parameters
# sigma = 0.8	# variance or sill
# L = 5000 # correlation length

# # grid parameters
# nrow = 109
# ncol = 63
# Lx = 1000.
# Ly = 900.
# ==================================
# python (with c++ kernel function) took  102.93483901 seconds for SVD and 
# 36.87987113 seconds to build C
# ==================================
# ==================================
# C 'not in memory' SVD took  87.4241662025 seconds and 
# 8.08238983154e-05 seconds to set up
# no memory approach seems to match pure python implementation
# ==================================
