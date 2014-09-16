import time
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh

# correlation kernel/function parameters
sigma = 0.8	# variance or sill
L = 5000 # correlation length

# grid parameters
nrow = 19
ncol = 63
Lx = 1000.
Ly = 900.

num_eig = 10

# initialise swig wrapped cpp class
import correlation
C_matrix = correlation.C_matrix(sigma, L)


C_matrix.set_dims(nrow, ncol, Lx, Ly)

out_vec = np.zeros((nrow*ncol), 'd')
def Av(v):
	C_matrix.av_no_C(nrow*ncol, v, out_vec)
	return out_vec

A = LinearOperator((nrow*ncol,nrow*ncol), matvec=Av, dtype='d')
# A.matvect(scipy.ones(nrow*ncol))

t0 = time.time()
eig_vals, eig_vecs = eigsh( A, k=num_eig)
t1 = time.time()

# save 
np.savez('kle_eigen.npz', eig_vals=eig_vals, eig_vecs=eig_vecs)

# load
npzfile = np.load('kle_eigen.npz')
eig_vals = npzfile['eig_vals']
eig_vecs = npzfile['eig_vecs']


# how you'd use it 
##########################################

# mu = np.zeros(nrow*ncol) #np.mean([ p['hk'] for p in observed_points ]) * np.ones(len(C))  
# scale = 1.
# # evaluate KLE with given modes
# # other parameters (other than modes) are mu, scale, eig_vecs, eig_vals (or really sigma and L on which they depend)
# def KLE(modes):
#     coefs = np.sqrt(eig_vals) * modes # elementwise
#     truncated_M = mu + scale * np.dot( eig_vecs, coefs)
#     unflattened = truncated_M.reshape(nrow,ncol)
#     return unflattened

# modes = np.ones(num_eig)
# print KLE(modes)

############################################


print "=================================="
print "done in ", t1-t0, "seconds"
print "=================================="


for i in range(num_eig):
	assert np.allclose(A.matvec(eig_vecs[:,i]), eig_vals[i]*eig_vecs[:,i]) 

print "=================================="
print "they are indeed eigenvectors"
print "=================================="
