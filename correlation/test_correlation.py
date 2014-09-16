import time
import numpy as np

# correlation kernel/function parameters
sigma = 0.8	# variance or sill
L = 5000 # correlation length

# grid parameters
nrow = 19
ncol = 33
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

t1= time.time()

C = np.zeros((nrow*ncol, nrow*ncol))
for i in range(nrow*ncol):
    for j in range(i, nrow*ncol):
        C[i,j] = py_correlation(coordinates[i], coordinates[j])
        C[j,i] = C[i,j]

t2 = time.time()

assert np.allclose(C,Cs)

print "=================================="
print "python took ", t2-t1, "seconds and "
print t1-t0, "seconds with swig correlation"
print "=================================="


for i in range(100):
	x1 = 100*np.random.rand()*np.random.rand(2)
	x2 = 100*np.random.rand()*np.random.rand(2)
	# print x1, x2,py_correlation(x1,x2),C_matrix.correlation(x1,x2)
	assert(py_correlation(x1,x2) ==  C_matrix.correlation(x1,x2))

print "=================================="
print "C_matrix.correlation seems to match pure python implementation"
print "=================================="




test_vecs = np.random.rand(nrow*ncol, 10)
test_outs = np.zeros((nrow*ncol, 10), 'd')
for i in range(np.shape(test_vecs)[1]):
	test_outs[:,i] = C.dot(test_vecs[:,i])



t0 = time.time()

C_matrix.make_C(nrow, ncol, Lx, Ly)

t1 = time.time()

out_vec = np.empty((nrow*ncol), 'd')
for i in range(np.shape(test_vecs)[1]):
	in_vec = test_vecs[:,i]
	C_matrix.av(nrow*ncol, in_vec, out_vec) # remember swig adds dim arg
	# print out_vec,test_outs[:,i]
	assert np.allclose(out_vec, test_outs[:,i])

t2 = time.time()

print "=================================="
print "sparse approach took ", t2-t1, "seconds and "
print t1-t0, "to set up"
print "sparse approach seems to match pure python implementation"
print "=================================="


t0 = time.time()

C_matrix.set_dims(nrow, ncol, Lx, Ly)

t1 = time.time()

out_vec = np.empty((nrow*ncol), 'd')
for i in range(np.shape(test_vecs)[1]):
	in_vec = test_vecs[:,i]
	C_matrix.av_no_C(nrow*ncol, in_vec, out_vec)
	assert np.allclose(out_vec,test_outs[:,i])

t2 = time.time()

print "=================================="
print "C not in memory approach took ", t2-t1, "seconds and "
print t1-t0, "to set up"
print "no memory approach seems to match pure python implementation"
print "=================================="

# first_significant = len(C)-10
# eig_vals, eig_vecs = scipy.linalg.eigh( C, turbo = True, eigvals = (first_significant, len(C)-1) )
