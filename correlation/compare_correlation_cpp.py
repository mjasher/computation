import time
import numpy as np

# correlation kernel/function parameters
sigma = 0.8	# variance or sill
L = 5000 # correlation length

# grid parameters
nrow = 300
ncol = 30
Lx = 1000.
Ly = 1000.

# initialise swig wrapped cpp class
import correlation
C_matrix = correlation.C_matrix(sigma, L)

test_vecs = np.random.rand(nrow*ncol, 10)
out_vec = np.empty((nrow*ncol), 'd')


t0 = time.time()

C_matrix.set_dims(nrow, ncol, Lx, Ly)

t1 = time.time()

for i in range(np.shape(test_vecs)[1]):
	in_vec = test_vecs[:,i]
	C_matrix.av_no_C(nrow*ncol, in_vec, out_vec)

t2 = time.time()

print "=================================="
print "C not in memory approach took ", t2-t1, "seconds and "
print t1-t0, "to set up"
print "=================================="

# t0 = time.time()
#
# C_matrix.make_C(nrow, ncol, Lx, Ly)
#
# t1 = time.time()
#
# for i in range(np.shape(test_vecs)[1]):
# 	in_vec = test_vecs[:,i]
# 	C_matrix.av(nrow*ncol, in_vec, out_vec) # remember swig adds dim arg
# 	# print out_vec,test_outs[:,i]
#
# t2 = time.time()
#
# print "=================================="
# print "sparse approach took ", t2-t1, "seconds and "
# print t1-t0, "to set up"
# print "=================================="



# exmample results
# ==============================
# sigma = 0.8	# variance or sill
# L = 5000 # correlation length
#
# # grid parameters
# nrow = 300
# ncol = 300
# Lx = 1000.
# Ly = 1000.
#
# no openMP
# ==================================
# C not in memory approach took  1081.95444584 seconds (20min) and
# 1.21593475342e-05 to set up
# ==================================
# ==================================
# sparse approach took  44738.708575 seconds and
# 115.018224955 to set up
# ==================================

# with openMP on a 4 core desktop
# ==================================
# C not in memory approach took  416.903661966 seconds and
# 1.09672546387e-05 to set up
# ==================================
