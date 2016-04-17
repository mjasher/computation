"""
preprocess, compile, and link
"""

import subprocess
# compile_cmd = 'gfortran -shared -fPIC  -g -o local.so local.f90'
compile_cmd = 'gcc -shared -fPIC  -g -o local.so local.c'
subprocess.check_call(compile_cmd.split())


from ctypes import c_char_p, c_long, c_float, c_double, byref, CDLL, POINTER
local = CDLL('./local.so')

# test.prnt_.argtypes = [c_char_p, c_long]

# local.sin_degrees.argtypes = [POINTER(c_double)]
# local.sin_degrees.restype = POINTER(c_double)
# sin_x = local.sin_degrees(byref(x))


local.sin_degrees.argtypes = [c_double]
local.sin_degrees.restype = c_double

x = c_double(4.56)
sin_x = local.sin_degrees(x)
print ""
print sin_x


# test.qref_.argtypes = [POINTER(c_float), c_float, c_long, c_long, c_long]
# test.qref_.restype = c_float

# import numpy 

# x = 4
# y = 3
# z = 2
# a = numpy.random.random((x,y,z))
# print a
# print "SUM", numpy.sum(a)

# a_c_types = (c_float*(x*y*z))()
# a_c_types[:] = a.flatten()
# tot = c_float()


# test.qref_( byref(a_c_types), byref(tot), byref(c_long(z)), byref(c_long(y)), byref(c_long(x)) )
# print "SUM", tot

# # c =[0.501212060,8.64220131E-03,0.589656591,0.670281947,0.714529157,0.222474501,0.289629757,0.858183444,2.30114143E-02,0.602795005,0.404510617,0.784162641,0.283252627,6.63019419E-02,0.695632756,0.633914351,3.02934125E-02,0.368286878,4.63037528E-02,0.342094421,0.849391580,0.408671349,0.758220494,0.166015476]



# import numpy as np


# ncol = 4
# nrow = 3
# nlay = 1
# ibound = np.ones((ncol,nrow,nlay), dtype=np.long, order="F")
# strt = 10.0*np.ones((ncol,nrow,nlay), dtype=np.float64, order="F")
# hnoflo = -9999.0
# name = '../tutorial2/tutorial2'

# strt[:,:1,:] = 5.

# # np.float64 corresponds with REAL(8)
# # np.float32 corresponds with REAL
# # x.ctypes.data_as(ctypes.POINTER(ctypes.c_long))

# test.halfhalf_.argtypes = [
#                             POINTER(c_float), 
#                             POINTER(c_long), POINTER(c_long), POINTER(c_long), 
#                             np.ctypeslib.ndpointer(dtype=np.float64,
#                                                 # ndim=3,
#                                                  shape=(ncol,nrow,nlay),
#                                                  flags='F_CONTIGUOUS'),
#                             np.ctypeslib.ndpointer(dtype=np.long,
#                                                 # ndim=3,
#                                                  shape=(ncol,nrow,nlay),
#                                                  flags='F_CONTIGUOUS'),
#                             c_char_p, c_long
#                             ]

# test.halfhalf_(
#     c_float(hnoflo), 
#     c_long(ncol), c_long(nrow), c_long(nlay), 
#     strt,
#     ibound,
#     name, len(name)
#     )
