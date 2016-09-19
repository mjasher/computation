import numpy
import subprocess

import ctypes as C

"""
Pre-process, compile, and link
Fortran can be used similarly 
    - `gfortran -shared -fPIC  -g -o local.so local.f90`
    - add _ to end of function names, `local.sin_x_`
    - use `[numpy.ctypeslib.ndpointer(dtype=numpy.float64, shape=(n,), flags='F_CONTIGUOUS'), x.astype(dtype=numpy.float64, order="F")]`
"""

compile_cmd = 'gcc -shared -fPIC  -g -o local.so local.c'  
subprocess.check_call(compile_cmd.split())

local = C.CDLL('./local.so')

"""
simplest example, args and ret by value
"""

local.sin_degrees.argtypes = [C.c_double]
local.sin_degrees.restype = C.c_double
x = 4.56
sin_x = local.sin_degrees(C.c_double(x))
assert abs(sin_x - numpy.sin(x * numpy.pi / 180.)) < 1e-10
print 'c:        ', sin_x
print 'numpy:    ', numpy.sin(x * numpy.pi / 180.)

"""
args and ret by reference
"""

local.sin_degrees_p.argtypes = [C.POINTER(C.c_double)]  # redundant?
local.sin_degrees_p.restype = C.POINTER(C.c_double)
sin_x_p = local.sin_degrees_p(C.byref(C.c_double(x)))
assert abs(sin_x_p[0] - numpy.sin(x * numpy.pi / 180.)) < 1e-10
print 'c pointer:', sin_x


"""
numpy array args

http://scipy.github.io/old-wiki/pages/Cookbook/Ctypes
"""
ncol = 4
nrow = 3
rand_array = 10.0 * numpy.random.random((ncol, nrow))

local.sum.restype = C.c_double

# complex way
"""
sum_types, sum_args = zip(*[
                        [numpy.ctypeslib.ndpointer(dtype=numpy.float64, shape=rand_array.shape), 
                        rand_array.astype(dtype=numpy.float64)],
                        [C.c_int, C.c_int(ncol)], 
                        [C.c_int, C.c_int(nrow)], 
                        # [C.C_FUNC, C.C_FUNC(f)]
                        ])

local.sum.argtypes = sum_types
c_sum = local.sum(* sum_args)
"""

# easy way
c_sum = local.sum(rand_array.ctypes.data_as(C.POINTER(C.c_double)), C.c_int(ncol), C.c_int(nrow))

print "(1,2) = ", rand_array[1,2]
assert abs(c_sum - numpy.sum(rand_array)) < 1e-10



