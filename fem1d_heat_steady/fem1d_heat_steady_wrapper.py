"""

Wrapper for [John Burkardt's 1D FEM steady state heat equation code](https://people.sc.fsu.edu/~jburkardt/f_src/fem1d_heat_steady/fem1d_heat_steady.html)

Serves as a demo of using python ctypes with fortran.

Note that conductivity k and forcing function f specified in python must treat arguments as pointers/arrays.

michael.james.asher@gmail.com
August 2015


"""

import numpy 
import subprocess
import os
import ctypes as C

"""
compile
"""

base = os.path.abspath(os.path.join(os.path.dirname(__file__), 'fem1d_heat_steady'))
compile_cmd = 'gfortran -shared -fPIC  -g -o '+base+'.so '+base+'.f90'
subprocess.check_call(compile_cmd.split())

"""
wrapping function
"""

fem1d_heat_steady_so = C.CDLL(base+'.so')

def fem1d_heat_steady(	n, 					# number of mesh points
						a, b, ua, ub, 		# u(a)=ua on left boundary, u(b)=ub on right boundary
						k, 					# conductivity function 
						f, 					# forcing function
						x 					# mesh points
					 ):	

	C_FUNC = C.CFUNCTYPE(C.c_double, C.POINTER(C.c_double))

	u = numpy.zeros((n), dtype=numpy.float64, order="F")

	fem1d_heat_steady_types, fem1d_heat_steady_args = zip(*[
	 										[C.POINTER(C.c_int), C.c_int(n)], 
											[C.POINTER(C.c_double), C.c_double(a)], 
											[C.POINTER(C.c_double), C.c_double(b)], 
											[C.POINTER(C.c_double), C.c_double(ua)], 
											[C.POINTER(C.c_double), C.c_double(ub)], 
											[C_FUNC, C_FUNC(k)], 
											[C_FUNC, C_FUNC(f)], 
											[numpy.ctypeslib.ndpointer(dtype=numpy.float64, shape=(n,), flags='F_CONTIGUOUS'), x.astype(dtype=numpy.float64, order="F")],
											[numpy.ctypeslib.ndpointer(dtype=numpy.float64, shape=(n,), flags='F_CONTIGUOUS'), u],
										])

	fem1d_heat_steady_so.fem1d_heat_steady_.argtypes = fem1d_heat_steady_types

	fem1d_heat_steady_so.fem1d_heat_steady_( * fem1d_heat_steady_args )

	return u

"""
demo
"""

if __name__ == '__main__':

	n = 2**7
	n_lf = 2**3
	a = 0.
	b = 1.
	ua = 0.
	ub = 0.

	a_par = 0.5
	def k(x):
		py_x = x[0]
		return 1. + a_par*py_x	

	Z_par = 0.9
	def f(x):
		py_x = x[0]
		return 50. * Z_par**2

	x = numpy.linspace(a, b, n)
	x_lf = numpy.linspace(a, b, n_lf)

	u = fem1d_heat_steady(n, a, b, ua, ub, k, f, x)
	u_lf = fem1d_heat_steady(n_lf, a, b, ua, ub, k, f, x_lf)

	import matplotlib.pylab as plt

	plt.plot(x, u, '-o')
	plt.plot(x_lf, u_lf)
	plt.show()