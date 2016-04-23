"""
FEM and Chebyshev collocation solutions to 1D diffusion equation

John Jakeman (jdjakem@sandia.gov) and Michael Asher (michael.james.asher@gmail.com)
September 2015

"""

import numpy 
import scipy.interpolate

"""
Wrapper for [John Burkardt's 1D FEM steady state heat equation code](https://people.sc.fsu.edu/~jburkardt/f_src/fem1d_heat_steady/fem1d_heat_steady.html)
Note that conductivity k and forcing function f specified in python must treat arguments as pointers/arrays.
"""

from fem1d_heat_steady_wrapper import fem1d_heat_steady

# solve 1d diffusion using FEM
def run_fem1d(Zs, x, forcing, diffusivity):
	n = len(x) # 2**7
	a = x[0] # 0.
	b = x[-1] # 1.
	ua = 0.
	ub = 0.

	num_samples = Zs.shape[1]
	solutions = numpy.empty((order, num_samples))
	for i in range(num_samples):

		def k(x_pointer):
			return diffusivity(Zs[:,i], x_pointer[0])

		def f(x_pointer):
			return forcing(x_pointer[0])

		# x = numpy.linspace(a, b, n)

		u = fem1d_heat_steady(n, a, b, ua, ub, k, f, x)
		solutions[:, i] = u

	return solutions

'''
Chebyshev collocation for 1D diffusion equation
functions are from John Jakeman's heat/pyheat/models/spectral_collocation/diffusion_1d.py
'''

def chebyshev_derivative_matrix( order ):
    if order == 0:
        pts = numpy.array( [1], float )
        derivative_matrix = numpy.array( [0], float )
    else:
        # this is reverse order used by matlab cheb function
        pts = -numpy.cos( numpy.linspace(0.,numpy.pi,order+1 ) )
        # mja - scale from [-1,1] to [0,1]
        pts = (pts+1.)/2.
        scalars = numpy.ones( ( order+1 ), float )
        scalars[0] = 2.; scalars[order] = 2.
        scalars[1:order+1:2] *= -1
        derivative_matrix = numpy.empty( ( order+1, order+1 ), float )
        for i in xrange( order+1 ):
            row_sum = 0.
            for j in xrange( order+1 ):
                if (i==j): denominator = 1.
                else: denominator = pts[i]-pts[j]
                numerator = scalars[i] / scalars[j]
                derivative_matrix[i,j] = numerator / denominator;
                row_sum += derivative_matrix[i,j]
            derivative_matrix[i,i] -= row_sum

    # I return points and calculate derivatives using reverse order of points
    # compared to what is used by Matlab cheb function thus the 
    # derivative matrix I return will be the negative of the matlab version
    return pts, derivative_matrix


def form_collocation_matrix(derivative_matrix, diagonal ):
    scaled_matrix = numpy.empty( derivative_matrix.shape )
    for i in xrange( scaled_matrix.shape[0] ):
        scaled_matrix[i,:] = derivative_matrix[i,:] * diagonal[i]
    matrix = numpy.dot( derivative_matrix, scaled_matrix )
    return matrix


def apply_boundary_conditions(matrix, forcing ):
    matrix[0,:] = 0; matrix[-1,:] = 0
    matrix[0,0] = 1; matrix[-1,-1] = 1
    forcing[0] = 0; forcing[-1] = 0;
    return matrix, forcing


# solve 1d diffusion using Chebyshev collocation
def run_chebyshev_collocation(Zs, order, forcing, diffusivity, pts, D):
	num_samples = Zs.shape[1]
	solutions = numpy.empty((order+1, num_samples))
	for i in range(num_samples):
		diagonal = [diffusivity(Zs[:,i], x) for x in pts]
		collocation_matrix = form_collocation_matrix(-1*D, diagonal)
		forcing_vector = [forcing(x) for x in pts]
		m, f = apply_boundary_conditions(collocation_matrix, forcing_vector)
		# m, f = (collocation_matrix, forcing_vector)
		solutions[:, i] = numpy.linalg.solve(m,f)
	return solutions


'''
Barycentric Lagrange interpolation (JJ and scipy implementations)
	arguments as numpy.interp
'''


def lagrange_interpolation_scipy( x, xp, fp ):
    interpolating_polynomial = scipy.interpolate.BarycentricInterpolator(xp, fp)
    return interpolating_polynomial(x)




# from John Jakeman's heat/pyheat/models/spectral_collocation/diffusion_1d.py
# import sys
# sys.path.append("/home/mikey/Desktop/keep/heat/cppheat/build/swig/")
# sys.path.append("/home/mikey/Desktop/keep/heat/pyheat/")
# from interpolation_cpp import compute_barycentric_weights_1d, multivariate_barycentric_lagrange_interpolation

# def lagrange_interpolation_JJ( x, xp, fp ):
#     if x.ndim == 1:
#         x = x.reshape( ( 1, x.shape[0] ) )
#     num_dims = x.shape[0]
#     weights_1d = compute_barycentric_weights_1d( xp )
#     interp_vals = multivariate_barycentric_lagrange_interpolation( 
#         x,
#         [xp]*num_dims, 
#         [weights_1d]*num_dims,
#         fp.reshape( ( 1, fp.shape[0] ) ),
#         # cannot use numpy.int (c++ long) for IntVectors
#         # must use numpy.int32 (c++ int)
#         numpy.arange( num_dims, dtype=numpy.int32 ) ).squeeze()
#     return interp_vals


if __name__ == '__main__':

	'''
	compare FEM and Chebyshev Collocation
	'''

	# equation (4.10) from http://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1048&context=prism
	# also (5.5) from "Computational Aspects of Stochastic Collocation with Multifidelity Models"
	def diffusivity_function(Z, x):
		sigma = .5
		return 1. + sigma * numpy.sum([Z[k-1]*numpy.cos(2*numpy.pi*k*x)/(k**2 * numpy.pi**2) for k in range(1, len(Z)+1)])

	def forcing_function(x):
		return 1.

	num_dims = 10
	num_samples = 10
	sample_to_plot = 0
	Zs = numpy.random.RandomState().uniform(-1.,1.,(num_dims, num_samples))

	order = 2**7
	x_cheb, D = chebyshev_derivative_matrix(order-1)

	u_cheb_pre = run_chebyshev_collocation(Zs, order-1, forcing_function, diffusivity_function, x_cheb, D)
	u_cheb_pre *= -1.
	u_cheb_pre = u_cheb_pre[:,sample_to_plot]

	# either solve FEM on grid, then interpolate Cheby to that grid
	#	or just solve at Cheby's nodes/locations/abscissae 
	interpolate = True

	if interpolate:
		x_fem = numpy.linspace(0., 1., order)
		u_fem = run_fem1d(Zs, x_fem, forcing_function, diffusivity_function)
		u_fem = u_fem[:,sample_to_plot]
		u_cheb = lagrange_interpolation_scipy(x_fem, x_cheb, u_cheb_pre)

		# u_cheb_jj = lagrange_interpolation_JJ(x_fem, x_cheb, u_cheb_pre)
		# assert numpy.allclose(u_cheb_jj, u_cheb)

	else:	
		x_fem = x_cheb
		u_fem = run_fem1d(Zs, x_fem, forcing_function, diffusivity_function)
		u_fem = u_fem[:,sample_to_plot]
		u_cheb = u_cheb_pre

	assert numpy.allclose(u_cheb, u_fem)

	import matplotlib.pylab as plt

	plt.plot(x_cheb, u_cheb,'x', label="Cheby")
	plt.plot(x_cheb, u_fem, label="FEM")
	plt.legend()
	plt.show()

	plt.plot(u_fem - u_cheb)
	plt.title('Difference between FEM and Chebyshev')
	plt.show()


	'''
	check Chebyshev Collocation interpolation (for comparing multifidelity)
	'''

	lf_n = 2**6-1
	lf_x, lf_D = chebyshev_derivative_matrix(lf_n)
	lf_u = run_chebyshev_collocation(Zs, lf_n, forcing_function, diffusivity_function, lf_x, lf_D)
	hf_n = 2**7-1
	hf_x, hf_D = chebyshev_derivative_matrix(hf_n)
	hf_u = run_chebyshev_collocation(Zs, hf_n, forcing_function, diffusivity_function, hf_x, hf_D)

	interpolated_lf_u = lagrange_interpolation_scipy( hf_x, lf_x, lf_u )

	# assert numpy.allclose( hf_u,  interpolated_lf_u)
	print "L2 error between HF and LF", numpy.linalg.norm(hf_u - interpolated_lf_u, ord=2)

	plt.plot(hf_x, hf_u[:, sample_to_plot], '.', label="hf")
	plt.plot(hf_x, interpolated_lf_u[:, sample_to_plot], label="lf")
	plt.legend()
	plt.show()