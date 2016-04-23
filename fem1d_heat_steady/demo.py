import numpy
import scipy.stats
from fem1d_heat_steady_wrapper import fem1d_heat_steady

a = 0.
b = 1.
ua = 0.
ub = 0.

n_hf = 2**9
x_hf = numpy.linspace(a, b, n_hf)

# number of random paramters/inputs
num_dims = 11
# number of quantities of interest/outputs
num_QOI = n_hf
# num_QOI = n_hf+1
# number of initial candidates/snapshots for low-fidelity model
num_lf_candidates = 1e2
# number of interpolations nodes/high-fidelity runs
# num_hf_runs = 2
# number of low-fidelity, multi-fidelity, and high-fidelity runs to do for testing
num_test_samples = 1e2


def post_process(uf):
	nf = len(uf)
	xf = numpy.linspace(a, b, nf)
	return numpy.interp(x_hf, xf, uf)


# should return (num_dims, num_samples)
def sample_inputs(num_samples):
	# KLE coefficients
    input_samples = numpy.random.RandomState().uniform(-1,1.,(num_dims, num_samples))
    # forcing
    input_samples[10,:] = (input_samples[10,:]+1.)/2.
    input_samples[10,:] = scipy.stats.distributions.norm(loc=0., scale=1.).ppf(input_samples[10,:])
    return input_samples


def build(n):

	def run(z):
		Z = z[:-1]
		f_par = z[-1]


		def k(x_pointer):
			x = x_pointer[0]
			# equation (4.10) from http://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1048&context=prism
			sigma = 1.
			# Z = numpy.random.uniform(low=-1., high=1., size=10)
			return 1. + sigma * numpy.sum([Z[k-1]*numpy.cos(2*numpy.pi*k*x)/(k**2 * numpy.pi**2) for k in range(1, len(Z)+1)])

		def f(x):
			py_x = x[0]
			return 50. * f_par**2

		x = numpy.linspace(a, b, n)

		u = fem1d_heat_steady(n, a, b, ua, ub, k, f, x)

		return post_process(u)

	return run

if __name__ == '__main__':
	
	import matplotlib.pylab as plt

	n = 2**4
	low_fidelity = build(n)
	high_fidelity = build(n_hf)

	zs = sample_inputs(5)
	for i in range(zs.shape[1]):
		u = low_fidelity(zs[:,i])
		plt.plot(x_hf,u, '-', label='lf')
		u_hf = high_fidelity(zs[:,i])
		plt.plot(x_hf,u_hf, '-', label='hf')
		plt.legend()
		plt.show()
