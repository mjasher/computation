import numpy

from chebyshev_collocation_fem1d_heat_steady import run_chebyshev_collocation

num_dims = 13
num_samples = 1
input_samples = numpy.random.RandomState().uniform(-1.,1.,(num_dims, num_samples))

Z = input_samples[:, 0]

order = 2**7
x_cheb, u_cheb = run_chebyshev_collocation(Z, order-1)


order = 2**4
x_cheb_low, u_cheb_low = run_chebyshev_collocation(Z, order-1)

import matplotlib.pylab as plt

plt.plot(x_cheb, u_cheb, '.')
plt.plot(x_cheb_low, u_cheb_low, 'x')

print len(x_cheb)
print len(x_cheb_low)

plt.show()
