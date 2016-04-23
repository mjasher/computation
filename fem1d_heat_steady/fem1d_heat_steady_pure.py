import numpy

def fem1d_heat_steady( n, a, b, ua, ub, k, f, x ):

  #*****************************************************************************80
  #
  # FEM1D_HEAT_STEADY solves the steady 1D heat equation with finite elements.
  #
  #  Discussion:
  #
  #    The program uses the finite element method, with piecewise linear basis
  #    functions to solve the steady state heat equation in one dimension.
  #
  #    The problem is defined on the region A <= x <= B.
  #
  #    The following differential equation is imposed between A and B:
  #
  #      - d/dx k(x) du/dx = f(x)
  #
  #    where k(x) and f(x) are given functions.
  #
  #    At the boundaries, the following conditions are applied:
  #
  #      u(A) = UA
  #      u(B) = UB
  #
  #    A set of N equally spaced nodes is defined on this
  #    interval, with A = X(1) < X(2) < ... < X(N) = B.
  #
  #    At each node I, we associate a piecewise linear basis function V(I,X),
  #    which is 0 at all nodes except node I.  This implies that V(I,X) is
  #    everywhere 0 except that
  #
  #    for X(I-1) <= X <= X(I):
  #
  #      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
  #
  #    for X(I) <= X <= X(I+1):
  #
  #      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
  #
  #    We now assume that the solution U(X) can be written as a linear
  #    sum of these basis functions:
  #
  #      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
  #
  #    where U(X) on the left is the function of X, but on the right,
  #    is meant to indicate the coefficients of the basis functions.
  #
  #    To determine the coefficient U(J), we multiply the original
  #    differential equation by the basis function V(J,X), and use
  #    integration by parts, to arrive at the I-th finite element equation:
  #
  #        Integral K(X) * U'(X) * V'(I,X) dx = Integral F(X) * V(I,X) dx
  #
  #    We note that the functions U(X) and U'(X) can be replaced by
  #    the finite element form involving the linear sum of basis functions,
  #    but we also note that the resulting integrand will only be nonzero
  #    for terms where J = I - 1, I, or I + 1.
  #
  #    By writing this equation for basis functions I = 2 through N - 1,
  #    and using the boundary conditions, we have N linear equations
  #    for the N unknown coefficients U(1) through U(N), which can
  #    be easily solved.
  #
  #  Licensing:
  #
  #    This code is distributed under the GNU LGPL license.
  #
  #  Modified:
  #
  #    07 April 2011
  #
  #  Author:
  #
  #    John Burkardt
  #
  #  Parameters:
  #
  #    Input, integer ( kind = 4 ) N, the number of nodes.
  #
  #    Input, real ( kind = 8 ) A, B, the left and right endpoints.
  #
  #    Input, real ( kind = 8 ) UA, UB, the prescribed value of U at A and B.
  #
  #    Input, external K, a function which evaluates k(x);
  #
  #    Input, external F, a function which evaluates f(x);
  #
  #    Input, real ( kind = 8 ) X(N), the mesh points.
  #
  #    Output, real ( kind = 8 ) U(N), the finite element coefficients, which 
  #    are also the value of the computed solution at the mesh points.
  #

  #
  #  Define a quadrature rule on the interval [-1,+1].
  #
  quad_num = 2
  abscissa = numpy.empty((quad_num), dtype=numpy.float64)    
  abscissa[0] = -0.577350269189625764509148780502E+00
  abscissa[1] = +0.577350269189625764509148780502E+00
  weight = numpy.empty((quad_num), dtype=numpy.float64)    
  weight[0] = 1.0E+00
  weight[1] = 1.0E+00
  #
  #  Zero out the matrix and right hand side.
  #
  amat = numpy.zeros((n,n), dtype=numpy.float64)
  bvec = numpy.zeros((n), dtype=numpy.float64)
  #
  #  Equation 1 is the left boundary condition, U(A) = UA;
  #
  amat[0,0] = 1.0E+00
  bvec[0] = ua
  #
  #  Equation I involves the basis function at node I.
  #  This basis function is nonzero from X(I-1) to X(I+1).
  #  Equation I looks like this:
  #
  #    Integral K(X) U'(X) V'(I,X) 
  #           + C(X) * U(X) V(I,X) dx 
  #  = Integral F(X) V(I,X) dx
  #
  #  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
  #  (U(X) means the function; U(J) is the coefficient of V(J,X) ).
  #
  #  The only V functions that are nonzero when V(I,X) is nonzero are
  #  V(I-1,X) and V(I+1,X). 
  #
  #  Let's use the shorthand 
  #
  #    VL(X) = V(I-1,X)
  #    VM(X) = V(I,X)
  #    VR(X) = V(I+1,X)
  #
  #  So our equation becomes
  #
  #    Integral K(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X) dx
  #  = Integral F(X) VM(X) dx.
  #
  #  
  #
  #  This is actually a set of N-2 linear equations for the N coefficients U.
  #
  #  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
  #  and so on.
  #
  for i in range(1, n-1):

    #
    #  Get the left, right and middle coordinates.
    #
    xl = x[i-1]
    xm = x[i]
    xr = x[i+1]
    #
    #  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
    #
    al = 0.0E+00
    am = 0.0E+00
    ar = 0.0E+00
    bm = 0.0E+00
    #
    #  We approximate the integrals by using a weighted sum of
    #  the integrand values at quadrature points.
    #
    
    for q in range(quad_num):

      #
      #  Integrate over the LEFT interval, between XL and XM, where:
      #
      #  VL(X) = ( XM - X       ) / ( XM - XL )
      #  VM(X) = (      X  - XL ) / ( XM - XL )
      #  VR(X) = 0
      #
      #  VL'(X) =             - 1 / ( XM - XL )
      #  VM'(X) =             + 1 / ( XM - XL ) 
      #  VR'(X) = 0
      #
      xq = ( ( 1.0E+00 - abscissa[q] ) * xl \
           + ( 1.0E+00 + abscissa[q] ) * xm ) \
           /   2.0E+00

      wq = weight[q] * ( xm - xl ) / 2.0E+00

      vl =  ( xm - xq ) / ( xm - xl )
      vlp =  - 1.0E+00  / ( xm - xl )

      vm =  ( xq - xl ) / ( xm - xl )
      vmp =  + 1.0E+00  / ( xm - xl )

      vr =  0.0E+00
      vrp = 0.0E+00

      kxq = k( xq )
      fxq = f( xq )

      al = al + wq * ( kxq * vlp * vmp )
      am = am + wq * ( kxq * vmp * vmp )
      ar = ar + wq * ( kxq * vrp * vmp )
      bm = bm + wq * ( fxq * vm )
      #
      #  Integrate over the RIGHT interval, between XM and XR, where:
      #
      #  VL(X) = 0
      #  VM(X) = ( XR - X       ) / ( XR - XM )
      #  VR(X) = (      X  - XM ) / ( XR - XM )
      #
      #  VL'(X) = 0
      #  VM'(X) =             - 1 / ( XR - XM )
      #  VR'(X) =             + 1 / ( XR - XM ) 
      #
      xq = ( ( 1.0E+00 - abscissa[q] ) * xm   \
           + ( 1.0E+00 + abscissa[q] ) * xr ) \
           /   2.0E+00

      wq = weight[q] * ( xr - xm ) / 2.0E+00

      vl = 0.0E+00
      vlp = 0.0E+00

      vm = ( xr - xq ) / ( xr - xm )
      vmp = - 1.0E+00  / ( xr - xm )

      vr = ( xq - xm ) / ( xr - xm )
      vrp =  1.0E+00   / ( xr - xm )

      kxq = k ( xq )
      kxq = -50.0 * 0.2**2
      fxq = f ( xq )
      fxq = 1.0 + 0.5*xq

      al = al + wq * ( kxq * vlp * vmp )
      am = am + wq * ( kxq * vmp * vmp )
      ar = ar + wq * ( kxq * vrp * vmp )
      bm = bm + wq * ( fxq * vm )

    amat[i,i-1] = al
    amat[i,i]   = am
    amat[i,i+1] = ar
    bvec[i]     = bm

#
#  Equation N is the right boundary condition, U(B) = UB;
#
  amat[n-1,n-1] = 1.0E+00
  bvec[n-1] = ub
#
#  Solve the linear system.
#
  # call r8mat_print ( n, n, amat, '  Matrix A:' )

  # call r8mat_solve2 ( n, amat, bvec, u, ierror )
  print amat
  print bvec
  u = numpy.linalg.solve(amat, bvec)

  return u


if __name__ == '__main__':
    
  n = 128
  a = 0.
  b = 1.
  ua = 0.
  ub = 0.

  a_par = 0.5
  def k(x):
    return 1. + a_par*x  

  Z_par = 0.2
  def f(x):
    return -50. * Z_par**2

  x = numpy.linspace(a, b, n)
  print "X", len(x)

  u = fem1d_heat_steady(n, a, b, ua, ub, k, f, x)

  print "ANSWER"
  print u

  import matplotlib.pylab as plt

  plt.plot(u, '-o')
  # plt.plot(x)
  plt.show()