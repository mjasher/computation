computation
===========

Computational methods I've implemented or collected.

[cache](cache)
-----------
A succinct and fairly flexible way to cache a given python function.
Easy as 
```
@cache.money(key_func, cache_directory)
def slow_function(a, b, c):
    return y
```

[ctypes_boiler](ctypes_boiler)
-----------
Python performance can be dramatically improved by writing slow functions in C. [Ctypes](https://docs.python.org/2/library/ctypes.html) allows C functions to be called from python. Here are some simple examples, including some using numpy arrays.


[fem1d_heat_steady](fem1d_heat_steady)
-----------
Ctypes can also be used to wrap Fortran in python! Excellent if you want access to troves of well written numerical code such as [John Burkardt's collection](http://people.sc.fsu.edu/~jburkardt/). 
[fem1d_heat_steady](fem1d_heat_steady/fem1d_heat_steady_wrapper.py) for example, wraps [John Burkardt's 1D FEM steady state heat equation code](https://people.sc.fsu.edu/~jburkardt/f_src/fem1d_heat_steady/fem1d_heat_steady.html).


[correlation](correlation)
-----------
To call c++ from python, SWIG can be used. Here is a simple example for computing a correlation matrix used in applications such as Karhunen-Loève expansions.

