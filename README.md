computation
===========

I enjoy maths, python, and fast code.
Here are some methods like minded people may find useful.

[b64_array_js_py](b64_array_js_py)
-----------
Demo of how to transfer a float32 array javascript <--> python as a base64 encoded byte string.


[cache](cache)
-----------
A succinct and fairly flexible way to cache a given python function.
Easy as 
```
@cache.money(key_func, cache_directory)
def slow_function(a, b, c):
    return y
```

Here's a [demo](cache/cache.py#L135).

[ctypes_boiler](ctypes_boiler)
-----------
Python performance can be dramatically improved by writing slow functions in C. [Ctypes](https://docs.python.org/2/library/ctypes.html) allows C functions to be called from python. 
Here's a [demo](ctypes_boiler/local.py) which compiles and runs some trivial examples, including some using numpy arrays.


[fem1d_heat_steady](fem1d_heat_steady)
-----------
Ctypes can also be used to wrap Fortran in python! Excellent if you want access to troves of well written numerical code such as [John Burkardt's collection](http://people.sc.fsu.edu/~jburkardt/). 
[This example](fem1d_heat_steady/fem1d_heat_steady_wrapper.py) wraps [John Burkardt's 1D FEM steady state heat equation code](https://people.sc.fsu.edu/~jburkardt/f_src/fem1d_heat_steady/fem1d_heat_steady.html).

This also contains [a comparison of Chebyshev Collocation and Finite Element solutions of the 1d heat equation](fem1d_heat_steady/chebyshev_collocation_fem1d_heat_steady.py).


[correlation](correlation)
-----------
To call c++ from python, SWIG can be used. Here is a simple example for computing a correlation matrix used in applications such as Karhunen-Loève expansions.

