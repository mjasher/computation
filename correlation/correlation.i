%module correlation

%{
    #define SWIG_FILE_WITH_INIT
    #include "correlation.hpp"
%}

%include "numpy.i"

%init %{
    import_array();
%}

// TODO very curious what the logic behind this is
// see http://docs.scipy.org/doc/numpy/reference/swig.interface-file.html
%apply (double* IN_ARRAY1, int DIM1) {(double* x1, int l1), (double* x2, int l2), (double *in, int nin)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *out, int nout)}
%include "correlation.hpp"