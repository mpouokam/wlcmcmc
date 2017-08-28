// DON'T FORGET TO CHANGE THIS LINE HERE FOR THE NAME OF THE MODULE!
%module wlcmcmclib

%{
    #define SWIG_FILE_WITH_INIT
    // Here make sure you include the header file as in the next line:
    #include "wlcmcmclib.h"
%}

// This file is what SWIG uses to reinterpret numpy arrays as C data types for C functions
%include "numpy.i"

// This is necessary for numpy.i thing to work.
%init %{
    import_array();
%}

%pythoncode %{
%}

// Here list functions that will be made available to Python.
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* coords, int numberofvertices, int ibe3)};
%include "wlcmcmclib.h"
//TODO: write your own typemap to take care of these guys below and get rid of wrapper functions in the cpp filenum
