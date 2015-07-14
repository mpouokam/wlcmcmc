// DON'T FORGET TO CHANGE THIS LINE HERE FOR THE NAME OF THE MODULE!
%module wlcmclib

%{
    #define SWIG_FILE_WITH_INIT
    // Here make sure you include the header file as in the next line:
    // #include "blah.h"
    // Or do as I did and just put the function prototype here directly as I do here:
    // This is the part that is used for compiling
    #include "wlcmclib.h"
    // TODO: nothing calls these. perhaps we can delete these lines above

%}

// This file is what SWIG uses to reinterpret NumPY arrays as C data types for C functions
%include "numpy.i"

// This I think is necessary for the whole numpy.i Thing to work.
%init %{
    import_array();
%}

// This is the part that converts the call from Python with an array as an argument as
// C arguments, courtesy of numpy.i



%pythoncode %{
%}

// Here list functions that will be made available to Python. You could do
// #include "blah"
// if you wanted... but I won't.

//TODO: write your own typemap to take care of these guys below and get rid of wrapper functions in the cpp filenum
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* coords, int numberofvertices, int ibe3)};
%include "wlcmclib.h"
// void docrankshaft(double* coords, int numberofvertices, int ibe3, int vert1, int numtonextvert, double theta);
// void dorandomcrankshaft(double* coords, int numberofvertices, int ibe3, double thetamax);
// double getwlcenergy(double* coords, int numberofvertices, int ibe3, double bendingrigidityconstant);
// int domontecarlosteps(double* coords, int numberofvertices, int ibe3, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength);
// double distancebetweenedges(double* coords, int numberofvertices, int ibe3, int _a1, int _a2, int _b1, int _b2);
// bool iscollision(double* coords, int numberofvertices, int ibe3, double diameter, double edgelength);
// void setrandomseedtoclocktime(int anumber);
// void set_sRand_seed_to_clocktime();

// %rename (foo) foo_safe;
// %inline %{
//   void foo_safe() {
//     std::cout << "Hello world" << std::endl;
//     foo(); // Calls the foo() from test.h, as you'd hope
//   }
// %}
