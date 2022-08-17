// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef johaTools_H
#define johaTools_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
mat M(mat X,mat Y);
mat S(mat X, mat Y);
mat M_perp(mat M);
mat M_bar(mat M);
// This is the end of the header guard
#endif
