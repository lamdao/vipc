//--------------------------------------------------------------------------
// Eigen decomposition code for symmetric 3x3 matrices,
// copied from the public domain Java Matrix library JAMA.
//--------------------------------------------------------------------------
#ifndef __EIGEN3x3_H
#define __EIGEN3x3_H
//--------------------------------------------------------------------------
// Symmetric matrix Mx => eigenvectors in columns of V,
//                        corresponding eigenvalues in d
//--------------------------------------------------------------------------
void eigen_decomposition(double *Mx, double *V, double *d);
//--------------------------------------------------------------------------
#endif
