#ifndef _diagonalizer_h_
#define _diagonalizer_h_
// #include "diagonalizer.decl.h"

struct diagData_t {
  int dim;  // dim x dim matrix
  int qindex;

  int num_handoff = 8;  // Number of diagonalizations to be performed

  double* input;  // Input matrix (dim x dim)
  double* eig_v;  // Eigenvectors (dim x dim)
  double* eig_e;  // Eigenvalues (dim x 1)


  // int numStateddsOA;
  // int grainSizeOrtho;
  // int numEOrthosPerDim = 1;
  // int numOrthosPerDim = 1;
  int inputsize = 0;
  int row_size = 0;
  int col_size = 0;
  int nprow = 2;
  int npcol = 2;
  int n = 8;
  int nb = 16;
};


#endif  // _diagonalizer_h_
