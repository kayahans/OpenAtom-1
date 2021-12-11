#ifndef _diagonalizer_h_
#define _diagonalizer_h_
// #include "diagonalizer.decl.h"

struct diagData_t {
  int qindex;

  int num_handoffs;  // Number of diagonalizations to be performed = num of q points

  // Dimensions of data stored here
  int inputsize = 0;
  int row_size = 0;
  int col_size = 0;
  int nprow = 0; // number of processors in row of processor grid
  int npcol = 0;
  int n = 0; // size of the original epsilon matrix
  int nb = 0; // block size
  
  double* input;  // Input matrix (row_size x col_size = inputsize)
  double* eig_v;  // Eigenvectors (n x n)
  double* eig_e;  // Eigenvalues (n x 1)

};

void restartCharm();
#endif  // _diagonalizer_h_
