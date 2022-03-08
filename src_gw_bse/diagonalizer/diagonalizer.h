#include <complex>
#ifndef _diagonalizer_h_
#define _diagonalizer_h_
struct matel {
  int x;
  int y;
  std::complex<double> val;
};

struct diagData_t {
  int qindex;

  int num_handoffs;  // Number of diagonalizations to be performed = num of q points

  // Dimensions of data stored here
  int inputsize = 0;
  int row_size = 0;
  int col_size = 0;
  int nprow = 0; // number of processors in row of processor grid
  int npcol = 0;
  int pe_row;
  int pe_col;
  int n = 0; // size of the original epsilon matrix
  int nb = 0; // block size
  
  std::complex<double>* input;  // Input matrix (row_size x col_size = inputsize)
  std::complex<double>* eig_v;  // Input matrix (row_size x col_size = inputsize)
  matel* eig_v2;  // Eigenvectors (n x n)
  double* eig_e;  // Eigenvalues (n x 1)

};

void restartCharm();
#endif  // _diagonalizer_h_
