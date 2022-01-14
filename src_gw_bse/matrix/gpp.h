
#include "matrix.h"
#include "gpp.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

class DiagBridge : public CBase_DiagBridge {
  DiagBridge_SDAG_CODE
  
  public:
    unsigned int row_size, col_size;
    unsigned int totaldata; // Total data in the diagData pointer
    unsigned int numBlocks; // Number of blocks MB X MB block grid
    unsigned int proc_rows, proc_cols; // processor grid
    double* data;
    DiagBridge() {};

    DiagBridge(int totaldata) : totaldata(totaldata) {
      data = new double[totaldata];
    }
    void prepareData(int qindex, int size, int num_qpts);
};


class Gpp : public CBase_Gpp {
  Gpp_SDAG_CODE
  private:
    // Structural basic input 
    // used in setQindex
    int qindex;
    double* qcryst;  // q-point in crystalline coord.
    double qcart[3]; // q-point in cartesian coord
    double vol;   // volume of the cell
    double alat;   // Alat from QE
    int nkpt;     // number of k-points
    double* b1; // reciprocal lattice parameters b1-b3 in cart. coord. [3] array
    double* b2;
    double* b3;
    int* nfft;
    int ndata;

    // Strictly GPP related input below
    // used in readInputFile
    bool qespresso;      // if rho and state data are from Quantum Espresso(true/1) or not(false/0)
    char rhoFile[1000];  // file name that reads density data
    int num_q;              // number of q vectors
    int num_w;       // number of points for w that we want to calculate.
    double *w;      // values to be evaluated. read in
    // used in readRho
    int nr[3];           // number of (real-space) data points in each direction for rhofile
    int *ga, *gb, *gc;   // g index for density
    std::complex<double>* rhoData;
    double Wpl2;         // Plasmon frequency squared
    // used in calculate_vc
    double* vcoulb;
    int ng; // check if this is filtered?

    double* omsq;
    double* eigval;

    CLA_Matrix_interface matrix;
    unsigned data_received;
    double total_time;
    unsigned int blockSize, numBlocks, block;
    
  public:
    bool gpp_is_on = false;      // if gpp is on(true/1) or not(false/0)
    Gpp();
    Gpp(MatrixConfig config);
    void readRho();
    void readInputFile();
    void setQIndex(int qindex);
    void find_rho_gdiff(int gdiff[3], int &gindex, bool &gdiffTrue);
    void calc_Omsq();
    void calculate_vc();
    void recvCopy(std::vector<complex> new_data);
    DiagMessage* receiveDataSimple(DiagMessage* msg);
    DiagMessage* sendDataSimple(DiagMessage* msg);
    
};
