#include "matrix.h"
#include "gpp.decl.h"


#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

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
    bool qespresso = true;      // if rho and state data are from Quantum Espresso(true/1) or not(false/0)
    int num_q = 8;              // number of q vectors
    int num_w = 1;       // number of points for w that we want to calculate.
    double *w;      // values to be evaluated. read in
    // int nr[3];           // number of (real-space) data points in each direction for rhofile
    
    // used in calculate_vc
    double* vcoulb;
    unsigned ng; // check if this is filtered?
    
    double factor;
    double* omsq;
    double* eigval;
    complex* rhoData;
    CLA_Matrix_interface matrix;
    unsigned data_received;
    double total_time;
    unsigned int blockSize, numBlocks, block;
    // complex* rhoDataG;
    
  public:
    bool gpp_is_on = false;      // if gpp is on(true/1) or not(false/0)
    Gpp();
    Gpp(MatrixConfig config);
    void RtoG(int qindex);
    void calc_omsq();
    void debug();
    void print(int qindex, int fnum);
    void setQIndex(int qindex);
    // void find_rho_gdiff(int gdiff[3], int &gindex, bool &gdiffTrue);
    void calc_M0();
    PsiMessage* send_data(PsiMessage* msg);
    void print_col(int num);
    void calculate_vc();
    void recvCopy(std::vector<complex> new_data);
    void recv_eig(std::vector<double> new_data);
    void sendToCacheV(int size);
    void sendToCacheE(int size);
    void sendToCacheO(int size);
};
