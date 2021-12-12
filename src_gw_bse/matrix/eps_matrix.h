#include "matrix.h"
#include "eps_matrix.decl.h"

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

class EpsMatrix : public CBase_EpsMatrix {
  EpsMatrix_SDAG_CODE
  private:
    unsigned K, L; // Number of occupied psis
    int* nfft; // number of fft grids in each direction
    unsigned qindex;

    unsigned data_received;
    double total_time;
    CLA_Matrix_interface matrix;
    unsigned int blockSize, numBlocks, block;

  public:
    EpsMatrix();
    EpsMatrix(MatrixConfig config);
    EpsMatrix(CkMigrateMessage* msg) {}

    void checkReady();
    void createTranspose(CProxy_EpsMatrix other, bool todo);
    void createConjugate(CProxy_EpsMatrix other);
    void receiveTranspose(std::vector<complex> incoming);
    void receiveConjugate(std::vector<complex> incoming);
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);
    void receiveFs(Phase3Message* msg);
    void multiply(double alpha, double beta);
    void round_done(void);
    void findAlpha(void);
    void screenedExchange();
    void bareExchange();
    void coh();
    void scalar_multiply(double alpha);
    void convergence_check(CProxy_EpsMatrix cmp_proxy);
    void add_compl_two();
    void multiply_coulb();
    void createCopy(CProxy_EpsMatrix other, bool todo);
    void recvCopy(std::vector<complex> new_data);
    void setI(CLA_Matrix_interface mat, bool clean);
    void receiveConvCheck(std::vector<complex> incoming);
    DiagMessage* receiveDataSimple(DiagMessage* msg);
    static void done_cb(void *obj){
     ((EpsMatrix*) obj)->round_done();
    }
};

