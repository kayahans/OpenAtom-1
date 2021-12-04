
#include "matrix.h"
#include "eps_matrix.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

class DiagBridge : public CBase_DiagBridge {
  DiagBridge_SDAG_CODE
  private:
   unsigned int data_received = 0;
   unsigned total_data = 0;
   bool busy, result;
   CthThread t;
   void (DiagBridge::*dataHandler)(DiagMessage*);

  public:
    int totaldata;
    int row_size;
    int col_size;
    int total_transfers = 49;
    double* data;
    int x;
    int y;
    int eps_pe;


    DiagBridge() {};

    DiagBridge(int totaldata) : totaldata(totaldata) {
      data = new double[totaldata];
    }

    void pup(PUP::er &p) {
      p|totaldata;
      p|row_size;
      p|col_size;
      p|x;
      p|y;
      p|eps_pe;
      if (p.isUnpacking())
        data = new double[totaldata];
      PUParray(p, data, totaldata);
    }
    void prepareData(int qindex, int size);
    void receiveDataDSimple(DiagMessage* msg);
    void waitFor();
    // void receiveHeapDSimple(const DiagBridge &inData) {
      // // totaldata = inData.totaldata;
      // // data = new double[inData.totaldata];
      // int mype = CkMyPe();
      // for (int i = 0; i < 1; i++) {
      //   // data[i] = inData.data[i];
      //   if (i < 1) {
      //     CkPrintf("[DIAGBRIDGE] x %d y %d eps_pe %d diag_pe %d i %d value %.6e\n",
      //       inData.x, inData.y, inData.eps_pe, mype, i, inData.data[i]);
      //   }
      // }
      // thisProxy[thisIndex].waitForData()
    // }
    
};

class EpsMatrix : public CBase_EpsMatrix {
  EpsMatrix_SDAG_CODE
  private:
    unsigned K, L; // Number of occupied psis
    int* nfft; // number of fft grids in each direction
    unsigned qindex;

    unsigned data_received;
    unsigned proc_rows;
    unsigned proc_cols;
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
    void done(int result);
    void createCopy(CProxy_EpsMatrix other, bool todo);
    void recvCopy(std::vector<complex> new_data);
    void setI(CLA_Matrix_interface mat, bool clean);
    void receiveConvCheck(std::vector<complex> incoming);
    DiagMessage* receiveDataSimple();
    void receiveHeapSimple();
    // void copyToMPI(int qindex, int epsilon_size, CProxy_DiagBridge diag_proxy);
    static void done_cb(void *obj){
     ((EpsMatrix*) obj)->round_done();
    }
};

