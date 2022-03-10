#ifndef PMATRIX_H
#define PMATRIX_H

#include "matrix.h"
#include "pmatrix.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"
#include "../utils/sigma_util.h"
#include "../utils/windowing.h"
class FFTController;

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix(MatrixConfig config);
    PMatrix(CkMigrateMessage* msg) {}

    void fftRows(int);
    void generateEpsilon(CProxy_EpsMatrix proxy, std::vector<int> accept);
    void applyFs();
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);
    void computeP();
    void sigma();
    void sigma_init(std::vector<int> accept);

    void reportPTime();
    void generateEpsilon(std::vector<double> vcoulb, std::vector<int> accept, int inew, int jnew, int size, int max_inew, int max_jnew);   
    void registerTileSections();
 
  private:
    unsigned L; // Number of occupied psis
    unsigned trans_count, num_chares; // SDAG variables
    unsigned completed_chunks;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    int psi_ndata_local;
    int tile_size;
    int region_ridx;
    int region_cidx;
    FFTController* fft_controller;

    unsigned start_index, end_index;
    unsigned chare_chunk;
    unsigned send_count, recv_count;
    unsigned iteration;
    unsigned max_iterations;

    int arrival_counter;
    int receive_counter;

    complex *P_m; // N3 P local tile
    complex *F_m; // N3 Sigma local psi tile
    complex *B_m; // N3 Sigma local gpp tile
    complex *sigma_m; // N3 Sigma tile
    std::vector<int> accept;

    complex total[144];
    void kqIndex(unsigned, unsigned&, int*);
    complex* umklapp_factor;
    void getUmklappFactor(complex*, int[3]);
    void cubic_sigma_per_window(const int is,
                                    const int ik,
                                    const int iq,
                                    const int ikq,
                                    const WINPAIR& winpair,
                                    const std::vector<double>& state_e,
                                    const std::vector<int>& state_idx,
                                    const std::vector<double>& pp_e,
                                    const std::vector<int>& pp_idx,
                                    const int uklpp[3]);
    void PerformPPEnergySumThisNode(const int is,
                                        const int iq,
                                        const std::vector<double>& pp_wp,
                                        const std::vector<int>& pp_idx,
                                        const WINPAIR& winpair,
                                        const int inode,
                                        const bool IsHGL,
                                        const bool IsFirstHGL);
    void PerformStateSumThisNode(const int& is,
                                    const int& ik,
                                    const int& ikq,
                                    const std::vector<double>& state_e,
                                    const std::vector<int>& state_idx,
                                    const WINPAIR& winpair,
                                    const int& inode,
                                    const bool& IsHGL,
                                    const bool& IsFirstHGL,
                                    const int uklpp[3]);

    void compute_fr(complex* fr, complex* psikq, const int uklpp[3]);
    void sigma_cubic_main(std::complex<double>* sigma,
                              const int& is,
                              const int& ik,
                              const double& w,
                              const std::vector<std::pair<int, int>>& n12,
                              const bool& bIsOccupied);
    void cubicSigma(const int& is,
                    const int& ik,
                    const SIGMAINDICES& iwn12);
    double total_time;
};

extern /* readonly */ CProxy_PMatrix pmatrix2D_proxy;
extern /* readonly */ CProxy_PMatrix pmatrix1D_proxy;
#endif
