#include "standard_include.h"
#include "allclass_gwbse.h"
#include "gpp.h"
#include "eps_matrix.h"
#include "messages.h"
#include "pmatrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"
#include "../diagonalizer/diagonalizer.h"

#include <cstring> // for memcpy
using std::memcpy;
#include <cmath>

extern diagData_t* diagData;

#define eps_rows 20
#define eps_cols 20
#define IDX_eps(r,c) ((r)*eps_cols + (c))

void example_dgemm(int M, int N, int K, double alpha,
                   complex *A, complex *B, complex *C) {
  /* multiply */
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      complex sum = 0.0;
      for (int k = 0; k < K; ++k) {
        sum += A[i*K + k] * B[k*N + j];
      }
      C[N*i + j] = C[N*i + j] + alpha*sum;
    }
  }
}

EpsMatrix::EpsMatrix() {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = 0; // The controller will set this
  CkPrintf("\nWarning!! The controller should set this q-point!\n");

  data_received = 0;
  total_time = 0.0;
}

EpsMatrix::EpsMatrix(MatrixConfig config) : CBase_EpsMatrix(config) {
  GWBSE* gwbse = GWBSE::get();
  blockSize = config.tile_rows;
  numBlocks = config.chareRows();

  // Set some constants
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = config.qindex; // The controller sets this

  total_time = 0.0;
  data_received = 0;
}

void EpsMatrix::setI(CLA_Matrix_interface mat, bool clean){
  matrix = mat;
  if (clean) {
    delete [] data;
    initialize();
  }
}

void EpsMatrix::receiveFs(Phase3Message* msg) {
  int n = 0;
  // TODO: memcpy
  for(int i=msg->start_i;i<=msg->end_i;i++)
    for(int j=msg->start_j;j<=msg->end_j;j++)
        data[IDX_eps(i,j)] = msg->data[n++];

  data_received+=n;
  if(data_received == eps_cols*eps_rows) {
    CkCallback cb(CkReductionTarget(Controller, epsilon_created), controller_proxy);
    contribute(cb);
  }
}

void EpsMatrix::multiply(double alpha, double beta) {
  matrix.multiply(alpha, beta, data, EpsMatrix::done_cb, (void*) this,
       thisIndex.x, thisIndex.y);
}

void EpsMatrix::add_compl_two() {
  int i = 0;
  complex compl_two(2.0, 0.0);
  if(thisIndex.x==thisIndex.y)
  for(int i=0;i<eps_rows;i++)
    data[IDX_eps(i,i)] += compl_two;

  CkCallback cb(CkReductionTarget(Controller, complement_multiplied), controller_proxy);
  contribute(cb);
}

void inline EpsMatrix::round_done(void) {
  CmiMemoryCheck();
  CkCallback cb(CkReductionTarget(Controller, m_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::scalar_multiply(double alpha) {
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = alpha*data[IDX_eps(i,j)]; 

  CkCallback cb(CkReductionTarget(Controller, scalar_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::screenedExchange() {

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex total_contribution(0.0,0.0);
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;

  for (int k = 0; k < K; k++) {
    for (int i = 0; i < f_cache->getNSize(); i++) {
        complex contribution(0.0,0.0);
  //      for (int j = 0; j < f_cache->getNSize(); j++) { //Performs only <n|Sigma|n> as does the fortran code
  // Uncommenting above loop will perform <n|Sigma|n’>
        for (int l = 0; l < L; l++) {
          complex* fi = f_cache->getFVec(k, i, l, thisIndex.x, eps_rows);
          complex* fj = f_cache->getFVec(k, i, l, thisIndex.y, eps_cols);
          for (int r = 0; r < config.tile_rows; r++) {
            for (int c = 0; c < config.tile_cols; c++) {
              contribution += fi[r]*fj[c].conj()*data[IDX_eps(r,c)];
            }
          }
        }
        contrib_data[ik] = contribution * -1.0;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution * -1.0;
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::screenedExchangeComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
}

void EpsMatrix::screenedExchangeGPP() {
  // TODO (kayahans) will be modified, for now just copy/paste from static
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex total_contribution(0.0,0.0);
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;

  for (int k = 0; k < K; k++) {
    for (int i = 0; i < f_cache->getNSize(); i++) {
        complex contribution(0.0,0.0);
  //      for (int j = 0; j < f_cache->getNSize(); j++) { //Performs only <n|Sigma|n> as does the fortran code
  // Uncommenting above loop will perform <n|Sigma|n’>
        for (int l = 0; l < L; l++) {
          complex* fi = f_cache->getFVec(k, i, l, thisIndex.x, eps_rows);
          complex* fj = f_cache->getFVec(k, i, l, thisIndex.y, eps_cols);
          for (int r = 0; r < config.tile_rows; r++) {
            for (int c = 0; c < config.tile_cols; c++) {
              contribution += fi[r]*fj[c].conj()*data[IDX_eps(r,c)];
            }
          }
        }
        contrib_data[ik] = contribution * -1.0;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution * -1.0;
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::screenedExchangeComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
}

void EpsMatrix::cohGPP(){

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

  complex* states = psi_cache->getStates();
  int psi_size = nfft[0]*nfft[1]*nfft[2];
  std::vector<int> accept_v = f_cache->getAcceptVector();

  int ga[psi_size];
  int gb[psi_size];
  int gc[psi_size];
  fftidx_to_gidx(ga,gb,gc,nfft);

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  complex total_contribution = (0.0,0.0);

  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  int n_np = gw_sigma->num_sig_matels;
  int *n_list = gw_sigma->n_list_sig_matels;
  int *np_list = gw_sigma->np_list_sig_matels;

  complex *f = new complex[n_np*psi_size];
  std::vector<int> map(psi_size);
  for (int k = 0; k < K; k++) {
    int epsilon_size = 0;

    for(int g=0;g<psi_size;g++){
      if(accept_v[g]){
        map[epsilon_size] = g;
        epsilon_size++;
      }
    }
    map.resize(epsilon_size);

    int base_index = k*2*n*psi_size;


//This could probably be done once per node and cached
    for (int i = 0; i < f_cache->getNSize(); i++){
      int i_index = n_list[i]-1;
      int j_index = np_list[i]-1;

      int state_index = i*2*psi_size;
      int f_base = i*psi_size;

      for(int g=0; g < psi_size; g++){
        f[f_base+g] = states[base_index + state_index + g].conj() * states[base_index + state_index + g];
      }


      fft_controller->setup_fftw_3d(nfft, -1);
      fftw_complex* in_pointer = fft_controller->get_in_pointer();
      fftw_complex* out_pointer = fft_controller->get_out_pointer();
      // Pack our data, do the fft, then get the output
      put_into_fftbox(nfft, &f[f_base], in_pointer);
      fft_controller->do_fftw();
      fftbox_to_array(psi_size, out_pointer, &f[f_base], 1);
    }

    for (int i = 0; i < f_cache->getNSize(); i++) {
      int f_base = i*psi_size;
      complex contribution = (0.0,0.0);
      for (int r = 0; r < config.tile_rows; r++) {
        int g1 = thisIndex.x*eps_rows+r;
        for (int c = 0; c < config.tile_cols; c++) {
          int g2 = thisIndex.y*eps_cols+c;
          if(g1>=epsilon_size || g2>=epsilon_size) continue;

          int gdiff[3];
          gdiff[0] = ga[map[g1]]-ga[map[g2]];
          gdiff[1] = gb[map[g1]]-gb[map[g2]];
          gdiff[2] = gc[map[g1]]-gc[map[g2]];
          // flip the value and
          // set back to gdiff values

          for (int ii=0; ii<3; ii++){
            if (gdiff[ii] < -nfft[ii]/2){
              gdiff[ii] += nfft[ii];
            }
            if (gdiff[ii] >= nfft[ii]/2){
              gdiff[ii] -= nfft[ii]/2;
            }
          }

          int gdiffIndex = -1;
          for (int ii=0; ii<psi_size; ii++){
            if (gdiff[0]==ga[ii] && gdiff[1]==gb[ii] && gdiff[2]==gc[ii]){
              gdiffIndex = ii;
              break;
            }
          }

          contribution += f[f_base+gdiffIndex]*data[IDX_eps(r,c)];
        }
      }
      contrib_data[ik] = contribution;
      tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
      ik++;
      total_contribution += contribution;
    }
  } //end of K loop

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::cohComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
  delete[] f;
}


void EpsMatrix::bareExchange() {
  complex total_contribution = (0.0,0.0);
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  std::vector<double> vcoulb = psi_cache->getVCoulb();

  if(qindex==0)
    vcoulb[0] = psi_cache->getVCoulb0();

  if(thisIndex.x == thisIndex.y) {
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < f_cache->getNSize(); i++) {//ib = 5 to 8 actually, map to a number from 0
        complex contribution = (0.0,0.0);
        for (int l = 0; l < L; l++) {
          complex* f = f_cache->getFVec(k, i, l, thisIndex.x, eps_rows);
          for(int ii=0; ii < config.tile_rows; ii++) {
            int g = thisIndex.x*eps_rows+ii;
            if(g < vcoulb.size())
              contribution += f[ii]*f[ii].conj()*vcoulb[g];
          }
        }
        contrib_data[ik] = contribution * -1.0;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution * -1.0;
      }
    }
  }
  else{
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < f_cache->getNSize(); i++) {
        complex contribution = (0.0,0.0);
        tuple_reduction[ik++] =  CkReduction::tupleElement(sizeof(complex), &contribution, CkReduction::sum_double);
      }
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::bareExchangeComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
}

void EpsMatrix::coh(){

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

  complex* states = psi_cache->getStates();
  int psi_size = nfft[0]*nfft[1]*nfft[2];
  std::vector<int> accept_v = f_cache->getAcceptVector();

  int ga[psi_size];
  int gb[psi_size];
  int gc[psi_size];
  fftidx_to_gidx(ga,gb,gc,nfft);

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  complex total_contribution = (0.0,0.0);

  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  int n_np = gw_sigma->num_sig_matels;
  int *n_list = gw_sigma->n_list_sig_matels;
  int *np_list = gw_sigma->np_list_sig_matels;

  complex *f = new complex[n_np*psi_size];
  std::vector<int> map(psi_size);
  for (int k = 0; k < K; k++) {
    int epsilon_size = 0;

    for(int g=0;g<psi_size;g++){
      if(accept_v[g]){
        map[epsilon_size] = g;
        epsilon_size++;
      }
    }
    map.resize(epsilon_size);

    int base_index = k*2*n*psi_size;


//This could probably be done once per node and cached
    for (int i = 0; i < f_cache->getNSize(); i++){
      int i_index = n_list[i]-1;
      int j_index = np_list[i]-1;

      int state_index = i*2*psi_size;
      int f_base = i*psi_size;

      for(int g=0; g < psi_size; g++){
        f[f_base+g] = states[base_index + state_index + g].conj() * states[base_index + state_index + g];
      }


      fft_controller->setup_fftw_3d(nfft, -1);
      fftw_complex* in_pointer = fft_controller->get_in_pointer();
      fftw_complex* out_pointer = fft_controller->get_out_pointer();
      // Pack our data, do the fft, then get the output
      put_into_fftbox(nfft, &f[f_base], in_pointer);
      fft_controller->do_fftw();
      fftbox_to_array(psi_size, out_pointer, &f[f_base], 1);
    }

    for (int i = 0; i < f_cache->getNSize(); i++) {
      int f_base = i*psi_size;
      complex contribution = (0.0,0.0);
      for (int r = 0; r < config.tile_rows; r++) {
        int g1 = thisIndex.x*eps_rows+r;
        for (int c = 0; c < config.tile_cols; c++) {
          int g2 = thisIndex.y*eps_cols+c;
          if(g1>=epsilon_size || g2>=epsilon_size) continue;

          int gdiff[3];
          gdiff[0] = ga[map[g1]]-ga[map[g2]];
          gdiff[1] = gb[map[g1]]-gb[map[g2]];
          gdiff[2] = gc[map[g1]]-gc[map[g2]];
          // flip the value and
          // set back to gdiff values

          for (int ii=0; ii<3; ii++){
            if (gdiff[ii] < -nfft[ii]/2){
              gdiff[ii] += nfft[ii];
            }
            if (gdiff[ii] >= nfft[ii]/2){
              gdiff[ii] -= nfft[ii]/2;
            }
          }

          int gdiffIndex = -1;
          for (int ii=0; ii<psi_size; ii++){
            if (gdiff[0]==ga[ii] && gdiff[1]==gb[ii] && gdiff[2]==gc[ii]){
              gdiffIndex = ii;
              break;
            }
          }

          contribution += f[f_base+gdiffIndex]*data[IDX_eps(r,c)];
        }
      }
      contrib_data[ik] = contribution;
      tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
      ik++;
      total_contribution += contribution;
    }
  } //end of K loop

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::cohComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
  delete[] f;
}

void EpsMatrix::findAlpha() {
  if (config.chareCols() != 1) {
    CkAbort("findAlpha() only implemented for 1D decompositions\n");
  }
  double R = 0;
  for(int i = 0; i < config.tile_cols; i++) {
    R += abs(data[i]);
  }

  CkCallback cb(CkReductionTarget(Controller, found_alpha), controller_proxy);
  contribute(sizeof(long double), &R, CkReduction::max_double, cb);
}

void EpsMatrix::convergence_check(CProxy_EpsMatrix cproxy){
    // TODO: memcpy
    std::vector<complex> data_out(total_data);
    for(int i=0;i<total_data;i++)
      data_out[i] = data[i];

    cproxy(thisIndex.x, thisIndex.y).receiveConvCheck(data_out);
}

void EpsMatrix::receiveConvCheck(std::vector<complex> data_in) {
  double Rmax=0;  // the largest element
  double tmp;
  for(int i=0; i<total_data; i++) {
    tmp = abs(data[i] - data_in[i]);
    if( tmp > Rmax ){ Rmax = tmp; }
  }
  contribute(sizeof(complex), &Rmax, CkReduction::max_double,
      CkCallback(CkReductionTarget(Controller, converge_results), controller_proxy));
}

void EpsMatrix::createTranspose(CProxy_EpsMatrix other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
      if (todo) {
        complex tranpose = data[IDX_eps(i,j)];
        tranpose.im *= -1;
        incoming.push_back(tranpose);
      } else {
        incoming.push_back(data[IDX_eps(i,j)]);
      }
    }
  }
  if(todo) {
    other(thisIndex.y, thisIndex.x).receiveTranspose(incoming);
  } else {
    other(thisIndex.x, thisIndex.y).receiveTranspose(incoming);
  }
}

void EpsMatrix::receiveTranspose(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, transpose_complete), controller_proxy));
}

void EpsMatrix::transferToGpp(CProxy_Gpp other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
        incoming.push_back(data[IDX_eps(i,j)]);
    }
  }
    other(thisIndex.x, thisIndex.y).recvCopy(incoming);
}

void EpsMatrix::transferFromGpp(CProxy_Gpp other, bool todo) {
  PsiMessage* msg;
  int tile_size = eps_rows*eps_cols;
  complex* tile_data = new complex[tile_size];
  msg = new (tile_size) PsiMessage(tile_size, tile_data);
  msg = other(thisIndex.x, thisIndex.y).send_data(msg);
  int idx = 0;
  // CkPrintf("GPP transfer 2 tile %d %d %d %d %d %d\n", CKMYPE(), msg->k_index, msg->spin_index, msg->state_index, config.tile_rows, config.tile_cols);
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
      data[IDX_eps(i,j)] = msg->psi[idx];
      idx++;
    }
  }
  contribute(CkCallback(CkReductionTarget(Controller, transferfromGpp_complete), controller_proxy));
    // other(thisIndex.x, thisIndex.y).recvCopy(incoming);
}

void EpsMatrix::createCopy(CProxy_EpsMatrix other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
        incoming.push_back(data[IDX_eps(i,j)]);
    }
  }
    other(thisIndex.x, thisIndex.y).recvCopy(incoming);
}

void EpsMatrix::recvCopy(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, copy_complete), controller_proxy));
}

void EpsMatrix::createConjugate(CProxy_EpsMatrix other){
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++)
    for(int j=0; j < config.tile_cols; j++){
      complex conj = data[IDX_eps(i,j)].conj();
      incoming.push_back(conj);
    }
  other(thisIndex.x, thisIndex.y).receiveConjugate(incoming);
}

void EpsMatrix::receiveConjugate(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];
    
  contribute(CkCallback(CkReductionTarget(Controller, conjugateComplete), controller_proxy));
}

void EpsMatrix::multiply_coulb(){
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  std::vector<double> coulb = psi_cache->getVCoulb();
  if(qindex==0)
    coulb[0] = psi_cache->getVCoulb0();

  for(int i=0;i<config.tile_rows;i++){
    for(int j=0;j<config.tile_cols;j++){
      int g = thisIndex.x*config.tile_rows+i;
      int gp = thisIndex.y*config.tile_cols+j;
      if(g==gp && g<coulb.size())
        data[IDX_eps(i,j)] -= 1.0;
      if(g<coulb.size() && gp<coulb.size())
        data[IDX_eps(i,j)] *= sqrt(coulb[g])*sqrt(coulb[gp]);
    }
  }

  contribute(CkCallback(CkReductionTarget(Controller, s_ready), controller_proxy));
}


void EpsMatrix::print_col(int num) {
  char fname[100];
  
  int x = thisIndex.x;
  int y = thisIndex.y;
  sprintf(fname,"debug/S_x%d_y%d_c%d.dat", x, y, num); 
  int global_row, global_col;
  int i = 0;
  FILE* fp = fopen(fname, "w");
  for (int r = 0; r < config.tile_rows; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
        global_row = x*config.tile_rows + r;
        global_col = y*config.tile_cols + c;
        if ( global_col == num ) {
          double re = data[c*config.tile_rows + r].re;
          double im = data[c*config.tile_rows + r].im;
          fprintf(fp, "%d %d %d %d %.8e %.8e\n", x, y, global_row, global_col, re, im);
          // CkPrintf("[EPSMATRIX] x/y %d %d global_r/c %d %d val %.8e\n", x, y, global_row, global_col, val);
        }
      i++;
    }
  }
  fclose(fp);
}
void EpsMatrix::print(int qindex, int fnum) {
  char fname[100];
  int x = thisIndex.x;
  int y = thisIndex.y;
  sprintf(fname,"./debug/EPS_%d_q%d_x%d_y%d.dat", fnum, qindex, x, y); 
  int global_row, global_col;
  int i = 0;
  FILE* fp = fopen(fname, "w");
  for (int r = 0; r < config.tile_rows; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
        global_row = x*config.tile_rows + r;
        global_col = y*config.tile_cols + c;
        // double re = data[c*config.tile_rows + r].re;
        // double im = data[c*config.tile_rows + r].im;
        double re = data[IDX_eps(r,c)].re;
        double im = data[IDX_eps(r,c)].im;
        fprintf(fp, "%d %d %d %d %.8e %.8e\n", x, y, global_row, global_col, re, im);
        // CkPrintf("[EPSMATRIX] x/y %d %d global_r/c %d %d val %.8e\n", x, y, global_row, global_col, val);
      i++;
    }
  }
  fclose(fp);
}

void EpsMatrix::print_row(int num) { 
  int x = thisIndex.x;
  int y = thisIndex.y;
  int global_row, global_col;
  for (int r = 0; r < config.tile_rows; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
        global_row = x*config.tile_rows + r;
        global_col = y*config.tile_cols + c;
        if ( global_row == num ) {
          double val = data[c*config.tile_rows + r].re;
          //kayahan debug
          CkPrintf("[EPSMATRIX] x/y %d %d global_r/c %d %d val %.8e\n",x, y, global_row, global_col, val);
        }
    }
  }
}

DiagMessage* EpsMatrix::sendDataSimple(DiagMessage* msg) {
  
  int msg_cols = msg->cols;
  int msg_rows = msg->rows;
  // int eps_local_cols = msg_rows;
  // int eps_local_rows = msg_cols;
  int eps_local_cols = msg_cols;
  int eps_local_rows = msg_rows;
  
  int idx_col, idx_row;
  int i = 0;

  int x = thisIndex.x;
  int y = thisIndex.y;
  int global_row, global_col;
  int global_index;

  for (int r = 0; r < config.tile_rows ; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
      // global_index = c*config.tile_rows + r;
        if (r < msg_rows && c < msg_cols) {
          // data[global_index] = msg->data[i];
          data[IDX_eps(r,c)] = msg->data[i];
          // data[r*config.tile_cols + c] = msg->data[c*msg_rows + r];
          i++;
        }
        else {
          data[global_index] = 0.0;
          // data[c*config.tile_rows + r] = msg->data[i];
        }
        // global_row = x*config.tile_rows + r;
        // global_col = y*config.tile_cols + c;
        // // if ( global_row == 17 ) {
        //   CkPrintf("[DIAGMESSAGE] rc %d %d msg_rc %d %d config_rc %d %d global_rc %d %d start_rc %d %d val %.8e valmsg %.8e\n",
        //    r, c, 
        //    msg_rows, msg_cols, 
        //    config.tile_rows, config.tile_cols,  
        //    global_row, global_col, 
        //    start_row, start_col, 
        //    data[c*config.tile_rows + r].re, msg->data[i-1].re);
        // // }
    }
  }
  return msg;
}

DiagMessage* EpsMatrix::receiveDataSimple(DiagMessage* msg) {
  msg->eps_pe = CkMyPe();
  msg->x = thisIndex.x;
  msg->y = thisIndex.y;

  bool borderX = false;
  bool borderY = false;
  if (thisIndex.x + 1 == numBlocks) {
    borderX = true;
  }
  if (thisIndex.y + 1 == numBlocks) {
    borderY = true;
  }
  int real_epsilon_size = msg->eps_size;
  unsigned int dataSize = 0;
  int remElems2 = real_epsilon_size % eps_rows;  // eps_rows = eps_col square matrix
  int stdElems = eps_rows * eps_cols;
  int remElems = remElems2 * eps_rows;
  int cornerElems = remElems2 * remElems2;

  int rows = 0;
  int cols = 0;
  if (borderX && !borderY) {
    dataSize = remElems;
    rows = remElems2;
    cols = eps_rows;
  } else if (!borderX && borderY) {
    dataSize = remElems;
    rows = eps_rows;
    cols = remElems2;
  } else if (borderX && borderY) {
    dataSize = cornerElems;
    rows = remElems2;
    cols = remElems2;
  } else {
    dataSize = stdElems;
    rows = eps_rows;
    cols = eps_rows;
  }
  msg->size = dataSize;
  msg->rows = rows;
  msg->cols = cols;

  int idx_col, idx_row;
  int i = 0;

  // Transfer data to column major
  for (int c = 0; c < cols; c++) {
    for (int r = 0; r < rows; r++) {
      idx_row = start_row + r;
      idx_col = start_col + c;
      
      // CHARM++ data is row-major
      msg->data[i] = data[r*config.tile_cols + c];
      // msg->data[i].imag(data[r*config.tile_cols + c].im);
      i++;
    }
  }
  
  return msg;
}

void DiagBridge::prepareData(int qindex, int eps_size, int num_qpts) {
  // This routine sets up the containers (diagData) for each element of process grid
  // 2D Block cyclic mapping
  // https://www.netlib.org/scalapack/slug/node76.html, Fig. 4.6
  // Blocks are square and block sizes are equal to EpsMatrix tile size in charm++
  // eps_size: rank of the matrix for an N x N invertable matrix size is N
  // qindex: q index of the epsilon. Not really needed here, but used in debugging

  int mype = CkMyPe();
  
  // Possible Charm++ tile sizes 
  int remElems2 = eps_size % eps_rows;  
  int stdElems = eps_rows * eps_rows; // square matrix 
  int remElems = remElems2 * eps_rows; // sides of eps
  int cornerElems = remElems2 * remElems2; // corner of eps

  // Proc distribution to be used for diagonalization
  GWBSE* gwbse = GWBSE::get();
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  // TODO (kayahans): read from input
  proc_rows = gw_sigma->proc_rows;
  proc_cols = gw_sigma->proc_cols;

  // Number of blocks
  numBlocks = eps_size / eps_rows + 1;
  // Eps is a square matrix
  int numx = numBlocks;
  int numy = numBlocks;

  // Dimensions of the matrix stored on process grid
  row_size = 0;
  col_size = 0;
  totaldata = 0; // totaldata = row_size * col_size

  // These are used to count the number of charm++ tiles in a process grid block
  int x_prev = -1;
  int y_prev = -1;
  int x_num_tiles = 0; // num of charm++ tiles in row dimension
  int y_num_tiles = 0; // num of charm++ tiles in col dimension

  
  for (int x=0; x < numx; x++) {
    for (int y=0; y < numy; y++) {
      int dest_pe_row = x%proc_rows;
      int dest_pe_col = y%proc_cols;
      int dest_pe = dest_pe_row*proc_cols + dest_pe_col;

      if (dest_pe == mype) {
        if (x > x_prev) {
          x_prev = x;
          x_num_tiles += 1;
        }
        if (y > y_prev) {
          y_prev = y;
          y_num_tiles +=1;
        }
        bool borderX = false;
        bool borderY = false;
        if (x + 1 == numBlocks) {
          borderX = true;
        }

        if (y + 1 == numBlocks) {
          borderY = true;
        }

        int xy_datasize = 0;
        int rows = 0;
        int cols = 0;
        if (borderX && !borderY) {
          xy_datasize = remElems;
          rows = remElems2;
          cols = eps_rows;
        } else if (!borderX && borderY) {
          xy_datasize = remElems;
          rows = eps_rows;
          cols = remElems2;
        } else if (borderX && borderY) {
          xy_datasize = cornerElems;
          rows = remElems2;
          cols = remElems2;
        } else {
          xy_datasize = stdElems;
          rows = eps_rows;
          cols = eps_rows;
        }
        totaldata += xy_datasize;
        row_size  += rows;
        col_size  += cols;
      } // end if
    } // end for
  } // end for

  row_size = row_size / y_num_tiles;
  col_size = col_size / x_num_tiles;

  // Setup the container to be transferred to MPI
  diagData = new diagData_t;
  diagData->qindex = qindex;
  diagData->num_handoffs = num_qpts;
  
  diagData->inputsize = totaldata; 
  diagData->row_size = row_size;
  diagData->col_size = col_size;
  diagData->nprow = proc_rows;
  diagData->npcol = proc_cols;
  diagData->nb = eps_rows;
  diagData->n = eps_size;

  diagData->input = new std::complex<double>[totaldata];
  // TODO (kayahans): Not sure which pe gets the final result, for now allocate this in all 
  // Later when we decide how to distribute eigenvectors/values, we can make this smarter
  // diagData->eig_e = new std::complex<double>[eps_size];
  diagData->eig_e = new double[eps_size];
  // diagData->eig_v = new std::complex<double>[eps_size*eps_size];  
  // diagData->eig_v = new matel[eps_size*eps_size];  
  // kayahan debug
  // CkPrintf("[DIAGONALIZER] Created a pointer with totalsize %d numblocks %d for S matrix global dim %d at pe %d for qindex %d\n", totaldata, numBlocks, eps_size, CkMyPe(), qindex);
  contribute(CkCallback(CkReductionTarget(Controller, diag_setup), controller_proxy));
}



#include "eps_matrix.def.h"
