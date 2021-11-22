#include "standard_include.h"
#include "allclass_gwbse.h"
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
// diagData_t* diagData;

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

void EpsMatrix::packMsg(DiagMessage* msg) {
  // int msg_idx = 0;
  // int local_idx = (msg->row - start_row) * config.tile_cols + (msg->col - start_col);
  // for (int r = 0; r < msg->num_rows; r++) {
  //   memcpy(&msg->data[msg_idx],&data[local_idx],sizeof(complex)*msg->num_cols);
  //   msg_idx += msg->num_cols;
  //   local_idx += config.tile_cols;
  // }
}

void EpsMatrix::sendDiagData(DiagMessage* msg) {
  // printf("eps x y %d %d \n", thisIndex.x, thisIndex.y);

  int i = 0;
  int idx_row;
  int idx_col;

  for (int r = 0; r < eps_rows; r++) {
      for (int c = 0; c < eps_cols; c++) {
        idx_row = start_row + r;
        idx_col = start_col + c;
        msg->data[i] = data[r*config.tile_cols + c].re;
        printf("x %d y %d r %d c %d i %d %.6e\n", thisIndex.x, thisIndex.y, r, c, i, data[r*config.tile_cols + c].re);
        i++;
      }
  }
  printf("Before copy\n");
  for (int i = 0; i < 400; i++) {
      printf("i %d val %.6e\n", i, msg->data[i]);
  }
}

void EpsMatrix::copyToMPI(int qindex, int epsilon_size) {
  int x = 2;
  int y = 5;
  int i = 0;
  int idx_row;
  int idx_col;
  int numblocks = 7;
  int size = 400;

  // if (thisIndex.x == x && thisIndex.y == y) {
  //   printf("%.6e\n", thisProxy(2, 5).data[0]);
  // }
  // printf("pe_eps %d x %d y %d \n", CKMYPE(), x,y);
  // if (thisIndex.x == 0 && thisIndex.y == 0) {
  //   i = 0;
  //   DiagMessage* msg; 
  //   msg = new DiagMessage(size);
  //   thisProxy(x,y).sendDiagData(msg);
  //   for (int i = 0; i < 20; i++) {
  //     printf("isend %d val %.6e\n", i, msg->data[i]);
  //   }
  // }

  // if (thisIndex.x == x && thisIndex.y == y) {
  //   i = 0;
  //   for (int r = 0; r < eps_rows; r++) {
  //     for (int c = 0; c < eps_cols; c++) {
  //       idx_row = start_row + r;
  //       idx_col = start_col + c;
  //       if (i < 20) {
  //         printf("i %d val %.6e\n", i, data[r*config.tile_cols + c].re);
  //       }
  //       i++;
  //     }
  //   }
  // }
  if (thisIndex.x == x && thisIndex.y == y) {
    for (int r = 0; r < eps_rows; r++) {
      for (int c = 0; c < eps_cols; c++) {
        idx_row = start_row + r;
        idx_col = start_col + c;
        if (i < 20) {
          printf("pe_eps %d x %d y %d r %d c %d val %.6e\n", CKMYPE(), x,y,r, c, data[r*config.tile_cols + c].re);
        }
        i++;
      }
    }
  }
  if (thisIndex.x == x && thisIndex.y == y) {
    
    DiagMessage* msg;
    int size = 400;
    msg = new DiagMessage(size);
    msg->x = x;
    msg->y = y;
    int dest_pe_row = x%proc_rows;
    int dest_pe_col = y%proc_cols;
    int eps_dest_pe = dest_pe_row*proc_cols + dest_pe_col;
    // msg->dest_pe_row = dest_pe_row;
    // msg->dest_pe_col = dest_pe_col;
    msg->dest_pe     = eps_dest_pe;
    // msg->total_data  = dataSize;
    msg->eps_source_pe = CKMYPE();
    msg->eps_source_pe2 = x*numblocks+y;
    i = 0;
    for (int r = 0; r < eps_rows; r++) {
      for (int c = 0; c < eps_cols; c++) {
        idx_row = start_row + r;
        idx_col = start_col + c;
        msg->data[i] = data[r*config.tile_cols + c].re;
        i++;
      }
    }
    for (int i = 0; i < 20; i++) {
      printf("isend %d val %.6e\n", i, msg->data[i]);
    }
    // packMsg(msg);
    // diag_bridge_proxy.ckLocalBranch()->receiveDataSimple(msg);
    diag_bridge_proxy[eps_dest_pe].receiveDataSimple(msg);
  }
    // msg = new DiagMessage(size);
    // thisProxy(0, 0).sendDiagData(msg);
    // printf("Other routine\n");
    // for (int i = 0; i < size; i++) {
    //   printf("i %d val %.6e\n", i, msg->data[i]);
    // }
    // for (int r = 0; r < eps_rows; r++) {
    //   for (int c = 0; c < eps_cols; c++) {
    //     idx_row = 0 + r;
    //     idx_col = 0 + c;
    //     DiagMessage* msg;
    //     if (idx_row < epsilon_size && idx_col < epsilon_size) {
    //       thisProxy(0, 1).sendDiagData();
    //       // printf("x %d y %d r %d c %d val %.6e\n", x,y,r, c, thisProxy(0,0)->data[r*config.tile_cols + c].re);
    //       i++;
    //     }
    //   }
    // }
  // }
  // if (thisIndex.x == 0 && thisIndex.y == 0) {
  //   int numblocks = 7;
  //   for (int x = 0; x < numblocks; x++) {
  //     for (int y = 0; y < numblocks; y++) {
  //       GWBSE* gwbse = GWBSE::get();
  //       proc_rows = 2; // gwbse->gw_parallel.proc_rows;
  //       proc_cols = 2; // gwbse->gw_parallel.proc_cols;
  //       int dest_pe_row = x%proc_rows;
  //       int dest_pe_col = y%proc_cols;
  //       int eps_dest_pe = dest_pe_row*proc_cols + dest_pe_col;

  //       bool borderX = false;
  //       bool borderY = false;

  //       if (x + 1 == numBlocks) {
  //         borderX = true;
  //       }

  //       if (y + 1 == numBlocks) {
  //         borderY = true;
  //       }

  //       int dataSize = -1;

  //       int remElems2 = epsilon_size % eps_rows;  // eps_rows = eps_col square matrix
  //       int stdElems = eps_rows * eps_cols;
  //       int remElems = remElems2 * eps_rows;
  //       int cornerElems = remElems2 * remElems2;

  //       int rows = 0;
  //       int cols = 0;
  //       if (borderX && !borderY) {
  //         dataSize = remElems;
  //         rows = remElems2;
  //         cols = eps_rows;
  //       } else if (!borderX && borderY) {
  //         dataSize = remElems;
  //         rows = eps_rows;
  //         cols = remElems2;
  //       } else if (borderX && borderY) {
  //         dataSize = cornerElems;
  //         rows = remElems2;
  //         cols = remElems2;
  //       } else {
  //         dataSize = stdElems;
  //         rows = eps_rows;
  //         cols = eps_rows;
  //       }

  //       DiagMessage* msg = new DiagMessage();
  //       msg->x = x;
  //       msg->y = y;
  //       msg->dest_pe_row = dest_pe_row;
  //       msg->dest_pe_col = dest_pe_col;
  //       msg->dest_pe     = eps_dest_pe;
  //       msg->total_data  = dataSize;
  //       msg->eps_source_pe = CKMYPE();
  //       msg->eps_source_pe2 = x*numblocks+y;

  //       packMsg(msg);
  //       diag_bridge_proxy[eps_dest_pe].receiveDataSimple(msg);
  //     }
  //   }
    
  // }
  // GWBSE* gwbse = GWBSE::get();
  // proc_rows = 2; // gwbse->gw_parallel.proc_rows;
  // proc_cols = 2; // gwbse->gw_parallel.proc_cols;
  // int dest_pe_row = thisIndex.x%proc_rows;
  // int dest_pe_col = thisIndex.y%proc_cols;
  // int eps_dest_pe = dest_pe_row*proc_cols + dest_pe_col;


  // int x = thisIndex.x;
  // int y = thisIndex.y;

  // bool borderX = false;
  // bool borderY = false;

  // if (x + 1 == numBlocks) {
  //   borderX = true;
  // }

  // if (y + 1 == numBlocks) {
  //   borderY = true;
  // }

  // int dataSize = -1;

  // int remElems2 = epsilon_size % eps_rows;  // eps_rows = eps_col square matrix
  // int stdElems = eps_rows * eps_cols;
  // int remElems = remElems2 * eps_rows;
  // int cornerElems = remElems2 * remElems2;

  // int rows = 0;
  // int cols = 0;
  // if (borderX && !borderY) {
  //   dataSize = remElems;
  //   rows = remElems2;
  //   cols = eps_rows;
  // } else if (!borderX && borderY) {
  //   dataSize = remElems;
  //   rows = eps_rows;
  //   cols = remElems2;
  // } else if (borderX && borderY) {
  //   dataSize = cornerElems;
  //   rows = remElems2;
  //   cols = remElems2;
  // } else {
  //   dataSize = stdElems;
  //   rows = eps_rows;
  //   cols = eps_rows;
  // }


  // std::vector<complex> data_out(dataSize);

  // int idx_row = 0;
  // int idx_col = 0;
  // int i = 0;
  // for (int r = 0; r < eps_rows; r++) {
  //   for (int c = 0; c < eps_cols; c++) {
  //     idx_row = start_row + r;
  //     idx_col = start_col + c;
  //     if (idx_row < epsilon_size && idx_col < epsilon_size) {
  //       data_out[i] = data[r*config.tile_cols + c].re;  // True?
  //       i++;
  //     }
  //   }
  // }
  // std::vectoqr<complex> data_out(total_data);
  // for(int i=0;i<total_data;i++)
  //   data_out[i] = data[i];
  // CkPrintf("qindex %d, eps_size %d \n", qindex, epsilon_size);
  // CkPrintf("dest_pe %d thisIndex.x %d thisIndex.y %d rows %d cols %d total_data %d \n", dest_pe, thisIndex.x, thisIndex.y, ng, ng, dataSize);
  // CkPrintf("dest_pe %d thisIndex.x %d thisIndex.y %d total_data %d \n", dest_pe, thisIndex.x, thisIndex.y, dataSize);
  // CkPrintf(" cols %d eps %d %d %d %d %d %d\n", config.tile_cols, eps_rows, blockSize, numBlocks, block, start_row, start_col);
  // diag_bridge_proxy[dest_pe].prepareData();
  // int eps_source_pe = CKMYPE();
  // CkPrintf("[EPSMATRIX] Eps chare (%d,%d) will copy global (%d, %d)-(%d, %d) "
  //   "S=%d eps_source_pe %d eps_dest_pe %d\n",
  //   x, y, start_row, start_col, idx_row, idx_col, dataSize, eps_source_pe, eps_dest_pe);

  
  // diag_bridge_proxy.receiveData(x, y, data_out, dataSize, rows, cols, dest_pe, sending_pe);
  // diag_bridge_proxy[eps_dest_pe].receiveDataSimple(x, y, eps_source_pe, eps_dest_pe, start_row, start_col, idx_row, idx_col );

  // CkCallback cb(CkIndex_EpsMatrix::all_copied(NULL), thisProxy);
  // int contribution = 1;
  // contribute(sizeof(int), &contribution, CkReduction::sum_int, cb);
}

void EpsMatrix::all_copied(CkReductionMsg *msg) {
  int total = *(int *)msg->getData();
  if (total == 49) {
    CkPrintf("[EPSMATRIX] All chares are copied\n");
    // contribute(CkCallback(CkReductionTarget(Controller, mpi_copy_complete), controller_proxy));
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

void EpsMatrix::prep_diag(){
  CkPrintf("\nDIAG data %d \n", CKMYPE());
  // diagData = new diagData_t();
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

#include "eps_matrix.def.h"
