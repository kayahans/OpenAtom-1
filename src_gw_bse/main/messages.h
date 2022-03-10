#ifndef __MESSAGES_H__
#define __MESSAGES_H__

#include <cstdlib>
#include "ckcomplex.h"
#include "diagonalizer.h"
#include "messages.decl.h"

// Message sent from psi used to compute f. Some of these are cached on each
// node, and others are streamed in to the PMatrix as needed.
class PsiMessage : public CMessage_PsiMessage {
  public:
    PsiMessage(unsigned s, complex* p) : size(s) {
      std::copy(p, p+size, psi);
    }
    unsigned spin_index, k_index, state_index, size;
    bool shifted;
    bool sigma;
    complex* psi;
};

class GppVMessage : public CMessage_GppVMessage {
  public:
    GppVMessage(unsigned s, complex* _eigv) : size(s) {
      std::copy(_eigv, _eigv+size, eigv);
      // std::copy(_eige, _eige+size, eige);
    }
    unsigned spin_index, q_index, alpha_idx, size;
    complex* eigv;
    // double* eige;
};

class GppEMessage : public CMessage_GppEMessage {
  public:
    GppEMessage(unsigned s, double* _eige) : size(s) {
      // std::copy(_eigv, _eigv+size, eigv);
      std::copy(_eige, _eige+size, eige);
    }
    unsigned spin_index, q_index, size;
    // complex* eigv;
    double* eige;
};

// Message sent between PMatrix chares to exchange data during a transpose.
class TransposeMessage : public CMessage_TransposeMessage {
  public:
    complex* data;
    unsigned start_col;
};


class Phase2Message : public CMessage_Phase2Message {
  public:
    complex* data;
    int size;
    int global_x;
    int global_y;
};

class Phase3Message : public CMessage_Phase3Message {
  public:
    complex* data;
    int size;
    int global_x;
    int global_y;
    int start_i;
    int start_j;
    int end_i;
    int end_j;
};

class Phase4Message : public CMessage_Phase3Message {
  public:
    complex* data;
    int n;
};

class DiagMessage : public CMessage_DiagMessage {
  public:
  // TODO (kayahans): Using fixed size messages (simpler) 
  //                  400 is 20 x 20 that is eps_cols x eps_rows (global variables)
  //                  In any way, the actual data transferred will be either
  //                  equal to or smaller than this, later here can be generalized
    DiagMessage() {
    }
    complex data[400];
    double eigenvalues[20];
    int size;
    int eps_size;
    int x;
    int y;
    int eps_pe;
    int rows;
    int cols;
    
};
#endif
