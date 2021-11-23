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
    DiagMessage(int s) : size(s) {
      data = new double[size];
    }
    double* data;
    int size;
    int x;
    int y;
    int dest_pe_row;
    int dest_pe_col;
    int dest_pe;
    int total_data;
    int eps_source_pe;
    int eps_source_pe2;
};
#endif
