//////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2020 OpenAtom developers.
//
// File developed by: Kayahan Saritas, saritaskayahan@gmail.com, Yale University
//
//
// File created by: Kayahan Saritas, saritaskayahan@gmail.com, Yale University
//////////////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef UTILS_WINDOWING_H_
#define UTILS_WINDOWING_H_

#include <vector>
#include <string>
#include <cassert>
#include "./vector_util.h"
#include "./gl.h"
#include "standard_include.h"
#include "allclass_gwbse.h"

struct WINPAIR {
  double gap;
  double bw;
  double w1[2];  // First window slice
  double w2[2];  // Second window slice
  double zeta;
  double gamma;
  double errfrac;
  bool isHgl;
  bool sigma;  // is this a window pair for sigma calculations
  int sigma_index;  // Sigma has + and - calculations as in Eq. 29 in PRB 101 035139 (2020)
                    // 0 means positive (valence) 1 means negative (conduction)
                    // These are default to false and zero

  int ngl;
  std::vector<double> nodes;
  std::vector<double> weights;
  WINPAIR(double a_min, double a_max, double b_min, double b_max, double _gap, double _bw, int _ngl, double _errfrac, bool _isHgl, bool _sigma = false, int _sigma_index = 0) {
    GL gl;
    assert(a_max > a_min);
    assert(b_max > b_min);
    // assert(_errfrac > 0.);
    // (kayahans) Check validity for overlapping windows
    gap   = _gap;
    bw    = _bw;
    w1[0] = a_min;
    w1[1] = a_max;
    w2[0] = b_min;
    w2[1] = b_max;
    isHgl = _isHgl;
    if (isHgl) {
      zeta = 0.0367493;  // FIXME kayahans hard coded!
    } else {
      zeta  = 1./sqrt(gap*bw);
    }
    
    errfrac = _errfrac;
    ngl   = _ngl;
    gamma = 0.0; // FIXME

    sigma = _sigma;
    sigma_index = _sigma_index;

    if (ngl > gl.nptmax) {
      printf("For the window pair, number of required quadrature points is %d, but max available is %d. Max will be used!\n", ngl, gl.nptmax);
      ngl = gl.nptmax;
    }
    nodes.clear();
    weights.clear();

    for (int i = 0; i < ngl; i++) {
        nodes.push_back(gl.n[ngl-1][i]);
        weights.push_back(gl.w[ngl-1][i]);
    }
    // printf("Fitted G=%f, W=%f, nnodes=%d, zeta=%f\n", gap, bw, ngl, zeta);
  }
  bool in_window(const double &energy, const int& window_number, const int& sigma_index = 0);  // If energy is in the window pair
  double find_zeta() const;
  void strip(const std::vector<double>& A, const std::vector<double>& B);
};

class WINDOWING {
 public:
  int nspin = 0;
  int nq = 0;
  int nkpt = 0;
  int nv = 0;
  int nc = 0;
  int counter = 0;
  double evmin, evmax;
  double ecmin, ecmax;
  double gap;
  double wppmax, wppmin;

  bool sigma;   // if true, then it means there are two sets of windows like Eq. 29 in PHYSICAL REVIEW B 101, 035139 (2020)
  int w_sets[2];

  int max_windows[2];
  std::vector<double> opt_awins_list[2];
  std::vector<double> opt_bwins_list[2];
  int opt_num_windows[2];
  std::vector<WINPAIR> winpairs;
  
  WINDOWING() {};
  WINDOWING(double*** e_occ, double*** e_unocc) {
    GWBSE *gwbse = GWBSE::get();
    int nocc = gwbse->gw_parallel.L;
    int nunocc = gwbse->gw_parallel.M;
    int nspin = 1; 
    int nkpt = gwbse->gw_parallel.K;
    printf("WIN %d %d %d %d\n", nocc, nunocc, nspin, nkpt);
    initialize(e_occ, e_unocc, nocc,nunocc, nspin, nkpt);
  }

  void sigma_win(const double w, double*** const omsq, int* const ng);
  void read_from_file();
  void initialize(double*** e_occ, double*** e_unocc, int nocc, int nunocc, int nspin, int nkpt);
  void printenergies() const;
  void printparameters() const;
  void searchwins(char*);
  void from_input(std::vector<double>, std::vector<double>);

 private:
  // Available options
  // Static polarizability: P0
  // Dynamic polarizability: Pw
  // Sigma
  std::string opt;
  static const int num_opt = 3;
  std::string available_opt[num_opt] = {"P0", "Pw", "Sigma"};
  
  std::string opt_option;  
  static const int num_opt_options = 2;
  std::string available_opt_options[num_opt_options] = {"constant_dos", "read_dos"};

  double omega;
  double ecut;
  
  // optimizer parameters
  double errfrac;
  double ptol;  // errfrac = ptol/100
  double min_gap;
  
  std::vector<double> eigs;
  std::vector<double> eigs_v;
  std::vector<double> eigs_c;
  std::vector<double> wpp;
  
  // void loadstates(const STATES* const* const* _psi, const SYSINFO& _sys);
  // void loadgpp(double*** omsq, int* ng, const SYSINFO& _sys);
  double fixedwin_costfun(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
  double fixedwin_fixeddos_costfun(const std::vector<double>, const std::vector<double>);
  double costfun(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
  double fixed_win_optimizer(std::vector<double>&, std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
  double Nlm_cost(double, double, double, double);
  double quad_cost(double, double);
  double hgl_quad_cost(double);
  double win_P0_cost(int, int, std::vector<double>(&)[2], std::vector<double>(&)[2]);
  double win_Pw_cost(int, int, std::vector<double>(&)[2], std::vector<double>(&)[2]);
  double win_Sigma_cost(int, int, std::vector<double>(&)[2], std::vector<double>(&)[2]);
  void add_padding(std::vector<double> (&)[2]);
  void set_winpairs();
};
#endif  // UTILS_WINDOWING_H_
