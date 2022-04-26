//////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2020 OpenAtom developers.
//
// File developed by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//
//
// File created by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "./windowing.h"
#include "./vector_util.h"
#include "./gl.h"
#include <math.h>
#include "standard_include.h"
#include "allclass_gwbse.h"

bool WINPAIR::in_window(const double &a, const int& w_number, const int& _sigma_index) {
  bool _inside = false;
  if (_sigma_index == sigma_index) {
    if (w_number == 0) {
      if ((a >= w1[0] && a <= w1[1])) {
        _inside = true;
      }
    } else if (w_number == 1) {
      if ((a >= w2[0] && a <= w2[1])) {
        _inside = true;
      }
    } else {
      CkPrintf("Window number error!\n");
    }
  }
  return _inside;
}

double WINPAIR::find_zeta() const {
  if (!isHgl) {

    double G = gap;     // min(Ec-Ev)
    double W = bw;   // max(Ec-Ev)
    int ng = nodes.size();      // number of Gauss-Laguerre quadrature node
    const GL gl;

    // list of u: logarithmic grid from 10^-2 to 10^2 with 10,000 points
    int nmax = 10000;
    std::vector<double> ulist = linspace(-2, 2, nmax);
    ulist = pow(10, ulist);
    
    double exact = 1;  // exact value of the integrand
    std::vector<double> errdiff(nmax);
    double sum_u, sum_uWG;
    for (int i = 0; i < nmax; i++) {
      // 1) integrand of u (=E/a)
      sum_u = 0;
      for (int ig = 0; ig < ng; ig++) {
        sum_u += gl.w[ng-1][ig]*exp(-(ulist[i]-1)*gl.n[ng-1][ig]);
      }
      sum_u *= ulist[i];
      // 2) integrand of u*W/G
      sum_uWG = 0;
      for (int ig = 0; ig < ng; ig++) {
        sum_uWG += gl.w[ng-1][ig]*exp(-(ulist[i]*W/G-1)*gl.n[ng-1][ig]);
      }
      sum_uWG *= ulist[i]*W/G;
      
      errdiff[i] = abs(sum_u - sum_uWG);
    }
    
    // find minimum error index
    int idx = 0;
    double minerrdiff = 1;
    for (int i = 0; i < nmax; i++) {
      // search only if ulist < 1
      if (ulist[i] < 1) {
        if (errdiff[i] < minerrdiff) {
          minerrdiff = errdiff[i];
          idx = i;
        }
      }
    }

    // set zeta
    return ulist[idx]/G;
  
  } else {
    CkPrintf("Only for GL quadrature!\n");
  }
}

void WINPAIR::strip(const std::vector<double>& A, const std::vector<double>& B) {
  // Update the gap and bandwidth values for a window pair after the A and B values are accumulated
  // This is not very useful at the moment because, we dont update the number of quadrature points 
  // as a result, but this is what Minjung did, so lets comply with it now.
  double Amin = *std::min_element(A.begin(),A.end());
  w1[0] = Amin;
  
  double Amax = *std::max_element(A.begin(),A.end());
  w1[1] = Amax;

  double Bmin = *std::min_element(B.begin(),B.end());
  w2[0] = Bmin;
  double Bmax = *std::max_element(B.begin(),B.end());
  w2[1] = Bmax;
  if (Bmin > Amax) {
    gap = Bmin-Amax;
    bw  = Bmax-Amin;
  } else if (Amin > Bmax) {
    gap = Amin-Bmax;
    bw  = Amax-Bmin;
  } else {
    gap = 0;
    std::vector<double> win_maxmins{Amin, Amax, Bmin, Bmax};
    double win_min = *std::min_element(win_maxmins.begin(),win_maxmins.end());
    double win_max = *std::max_element(win_maxmins.begin(),win_maxmins.end());
    bw = win_max - win_min;
  }
  // After bw and gap are updated, update zeta!
  if (!isHgl) 
    zeta = find_zeta();
  // FIXME add gamma here as well 
}  // end reoptimize

void WINDOWING::set_winpairs() {
    // (kayahans) valid for only P0
    double _gap, _bw;
    int _ngl;  // Number of GL nodes
    bool _isHgl;

    winpairs.clear();
    if (opt == "P0") {
        for (int i = 0; i < opt_num_windows[0]; i++) {
            for (int j = 0; j < opt_num_windows[1]; j++) {
                _gap = opt_bwins_list[0][j]   - opt_awins_list[0][i+1];
                _bw  = opt_bwins_list[0][j+1] - opt_awins_list[0][i];
                if (_gap > min_gap) {
                    _ngl = static_cast<int> (quad_cost(_bw, _gap));
                    _isHgl = false;
                } else {
                    _ngl = static_cast<int> (hgl_quad_cost(_bw));
                    _isHgl = true;
                }
                WINPAIR _pair(opt_awins_list[0][i],
                              opt_awins_list[0][i+1],
                              opt_bwins_list[0][j],
                              opt_bwins_list[0][j+1],
                              _gap, _bw, _ngl, errfrac,
                              _isHgl);
                winpairs.push_back(_pair);
            }
        }
    } else if (opt == "Sigma") {
        for (int cond = 0; cond < 2; cond++) {
            // If cond == 0, scanning over valence states and pp
            // if cond == 1, scanning over conduction states and pp
            for (int i = 0; i < opt_awins_list[cond].size() - 1 ; i++) {
                for (int j = 0; j < opt_bwins_list[cond].size() - 1; j++) {
                    _gap = opt_bwins_list[cond][j]   - opt_awins_list[cond][i+1];
                    _bw  = opt_bwins_list[cond][j+1] - opt_awins_list[cond][i];
                    if (abs(_gap) > abs(_bw)) {
                        _gap = abs(_gap);
                        _bw  = abs(_bw);
                        std::swap(_gap, _bw);
                    }
                    if (_gap > min_gap) {
                        _ngl = static_cast<int>(floor(quad_cost(_bw, _gap)));
                        _isHgl = false;
                    } else {
                        _ngl = static_cast<int>(floor(hgl_quad_cost(_bw)));
                        _isHgl = true;
                    }
                    WINPAIR _pair(opt_awins_list[cond][i],
                                  opt_awins_list[cond][i+1],
                                  opt_bwins_list[cond][j],
                                  opt_bwins_list[cond][j+1],
                                  _gap, _bw, _ngl, errfrac,
                                  _isHgl, true, cond);
                    winpairs.push_back(_pair);
                }
        }

        }
    }
}  // end set_pairs

void WINDOWING::initialize(double*** e_occ, double*** e_unocc, int _nocc, int _nunocc, int _nspin, int _nkpt) {
    nspin   = _nspin;
    nq      = _nkpt;  // FIXME
    nkpt    = _nkpt;
    double _eV       = 1./27.2114;
    ptol            = 10;                  // (hard coded)
    errfrac         = ptol/100;             // Percent tolerance converted to errfrac
    min_gap         = 0.1*_eV;              // (hard coded) Minimum gap to be considered metallic (0.1 eV)
    max_windows[0]  = 5;                    // (hard coded)
    max_windows[1]  = 10;                   // (hard coded)
    omega           = 0.23161204729081489;  // (hard coded) For Si!

    // below assuming insulator
    nv          = _nocc   * _nkpt;
    nc          = _nunocc * _nkpt;

    std::vector<double> _eigs;
    std::vector<double> _eigs_v;
    std::vector<double> _eigs_c;
    for (int is = 0; is < nspin; is++) {
        for (int ik = 0; ik < nkpt; ik++) {
            double* _eig_occ = e_occ[is][ik];
            double* _eig_unocc = e_unocc[is][ik];
            
            int neig = _nocc + _nunocc;
            int nocc = _nocc;
            int nunocc = _nunocc;  // FIXME these should be Psi occ and nunocc
            
            
            _eigs.insert(_eigs.end(),     _eig_occ,       _eig_occ+nocc);
            _eigs.insert(_eigs.end(),     _eig_unocc,       _eig_unocc+nunocc);
            _eigs_v.insert(_eigs_v.end(), _eig_occ,       _eig_occ+nocc);
            _eigs_c.insert(_eigs_c.end(), _eig_unocc,       _eig_unocc+nunocc);
        }
    }
    std::sort(_eigs.begin(), _eigs.end());
    std::sort(_eigs_v.begin(), _eigs_v.end());
    std::sort(_eigs_c.begin(), _eigs_c.end());
    eigs = _eigs;
    eigs_v = _eigs_v;
    eigs_c = _eigs_c;
    
    // Number of valence bands should be read using "loadparameters" first!
    assert(nv != 0);

    evmin       = eigs.front();
    evmax       = eigs[nv-1];
    ecmin       = eigs[nv];
    ecmax       = eigs.back();
    gap         = ecmin - evmax;
}

// void WINDOWING::loadgpp(double*** omsq, int* ng, const SYSINFO& sys) {
//     const int _nspin   = sys.nspin;
//     const int _nq      = sys.nkpt;  // FIXME 
//     std::vector<double> _wppsq;
//     for (int is = 0; is < _nspin; is++) {
//         for (int iq = 0; iq < _nq; iq++) {
//             double* _wpp = omsq[is][iq];
//             int _ng = ng[iq];
//             _wppsq.insert(_wppsq.end(), _wpp, _wpp+_ng);
//         }
//     }
//     std::vector<double> _wpp_pos;
//     for (double i : _wppsq) {
//         if (i > 0) {
//             _wpp_pos.push_back(sqrt(i));
//         }
//     }
//     std::sort(_wpp_pos.begin(), _wpp_pos.end());  // sorted
//     // Use only positive values, negative values are unphysical

//     wpp = _wpp_pos;
//     wppmin = wpp.front();  // first element
//     wppmax = wpp.back();   // last element
// }
/**
 * @brief Load windowing parameters for a matrix element Sigma_nn'
 * 
 * @param w evaluation frequency (E_nk for now, interpolation later)
 * @param omsq GPP mode strength squares
 * @param ng number of G points for each iq
 */
void WINDOWING::sigma_win(const double _w, double* const omsq, int const ng) {
    omega = _w;
    // printf("omega sigmawin %f\n", omega);
    std::vector<double> _wppsq;
    // for (int is = 0; is < nspin; is++) {
    //     for (int iq = 0; iq < nq; iq++) {
            // double* _wpp = omsq[is][iq];
            double* _wpp = omsq;
            // int _ng = ng[iq];
            int _ng = ng;
            _wppsq.insert(_wppsq.end(), _wpp, _wpp+_ng);
    //     }
    // }
    // print_vector(_wppsq);
    // FIXME using only positive ones here!
    std::vector<double> _wpp_pos;
    for (double i : _wppsq) {
        if (i > 0) {
            _wpp_pos.push_back(sqrt(i));
        }
    }
    std::sort(_wpp_pos.begin(), _wpp_pos.end());  // sorted
    // Use only positive values, negative values are unphysical
    // print_vector(_wpp_pos);
    wpp = _wpp_pos;
    wppmin = wpp.front();  // first element
    wppmax = wpp.back();   // last element
}

void WINDOWING::read_from_file() {
  std::string _opt("Sigma");  // Convert char to string
  // printf("reading windows\n");
  if (_opt == "Sigma") {
    bool sigma_windows = true;
    int sigma_window_index = 0;
    FILE *fp = fopen("sigmaHGLinfo.dat", "r");
    
    double omega(0);
    int nwocc(0);
    int nwPPocc(0);
    // read omega
    fscanf(fp, "%lf\n", &omega);
    // printf("omega: %f\n",omega);
    // read data for occupied states
    fscanf(fp, "%d %d\n", &nwocc, &nwPPocc);
    double _occwin_i(0);
    double _PPoccwin_i(0);
    for (int i = 0; i < nwocc+1; i++) {
      fscanf(fp, "%lf", &_occwin_i);
      opt_awins_list[0].push_back(_occwin_i);
    }
    for (int i = 0; i < nwPPocc+1; i++) {
      fscanf(fp, "%lf", &_PPoccwin_i);
      opt_bwins_list[0].push_back(_PPoccwin_i);
    }
    std::vector<int> occ_nnodes;
    int nnodes(0);

    for (int i = 0; i < nwocc; i++) {
      for (int j = 0; j < nwPPocc; j++) {
        fscanf(fp, "%d", &nnodes);
        occ_nnodes.push_back(nnodes);
      }
    }

    std::vector<int> hgl_list;
    int is_hgl(0);

    for (int i = 0; i < nwocc; i++) {
      for (int j = 0; j < nwPPocc; j++) {
        fscanf(fp, "%d", &is_hgl);
        hgl_list.push_back(is_hgl);
      }
    }

    for (int i = 0; i < nwocc; i++) {
      for (int j = 0; j < nwPPocc; j++) {
        double _a_min = opt_awins_list[0][i];
        double _a_max = opt_awins_list[0][i+1];

        double _b_min = opt_bwins_list[0][j];
        double _b_max = opt_bwins_list[0][j+1];

        double _gap = _b_min - _a_max;
        double _bw  = _b_max - _a_min;

        int _ngl = occ_nnodes[i*j + j];  // FIXME: Check if this is actually the case
        double _errfrac = 0;
        bool _is_hgl = hgl_list[i*j +j];
        
        WINPAIR winp(_a_min, _a_max, _b_min, _b_max, _gap, _bw, _ngl, _errfrac, is_hgl, sigma_windows, sigma_window_index);
        winpairs.push_back(winp);
      }
    }
    // Completed reading occupied states
    //
    // Start reading unoccupied states
    sigma_window_index = 1;
    int nwunocc(0);
    int nwPPunocc(0);
    fscanf(fp, "%d %d\n", &nwunocc, &nwPPunocc);
    double _unoccwin_i(0);
    double _PPunoccwin_i(0);
    for (int i = 0; i < nwunocc+1; i++) {
      fscanf(fp, "%lf", &_unoccwin_i);
      opt_awins_list[1].push_back(_unoccwin_i);
    }
    for (int i = 0; i < nwPPunocc+1; i++) {
      fscanf(fp, "%lf", &_PPunoccwin_i);
      opt_bwins_list[1].push_back(_PPunoccwin_i);
    }
    std::vector<int> unocc_nnodes;
    nnodes = 0;

    for (int i = 0; i < nwunocc; i++) {
      for (int j = 0; j < nwPPunocc; j++) {
        fscanf(fp, "%d", &nnodes);
        unocc_nnodes.push_back(nnodes);
      }
    }
    std::vector<int> hgl_list_unocc;
    is_hgl = 0;

    for (int i = 0; i < nwocc; i++) {
      for (int j = 0; j < nwPPocc; j++) {
        fscanf(fp, "%d", &is_hgl);
        hgl_list_unocc.push_back(is_hgl);
      }
    }

    for (int i = 0; i < nwunocc; i++) {
      for (int j = 0; j < nwPPunocc; j++) {
        double _a_min = opt_awins_list[1][i];
        double _a_max = opt_awins_list[1][i+1];

        double _b_min = opt_bwins_list[1][j];
        double _b_max = opt_bwins_list[1][j+1];

        double _gap = _b_min - _a_max;
        double _bw  = _b_max - _a_min;

        int _ngl = unocc_nnodes[i*j + j];
        double _errfrac = 0;
        bool _is_hgl = hgl_list_unocc[i*j +j];
        WINPAIR winp(_a_min, _a_max, _b_min, _b_max, _gap, _bw, _ngl, _errfrac, is_hgl, sigma, sigma_window_index);
        winpairs.push_back(winp);
      }
    }
    double gamma(0);  // FIXME
    // read gamma
    fscanf(fp, "%lg", &gamma);
    fclose(fp);
    // printf("read_windows done\n");
  }  // end if
}
// void WINDOWING::loadparameters(const USRINPUT usrin) {
//     // defaults
//     // (kayahans) keep default for now
//     // double check later if it works for metal
//     double _eV       = 1./27.2114;
//     ptol            = usrin.ptol;
//     errfrac         = ptol/100;             // Percent tolerance converted to errfrac
//     min_gap         = 0.1*_eV;              // (hard coded) Minimum gap to be considered metallic (0.1 eV)
//     max_windows[0]  = 5;                    // (hard coded)
//     max_windows[1]  = 10;                   // (hard coded)
//     omega           = 0.23161204729081489;  // (hard coded) For Si!

//     int _nkpt   = usrin.nkpt;
//     int _nstate = usrin.nstate * _nkpt;
//     nv          = usrin.nocc   * _nkpt;
//     nc          = usrin.nunocc * _nkpt;
//   }

void WINDOWING::printenergies() const {
    for (int i = 0; i < eigs.size(); i++) {
        printf("%d %f\n", i, eigs[i]);
    }
}

void WINDOWING::printparameters() const {
    printf("Windowing method parameters \n");
    printf("\tNumber of valence bands over all k-points: %d \n", nv);
    printf("\tNumber of conduction bands over all k-points: %d \n", nc);
    printf("\tValence band range: from %f to %f Ha\n", evmin, evmax);
    printf("\tConduction band range: from %f to %f Ha\n", ecmin, ecmax);
    printf("\tHOMO - LUMO gap: %f Ha\n", ecmin - evmax);
    printf("\tErrfrac : %f (equal to ptol = %f)\n", errfrac, ptol);
    printf("\tMinimum gap: %f Ha \n", min_gap);
    printf("\tMax valence and conduction windows : (%d, %d) \n", max_windows[0], max_windows[1]);
    // printf("\t w1_d     w1_u     w2_d     w2_u     n ns nb ng\n");
}

void WINDOWING::searchwins(char* _option) {
    std::string _opt(_option);  // Convert char to string
    opt = _opt;
    double _cost_matrix[max_windows[0]][max_windows[1]];  // Store cost_matrix
    int _min_a_windows = 0;
    int _min_b_windows = 0;
    double _max_cost = 1E9;
    double _min_cost = _max_cost;
    std::vector<double> _eawins[2];
    std::vector<double> _ebwins[2];
    for (int Naw : range(1, max_windows[0]+1)) {
        for (int Nbw : range(1, max_windows[1]+1)) {
            if (_opt == "P0") {
                _cost_matrix[Naw-1][Nbw-1] = win_P0_cost(Naw, Nbw, _eawins, _ebwins);
                // print_vector(eawins[0]);
                // print_vector(ebwins[0]);
            } else if (_opt == "Pw") {
                // New, Npw == Naw, Nbw
                _cost_matrix[Naw-1][Nbw-1] = win_Pw_cost(Naw, Nbw, _eawins, _ebwins);
            } else if (_opt == "Sigma") {
                // New, Npw == Naw, Nbw
                if (Naw > 1) { 
                    _cost_matrix[Naw-1][Nbw-1] = win_Sigma_cost(Naw, Nbw, _eawins, _ebwins);
                } else {
                    // When Naw == 1, for sigma case this doesnt make sense 
                    // because valence and conduction bands are all in one window, therefore 
                    // Naw >= 2. Therefore, the cost is set to an arbitrarily high value.  
                    _cost_matrix[Naw-1][Nbw-1] = _max_cost;
                }
                
            } else {
                std::string _err_msg = _opt +" is not in options. Available options are: ";
                for (int i = 0; i < num_opt; i++) {
                    _err_msg += " "+available_opt[i];
                }
                _err_msg += "!\n";
                printf(_err_msg.c_str());
            }
            if (_cost_matrix[Naw-1][Nbw-1] < _min_cost) {
                _min_a_windows   = Naw;
                _min_b_windows   = Nbw;
                _min_cost  = _cost_matrix[Naw-1][Nbw-1];
                opt_awins_list[0] = _eawins[0];
                opt_bwins_list[0] = _ebwins[0];
                opt_awins_list[1] = _eawins[1];
                opt_bwins_list[1] = _ebwins[1];
            }
#ifdef DEBUG
                // printf("Naw=%d, Nbw=%d, cost=%f\n", Naw, Nbw, _cost_matrix[Naw-1][Nbw-1]);
#endif
        }
    }
    // printf("============\n");
    // printf("Min windowing cost for %s calculation Naw=%d, Nbw=%d, cost=%f\n", _option, _min_a_windows, _min_b_windows, _min_cost);
    opt_num_windows[0] = _min_a_windows;
    opt_num_windows[1] = _min_b_windows;
    // add_padding(opt_awins);

    // printf("Window A\n");
    // print_vector(opt_awins_list[0]);
    // if (!opt_awins_list[1].empty()) {
    //     printf("Window A-2\n");
    //     print_vector(opt_awins_list[1]);
    // }
    // printf("Window B\n");
    // print_vector(opt_bwins_list[0]);
    // if (!opt_bwins_list[1].empty()) {
    //     printf("Window B-2\n");
    //     print_vector(opt_bwins_list[1]);
    // }
    // printf("============\n");
    set_winpairs();
}

double WINDOWING::win_P0_cost(int Naw, int Nbw, std::vector<double> (&eawins)[2], std::vector<double> (&ebwins)[2]) {
    // Static polarizability windowing
    // Initial conditions
    double padding = 0.001;
    eawins[0] = linspace(evmin-padding, evmax+padding, Naw+1);
    ebwins[0] = linspace(ecmin-padding, ecmax+padding, Nbw+1);
    // Optimizer
    // Dos is not used to optimize P0, using fixed, flat DOS
    opt_option = "read_dos";
    double cost = fixed_win_optimizer(eawins[0], ebwins[0], eigs, eigs);
    return cost;
}

double WINDOWING::win_Pw_cost(int Naw, int Nbw, std::vector<double> (&eawins)[2], std::vector<double> (&ebwins)[2]) {
    /** Dynamic polarizability windowing
     \sum_{c,v} ( 1/(w-(Ec-Ev)) -1/(w+(Ec-Ev)) )
     Rewrite this as ( 1/(-(Ec-w/2-(Ev+w/2)) -1/(Ec+w/2-(Ev-w/2))) )
     With w > 0 the first term can have 0 in the denominator, but 
     the second term is always positive
    */
    // Initial conditions
    double w_abs = abs(omega);
    double cost = 0.;
    
    std::vector<double> eawins_1 = linspace(evmin+w_abs/2, evmax+w_abs/2, Naw+1);
    std::vector<double> eawins_2 = linspace(evmin-w_abs/2, evmax-w_abs/2, Naw+1);
    std::vector<double> ebwins_1 = linspace(ecmin-w_abs/2, ecmax-w_abs/2, Nbw+1);
    std::vector<double> ebwins_2 = linspace(ecmin+w_abs/2, ecmax+w_abs/2, Nbw+1);

    // Optimizer
    opt_option = "read_dos";
    cost  = fixed_win_optimizer(eawins_1, ebwins_1, eigs, wpp);
    cost += fixed_win_optimizer(eawins_2, ebwins_2, eigs, wpp);
    // Results
    eawins[0] = eawins_1;
    eawins[1] = eawins_2;
    ebwins[0] = ebwins_1;
    ebwins[1] = ebwins_2;
    return cost;
}
/**
 * @brief Calculate windowing cost for a Sigma calculation
 * 
 * Eq. 30 and 31 of PRB 101, 035139 (2020)
 * Since there is Sigma+ and Sigma-, you need two windowings for eq. 30 and 31 separately
 * Here, New = Nvw + Ncw
 * Npw, Nvw and Ncw are as defined in the equations
 * 
 * @param New number of windows over all bands
 * @param Npw number of windows for poles
 * @param eawins window boundaries of valence and conduction bands
 * @param ebwins window boundaries of pole strengths
 * @return double 
 */
double WINDOWING::win_Sigma_cost(int New, int Npw, std::vector<double> (&eawins)[2], std::vector<double> (&ebwins)[2]) {

  // Initial conditions
  double padding = 0.001;
  std::vector<double> a_v_wins, a_c_wins;  // valence and conduction band windows

  std::vector<double> pp_wins = linspace(log(wppmin-padding), log(wppmax+padding), Npw+1);
  std::vector<double> b_pp_v_wins = exp(pp_wins); // +ve pp windows (with valence bands)
  
  // std::vector<double> pp_wins = linspace(wppmin-padding, wppmax+padding, Npw+1);
  // std::vector<double> b_pp_v_wins = pp_wins; // +ve pp windows (with valence bands)
  std::vector<double> b_pp_c_wins = b_pp_v_wins; // -ve pp windows (with conduction bands)
  
  // always sort windows when you initialize
  std::sort(b_pp_v_wins.begin(), b_pp_v_wins.end());
  std::sort(b_pp_c_wins.begin(), b_pp_c_wins.end());
  
  int Ncw;
  double cost = 1E9;
  double new_cost;
  opt_option = "read_dos";
  // printf("omega %f\n", omega);
  for (int Nvw = 1; Nvw < New; Nvw++) {
    new_cost = 0;
    Ncw = New - Nvw;  // Number of conduction band windows are implicitly defined
    // total number of windows is fixed (New), number of windows over
    //
    // -1/((E_v - w) - wpp)
    // 1/(w - E_c - wpp)
    std::vector<double> ea_v = eigs_v - omega;  // E_v - w (for Sigma+)
    std::vector<double> ea_c = omega - eigs_c;  // w - E_c (for Sigma-)

    std::vector<double> eb_pp_v = wpp;          //  wpp (for Sigma+)
    std::vector<double> eb_pp_c = wpp;          //  wpp (for Sigma-)

    // Not sure if this is really needed here, maybe pre-sorted, but doesn't matter much
    std::sort(ea_v.begin(), ea_v.end());
    std::sort(ea_c.begin(), ea_c.end());
    std::sort(eb_pp_v.begin(), eb_pp_v.end());
    std::sort(eb_pp_c.begin(), eb_pp_c.end());

    a_v_wins = linspace(ea_v.front() - padding, ea_v.back() + padding, Nvw+1);
    a_c_wins = linspace(ea_c.front() - padding, ea_c.back() + padding, Ncw+1);

    new_cost += fixed_win_optimizer(a_v_wins, b_pp_v_wins, ea_v, eb_pp_v);
    new_cost += fixed_win_optimizer(a_c_wins, b_pp_c_wins, ea_c, eb_pp_c);

    if (new_cost < cost) {
        cost = new_cost;
        // Results
        eawins[0] = a_v_wins;
        eawins[1] = a_c_wins;
        ebwins[0] = b_pp_v_wins;
        ebwins[1] = b_pp_c_wins;
    }
  }
  return cost;
}

double WINDOWING::fixed_win_optimizer(std::vector<double>& eawins,
                                    std::vector<double>& ebwins,
                                    const std::vector<double>& ea,
                                    const std::vector<double>& eb) {
    int Naw = eawins.size() - 1;
    int Nbw = ebwins.size() - 1;
    double dE_a = eawins[1] - eawins[0];
    double dE_b = ebwins[1] - ebwins[0];  // Since this was initialized in log scale, this is the smallest step size 
    double dEaderiv = dE_a/100;
    double dEbderiv = dE_b/100;
    // printf("dE_a=%f, dE_b=%f\n", dE_a, dE_b);
    // Define cost function
    // Initial cost
    double cost = costfun(eawins, ebwins, ea, eb);
    // Cost of windows (l, m)
    double costl, costm;
    std::vector<double> dcostl(eawins.size(), 0.0);
    std::vector<double> dcostm(ebwins.size(), 0.0);
    
    std::vector<double> step_eawins(eawins.size(), 0.0);
    std::vector<double> step_ebwins(ebwins.size(), 0.0);
    
    static double fracquit = 1E-6;
    double dcost = 100*cost;
    double cost_new = 1E9;
    double conv     = 1E9;
    double gamma = fracquit*100;
    int num_break = 0;
    // If number of windows is larger than 1, then optimize window boundaries
    if (Naw > 1 || Nbw > 1) {
        while (conv > fracquit && num_break < 3) {
            for (int l = 1; l < Naw; l++) {
                eawins[l] += dEaderiv;
                costl      = costfun(eawins, ebwins, ea, eb);
                // FIXME returned cost above shouldnt be integer, then derivative cant be calculated!
                // Store number of quadrature points at the end once more
                eawins[l] -= dEaderiv;
                dcostl[l]  = (costl-cost)/dEaderiv;
                if (dcostl[l] != 0.0) {
                    step_eawins[l]  = -cost/dcostl[l]*gamma;
                }
            }
            for (int m = 1; m < Nbw; m++) {
                // printf("Before\n");
                // print_vector(ebwins);
                ebwins[m] += dEbderiv;
                // printf("After\n");
                // print_vector(ebwins);
                costm      = costfun(eawins, ebwins, ea, eb);
                ebwins[m] -= dEbderiv;
                dcostm[m]  = (costm-cost)/dEbderiv;
                // print_vector(dcostm);
                if (dcostm[m] != 0.0) {
                    step_ebwins[m]  = -cost/dcostm[m]*gamma;
                }
            }
            // print_vector(ebwins);
            // print_vector(step_ebwins);
            eawins   = eawins + step_eawins;
            ebwins   = ebwins + step_ebwins;
            cost_new = costfun(eawins, ebwins, ea, eb);
            dcost    = cost_new - cost;
            conv     = dcost/cost;
            if (conv < 0) {
                // Successful step
                // printf("1 Naw=%d, Nbw=%d, gamma=%f, conv=%f, cost=%f, cost_new=%f\n", Naw, Nbw, gamma, conv, cost, cost_new);
                cost      = cost_new;
                num_break = 0;
                conv      = std::abs(conv);
                gamma    *= 1.1;
            } else {
                // If consecutive three steps fail, just report that cost
                // printf("2  Naw=%d, Nbw=%d, gamma=%f, conv=%f, cost=%f, cost_new=%f\n", Naw, Nbw, gamma, conv, cost, cost_new);
                num_break += 1;
                eawins     = eawins - step_eawins;
                ebwins     = ebwins - step_ebwins;
                gamma     *= 0.1;
                dcost      = 1E9;
                cost_new   = 1E9;
                conv       = dcost/cost_new;
            }
            // printf("=======\n");
        }
    }
    // exit(1);
    return cost;
}

double WINDOWING::costfun(const std::vector<double>& eawins, 
                        const std::vector<double>& ebwins, 
                        const std::vector<double>& ea, 
                        const std::vector<double>& eb) {
    double cost;
    if (opt_option == "constant_dos") {
        // constant dos (flat)
        cost = fixedwin_fixeddos_costfun(eawins, ebwins);
    } else if (opt_option == "read_dos") {
        // read dos
        cost = fixedwin_costfun(eawins, ebwins, ea, eb);
    } else {
        std::string err_msg = opt_option +" is not in dos reading options";
        for (int i = 0; i < num_opt; i++) {
            err_msg += " "+available_opt_options[i];
        }
        err_msg += "!\n";
        printf(err_msg.c_str());
    }
    return cost;
}

double WINDOWING::fixedwin_costfun(const std::vector<double>& eawins,
                                const std::vector<double>& ebwins,
                                const std::vector<double>& ea,
                                const std::vector<double>& eb) {
    int Naw = eawins.size() - 1;
    int Nbw = ebwins.size() - 1;
    double cost      = 0.;  // Total cost
    double N_lm      = 0.;  // Cost to calculate laplace transform
    int l_la         = 0;   // number of states in a
    int l_mb         = 0;   // number of states in b
    double a_bw      = 0.;  // a bandwidth
    double b_bw      = 0.;  // b bandwidth
    double ab_bw     = 0.;  // ab bandwidth
    double ab_gap    = 0.;  // ab gap
    for (int l = 0; l < Naw; l++) {
        a_bw  = abs(eawins[l+1] - eawins[l]);
        l_la  = std::upper_bound(ea.begin(), ea.end(), eawins[l+1])
              - std::lower_bound(ea.begin(), ea.end(), eawins[l]);
        for (int m = 0; m < Nbw; m++) {
            b_bw    = abs(ebwins[m+1] - ebwins[m]);
            l_mb    = std::upper_bound(eb.begin(), eb.end(), ebwins[m+1])
                    - std::lower_bound(eb.begin(), eb.end(), ebwins[m]);
            ab_bw   = ebwins[m+1] - eawins[l];
            ab_gap  = ebwins[m] - eawins[l+1];
            N_lm    = floor(Nlm_cost(a_bw, b_bw, ab_bw, ab_gap)); //FIXME
            cost   += N_lm *(l_la + l_mb);
            // printf("a_bw: %f, a_gap: %f, Naw: %d, Nbw: %d, l: %d, m: %d, N_lm: %f, l_la: %d, l_mb: %d\n", ab_bw, ab_gap, Naw, Nbw, l, m, N_lm, l_la, l_mb);
        }
    }
    return cost;
}

double WINDOWING::fixedwin_fixeddos_costfun(const std::vector<double> eawins, 
                                            const std::vector<double> ebwins) {
    int Naw = eawins.size() - 1;
    int Nbw = ebwins.size() - 1;
    double cost      = 0.;  // Total cost
    double N_lm      = 0.;  // Cost to calculate laplace transform
    double a_bw      = 0.;  // a bandwidth
    double b_bw      = 0.;  // b bandwidth
    double ab_bw     = 0.;  // ab bandwidth
    double ab_gap    = 0.;  // ab gap
    for (int l = 0; l < Naw; l++) {
        a_bw  = eawins[l+1] - eawins[l];
        for (int m = 0; m < Nbw; m++) {
            b_bw    = ebwins[m+1] - ebwins[m];
            ab_bw   = ebwins[m+1] - eawins[l];
            ab_gap  = ebwins[m] - eawins[l+1];
            N_lm    = Nlm_cost(a_bw, b_bw, ab_bw, ab_gap);
            cost   += N_lm * ((eawins[l+1]- eawins[l]) + (ebwins[m+1]- ebwins[m]));
            // printf("ab_bw: %f, ab_gap: %f, Naw: %d, Nbw: %d, l: %d, m: %d, N_lm: %f\n", ab_bw, ab_gap, Naw, Nbw, l, m, N_lm);
            // if (ab_gap < 0) {
            //     print_vector(eawins);
            //     print_vector(ebwins);
            // }
        }
    }
    return cost;
}





double WINDOWING::Nlm_cost(double a_bw, double b_bw, double ab_bw, double ab_gap) {
    double N_lm = 0.;  // Number of quadrature points for window pair (l,m)
    if (ab_bw > a_bw + b_bw) {
        if (ab_gap < min_gap) { // For gapless systems, gap is set to 1/beta via scissors Sec. D.
            double nquad_o = hgl_quad_cost(ab_bw);
            double nquad_n = quad_cost(ab_bw, ab_gap);
            N_lm           = ceil(nquad_n + (nquad_o-nquad_n)*1./(exp((ab_gap-min_gap/2)/(min_gap/10))+1));
        } else { // Non-overlapping window
            N_lm = quad_cost(ab_bw, ab_gap);
        }
    } else {  // Overlapping windows
        N_lm = hgl_quad_cost(ab_bw);
    }
    return N_lm;
}

double WINDOWING::quad_cost(double _ab_bw, double _ab_gap) {
    // Eq. 25 in PHYSICAL REVIEW B 101, 035139 (2020)
    double _alpha = sqrt(_ab_bw/_ab_gap);
    double _nq    = (-0.27*log(errfrac) + 0.4) * _alpha;
    // (kayahans) in the paper -0.49 doesn't exist
    // printf("%f, %f, %f, %f, %f\n", nq, ceil(nq), alpha, ab_bw, ab_gap);
    return _nq;  // (kayahans) Ceil would be fail safe but floor is good enough
}

double WINDOWING::hgl_quad_cost(double ab_bw) {
    // Eq. D3 in PHYSICAL REVIEW B 101, 035139 (2020)
    double c2 = -0.0036*log(errfrac)+0.11;
    double c1 = -0.0043*pow(log(errfrac), 2)-0.13*log(errfrac)+0.54;
    double c0 = -0.204*log(errfrac)-0.29;

    double nq = c2*pow(ab_bw, 2) + c1*ab_bw + c0;
    return floor(nq);  // (kayahans) Ceil would be fail safe but floor is good enough
}

double WINDOWING::get_omega() {
    return omega; 
}

// void WINDOWING::from_input(std::vector<double> wa, std::vector<double> wb) {
//     // Kept for debugging
//     opt_windows[0]= wa.size() - 1;
//     opt_windows[1]= wb.size() - 1;
//     opt_awins[0] = wa;
//     opt_bwins[0] = wb;
//     set_winpairs();
// }

// void WINDOWING::add_padding(std::vector<double> (&win)[2]) {
//     double padding = 0.001;
//     win[0][0] -= padding;
//     win[0].back() += padding;

//     if (!win[1].empty()) {
//         win[1][0] -= padding;
//         win[1].back() += padding;
//     }
// }


// void WINPAIR::set_numnodes() {
//     const GL gl;  // Gauss-Laguerre library, stores all nodes and weights
//     double exact, estimated, errorG, errorW;
//     const int maxnnode = gl.nptmax;
//     const double G = gap;
//     const double W = bw;
//     for (int nnode = 1; nnode < maxnnode+1; nnode++) {
//         func_opta(nnode, a);
//         // 1. Ec-Ev = gap
//         exact = 1./G;
//         estimated = 0;
//         for (int k = 0; k < nnode; k++) {
//             estimated += gl.w[nnode-1][k] * exp(-1.0 * gl.n[nnode-1][k] * (G/a-1));
//         }
//         estimated *= 1.0/a;
//         errorG = abs(exact - estimated) * G;
//         // 2. Ec-Ev = bandwidth
//         exact = 1./bw;
//         estimated = 0;
//         for (int k = 0; k < nnode; k++) {
//             estimated += gl.w[nnode-1][k] * exp(-1.0 * gl.n[nnode-1][k] * (W/a-1));
//         }
//         estimated *= 1.0/a;
//         errorW = abs(exact - estimated) * W;

//         if ( errorG < errfrac &&  errorW < errfrac ) {
//             gl_nnodes = nnode;
//             for (int i = 0; i < nnode; i++) {
//                 gl_nodes.push_back(gl.n[nnode-1][i]);
//                 gl_weights.push_back(gl.w[nnode-1][i]);
//             }
//             // printf(" Kayahan G=%f, W=%f, nnodes %d, opta %f, errorG %f, errorW %f\n", G, W, nnode, a, errorG, errorW);
//             break;
//         } else if (nnode == maxnnode) {
//             printf("Convergence errfrac cannot be reached using Gauss-Laguerre quadrature!");
//         }
//         // else {
//         //     printf(" G=%f, W=%f, nnodes %d, opta %f, errorG %f, errorW %f\n", G, W, nnode, a, errorG, errorW);
//         // }
//     }
// }

// void WINPAIR::set_numnodes2() {
//     const GL gl;  // Gauss-Laguerre library, stores all nodes and weights
//     const double G = gap;
//     const double W = bw;
//     a = sqrt(G*W);
//     int nnode = quad_cost(W, G);
//     double exact, estimated, errorG, errorW;
//     exact = 1./G;
//     estimated = 0;
//     for (int k = 0; k < nnode; k++) {
//         estimated += gl.w[nnode-1][k] * exp(-1.0 * gl.n[nnode-1][k] * (G/a-1));
//     }
//     estimated *= 1.0/a;
//     errorG = abs(exact - estimated) * G * 100;
//     // 2. Ec-Ev = bandwidth
//     exact = 1./bw;
//     estimated = 0;
//     for (int k = 0; k < nnode; k++) {
//         estimated += gl.w[nnode-1][k] * exp(-1.0 * gl.n[nnode-1][k] * (W/a-1));
//     }
//     estimated *= 1.0/a;
//     errorW = abs(exact - estimated) * W * 100;
//     for (int i = 0; i < nnode; i++) {
//         gl_nodes.push_back(gl.n[nnode-1][i]);
//         gl_weights.push_back(gl.w[nnode-1][i]);
//     }
//     printf(" New Fitted    G=%f, W=%f, nnodes %d, opta %f, errorG %f, errorW %f\n", G, W, nnode, a, errorG, errorW);
// }

// double WINPAIR::quad_cost(double ab_bw, double ab_gap) {
//     double alpha = sqrt(ab_bw/ab_gap);
//     double nq    = (-0.2788*log(errfrac) + 0.4065) * alpha- 0.49;  // (kayahans) in the paper -0.49 doesn't exist
//     return ceil(nq);
// }


// double WINDOWING::overlapping_quad_cost(double ab_bw) {
    // Below here is from MATLAB code
    // double x     = ab_bw/delta;
    // int maxiter  = 10000;
    // double q     = errfrac;
    // double n     = 0.1 * pow(x, 2);
    // double q_tmp = tanh(pow(x, 2*n))*exp(-(1+3.3*n)*exp(-0.68*pow(x, 2)/n));
    // for (int i = 0; i < maxiter; i++) {
    //     if (q_tmp > q) {
    //         n = n * 1.1;
    //     } else {
    //         n = n * 0.9;
    //     }
    //     q_tmp = tanh(pow(x, 2*n))*exp(-(1+3.3*n)*exp(-0.68*pow(x, 2)/n));

    //     if (abs(q_tmp - q)/q < 0.01) {
    //         i = maxiter;
    //     }
    // }
    // double q_NR = tanh(pow(x, 2*n))*exp(-(1+3.3*n)*exp(-0.68*pow(x, 2)/n));
    // if (abs(q_NR - q) > 0.1 * q) {
    //     printf("n does not produce proper error");
    // }
    // double nq = n;
    // return nq;
// }

// double WINDOWING::win_GPP_cost(int New, int Nbw) {
//     /** 
//      * Fixed window cost for Sigma
//      * New: number of windows over all bands
//      * Nvw: number of windows over filled states
//      * Ncw: number of windows over empty states
//      * Nabw: Naw + Nbw
//      * Npw: number of windows for poles
//      */
    
//     std::vector<double> vwins, cwins;
//     std::vector<double> pwins = linspace(log(wppmin), log(wppmax), Npw+1);
//     std::vector<double> pvwins = exp(pwins);
//     std::vector<double> pcwins = -exp(pwins);
    
//     double N_lm      = 0.;  // Cost to calculate laplace transform
//     int l_m          = 0;   // number of states in a
//     int l_l          = 0;   // number of states in b
//     double mv_bw      = 0.;  // valence window bandwidth
//     double mc_bw      = 0.;  // conduction window bandwidth
//     double l_bw      = 0.;  // wpp bandwidth
//     double ml_bw     = 0.;  // ab bandwidth
//     double ml_gap    = 0.;  // ab gap
//     int opt_Nvw      = 0;
//     int opt_Ncw      = 0;
//     int opt_Npw      = 0;
//     int Ncw;
//     double cost = 10E6;

//     for (int Nvw = 1; Nvw < New; Nvw++) {
//         Ncw = New - Nvw;
//         // total number of windows is fixed (New), number of windows over
//         // valence band ()
//         vwins = linspace(omega - evmin, omega - evmax, Nvw+1);
//         cwins = linspace(omega - ecmin, omega - ecmax, Ncw+1);
//         for (int mv = 0; mv < Nvw - 1; mv++) {
//             for (int mc = 0; mc < Ncw - 1; mc++) {
//                 for (int l = 0; l < Npw - 1; l++) {
//                     // Filled states
//                     // x = w - E_v + w_p
//                     mv_bw = vwins[mv+1] - vwins[mv];
//                     l_bw = pvwins[l+1] - pvwins[l];
//                     ml_bw   = pvwins[l+1] - vwins[mv];
//                     ml_gap  = pvwins[l] - vwins[mv+1];
                    
//                     N_lm = Nlm_cost(mv_bw, l_bw, ml_bw, ml_gap);
                    
//                     // Empty states
//                     // x = w - E_c - w_p
//                     mc_bw = cwins[mc+1] - cwins[mc];
//                     ml_bw   = pcwins[l+1] - cwins[mc];
//                     ml_gap  = pcwins[l] - cwins[mc+1];
//                     N_lm += Nlm_cost(mc_bw, l_bw, ml_bw, ml_gap);

//                     if (N_lm < cost) {
//                         // cost = Densities and the factor of 2 in front is
//                         // removed as in the cost function in Sec. E. 4.
//                         cost = N_lm;
//                         opt_Nvw = mv;
//                         opt_Ncw = mc;
//                         opt_Npw = l;
//                     }
//                 }
//             }
//         }
//     }
//     return cost;
// }

// double WINDOWING::win_Pw_cost(int Naw, int Nbw) {
//     /** Dynamic polarizability windowing
//      \sum_{c,v} ( 1/(w-(Ec-Ev)) -1/(w+(Ec-Ev)) )
//      Rewrite this as ( 1/(-(Ec-w/2-(Ev+w/2)) -1/(Ec+w/2-(Ev-w/2))) )
//      With w > 0 the first term can have 0 in the denominator, but 
//      the second term is always positive
//     */
//     double w_abs = abs(omega);

//     std::vector<double> eawins_1 = linspace(evmin+w_abs/2, evmax+w_abs/2, Naw+1);
//     std::vector<double> eawins_2 = linspace(evmin-w_abs/2, evmax-w_abs/2, Naw+1);
//     std::vector<double> ebwins_1 = linspace(ecmin-w_abs/2, ecmax-w_abs/2, Nbw+1);
//     std::vector<double> ebwins_2 = linspace(ecmin+w_abs/2, ecmax+w_abs/2, Nbw+1);
    
    // double cost      = 0.;  // Total cost
    // double N_lm      = 0.;  // Cost to calculate laplace transform
    // int l_la         = 0;   // number of states in a
    // int l_mb         = 0;   // number of states in b
    // double a_bw_1    = 0.;  // a bandwidth for the first term
    // double b_bw_1    = 0.;  // b bandwidth for the first term
    // double a_bw_2    = 0.;  // a bandwidth for the second term
    // double b_bw_2    = 0.;  // b bandwidth for the second term
    // double ab_bw_1   = 0.;  // ab bandwidth for the first term
    // double ab_gap_1  = 0.;  // ab gap for the first term
    // double ab_bw_2   = 0.;  // ab bandwidth for the second term
    // double ab_gap_2  = 0.;  // ab gap for the second term

    // for (int l = 0; l < Naw; l++) {
    //     a_bw_1  = eawins_1[l+1] - eawins_1[l];
    //     a_bw_2 = 0;  // No need to calculate bandwidths for the second term
    //     l_la  = std::upper_bound(eigs.begin(), eigs.end(), eawins_1[l+1])
    //           - std::lower_bound(eigs.begin(), eigs.end(), eawins_1[l]);
    //     for (int m = 0; m < Nbw; m++) {
    //         b_bw_1    = ebwins_1[m+1] - ebwins_1[m];
    //         b_bw_2    = 0;  // No need to calculate bandwidths for the second term
    //         l_mb    = std::upper_bound(eigs.begin(), eigs.end(), ebwins_1[m+1])
    //                 - std::lower_bound(eigs.begin(), eigs.end(), ebwins_1[m]);
    //         ab_bw_1   = ebwins_1[m+1] - eawins_1[l];
    //         ab_gap_1  = ebwins_1[m] - eawins_1[l+1];
    //         ab_gap_2  = 1.0;
    //         N_lm    = Nlm_cost(a_bw_1, b_bw_1, ab_bw_1, ab_gap_1);
    //         N_lm   += Nlm_cost(a_bw_2, b_bw_2, ab_bw_2, ab_gap_2);
    //         cost   += N_lm *(l_la + l_mb);
    //         // printf("a_bw: %f, a_gap: %f\n", ab_bw, ab_gap);
    //         // printf("Naw: %d, Nbw: %d, l: %d, m: %d, N_lm: %f, l_la: %d, l_mb: %d\n", Naw, Nbw, l, m, N_lm, l_la, l_mb);
    //     }
    // }
    // return cost;
// }

