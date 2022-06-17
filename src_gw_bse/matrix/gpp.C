#include "standard_include.h"
#include "allclass_gwbse.h"
#include "gpp.h"
#include "controller.h"
// #include "eps_matrix.h"
#include "messages.h"
// #include "pmatrix.h"
// #include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"
#include "../diagonalizer/diagonalizer.h"
#include "ckcomplex.h"
#include "fftw3.h"
// #include "constant.h"

#include <cstring> // for memcpy
using std::memcpy;
#include <cmath>
#include <stdio.h>

#define eps_rows 20
#define eps_cols 20
#define IDX_eps(r,c) ((r)*eps_cols + (c))

//------------------------------
// dot product of two vectors (double)
//------------------------------
double dot_product(double v1[3], double v2[3]){
  double result=0.;
  for(int i=0; i<3; i++){
    result += v1[i]*v2[i];
  }
  return result;
}//-----end function

double calc_vola(double* a1, double* a2, double* a3){
  double a[3][3];
  double m[3][3];
  double vol;

  for (int i=0; i<3; i++){
      m[0][i] = a1[i];
      m[1][i] = a2[i];
      m[2][i] = a3[i];
  }

  /* compute matrix of cofactors */
  a[0][0] =  m[1][1]*m[2][2] - m[1][2]*m[2][1];
  a[1][0] = -m[1][0]*m[2][2] + m[1][2]*m[2][0];
  a[2][0] =  m[1][0]*m[2][1] - m[1][1]*m[2][0];
  a[0][1] = -m[0][1]*m[2][2] + m[0][2]*m[2][1];
  a[1][1] =  m[0][0]*m[2][2] - m[0][2]*m[2][0];
  a[2][1] = -m[0][0]*m[2][1] + m[0][1]*m[2][0];
  a[0][2] =  m[0][1]*m[1][2] - m[0][2]*m[1][1];
  a[1][2] = -m[0][0]*m[1][2] + m[0][2]*m[1][0];
  a[2][2] =  m[0][0]*m[1][1] - m[0][1]*m[1][0];

  vol = m[0][0]*a[0][0] + m[0][1]*a[1][0] + m[0][2]*a[2][0];
  return vol;
}

Gpp::Gpp () {
  // readInputFile();
  // if (1); // TODO(kayahans): correct later (gpp_is_on)
    // readRho(); // read rho data
}

Gpp::Gpp(MatrixConfig config) : CBase_Gpp(config) {
  GWBSE* gwbse = GWBSE::get();
  blockSize = config.tile_rows;
  numBlocks = config.chareRows();

  // Set some constants
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = config.qindex; // The controller sets this

  total_time = 0.0;
  data_received = 0;
  calculate_vc();
  // readInputFile();
  // if (1); // TODO(kayahans): correct later (gpp_is_on)
  //   readRho(); // read rho data
}

// read input files
// void Gpp::readInputFile(){
//   // TODO replace this using the regular input file later
  
//   if (thisIndex.x == 0 and thisIndex.y == 0) {
//     char fileName[1000];
//     int mype = CKMYPE();
//     sprintf(fileName, "gpp.in");
//     FILE *fp;
//     if ((fp = fopen(fileName,"r"))) {
//       fscanf(fp, "%d", &gpp_is_on);
//       fscanf(fp, "%d", &qespresso);
//       fscanf(fp, "%s", &rhoFile);
//       fscanf(fp, "%d", &num_q);
//       fscanf(fp, "%d", &num_w);
//       w = new double [num_w];
//       for(int i=0; i<num_w; i++){
//         fscanf(fp,"%lg",&w[i]);
//       }
//       gpp_is_on = true;
//       CkPrintf("[GPP] GPP is on!\n");
//     } else {
//       gpp_is_on = false;
//       CkPrintf("[GPP] No GPP requested!\n");
//     }
//     fclose(fp);
//   }
// }

void Gpp::setQIndex(int _qindex) {
  qindex = _qindex;
  GWBSE *gwbse = GWBSE::get();
  
  qcryst = gwbse->gwbseopts.qvec[qindex];

  b1 = gwbse->gwbseopts.b1;
  b2 = gwbse->gwbseopts.b2;
  b3 = gwbse->gwbseopts.b3;
  for(int i=0; i<3; i++){
    qcart[i] = b1[i]*qcryst[0] + b2[i]*qcryst[1] + b3[i]*qcryst[2];
    qcart[i] *= 2.0 * M_PI / alat;
  }

  double* a1 = gwbse->gwbseopts.a1;
  double* a2 = gwbse->gwbseopts.a2;
  double* a3 = gwbse->gwbseopts.a3;
  vol = calc_vola(a1, a2, a3);
  alat = 10.261200; // TODO (kayahans) hardcoded
  nkpt = gwbse->gwbseopts.nkpt;

  nfft = gwbse->gw_parallel.fft_nelems;
  ndata = nfft[0]*nfft[1]*nfft[2];
  
}

// void find_rho_gdiff(int gdiff[3], int &gindex, bool &gdiffTrue int ga[3], int gb[3], int ndata){
  
//   // int ndata = nr[0]*nr[1]*nr[2];
//   gdiffTrue = false;

//   for(int i=0; i<ndata; i++){
//     //we need to compare gdiff with density g index
//     if( gdiff[0] == ga[i] && gdiff[1] == gb[i] && gdiff[2] == gc[i] ){
//       gindex = i;
//       gdiffTrue = true;
//       printf("Inside %d %d %d %d\n", gdiff[0], gdiff[1], gdiff[2], gindex);
//       break;
//     }
//   }
// }


void Gpp::recvCopy(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, transferToGpp_complete), controller_proxy));
}

void Gpp::recv_eig(std::vector<double> new_data) {
  int x = thisIndex.x;
  int start_col = x * eps_rows;
  eigval = new double[eps_rows];
  for(int i=0;i<eps_rows;i++) {
    eigval[i] = new_data[i+start_col];
    // printf("eig %d %.8e\n", i+start_col, eigval[i]);
  }
}

void Gpp::calculate_vc() {
  // TODO (kayahans) For now just use 1/q2 regular vc, but change to averages in the future
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  std::vector<double> vcoulbp = psi_cache->getVCoulb();
  ng = vcoulbp.size();
  vcoulb = new double [ng];
  std::copy(vcoulbp.begin(), vcoulbp.end(), vcoulb);
  if(qindex==0)
      vcoulb[0] = psi_cache->getVCoulb0();
  
  // omsq values are stored only on the diagonal tiles
  if (thisIndex.x == thisIndex.y) {
    omsq = new double [eps_cols];
  }
  contribute(CkCallback(CkReductionTarget(Controller, gpp_vc_complete), controller_proxy));
  // Check for garbage collector?
}

void Gpp::fft_RtoG(int qindex) {

  if (thisIndex.x == 0 && thisIndex.y == 0) {
    FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();
    PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
    complex* rhoData = psi_cache->getRhoData();
    int* nr = psi_cache->getRhosize();
    int nrsize = nr[0]*nr[1]*nr[2];
    // for (int i = 0; i<100; i++) {
    //   printf("rho %d %f\n", i, rhoData[i].re);
    // }
    // CkPrintf("Rhodata R %f\n", rhoData[0].re);
    // complex* rhoDataG;
    // int nfft[3];
    // nfft[0] = 10;
    // nfft[1] = 10;
    // nfft[2] = 10;
    // int ndata = nfft[0]*nfft[1]*nfft[2];
    // rhoDataG = new complex[nrsize];
    
    // TODO: kayahans hardcoded, need to check here more carefully!
    // // int ga[nfft[0]];
    // // int gb[nfft[1]];
    // // int gc[nfft[2]];
    // // fftidx_to_gidx(ga,gb,gc,nr);
    int forward = -1;
    fft_controller->setup_fftw_3d(nr,forward);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();
    put_into_fftbox(nr, rhoData, in_pointer);
    fft_controller->do_fftw();
    // delete [] rhoData;
    // rhoData = new complex[ndata];
    fftbox_to_array(nrsize, out_pointer, rhoData, 1.0);
    // 
    
    PsiMessage* msg = new (nrsize) PsiMessage(nrsize, rhoData);
    msg->spin_index = -1;
    msg->k_index = qindex;
    msg->state_index = -1;
    msg->shifted = false;
    psi_cache_proxy.send_rhodata(msg);
  }
  
  // // contribute(CkCallback(CkReductionTarget(Controller, gpp_fft_complete), controller_proxy));
}

void Gpp::ifft_GtoR(int qindex, int real_epsilon_size, std::vector<int> accept) {
  if (config.chareCols() != 1) {
    CkAbort("FFT not supported for 2D decompositions\n");
  }
  // After the GPP eigenvalues are 1D decomposed in G-space
  // Do FFT from G to R space on each eigenvector

  GWBSE *gwbse = GWBSE::get();
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  int ndata = nfft[0]*nfft[1]*nfft[2];
  int backward = 1;
  int numCoeff = real_epsilon_size;

  // First set up the data structures in the FFTController
  fft_controller->setup_fftw_3d(nfft, backward);
  fftw_complex* in_pointer = fft_controller->get_in_pointer();
  fftw_complex* out_pointer = fft_controller->get_out_pointer();

  // we need to setup fftidx
  int *g[3]; // put_into_fftbox routine takes 2D g array, so we need to do this
  std::vector<int> geps_x = f_cache->getGepsXVector();
  std::vector<int> geps_y = f_cache->getGepsYVector();
  std::vector<int> geps_z = f_cache->getGepsZVector();

  int *ga = new int[real_epsilon_size];
  int *gb = new int[real_epsilon_size];
  int *gc = new int[real_epsilon_size];
  for (int i = 0; i < real_epsilon_size; i++) {
    ga[i] = geps_x[i];
    gb[i] = geps_y[i];
    gc[i] = geps_z[i];
  }
  g[0] = ga;
  g[1] = gb;
  g[2] = gc;
  // int index_deb = 1;
  // if (thisIndex.x == index_deb) {
  //   printf("xy %d,%d\n", thisIndex.x, thisIndex.y);
  //   // printf("%d eig %f omsq %f\n", index_deb, eigval[index_deb], omsq[index_deb]);
  //   for (int i=0; i<real_epsilon_size; i++) {
  //     printf("rhop %d %f %f \n", i, data[i].re, data[i].im);
  //     // i = real_epsilon_size - 1;
  //     // printf("rhop %d %f %f \n", i, data[i].re, data[i].re * data[i].re + data[i].im * data[i].im);
  //   }
  //   CkExit();
  // }
  
  int **fftidx;
  fftidx = new int *[numCoeff];
  for(int i=0; i<numCoeff;i++){ fftidx[i] = new int [3]; }

  // this routine changes negative g index to be a positive numbers
  // since it is origianlly written with Fortran, fftidx has fortran counting,
  // i.e., if gidx is (0,0,0), then (1,1,1) in fftidx
  gidx_to_fftidx(numCoeff, g, nfft, fftidx);
  B_r = new complex [ndata];
  
  // int gidx = 0;
  // for (int inr=0; inr < ndata; inr++) {
  //   B_r[inr] = (0, 0);
  //   if (accept[inr]) {
  //     B_r[gidx] = data[gidx];
  //     gidx++;
  //   }
  // }
  // if (thisIndex.x==0) {
  //   for (int i = 0; i < numCoeff; i++) {
  //     // data[i].re = 1.0;
  //     // data[i].im = 0.0;
  //     printf("input %d %f %f\n", i, data[i].re, data[i].im);
  //   } 
  //   // for (int i = 0; i < numCoeff; i++) {
  //   //   printf("idx %d %d %d %d\n", i, fftidx[i][0], fftidx[i][1], fftidx[i][2]);
  //   // } 
  // }  
  // state coefficients are copied to in_pointer
  // put_into_fftbox was originally written for doublePack = 0 (false)
  // double pack is hard coded. FIXME
  bool doublePack = false;
  put_into_fftbox(numCoeff, data, fftidx, nfft, in_pointer, doublePack);
  fft_controller->do_fftw();
  // if (thisIndex.x == 0 && thisIndex.y == 0) {
  //   for (int i=0; i<ndata; i++)
  //     printf("rho %d %.10f \n", i, in_pointer[i][0]);
  //   // i = real_epsilon_size-1;
  //   // printf("rho %d %.10f %.10f \n", i, data[i].re, in_pointer[i][0]);
  // }  

  // CkExit(1);
  // delete[] data;
  // data = new complex[ndata];
  
  double scale = 1.0; // sqrt(1.0/double(ndata)); // IFFT requires normalization?
  fftbox_to_array(ndata, out_pointer, B_r, scale);
  // if (thisIndex.x == 0 && thisIndex.y == 0) {
  //   // for (int i=0; i<ndata; i++)
  //   int i =0;
  //   printf("rho %d %.10f %.10f \n", i, B_r[i].re, out_pointer[i][0]);
  //   i = ndata-1;
  //   printf("rho %d %.10f %.10f \n", i, B_r[i].re, out_pointer[i][0]);
  // }
  // delete[] ga;
  // delete[] gb;
  // delete[] gc;
  // if (thisIndex.x==0) {
  //   for (int i = 0; i < ndata; i++) {
  //     printf("output %d %f %f\n", i, B_r[i].re, B_r[i].im);
  //   } 
  // }  
  for(int i=0; i<numCoeff;i++)
    delete[] fftidx[i];
  
  delete[] fftidx;
  //   if (thisIndex.x==0) {
  //   for (int i = 0; i < ndata; i++) {
  //     // data[i].re = 1.0;
  //     // data[i].im = 0.0;
  //     printf("output %d %f %f\n", i, B_r[i].re, B_r[i].im);
  //   } 
  //   // for (int i = 0; i < numCoeff; i++) {
  //   //   printf("idx %d %d %d %d\n", i, fftidx[i][0], fftidx[i][1], fftidx[i][2]);
  //   // } 
  // }  
  // printf("%d %d %d\n", thisIndex.x, thisIndex.y, ndata);
  contribute(CkCallback(CkReductionTarget(Controller, gpp_ifft_complete), controller_proxy));
}

// void Gpp::ifft_GtoR(int qindex, int real_epsilon_size, std::vector<int> accept) {
//   if (config.chareCols() != 1) {
//     CkAbort("FFT not supported for 2D decompositions\n");
//   }
//   // After the GPP eigenvalues are 1D decomposed in G-space
//   // Do FFT from G to R space on each eigenvector

//   GWBSE *gwbse = GWBSE::get();
//   FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();
//   FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
//   int* nfft;
//   nfft = gwbse->gw_parallel.fft_nelems;
//   int ndata = nfft[0]*nfft[1]*nfft[2];
//   int ng = real_epsilon_size;
//   int forward = -1;
//   B_r = new complex [ndata];
//   int gidx = 0;
//   for (int inr=0; inr < ndata; inr++) {
//     // B_r[inr] = (0, inr*1.0/10);
//     if (accept[inr]) {
//       B_r[gidx] = data[gidx];
//       gidx++;
//     }
//   }

//   // First set up the data structures in the FFTController
//   fft_controller->setup_fftw_3d(nfft, forward);
//   fftw_complex* in_pointer = fft_controller->get_in_pointer();
//   fftw_complex* out_pointer = fft_controller->get_out_pointer();

//   // we need to setup fftidx
//   int *g[3]; // put_into_fftbox routine takes 2D g array, so we need to do this
//   std::vector<int> geps_x = f_cache->getGepsXVector();
//   std::vector<int> geps_y = f_cache->getGepsYVector();
//   std::vector<int> geps_z = f_cache->getGepsZVector();

//   int **ag;
//   ag = new int *[ng];
//   for (int ig=0; ig < ng; ig++)  {
//     ag[ig] = new int[3];
//     ag[ig][0] = geps_x[ig];
//     ag[ig][1] = geps_y[ig];
//     ag[ig][2] = geps_z[ig];
//     printf("%d %d %d %d\n", ig, ag[ig][0], ag[ig][1], ag[ig][2]);
//   }
//   int **fftidx;
//   fftidx = new int *[ng];
//   for (int ig=0; ig < ng; ig++) {
//     fftidx[ig] = new int[3];
//   }
//   gidx_to_fftidx(ng, ag, nfft, fftidx);

//   if (thisIndex.x==0) {
//     for (int i = 0; i < ng; i++) {
//       printf("input %d %f %f\n", i, data[i].re, data[i].im);
//     } 
//   }
//   // printf("acceptsize %d\n", accept.size());
//   bool doublePack = false;
//   put_into_fftbox(ng, B_r, fftidx, nfft, in_pointer, doublePack);
//   fft_controller->do_fftw();
//   // if (thisIndex.x == 0 && thisIndex.y == 0) {
//   //   // for (int i=0; i<ndata; i++)
//   //   int i =0;
//   //   printf("rho1 %d %.10f %.10f \n", i, data[i].re, in_pointer[i][0]);
//   //   i = real_epsilon_size-1;
//   //   printf("rho %d %.10f %.10f \n", i, data[i].re, in_pointer[i][0]);
//   // }  
//   // delete[] data;
//   // data = new complex[ndata];
  
//   double scale = 1.0; // sqrt(1.0/double(ndata)); // IFFT requires normalization?
//   fftbox_to_array(ndata, out_pointer, B_r, scale);
//   // if (thisIndex.x == 0 && thisIndex.y == 0) {
//   //   // for (int i=0; i<ndata; i++)
//   //   int i =0;
//   //   printf("rho %d %.10f %.10f \n", i, B_r[i].re, out_pointer[i][0]);
//   //   i = ndata-1;
//   //   printf("rho %d %.10f %.10f \n", i, B_r[i].re, out_pointer[i][0]);
//   // }
//   // delete[] ga;
//   // delete[] gb;
//   // delete[] gc;
//   // if (thisIndex.x==0) {
//   //   int vsum = 0;
//   //   for (int i = 0; i < ndata; i++) {
//   //     printf("output %d %f %f %d\n", i, B_r[i].re, B_r[i].im, accept[i]);
//   //     vsum += accept[i];
//   //   } 
//   //   printf("vsum %d\n", vsum);
//   // }  
//   for(int i=0; i<ng;i++)
//     delete[] fftidx[i];
  
//   delete[] fftidx;
  
//   contribute(CkCallback(CkReductionTarget(Controller, gpp_ifft_complete), controller_proxy));
// }


void Gpp::fft_skip() {
  contribute(CkCallback(CkReductionTarget(Controller, fft_skipped), controller_proxy));
}


void Gpp::print_col(int num) { 
  int x = thisIndex.x;
  int y = thisIndex.y;
  int global_row, global_col;
  int i = 0;
  for (int r = 0; r < eps_rows; r++) {
    for (int c = 0; c < eps_cols; c++) {
        global_row = x*eps_rows + r;
        global_col = y*eps_cols + c;
        if ( global_col == num ) {
          double val = data[c*eps_rows + r].re;
          // kayahan debug
          // CkPrintf("[GPP] x/y %d %d global_r/c %d %d val %.8e\n", x, y, global_row, global_col, val);
        }
      i++;
    }
  }
}

void Gpp::calc_M0() {
  // M0 = (q+G)*(q+G')/(|q+G|^2*|q+G'|^2)*rho(G-G')/rho(0)
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();
  std::vector<double> vcoulb = psi_cache->getVCoulb();
  std::vector<double> ga = psi_cache->get_ga();
  std::vector<double> gb = psi_cache->get_gb();
  std::vector<double> gc = psi_cache->get_gc();
  ng = psi_cache->get_ng();
  
  if(qindex==0)
      vcoulb[0] = psi_cache->getVCoulb0(); 
  // if (thisIndex.x==0 and thisIndex.y==0) {
  //   for (int i=0; i < ng; i++) {
  //     printf("%d %f\n", i, vcoulb[i]);
  //   }
  // }
  // CkExit(0);
  complex* rhoData = psi_cache->getRhoData();
  int* nr = psi_cache->getRhosize();
  int nrsize = nr[0] * nr[1] * nr[2];
  int *rhoga, *rhogb, *rhogc;
  rhoga = new int[nrsize];
  rhogb = new int[nrsize];
  rhogc = new int[nrsize];
  fftidx_to_gidx(rhoga, rhogb, rhogc, nr);
  double Wpl2 = double(4) * M_PI * rhoData[0].re/vol;  
  factor = -Wpl2 * 4.0 * M_PI / (vol * nkpt);
  int ig1_start = thisIndex.x * eps_rows;
  int ig2_start = thisIndex.y * eps_cols;
  int ig1_end = (ig1_start + eps_rows > ng) ? ng : ig1_start + eps_rows;
  int ig2_end = (ig2_start + eps_cols > ng) ? ng : ig2_start + eps_cols;
  // CkPrintf("pe %d tile %d %d %d %d ng %d %f %f %f %d\n", CKMYPE(), ig1_start, ig2_start, ig1_end, ig2_end, ng, ga[0],gb[0], gc[0], ga.size());
  // for (int ig1 = ig1_start; ig1 < ig1_end; ig1++) {
  //   for (int ig2 = ig2_start; ig2 < ig2_end; ig2++) {
  int ig1_glob, ig2_glob;
  for (int ig1 = 0; ig1 < ig1_end-ig1_start; ig1++) {
    ig1_glob = ig1 + ig1_start;
    for (int ig2 = 0; ig2 < ig2_end-ig2_start; ig2++) {
      // g vector at ig1, ig2
      ig2_glob = ig2 + ig2_start;
      double g1cryst[3], g2cryst[3];
      
      // g vector in cartesian coordinates at ig1, ig2
      // vc routine adds q to each G when vc routine is called
      g1cryst[0] = double(ga[ig1_glob]);
      g1cryst[1] = double(gb[ig1_glob]);
      g1cryst[2] = double(gc[ig1_glob]);
      g2cryst[0] = double(ga[ig2_glob]);
      g2cryst[1] = double(gb[ig2_glob]);
      g2cryst[2] = double(gc[ig2_glob]);
      // CkPrintf("G1 %d %d %d G2 %d %d %d \n", ga[ig1], gb[ig1], gc[ig1], ga[ig2], gb[ig2], gc[ig2]);
      // change g1 and g2 to cartesian coordiates
      double g1[3], g2[3];
      for (int i = 0; i < 3; i++) {
        g1[i] = b1[i]*g1cryst[0] + b2[i]*g1cryst[1] + b3[i]*g1cryst[2];
        g2[i] = b1[i]*g2cryst[0] + b2[i]*g2cryst[1] + b3[i]*g2cryst[2];
        g1[i] *= 2*M_PI/alat;
        g2[i] *= 2*M_PI/alat;
      }
      int gdiff[3];
      gdiff[0] = int(ga[ig1_glob] - ga[ig2_glob]);
      gdiff[1] = int(gb[ig1_glob] - gb[ig2_glob]);
      gdiff[2] = int(gc[ig1_glob] - gc[ig2_glob]);
      int gdiffIndex = -1;
      bool gdiffTrue = false;
      // find the gdiffIndex
      for(int idx=0; idx<nrsize; idx++){
        //we need to compare gdiff with density g index
        // printf("gdiff %d %d %d gabc %d %d %d \n", gdiff[0], gdiff[1], gdiff[2], ga[idx], gb[idx], gc[idx]);
        if( gdiff[0] == rhoga[idx] && gdiff[1] == rhogb[idx] && gdiff[2] == rhogc[idx] ){
        // if( !gdiffTrue && gdiff[0] == geps->ig[idx] && gdiff[1] == geps->jg[idx] && gdiff[2] == geps->kg[idx] ){  
          gdiffIndex = idx;
          gdiffTrue = true;
          break;
        }
      }
      // calculate M_GG' matrix element
      // CkPrintf("gdiff %d %d\n", gdiffIndex, ndata);
      if( gdiffTrue ) {
        // if (thisIndex.x==0 and thisIndex.y==0) {
        // CkPrintf("Returned %d %d %d %d %d %d \n", ig1, ig2, gdiff[0], gdiff[1], gdiff[2], gdiffIndex);
        // }
        double dp12 = dot_product(g1,g2);
        double dp11 = dot_product(g1,g1);
        double dp22 = dot_product(g2,g2);
        double vcqg1 = vcoulb[ig1_glob]*nkpt*vol/(4.0*M_PI); 
        double vcqg2 = vcoulb[ig2_glob]*nkpt*vol/(4.0*M_PI);
        
        complex rdg0;
        std::complex<double> rdg(rhoData[gdiffIndex].re, rhoData[gdiffIndex].im);
        std::complex<double> rd0(rhoData[0].re, rhoData[0].im);
        std::complex<double> rdgc0;
        rdgc0 = rdg/rd0; 
        rdg0.re = rdgc0.real();
        rdg0.im = rdgc0.imag();
        if (dp11 < 1E-12 && dp22 < 1E-12){
          data[IDX_eps(ig1, ig2)] = vcqg1 * rdg0;
        } else if (dp11 < 1E-12 && dp22 > 1E-12) {
          data[IDX_eps(ig1, ig2)] = 0.0;
        } else if (dp11 > 1E-12 && dp22 < 1E-12) {
          data[IDX_eps(ig1, ig2)] = 0.0;
        } else {
          data[IDX_eps(ig1, ig2)] = sqrt(vcqg1*vcqg2) * dp12/(dp11*dp22) * rdg0;
        }  //end if
        // if (ig1_glob == ig2_glob) {
        //   // CkPrintf("g %d %d dp11 %.12f dp12 %.12f dp22 %.12f vg1 %.12f vg2 %.12f rdg0 %.12f mgpp %.12f\n", ig1_glob, ig2_glob, dp11, dp12, dp22, vcqg1, vcqg2, rdg0.re, data[IDX_eps(ig1, ig2)].re);
        //   CkPrintf("%d %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", ig1_glob, ig2_glob, dp11, dp12, dp22, vcqg1, vcqg2, rdg0.re, data[IDX_eps(ig1, ig2)].re, data[IDX_eps(ig1, ig2)].im);
        //   FFLUSH(stdout);
        // //   // CkPrintf("%f %f %f %f %f %f %d %d %d \n", g1[0], g1[1], g1[2],g2[0], g2[1], g2[2], gdiff[0], gdiff[1], gdiff[2]);
        // }
      } else{
          // Mggp[ig1][ig2] = 0.;
          data[IDX_eps(ig1, ig2)] = 0.;
      } //end if
    } //end ig2
  } //end ig1
  // if (thisIndex.x==0 and thisIndex.y==0) {
  //   for (int i=0; i<10;i++) {
  //     for (int j=0; j<10;j++) {
  //       CkPrintf("M0 %d %d %f\n",i, j, data[IDX_eps(i, j)].re);
  //     }
  //   }
  // }
  contribute(CkCallback(CkReductionTarget(Controller, gpp_M0_complete), controller_proxy));
}


void Gpp::calc_omsq() {
  // omsq values are stored only on the diagonal tile elements
  if (thisIndex.x == thisIndex.y) {
    int start_index = thisIndex.x * eps_cols;
    int end_index = (thisIndex.x + 1) * eps_cols;
    end_index = ( end_index < ng) ? end_index : ng;
    for (int i = 0; i < end_index-start_index; i++ ) {
      omsq[i] = data[IDX_eps(i, i)].re * factor / eigval[i];
      // if (i == 0) {
        // printf("i %d omsq %f %.8e factor %f data %f\n", i+start_index , omsq[i], eigval[i], factor, data[IDX_eps(i, i)].re);
        // fflush(stdout);
      // }
    }
  }

  contribute(CkCallback(CkReductionTarget(Controller, gpp_omsq_complete), controller_proxy));
}

PsiMessage* Gpp::send_data(PsiMessage* msg) {
  // CkPrintf("GPP transfer 1 pe %d tile %d %d %d %d ng %d\n", CKMYPE(), ig1_start, ig2_start, ig1_end, ig2_end, ng);
  int idx = 0;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
      msg->psi[idx] = data[IDX_eps(i,j)];
      idx++;
    }
  }
  // msg->spin_index = ig1_start;
  // msg->state_index = ig2_start;
  msg->k_index = CKMYPE();
  msg->shifted = false;
  return msg;
}

void Gpp::print(int qindex, int fnum) {
  char fname[100];
  int x = thisIndex.x;
  int y = thisIndex.y;
  sprintf(fname,"./debug/GPP_%d_q%d_x%d_y%d.dat", fnum, qindex, x, y); 
  int global_row, global_col;
  int i = 0;
  FILE* fp = fopen(fname, "w");
  for (int r = 0; r < config.tile_rows; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
        global_row = x*config.tile_rows + r;
        global_col = y*config.tile_cols + c;
        double re = data[IDX_eps(r,c)].re;
        double im = data[IDX_eps(r,c)].im;
        fprintf(fp, "%d %d %d %d %.8e %.8e\n", x, y, global_row, global_col, re, im);
      i++;
    }
  }
  fclose(fp);
}

void Gpp::debug() {
  int mype = CkMyPe();
  int r = config.tile_rows;
  int c = config.tile_cols;
  if (thisIndex.x == thisIndex.y) {
    int start_index = 0;
    int end_index = 10;
    for (int i = start_index; i < end_index; i++ ) {
      for (int j = start_index; j < end_index; j++ ) {
        printf("GPP Print %f\n", data[IDX_eps(i, j)].re);
      }
    }
  }
  CkPrintf("[GPP] pe %d r %d c %d\n", mype, r, c);
  // contribute(CkCallback(CkReductionTarget(Controller, gpp_debug_complete), controller_proxy));
}

void Gpp::sendToCacheV(int size) {
  if (config.chareCols() != 1) {
    CkAbort("Eigenvectors are stored in 1D chares\n");
  }  
  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  int r_size = nfft[0]*nfft[1]*nfft[2];
  if (thisIndex.x < size) {
    // complex* new_data;
    // new_data = new complex[r_size];
    int row_idx = thisIndex.x;
    // for(int j=0;j<r_size;j++) {
    //   new_data[j] = data[j];
    //   // printf("col_idx %d idx %d data %f\n", row_idx, j, data[j].re);
    // }
    
    GppVMessage* msg;
    msg = new (r_size) GppVMessage(r_size, B_r);
    msg->spin_index = 0; // TODO 
    msg->q_index = qindex;
    msg->alpha_idx = row_idx;
    msg->size = r_size;
    msg->tot_alpha = size;
    psi_cache_proxy.receiveGppV(msg);
    delete[] B_r;
  }
}

void Gpp::sendToCacheE(int total_size) {
  int start_index = thisIndex.x * eps_rows;
  int end_index = (thisIndex.x + 1) * eps_rows;
  end_index = ( end_index < total_size) ? end_index : ng;
  int size = end_index - start_index;
  // printf("GPP2d %d %d %d %u\n", thisIndex.x, thisIndex.y, config.tile_rows, size);
  // TODO keep it like this for now since GPP array size should always be larger than num nodes
  if (thisIndex.y == 0) {
    GppEMessage* msge;
    // printf("Gpp xy %d %d %f %f %f %f\n", thisIndex.x, thisIndex.y, eigval[0], eigval[end_index-1], omsq[0], omsq[end_index-1]);
    msge = new (size) GppEMessage(size, eigval);
    msge->spin_index = 0; // TODO 
    msge->q_index = qindex;
    msge->start_idx = start_index;
    msge->tot_sent = int(total_size/eps_rows) + (total_size % eps_rows != 0);
    psi_cache_proxy.receiveGppE(msge);
  }
}

void Gpp::sendToCacheO(int total_size) {
  int start_index = thisIndex.x * eps_rows;
  int end_index = (thisIndex.x + 1) * eps_rows;
  end_index = ( end_index < total_size) ? end_index : ng;
  int size = end_index - start_index;
  // printf("GPP2d %d %d %d %u\n", thisIndex.x, thisIndex.y, config.tile_rows, size);
  // TODO keep it like this for now since GPP array size should always be larger than num nodes
  if (thisIndex.x == thisIndex.y) {
    GppEMessage* msgo;
    msgo = new (size) GppEMessage(size, omsq);
    msgo->spin_index = 0; // TODO 
    msgo->q_index = qindex;
    msgo->start_idx = start_index;
    msgo->tot_sent = int(total_size/eps_rows) + (total_size % eps_rows != 0);
    psi_cache_proxy.receiveGppO(msgo);
  }
}


#include "gpp.def.h"
