#include "standard_include.h"
#include "allclass_gwbse.h"
#include "gpp.h"
#include "controller.h"
// #include "eps_matrix.h"
#include "messages.h"
// #include "pmatrix.h"
// #include "states.h"
// #include "fft_routines.h"
// #include "fft_controller.h"
#include "../diagonalizer/diagonalizer.h"
// #include "ckcomplex.h"
// #include "constant.h"

#include <cstring> // for memcpy
using std::memcpy;
#include <cmath>
#include <stdio.h>

extern diagData_t* diagData;

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
  readInputFile();
  if (gpp_is_on)
    readRho(); // read rho data
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
  readInputFile();
  if (gpp_is_on)
    readRho(); // read rho data
}

// read input files
void Gpp::readInputFile(){
  // TODO replace this using the regular input file later
  char fileName[1000];
  int mype = CKMYPE();
  
  sprintf(fileName, "gpp.in");
  FILE *fp;
  if ((fp = fopen(fileName,"r"))) {
    fscanf(fp, "%d", &gpp_is_on);
    fscanf(fp, "%d", &qespresso);
    fscanf(fp, "%s", &rhoFile);
    fscanf(fp, "%d", &num_q);
    fscanf(fp, "%d", &num_w);
    w = new double [num_w];
    for(int i=0; i<num_w; i++){
      fscanf(fp,"%lg",&w[i]);
    }
    gpp_is_on = true;
    if (thisIndex.x == 0 and thisIndex.y == 0)
      CkPrintf("[GPP] GPP is on!\n");
  } else {
    gpp_is_on = false;
    if (thisIndex.x == 0 and thisIndex.y == 0)
      CkPrintf("[GPP] No GPP requested!\n");
  }
  fclose(fp);
  
}

void Gpp::readRho(){
  FILE *fp;
  if ((fp = fopen(rhoFile,"r"))) { 
    // read number of grids 
    fscanf(fp, "%d %d %d", &nr[0], &nr[1], &nr[2]);
    int ndata = nr[0]*nr[1]*nr[2];

    // malloc
    ga = new int [ndata];
    gb = new int [ndata];
    gc = new int [ndata];
    rhoData = new std::complex<double> [ndata];

    // scale rhoData if qespresso is true
    double scale;
    if(qespresso)
      scale = vol/double(ndata); 
    else
      scale = 1.0; 

    // read data
    int counter=0;
    for(int i=0; i<nr[0]; i++){
        for(int j=0; j<nr[1]; j++){
        for(int k=0; k<nr[2]; k++){
            double rho_tmp;
            fscanf(fp, "%lg", &(rho_tmp));
            if( qespresso ){ rho_tmp *= scale; }
            rhoData[counter] = rho_tmp;
            counter += 1;
        }
        }
    }
    fclose(fp);
  } else {
    CkPrintf("rho.dat file couldn't be read, exiting!\n");
    CkExit();
  }
  Wpl2 = double(4) * M_PI * rhoData[0].real()/vol;
}

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

void Gpp::find_rho_gdiff(int gdiff[3], int &gindex, bool &gdiffTrue){

  int ndata = nr[0]*nr[1]*nr[2];
  gdiffTrue = false;
  
  for(int i=0; i<ndata; i++){
    //we need to compare gdiff with density g index
    if( gdiff[0] == ga[i] && gdiff[1] == gb[i] && gdiff[2] == gc[i] ){
      gindex = i;
      gdiffTrue = true;
      break;
    }
  }
}


void Gpp::recvCopy(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, eps_to_gpp_complete), controller_proxy));
}



void Gpp::calculate_vc() {
    // TODO (kayahans) For now just use 1/q2 regular vc, but change to averages in the future
    FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
    PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
    std::vector<double> vcoulb = psi_cache->getVCoulb();
    if(qindex==0)
        vcoulb[0] = psi_cache->getVCoulb0();
    ng = vcoulb.size();
}

// void Gpp::calc_Omsq(){
//   PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
//   std::vector<double> coulb = psi_cache->getVCoulb();
//   double Wpl2 = double(4) * M_PI * opt->rhoData[0].real()/opt->sys->vol;
//   if(qindex==0)
//     coulb[0] = psi_cache->getVCoulb0();

//   for(int i=0;i<config.tile_rows;i++){
//     for(int j=0;j<config.tile_cols;j++){
//       int g = thisIndex.x*config.tile_rows+i;
//       int gp = thisIndex.y*config.tile_cols+j;
//       if(g==gp && g<coulb.size())
//         data[IDX_eps(i,j)] -= 1.0;
//       if(g<coulb.size() && gp<coulb.size())
//         data[IDX_eps(i,j)] *= sqrt(coulb[g])*sqrt(coulb[gp]);
//     }
//   }

//   contribute(CkCallback(CkReductionTarget(Controller, s_ready), controller_proxy));
// }

void Gpp::calc_Omsq() {
  double factor = -Wpl2 * 4.0 * M_PI / (vol * nkpt);
  

  std::complex<double>** Mggp;
  Mggp = new std::complex<double>*[ng];
  for (int iq = 0; iq < ng; iq++) {
    Mggp[iq] = new std::complex<double>[ng];
  }
  // say Mggp = data
  // save (q+G)*(q+G')/(|q+G|^2*|q+G'|^2)*rho(G-G')/rho(0)
  int ig1_start = thisIndex.x * eps_rows;
  int ig2_start = thisIndex.y * eps_cols;
  for (int ig1 = ig1_start; ig1 < eps_rows; ig1++) {
    for (int ig2 = ig2_start; ig2 < eps_cols; ig2++) {
      // g vector at ig1, ig2
      double g1cryst[3], g2cryst[3];

      // g vector in cartesian coordinates at ig1, ig2
      // Note(10/22/18): vcoulb.C changed
      // vc routine adds q to each G when vc routine is called
      double g1[3], g2[3];
      g1cryst[0] = double(ga[ig1]);
      g1cryst[1] = double(gb[ig1]);
      g1cryst[2] = double(gc[ig1]);
      g2cryst[0] = double(ga[ig2]);
      g2cryst[1] = double(gb[ig2]);
      g2cryst[2] = double(gc[ig2]);

      // change g1 and g2 to cartesian coordiates
      for (int i = 0; i < 3; i++) {
        g1[i] = b1[i]*g1cryst[0] + b2[i]*g1cryst[1] + b3[i]*g1cryst[2];
        g2[i] = b1[i]*g2cryst[0] + b2[i]*g2cryst[1] + b3[i]*g2cryst[2];
        g1[i] *= 2*M_PI/alat;
        g2[i] *= 2*M_PI/alat;
      }
      int gdiff[3];
      gdiff[0] = int( ga[ig1] - ga[ig2] );
      gdiff[1] = int( gb[ig1] - gb[ig2] );
      gdiff[2] = int( gc[ig1] - gc[ig2] );
      int gdiffIndex;
      bool gdiffTrue;
      // find the gdiffIndex
      find_rho_gdiff(gdiff, gdiffIndex, gdiffTrue);

      // calculate M_GG' matrix element
      if( gdiffTrue ){
        double dp12 = dot_product(g1,g2);
        double dp11 = dot_product(g1,g1);
        double dp22 = dot_product(g2,g2);
        double vcqg1 = vcoulb[ig1]*nkpt*vol/(4.0*M_PI); 
        double vcqg2 = vcoulb[ig2]*nkpt*vol/(4.0*M_PI); 
        if (dp11 < 1E-12 && dp22 < 1E-12){
          // Mggp[ig1][ig2] = vcqg1 * rhoData[gdiffIndex]/ rhoData[0];
          data[IDX_eps[ig1, ig2]] = vcqg1 * rhoData[gdiffIndex]/ rhoData[0];
        } else if (dp11 < 1E-12 && dp22 > 1E-12) {
          // Mggp[ig1][ig2] = 0.0;
          data[IDX_eps[ig1, ig2]] = 0.0;
        } else if (dp11 > 1E-12 && dp22 < 1E-12) {          
          // Mggp[ig1][ig2] = 0.0; 
          data[IDX_eps[ig1, ig2]] = 0.0;
        } else {
          // Mggp[ig1][ig2] = dp12*(vcqg1*vcqg2) * rhoData[gdiffIndex]/ rhoData[0];
          data[IDX_eps[ig1, ig2]] = dp12*(vcqg1*vcqg2) * rhoData[gdiffIndex]/ rhoData[0];
        }  
      } else{
          // Mggp[ig1][ig2] = 0.;
          data[IDX_eps[ig1, ig2]] = 0.;
      } //end if
    } //end ig2
  } //end ig1
  // do matrix multiplication
  // omega^2 = V^-1 * Mggp * V
  // eigenvectors are stored in S (columns are eigenvector)

  // call LAPACK

  char transformT = 'T'; // transpose
  char transform = 'N';  // do nothing
  double Lalpha = double(1.0);  // scale A*B by this factor
  double Lbeta = double(0.0); // scale initial value of C by this factor
  std::complex<double>** M1;
  std::complex<double>** S;

  

  M1 = new std::complex<double>*[ng];
  S = new std::complex<double>*[ng];
  for (int iq = 0; iq < ng; iq++) {
    M1[iq] = new std::complex<double>[ng];
    S[iq] = new std::complex<double>[ng];
  }
  
  

  // 1. M1 = Mggp * V
  // S[is][iq]->transpose();
  // myGEMM( &transformT, &transform, &ng, &ng, &ng, &Lalpha, Mggp.m, &ng, S->m, &ng, &Lbeta, M1.m, &ng);

  // 2. Mggp = V^-1 * M1
  // S[is][iq]->ctranspose();
  // myGEMM( &transform, &transform, &ng, &ng, &ng, &Lalpha, S->m, &ng, M1.m, &ng, &Lbeta, Mggp.m, &ng);

  // std::ofstream f;
  // std::ostringstream oss;
  // oss << "debug/omsq-iq-" << iq << ".dat";
  // std::string fname_omsq = oss.str();
  // f.open(fname_omsq);
  for (int i = 0; i < ng; i++) {
    omsq[i] = Mggp[i][i].real() * factor / eigval[i];

    // below block is for debugging
    // print omega^2
    // if (i == 0) {
    //   // printf("Printing Omsq and eigval values iq=%d\n", iq);
    //   f << "Printing Omsq and eigval values iq=" << iq << std::endl;
    // }
    // printf("%lg %g\n", omsq[is][iq][i], eigval[is][iq][i]);
    // f << omsq[i] << eigval[i] << std::endl;
  }
  // f.close();
  // FIXME Revert S matrix to normal
  // S[is][iq]->ctranspose();
  // S[is][iq]->transpose();
  // char fname[100];
  // sprintf(fname, "debug/S_s%d_q%d.dat", is, iq);
  // S[is][iq]->printallG(fname,geps);
  // sprintf(fname, "debug/Mggp_s%d_q%d.dat", is, iq);
  // Mggp.printallG(fname, geps);
}

DiagMessage* Gpp::sendDataSimple(DiagMessage* msg) {
  
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

  
  for (int r = 0; r < config.tile_rows ; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
        if (r < msg_rows && c < msg_cols) {
          // data[c*config.tile_rows + r] = msg->data[i];
          data[r*config.tile_cols + c] = msg->data[c*msg_rows + r];
          i++;
        }
        else {
          data[r*config.tile_cols + c].re = 0.0;
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

DiagMessage* Gpp::receiveDataSimple(DiagMessage* msg) {
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
  // TODO (kayahans): read from input
  proc_rows = 2; // gwbse->gw_parallel.proc_rows;
  proc_cols = 2; // gwbse->gw_parallel.proc_cols;

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
  diagData->eig_e = new std::complex<double>[eps_size];
  diagData->eig_v = new std::complex<double>[eps_size*eps_size];  
  CkPrintf("[DIAGONALIZER] Created a pointer with totalsize %d numblocks %d for S matrix global dim %d at pe %d for qindex %d\n", totaldata, numBlocks, eps_size, CkMyPe(), qindex);
  contribute(CkCallback(CkReductionTarget(Controller, diag_setup), controller_proxy));
}



#include "gpp.def.h"
