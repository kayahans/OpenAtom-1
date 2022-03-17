#include "standard_include.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "eps_matrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"
#include "../utils/sigma_util.h"

#define CHARE_NUM 10
#define IDX(r,c) ((r)*config.tile_cols + (c))
#define IDX_eps(r,c) ((r)*eps_cols + (c))
#define eps_chares_x 10
#define eps_chares_y 10
#define eps_rows 20
#define eps_cols 20

PMatrix::PMatrix(MatrixConfig config) : CBase_PMatrix(config) {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = qindex = config.qindex; // The controller will sets this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?
  fft_controller = fft_controller_proxy.ckLocalBranch();

  max_iterations = gwbse->gw_parallel.transpose_stages;
  iteration = 0;
  arrival_counter = 0;
  receive_counter = 0;
  completed_chunks = 0;
  total_time = 0.0;
}

void PMatrix::computeP() {
//if(!(thisIndex.x==0 && thisIndex.y==0)) return;
  GWBSE *gwbse = GWBSE::get();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  LAPLACE *lp = psi_cache->getLP();
  int counter = 0;
  int nocc = gwbse->gw_parallel.L;
  int nwocc = lp->nwocc;
  int p_metal = 0;

  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;
  double*** occ_occ = gwbse->gw_epsilon.Occ_occ;
  double*** occ_unocc = gwbse->gw_epsilon.Occ_unocc;
  int M = gwbse->gw_parallel.M;

  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];
  ndata = config.tile_rows;//ndata/config.tile_rows;
  P_m = new complex[ndata*ndata];

for(int ik=0;ik<gwbse->gw_parallel.K;ik++) {

  int indexes[3]; // This is for metal calculation
  indexes[0] = 0;//is;
  unsigned ikq;
  int umklapp[3];
  kqIndex(ik, ikq, umklapp);
  indexes[1] = ik;
  indexes[2] = ikq;

  // loop over window : occupied states
  for (int iwocc=0; iwocc<nwocc; iwocc++){

    // Occupied state energy in this window
    // Find number of occupied states in this window  
    int this_nocc=0;
#if 0
    if (p_metal){
      int kindex = sys.nkpt * indexes[0] + indexes[2];
      nocc = L.nvk[kindex];
      printf("\nL.nvk[%d] = %d\n", kindex, L.nvk[kindex]);
    }
#endif

    for (int iv=0; iv<nocc; iv++){
      if (lp->bin_occ[iwocc] <= e_occ[0][ikq][iv]  &&  e_occ[0][ikq][iv] < lp->bin_occ[iwocc+1])
        this_nocc += 1;
    } // done calculating occupied states in this window
    // now we save eigenvalues in this window to Eocc

    double Eocc[this_nocc];
    double Occp_occ[this_nocc]; // occupation
    int Eoccidx[this_nocc];
    counter = 0;
    for (int iv=0; iv<nocc; iv++){
      if (lp->bin_occ[iwocc] <= e_occ[0][ikq][iv]  &&  e_occ[0][ikq][iv] < lp->bin_occ[iwocc+1])
        {
 //         if(counter >= this_nocc) CkPrintf("\nNode-%d Error", CkMyNode());
#if 1
          Eocc[counter] = e_occ[0][ikq][iv];
          Occp_occ[counter] = psi_cache->get_OccOcc(ikq, iv);//gwbse->gw_epsilon.Occ_occ[0][ikq][iv];//occ_occ[0][ikq][iv];
#endif
#if 1
        // Reset occupation for metal to prevent some crazy occupation numbers (like 1E-123)
        if(p_metal){
          if ( abs ( Occp_occ[counter]-1e-0 ) < 1e-6){
            Occp_occ[counter] = 1.0;
          }
          else if( abs ( Occp_occ[counter] ) < 1e-6 ){
            Occp_occ[counter] = 0.0;
          }
        }
#endif
        Eoccidx[counter] = iv;
        counter += 1;
      }
    }

    // loop over window : unoccupied states
    int nwunocc = lp->nwunocc;
    for (int iwunocc=0; iwunocc<nwunocc; iwunocc++){
      int nunocc = gwbse->gw_parallel.M;
      int nocc = gwbse->gw_parallel.L;

      // Unoccupied state energy in this window
      // Find number of unoccupied states in this window
      int this_nunocc=0;
      // if metal, reset nunocc (and also need to re-reset nocc for this k
      if (p_metal){
        int kindex = gwbse->gw_parallel.K * indexes[0] + indexes[1];
        nunocc = lp->nck[kindex];
        nocc = lp->nvk[kindex];
      }

      for (int ic=0; ic<nunocc; ic++){
        if(lp->bin_unocc[iwunocc] <= e_unocc[0][ik][ic]  &&  e_unocc[0][ik][ic] < lp->bin_unocc[iwunocc+1])
          this_nunocc += 1;
      } // done calculating unoccupied states in this window
      // now we save eigenvalues in this window to Eunocc
      double Eunocc[this_nunocc];
      double Occp_unocc[this_nunocc];
      int Eunoccidx[this_nunocc];
      counter = 0;
      for (int ic=0; ic<nunocc; ic++){
        if(lp->bin_unocc[iwunocc] <= e_unocc[0][ik][ic] &&  e_unocc[0][ik][ic] < lp->bin_unocc[iwunocc+1]){
          Eunocc[counter] = e_unocc[0][ik][ic];
          // Reset occupation for metal to prevent some crazy occupation numbers (like 1E-123)
          Occp_unocc[counter] = occ_occ[0][ik][nocc+ic];
          if(p_metal){
            if ( abs ( Occp_unocc[counter]-1.0 ) < 1e-6){
              Occp_unocc[counter] = 1.0;
            }
            else if( abs ( Occp_unocc[counter] ) < 1e-6 ){
              Occp_unocc[counter] = 0.0;
            }
          }
          Eunoccidx[counter] = nocc+ic;
          counter += 1;
        }
      }

      // debugging message
#ifdef DEBUG
      printf("\n unoccupied states window %d , number of unoccupied states in this window: %d\n",iwunocc, this_nunocc);
      for (int i=0;i<this_nunocc;i++){
  printf("This unoccupied state: %d  Eunocc: %lg (%lg) , Eunoccidx: %d\n",i,Eunocc[i],Occp_unocc[i],Eunoccidx[i]);
      }
#endif

      // number of nodes and optimized a at this window pair (iwocc,iwunocc)
      int Ng = lp->Nnodes[iwocc][iwunocc];
      double opta = lp->opta[iwocc][iwunocc];

      //----------------------------------
      // Gaussian quadrature starts here
      //----------------------------------
      // loop over the GL nodes
      double gap = lp->gap;
      for (int n=0; n<Ng; n++){
        double Wn = lp->w[Ng-1][n];
        double Xn = lp->n[Ng-1][n];
#if DEBUG4
        printf("\nL.w[%d][%d] = %lf, L.n[Ng-1][n] = %lf gap = %lf\n", Ng-1, n, Lw, Ln, gap);
#endif
      double factor =  (-4.0/opta/1) * Wn * exp( -( gap/opta -1.0 ) * Xn );

        // adding p_metal calculation here
        // set internal counter for p_metal calculations 
        // if p_metal is true, nloop = 2
        // if not p_metal, nloop = 1
        int nloop = 1;
        // now let's check if we have to include f(Ec) or f(Ev) below. 
        // if occupancies are not 1 or 0, then we perform metal calculations for this window pair
        if( p_metal ){
          double one = 1.0;
          double zero = 0.0;
          for (int i=0; i<this_nocc; i++){
            if (Occp_occ[i] != zero && Occp_occ[i] != one){
              nloop = 2;
              // report only once
            if( n==0 ){
                printf("Occp_occ[%d]=%f\n",i,Occp_occ[i]);
              }
              break;
            }
          }
          for (int i=0; i<this_nunocc; i++){
            if( Occp_unocc[i] != zero && Occp_unocc[i] != one){
              nloop = 2;
              // report only once
              if ( n==0 ){
                printf("Occp_unocc[%d]=%f\n",i,Occp_unocc[i]);
              }
              break;
            }
          }
        }

        // if not p_metal, then nloop is 1 (and assume the system is gapped)
        for ( int ploop=0; ploop<nloop; ploop++ ){
          complex *focc; //rho matrix
          complex *funocc; //rho-bar matrix
          focc = new complex[ndata*ndata];
          funocc = new complex[ndata*ndata];

          // Occupied states
          complex *_psis_occ1 = new complex[this_nocc*ndata];
          complex *_psis_occ2 = new complex[this_nocc*ndata];
          complex *_psis_unocc1 = new complex[this_nunocc*ndata];
          complex *_psis_unocc2 = new complex[this_nunocc*ndata];

          int region_ridx = 0;
          for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
            if(start_row == psi_cache->regions[r_i].second) {
              region_ridx = r_i;
              break;
            }

          int region_cidx = 0;
          for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
            if(start_col == psi_cache->regions[r_i].second) {
              region_cidx = r_i;
              break;
            }
          // printf("%d %d %d %d\n", thisIndex.x, thisIndex.y, region_ridx, region_cidx);
          for (int iv=0; iv<this_nocc; iv++){

            // copy of the occupied state wavefunction as it may change due to U process
            complex *psi_occ1 = new complex[ndata];
            complex *psi_occ2 = new complex[ndata];

            

            for (int i=0; i<ndata; i++){
              psi_occ1[i] = psi_cache->psis[ikq][Eoccidx[iv]][region_ridx*ndata+i];
              psi_occ2[i] = psi_cache->psis[ikq][Eoccidx[iv]][region_cidx*ndata+i];
            }
            // printf("%d %d %d %d\n", thisIndex.x, thisIndex.y, region_ridx, region_cidx);
#if 0
            // modify wavefunction if umklapp scattering applies
            int uklsq = uklpp[0]*uklpp[0] + uklpp[1]*uklpp[1] + uklpp[2]*uklpp[2];
            if ( uklsq != 0 ){
              modify_state_Uproc(psi_occ, uklpp, nfft, sys);
            }
#endif

            for (int i=0; i<ndata; i++){
              // for gapped system
              if (nloop==1){
                _psis_occ1[iv*ndata+i] = psi_occ1[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) );
                _psis_occ2[iv*ndata+i] = psi_occ2[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) );
              }
              // for metal
              if (nloop==2){
                if(ploop==0){
                  _psis_occ1[iv*ndata+i] = psi_occ1[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) ) * sqrt( Occp_occ[iv] );
                  _psis_occ2[iv*ndata+i] = psi_occ2[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) ) * sqrt( Occp_occ[iv] );
                }
                else if(ploop==1){
                  _psis_occ1[iv*ndata+i] = psi_occ1[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) );
                  _psis_occ2[iv*ndata+i] = psi_occ2[i] * sqrt( exp(-(lp->maxEocc-Eocc[iv])/opta*Xn) );
                }
              }
            }
            delete[] psi_occ1;
            delete[] psi_occ2;
          } // end occupied states
#ifdef USE_LAPACK
    char transformT = 'C'; // conjugate transpose (Hermitian conjugate)
    char transform = 'N'; // do nothing
    double Lalpha = double(1.0); // scale A*B by this factor
    double Lbeta = double(1.0); // scale initial value of C by this factor
    // LAPACK GEMM: C = alpha*A*B + beta*C
    myGEMM( &transform, &transformT, &ndata, &ndata, &this_nocc, &Lalpha, _psis_occ2, &ndata, _psis_occ1, &ndata, &Lbeta, focc, &ndata);
#else
    Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif

    // Unoccupied states

    for (int ic=0; ic<this_nunocc; ic++){
      // copy unoccupied states into the temporary array
      complex psi_unocc1[ndata];
      complex psi_unocc2[ndata];
      int region_ridx = 0;
      for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
        if(start_row == psi_cache->regions[r_i].second) {
          region_ridx = r_i;
          break;
        }

      int region_cidx = 0;
      for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
        if(start_col == psi_cache->regions[r_i].second) {
          region_cidx = r_i;
          break;
        }
      for (int i=0; i<ndata; i++){
        psi_unocc1[i] = psi_cache->psis[ik][Eunoccidx[ic]][region_ridx*ndata+i].conj();
        psi_unocc2[i] = psi_cache->psis[ik][Eunoccidx[ic]][region_cidx*ndata+i].conj();
      }

      for (int i=0; i<ndata; i++){
        // for gapped system
        if (nloop==1){
          _psis_unocc1[ic*ndata+i] = psi_unocc1[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta * Xn) );
          _psis_unocc2[ic*ndata+i] = psi_unocc2[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta * Xn) );
        }
        // for metal
        if (nloop==2){
          if(ploop==0){
            _psis_unocc1[ic*ndata+i] = psi_unocc1[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta*Xn) );
            _psis_unocc2[ic*ndata+i] = psi_unocc2[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta*Xn) );
          }
          else if(ploop==1){
            _psis_unocc1[ic*ndata+i] = psi_unocc1[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta*Xn) );
            _psis_unocc2[ic*ndata+i] = psi_unocc2[i] * sqrt( exp(-(Eunocc[ic]-lp->minEunocc)/opta*Xn) );
          }
        }
      }
    }// end unoccupied state
#ifdef USE_LAPACK
    myGEMM( &transform, &transformT, &ndata, &ndata, &this_nunocc, &Lalpha, _psis_unocc2, &ndata, _psis_unocc1, &ndata, &Lbeta, funocc, &ndata);
#else
    Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif

    // Update P
    for (int i=0; i<ndata*ndata; i++){
      // for the first ploop, the sign is "+"  (either gapped system or metal first loop)
      if(ploop==0){
        P_m[i] += ((factor * focc[i] * funocc[i])*100000)/100000;
      }
      // for the second ploop, the sign is "-"
      if(ploop==1){
        P_m[i] -= ((factor * focc[i] * funocc[i])*100000)/100000;
      }
    }
#if DEBUG4
//    printf("\nfocc[10] = %lf, %lf, funocc[0] = %lf, %lf\n", focc[10].re, focc[10].im, funocc[0].re, funocc[0].im);
    printf("\n[%d,%d, %d]P->m[0] = %lf, %lf factor = %lf\n", iwocc, iwunocc, n, P_m[0].re, P_m[0].im, factor);
#endif
    delete[] focc;
    delete[] funocc;
    delete[] _psis_occ1;
    delete[] _psis_occ2;
    delete[] _psis_unocc1;
    delete[] _psis_unocc2;
        }

      }

    }
  }

}
  if(thisIndex.x==0 && thisIndex.y==0)
  {
    printf("\n[Node-%d][%d,%d]P->m[0] = %lf %lf\n", CkMyNode(), thisIndex.x, thisIndex.y,P_m[0].re, P_m[0].im);
    fflush(stdout);
  }

  psi_cache->elements++;
  int node_elems = psi_cache->total_elements/CkNumNodes();
  int tens = node_elems/10;
  if(psi_cache->elements%tens == 0)
    CkPrintf("\n%d%% p-matrix done", 10*psi_cache->elements/tens);
for(int i=0;i<ndata*ndata;i++)
  data[i] = P_m[i].conj();
//CkPrintf("\nContrib %d,%d\n", thisIndex.x, thisIndex.y);
  // contribute(CkCallback(CkReductionTarget(Controller, newPMatrixComplete), controller_proxy)); // enable for N3 P
}


void PMatrix::sigma_init(std::vector<int> _accept) {
  accept = _accept;
  // printf("sigma %d %d %f\n", thisIndex.x, thisIndex.y, accept[0]);
  contribute(CkCallback(CkReductionTarget(Controller, sigma_n3_initialized), controller_proxy));
}

// atakan
void PMatrix::sigma() {
  GWBSE *gwbse = GWBSE::get();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  WINDOWING* WIN = psi_cache->getWin();
  // WINDOWING WIN;
  // WIN->read_from_file();
  // WINDOWING* WIN;
  // System energies, parameters
  int nq = 1; //gwbse->gw_parallel.K; // for now nq = nkpt TODO kayahans hardcoded
  int nocc = gwbse->gw_parallel.L;
  int nunocc = gwbse->gw_parallel.M;
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;
  int nkpt = gwbse->gw_parallel.K;
  int nspin = 1; // TODO (kayahans) hardcoded
  // WINDOWING WIN(e_occ, e_unocc); // = WINDOWING(gwbse->gw_epsilon.Eocc, gwbse->gw_epsilon.Eunocc);
  // WIN.read_from_file(); 
  // if (thisIndex.x == 0 && thisIndex.y ==0) {
  //   WIN.printparameters();
  // }
  
  // Set up regioning data
  psi_ndata_local = config.tile_rows;
  tile_size = config.tile_rows * config.tile_cols;
  region_ridx = 0;
  // for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
  //   if(start_row == psi_cache->regions[r_i].second) {
  //     region_ridx = r_i;
  //     break;
  //   }

  region_cidx = 0;
  // for(int r_i=0;r_i<psi_cache->regions.size();r_i++)
  //   if(start_col == psi_cache->regions[r_i].second) {
  //     region_cidx = r_i;
  //     break;
  //   }
  // printf("xy %d %d %d %d %d %d\n", thisIndex.x, thisIndex.y, region_ridx, region_cidx, psi_ndata_local, tile_size);
  sigma_m = new complex[tile_size]; // sigma rr' tile
  for (int idata = 0; idata < tile_size; idata++) {
    sigma_m[idata] = 0.0;
  }  
  // // Sigma mtxels, frequencies
  double w = 0.23;
  std::vector<SIGMAINDICES> sigma_index_list;
  SIGMAINDICES n44(w, 4,4);
  SIGMAINDICES n55(w, 5,5);
  sigma_index_list.push_back(n44);
  sigma_index_list.push_back(n55);
  int ik = 0;
  int is = 0;
  

  
  // 1. Loop over occ/unocc
  for (int bloop=0; bloop < 2; bloop++) {
    bool bIsOccupied;
    if (bloop == 0)
      bIsOccupied = true;
    else
      bIsOccupied = false;
    
    // 2. Loop over q-points
    // printf("Nq sigma %d\n", nq);
    for (int iq = 0; iq < nq; iq++) {
      unsigned ikq;
      int umklapp[3];
      kqIndex(ik, ikq, umklapp);
      double* gpp_eig = psi_cache->get_gpp_eige(iq);
      double* gpp_omsq = psi_cache->get_gpp_omsq(iq);
      int ng = psi_cache->get_ng();
      // WIN->sigma_win(w, gpp_omsq, ng);
      // WIN->searchwins("Sigma"); // on the fly windowing optimization
      if (thisIndex.x == 0 && thisIndex.y == 0) {
        WIN->printparameters();
      }      
      // printf("xy %d %d %f %f %f %f\n", thisIndex.x, thisIndex.y, gpp_eig[0], gpp_eig[130], gpp_omsq[0], gpp_omsq[130]);
      std::vector<double> psi_eig;// shifted state eigenvalues
      double* _eig_occ = e_occ[is][ikq];
      double* _eig_unocc = e_unocc[is][ikq];
      std::copy(_eig_occ, _eig_occ+nocc, back_inserter(psi_eig));
      std::copy(_eig_unocc, _eig_unocc+nunocc, back_inserter(psi_eig));

      int state_start_index = 0;
      int state_end_index = nocc;
      if (!bIsOccupied) {
        state_start_index = nocc;
        state_end_index   = nunocc; // opts.nstateCh;  // TODO (kayahans)
      }
      // 3. loop over window pairs
      for (WINPAIR winpair : WIN->winpairs) {
        std::vector<double> state_e;  // E_v - w or w - E_c
        std::vector<double> pp_wp;     // PP energies w
        std::vector<int> state_idx;   // indexes of the state psi in window
        std::vector<int> pp_idx;      // indexes of the B in window
        double shiftedStateEnergy = 0.;
        int window_index = 0;
        int sigma_window_index;
        // printf("%d %d %f %f %f %f\n", thisIndex.x, thisIndex.y, winpair.w1[0], winpair.w1[1], winpair.w2[0], winpair.w2[1]);
        // istate loop: start accumulating states in the window pair
        for (int istate = state_start_index; istate < state_end_index; istate++) {
          // printf("%d %d %d %f %f\n", thisIndex.x, thisIndex.y, istate, psi_eig[istate] - w, w);
          if (bIsOccupied) {
            shiftedStateEnergy = psi_eig[istate] - w;  // E_v - w
            sigma_window_index = 0;
          } else {
            shiftedStateEnergy = w - psi_eig[istate];  // w - E_c
            sigma_window_index = 1;
          }
          if (winpair.in_window(shiftedStateEnergy, window_index, sigma_window_index)) {
            state_e.push_back(shiftedStateEnergy);
            state_idx.push_back(istate);
          }
        }  // end istate loop: accumulate states in state_e and state_idx
        
        double wpsq = 0;
        double wp = 0;
        int ipp = 0;
        window_index = 1;
        nfft = gwbse->gw_parallel.fft_nelems;
        int ndata = nfft[0] * nfft[1] * nfft[2];
        // i_ndata loop: start accumulate PP in the window pair
        for (int i_ndata = 0; i_ndata < ndata; i_ndata++) {
          if (accept[i_ndata]) {
            wpsq = gpp_omsq[ipp];
            wp = sqrt(wpsq);
            // Select w2>0 windows only!
            // printf("ipp %d wp %f \n", ipp, wp);
            if (wpsq > 0.0) {
              if (winpair.in_window(wp, window_index, sigma_window_index)) {
                pp_wp.push_back(wp);
                pp_idx.push_back(ipp);
              }
            }
            ipp++;
          }
        }  // end i_ndata accumulate pp loop
        
        // Important!
        // We want the denominator to be always negative or positive.
        // In this code it is assumed to be always negative for GL nodes.
        // Ev - w  and w - Ec are always negative
        // Normally the terms in the denominator are 1/(w - Ev + wp) and 1/(w - Ec - wp)
        // To keep the 1/(ai-bj) form, we rewrite the as:
        // -1/((Ev-w) - wp)
        // When 1/(ai-bj) is always negative, it means that zeta should be negative number.
        // Hence, in the succeeding functions:
        //  cubic_sigma_per_window
        //  PerformPPEnergySumThisNode
        //  PerformStateSumThisNode
        // we use zeta = - winpair.zeta;
        if (thisIndex.x == 0 && thisIndex.y==0) {
          printf("Window %f %f %f %f %d %d %d %d \n", winpair.w1[0], winpair.w1[1], winpair.w2[0], winpair.w2[1], winpair.nodes.size(), state_e.size(), pp_wp.size(), ipp);
        }
        
        // 4. If any states/PP are in this window pair
        if (state_e.size() > 0 && pp_wp.size() > 0) {
          // cubic_sigma_per_window(is, ik, iq, ikq, winpair, state_e, state_idx, pp_wp, pp_idx, umklapp);
          // for (int idata = 0; idata < ndata*ndata; idata++) {
          //   sigmat->m[idata] = 0.0;
          // }
          // 5. loop over GL quadrature
          for (int inode = 0; inode < winpair.nodes.size(); inode++) {
            // Quadrature constants
            int nloops;
            bool bIsFirstLoop;
            double iQuadFactor;
            const bool bIsHGL         = winpair.isHgl;
            const double tau_u        = winpair.nodes[inode];
            const double w_u          = winpair.weights[inode];
            const double zeta         = winpair.zeta;
            const double gap          = winpair.gap;
            const double gamma        = winpair.gamma;
            const double state_emax   = winpair.w1[1];
            const double wp_min       = winpair.w2[0];
            const int Na              = state_idx.size();
            const int Nb              = pp_idx.size();
            // printf("Pindex %d %d %f %f %f %f %d %f %f %f %f %f\n", thisIndex.x, thisIndex.y, winpair.w1[0], winpair.w1[1], winpair.w2[0], winpair.w2[1], inode, w_u, tau_u, zeta, gap, gamma);
            // Hermite GL or regular GL?
            if (bIsHGL) {
              iQuadFactor = w_u * winpair.gamma;  // MJ used -gamma here, check why and how it differs from zeta
              nloops = 2;
              bIsFirstLoop = true;
            } else {
              iQuadFactor = w_u * zeta * exp(-tau_u*(gap*zeta - 1.0));
              nloops = 1;
              bIsFirstLoop = false;
            }
            // 6. iLoop : HGL needs to calcualte B_p and f_n twice (Eq. 35)
            
            for (int iLoop=0; iLoop < nloops; iLoop++) {
              
              F_m = new complex[tile_size]; // psi rr' tile
              B_m = new complex[tile_size]; // gpp rr' tile

              // start B^p_{rr'}  Eq. 35
              // PerformPPEnergySumThisNode(is, iq, pp_e, pp_idx, winpair, inode, bIsHGL, bIsFirstLoop);
              // 7a. Loop over PP moves with energies in window pair
              for (int j = 0; j < Nb; j++) { 
                double omsq_j = gpp_omsq[pp_idx[j]];
                // if (1) {
                if (omsq_j > 0) {
                  double sigma_j = gpp_eig[pp_idx[j]];
                  double omega_j = sqrt(omsq_j);
                  // double omega_j = 1000000;

                  double factor  = sigma_j*omega_j/2.0;
                  // double factor  = sigma_j/2.0;
                  // printf("omega %f pp_wp %f \n", omega_j, pp_wp[j]);
                  if (!bIsHGL) {  // GL factor
                    factor *= exp(-zeta*(omega_j - wp_min)*tau_u);
                  } else {  // HGL factor
                    if (bIsFirstLoop) {
                      factor *= sin(omega_j*tau_u*gamma);
                    } else {
                      factor *= cos(omega_j*tau_u*gamma);
                    }  // end if
                  }  // end if

                  complex* B_r = psi_cache->get_gpp_eigv(iq, pp_idx[j]);

#ifdef USE_LAPACK
                  char transformT = 'C';
                  char transform  = 'N';
                  int Lone        = 1;
                  double one    = 1.0;
                  // myGEMM( &transform, &transformT, &ndata, &ndata, &this_nocc, &Lalpha, _psis_occ2, &ndata, _psis_occ1, &ndata, &Lbeta, focc, &ndata);
                  myGEMM(&transform, &transformT, &psi_ndata_local, &psi_ndata_local, &Lone,  &factor, B_r, &psi_ndata_local, B_r, &psi_ndata_local, &one, B_m, &psi_ndata_local);
#else
                Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif
                }
              }  // end B^p_{rr'}  Eq. 35 (7a)
  
              // A_{r,r'} Eq. 35
              // PerformStateSumThisNode(is, ik, ikq, state_e, state_idx, winpair, inode, bIsHGL, bIsFirstLoop, uklpp);
              // 7b. Loop over state energies in window pair
              for (int i_s = 0; i_s < Na; i_s++) {
                complex* psiRkq_i = new complex[psi_ndata_local];
                complex *fr1 = new complex[psi_ndata_local]; //rho-bar matrix
                for (int idx=0; idx<psi_ndata_local; idx++){
                  psiRkq_i[idx] = psi_cache->psis[ikq][state_idx[i_s]][region_ridx*psi_ndata_local+idx];
                  fr1[idx] = psi_cache->psis[ikq][state_idx[i_s]][region_ridx*psi_ndata_local+idx];
                }
                
                // TODO later
                // Need to modify psikq if umklapp vector is not zero
                // if( uklpp[0]!=0 || uklpp[1]!=0 || uklpp[2]!=0 ){
                //   complex* psikq_tmp[tile_size];
                //   // TODO modify_state_Uproc(fr, uklpp, nfft, sys);
                // }

                double factor;
                if (!bIsHGL) {
                  factor = exp(-zeta * (state_emax - state_e[i_s])*tau_u);
                  // factor = exp(-zeta * (state_emax - state_e[i])*tau_u);
                } else {  // HGL
                  if (bIsFirstLoop) {
                    factor = cos(tau_u*state_e[i_s]*gamma);
                  } else {
                    factor = sin(tau_u*state_e[i_s]*gamma);
                  }
                }
                // printf("i: %d, inode: %d, State factor: %lf, tau_u: %lf, zeta %lf\n", i, inode, factor, tau_u, zeta);
                // LAPACK
                char transformT = 'C';
                char transform  = 'N';
                int Lone        = 1;
                double beta    = 1.0;
#ifdef USE_LAPACK
                myGEMM(&transform, &transformT, &psi_ndata_local, &psi_ndata_local, &Lone,  &factor, fr1, &psi_ndata_local, fr1, &psi_ndata_local, &beta, F_m, &psi_ndata_local);
#else
                Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif    
                // delete [] psiRkq_i;
                delete [] fr1;
              } // end A^p_{rr'}  Eq. 35 (7b)

              // for (int idata = 0; idata < tile_size; idata++) {
              //     sigma_m[idata] += iQuadFactor * F_m[idata] * B_m[idata];
              // }
              
              delete [] B_m;
              delete [] F_m;
            } // end iLoop (6)
          }  // end quadrature loop (5)
        } // end if (4)
      }  // end winpair loop (3)
    }  // end iloop (2)
  } // end bloop (1)


  // Now multiply with state vectors <psi^*|Sigma|psi> to get sigma energies
  int i = 0;
  for (SIGMAINDICES iwn12 : sigma_index_list) {
    for (std::pair<int, int> in12 : iwn12.n12) {
      complex psiRk_n1[psi_ndata_local];
      complex psiRk_n2[psi_ndata_local];
      for (int i=0; i<psi_ndata_local; i++){
        psiRk_n1[i] = psi_cache->psis[ik][in12.first][region_ridx*psi_ndata_local + i];
        psiRk_n2[i] = psi_cache->psis[ik][in12.second][region_cidx*psi_ndata_local + i];
      }
      char transformT = 'C';
      char transform  = 'N';
      int    Lone   = 1;
      double beta   = 0.0;
      double factor = 1.0;
      complex sigma_n2[tile_size] = {0.0};
      complex n1_sigma_n2[1]  = {0.0};
      myGEMM(&transform,  &transform, &psi_ndata_local, &Lone, &psi_ndata_local, &factor, sigma_m, &psi_ndata_local, psiRk_n2, &psi_ndata_local, &beta, sigma_n2, &psi_ndata_local);
      myGEMM(&transformT, &transform, &Lone,  &Lone, &psi_ndata_local, &factor, psiRk_n1, &psi_ndata_local, sigma_n2, &psi_ndata_local, &beta, n1_sigma_n2, &Lone);
      sigma_m[i] = n1_sigma_n2[0].re;
      i++;
    }
  }
  contribute(CkCallback(CkReductionTarget(Controller, sigma_n3_complete), controller_proxy));
}

void PMatrix::compute_fr(complex* fr, complex* psikq, const int uklpp[3]){
  // set f_r
  for(int i=0; i<tile_size; i++){
      fr[i] = psikq[i];
  }
  // Need to modify psikq if umklapp vector is not zero
  if( uklpp[0]!=0 || uklpp[1]!=0 || uklpp[2]!=0 ){
    complex* psikq_tmp[tile_size];
    // TODO modify_state_Uproc(fr, uklpp, nfft, sys);
  }
}


void PMatrix::generateEpsilon(CProxy_EpsMatrix proxy, std::vector<int> accept){
  int inew = 0;
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  std::vector<double> vcoulb = psi_cache->getVCoulb();
  for(int i=0;i<thisIndex.x;i++)
    if(accept[i])
      inew++;

  if(accept[thisIndex.x]){
    int jnew_local = 0;
    int jnew = 0;
    Phase3Message *msg;
    msg = new(eps_cols)Phase3Message();
    for(int j=0;j<config.tile_cols;j++){

      if(accept[j]){
        msg->data[jnew_local] = data[j];
        msg->data[jnew_local] *= -1 * sqrt(vcoulb[inew]) * sqrt( vcoulb[jnew]);
        if ( inew == jnew )
          msg->data[jnew_local] += double(1);
        jnew_local++;
        jnew++;
        if(jnew_local == eps_cols){
          int dest_chare_x = inew/eps_rows;
          int dest_chare_y = (jnew/eps_cols)-1;
          msg->start_i = inew%eps_rows;
          msg->start_j = 0;
          msg->end_i = inew%eps_rows;
          msg->end_j = (jnew-1)%eps_cols;
          proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
          jnew_local = 0;
          msg = new(eps_cols)Phase3Message();  
        }
      }
    } 
    if(jnew_local != eps_cols){
      for(int i=jnew_local; i<eps_cols; i++){
        msg->data[i] = 0;
        jnew++;
      }
        
      int dest_chare_x = inew/eps_rows;
      int dest_chare_y = (jnew/eps_cols)-1;
      msg->start_i = inew%eps_rows;
      msg->start_j = 0;
      msg->end_i = inew%eps_rows;
      msg->end_j = (jnew-1)%eps_cols;
      proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
    }
  }
 
  int counter = 0;
  if(thisIndex.x == config.tile_cols-1){
    
    int padded_send_size = inew + (eps_cols - (inew%eps_cols));
    int inew_start = inew;
    if(accept[thisIndex.x]) inew_start+=1;
    int remainder = eps_cols - inew_start%eps_cols;
//    CkPrintf("\nSending to i=%d,j=%d\n", remainder, padded_send_size);
    for(int i=inew_start;counter<remainder;i++){
      counter++;
      Phase3Message *msg;
      msg = new(eps_cols)Phase3Message();
      int jnew_local = 0;
      int jnew = 0;
      msg->start_i = i%eps_rows;

      for(int j=0;j<padded_send_size;j++){
        msg->data[jnew_local++] = 0;
        jnew++;
        if(jnew_local == eps_cols){
          int dest_chare_x = i/eps_rows;
          int dest_chare_y = (jnew/eps_cols)-1;
          msg->start_i = i%eps_rows;
          msg->start_j = 0;
          msg->end_i = i%eps_rows;
          msg->end_j = (jnew_local-1)%eps_cols;
          proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
          jnew_local = 0;
          msg = new(eps_cols)Phase3Message();
        }
      }
    }
  }
}

void PMatrix::reportPTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::sum_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportPTime(NULL), controller_proxy));
  contribute(msg);
}

void PMatrix::applyFs() {
  double start = CmiWallTimer();

  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();

#ifdef USE_LAPACK
  // Common variables for both ZGERC and ZGEMM
  int M = config.tile_rows, N = config.tile_cols;
  complex alpha = -1.0;
#ifdef USE_ZGEMM
  int K = L; // If using ZGEMM, we compute all outer products with one call
  int LDF = config.mat_rows; // Leading dimension of fs
  complex beta = 1.0;
  char opA = 'N', opB = 'C';
  complex* fs = psi_cache->getF(0,completed_chunks);
  ZGEMM(&opA, &opB, &N, &M, &K,
    &alpha, &(fs[start_col]), &LDF,
    &(fs[start_row]), &LDF,
    &beta, data, &N);
#else
  int K = 1; // If using ZGERC, we compute each outer product one at a time
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l,completed_chunks);
    ZGERC(&N, &M, &alpha, &(f[start_col]), &K, &(f[start_row]), &K, data, &N);
  }
#endif // endif for ifdef USE_ZGEMM
#else
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l,completed_chunks);
    for (int r = 0; r < config.tile_rows; r++) {
      for (int c = 0; c < config.tile_cols; c++) {
        data[IDX(r,c)] += f[r+start_row]*f[c+start_col].conj() * -1.0;
      }
    }
  }
#endif // endif for ifdef USE_LAPACK
  completed_chunks++;
  contribute(CkCallback(CkReductionTarget(Controller, psiComplete), controller_proxy));
  total_time += CmiWallTimer() - start;
}

void PMatrix::fftRows(int direction) {
  if (config.chareCols() != 1) {
    CkAbort("FFT not supported for 2D decompositions\n");
  }
  // FFT each row stored in this chare
  for (int i=0; i < config.tile_rows; i++){
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();

    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &data[IDX(i,0)], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(config.tile_cols, out_pointer, &data[IDX(i,0)], 1);
  }

  // Let the controller know we have completed the fft
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
}

// TODO: These methods shouldn't be part of PMatrix, and should also just be
// computed once at startup.
void PMatrix::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
  GWBSE* gwbse = GWBSE::get();

  // temporary space to save k/q/k+q vectors
  double *this_k, *this_q;
  double k_plus_q[3], k_plus_q_orig[3];

  this_k = gwbse->gwbseopts.kvec[ikpt];
  this_q = gwbse->gwbseopts.qvec[qindex];

  for (int i=0; i<3; i++) {
    // calculate k+q vector
    k_plus_q[i] = this_k[i] + this_q[i]; // k+q vector
    k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp process
    // if not 0 =< k+q [i] <1, adjust k+q so that k+q[i] is in the Brillouine zone
    if ( k_plus_q[i] >= 1 ) {
      k_plus_q[i] -= 1;
    }
    else if( k_plus_q[i] < 0 ){
      k_plus_q[i] += 1;
    }
  }

  // find k+q vector index
  for (int kk=0; kk < gwbse->gwbseopts.nkpt; kk++) {
    bool match = true;
    this_k = gwbse->gwbseopts.kvec[kk];
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      if (this_k[i] != k_plus_q[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ikq = kk;
      break;
    }
  }
  // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklapp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }

}

void PMatrix::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

  if (uklpp[0]==0 && uklpp[1]==0 && uklpp[2]==0){
    // do nothing
  }
  else{
    GWBSE *gwbse = GWBSE::get();
    int* nfft;
    nfft = gwbse->gw_parallel.fft_nelems;
    double *a1, *a2, *a3, *b1, *b2, *b3;
    a1 = gwbse->gwbseopts.a1;
    a2 = gwbse->gwbseopts.a2;
    a3 = gwbse->gwbseopts.a3;
    b1 = gwbse->gwbseopts.b1;
    b2 = gwbse->gwbseopts.b2;
    b3 = gwbse->gwbseopts.b3;
    double lattconst = gwbse->gwbseopts.latt;

    double rijk, G0, phase;
    unsigned counter = 0;
    for(int i=0; i<nfft[0]; i++){
      for(int j=0; j<nfft[1]; j++){
        for(int k=0; k<nfft[2]; k++){
          phase = 0;
          for (int l=0; l<3; l++){
            rijk = a1[l]*i/nfft[0] + a2[l]*j/nfft[1] + a3[l]*k/nfft[2];
            G0 = b1[l]*uklpp[0] + b2[l]*uklpp[1] + b3[l]*uklpp[2];
            G0 *= -2*M_PI/lattconst;
            phase += rijk*G0;
          }
          umklapp_factor[counter].re = cos(phase);
          umklapp_factor[counter].im = sin(phase);
          counter += 1;
        }// end k loop
      }// end j loop
    }// end i loop
  }//end if-else statement

}//end function

//each chare in the matrix gives its row and column ranges to its local PsiCache
void PMatrix::registerTileSections() {
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  psi_cache->setRegionData(this, start_row, start_col, config.tile_rows, config.tile_cols);

  contribute(CkCallback(CkReductionTarget(Controller,registrationComplete), controller_proxy));
}//end function

#include "pmatrix.def.h"
