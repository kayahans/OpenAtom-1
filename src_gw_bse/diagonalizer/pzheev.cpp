#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
//header file from Charm needed for Interoperation
#include "mpi-interoperate.h"
#include "diagonalizer.h"
using namespace std;

extern "C" {
// BLACS
  void Cblacs_pinfo(int* mypnum, int* nprocs);
  void Cblacs_get(int context, int request, int* value);
  int Cblacs_gridinit(int* context, const char * order, int np_row, int np_col);
  void Cblacs_gridinfo(int context, int* np_row, int* np_col, int* my_row, int* my_col);
  void Cblacs_gridexit(int context);
  void Cblacs_exit(int error_code);

// SCALAPACK
  int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
  void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                 int *ictxt, int *lld, int *info);
  void pzheev_(char *jobz, char *uplo, int *n, std::complex<double> *a, int *ia, int *ja,
                 int *desca, double *w, std::complex<double> *z, int *iz, int *jz, int *descz,
                 std::complex<double> *work, int *lwork, double *rwork, int *lrwork, int *info);
  void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja,
                  int *desca, double *w, double *z, int *iz, int *jz, int *descz,
                  double *work, int *lwork, int *info);
}

static int max(int a, int b) {
  if (a>b) return(a); else return(b);
}

static int min(int a, int b) {
  if (a<b) return(a); else return(b);
}

extern diagData_t* diagData;


int main(int argc, char **argv){
  printf("\nCalling MPI code\n");

  // MPI parameters
  int peid, numpes;
  MPI_Comm newComm;
  int iam, nprocs;
  int ictxt, nprow, npcol, myrow, mycol;
  int np, nq, nb, n;
  int mpA, nqA;
  int izero=0,ione=1;
  double mone=(-1.0e0),pone=(1.0e0),dzero=(0.0e0);
  double MPIt1, MPIt2, MPIelapsed;

  // SCALAPACK
  // http://www.netlib.org/scalapack/explore-html/d3/d05/pzheev_8f_source.html
  char jobz, uplo;
  int ia, info, iz, ja, jz, lrwork, lwork, min_mn;
  int descA[9], descZ[9];
  double *rwork, *W;
  std::complex<double> *A, *Z, *work;
  
  //basic MPI initilization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &peid);
  MPI_Comm_size(MPI_COMM_WORLD, &numpes);

  if(numpes % 2 != 0){
    if(peid==0){
      printf("This test program must be run with number of procs = 2x\n");
    }
    MPI_Finalize();
    return 1;
  }

  // Currently, only MPI builds support this style of interop
  // #if CMK_CONVERSE_MPI
  MPI_Comm_split(MPI_COMM_WORLD, 1, peid, &newComm);

  // Initial call to Charm++
  MPI_Barrier(MPI_COMM_WORLD);
  CharmLibInit(newComm, argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);


  // Return to MPI
  const int total_iter  = diagData->num_handoffs;
  int iternum = 1;

  for ( ; iternum <= total_iter; iternum++) {
    printf("Handoff to MPI for diagonalization at PE %d\n", CkMyPe());
    nprow   =  diagData->nprow;
    npcol   =  diagData->npcol;
    n       =  diagData->n;
    nb      =  diagData->nb;
    jobz    = 'V';
    uplo    = 'U';
    // printf("nprow %d npcol %d n %d nb %d pe %d \n", )
    if (nprow * npcol > numpes) {
        if (peid == 0) {
          printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
        }
        printf("%d %d\n", nprow, npcol);
        printf(" **** Bye-bye ***\n");
        MPI_Finalize();
        exit(1);
    }

    Cblacs_pinfo(&iam, &nprocs);
    Cblacs_get(0, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "R", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    int size = diagData->inputsize;
    printf("[DIAGONALIZER] on proc %dx%d of %dx%d inputsize %d nb %d n %d iternum %d\n",myrow+1, mycol+1, nprow, npcol, size, nb, n, iternum);
    // for (int i=0; i < 10 ; i++) {
    //   printf("[DIAGONALIZER] on proc %dx%d of %dx%d inputsize %d nb %d n %d iternum %d, %.8e + %.8e\n",myrow+1, mycol+1, nprow, npcol, size, nb, n, iternum, diagData->input[i].real(), diagData->input[i].imag());
    // }
    mpA = numroc_(&n, &nb, &myrow, &izero, &nprow);
    nqA = numroc_(&n, &nb, &mycol, &izero, &npcol);
    printf("pa qa %d, %d %dx%d \n", mpA, nqA, myrow+1, mycol+1);
    min_mn = n;

    A = (std::complex<double> *)calloc(mpA*nqA,sizeof(std::complex<double>));
    if (A == NULL) { printf("error of memory allocation A on proc %dx%d\n",
      myrow, mycol);
      exit(0); }
    Z = (std::complex<double> *)calloc(mpA*nqA,sizeof(std::complex<double>));
    if (Z == NULL) { printf("error of memory allocation VT on proc %dx%d\n",
      myrow, mycol);
      exit(0);}
    W = (double *)calloc(min_mn,sizeof(double));
    if (W == NULL) { printf("error of memory allocation S on proc %dx%d\n",
      myrow, mycol);
      exit(0); }
    
    {
      // copy to A, maybe memcpy will be faster
      int k = 0;
      for (int i = 0; i < mpA; i++) {
          for (int j = 0; j < nqA; j++) {
              A[k] = diagData->input[(i*nqA)+j];
              k++;
          }
      }
    }

    int itemp = max(1, mpA);
    descinit_(descA,  &n, &n, &nb, &nb,
      &izero, &izero, &ictxt, &itemp, &info);
    descinit_(descZ,  &n, &n, &nb, &nb,
      &izero, &izero, &ictxt, &itemp, &info);
    
    // Initialize work and rwork to get work sizes
    work = (std::complex<double> *)calloc(2,sizeof(std::complex<double>));
    if (work == NULL) {
      printf("error of memory allocation for work on proc %dx%d (1st time)\n",
      myrow, mycol);
      exit(0);
    }
    rwork = (double *)calloc(2,sizeof(double));
    if (rwork == NULL) {
      printf("error of memory allocation for rwork on proc %dx%d (1st time)\n",
      myrow, mycol);
      exit(0);
    }

    lwork = lrwork = -1;
    pzheev_(&jobz, &uplo, &n, A, &ione, &ione, descA, W, Z,
            &ione, &ione, descZ, work, &lwork, rwork, &lrwork, &info);
    
    // Resize work array 
    lwork = static_cast<int>(work[0].real()); 
    // printf("Work size on proc %dx%d: %d\n", myrow,mycol, lwork);
    free(work);
    work = (std::complex<double> *)calloc(lwork, sizeof(std::complex<double>));
    if (work == NULL) { printf("error of memory allocation work "
      "on proc %dx%d\n",myrow,mycol); exit(0);}

    // Resize rwork array 
    lrwork = static_cast<int>(rwork[0]) * 2; // Bug in Scalapack!
    // printf("Rwork size on proc %dx%d: %d\n", myrow,mycol, lrwork);
    free(rwork);
    rwork = (double *)calloc(lrwork, sizeof(double));
    if (rwork == NULL) { printf("error of memory allocation rwork "
      "on proc %dx%d\n",myrow,mycol); exit(0);}

    MPIt1 = MPI_Wtime();
    pzheev_(&jobz, &uplo, &n, A, &ione, &ione, descA, W, Z,
            &ione, &ione, descZ, work, &lwork, rwork, &lrwork, &info);
    // fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    MPIt2 = MPI_Wtime();
    MPIelapsed = MPIt2 - MPIt1;

    if (iam == 1) {
      int looptot = mpA * nqA;
      // diagData->eig_v = new double[looptot];
      // for (i = 0; i < mpA; i++) {
      //     for (j = 0; j < nqA; j++) {
      //         // diagData->eig_v[(i*nqA)+j] = Z[(i*nqA)+j];
      //         printf("[DIAGONALIZER] eig_v (%d,%d) %.8e %d, %d, %d\n",
      //           i, j, Z[(i*nqA)+j].real(), myrow+1, mycol+1, iam);
      //     }
      // }
      // diagData->eig_e = new double[mpA];
      for (int i = 0; i < min_mn; i++) {
        // diagData->eig_e[i] = A[i];
        printf("[DIAGONALIZER] eig_e %.8e %f %d, %d, %d\n",
          i, W[i], myrow+1, mycol+1, iam);
      }
      printf("n=%d\t(%d,%d)\t%d\tjobz=%c\t%8.2fs \n",
        n, nprow, npcol, nb, jobz, MPIelapsed);
    }
    free(work);
    free(rwork);
    free(W);
    free(Z);
    free(A);
    MPI_Barrier(MPI_COMM_WORLD);
    restartCharm();
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(!peid)
    printf("Invoke charmlib exit\n");

  CharmLibExit();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;  
}
