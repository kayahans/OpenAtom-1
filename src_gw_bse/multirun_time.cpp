/*  Example code to demonstrate time sharing interoperability between MPI and Charm
    Author - Nikhil Jain
    Contact - nikhil@illinois.edu
 */

//standard header files
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//header files for libraries in Charm I wish to use with MPI
#include "main/hi.h"
//header file from Charm needed for Interoperation
#include "mpi-interoperate.h"
#include "diagonalizer/diagonalizer.h"
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

  int peid, numpes;
  
  MPI_Comm newComm;

  int iam, nprocs;
  int ictxt, nprow, npcol, myrow, mycol;
  int np, nq, nb, n;
  int mpA, nqA;
  int i, j, k, info, itemp, seed, lwork, min_mn;
  int descA[9], descZ[9];
  // double *A, *Z, *work, *W;
  int izero=0,ione=1;
  double mone=(-1.0e0),pone=(1.0e0),dzero=(0.0e0);
  double MPIt1, MPIt2, MPIelapsed;
  char jobz, uplo;

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
    const int total_iter  = diagData->num_handoff;
    int iternum = 1;
    
    for(; iternum <= total_iter ; iternum++) {
        nprow   = diagData->nprow;
        npcol   = diagData->npcol;
        n       = diagData->n;
        nb      = diagData->nb;
        jobz    = 'V';
        uplo    = 'U';

        if (nprow*npcol>numpes) {
            if (peid==0) {
            printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
            }
            printf("%d %d\n", nprow, npcol);
            printf(" **** Bye-bye ***\n");
            MPI_Finalize();
            exit(1);
        }        
        // printf("Handoff to MPI for diagonalization %d\n", CkMyPe());
        Cblacs_pinfo(&iam, &nprocs);
        Cblacs_get(0, 0, &ictxt);
        Cblacs_gridinit(&ictxt, "R", nprow, npcol);
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        printf("[DIAGONALIZER] on proc %dx%d of %dx%d nb %d n %d iternum %d\n",myrow+1, mycol+1, nprow, npcol, nb, n, iternum);

        MPI_Barrier(MPI_COMM_WORLD);
        StartHi(16);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(!peid)
      printf("Invoke charmlib exit\n");
    CharmLibExit();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();



  // //do some MPI work
  // for(int q=0; q<8; q++) {
  //   for(int i=0; i<5; i++) {
  //     if(peid % 2 == 0) {    
  //       MPI_Send(&peid, 1, MPI_INT, peid+1, 808, MPI_COMM_WORLD);
  //     } else {
  //       int recvid = 0;
  //       MPI_Status sts;
  //       MPI_Recv(&recvid, 1, MPI_INT, peid-1, 808, MPI_COMM_WORLD, &sts);
  //     }
  //   }

  //   //Hello
  // //  HelloStart(5);
  //   MPI_Barrier(newComm);

  //   for(int i=0; i<5; i++) {
  //     if(peid % 2 == 1) {    
  //       MPI_Send(&peid, 1, MPI_INT, peid-1, 808, MPI_COMM_WORLD);
  //     }  else {
  //       int recvid = 0;
  //       MPI_Status sts;
  //       MPI_Recv(&recvid, 1, MPI_INT, peid+1, 808, MPI_COMM_WORLD, &sts);
  //     }
  //   }

  //   StartHi(16);
  //   MPI_Barrier(newComm);
  // }
  // if(!peid)
  //   printf("Invoke charmlib exit\n");

  // CharmLibExit();

  // //final synchronization
  // MPI_Barrier(MPI_COMM_WORLD);
// #else
//   if (peid == 0) {
//     printf("This test program is a no-op with non-MPI builds\n");
//   }
// #endif

  // MPI_Finalize();
  return 0;  
}
