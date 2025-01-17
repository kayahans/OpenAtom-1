      INTEGER FUNCTION ILAQLC( M, N, A, LDA )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2.2)                        --
*
*  -- June 2010                                                       --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      REAL*16  A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ILAQLC scans A for its last non-zero column.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input) REAL*16 array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*16 ZERO
      PARAMETER ( ZERO = 0.0Q+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILAQLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAQLC = N
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ILAQLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILAQLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
