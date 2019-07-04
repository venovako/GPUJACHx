PROGRAM PP80
  USE HDF5
  USE H5LT

  IMPLICIT NONE

  REAL(10), PARAMETER :: ZERO =  0.0_10
  REAL(10), PARAMETER :: ONE  =  1.0_10
  REAL(10), PARAMETER :: TWO  =  2.0_10
  REAL(10), PARAMETER :: MONE = -1.0_10

  INTEGER, PARAMETER :: INRMSL = 11
  INTEGER, PARAMETER :: IDDIML =  4

  CHARACTER(LEN=256) :: FH5
  CHARACTER(LEN=256) :: GH5
  CHARACTER(LEN=256) :: RH5

  INTEGER(HID_T) :: FID
  INTEGER(HID_T) :: GID
  INTEGER(HID_T) :: DID
  INTEGER(HSIZE_T) :: DIMS1(1)
  INTEGER(HSIZE_T) :: DIMS2(2)

  INTEGER :: IDADIM(IDDIML)
  INTEGER :: LDA
  INTEGER :: M
  INTEGER :: N
  INTEGER :: P
  EQUIVALENCE (IDADIM(1),LDA), (IDADIM(2),M), (IDADIM(3),N), (IDADIM(4),P)

  LOGICAL :: FEXIST
  INTEGER :: I
  INTEGER :: J
  INTEGER :: K
  INTEGER :: INFO
  DOUBLE PRECISION :: SCL
  DOUBLE PRECISION :: NRMS(INRMSL)
  REAL(10) :: QTMP0, QMAXA, QSUMSQ

  REAL(10), ALLOCATABLE :: D(:)
  REAL(10), ALLOCATABLE :: U(:,:)
  REAL(10), ALLOCATABLE :: V(:,:)
  REAL(10), ALLOCATABLE :: MM(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: DD(:)
  DOUBLE PRECISION, ALLOCATABLE :: LAMBDA(:)
  DOUBLE PRECISION, ALLOCATABLE :: ABSERR(:)
  DOUBLE PRECISION, ALLOCATABLE :: RELERR(:)
  DOUBLE PRECISION, ALLOCATABLE :: DTMP2(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:)

#ifdef USE_FMA
  INTERFACE
     REAL(10) FUNCTION C_FMA(x, y, z) BIND(C,NAME='fmal')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
       REAL(10), VALUE :: x, y, z
     END FUNCTION C_FMA
  END INTERFACE
#endif

  EXTERNAL :: DLASRT

  IF (COMMAND_ARGUMENT_COUNT() .NE. 3) THEN
#ifdef MAGMA
     STOP 'pp80M.exe FH5 GH5 RH5'
#else
     STOP 'pp80J.exe FH5 GH5 RH5'
#endif
  END IF

  CALL GET_COMMAND_ARGUMENT(1,FH5)
  CALL GET_COMMAND_ARGUMENT(2,GH5)
  CALL GET_COMMAND_ARGUMENT(3,RH5)

  INQUIRE(FILE=FH5,EXIST=FEXIST)
  IF (.NOT. FEXIST) STOP 'FH5 nonexistent!'
  INQUIRE(FILE=RH5,EXIST=FEXIST)
  IF (.NOT. FEXIST) STOP 'RH5 nonexistent!'

  CALL h5open_f(INFO)
  IF (INFO .NE. 0) STOP 'Error initializing HDF5!'

  CALL h5fopen_f(RH5, H5F_ACC_RDONLY_F, FID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening RH5!'

  CALL h5gopen_f(FID, GH5, GID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening RH5/GH5!'

  DIMS1(1) = IDDIML
  CALL h5ltread_dataset_int_f(GID, 'IDADIM', IDADIM, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading IDADIM!'

  IF (LDA .LT. 0) STOP 'LDA < 0'
  IF (M .LT. 0) STOP 'M < 0'
  IF (N .LT. 0) STOP 'N < 0'
  IF (P .LT. 0) STOP 'N < 0'

  IF (LDA .LT. M) STOP 'LDA < M'
  IF (M .LT. N) STOP 'M < N'
  IF (N .LT. P) STOP 'N < P'
  IF (MOD(N,2) .NE. 0) STOP 'N ODD'

#ifdef MAGMA
  IF (P .NE. N) STOP 'P <> N'
#endif

  ALLOCATE(D(N))

  ALLOCATE(DD(N))
  DIMS1(1) = N
#ifdef MAGMA
  CALL h5ltread_dataset_double_f(GID, 'SIGMA', DD, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading SIGMA!'
  !$OMP PARALLEL DO PRIVATE(QTMP0)
  DO I = 1, N
     QTMP0 = REAL(DD(I),10)
     D(I) = QTMP0 * QTMP0
     DD(I) = DBLE(D(I))
  END DO
#else
  CALL h5ltread_dataset_double_f(GID, 'D', DD, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading D!'
  !$OMP PARALLEL DO
  DO I = 1, N
     D(I) = REAL(DD(I),10)
  END DO
#endif

  ALLOCATE(U(M,N))

  ALLOCATE(DTMP2(M,N))
  DIMS2(1) = M
  DIMS2(2) = N
#ifdef MAGMA
  CALL h5ltread_dataset_double_f(GID, 'U', DTMP2, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading U!'
#else
  CALL h5ltread_dataset_double_f(GID, 'G', DTMP2, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading G!'
#endif
  !$OMP PARALLEL DO PRIVATE(I)
  DO J = 1, N
     DO I = 1, M
        U(I,J) = REAL(DTMP2(I,J),10)
     END DO
  END DO
  DEALLOCATE(DTMP2)

  ALLOCATE(MM(M,M))

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, N
     MM(J,J) = MONE
     DO K = 1, M
#ifdef USE_FMA
        MM(J,J) = C_FMA(U(K,J), U(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) + U(K,J) * U(K,J)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, M
#ifdef USE_FMA
           MM(I,J) = C_FMA(U(K,J), U(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) + U(K,J) * U(K,I)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = J, N
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(1) = DBLE(QMAXA)
  WRITE (*,1) '|| U^T U - I ||_M', NRMS(1)

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, N
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, N
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(2) = DBLE(SQRT(QSUMSQ))
  WRITE (*,1) '|| U^T U - I ||_F', NRMS(2)

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, M
     MM(J,J) = MONE
     DO K = 1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(U(J,K), U(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) + U(J,K) * U(J,K)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, M
        MM(I,J) = ZERO
        DO K = 1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(U(J,K), U(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) + U(J,K) * U(I,K)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, M
     DO I = J, M
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(3) = DBLE(QMAXA)
  WRITE (*,1) '|| U U^T - I ||_M', NRMS(3)

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, M
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, M
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(4) = DBLE(SQRT(QSUMSQ))
  WRITE (*,1) '|| U U^T - I ||_F', NRMS(4)

  ALLOCATE(V(M,N))

  ALLOCATE(DTMP2(N,N))
  DIMS2(1) = N
  DIMS2(2) = N
#ifdef MAGMA
  CALL h5ltread_dataset_double_f(GID, 'VT', DTMP2, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading VT!'
#else
  CALL h5ltread_dataset_double_f(GID, 'V', DTMP2, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading V!'
#endif
  !$OMP PARALLEL DO PRIVATE(I)
  DO J = 1, N
     DO I = 1, N
        V(I,J) = REAL(DTMP2(I,J),10)
     END DO
  END DO
  DEALLOCATE(DTMP2)

  CALL h5gclose_f(GID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing RH5/GH5!'

  CALL h5fclose_f(FID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing RH5!'

#ifdef MAGMA
  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, N
     MM(J,J) = MONE
     DO K = 1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(J,K), V(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(J,K) * V(J,K)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(J,K), V(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(J,K) * V(I,K)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = J, N
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(5) = DBLE(QMAXA)
  WRITE (*,1) '|| V^T V - I ||_M', NRMS(5)

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, M
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, M
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(6) = DBLE(SQRT(QSUMSQ))
  WRITE (*,1) '|| V^T V - I ||_F', NRMS(6)

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, N
     MM(J,J) = MONE
     DO K = 1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(K,J), V(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(K,J) * V(K,J)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(K,J), V(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(K,J) * V(K,I)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = J, N
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(7) = DBLE(QMAXA)
  WRITE (*,1) '|| V V^T - I ||_M', NRMS(7)

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, N
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, N
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(8) = DBLE(SQRT(QSUMSQ))
  WRITE (*,1) '|| V V^T - I ||_F', NRMS(8)
#else
  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, P
     MM(J,J) = MONE
     DO K = 1, P
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(K,J), V(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(K,J) * V(K,J)
#endif
     END DO
     DO K = P+1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(-V(K,J), V(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) - V(K,J) * V(K,J)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, P
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(K,J), V(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(K,J) * V(K,I)
#endif
        END DO
        DO K = P+1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(-V(K,J), V(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) - V(K,J) * V(K,I)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = P+1, N
     MM(J,J) = ONE
     DO K = 1, P
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(K,J), V(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(K,J) * V(K,J)
#endif
     END DO
     DO K = P+1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(-V(K,J), V(K,J), MM(J,J))
#else
        MM(J,J) = MM(J,J) - V(K,J) * V(K,J)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, P
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(K,J), V(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(K,J) * V(K,I)
#endif
        END DO
        DO K = P+1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(-V(K,J), V(K,I), MM(I,J))
#else
           MM(I,J) = MM(I,J) - V(K,J) * V(K,I)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = J, N
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(5) = DBLE(QMAXA)
  IF (P .EQ. N) THEN
     WRITE (*,1) '|| V^T V - I ||_M', NRMS(5)
  ELSE
     WRITE (*,1) '|| V^T J V - J ||_M', NRMS(5)
  END IF

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, N
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, N
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(6) = DBLE(SQRT(QSUMSQ))
  IF (P .EQ. N) THEN
     WRITE (*,1) '|| V^T V - I ||_F', NRMS(6)
  ELSE
     WRITE (*,1) '|| V^T J V - J ||_F', NRMS(6)
  END IF

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, P
     MM(J,J) = MONE
     DO K = 1, P
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(J,K), V(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(J,K) * V(J,K)
#endif
     END DO
     DO K = P+1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(-V(J,K), V(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) - V(J,K) * V(J,K)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, P
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(J,K), V(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(J,K) * V(I,K)
#endif
        END DO
        DO K = P+1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(-V(J,K), V(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) - V(J,K) * V(I,K)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = P+1, N
     MM(J,J) = ONE
     DO K = 1, P
#ifdef USE_FMA
        MM(J,J) = C_FMA(V(J,K), V(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) + V(J,K) * V(J,K)
#endif
     END DO
     DO K = P+1, N
#ifdef USE_FMA
        MM(J,J) = C_FMA(-V(J,K), V(J,K), MM(J,J))
#else
        MM(J,J) = MM(J,J) - V(J,K) * V(J,K)
#endif
     END DO
     MM(J,J) = ABS(MM(J,J))
     DO I = J+1, N
        MM(I,J) = ZERO
        DO K = 1, P
#ifdef USE_FMA
           MM(I,J) = C_FMA(V(J,K), V(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) + V(J,K) * V(I,K)
#endif
        END DO
        DO K = P+1, N
#ifdef USE_FMA
           MM(I,J) = C_FMA(-V(J,K), V(I,K), MM(I,J))
#else
           MM(I,J) = MM(I,J) - V(J,K) * V(I,K)
#endif
        END DO
        MM(I,J) = ABS(MM(I,J))
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = J, N
        IF (MM(I,J) .GT. QMAXA) QMAXA = MM(I,J)
     END DO
  END DO

  NRMS(7) = DBLE(QMAXA)
  IF (P .EQ. N) THEN
     WRITE (*,1) '|| V V^T - I ||_M', NRMS(7)
  ELSE
     WRITE (*,1) '|| V J V^T - J ||_M', NRMS(7)
  END IF

  QTMP0 = ZERO
  QSUMSQ = ZERO
  DO J = 1, N
#ifdef USE_FMA
     QSUMSQ = C_FMA(MM(J,J), MM(J,J), QSUMSQ)
#else
     QSUMSQ = QSUMSQ + MM(J,J) * MM(J,J)
#endif
     DO I = J+1, N
#ifdef USE_FMA
        QTMP0 = C_FMA(MM(I,J), MM(I,J), QTMP0)
#else
        QTMP0 = QTMP0 + MM(I,J) * MM(I,J)
#endif
     END DO
  END DO
#ifdef USE_FMA
  QSUMSQ = C_FMA(TWO, QTMP0, QSUMSQ)
#else
  QSUMSQ = TWO * QTMP0 + QSUMSQ
#endif

  NRMS(8) = DBLE(SQRT(QSUMSQ))
  IF (P .EQ. N) THEN
     WRITE (*,1) '|| V V^T - I ||_F', NRMS(8)
  ELSE
     WRITE (*,1) '|| V J V^T - J ||_F', NRMS(8)
  END IF
#endif

  CALL h5fopen_f(FH5, H5F_ACC_RDONLY_F, FID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening FH5!'

  CALL h5gopen_f(FID, GH5, GID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening FH5/GH5!'

  ALLOCATE(DTMP2(M,M))
  DIMS2(1) = M
  DIMS2(2) = M
  CALL h5ltread_dataset_double_f(GID, 'A', DTMP2, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading A!'
  !$OMP PARALLEL DO PRIVATE(I)
  DO J = 1, M
     DO I = 1, M
        MM(I,J) = REAL(DTMP2(I,J),10)
     END DO
  END DO
  DEALLOCATE(DTMP2)

  ALLOCATE(IPIV(N))
  DIMS1(1) = N
  CALL h5ltread_dataset_int_f(GID, 'IPIV', IPIV, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading IPIV!'
  CALL QBAKCP(M, N, U, M, IPIV)
  DEALLOCATE(IPIV)

  !$OMP PARALLEL DO PRIVATE(I)
  DO J = 1, N
     DO I = 1, M
        V(I,J) = U(I,J) * D(J)
     END DO
  END DO

  !$OMP PARALLEL DO PRIVATE(I,K)
  DO J = 1, N
     DO I = 1, M
        DO K = 1, M
#ifdef USE_FMA
           V(I,J) = C_FMA(-MM(I,K), U(K,J), V(I,J))
#else
           V(I,J) = V(I,J) - MM(I,K) * U(K,J)
#endif
        END DO
        V(I,J) = ABS(V(I,J))
     END DO
  END DO

  DEALLOCATE(MM)
  DEALLOCATE(U)

  QMAXA = ZERO
  DO J = 1, N
     DO I = 1, M
        IF (V(I,J) .GT. QMAXA) QMAXA = V(I,J)
     END DO
  END DO

  NRMS(9) = DBLE(QMAXA)
  WRITE (*,1) '|| A U - U D ||_M', NRMS(9)

  QSUMSQ = ZERO
  DO J = 1, N
     DO I = 1, M
#ifdef USE_FMA
        QSUMSQ = C_FMA(V(I,J), V(I,J), QSUMSQ)
#else
        QSUMSQ = QSUMSQ + V(I,J) * V(I,J)
#endif
     END DO
  END DO

  NRMS(10) = DBLE(SQRT(QSUMSQ))
  WRITE (*,1) '|| A U - U D ||_F', NRMS(10)

  !$OMP PARALLEL DO PRIVATE(I,QTMP0)
  DO J = 1, N
#ifdef MAGMA
     QTMP0 = SQRT(D(J))
#else
     QTMP0 = SQRT(ABS(D(J)))
#endif
     DO I = 1, N
#ifdef MAGMA
        V(I,J) = V(I,J) / (SQRT(D(I)) * QTMP0)
#else
        V(I,J) = V(I,J) / (SQRT(ABS(D(I))) * QTMP0)
#endif
     END DO
  END DO

  QMAXA = ZERO
  DO J = 1, N
     DO I = 1, M
        IF (V(I,J) .GT. QMAXA) QMAXA = V(I,J)
     END DO
  END DO

  NRMS(11) = DBLE(QMAXA)
  WRITE (*,1) '|| A U - U D ||_R', NRMS(11)

  DEALLOCATE(V)

  ALLOCATE(LAMBDA(N))

  DIMS1(1) = N
  CALL h5ltread_dataset_double_f(GID, 'LAMBDA', LAMBDA, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading LAMBDA!'

  CALL h5gclose_f(GID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing FH5/GH5!'

  CALL h5fclose_f(FID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing FH5!'

  CALL DLASRT('D', N, DD, INFO)
  IF (INFO .NE. 0) STOP 'Error sorting D!'

  CALL h5fopen_f(RH5, H5F_ACC_RDWR_F, FID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening RH5 for writing!'

  CALL h5gopen_f(FID, GH5, GID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening RH5/GH5 for writing!'

  DIMS1(1) = N
  INFO = h5ltfind_dataset_f(GID, 'LAMBDA')
  IF (INFO .EQ. 0) THEN
     CALL h5ltmake_dataset_double_f(GID, 'LAMBDA', 1, DIMS1, LAMBDA, INFO)
     IF (INFO .NE. 0) STOP 'Error writing LAMBDA!'
  ELSE
     CALL h5dopen_f(GID, 'LAMBDA', DID, INFO)
     IF (INFO .NE. 0) STOP 'Error opening LAMBDA!'
     CALL h5dwrite_f(DID, H5T_NATIVE_DOUBLE, LAMBDA, DIMS1, INFO)
     IF (INFO .NE. 0) STOP 'Error writing to LAMBDA!'
     CALL h5dclose_f(DID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing LAMBDA!'
  END IF

  DIMS1(1) = INRMSL
  INFO = h5ltfind_dataset_f(GID, 'NRMS')
  IF (INFO .EQ. 0) THEN
     CALL h5ltmake_dataset_double_f(GID, 'NRMS', 1, DIMS1, NRMS, INFO)
     IF (INFO .NE. 0) STOP 'Error writing NRMS!'
  ELSE
     CALL h5dopen_f(GID, 'NRMS', DID, INFO)
     IF (INFO .NE. 0) STOP 'Error opening NRMS!'
     CALL h5dwrite_f(DID, H5T_NATIVE_DOUBLE, NRMS, DIMS1, INFO)
     IF (INFO .NE. 0) STOP 'Error writing to NRMS!'
     CALL h5dclose_f(DID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing NRMS!'
  END IF

  CALL DLASRT('D', N, LAMBDA, INFO)
  IF (INFO .NE. 0) STOP 'Error sorting LAMBDA!'

  ALLOCATE(ABSERR(N+1))

  SCL = ZERO
  ! Not so absoulte, keep error directions.
  DO I = 1, N
     ABSERR(I) = LAMBDA(I) - DD(I)
     IF (ABS(SCL) .LT. ABS(ABSERR(I))) SCL = ABSERR(I)
  END DO
  ABSERR(N+1) = SCL
  WRITE (*,1) 'max abs err', SCL

  DIMS1(1) = N+1
  INFO = h5ltfind_dataset_f(GID, 'ABSERR')
  IF (INFO .EQ. 0) THEN
     CALL h5ltmake_dataset_double_f(GID, 'ABSERR', 1, DIMS1, ABSERR, INFO)
     IF (INFO .NE. 0) STOP 'Error writing ABSERR!'
  ELSE
     CALL h5dopen_f(GID, 'ABSERR', DID, INFO)
     IF (INFO .NE. 0) STOP 'Error opening ABSERR!'
     CALL h5dwrite_f(DID, H5T_NATIVE_DOUBLE, ABSERR, DIMS1, INFO)
     IF (INFO .NE. 0) STOP 'Error writing to ABSERR!'
     CALL h5dclose_f(DID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing ABSERR!'
  END IF

  ALLOCATE(RELERR(N+1))

  SCL = ZERO
  DO I = 1, N
     RELERR(I) = ABSERR(I) / LAMBDA(I)
     IF (ABS(SCL) .LT. ABS(RELERR(I))) SCL = RELERR(I)
  END DO
  RELERR(N+1) = SCL
  WRITE (*,1) 'max rel err', SCL

  DIMS1(1) = N+1
  INFO = h5ltfind_dataset_f(GID, 'RELERR')
  IF (INFO .EQ. 0) THEN
     CALL h5ltmake_dataset_double_f(GID, 'RELERR', 1, DIMS1, RELERR, INFO)
     IF (INFO .NE. 0) STOP 'Error writing RELERR!'
  ELSE
     CALL h5dopen_f(GID, 'RELERR', DID, INFO)
     IF (INFO .NE. 0) STOP 'Error opening RELERR!'
     CALL h5dwrite_f(DID, H5T_NATIVE_DOUBLE, RELERR, DIMS1, INFO)
     IF (INFO .NE. 0) STOP 'Error writing to RELERR!'
     CALL h5dclose_f(DID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing RELERR!'
  END IF

  CALL h5gclose_f(GID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing RH5/GH5 for writing!'

  CALL h5fclose_f(FID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing RH5 for writing!'

  CALL h5close_f(INFO)
  IF (INFO .NE. 0) STOP 'Error closing HDF5!'

  DEALLOCATE(RELERR)
  DEALLOCATE(ABSERR)
  DEALLOCATE(LAMBDA)
  DEALLOCATE(DD)
  DEALLOCATE(D)

1 FORMAT(A,' = ',ES25.17E3)

CONTAINS

  SUBROUTINE QBAKCP(M, N, V, LDV, IPIV)
    ! Back-permutes the rows of V according to information stored in IPIV.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: LDV

    INTEGER, INTENT(IN) :: IPIV(N)
    REAL(10), INTENT(INOUT) :: V(LDV,N)

    REAL(10) :: QTMP
    INTEGER :: K
    INTEGER :: KP
    INTEGER :: J

    DO K = N, 1, -1
       KP = ABS(IPIV(K))
       IF (KP .NE. K) THEN
          DO J = 1, N
             QTMP = V(K,J)
             V(K,J) = V(KP,J)
             V(KP,J) = QTMP
          END DO
       END IF
    END DO
  END SUBROUTINE QBAKCP

END PROGRAM PP80
