MODULE modulozheevx
    IMPLICIT NONE
!
    CONTAINS
        SUBROUTINE zheevx_(MATRIXA,N,IU,IL)
        IMPLICIT NONE
        INCLUDE 'mkl_lapack.fi'
!        
        INTEGER,INTENT(IN) :: N,IU,IL   
        DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: MATRIXA
        DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: EV
        INTEGER :: Y,R      
!    
        CHARACTER,PARAMETER :: JOBZ = 'V',RANGE = 'I',UPLO = 'U'
        DOUBLE PRECISION :: VL = 0,VU = 10
        INTEGER :: UB
        DOUBLE PRECISION :: ABSTOL
        INTEGER :: INFO, LDZ, LDA, LWORK, M
        DOUBLE PRECISION,DIMENSION(1:N) :: W
        DOUBLE COMPLEX,DIMENSION(:,:),ALLOCATABLE :: Z
        DOUBLE COMPLEX,DIMENSION(1:2*N) :: WORK
        INTEGER,DIMENSION(1:5*N) :: IWORK
        DOUBLE PRECISION,DIMENSION(1:7*N) :: RWORK
        INTEGER,DIMENSION(1:N) :: IFAIL
!
        CHARACTER(LEN=*),PARAMETER :: format='(E15.5,200(",",E15.5))'
!
! Utilizzo della subroutine Zheevx 
!
        LWORK = 2*N
        LDA = N
        LDZ = N
        UB = IU + 5     !upper bound
        ABSTOL = 2 * DLAMCH('S')
        ALLOCATE(Z(LDZ,1:UB))
        CALL zheevx(JOBZ,RANGE,UPLO,N,MATRIXA,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

        ALLOCATE(EV(1:M))
        EV = W(1:M)
        
!
! Scrittura dei file di output.
!
        OPEN(23,FILE='autovalori.out')
        WRITE(23,*) EV
        CLOSE(23)
!       

        OPEN(27,FILE='autovettori.out')
        DO Y=1,LDZ
            WRITE(27,FMT=format) (REAL(Z(Y,R)), R=1,M)
        END DO
        CLOSE(27)
!
        OPEN(28,FILE='imm.out')
        DO Y=1,LDZ
            WRITE(28,FMT=format) (AIMAG(Z(Y,R)), R=1,M)
        END DO
        CLOSE(28)
!
        END SUBROUTINE zheevx_
!
END MODULE modulozheevx
