PROGRAM taskA
IMPLICIT NONE
INCLUDE 'mkl_lapack.fi'
!
DOUBLE PRECISION:: A,B,a1,b1,c1,d1
INTEGER :: N
INTEGER :: K,I,J,T,C,Q,R
DOUBLE PRECISION :: L,h,S
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: D,X
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: E
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: V,v1
!
CHARACTER(LEN=*),PARAMETER :: format='(E15.5,200(",",E15.5))'
CHARACTER,PARAMETER :: JOBZ = 'V',RANGE = 'I'


DOUBLE PRECISION :: VL = 0,VU = 20
INTEGER :: IL,IU,UB
DOUBLE PRECISION :: ABSTOL
INTEGER :: M, INFO, LDZ
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: W
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: Z
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: WORK
INTEGER,DIMENSION(:),ALLOCATABLE :: IWORK
INTEGER,DIMENSION(:),ALLOCATABLE :: IFAIL

CHARACTER(LEN=*), PARAMETER :: inputA = 'inputA.txt'
OPEN(UNIT=7, FILE=inputA)
    READ(7,*) A,B,N,IL,IU,a1,b1,c1,d1
    CLOSE(7)
!
!
! definizione dei vettori che verranno utilizzati nella subroutine.
!
!
ALLOCATE(D(1:N))
ALLOCATE(X(1:N))
ALLOCATE(E(1:N-1))
!
L = B-A        !larghezza della doppia buca di potenziale
h = L/dble(N-1)        !distanza tra i punti della griglia
S = A
X(1) = A
X(N) = B
DO K=2,N-1        !creazione della griglia 
    S = S+h
    X(K) = S
END DO
DO I=1,N        !creazione della diagonale principale della matrice
    D(I) = 1/h**2+(a1 + b1*X(I) + c1*(X(I)**2) + d1*(X(I)**4))/2
END DO
    E = (/(-dble(0.5)/h**2,J=1,N-1)/)        !diagonale secondaria della matrice
!
!
! Utilizzo della subroutine per il calcolo di autovalori e autovettori.
!
!
M=IU-IL+1
LDZ = N
UB = IU + 5     !Per avere un upper bound da definire in precedenza
!ALLOCATE(Z(LDZ,1:UB))
ALLOCATE(Z(LDZ,max(1,M)))
ALLOCATE(W(1:N))
ALLOCATE(WORK(1:5*N))
ALLOCATE(IWORK(1:5*N))
ALLOCATE(IFAIL(1:N))
ABSTOL = 2 * DLAMCH('S')
CALL dstevx(JOBZ,RANGE,N,D,E,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,IWORK,IFAIL,INFO)
ALLOCATE(V(1:M))
V = (/(W(C),C=1,M)/)        !autovalori
IF(M .NE. (IU-IL+1)) THEN
    PRINT *, "Il numero di valori trovati non corrisponde a quello cercato."
END IF
!
!
!   i file di output del programma.
!
!
OPEN(23,FILE='autovalori.out')
WRITE(23,*) V
CLOSE(23)

!
OPEN(27,FILE='autovettoritaskA.out')

DO Q=1,LDZ
WRITE(27,FMT=format) (Z(Q,R), R=1,M)
end do
CLOSE(27)
END PROGRAM taskA
