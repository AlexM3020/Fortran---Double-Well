PROGRAM oscarmB
   USE fftw
    USE modulozheevx
    IMPLICIT NONE
!
    DOUBLE PRECISION :: A,B,L,deltax,deltaq,O,h
    INTEGER :: IL,IU,N,K,i,T,C,S,P,g
    DOUBLE PRECISION, PARAMETER :: pi=ACOS(-1.0)
    DOUBLE COMPLEX,DIMENSION(:), ALLOCATABLE :: V,V_T
    DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: X,Q,W
    DOUBLE COMPLEX,DIMENSION(:,:), ALLOCATABLE :: MATRIXA
    CHARACTER(LEN=*), PARAMETER :: inputB1 = 'inputB1.txt'
    CHARACTER(LEN=*),PARAMETER :: format='(E15.5,200(",",E15.5))'
  
    OPEN(UNIT=35, FILE=inputB1) ! Lettura del file di input.


    READ(35,*) A,B,N,IL,IU,K
    CLOSE(35)
        ALLOCATE(X(1:N+1))
        ALLOCATE(Q(1:N+1))
        ALLOCATE(V(1:N+1))
        ALLOCATE(V_T(1:N+1))
      
!
! Definizione dei vettori utili e del potenziale da utilizzare.
!
        L = B-A
        deltax = L/dble(N)
        deltaq = 2*pi/L
        X(1) = -L/2
        Q(1) = -(pi*N)/L
        DO i = 1,N
            X(i+1) = X(i) + deltax
        END DO
        DO g = 1,N
            Q(g+1) = Q(g) + deltaq
        END DO
        V = X**2

        CALL fttw_(L,N,deltax,V,V_T)
        OPEN(36,FILE='V_T.txt')
        OPEN(37,FILE='DeltaK.txt')
        DO P=1,N+1
            WRITE(36,*) DREAL(V_T(P))    
            WRITE(37,*) Q(P)
        END DO
        CLOSE(36)
        CLOSE(37)


! Costruzione della matrice da diagonalizzare
        ALLOCATE(MATRIXA(1:2*K+1,1:2*K+1)) 
        MATRIXA = 0.0
        DO S=-K,K
            DO C=S,K
                MATRIXA(K+S+1,K+C+1) = V_T(N/2+1+S-C)*0.5
            END DO
            MATRIXA(K+S+1,K+S+1) = MATRIXA(K+S+1,K+S+1) + 0.5*(S*2*pi/L)**2
        END DO

        
       CALL zheevx_(MATRIXA,2*K+1,IU,IL)

   
END PROGRAM oscarmB

        

            
