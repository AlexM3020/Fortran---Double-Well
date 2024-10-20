MODULE fftw
    IMPLICIT NONE
!
    CONTAINS
!
        SUBROUTINE fttw_(L,N,deltax,A,ff3)
            IMPLICIT NONE
            INCLUDE 'fftw3.f'
            DOUBLE COMPLEX, DIMENSION(1:N+1), INTENT(IN) :: A
            DOUBLE COMPLEX, DIMENSION(1:N+1) :: ff1,ff2,ff4,ff5,ff6
            DOUBLE COMPLEX, DIMENSION(1:N+1), INTENT(OUT) :: ff3
          
            INTEGER(KIND=8) :: plan1,plan2
            DOUBLE PRECISION, INTENT(IN) :: L,deltax
            INTEGER, INTENT(IN) :: N
!
                ff1(N/2+2:) = A(1:N/2)
                ff1(1:N/2+1) = A(N/2+1:)
!
                CALL dfftw_plan_dft_1d(plan1,N+1,ff1,ff2,FFTW_forward,fftw_estimate)
                CALL dfftw_execute(plan1,ff1,ff2)
                CALL dfftw_destroy_plan(plan1)
                ff2 = ff2/dble(N+1)
                ff3(N/2+1:) = ff2(1:N/2+1)
                ff3(1:N/2) = ff2(N/2+2:)
!
                ff4(N/2+2:) = ff3(1:N/2)
                ff4(1:N/2+1) = ff3(N/2+1:)

                CALL dfftw_plan_dft_1d(plan2,N+1,ff4,ff5,FFTW_backward,fftw_estimate)
                CALL dfftw_execute(plan2,ff4,ff5)
                CALL dfftw_destroy_plan(plan2)
                ff6(N/2+1:) = ff5(1:N/2+1)
                ff6(1:N/2) = ff5(N/2+2:)
              
                   
        END SUBROUTINE fttw_
!
END MODULE fftw
