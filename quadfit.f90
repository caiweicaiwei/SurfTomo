!***********************************************************************************************************************************
!  QUADFIT
!
!  Inputs:
!    N      Array size
!    XARR   1-D array of X values (N elements)
!    YARR   1-D array of Y values (N elements)
!  Outputs:
!    A, B, C    Quadratic coefficients of equation Y = A X^2 + B X + C
!
!  Ref:  Jean Meeus, Astronomial Algorithms (2nd ed.), Chapter 4.
!***********************************************************************************************************************************

      SUBROUTINE QUADFIT (N, XARR, YARR, A, B, C)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N) :: XARR, YARR
      REAL, INTENT(OUT) :: A, B, C
      REAL :: P, Q, R, S, T, U, V, D


      P = SUM(XARR)
      Q = SUM(XARR**2)
      R = SUM(XARR**3)
      S = SUM(XARR**4)
      T = SUM(YARR)
      U = SUM(XARR*YARR)
      V = SUM(XARR**2*YARR)

      D = N*Q*S + 2.0D0*P*Q*R - Q**3 - P**2*S - N*R**2

      A = (N*Q*V + P*R*T + P*Q*U - Q**2*T - P**2*V - N*R*U)/D
      B = (N*S*U + P*Q*V + Q*R*T - Q**2*U - P*S*T - N*R*V)/D
      C = (Q*S*T + Q*R*U + P*R*V - Q**2*V - P*S*U - R**2*T)/D

      RETURN

      END SUBROUTINE QUADFIT
