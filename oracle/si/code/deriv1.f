
      SUBROUTINE DERIV1(X,F,DFDX,N,NS)
      IMPLICIT REAL*4 (A-H,O-Z)
C     ASSUMES THAT ANY ZERO IN X OCCURS AT A ENDPOINT
      DIMENSION X(N),F(N),DFDX(N)
      DFDX(1)=(F(2)-F(1))/(X(2)-X(1))
      N1=N-1
      DFDX(N)=(F(N)-F(N1))/(X(N)-X(N1))
      IF(N.EQ.2)RETURN
      S=ABS(X(2)-X(1))/(X(2)-X(1))
      DO 1 J=MAX(2,NS),N1
       SCALE=AMAX1(ABS(F(J-1)),ABS(F(J)),ABS(F(J+1)))/ABS(X(J))
       IF(SCALE.EQ.0.)SCALE=1.
       D1=(F(J+1)-F(J))/(X(J+1)-X(J))/SCALE
       D=(F(J)-F(J-1))/(X(J)-X(J-1))/SCALE
       TAN1=D1/(S*SQRT(1.+D1**2)+1.)
       TAN=D/(S*SQRT(1.+D**2)+1.)
    1 DFDX(J)=(TAN1+TAN)/(1.-TAN1*TAN)*SCALE
      RETURN
      END
