

      SUBROUTINE EQUI_SPLINE(IBUF,A,B,C,N)
      INTEGER*2 IBUF(1)
      REAL A(N),B(N),C(N)
C
C-----CUBIC SPLINE COEFFICIENTS FOR EQUIDISTANT INTERVALS
      DO 10 I=2,N-1
      B(I)=.5
   10 C(I)=3.*(IBUF(I+1)-2.*IBUF(I)+IBUF(I-1))
      B(2)=2.
C
      DO 20 I=3,N-1
      X=B(I)/B(I-1)
      B(I)=2.-.5*X
   20 C(I)=C(I)-X*C(I-1)
      C(N-1)=C(N-1)/B(N-1)
      DO 30 I=N-2,2,-1
   30 C(I)=(C(I)-.5*C(I+1))/B(I)
      C(1)=.0
      C(N)=.0
C
      DO 40 I=1,N-1
      A(I)=IBUF(I+1)-IBUF(I)-(2.*C(I)+C(I+1))/6.
      B(I)=.5*C(I)
   40 C(I)=(C(I+1)-C(I))/6.
C
      RETURN
      END
