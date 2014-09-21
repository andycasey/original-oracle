
      SUBROUTINE DSPLINE(X,Y,A,B,C,N)                ! Double Precision
      IMPLICIT REAL*8 (A-H),REAL*8 (O-Z)
      DIMENSION X(N),Y(N),A(N),B(N),C(N)
C
C-----COMPUTATION OF CUBIC SPLINE COEFFICIENTS A,B,C
C-----FOR APPROXIMATION OF A DATA ARRAY (X,Y) SUCH THAT
C-----F(X)=Y(I)+DEL*[A(I)+DEL*(B(I)+C(I)*DEL)] IN THE
C-----INTERVAL X(I)<X<X(I+1), WHERE DEL=X-X(I).
C-----VANISHING CURVATURE AT THE BOUNDARIES IS ASSUMED:
C-----F"[X(1)]=F"[X(N)]=0
      X1=X(1)
      DO 10 I=2,N
      Z1=X(I)
      X(I)=Z1-X1
   10 X1=Z1
C
      DO 20 I=2,N-1
      A(I)=X(I+1)/(X(I)+X(I+1))
      B(I)=1.-A(I)
      C(I)=6.*((Y(I+1)-Y(I))/X(I+1)-(Y(I)-Y(I-1))/X(I))
     *  /(X(I)+X(I+1))
   20 CONTINUE
C
      B(2)=2.
      DO 30 I=3,N-1
      X1=B(I)/B(I-1)
      B(I)=2.-X1*A(I-1)
   30 C(I)=C(I)-X1*C(I-1)
      C(N-1)=C(N-1)/B(N-1)
C
      DO 40 I=N-2,2,-1
   40 C(I)=(C(I)-A(I)*C(I+1))/B(I)
      C(1)=.0
      C(N)=.0
C
      DO 50 I=1,N-1
      A(I)=(Y(I+1)-Y(I))/X(I+1)-X(I+1)*(2.*C(I)+C(I+1))/6.
      B(I)=.5*C(I)
   50 C(I)=(C(I+1)-C(I))/X(I+1)/6.
      DO 60 I=2,N
   60 X(I)=X(I-1)+X(I)
C
      RETURN
      END
