

      FUNCTION SPOL(Z,X,Y,A,B,C,N)
      REAL*4 X(N),Y(N),A(N),B(N),C(N)
      INTEGER*4 N
      INTEGER*4 ISEARCH
C
C-----CUBIC SPLINE INTERPOLATION
      K=ISEARCH(Z,X,N)
      DEL=Z-X(K)
      SPOL=Y(K)+DEL*(A(K)+DEL*(B(K)+DEL*C(K)))
      RETURN
      END
