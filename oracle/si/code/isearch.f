
      INTEGER*4 FUNCTION ISEARCH(X,A,N)
      REAL A(N)
C
C-----INDEX SEARCH
      IF (X.LE.A(1)) THEN
      ISEARCH=1
      ELSE IF (X.GE.A(N)) THEN
      ISEARCH=N-1
      ELSE
      K=1
      J=N
   10 IF (K.EQ.J-1) GO TO 20
      I=(K+J)/2
      IF (X.LE.A(I)) THEN
      J=I
      ELSE
      K=I
      END IF
      GO TO 10
   20 ISEARCH=K
      END IF
C
      RETURN
      END
