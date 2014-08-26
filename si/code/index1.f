
      SUBROUTINE INDEX1(XVAR,X,NDIMX,K,XMIN)
      IMPLICIT REAL*4 (A-H,O-Z)
      REAL*4 X(NDIMX) 
      IF ( XVAR .LT. XMIN ) GOTO 300
 100  K = K + 1
      IF ( ( K .GE. NDIMX ) .OR. ( XVAR .LE. X(K) ) ) GOTO 200
      GOTO 100
 200  K = K - 1
      GOTO 500
 300  K = K + 1
 400  K = K - 1
      IF ( ( K .LE. 1 ) .OR. ( XVAR .GE. X(K) ) ) GOTO 500
      GOTO 400
 500  CONTINUE
      RETURN
      END 
