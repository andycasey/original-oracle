
      SUBROUTINE SPLINE_CENTER(IBUF,LONG,NN,INCR,CEN,EXTR,NUM,*)
      INTEGER*2 IBUF(1),IMAX,IA,IB,IC,ID,IS
C
C-----LOCAL EXTREMUM OF IBUF COMPUTED FROM CUBIC SPLINES:
C-----ON ENTRY N DEFINES AN INITIAL GUESS OF POSITION,
C-----ON EXIT N CONTAINS THE NEXT LOWER PIXEL POSITION.
C-----FOR MULTIDIMENSIONAL CENTERING, INCR DEFINES THE
C-----DIRECTION OF THE SEARCH.
C-----CEN OBTAINS THE FRACTIONAL DISTANCE FROM N, AND
C-----EXTR IS THE INTERPOLATED MAXIMUM (NUM>0) OR MINIMUM
C-----(NUM<0). ABS(NUM) DETERMINES THE INCREMENTAL EXTENSION
C-----OF THE SEARCH.
      IS=ISIGN(1,NUM)
      IMAX=-30000
      NDIF=IABS(NUM)*INCR
      N=NN
C
   10 I1=N-NDIF
      I2=N+NDIF
      IF (I1.LT.1 .OR. I2.GT.LONG) RETURN 1
      DO 20 I=I1,I2,INCR
      IB=IS*IBUF(I)
      IF (IB.GT.IMAX) THEN
      N=I
      IMAX=IB
      END IF
   20 CONTINUE
      IF (N.EQ.I1 .OR. N.EQ.I2) RETURN 1
C
      IF (ISIGN(1,IBUF(N-INCR)-IBUF(N+INCR)).EQ.IS) N=N-INCR
      IA=IBUF(N-INCR)
      IB=IBUF(N)
      IC=IBUF(N+INCR)
      ID=IBUF(N+INCR+INCR)
      A=IS*(12.*IC-2.*ID-3.*IB-7.*IA)/15.
      B=IS*(6.*IC-ID-9.*IB+4.*IA)/5.
      C=IS*(IB-IC+(ID-IA)/3.)
      IF (C.EQ..0) THEN
      CEN=-.5*A/B
      ELSE IF (B.EQ..0) THEN
      CEN=SQRT(-A/C/3.)
      ELSE
      D=-B/C/3.
      E=SQRT(D*D-A/C/3.)
      CEN=D+E
      IF (CEN.GT.1.) CEN=D-E
      END IF
      EXTR=IB+IS*CEN*(A+CEN*(B+CEN*C))
      NN=N
C
      RETURN
      END
