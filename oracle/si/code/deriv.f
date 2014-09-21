

      SUBROUTINE DERIV(nmax,nstart,x,f,d,mode)   
      INCLUDE   'lf_param.inc'                             ! nxdepth,taux_lam,coefj(k,l)
      INTEGER*4 nmax,nstart,mode
      REAL*4 x(nmax),f(nmax),d(nmax),y(NDMAX)
      CALL DERIV1(x,f,d,nmax,nstart)
      IF (mode.EQ.2) THEN
        DO i = 1,nmax
           y(i) = d(i)
        ENDDO
        CALL DERIV1(x,y,d,nmax,nstart)
      ENDIF
      RETURN
      END
