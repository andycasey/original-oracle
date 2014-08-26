C --- ------------------------- CHECK_STEPWIDTH_CRIT ---------------------------
      LOGICAL*4 FUNCTION CHECK_STEPWIDTH_CRIT(wsp,fsp,sw_crit)
C     Stepwidth control checking curvature of profile
      REAL*8    as(4),bs(4),cs(4)               ! Buffer for Spline Coefficients
      REAL*8    fsp(4)                          ! Fluxpoints
      REAL*8    wsp(4)                          ! wavelengthpoints
      REAL*8    z,dz
      REAL*8    f1,f2                           ! 1. and 2. derivation
      REAL*8    sw_crit                         ! stepwidth-criterium
      REAL*8    xbuf
      REAL*8    hex

      CALL DSPLINE(wsp,fsp,as,bs,cs,4)

c      dz = (wsp(4)-wsp(3))*0.5
c      f1 = (as(3) + dz*(2.D0*bs(3) + 3.D0*cs(3)*dz)) * 1.D-8
c      f2 = (2.D0*bs(3) + 6.D0*cs(3) * dz) * 1.D-16

      f1 = as(3) * 1.D-8
      f2 = bs(3) * 2.D-16
      hex= DABS(f2/dsqrt(1.D0+f1*f1)**3.D0)
      xbuf = DABS(f1) + hex
C      write(*,100) wsp(3)*1.D8,fsp(3),f1,f2,xbuf,hex
C100   FORMAT(' ',f8.3,' ',f7.4,' ',4(G13.4,' '))
      CHECK_STEPWIDTH_CRIT = sw_crit.LT.xbuf

      RETURN
      END
