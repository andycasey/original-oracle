
!     Convolution with Intrumental, Rotationprofile and/or Makroturbulence 

      SUBROUTINE convol(w,f,n,px,py,answ,fft,
     *                  vdop,gamma,expo,vsini,beta)
      INCLUDE 'physcnst.inc'
         
      REAL*4 vdop    ! Makroturbulence in cm/s
      REAL*4 expo    ! Makroturbulence in cm/s
      REAL*4 gamma   ! Gamma (Lorentzian) in cm
      REAL*4 vsini   ! Rotationprofile in cm/s
      REAL*4 beta    ! Limbdarkening Parameter (default is 1.5)

      REAL*8    xcen,x,x2
      REAL*8    w(n),px(n) , dw, a, width , fmean
      REAL*4    f(n),py(n) , offs
      REAL*4    answ(2*n)
      COMPLEX*8 fft(n)

      xcen    = (w(n)+w(1)) * .5D0 
      dw      = (w(n)-w(1)) / DBLE(n-1)


      px(1)   = 0.
      px(n)   = -dw
      DO i  = 1,n/2-1
         px(n-i) = px(n-i+1) - dw
         px(i+1) = px(i)     + dw
      END DO
      px(n/2+1)  = px(n/2) + dw

      offs  = 0.5*(f(1)+f(n))
      DO i=1,n
        f(i) = f(i) - offs
      ENDDO

      fmean = DBLE(f(1) + f(n))*0.5D0
      DO i=2,n-1
         fmean = fmean + DBLE(f(i))
      ENDDO
      fmean = fmean / DBLE(n-1)


! -------------------------------- rotation ----------------------------
      IF (vsini.GT.0.) THEN            
        IF (beta.LE.0.) beta = 1.5
        xnorm = C/(DBLE(vsini)*xcen)
        DO i    = 1,n/2+1
          x     = xnorm * px(i) 
          IF(x .LT. 1.) THEN
            x2    = x*x
            py(i) =REAL(dw*(2.D0/PI*SQRT(1.D0-x2)+(1.D0-x2) 
     *                 * beta * .5D0) * xnorm  / 
     *                  (1.D0+ 0.66666666667D0 * DBLE(beta)))
            IF (i.GT.1) py(n-i+2) = py(i)
          ELSE
            py(i)     = 0.
            py(n-i+1) = 0.
          ENDIF
        END DO

        CALL convlv(f,py,n,+1,answ,fft)
        DO i=1,n
          f(i) = answ(i)  + REAL(fmean) 
        ENDDO
      ENDIF
! -------------------------------- vdop ----------------------------
      IF (vdop.GT.0.) THEN            
        width = xcen * DBLE(vdop) / C
        IF (gamma.GT.0.) THEN                     ! and gamma = VOIGT
          DO i    = 1,n/2+1
            a     = gamma/ width
            py(i) = REAL(dw*VOIGT(a,px(i)/width) / (SQRT_PI*width))
            IF (i.GT.1) py(n-i+2) = py(i)
          END DO
        ELSE    
          DO i    = 1,n/2+1
            py(i) = REAL(dw*DEXP(-(px(i)/width)**2) / (SQRT_PI*width))
            IF (i.GT.1) py(n-i+2) = py(i)
          END DO
        ENDIF
        CALL convlv(f,py,n,+1,answ,fft)
        DO i=1,n
          f(i) = answ(i) + REAL(fmean) 
        ENDDO
      ENDIF

! -------------------------------- exponential ----------------------------
      IF (expo.GT.0.) THEN            
        width = DBLE(expo) / C  *  xcen 
        DO i  = 1,n/2+1
          py(i) = REAL(dw*DEXP( -px(i) / width)  / (2.D0 * width))
          IF (i.GT.1) py(n-i+2) = py(i)
        END DO

        CALL convlv(f,py,n,+1,answ,fft)

        DO i=1,n
          f(i) = answ(i) + REAL(fmean) 
        ENDDO
      ENDIF
! -------------------------------------------------------------------------

      DO i=1,n
        f(i) = f(i) + offs
      ENDDO

      RETURN
      END
