

!     JKR, 09.07.95
      SUBROUTINE resample(mode,rdst,wo,fo,so,no,go,w,f,g,s,n,k,xdst)
!     
      INCLUDE 'physcnst.inc'
      INTEGER*4 no,n,k  ! input number, output number, actual number
      INTEGER*4 ntmp
      INTEGER*4 mode    ! 0: exact rdst , 1: rdst is min. dw, 
!                         2: next power of 2 
      REAL*4    fo(no),so(no),go(no),f(n),s(n),g(n)   ! flux, sigma
      REAL*8    rdst    ! resmpling distance
      REAL*8    wo(no),w(n) ! wavelength
      REAL*8    xdst    ! new resampling distance
      REAL*8    dwmin   ! temporary min. resampling distance 

!     determination of xdst and establish new w scale
    
      xdst = 0.D0
      k    = INT( (wo(n)-wo(1)) / rdst + 0.5 )
      xdst = rdst
      IF     (mode.GE.1) THEN
        dwmin = 1.D36
        DO i=2,no
           IF ( (wo(i) - wo(i-1)) .LT. dwmin ) dwmin = wo(i)-wo(i-1)
        END DO
        IF (mode.EQ.2) THEN
         xdst = DMAX1(dwmin*0.9D0,rdst)
         k = 2**( 1 + INT(log10((wo(n)-wo(1))/xdst) / LOG2 ))
         xdst =  ( wo(n)-wo(1) ) / (DBLE(k - 1)) 
        ELSE
          xdst = DMAX1(dwmin*0.7D0,rdst)
          k    = INT( (wo(n)-wo(1)) / xdst + 0.5 )
        ENDIF
      ENDIF
     
!     Interpolation

      w(1) = wo(1)
      f(1) = fo(1)
      s(1) = so(1)
      g(1) = go(1)
      ntmp = 1
      DO i = 2,k
         w(i) = w(i-1) + xdst
         f(i) = QPOL1(w(i),wo,fo,no,ntmp)
         s(i) = QPOL1(w(i),wo,so,no,ntmp)
         g(i) = QPOL1(w(i),wo,go,no,ntmp)
      ENDDO

      RETURN
      END
