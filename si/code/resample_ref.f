
!     JKR, 09.07.95
      SUBROUTINE resample_ref(mode,rdst,ref,nref,k,xdst)
!     
      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_ref.inc'
      RECORD /ref_type/ ref(nref)      
      INTEGER*4 n,k
      INTEGER*4 ntmp
      INTEGER*4 mode    ! 0: exact rdst , 1: rdst is min. dw, 
!                         2: next power of 2 
      REAL*8    rdst          ! resmpling distance
      REAL*8    w(MAX_REF_NR) ! wavelength
      REAL*4    y(MAX_REF_NR) ! buffer
      REAL*8    xdst          ! new resampling distance
      REAL*8    dwmin         ! temporary min. resampling distance 
!     determination of xdst and establish new w scale
    
      DO i = 1,nref
         w(i) = ref(i).w
         y(i) = ref(i).f
      END DO

      xdst = 0.D0
      k    = INT( (w(nref)-w(1)) / rdst + 0.5 )
      xdst = rdst
      IF     (mode.GE.1) THEN
        dwmin = 1.D36
        DO i=2,nref
           IF (DABS(w(i) - w(i-1)).LT.dwmin) dwmin = DABS(w(i)-w(i-1))
        END DO
        IF (mode.EQ.2) THEN
          xdst = DMAX1(dwmin*0.9D0,rdst)
          k = 2**( 1 + INT(LOG10((w(nref)-w(1))/xdst) / LOG2 ))
          xdst =  ( w(nref)-w(1) ) / (DBLE(k - 1)) 
        ELSE
          xdst = DMAX1(dwmin*0.7D0,rdst)
          k    = INT( (w(nref)-w(1)) / xdst + 0.5 )
        ENDIF
      ENDIF
     
!     Interpolation

      ntmp = 1
      DO i = 2,k
         ref(i).w = ref(i-1).w + xdst
         ref(i).f = QPOL1(ref(i).w,w,y,nref,ntmp)
      ENDDO

      DO i=2,nref
        IF (ref(i).s.NE.ref(i-1).s) GOTO 111
      ENDDO
      DO i=nref,k
        ref(i).s = ref(1).s
      ENDDO
      GOTO 112
 111  CONTINUE
      DO i=1,nref
        y(i) = ref(i).s
      ENDDO
      DO i = 2,k         
         ref(i).s = QPOL1(ref(i).w,w,y,nref,ntmp)
      ENDDO
 112  CONTINUE
      DO i=2,nref
        IF (ref(i).g.NE.ref(i-1).g) GOTO 113
      ENDDO
      DO i=nref,k
        ref(i).g = ref(1).g
      ENDDO
      RETURN
 113  CONTINUE
      DO i=1,nref
        y(i) = ref(i).g
      ENDDO
      DO i = 2,k         
         ref(i).s = QPOL1(ref(i).w,w,y,nref,ntmp)
      ENDDO

      RETURN
      END
