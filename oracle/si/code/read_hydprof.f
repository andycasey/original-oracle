      SUBROUTINE READ_HYDPROF
!
!-----Input of tabulated VCS Balmer profiles (see also: HYDLINE)
!-----Electron densities are represented by a minimum value of 10^10 cm^-3
!-----Interpolation of temperature, electron (and hydrogen) densities
!-----for each depth point and all hydrogen lines in wavelength interval
!
C     THIS IS FUGs LINUX VERSION!


      INCLUDE 'physcnst.inc'                                ! Physical constants
      INCLUDE 'lf_decl0.inc'                               ! Common declarations
      INTEGER*4 idum,nnb,nwords                        ! Dummy integer variables
      REAL*4 nelog                                        ! log electron density
      REAL*4 tlog                                 ! log temperature (at depth n)
      REAL*4 nhlog                                        ! log hydrogen density
      REAL*4 grid_vcs(nnhw,nnhh,nnht,nnhe,nnhmax)                 ! VCS grid profiles
      REAL*4 vcs_ijkm(nnhh)                                    ! Temporary buffer
      REAL*4 vcs_ijm(nnht)                                     ! Temporary buffer
      REAL*4 vcs_im(nnhe)                                      ! Temporary buffer
      REAL*4 dw0,ddw,ne0,dne,nh0,dnh,lgt0,dlgt,dwalph
!
      maxh = 0
      minh = 1
      nhyd=1
      DO m = 1,nhmax-nhyd+1
        IF (wmin.LE.whyd(m,nhyd)+1.D-6)  maxh = m
      ENDDO
      DO m = nhmax-nhyd+1,1,-1
        IF (wmax.GE.whyd(m,nhyd)-1.D-6)  minh = m
      ENDDO

      env_len = IGETENV('SIU_MAIN',env)
      OPEN(UNIT=28,FILE=env(1:env_len)//'/data/vcs__balmer.dat',
     *  STATUS='UNKNOWN',
     *  FORM='FORMATTED',READONLY)
      READ(28,*) dw0,ddw,ne0,dne,nh0,dnh,lgt0,dlgt,
     *  (nlow(m),nup(m),nhhnr(m),((((grid_vcs(i,j,k,l,m),
     *  i=1,nnhw),j=1,nhhnr(m)),k=1,nnht),l=1,nnhe),m=1,maxh) 
      CLOSE(28)
     
      
!
! --- Calculate relevant Balmer line profiless using quadratic interpolation
!
      DO k = 1,nnht
        nt_vcs(k) = lgt0+(k-1)*dlgt
      ENDDO
      DO l = 1,nnhe
        ne_vcs(l) = ne0+(l-1)*dne
      ENDDO
      DO j = 1,nnhh
        nh_vcs(j) = nh0+(j-1)*dnh
      ENDDO
      DO i = 1,nnhw
        nw_vcs(i) = dw0+(i-1)*ddw
      ENDDO
!
      DO m = minh,maxh                              ! Begin of Balmer lines loop
        DO n = 1,atm.h.ndepth                     ! Begin of stratification loop
          tlog = LOG10(atm.t(n))
          nelog = LOG10(atm.pe(n))-tlog-LOG_K
          IF (nelog.LT.10.) nelog = 10.           ! Replace by lower table limit
          nhlog = LOG10(atm.pion(n,1,1))-tlog-LOG_K
          IF (m.EQ.minh) THEN           ! Check temperature and electron density
            IF ((tlog.LT.nt_vcs(1)) .OR. (tlog.GT.nt_vcs(nht))) THEN
              WRITE(chn,*) 'Temperature interpolation out of bounds'
              lf_status = JIBSET(lf_status,0)                      ! Fatal Error
              lf_status = JIBSET(lf_status,4)                      ! Fatal Error
              RETURN
            ENDIF
            IF (nelog.GT.ne_vcs(nnhe)) THEN
              WRITE(chn,*) 'Electron density interpolation out',
     &          ' of bounds => extrapolation'
              lf_status = JIBSET(lf_status,0)                      ! Fatal Error
              lf_status = JIBSET(lf_status,5)                      ! Fatal Error
              RETURN
            ENDIF
          ENDIF
!
!-----Second order interpolation of PHI(lambda) as a function of (NH,T,NE)
!
          DO i = 1,nnhw
            DO j = 1,nnhe
              DO k = 1,nnht
                DO l = 1,nhhnr(m)
                  vcs_ijkm(l) = grid_vcs(i,l,k,j,m)
                ENDDO
                IF (nhhnr(m).GT.1) THEN
                  vcs_ijm(k) = QPOL(nhlog,nh_vcs,vcs_ijkm,nnhh,idum)
                ELSE
                  vcs_ijm(k) = vcs_ijkm(1)             ! No resonance broadening
                ENDIF
              ENDDO
              vcs_im(j) = QPOL(tlog,nt_vcs,vcs_ijm,nnht,idum)
            ENDDO
            nphi_vcs(i,n,m) = QPOL(nelog,ne_vcs,vcs_im,nnhe,idum)
            !WRITE(*,*) "YYY",nphi_vcs(i,n,m),i,n,m
          ENDDO
        ENDDO                                              ! Stratification loop
      ENDDO                                                ! Hydrogen lines loop
!
      RETURN
      END
