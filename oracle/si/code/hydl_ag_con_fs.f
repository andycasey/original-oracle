      SUBROUTINE HYDL_AG_CON_FS(divtn,w,n,kappa_lam)
!
!-----Wavelength interpolation of Balmer line profiles from VCS data of
!-----Michael Lemke (1997, A&AS 122,285). Processed to include convolution with
!-----Doppler broadening, Ali & Griem resonance broadening, and radiation
!-----damping. Accounts for level fine structure in a single symmetric profile.
!-----Profiles for electron densities < 10^10 cm^-3 are dominated by thermal
!-----Doppler movements and are replaced by those for 10^10 cm^-3.
!-----kappa_lam returns the summed contributions of all Balmer lines at w.
!
!-----T. Gehren      25-Jun-05     Adapted from SIU
!     FRANKS LINUX VERSION fug 12.12.2005 
!
      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
!
      INTEGER*4 idum,istart                             ! Dummy integer variable
      REAL*4 phi,p                              ! Local sum of hydrogen profiles
      REAL*4    phi_vcs_tmp(NNHW)           ! temporary buffer for line profile
      REAL*4 exc(2),deltaw
      DATA exc/118353.,140270./            ! Excitation energies for n = 2 and 3
      DATA istart/0/
!
      phi = 0.
      scale = 2.*LOG10(w)+8.
      nhyd=1
     
      DO k = minh,maxh  !ALL BALMER LINES
        DO i_fs = 1,7   ! ALL FS COMPONENTS
           
           !deltaw = LOG10(REAL(MAX(ABS(w-whyd(k,nhyd)),0.01D-08)))+8.
           deltaw = 
     *          LOG10(REAL(MAX(ABS(w-1.0E-8*wfine(i_fs,k)),0.01D-08)))+8.
           !WRITE(*,*) w,1.0E-8*wfine(i_fs,k)
           DO i = 1,nnhw
              phi_vcs_tmp(i) = nphi_vcs(i,n,k)
                                !IF (n.eq.50) WRITE(*,*) "XXXXX",i,n,k,nphi_vcs(i,n,k)
           END DO
                                !p = QPOL(deltaw,nw_vcs,nphi_vcs(1,n,k),nnhw,idum)+scale
           p = QPOL(deltaw,nw_vcs,phi_vcs_tmp,nnhw,idum)+scale
           phi = phi+gfine(i_fs,k)*10**p
           IF (n.EQ.90) THEN
              IF (istart.EQ.0) THEN
                 WRITE(chn,"(10F9.4)") (nw_vcs(i),i=1,nnhw)
                 WRITE(chn,"(10F9.4)") (nphi_vcs(i,n,k),i=1,nnhw)
                 WRITE(chn,*) scale
                 istart = 1
              ENDIF
              WRITE(chn,"(2F11.4,2(1PE14.4))") 1.D8*w,deltaw,p,phi
              !WRITE(*,*) "P PHI ",p,phi
           ENDIF
        ENDDO
      ENDDO
      IF (phi.GT.0.) kappa_lam = fh(n)*phi*EXP(-exc(nhyd)*divtn)

      RETURN
      END
