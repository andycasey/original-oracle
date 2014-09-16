      SUBROUTINE HYDLINE(divtn,w,n,kappa_lam)
!
!-----Wavelength interpolation of Balmer line profiles from VCS data of
!-----Michael Lemke (1997, A&AS 122,285). Processed to include convolution with
!-----Doppler broadening, Ali & Griem resonance broadening, and radiation
!-----damping. Accounts for level fine structure in a single symmetric profile.
!-----Profiles for electron densities < 10^10 cm^-3 are dominated by thermal
!-----Doppler movements and are replaced by those for 10^10 cm^-3.
!-----kappa_lam returns the summed contributions of all Balmer lines at w.
!
!-----T. Gehren     25-Jun-05     Adapted from SIU
!
      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
!
      INTEGER*4 idum,istart                             ! Dummy integer variable
      REAL*4 phi,p                              ! Local sum of hydrogen profiles
      REAL*4 exc(2),deltaw
      DATA exc/118353.,140270./            ! Excitation energies for n = 2 and 3
	DATA istart/0/
!
      phi = 0.
      scale = 2.*LOG10(w)+8.
      DO k = minh,maxh
        deltaw = LOG10(REAL(MAX(ABS(w-whyd(k,nhyd)),0.01D-08)))+8.
        p = QPOL(deltaw,w_vcs,phi_vcs(1,n,k),nhw,idum)+scale
        phi = phi+gfhyd(k,nhyd)*10**p
	  IF (n.EQ.50) THEN
	    IF (istart.EQ.0) THEN
	      WRITE(chn,"(10F9.4)") (w_vcs(i),i=1,nhw)
	      WRITE(chn,"(10F9.4)") (phi_vcs(i,n,k),i=1,nhw)
	      WRITE(chn,*) scale
	      istart = 1
	    ENDIF
          WRITE(chn,"(2F11.4,2(1PE14.4))") 1.D8*w,deltaw,p,phi
        ENDIF
      ENDDO
      IF (phi.GT.0.) kappa_lam = fh(n)*phi*EXP(-exc(nhyd)*divtn)
!
      RETURN
      END
