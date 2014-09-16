C --- -------------------------------- HYDLINE  --------------------------------
      SUBROUTINE HYDLINE(divtn,w,n,kappa_lam)
C --- Calculation of Balmerlines according to the VCS-Theory
C --- DATA FOR LINE-ABSORPTION COEFFICIENTS

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'

      INTEGER*4 itmp                       ! temporary integer variable
      REAL*8    phi_lam                    ! additional broadened profile at lambda
      REAL*8    phi_lam_tmp                ! temporary buffer
      REAL*4    phi_vcs_tmp(NHW)           ! temporary buffer for line profile
      REAL*4    w_vcs_tmp(NHW)             ! temporary wavelengthpoints of vcs-profile
      REAL*4    dlam4                      ! dlam as REAL*4 variable
      PARAMETER (W_HBE = 3.64704 D-05)     ! wavelength of Balmer edge [cm]

      phi_lam   = 0.                       ! Init
      DO k = balmer_min, balmer_max
        dlam   = DABS(w - wbalmer(k))
        IF (dlam .LE. 2.e-6) THEN
          DO i = 1,nhw_nr(k)
            phi_vcs_tmp(i) = phi_vcs(i,n,k)
            w_vcs_tmp(i)   = w_vcs(i,k)
          END DO
          dlam4= dlam
          phi_lam_tmp = 
     *     10**QPOL(dlam4,w_vcs_tmp,phi_vcs_tmp,nhw_nr(k),itmp)
          IF (k.LE.5) THEN
            CALL RESON_RADIA (k,w,dlam,n,divtn,broad)
          ELSE
            broad = 0.
          ENDIF
          phi_lam = phi_lam + (broad + phi_lam_tmp) * C * balmer_gf(k)
        ENDIF
      END DO        
      IF (phi_lam.EQ.0.) RETURN 
      kappa_lam = fh(n) * phi_lam
      kappa_lam = kappa_lam * EXP(-(H_E2_DIV_K * divtn))
      RETURN
      END
