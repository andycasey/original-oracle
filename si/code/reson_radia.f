

C --- ------------------------------ RESON_RADIA -------------------------------
      SUBROUTINE RESON_RADIA (line_nr,w,dw,n,divtn,broad,rescale)
C --- Calculates
C ---  i) the resonance broadening effect on Balmer lines ( H-alpha ...
C ---     ... H-epsilon ) according to Ali and Griem, Phys.Rev.
C ---     140A,1044 (1965) and 144,366 (1966),
C ---     see also ---> Cayrel and Traving, ZfA 50,239 (1960) and
C ---     Furssow and Wlassow, Phys. Z. Sowjetunion 10,378 (1936) .
C ---     Resonance damping constants GAMMA_RES of upper and lower levels
C ---     are simply added.
C --- ii) broadening due to radiation damping ( H-alpha .... H-delta )
C     Routine adapted from Klaus Fuhrmann
 
      INCLUDE 'physcnst.inc'                     ! physical constants
      INCLUDE 'lf_decl0.inc'                     ! general declarations

      INTEGER*4 line_nr                          ! Nr of Balmerline (H-alpha=1)
                                  ! line_nr = 1,2,3,4,5 <=> H-alpha...H-epsilon
      REAL*4    gamma_res
      REAL*4    gamma_rad
      REAL*4    alpha
      REAL*4    alpha2
      REAL*4    w2
      REAL*4    aii(6,7)                     
      REAL*4    l1nf1n(6)                        ! lam(1,n) * f(1,n)
      REAL*4    c3
      REAL*8    rescale


      DATA l1nf1n/.5059619e-5,.811345e-6,.281938e-6,
     *            .132394 e-6,.73139 e-7,.44806 e-7/
      DATA aii/   0.,        0.,        0.,        0.,        0., 0.,
     *     4.699E+08,        0.,        0.,        0.,        0., 0.,
     *     5.575E+07, 4.410E+07,        0.,        0.,        0., 0.,
     *     1.278E+07, 8.419E+06, 8.986E+06,        0.,        0., 0.,
     *     4.125E+06, 2.530E+06, 2.201E+06, 2.699E+06,        0., 0.,
     *     1.644E+06, 9.732E+05, 7.783E+05, 7.711E+05, 1.010E+06, 0.,
     *            0.,        0.,        0.,        0.,        0., 0./ ! dummies
C            LYMAN     BALMER     PASCHEN    BRACKETT    PFUND

      gamma_res = 0.                      
      gamma_rad = 0.                      

      c3 = C3FAC * (l1nf1n(1) + l1nf1n(line_nr+1))
      gamma_res = PI2_4_DIV_K * c3 * divtn * 
     *                   atm.pion(n,1,1)/atm.uion(n,1,1)
      gamma_rad = aii(1,2)                

      DO i = 1,line_nr+1
        gamma_rad = gamma_rad + aii(i,line_nr+2) 
      END DO

      w2     = w * w
      alpha  = (gamma_res*rescale + gamma_rad) * w2 * DIV_4PIC
      alpha2 = alpha * alpha              
      broad  = alpha / ( (dw*dw) + alpha2 ) * DIV_PIC * w2

      RETURN
      END
