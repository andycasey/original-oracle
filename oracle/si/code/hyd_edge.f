      SUBROUTINE HYD_EDGE(divtn,w,n,kappa_lam)
C-----INCLUSION OF THE QUASI-CONTINUUM NEAR THE SERIES-EDGE W_EDGE(NL)
C-----CORRECTION FOR STIMULATED EMISSION IS NOT INCLUDED
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 W3
      REAL*8 GAUNT(2),W0(2,18),STARK(2,18)
     *       ,ZMN(2,18),ALPHA0(2,18),w_edge(2)
C-----DATA FOR LINE-ABSORPTION COEFFICIENTS

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'

      DATA GAUNT  /1.6678E-2, 2.2896E-3/
      DATA W_EDGE/911.76D-8,3647.04D-8/ 
  
C-----INITIALISATION
      KAPPA_LAM = 0.

C-----DETERMINATION OF THE LOWER PRINCIPAL QUANTUM-NUMBER : N = 1 - LYMAN
C-----                                                      N = 2 - BALMER
      NL = 1
      istart = 0
      IF ( W .GT. 3647.D-8 ) THEN
        nl=2
        istart = NSTART_ESW
      ENDIF
      IF ( ( W.LT.(W_EDGE(NL)+43.D-8)) .AND. (W .GT. W_EDGE(NL))) THEN
        kappa_lam = kappa_lam + GAUNT(NL) * (w**3) * 0.5    ! factor 0.5 introdiced by JKR
     *        * atm.pion(n,1,1)*divtn/(KB*atm.uion(n,1,1))
      ENDIF
C-----FINAL LINE-ABSORPTION COEFFICIENT 
      IF (NL .EQ. 2) 
     *    KAPPA_LAM = KAPPA_LAM * EXP(-(H_E2_DIV_K * DIVTN))
      RETURN
      END
