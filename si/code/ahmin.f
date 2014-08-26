
C --- -----------------  KAPPA ---------------------------------------  

      REAL*4 FUNCTION AHMIN(wlam)

C-----This function determines the photo-detachment cross-section for 
C-----HMINUS according to T.L John: Astron. Astrophys. 193,189 (1988)

      REAL*8 WLAM
      REAL*4 XMICLAM,XMICLAM3,CROSEC(6)

      DATA CROSEC/152.519,49.534,-118.858,92.536,-34.194,4.982/

      XMICLAM = wlam * 10000.           ! cm -> microns
      DIVLAM  = 1. / XMICLAM
      DIVLAM0 = 1. / 1.6419
      DIFFLAM = ( DIVLAM - DIVLAM0 ) ** .5
      XLAM    = 1.E-06 * XMICLAM
      XMICLAM3= XLAM * XLAM * XLAM

C-----Cross-section for 1250 A <= LAM <= 16419 A (=threshold wavelength)

      AHMIN = 0.
      DO  I = 6,1,-1
         AHMIN = CROSEC(I) + AHMIN * DIFFLAM
      ENDDO
      AHMIN = AHMIN * XMICLAM3 * (DIVLAM - DIVLAM0) ** 1.5

      RETURN
      END 
