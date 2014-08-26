  
  
  
C-----CALCULATION OF THE BOUND-FREE ABSORPTION COEFFICIENT
      REAL*4 FUNCTION AHYDBF(DIVTE,WLAM,Z,EMI)

      REAL*8 WLAM,WLAMZ
      REAL*4 GAUNT(7),GHBF1(5,3),GHBF2(5),DIVTE,Z,EMI
      DATA GAUNT/7*1./
      DATA GHBF1/
     *     1.105,1.101,1.101,1.102,1.099,
     *     -7.922E-05,-3.290E-05,-1.923E-05,-1.304E-05,-9.020E-06,
     *     4.536E-09,1.152E-09,5.110E-10,2.638E-10,1.367E-10/
  
      DATA GHBF2/2.02848E5,1.31854E4,2.06536E-1,9.36244E4,2.21980E3/

C-----DETERMINATION OF THE LOWER AND UPPER SUMMATION INDICES N1, N2 FOR THE
C-----THE FUNCTION F(N1,N2,1,WLAM,TE) AND THE GAUNT-FACTORS FOR LEVELS N = 1 - 15  
      ZZ = Z * Z
      WLAMZ = ZZ * (WLAM * 1.E8)             ! cm -> A
      N1 = 1. + SQRT( WLAMZ / 911.76 )
      N2 = MAX(N1,10)
      IF ( N1 .LE. 6 ) THEN
         DO 100 J = N1, 6
            IF ( J .EQ. 1 ) THEN
               GAUNT(J) = ( GHBF2(1) + WLAMZ * ( GHBF2(2) - WLAMZ ) ) * 
     *                      GHBF2(3) / ( GHBF2(4) 
     *                                  + WLAMZ * ( GHBF2(5) + WLAMZ ) )
            ELSE
               I = J - 1
               GAUNT(J) = GHBF1(I,1) + WLAMZ * ( GHBF1(I,2)
     *                                            + WLAMZ * GHBF1(I,3) )
            ENDIF
 100     CONTINUE
      ENDIF
C
C-----DETERMINATION OF THE FUNCTION F(N1,N2,1,WLAM,TE) 
C
      X1 = 1.57804E+05 * ZZ * DIVTE
      DIVX1 = 1. / X1
      IF ( N1 .LT. 10 ) THEN
         N3 = N2 + 1
         G = 1. / ( N3 * N3 ) 
         B0 = X1 * G
         EB0 = EXP(B0)
         F = .5 * ( EB0 - 1. ) * DIVX1
         DO 200  K = N1, N2
            KINDEX = MIN(K,7) 
            K2 = K * K
            K3 = K2 * K
            FK = GAUNT(KINDEX) * EXP( X1 / K2 ) / K3
            F = F + FK
 200     CONTINUE
      ELSE
         DIV = 1. / EMI
         F = .5 * DIVX1 * ( DIV - 1. )
      ENDIF
      FAC = 2.08952E-18 * WLAM * WLAMZ * WLAMZ * EXP(-X1)
      AHYDBF = FAC * F
      RETURN
      END 
