  
  
      REAL*4 FUNCTION GAHFF(TE,WLAM,ZZ)

C-----CALCULATION OF THE FREE-FREE GAUNT-FACTORS GHFF FOR A GIVEN TEMPERATURE 
C-----TE, WAVELENGTH WLAM AND IONIC CHARGE Z ( ZZ = Z * Z ) IN THE
C-----HYDROGENIC APPROXIMATION BY POLYNOMIAL INTERPOLATION
C-----ACCORDING TO TABLES BY KARZAS AND LATTER

      REAL*4 X,TE,T,GFF1(2,2),GFF2(3,4)
      REAL*8 WLAM        

      DATA GFF1/
     *     2.22E-10,6.7E-07,
     *     5.912E-06,1.0823E0/
  
      DATA GFF2/
     *     0.,2.351431E-18,3.4169006E-13,
     *     5.040843E-21,-7.4482864E-14,-5.7917786E-09,
     *     -2.4385858E-15,1.2822117E-09,2.6061249E-05,
     *     -3.0847031E-12,7.9357168E-07,1.070192E0/
C
      DIVZZ = 1. / ZZ
      X = ZZ * WLAM * 1.E8
      T = TE * DIVZZ
      GAHFF = 0.
      IF ( T .GT. 126010. ) T = 126010. 
      IF ( X .LT.  10000. ) GOTO 1000

C-----GAHFF, FOR X GREATER THAN 10000 A

      IF ( X .GT. 50000. )  X = 50000.
      DO K = 1, 2
         CT = 0.
         DO J = 1, 2
            CT = GFF1(J,K) + CT * T
         ENDDO
         GAHFF = CT + GAHFF * X 
      ENDDO
      RETURN

C-----GAHFF FOR X LESS THAN 10000 A

 1000 DO K = 1, 4
         CT = 0.
         DO J = 1, 3
            CT = GFF2(J,K) + CT * T
         ENDDO
         GAHFF = CT + GAHFF * X 
      ENDDO

      RETURN
      END 