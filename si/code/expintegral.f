
*                           
* --------------------------------------------------------------------------
* >>>
* NAME        : EXPINTEGRAL
* FILENAME    : EXPINTEGRAL.FOR
* DIRECTORY   : [REETZ.PROGLIB]
* 
* PURPOSE     : Exponential-Integralfunction
* MODUL       : FORTRAN - Subroutine
* CALLING SEQ.: EXPINTEGRAL(x,n)
* INPUTS      : x Argument
*               n order (1,2)
* OUTPUTS     : Return-Value: 
* AUTHOR      : Johannes Reetz
* DATE        :  9-OCT-1992 12:23:37.56
* <<<
* --------------------------------------------------------------------------
* 

      REAL*4 FUNCTION expintegral(x,n)
      REAL*4    x
      REAL*8    y
      INTEGER*4 n
      REAL*8 A1,A2,A3,A4,B1,B2,B3,B4,C0,C1,C2,C3,C4,C5

C --- Coefficients for rational approx. 5.1.56
      PARAMETER (A1 =  8.5733287401)
      PARAMETER (A2 = 18.0590169730)
      PARAMETER (A3 =  8.6347608925)
      PARAMETER (A4 =  0.2677737343)
      PARAMETER (B1 =  9.5733223454)
      PARAMETER (B2 = 25.6329561486)
      PARAMETER (B3 = 21.0996530827)
      PARAMETER (B4 =  3.9584969228)

C --- Coefficients for rational approx. 5.1.53
      PARAMETER (C0 = -0.57721566)
      PARAMETER (C1 =  0.99999193)
      PARAMETER (C2 = -0.24991055)
      PARAMETER (C3 =  0.05519968)
      PARAMETER (C4 = -0.00976004)
      PARAMETER (C5 =  0.00107857)

C --- Polynomial approximation 5.1.53  (eps < 2.e-7)

      GOTO (100,200), n
      WRITE(*,*) 'ERROR in EXPINTEGRAL: Ordernumber of E2-function > 2'
      expintegral = 0.      
      RETURN

100   CONTINUE
C --- Rational approximation 5.1.56    (eps < 2.e-8)
      IF (x.LT.1.) THEN 
        IF (x.GT.0.) expintegral=C0+x*(C1+x*(C2+x*(C3+x*
     *                        (C4+x*C5))))-log(x)
      ELSE
        expintegral = exp(-x)*(A4+x*(A3+x*(A2+x*(A1+x)))) / 
     *                   (x*(B4+x*(B3+x*(B2+x*(B1+x)))))
      ENDIF
      RETURN

200   CONTINUE
C --- Recurrence relation    5.1.14
      IF (x.LT.1) THEN
        expintegral = 1.
        IF (x.GT.0.) THEN 
            expintegral = exp(-x)-x*(C0+x*
     *                      (C1+x*(C2+x*(C3+x*(C4+x*C5))))-log(x))
        ENDIF
      ELSE                     
C --- Rational approximation 5.1.56    (eps < 2.e-8)
          expintegral = exp(-x)*(1 - (A4+x*(A3+x*(A2+x*(A1+x)))) / 
     *                          (B4+x*(B3+x*(B2+x*(B1+x))))  )
      ENDIF
      RETURN
      END 
