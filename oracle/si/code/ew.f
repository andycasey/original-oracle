
!     Calculation of equivalentwidth
      REAL*4 FUNCTION EW(mode,x,y,n,xb,xe)
!     mode=0 --> equivalentwidth   (ew > 0 : absorption  ew<0 : emission)
!     mode=1 --> integral
!     x,y    --> function
!     xb,xe  --> start and end value
      INTEGER*4   mode,n
      REAL*8 x(n),xb,xe, sum
      REAL*8 y(n)

      EW  = 0.D0
***      print*, x,y,n,xb,xe
      ib = 1
      ie = n
      DO i=2,n
       IF (x(i-1).LT.xb) THEN
          ib = i
       ELSE
          GOTO 10
       ENDIF
      ENDDO
 10   CONTINUE
      DO i=n-1,ib,-1
       IF (x(i+1).GT.xe) THEN
          ie = i
       ELSE
          GOTO 20
       ENDIF
      ENDDO
 20   CONTINUE
      
      sum = 0.D0
      IF (mode.EQ.0) THEN
         DO i=ib,ie-1
           sum = sum + DBLE(2. - ( y(i) + y(i+1) ) )  * 
     *                       0.5D0 * (x(i+1) - x(i))
         END DO
         EW = sum
         RETURN
      ELSEIF (mode.EQ.1) THEN
         DO i=ib,ie-1
           sum = sum + DBLE( y(i) + y(i+1) ) * 
     *                       0.5D0 * ( x(i+1) - x(i) )
         END DO
         EW = sum
         RETURN
      ELSE
         WRITE(*,*) ' EW: Unknown mode ! '
      ENDIF

      print*, 'ew ', ew

      RETURN
      END
