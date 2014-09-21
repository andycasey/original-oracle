

C --- Quadratic Interpolation ------------------------------------------------
      REAL*4 FUNCTION QPOL1(x,a,f,no,k)
      INTEGER*4  i,k,no,n
      REAL*4     f(no)
      REAL*8     u,v,w,x,y,z,a(no)

      k = no - 1
      IF (x .GT. a(no)) GOTO 50
      k = 1
      IF (x .LT .a(1) .OR. no .EQ. 2) GOTO 50
      n = no
   10 IF (k .EQ. n-1) GOTO 30
      i = (k + n) / 2
      IF (a(i).GE. x) GOTO 20
      k = i
      GOTO 10
   20 n = i
      GOTO 10
   30 IF (k .EQ. no-1) k = k - 1
      u = x - a(k)
      v = x - a(k+1)
      w = x - a(k+2)
      y = (f(k)  - f(k+1)) / (v-u)
      z = (f(k+1)- f(k+2)) / (w-v)
      IF (y*z .GT. .0) GOTO 40
      QPOL1 = f(k) + u * (y+v*(y-z)/(w-u))
      RETURN
   40 QPOL1 = f(k) + u*y / (1.-v*(y-z)/(f(k)-f(k+2)))
      RETURN
   50 QPOL1 = f(k) + (f(k+1)-f(k)) * (x-a(k))/(a(k+1)-a(k))
      RETURN                      
      END
