      SUBROUTINE SORT(nmax,ra,idx)
C --- see Heapsort: Numerical Recipes [Fortran] (1989,P.231)       

      INTEGER*4  nmax,idx(nmax),ibb
      REAL*8 raa,ra(nmax)

      n = nmax/2 + 1
      k = nmax
   10 CONTINUE
        IF (n.GT.1) THEN
           n = n - 1
           raa = ra(n)
           ibb = idx(n)
        ELSE
           raa = ra(k)
           ibb = idx(k)
           ra(k)  = ra(1)
           idx(k) = idx(1) 
           k        = k - 1
           IF (k.EQ.1) THEN
              ra(1) = raa
              idx(1)= ibb
              RETURN    
           ENDIF
        ENDIF
        i = n
        j = n + n
   20   IF (j.LE.k) THEN
          IF (j.LT.k) THEN
            IF ( ra(j) .LT. ra(j+1) ) j = j + 1            
          ENDIF
          IF ( raa .LT. ra(j) ) THEN
             ra(i) = ra(j)
             idx(i)= idx(j)
             i = j
             j = j + j
          ELSE
             j = k + 1
          ENDIF
        GOTO 20
        ENDIF
        ra(i) = raa
        idx(i)= ibb
      GOTO 10
      END
