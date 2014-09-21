c --- According to KURUCZ ATLAS9 !!!
c --- adapted by JKR
c --- correct of true absorption sourcefunction Slam for scattering term (beta*Jlam)
      SUBROUTINE ITERATE_SOURCE_FCT(ndepth,tau_lam,beta,slam,jlam,
     *                              tau_max)
      INTEGER*4 ndepth,n,i,k,itmp,k1,k2
      INTEGER*4 ITMAX
      INTEGER*4 n_max                ! for backward interpolation of tau-scale
      REAL*4    tau_lam(ndepth)                            ! monochromatic depthscale
      REAL*4    beta(ndepth)                               ! sigma / (kappa + sigma)
      REAL*4    slam(ndepth),jlam(ndepth)                  ! sourcefunction, mean intensity
      REAL*8    sum
      REAL*4    xold,dsmax,eps,tau_max
      PARAMETER (EPS   = 1.E-5)                            ! max. correction of Slam
      PARAMETER (EPS2  = 1.E-4)                            ! max. correction of Slam
      
      INCLUDE   'lf_decl3.inc'                             ! nxdepth,taux_lam,coefj(k,l)
      PARAMETER (ITMAX = 100)                              ! maximum nr of iterations

      REAL*8    strue(NDMAX1),xtrue(NDMAX1)                ! (1.-betax(n))*strue(n)
      REAL*4    sxlam(NXDEPTH),jxlam(NXDEPTH)              ! sourcefunction, mean intensity
      REAL*4    dsxlam(NXDEPTH)                            ! correction of sourcefunction
      REAL*4    betax(NXDEPTH)                             ! sigma / (kappa + sigma)
      REAL*8    div_diag(NXDEPTH)

      n_max = 0
      DO n=1,ndepth
         IF (tau_lam(n).LE.taux_lam(nxdepth)) n_max = n
      ENDDO

C --- 0) interpolate slam(tau_lam) ---> sxlam (taux_lam) --------------

C      k1 = MAP1(tau_lam,slam,ndepth,taux_lam,sxlam,nxdepth)
      DO n=1,nxdepth
         sxlam(n) = QPOL(taux_lam(n),tau_lam,slam,ndepth,k1)
         betax(n) = QPOL(taux_lam(n),tau_lam,beta,ndepth,k2)
         strue(n) = sxlam(n)
         div_diag(n)  = 1. / (1. - betax(n) * coefj(n,n))
         xtrue(n) = (1. - betax(n)) * strue(n) 
      ENDDO

C --- 1) starting with Slam_0 = slam (true absorption sourcefunction) -
     
      DO itmp=1,itmax
C --- 2) calculate dsxlam ---------------------------------------------
       dsmax = 0.
       xold  = -1.
       i     = 0
       DO n=1,nxdepth
         sum = 0.D0
         DO k=1,nxdepth
            sum = sum + coefj(n,k)*sxlam(k)
         ENDDO       
         jxlam(n)  = sum
         dsxlam(n) = div_diag(n)*(xtrue(n) - sxlam(n) + betax(n)*sum)
         dsmax = AMAX1(dsmax,ABS(dsxlam(n)/sxlam(n)))
         IF (xold.NE.dsmax) THEN
            xold = dsmax
            i    = n
         ENDIF
       ENDDO
       DO n=1,ndepth
         sxlam(n) = sxlam(n) + dsxlam(n)
       ENDDO

C --- 3) check if correction small enough -----------------------------
       IF (dsmax.LE.EPS) THEN
         k = itmp
         GOTO 10
       ENDIF  
      ENDDO

C --- if concergence not achieved -------------------------------------
      k = ITMAX
      WRITE(*,9999) i,sxlam(i),dsmax
 9999 FORMAT(' Convergence of S not achieved: n=',I3,
     *       '  (S,dS): ',2G13.6)
      DO n=1,nxdepth
        sxlam(n) = (1.-betax(n))*strue(n) + betax(n)*jxlam(n)
      ENDDO
C --- -----------------------------------------------------------------
   10 CONTINUE 
      DO n=MAX(1,n_max-2),ndepth
         strue(n) = slam(n) 
         xtrue(n) = (1. - beta(n))*strue(n)
      ENDDO
      
      DO n=1,n_max
         slam(n) = QPOL(tau_lam(n),taux_lam,sxlam,nxdepth,k1)
         jlam(n) = QPOL(tau_lam(n),taux_lam,jxlam,nxdepth,k2)
      ENDDO

C --- --- handle all depthpoints with tau_lam > MAX(taux_lam) ---------
      IF (tau_max.LT.taux_lam(NXDEPTH)) THEN
         DO n = n_max+1,ndepth
            slam(n) = strue(n) 
         ENDDO
         GOTO 90     ! finished      
      ENDIF
      i=0
   20 CONTINUE 
      CALL DERIV(ndepth,MAX(n_max-2,1),tau_lam,slam,dsxlam,2) 
      DO n=MAX(n_max,1),ndepth
         xold    = slam(n)         
         slam(n) = xtrue(n) + beta(n)*(dsxlam(n) + slam(n))
         dsmax = AMAX1(dsmax,ABS((xold-slam(n))/slam(n)))
      ENDDO
      i = i + 1
      IF (i.GT.ITMAX) THEN
        WRITE(*,9999) n_max+1,dsmax 
        GOTO 90
      ENDIF
      IF (dsmax.GT.EPS2) GOTO 20  
C --- --- -------------------------------------------------------------
 90   CONTINUE

      RETURN
      END
