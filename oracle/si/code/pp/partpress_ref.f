C--------------------------------------------------------------------------------------
      SUBROUTINE PARTPRESS
C
C     June,2011  MB     
C     currently works only for atoms, not bound in molecules
C--------------------------------------------------------------------------------------

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
      INCLUDE 'iondat.inc'
      INCLUDE 'ppcommon.inc'

C --- Calculation of Atomic Partition functions and their derivatives

      DO i=1,NATOM
         atm.h.nion(i) = nion(i)
      ENDDO
      CALL ATOM_PARTFUN(atm.h.ndepth,atm.t,atm.h.nel,atm.h.nion)
      
      k = 0
      DO i=1,NATOM
        nelpos(i) = 0
        IF (atm.h.nion(i).GT.0) THEN
           k = k + 1
           atm.h.idx_el(k) = i  
         ENDIF
      ENDDO
      atm.h.nel = ke
      DO i=1,ke
        nelpos(atm.h.idx_el(i)) = i
      ENDDO

C --- Change the first nelem elements which are forming molecules

      DO i=1,nelem 
         ii = atm.h.idx_el(i)
         atm.h.idx_el(i)          = atm.h.idx_el(elm_idx(i))
         atm.h.idx_el(elm_idx(i)) = ii
         nelpos(ii)               = nelpos(atm.h.idx_el(i))
         nelpos(atm.h.idx_el(i))  = i
      ENDDO
      atm.h.nelem = nelem

C --- compute Partial pressures for atoms

      CALL SAHA

      RETURN
      END

C---------------------------------------------------------
      SUBROUTINE ATOM_PARTFUN(ndepth,t,nel,nion,uion,chi)
C---------------------------------------------------------

C --- Calculation of Atomic(Ionic) Partition functions UION
C --- by linear interpolation of tables
C --- provided by KURUCZ

      PARAMETER (FAC=1.5664e-02 * 11607.5)
      INTEGER*4 ndepth,nion(nel)
      REAL*4	t(ndepth)

      INCLUDE 'lf_param.inc'
      INCLUDE 'ppcommon.inc'

      DO n = 1,ndepth
       DO m = 1,nel
	DO iz = 1,nion(m)
	  t_chi	     = FAC * chi(iz,m)
	  tk	     = t(n) / t_chi - 0.5
	  j	     = MAX(1,MIN(9,INT(tk)))
	  du	     = u_tab(j+1,iz,m) - u_tab(j,iz,m)
	  uion(n,iz,m) = u_tab(j,iz,m) + du * (tk - j)
	END DO
       END DO
      END DO
      RETURN
      END
C
C---------------------------------------------------------
      SUBROUTINE IONIS(te,divte,tp,pe,n)
C---------------------------------------------------------

C --- CALCULATION OF THE IONIZATION EQUILIBRIA
C --- according to the SAHA - EQUATION  

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
      INCLUDE 'iondat.inc'
      INCLUDE 'ppcommon.inc'

      tpe = tp / pe 
      nh  = nelpos(1)

C-----IONIZATION EQUILIBRIA FOR EVERY ELEMENT, EXCEPT HYDROGEN

      DO m = 1, atm.h.nel
         id     = atm.h.idx_el(m)
         nm     = atm.h.nion(id)
         s(m)   = 1.
         zs(m)  = 0. 
         zzs(m) = 0.
         fmz    = 1.
         IF ( nh .EQ. m ) nm = 2
         DO i   = 2, nm
            iz  = i - 1
            fmz = fmz * tpe * uion(n,i,id) / uion(n,iz,id) *
     *                                   EXP( - divte * chi(iz,id) )
            dfmz       = FLOAT(iz) * fmz
            sion(iz,m) = fmz
            s(m)       = s(m)   + fmz 
            zs(m)      = zs(m)  + dfmz
            zzs(m)     = zzs(m) + FLOAT(iz) * dfmz
         ENDDO
      ENDDO

C-----IONIZATION EQUILIBRIUM OF HYDROGEN INCLUDING H-MINUS

      sion2nh = EXP( divte * chi(3,1) ) / ( tpe * uion(n,1,1) )
      sion(2,nh) = sion2nh
      s(nh)   = s(NH)   + sion2nh 
      zs(mh)  = zs(NH)  - sion2nh
      zzs(mh) = zzs(NH) + sion2nh

      RETURN
      END 

C---------------------------------------------------------
      SUBROUTINE SAHA
C---------------------------------------------------------

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
      INCLUDE 'ppcommon.inc'

      REAL*4  te,tp           ! local kinetic temperature
      REAL*4  pe              ! electron pressure per depth
      REAL*4  pg              ! gas pressure          
      REAL*4  pp (NDMAX,NR_ION)       ! partial pressures of each ion class
      REAL*4  pk              ! stratif. of Pnuc      
      REAL*4  rho             ! mass density stratif. 

      DO n   = 1, atm.h.ndepth
       te    = atm.t(n)
       pg    = atm.pg(n)
       divte = 1. / te
       tp    = .66667 * te ** 2.5
       pe    = atm.pe(n)

       kma1 = atm.h.nelem+1       
       kma2 = atm.h.nelem+2       

       CALL PINIT(te,divte,tp,pg,n,pnuc,pe)

C-----NEWTON-RAPHSON ITERATION OF PARTIAL PRESSURES
C       it = 0
C 500   it = it + 1 
C       f = 0.
C       df = 0.

C --- use only elements from atm.h.idx_el

       DO k = kma1,atm.h.nel
         id = atm.h.idx_el(k)
         x1 = 1. /  s(id)
         x2 = zs(id) * x1
         f  = f  + eps(id) * x2 
         df = df + eps(id) * ( zzs(id) * x1 - x2 * x2 )
       ENDDO

       ge   = f * pnuc - pe
       DGEE = 0.
       DO k = 1, atm.h.nelem
         id   = atm.h.idx_el(k)
         epsi = eps(id)
         pitk = pit(k)
         fac1 = v(k,k) * pitk 
         sum1 = fac1
         DO l = 1, atm.h.nelem
            pitl = pit(l)
            sum1 = sum1 + v(k,l) * pitl 
            sum2 = v(l,l) * pitl
            DO j = 1,atm.h.nelem
               sum2 = sum2 + v(l,j) * pit(j)
            ENDDO
            dg(k,l) = pitk * v(k,l) - epsi * sum2 
         ENDDO
         hk   = sum1 + s(k)
         fac2 = zs(k) * pitk
         g(k) = epsi * pnuc - hk * pitk
         ge   = ge + fac2
         dg(k,k) = dg(K,K) + hk + fac1
         dg(k,kma1) = - fac2 / pe  + epsi
         dg(kma1,k) = - zs(k) - f * sum1
         dgee = dgee + zzs(k) * pitk
       ENDDO

       g(kma1) = ge
       dg(kma1,kma1) = 1. + f + ( dgee + pnuc * df ) / pe

C --- DETERMINATION OF THE PRESSURE CORRECTIONS 
C --- PREPARATIONS FOR THE NEXT ITERATION-STEP

      CALL LINSLV(dg,g,kma1,NMOLMAX)
  
      DO K = 1, kma1
         pit(k) = pit(k) + g(k)
      ENDDO
      pe = pit(kma1)
      pnuc = pg - pe
      DO k = 1,atm.h.nelem
         pitk = pit(k)
       DO l = 1, k
         pnuc = pnuc + pitk * v(k,l) * pit(l)
       ENDDO
      ENDDO

      CALL IONIS(te,divte,tp,pe,n)

C-----CONVERGENCE CHECK

C      DO k = 1, kma1
C         IF ( ABS(g(k) / pit(k)) .GT. .0001 ) GOTO 500
C      ENDDO

C-----FINAL PARTIAL PRESSURES OF THE ELEMENTS

      DO K = 1,atm.h.nel
         M = NELPOS(NEL(K))
         IF ( K .LE. KMOL ) THEN
            A = PIT(K)
         ELSE
            A = EPS(M) * PNUC / S(M)
         ENDIF
         PION(1,M,N) = A
C         PION_U(1,M,N) = PION(1,M,N) / U(1,M,N)
         DO I = 2, NION(M)
            II = I - 1
            PION(I,M,N) = A * SION(II,M)
C            PION_U(I,M,N) = PION(I,M,N) / U(I,M,N)
         ENDDO
      ENDDO

C-----FINAL ELECTRON AND NUCLEAR PRESSURE AS WELL AS MASS-DENSITY AND 
C-----OTHER RELATED PHYSICAL VARIABLES

C      PE(N) = PEL
C      PK(N) = PNUC
C      R = ATW * PNUC
C      RHO(N) = 7.243564E+15 * R * DIVTE 
C      RKT(N) = 1. / R
C
      NH = NELPOS(1)
      DO 800 I = 1, 2
         PP(I,N) = .0
 800  CONTINUE

      DO 810 K = 1, KEL
         IF (  K  .EQ. NH ) THEN
            NIONK = 2
         ELSE
            NIONK = NION(K)
         ENDIF
         DO 820 I = 2, NIONK
            II = I - 1
            FAC = PION(I,K,N) 
            PP(II,N) = PP(II,N) + FAC
 820     CONTINUE
 810  CONTINUE

      ENDDO
      RETURN
      END 
  
C---------------------------------------------------------  
      SUBROUTINE PINIT(TE,DIVTE,TP,PG,N,PNUC,PEO) 
C---------------------------------------------------------
C
C-----THIS SUBROUTINE CALCULATES THE INITIAL VALUES FOR THE ELECTRON, 
C-----NUCLEAR AND PARTIAL PRESSURES DISREGARDING MOLECULE FORMATION!

      PARAMETER (EVK=11607.5)
      COMMON/XEPS/ATW,KEL,NEL(92),NELPOS(92),NION(92),EPS(92),WEIGHT(92)
     *            ,EPS_SAVE(92),NEL_SAVE(92),EPS_ADDI(92)
     *      /XION/CHI(3,92),UION(10,3,92)
     *      /XMOL/NMOL,KMOL,KMA1,KMA2,NC(2,20),REDMA(20),D0(20),
     *            SYM(20),EMAX(20),DMOL(6,10,20)
     *      /WORK/Q(20),V(10,10),SION(2,92),S(92),ZS(92),ZZS(92),
     *            ZZZS(92),F,DF,PIT(11),G(12),DG(12,12)
C
C-----INITIAL VALUE OF THE ELECTRON PRESSURE, ASSUMING FIRST A PURE
C-----HYDROGEN PLASMA INCLUDING HMINUS AND THEN INCLUDING METALS WITH 
C-----LOW IONIZATION-ENERGIES
C
      FH0 = 2. * TP * EXP( - 8754.6 * DIVTE )
      FH1 = .5 * TP * EXP( - 157803.6 * DIVTE )
      X1 = FH0 * FH1
      X2 = FH0 + PG
      PEO = ( SQRT( X1 ) * SQRT ( X1 + X2 * PG  ) - X1 ) / X2

      CALL IONIS(TE,DIVTE,TP,PEO,N)

      NH = NELPOS(1)
      A = EPS(NH) * ZS(NH) / S(NH)
      DO 100 M = 1, KEL
         IF ( CHI(1,M) .LE. 8.*EVK ) THEN
            A = A + EPS(M) * SION(1,M) / ( 1. + SION(1,M) ) 
         ENDIF
 100  CONTINUE
      IF ( A .GT. 0. ) THEN
         PEO = A * PG / ( 1. + A )
         CALL IONIS(TE,DIVTE,TP,PEO,N)
      ENDIF
C
C-----INITIAL GUESSES FOR THE PARTIAL PRESSURES
C
      PNUC = PG - PEO
      DO 200 K = 1, KMOL
         M = NELPOS(NEL(K))
         PIT(K) = EPS(M) * PNUC / S(M)
 200  CONTINUE
      PIT(KMA1) = PEO
C
      RETURN

      END 
  
C -------------------------------------------------------------------   
      SUBROUTINE LINSLV(A,B,NO,NDIM)
C
C-----LINSLV, ADAPTED FROM MIHALAS, HEASLEY AND AUER (1975),
C-----SOLVES A SYSTEM OF LINEAR EQUATIONS ACCORDING TO A
C-----PARTIAL PIVOT ROUTINE
C
      DIMENSION A(NDIM,NDIM),B(NDIM),NUE(60),X(60)
      DOUBLE PRECISION DM,SUM,D(60)
C
      DO 70 I=1,NO
      DO 10 K=1,NO
   10 D(K)=A(K,I)
      I1=I-1
      IF (I1.LT.1) GO TO 40
      DO 30 J=1,I1
      II=NUE(J)
      A(J,I)=D(II)
      D(II)=D(J)
      J1=J+1
      DO 20 L=J1,NO                                                          
   20 D(L)=D(L)-A(L,J)*A(J,I) 
   30 CONTINUE
   40 DM=DABS(D(I)) 
      NUE(I)=I
      DO 50 L=I,NO
      IF (DM.GE.DABS(D(L))) GO TO 50
      NUE(I)=L
      DM=DABS(D(L)) 
   50 CONTINUE
      II=NUE(I)
      A(I,I)=D(II)
      D(II)=D(I)
      I1=I+1 
      IF (I1.GT.NO) GO TO 80
      DO 60 L=I1,NO 
   60 A(L,I)=D(L)/A(I,I)
   70 CONTINUE
C                            
   80 DO 100 L=1,NO 
      II=NUE(L)
      X(L)=B(II)
      B(II)=B(L)
      L1=L+1
      IF (L1.GT.NO) GO TO 110 
      DO 90 J=L1,NO 
   90 B(J)=B(J)-A(J,L)*X(L)
  100 CONTINUE
  110 DO 140 L=1,NO 
      K=NO-L+1
      SUM=0.0D+00
      K1=K+1               
      IF (K1.GT.NO) GO TO 130 
      DO 120 J=K1,NO
  120 SUM=SUM+A(K,J)*X(J)
  130 X(K)=(X(K)-SUM)/A(K,K)
  140 CONTINUE
      DO 150 I=1,NO 
  150 B(I)=X(I)
C
      RETURN
C
      END 

