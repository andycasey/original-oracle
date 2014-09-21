C--------------------------------------------------------------------------------------
      SUBROUTINE PARTPRESS
C
C     June,2011  MB     
C     currently works only for atoms, not locked in molecules
C     all PP's and PF's from the model atmosphere are kept fixed
C  
C     Note: the atomic number is here mapped to the equivalent array-element index,
C           thus also in a model atmosphere all atm.h.eps (element abundances) must 
C           be present, but set to zero if UION and PION are not available for these
C           elements
C
C     global variables from lf_param.inc are SUPERSCRIPT
C
C--------------------------------------------------------------------------------------
C
      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'        ! includes lf_param.inc, atm_def.inc

      INTEGER*2   chk
      INTEGER*4   nion_temp(NATOM)
      INTEGER*4   dummy(3,NATOM)
      REAL*4      idx_el_temp(NATOM), eps_temp(NATOM), ref_el(NATOM)
      CHARACTER*2 cel_temp(NATOM)
      REAL*4      uion_temp(NDMAX,NR_ION,NATOM)
      REAL*4      pion_temp(NDMAX,NR_ION,NATOM)
      REAL*4      uion(NDMAX,NR_ION,NATOM)
      REAL*4      pion(NDMAX,NR_ION,NATOM)
      
C reference index for each atom is its atomic number plus ionization stage

      chk = 1                       ! detailed output (for TESTING)
      ref_el = atm.h.idx_el         ! save model atmosphere atoms for further cross-check

      k = 0

      IF (chk) WRITE(*,*) 'Marcs: ',  atm.h.nel
      IF (chk) WRITE(*,*) 'Marcs: ',  atm.h.nion

      DO i = 1, NATOM               ! 92 (lf_param.inc) ==> use as ATOMIC NUMBER identifier
        idx_el_temp(i) = i          ! atomic mass
        nelpos(i)      = i          ! == index number
        cel_temp(i)    = cel(i)     ! Element-Label
        nion_temp(i)   = nion(i)    ! Nr of ionisation stages
        eps_temp(i)    = epss(i)    ! solar abundance from elemdata.inc

        DO m = 1, nion_temp(i)      ! ionization stage
          DO n = 1, atm.h.ndepth
            uion_temp(n,m,i) = 1.   ! not LOG
            pion_temp(n,m,i) = 0.   ! not LOG            
          ENDDO
        ENDDO

        DO j = 1, atm.h.nel              ! 16 in all MARCS models
        IF (k.GT.16) WRITE(*,*) 'Something wrong in partpress.f'

        IF (idx_el_temp(i).EQ.ref_el(j)) THEN
           k = k + 1
           eps_temp(i)  = atm.h.eps(i)   ! N(atm.h.eps)  = 63
           nion_temp(i) = atm.h.nion(j)  ! N(atm.h.nion) = 16
C
C uion and pion are arrays (ndepth, nion, atomic_number)!
C
          DO l = 1, atm.h.nion(j)              ! only one element
           DO n = 1, atm.h.ndepth
            uion_temp(n,l,i) = atm.uion(n,l,atm.h.idx_el(j))  ! assume equal N of depth points in the model atmosphere
            pion_temp(n,l,i) = atm.pion(n,l,atm.h.idx_el(j))  
           ENDDO 
          ENDDO
         ENDIF
        ENDDO

        DO m=1,nion_temp(i)
          IF (chk) dummy(m,i) = 100*i + m
        ENDDO

      ENDDO                              ! finish re-asssignment of arrays

      IF (chk) THEN
      WRITE(*,*) 'Model atmosphere partition functions'
      DO n = 1, atm.h.ndepth
           WRITE(*,"(I4,300F13.6)")
     *                 n,((uion_temp(n,m,k),
     *                 m=1,nion_temp(k)), k=1,NATOM)
      ENDDO
      ENDIF
 
      IF (chk) WRITE(*,99) NATOM-k 
 99   FORMAT('For ',I5,' elements UION and PION will be computed')
      IF (chk) WRITE(*,100) (eps_temp(ii),  ii = 1, NATOM)
 100  FORMAT(/,'New abundances',/,10(F5.2,2X))

      atm.h.nel  = NATOM            ! update total number of elements
      atm.h.nion = nion_temp        ! update ion stages
      atm.h.eps  = eps_temp         ! update abundance array
      atm.h.cel  = cel_temp         ! update element labels
      atm.h.idx_el = idx_el_temp    ! update atomic number array
      atm.h.eps(38) = -0.07         ! update abundance array

C --- load atomic partition functions

      IF (chk) WRITE(*,*)   atm.h.nel
      IF (chk) WRITE(*,*)   atm.h.nion
      IF (chk) WRITE(*,100) atm.h.eps
      IF (chk) WRITE(*,*)   atm.h.cel
      IF (chk) WRITE(*,*)   atm.h.idx_el
      IF (chk) WRITE(*,*)   atm.h.ndepth
      IF (chk) WRITE(*,*)   atm.t
 
      CALL ATOM_PARTFUN(chk,atm.h.ndepth,
     *                  atm.t,atm.h.nel,atm.h.nion,uion)

      DO i = 1, NATOM
        DO m = 1, atm.h.nion(i)      ! ionization stage
          DO n = 1, atm.h.ndepth
           IF (uion_temp(n,m,i).EQ.1.0) THEN
            atm.uion(n,m,i) = uion(n,m,i)       ! new
           ELSE  
            atm.uion(n,m,i) = uion_temp(n,m,i)  ! model atmosphere           
           ENDIF
          ENDDO
        ENDDO
      ENDDO                         ! finish re-assignment of arrays

      IF (chk) THEN 
      WRITE(*,*)'Final PFs: '
      WRITE(*,"(/,4X,300(I8,5X))")
     *         ((dummy(n,m),n=1,atm.h.nion(m)),m=1,atm.h.nel)
      DO n=1, atm.h.ndepth
           WRITE(*,"(I4,300F13.6)")
     *                 n,((atm.uion(n,m,k),
     *                 m=1,atm.h.nion(k)), k=1,atm.h.nel)
      ENDDO 
      ENDIF

C --- compute partial pressures for atoms

      CALL SAHA(chk,atm.uion,pion)

      DO i = 1, NATOM
        DO m = 1, atm.h.nion(i)      ! ionization stage
          DO n = 1, atm.h.ndepth
           IF (pion_temp(n,m,i).EQ.0.0) THEN
            atm.pion(n,m,i) = pion(n,m,i)       ! new
           ELSE  
            atm.pion(n,m,i) = pion_temp(n,m,i)  ! model atmosphere           
           ENDIF
          ENDDO
        ENDDO
      ENDDO                         ! finish re-assignment of arrays

      IF (chk) THEN 
      WRITE(*,*)'Final partial pressures: '
      WRITE(*,"(/,4X,300(I8,5X))")
     *         ((dummy(n,m),n=1,atm.h.nion(m)),m=1,atm.h.nel)
      DO n=1, atm.h.ndepth
           WRITE(*,"(I4,300F13.6)")
     *                 n,((log10(atm.pion(n,m,k)),
     *                 m=1,atm.h.nion(k)), k=1,atm.h.nel)
      ENDDO 
      ENDIF

      RETURN
      END

C---------------------------------------------------------
      SUBROUTINE ATOM_PARTFUN(chk,ndepth,t,nel,nion,uion)
C---------------------------------------------------------

C --- Calculation of Atomic(Ionic) Partition functions UION
C --- by linear interpolation of tables (iondat.inc)
C --- from KURUCZ

      INCLUDE 'lf_param.inc'
      INCLUDE 'iondat.inc'

      INTEGER*2     chk
      PARAMETER     (FAC=1.5664e-02 * 11607.5)
      INTEGER*4     ndepth, nion(nel)
      INTEGER*4     dummy(3,nel)
      REAL*4        t(ndepth), t_chi, tk, du
      REAL*4        uion(80,3,92)

      DO n = 1,ndepth
       DO m = 1,nel
	DO iz = 1,nion(m)
	  t_chi	       = FAC * chi(iz,m)
	  tk	       = t(n) / t_chi - 0.5
	  j	       = MAX(1,MIN(9,INT(tk)))
	  du	       = u_tab(j+1,iz,m) - u_tab(j,iz,m)
	  uion(n,iz,m) = u_tab(j,iz,m) + du * (tk - j)
          IF (chk.AND.n.EQ.1) dummy(iz,m) = 
     *        100*m + iz
	END DO
       END DO
      END DO

      IF (chk) THEN 
           WRITE(*,"('SIU internal partition functions')")
           WRITE(*,"(4X,300(I8,5X))")
     *                ((dummy(n,m),n=1,nion(m)),m=1,nel)
      DO n=1, ndepth
           WRITE(*,"(I4,300F13.6)")
     *                 n,((uion(n,m,k),
     *                 m=1,nion(k)), k=1,nel)
      ENDDO
      ENDIF

      RETURN
      END
C
C---------------------------------------------------------
      SUBROUTINE IONIS(chk,te,divte,tp,pe,n,uion,sion,s)
C---------------------------------------------------------

C --- CALCULATION OF THE IONIZATION EQUILIBRIA
C --- ACCORDING TO THE SAHA - EQUATION  

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
      INCLUDE 'iondat.inc'

      REAL*4   sion(2,NATOM),s(NATOM),zs(NATOM),zzs(NATOM)
      REAL*4   uion(NDMAX,NR_ION,NATOM)
      INTEGER*2 chk, dd

      tpe = tp / pe 
      dd  = n
      nh  = nelpos(1)

C-----IONIZATION EQUILIBRIA FOR EVERY ELEMENT, EXCEPT HYDROGEN

      DO m = 1, atm.h.nel
         id     = atm.h.idx_el(m)    ! atomic-numbers of used elements
         nm     = atm.h.nion(id)     ! Ionisation steps taken into account
         s(m)   = 1.
         fmz    = 1.
         DO i   = 2, nm
            iz  = i - 1
            fmz = fmz * tpe * uion(n,i,id) / uion(n,iz,id) *
     *            EXP( - divte * chi(iz,id) )
            sion(iz,m) = fmz            ! ionization fraction N(i+1)/N(i)
            s(m)       = s(m)   + fmz   ! 1+ N(i+1)/N(i)+ ...
         ENDDO
      ENDDO


C-----IONIZATION EQUILIBRIUM OF HYDROGEN INCLUDING H-MINUS

      sion2nh    = EXP( divte * chi(3,1) ) / ( tpe * uion(n,1,1) )
      sion(2,nh) = sion2nh
      s(nh)      = s(NH)   + sion2nh 


      IF (chk) WRITE(*,"(I4,300(A15,5X))"), dd,
     *          (atm.h.cel(m),m = 1, atm.h.nel)   
      IF (chk) WRITE(*,"(4X,300(F15.7,5X))")
     *          (log10(s(m)), m = 1, atm.h.nel)         ! N(tot)/N(i)
      IF (chk) WRITE(*,"(4X,300(F15.7,5X))")
     *       ((log10(sion(iz,m)), iz =1,atm.h.nion(m)-1),
     *         m = 1, atm.h.nel)

      RETURN
      END 

C---------------------------------------------------------
      SUBROUTINE SAHA(chk,uion,pion)
C---------------------------------------------------------

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'

      REAL*4  te,tp,divte             ! local kinetic temperature
      REAL*4  pe                      ! electron pressure per depth
      REAL*4  pg                      ! gas pressure          
      REAL*4  pk                      ! at pressure     
      REAL*4  uion(80,3,92)           ! partition functions
      REAL*4  pion(80,3,92)           ! partial pressure
      REAL*4  sion(2,NATOM),s(NATOM)  ! temporary arrays
      REAL*4  a, kbo
      INTEGER*2 chk                       

      kbo = 8.61734315D-5              ! k_Boltz in eV/K

      DO k = 1, atm.h.nel
        eps(k) = 10 ** (atm.h.eps(k) - 12.)
      ENDDO

      IF (chk) WRITE(*,*) (atm.h.eps(k), k = 1, atm.h.nel)

      IF (chk) WRITE(*,"('Ionization fractions')")

      DO n   = 1, atm.h.ndepth

         te    = atm.t(n)
         pg    = atm.pg(n)
         pe    = atm.pe(n)
         divte = 1. / (te*kbo)
         tp    = .66667 * te ** 2.5
         pk    = pg - pe

        CALL IONIS(chk,te,divte,tp,pe,n,uion,sion,s)

         DO k = 1, atm.h.nel
           a = eps(k) * pk/s(k)            ! total N for an element from IONIS
           pion(n,1,k) = a
           DO i = 2, atm.h.nion(k)
            ii = i - 1
            pion(n,i,k) = a * sion(ii,k)   ! sion - ionization fraction
           ENDDO
         ENDDO
      IF (chk) WRITE(*,"(300(E13.6,5X))")
     *      ((pion(n,i,k), i=1, atm.h.nion(k)),k=1, atm.h.nel)

      ENDDO
      
      RETURN
      END 
