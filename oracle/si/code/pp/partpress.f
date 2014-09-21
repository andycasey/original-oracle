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
      INTEGER*4   nion_temp(NATOM)
      REAL*4      idx_el_temp(NATOM), eps_temp(NATOM), ref_el(NATOM)
      CHARACTER*2 cel_temp(NATOM)
      REAL*4      uion_temp(NDMAX,NR_ION,NATOM)
      REAL*4      pion_temp(NDMAX,NR_ION,NATOM)

c      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'        ! includes lf_param.inc, atm_def.inc
c      INCLUDE 'iondat.inc'
c      INCLUDE 'ppcommon.inc'

C reference index for each atom is its atomic number plus ionization stage

      ref_el = atm.h.idx_el         ! save model atmosphere atoms for further cross-check

      k = 0
      DO i = 1, NATOM               ! 92 (lf_param.inc) ==> use as ATOMIC NUMBER identifier
        idx_el_temp(i) = i          ! atomic mass
        nelpos(idx_el_temp(i)) = i  ! == index number
        cel_temp(i)    = cel(i)     ! Element-Label
        nion_temp(i)   = nion(i)    ! Nr of ionisation stages
        eps_temp(i)    = epss(i)    ! solar abundance from elemdata.inc

        DO m = 1, nion_temp(i)      ! ionization stage
          DO n = 1, atm.h.ndepth
            uion_temp(n,m,i) = 0.
            pion_temp(n,m,i) = 0.              
          ENDDO
        ENDDO

        DO j = 1, atm.h.nel              ! 16 in all MARCS models
        IF (k.GT.16) WRITE(*,*) 'Something wrong in partpress.f'
        IF (atm.h.idx_el(j).EQ.idx_el_temp(i)) THEN
           k = k + 1
           eps_temp(i)  = atm.h.eps(i)   ! atm.h.neps = 63
           nion_temp(i) = atm.h.nion(j)
          DO l = 1, nion_temp(i)    ! only one element
           uion_temp(n,l,i) = atm.uion(n,atm.h.nion(j),j)
           pion_temp(n,l,i) = atm.pion(n,atm.h.nion(j),j)       
          ENDDO       
         ENDIF
        ENDDO

      ENDDO                              ! finish re-asssignment of arrays

      WRITE(*,*) 
     *   'For ',NATOM-k,' elements UION and PION will be computed'
      WRITE(*,100) (eps_temp(ii),  ii = 1, NATOM)
 100  FORMAT(/,'New abundances',10(F5.2,2X))

      atm.h.nel  = NATOM            ! update total number of elements
      atm.h.nion = nion_temp        ! update ion stages
      atm.h.eps  = eps_temp         ! update abundance array
      atm.h.cel  = cel_temp         ! update element labels
      atm.h.idx_el = idx_el_temp    ! update atomic number array

C --- load atomic partition functions

      CALL ATOM_PARTFUN(atm.h.ndepth,atm.t,atm.h.nel,atm.h.nion)

      DO i = 1, NATOM
        DO m = 1, nion_temp(i)      ! ionization stage
          DO n = 1, atm.h.ndepth
           IF (uion_temp(n,m,i).EQ.0.0.AND.pion_temp(n,m,i).EQ.0.0) THEN
            atm.uion(n,m,i) = uion(n,m,i)       ! new
           ELSE  
            atm.uion(n,m,i) = uion_temp(n,m,i)  ! model atmosphere           
          ENDDO
        ENDDO
      ENDDO                         ! finish re-asssignment of arrays

C --- compute partial pressures for atoms

      CALL SAHA

      RETURN
      END

C---------------------------------------------------------
      SUBROUTINE ATOM_PARTFUN(ndepth,t,nel,nion,uion,chi)
C---------------------------------------------------------

C --- Calculation of Atomic(Ionic) Partition functions UION
C --- by linear interpolation of tables (iondat.inc)
C --- from KURUCZ

      PARAMETER (FAC=1.5664e-02 * 11607.5)
      INTEGER*4 ndepth,nion(nel)
      REAL*4	t(ndepth)

      INCLUDE 'lf_param.inc'
      INCLUDE 'ppcommon.inc'

      DO n = 1,ndepth
       DO m = 1,nel
	DO iz = 1,nion(m)
	  t_chi	       = FAC * chi(iz,m)
	  tk	       = t(n) / t_chi - 0.5
	  j	       = MAX(1,MIN(9,INT(tk)))
	  du	       = u_tab(j+1,iz,m) - u_tab(j,iz,m)
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
C --- ACCORDING TO THE SAHA - EQUATION  

      INCLUDE 'physcnst.inc'
      INCLUDE 'lf_decl0.inc'
      INCLUDE 'iondat.inc'
      INCLUDE 'ppcommon.inc'

      tpe = tp / pe 
      nh  = nelpos(1)

C-----IONIZATION EQUILIBRIA FOR EVERY ELEMENT, EXCEPT HYDROGEN

      DO m = 1, atm.h.nel
         id     = atm.h.idx_el(m)    ! atomic-numbers of used elements
         nm     = atm.h.nion(id)     ! Ionisation steps taken into account
         s(m)   = 1.
         zs(m)  = 0. 
         zzs(m) = 0.
         fmz    = 1.
         IF ( nh .EQ. m ) nm = 2
         DO i   = 2, nm
            iz  = i - 1
            fmz = fmz * tpe * uion(n,i,id) / uion(n,iz,id) *
     *            EXP( - divte * chi(iz,id) )
            dfmz       = FLOAT(iz) * fmz
            sion(iz,m) = fmz            ! ionization fraction N(i+1)/N(i)
            s(m)       = s(m)   + fmz   ! total N
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

      REAL*4  te,tp,divte             ! local kinetic temperature
      REAL*4  pe                      ! electron pressure per depth
      REAL*4  pg                      ! gas pressure          
      REAL*4  pat                     ! at pressure     

      DO n   = 1, atm.h.ndepth

         te    = atm.t(n)
         pg    = atm.pg(n)
         divte = 1. / te
         tp    = .66667 * te ** 2.5
         pe    = atm.pe(n)
         pat   = pg - pe

        CALL IONIS(te,divte,tp,pe,n)

         DO k = 1, atm.h.nel
         m = nelpos(k)               ! index == atomic number
         a = eps(m) * pat/s(m)       ! total N for an element from IONIS
         pion(1,m,n) = a 
         DO i = 2, nion(m)
            ii = i - 1
            pion(i,m,n) = a * sion(ii,m)   !sion - ionization fraction
         ENDDO
         ENDDO

      ENDDO
      
      RETURN
      END 
