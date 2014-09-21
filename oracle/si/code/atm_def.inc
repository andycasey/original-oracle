* ----------------------------------------------------------------------
* >>>
* FILENAME  : ATM_DEF.INC
* DIRECTORY :
* PURPOSE   : Definition of ATM structure
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 22.03.96
* <<<
* ----------------------------------------------------------------------

      STRUCTURE /atm_head2_type/
        CHARACTER*30  date              ! atmosphere generation date
        CHARACTER*20  user              ! atmosphere calc. by user ...
        REAL*4        teff              ! effective temp.
        REAL*4        logg              ! log(g) (gravity)
        REAL*4        x                 ! [H]
        REAL*4        y                 ! [He/H]
        REAL*4        z                 ! [Fe/H] (log)
        REAL*4        xi                ! Microturb.
        REAL*4        alpha             ! l/Hp   Mixinglength
        INTEGER*2     status            ! Atm Status
        
        INTEGER*2     eps_nr            ! Nr of atm element abundances
        REAL*4        eps(NATOM)        ! atm element abundances
        REAL*4        eps_dev(NATOM)    ! temp. deviation of element abundances
        INTEGER*2     ref_nr            ! Nr of Ref.-character (atm)
        CHARACTER*1000 ref              ! atm references   
        INTEGER*2      cmt_nr           ! Nr of atm comment character
        CHARACTER*1000 cmt              ! atm comment

        INTEGER*4   ndepth              ! Nr of depthpoints
        INTEGER*4   nmol                ! Nr of different molecules
        INTEGER*4   nelem               ! Nr of different molecule formingelements
        INTEGER*4   nel                 ! Nr of different elements
        INTEGER*4   nbreak              ! Nr of breakpoints in tau-scale
        INTEGER*4   nscale(MAXBREAK)    ! breakpoints of tau-scale
        INTEGER*4   idx_el(NELMAX)      ! atomic-numbers of used elements
        INTEGER*4   nion  (NELMAX)      ! Ionisation steps taken into account
        INTEGER*4   eta_cont_nr         ! Amount of continuum opacities
        CHARACTER*2 cel (NELMAX)        ! Name of elements concerned
        CHARACTER*4 cmol(NMOLMAX)       ! name of molecules concerned
      END STRUCTURE

      STRUCTURE /atm_type/
        RECORD /atm_head2_type/ h
        REAL*4  tau(NDMAX)              ! tau(reference) depth scale
        REAL*4  logtau(NDMAX)           ! log(tau(reference))
        REAL*4  t  (NDMAX)              ! Temperature structure
        REAL*4  pe (NDMAX)              ! electron pressure per depth
        REAL*4  pg (NDMAX)              ! gas pressure          
        REAL*4  pp (NDMAX,NR_ION)       ! partial pressures of each ion class
        REAL*4  pk (NDMAX)              ! stratif. of Pnuc      
        REAL*4  rho(NDMAX)              ! mass density stratif. 
        REAL*4  pmol(NDMAX,NMOLMAX)     ! molecular partial pressures / partitionfunction
        REAL*4  umol(NDMAX,NMOLMAX)     ! molecular partial pressures / partitionfunction
        REAL*4  pion(NDMAX,NR_ION,NATOM)! atomic partial pressures    / partitionfunction
        REAL*4  uion(NDMAX,NR_ION,NATOM) ! atomic partial pressures    / partitionfunction
        REAL*4  kappa_ref(NDMAX)        ! reference kappa (i.e. kappa(5000A))
        REAL*4  etasig  (NDMAX,MAX_ETA_CONT)! scattering coefficient/kapref
        REAL*4  eta_cont(NDMAX,MAX_ETA_CONT)! continuum opacities
        REAL*8  eta_cont_w(MAX_ETA_CONT) ! Wavelength points of continuum opacities
      END STRUCTURE

      RECORD /atm_type/ atm

      COMMON /XATM/ atm