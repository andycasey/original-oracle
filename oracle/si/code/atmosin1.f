
C ---------------------------------  ATMOSIN1    ------------------------------
C --- ATMOSPHERE-INPUT(pion,uion,pmol,umol,.....) -----------------------------
C --- This program calculates the atmospheric stratification by linear 
C --- interpolation within a 3-dimensional precalculated (TEFF,LOGG,[M/H])-grid
C --- Interpolation Algorithm by K.Fuhrmann
C 
C                         
C            [T1,G1,M1]----------------------[T2,G1,M1]
C                /                               /
C               /                               /
C              /                               /
C             /                [T1,G1,M2]...../----------------[T2,G1,M2]
C            /                     .         /                     /
C           /                     .         /                     /
C          /                     .         /                     /
C         /                     .         /                     /
C    [T1,G2,M1]----------------------[T2,G2,M1]                /
C                             /                               /
C                            /                               /
C                           /                               /
C                      [T1,G2,M2]----------------------[T2,G2,M2]
C

      SUBROUTINE ATMOSIN1
           

      INCLUDE 'physcnst.inc'            ! physical constants
      INCLUDE 'atmgrid.inc'             ! atmosphere filenames
      INCLUDE 'lf_decl0.inc'            ! general declarations

      STRUCTURE /atmgrid_type/          ! map segments of atmosphere grid area
        INTEGER*4 idx_teff              ! all logarithmic (exept 't')
        INTEGER*4 idx_logg
        INTEGER*4 idx_z
        
        INTEGER*4 ndepth                ! Nr of depthpoints
        INTEGER*4 nmol                  ! Nr of different molecules
        INTEGER*4 nel                   ! Nr of different elements
        INTEGER*4 nbreak                ! Nr of breakpoints in tau-scale
        INTEGER*4 eta_cont_nr           ! Amount of continuum opacities
        INTEGER*4 nscale(MAXBREAK)      ! breakpoints of tau-scale
        INTEGER*4 idx_el(NELMAX)        ! atomic-numbers of used elements
        INTEGER*4 nion  (NELMAX)        ! Ionisation steps taken into account
        
        RECORD /atm_head_type/ h
      
        REAL*4  tau(NDMAX)               ! tau(reference) depth scale
        REAL*4  t  (NDMAX)               ! Temperature structure
        REAL*4  pe (NDMAX)               ! electron pressure per depth
        REAL*4  pg (NDMAX)               ! gas pressure          
        REAL*4  pk (NDMAX)               ! stratif. of Pnuc      
        REAL*4  rho(NDMAX)               ! mass density stratif. 
        REAL*4  pmol(NDMAX,NMOLMAX)      ! molecular partial pressures/partitionfunction
        REAL*4  umol(NDMAX,NMOLMAX)     ! molecular partial pressures/partitionfunction
        REAL*4  pion(NDMAX,NR_ION,NELMAX)! atomic partial pressures/partitionfunction
        REAL*4  uion(NDMAX,NR_ION,NELMAX)! atomic partial pressures/partitionfunction
        REAL*4  kappa_ref(NDMAX)         ! reference kappa (i.e. kappa(5000A))
        REAL*4  eta_cont(NDMAX,MAX_ETA_CONT)! background opacities
        REAL*4  etasig(NDMAX,MAX_ETA_CONT)  ! scattering coeff/kappa-ref
        REAL*4  eta_cont_w(MAX_ETA_CONT) ! wavelength continuum opacities
      END STRUCTURE                  

      INTEGER*2 ua_ptr
      INTEGER*2 uag_ptr

      RECORD  /atm_type/ ua(N_ATM_HOLD)     ! already used atm are stored  
      RECORD  /atmgrid_type/ uag(N_ATMGRID_HOLD) ! used atmgridpoints
      RECORD  /atmgrid_type/ atm_int(8)     ! gridpoints for interpolation 

      REAL*8    dtau                           ! delta tau
      REAL*4    x,y,z                          ! temp. variables for calc. of interp.coefficients
      REAL*4    co(8)                          ! interpolation coefficients
                            
      LOGICAL   atmgrid_exist                  ! result of request whether gridatmosphere files exist
      LOGICAL   atm_already_hold               ! result of request whether atmosphere files exist
      LOGICAL   exact_atm_exist,atm_exist      ! result of request whether exact atmosphere file exist
      LOGICAL   exact_grid                     ! result of request whether exact grid atmosphere exist

      INTEGER*4 it,ig,iz                       ! selected grid points
      INTEGER*4 it_mod,ig_mod,iz_mod           ! it+MOD(k,2),ig+MOD(k/2,2),.... 
      INTEGER*4 iteff,ilogg,iz_tmp             ! Integerbuffer for Teff, logg and z
      INTEGER*4 n_unit                         ! Variable for Unitnumbers

      INTEGER*4 idx_atmgrid(8)                 ! indices of uag if found
      
      REAL*4    sum                            ! sum of something
      INTEGER*4 nbreak,ndepth,nel,nmol,eta_cont_nr  ! temporary variables
      INTEGER*4 ic                             ! counter for coefficents co
      

      CHARACTER*255 atm_file_tmp
      CHARACTER*55 filename_tmp
      CHARACTER*25 atm_fn(8)

      COMMON/ATMOSIN_WORK/x,y,z,atmgrid_exist,atm_already_hold,
     *             exact_atm_exist,it,ig,iz,
     *             iteff,ilogg,iz_tmp,n_unit,dtau,atm_file_tmp,
     *             filename_tmp,atm_fn,atm_int,idx_atmgrid,
     *             ua,uag,ua_ptr,uag_ptr

      COMMON/ATMINT/co

      DATA ua_ptr  /0/
      DATA uag_ptr /0/

C --- Determination of grid filenames
      exact_atm_exist = .FALSE.
      exact_grid = .FALSE.

C --- -------------------------------------------------------
C --- CHECK if adapted atmosphere is already hold in UA
C --- -------------------------------------------------------

      atm_already_hold = .FALSE.
      IF (BJTEST(lf_ctrl,3)) THEN       ! IF search for exact atmosphere
       DO i=1,N_ATM_HOLD
        IF (ua(i).h.teff.EQ.0.) GOTO 2
        IF (atm.h.teff.EQ.ua(i).h.teff.AND.
     *      atm.h.logg.EQ.ua(i).h.logg.AND.
     *      atm.h.z   .EQ.ua(i).h.z        ) THEN
          atm = ua(i)                    ! assign existing atmosphere 
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *          'Atmosphere assigned from stack: ' //
     *           atm.h.ref(1:INDEX(atm.h.ref,'//'))
          GOTO 8000   ! jump to end of procedure
        ENDIF
       ENDDO
      ENDIF
 2    CONTINUE
      
C --- ---------------------------------------

      exact_atm_exist = .FALSE.

      IF (BJTEST(lf_ctrl,3).OR.BJTEST(lf_ctrl,22)) THEN  ! search first for exact atmosphere
       IF (BJTEST(lf_ctrl,22)) THEN
        exact_atm_exist = .TRUE.         
        GOTO 200                             ! jump to the begin of ATM-loop
       ELSE
         iteff  = INT((atm.h.teff + 5) / 10 )
         ilogg  = INT((atm.h.logg + 0.005) * 100)
         IF (atm.h.z .GE. 0) THEN 
           iz_tmp = NINT(100 * atm.h.z)              
         ELSE  
           iz_tmp = NINT(100 * atm.h.z)              
         ENDIF

         filename_tmp        = ' '
         filename_tmp(1:1)   = 't'
         ENCODE(3,101,filename_tmp(2: 4)) iteff         ! TEFF NAME
         ENCODE(3,101,filename_tmp(5: 7)) ilogg         ! LOGG NAME
         ENCODE(4,102,filename_tmp(8:11)) iz_tmp        ! Z    NAME
         filename_tmp(8:8)   = '-'
         IF (atm.h.z .GT. 0.) filename_tmp(8:8) = 'p'	
         IF (atm.h.z .GE. -0.59) THEN
         filename_tmp(12:15) = '.dat'
         ELSE
         filename_tmp(12:21) = '_alpha.dat'
         ENDIF
         lf_output_filename  = filename_tmp(1:11)
       
101      FORMAT(I3)
102      FORMAT(I4.3)                                   ! print min. 4 Character
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *             'first looking for ',filename_tmp(1:21)
C ---    Check if file does exist in the grid directory
c MB 02.02.2011
c       
         env_len = igetenv('SIU_MAIN',env)
         atm_file_tmp = env(1:env_len) 
     *       //'/atmospheres/special/'// filename_tmp
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *             ' ',atm_file_tmp
         INQUIRE(FILE = atm_file_tmp,EXIST = exact_atm_exist)
         IF (exact_atm_exist) THEN
          atm_fn(1) = filename_tmp
          GOTO 200                             ! jump to the begin of ATM-loop
         ENDIF
         env_len = igetenv('SIU_MAIN',env)
         atm_file_tmp = env(1:env_len) 
     *       //'/atmospheres/grid/'// filename_tmp
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *             ' ',atm_file_tmp
         INQUIRE(FILE = atm_file_tmp,EXIST = exact_grid)
         IF (exact_grid) THEN
          atm_fn(1) = filename_tmp
          GOTO 200                             ! jump to the begin of ATM-loop
         ENDIF

       ENDIF
      ENDIF

C ---- searching atmosphere grid edges -----------------------------------------

      IF (NTEFFS.GT.1) THEN
       DO i=1,NTEFFS-1
        IF (atm.h.teff.GE.atm_teff(i).AND.
     *      atm.h.teff.LE.atm_teff(i+1)) THEN
          it = i
          x  = (atm.h.teff-atm_teff(it))/(atm_teff(it+1)-atm_teff(it))
          GOTO 5
        ENDIF
       ENDDO

c - extrapolation - be carefull !!!
       it = 0
       IF (atm.h.teff.LT.atm_teff(1))      it = 1
       IF (atm.h.teff.GE.atm_teff(NTEFFS)) it = NTEFFS-1
       IF (it.GE.1)  THEN
        x  = (atm.h.teff-atm_teff(it))/(atm_teff(it+1)-atm_teff(it))
        WRITE(CHN,*) 'W A R N I N G !!! TEFF EXTRAPOLATION: DTeff= ',
     *               atm.h.teff-atm_teff(it)
        lf_status = JIBSET(lf_status,1) 
        GOTO 5
       ENDIF

c - extrapolation - be carefull !!!
      ELSE
        IF (ABS(atm.h.teff-atm_teff(1)).LE.1.) THEN
          it = 1
          x  = 0.
          GOTO 5
        ENDIF        
      ENDIF
      WRITE(CHN,*) 'NO ATMOSPHERE GRID POINT (TEFF) FOUND ! FATAL !'
      lf_status = JIBSET(lf_status,0)         ! Fatal Error
      lf_status = JIBSET(lf_status,1) 
      RETURN
 5    CONTINUE

      IF (NLOGGS.GT.1) THEN
       DO i=1,NLOGGS-1
        IF (atm.h.logg.GE.atm_logg(i).AND.
     *      atm.h.logg.LE.atm_logg(i+1)) THEN
         ig = i 
         y  = (atm.h.logg-atm_logg(ig))/(atm_logg(ig+1)-atm_logg(ig))
         IF ((atm_logg(ig) - atm_logg(ig+1)) .GT. .8) THEN
          WRITE(CHN,*) 'DISTANCE OF GRID POINTS (LOGG) > 0.8 ! FATAL !'
          lf_status = JIBSET(lf_status,0)         ! Fatal Error
          lf_status = JIBSET(lf_status,2) 
         RETURN
         ENDIF
         GOTO 6
        ENDIF
       ENDDO
c - extrapolation - be carefull !!!
        ig = 0
        IF (atm.h.logg.LT.atm_logg(1))      ig = 1
        IF (atm.h.logg.GT.atm_logg(NLOGGS)) ig = NLOGGS-1 
        IF (ig.GE.1) THEN
         y  = (atm.h.logg-atm_logg(ig))/(atm_logg(ig+1)-atm_logg(ig))
         WRITE(CHN,*) 'W A R N I N G !!! LOG G EXTRAPOLATION: Dlog g= ',
     *               atm.h.logg-atm_logg(ig)
         lf_status = JIBSET(lf_status,2) 
         GOTO 6
        ENDIF
c - extrapolation - be carefull !!!
      ELSE
        IF (ABS(atm.h.logg-atm_logg(1)).LE.0.01 ) THEN
         ig = 1 
         y  = 0.
         GOTO 6
        ENDIF       
      ENDIF                     
      WRITE(CHN,*) 'NO ATMOSPHERE GRID POINT (LOGG) FOUND ! FATAL !'
      lf_status = JIBSET(lf_status,0)         ! Fatal Error
      lf_status = JIBSET(lf_status,2) 
      RETURN
 6    CONTINUE

      IF (NZS.GT.1) THEN
       DO i=1,NZS-1
        IF (atm.h.z.GE.atm_z(i).AND.
     *      atm.h.z.LE.atm_z(i+1)) THEN
         iz = i
         z  = (atm.h.z-atm_z(iz))/(atm_z(iz+1)-atm_z(iz))
         IF ((atm_z(iz) - atm_z(iz+1)) .GT. 0.5) THEN
          WRITE(CHN,*) 'DISTANCE OF GRID POINTS (Z) > 0.5 ! FATAL !'
          lf_status = JIBSET(lf_status,0)         ! Fatal Error
          lf_status = JIBSET(lf_status,3) 
          RETURN
         ENDIF
         GOTO 7
        ENDIF
       ENDDO
c - extrapolation - be carefull !!!
        iz = 0
        IF (atm.h.z.LT.atm_z(1))      iz = 1
        IF (atm.h.z.GT.atm_z(NZS)) iz = NZS-1 
        IF (iz.GE.1) THEN
         y  = (atm.h.z-atm_z(iz))/(atm_z(iz+1)-atm_z(iz))
         WRITE(CHN,*) 'W A R N I N G !!! LOG Z EXTRAPOLATION: Dlog z= ',
     *               atm.h.z-atm_z(iz)
         lf_status = JIBSET(lf_status,3) 
         GOTO 7
        ENDIF
c - extrapolation - be carefull !!!
      ELSE
        IF (ABS( atm.h.z - atm_z(1) ).LT.0.01 ) THEN
         iz = 1
         z  = 0.
         GOTO 7
        ENDIF
      ENDIF
      WRITE(CHN,*) 'NO ATMOSPHERE GRID POINT (Z) FOUND ! FATAL !'
      lf_status = JIBSET(lf_status,0)         ! Fatal Error
      lf_status = JIBSET(lf_status,3) 
      RETURN
 7    CONTINUE                

      DO i = 0, 7
         it_mod = it + MIN( MOD(i,2)  ,NTEFFS - 1 )
         ig_mod = ig + MIN( MOD(i/2,2),NLOGGS - 1 )
         iz_mod = iz + MIN( MOD(i/4,2),NZS    - 1 )

         atm_fn(i+1) = atm_filename(it_mod,ig_mod,iz_mod)
         IF (atm_fn(i+1).LE.' ') THEN
          WRITE(CHN,111) atm.h.teff,atm.h.logg,atm.h.z 
 111      FORMAT('GRID POINT MISSING: Teff=',F7.1,'; log(g)=',F7.3,
     *           '; z=',F7.3)
          lf_status = JIBSET(lf_status,0)         ! Fatal Error
          RETURN
         ENDIF

C --- Assign Model Indices to ATM_INT(I)

         atm_int(i+1).idx_teff = it_mod
         atm_int(i+1).idx_logg = ig_mod
         atm_int(i+1).idx_z    = iz_mod

      END DO

C --- Calculate interpolation coefficients
        
      co(1) = (1. - x) * (1. - y) * (1. - z)    
      co(2) =     x    * (1. - y) * (1. - z)    
      co(3) = (1. - x) *     y    * (1. - z)    
      co(4) =     x    *     y    * (1. - z)    
      co(5) = (1. - x) * (1. - y) *     z       
      co(6) =     x    * (1. - y) *     z       
      co(7) = (1. - x) *     y    *     z       
      co(8) =     x    *     y    *     z       
    
C --- Logarithmic linear interpolation between 8 model atmosphere files
C --- first search gridpoints in uag ---------------------------------
      atmgrid_exist   = .TRUE.
      DO k=1,8
        idx_atmgrid(k) = 0.     ! nescessary for assigning gridpoints later
      ENDDO
      DO k=1,8
        IF (.NOT.atmgrid_exist) GOTO 200
        DO i=1,N_ATMGRID_HOLD
          IF (atm_int(k).idx_teff .EQ.uag(i).idx_teff.AND.
     *        atm_int(k).idx_logg .EQ.uag(i).idx_logg.AND.
     *        atm_int(k).idx_z    .EQ.uag(i).idx_z)       THEN
            idx_atmgrid(k) = i
            GOTO 151  
          ENDIF
c        WRITE(CHN,*) 'GRID STACK: ',i,uag(i).idx_teff,uag(i).idx_logg,
c     *                  uag(i).idx_z
        ENDDO 
        atmgrid_exist = .FALSE.     
 151    CONTINUE
      ENDDO

C --- Provide gridpoints from uag
      IF (atmgrid_exist) THEN
        DO i=1,8      
          atm_int(i) = uag(idx_atmgrid(i))
        ENDDO

        ndepth = atm_int(1).ndepth      
        nel    = atm_int(1).nel        
        nmol   = atm_int(1).nmol      
        nbreak = atm_int(1).nbreak
        eta_cont_nr = atm_int(1).eta_cont_nr

        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *             'Atmosphere gridpoints assigned from stack'
        GOTO 300
      ENDIF
c------------------------------------------------------------------------------
c
c      READ MODEL ATMOSPHERE FILE
c
c------------------------------------------------------------------------------

200   CONTINUE

      DO  i = 1,8                                              ! start ATM-loop
        n_unit = 20 + i
        IF (exact_atm_exist) THEN     ! assume model atm in SPECIAL directory
          IF (BJTEST(lf_ctrl,22)) THEN
           env_len = IGETENV('ATMOSPHERE',env)
           IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *        'Opening ',env(1:env_len)
           OPEN(UNIT = n_unit, FILE=env(1:env_len),STATUS='OLD',
     *         READONLY, FORM = 'UNFORMATTED' )
          ELSE
           env_len = IGETENV('SIU_MAIN',env)
           OPEN(UNIT = n_unit, FILE=env(1:env_len)
     *         //'/atmospheres/special/'//atm_fn(i),
     *         STATUS='OLD', READONLY, FORM = 'UNFORMATTED' )
          ENDIF                      
        ELSE                          ! assume model atm interpolated
          env_len = IGETENV('SIU_MAIN',env)
          OPEN(UNIT = n_unit, FILE=env(1:env_len)
     *        //'/atmospheres/grid/'//atm_fn(i),
     *        STATUS='OLD', READONLY, FORM = 'UNFORMATTED' )
        ENDIF

        modtp = INDEX(env,'/', back=.true.)

       IF (env(modtp+1:modtp+7).EQ.'siu_bin') THEN                         ! MARCS 
        WRITE(*,*) 'MARCS-type model'
        READ(n_unit)
     *     ndepth,nel,nmol,nbreak,
     *    (atm_int(i).idx_el(n)    ,n=1,nel   ),
     *    (atm_int(i).nion(n)      ,n=1,nel   ),
     *    (atm_int(i).nscale(n)    ,n=1,nbreak),
     *    (atm_int(i).tau(n)       ,n=1,ndepth),
     *    (atm_int(i).t  (n)       ,n=1,ndepth),
     *    (atm_int(i).kappa_ref(n) ,n=1,ndepth),
     *    (atm_int(i).pe(n)        ,n=1,ndepth),
     *    (atm_int(i).pg(n)        ,n=1,ndepth),
     *    (atm_int(i).pk(n)        ,n=1,ndepth),
     *    (atm_int(i).rho(n)       ,n=1,ndepth),
     *  (((atm_int(i).pion(n,m,k),n=1,ndepth),m=1,
     *                               atm_int(i).nion(k)),k=1,nel),
     *  (((atm_int(i).uion(n,m,k),n=1,ndepth),m=1,
     *                               atm_int(i).nion(k)),k=1,nel)

        atm_int(i).ndepth      = ndepth 
        atm_int(i).nel         = nel
        atm_int(i).nmol        = nmol
        atm_int(i).nbreak      = nbreak

        READ(n_unit)  atm_int(i).h.date,atm_int(i).h.user,
     *                atm_int(i).h.ref_nr,atm_int(i).h.teff,
     *                atm_int(i).h.logg,atm_int(i).h.x,atm_int(i).h.y,
     *                atm_int(i).h.z,atm_int(i).h.alpha,
     *                atm_int(i).h.eps_nr,atm_int(i).h.status,
     *                atm_int(i).h.cmt_nr,
     *                atm_int(i).h.ref(:atm_int(i).h.ref_nr),
     *               (atm_int(i).h.eps(k),k=1,atm_int(i).h.eps_nr),
     *                atm_int(i).h.cmt(:atm_int(i).h.cmt_nr)
        WRITE(*,*)  'The model was constructed by: ',
     *               atm_int(i).h.user

c ---------------- CONTROL OF INPUT MODEL ---------------------
c
c        WRITE(*,"('N depth points',I4)") ndepth
        WRITE(*,"('N of elements',I4)") nel
        WRITE(*,"('N of molecules',I4)") nmol
        WRITE(*,"('N of break points in tau scale',I4)") nbreak
         WRITE(*,"('Atomic numbers of elements')")
         WRITE(*,*) (atm_int(i).idx_el(n), n=1, nel)
         WRITE(*,"('Number of ionization ',
     *                     'stages for each element')")
         WRITE(*,*) (atm_int(i).nion(n), n=1, nel)
         WRITE(*,"('Logarithmic element ',
     *            'abundances, A = log(N_el/N_H) + 12')")
         WRITE(*,"(10(F5.2,2X))") (atm_int(i).h.eps(k),
     *                           k=1,atm_int(i).h.eps_nr)
         WRITE(*,"('Break poins in tau-scale', 10I4)")
     *             (atm_int(i).nscale(n), n=1, nbreak)
c
        DO n=1, ndepth
         WRITE(*,"(I4,14F13.6)") n, atm_int(i).tau(n),
     *                 10**(atm_int(i).t(n)), atm_int(i).kappa_ref(n),
     *                 atm_int(i).pe(n), atm_int(i).pg(n),
     *                 atm_int(i).pk(n), atm_int(i).rho(n)
        ENDDO
c
c        DO n=1, ndepth
c           WRITE(*,"(I4,300F13.6)")
c     *                 n,((atm_int(i).pion(n,m,k),
c     *                 m=1,atm_int(i).nion(k)), k=1,nel)
c        ENDDO
c
c        DO n=1, ndepth
c           WRITE(*,"(I4,300F13.6)")
c     *                 n,((atm_int(i).uion(n,m,k),
c     *                 m=1,atm_int(i).nion(k)), k=1,nel)
c        ENDDO

c ----------------
        ELSE           !MAFAGS models WITH molecular part functions
        WRITE(*,*) 'MAFAGS-type model'
        READ(n_unit)
     *     ndepth,nel,nmol,nbreak,
     *    (atm_int(i).idx_el(n)    ,n=1,nel   ),
     *    (atm_int(i).nion(n)      ,n=1,nel   ),
     *    (atm_int(i).nscale(n)    ,n=1,nbreak),
     *    (atm_int(i).tau(n)       ,n=1,ndepth),
     *    (atm_int(i).t  (n)       ,n=1,ndepth),
     *    (atm_int(i).kappa_ref(n) ,n=1,ndepth),
     *    (atm_int(i).pe(n)        ,n=1,ndepth),
     *    (atm_int(i).pg(n)        ,n=1,ndepth),
     *    (atm_int(i).pk(n)        ,n=1,ndepth),
     *    (atm_int(i).rho(n)       ,n=1,ndepth),
     *   ((atm_int(i).pmol(n,k)  ,n=1,ndepth),k=1,nmol),
     *   ((atm_int(i).umol(n,k)  ,n=1,ndepth),k=1,nmol),
     *  (((atm_int(i).pion(n,m,k),n=1,ndepth),m=1,
     *                               atm_int(i).nion(k)),k=1,nel),
     *  (((atm_int(i).uion(n,m,k),n=1,ndepth),m=1,
     *                               atm_int(i).nion(k)),k=1,nel)

        atm_int(i).ndepth      = ndepth 
        atm_int(i).nel         = nel
        atm_int(i).nmol        = nmol
        atm_int(i).nbreak      = nbreak

        READ(n_unit)  atm_int(i).h.date,atm_int(i).h.user,
     *                atm_int(i).h.ref_nr,atm_int(i).h.teff,
     *                atm_int(i).h.logg,atm_int(i).h.x,atm_int(i).h.y,
     *                atm_int(i).h.z,atm_int(i).h.alpha,
     *                atm_int(i).h.eps_nr,atm_int(i).h.status,
     *                atm_int(i).h.cmt_nr,
     *                atm_int(i).h.ref(:atm_int(i).h.ref_nr),
     *               (atm_int(i).h.eps(k),k=1,atm_int(i).h.eps_nr),
     *                atm_int(i).h.cmt(:atm_int(i).h.cmt_nr)
         

c        DO n=1, ndepth
c         WRITE(*,"(I4,14F13.6)") n, atm_int(i).tau(n),
c     *                 10**(atm_int(i).t(n)), atm_int(i).kappa_ref(n),
c     *                 atm_int(i).pe(n), atm_int(i).pg(n),
c     *                 atm_int(i).pk(n), atm_int(i).rho(n)
c        ENDDO
c

        ENDIF

        CLOSE(n_unit)
                                                           
        IF (exact_atm_exist.OR.exact_grid) THEN             ! exact atmosphere exist

          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *            'Exact Atmosphere found - no interpolation'
         
          atm.h.teff = atm_int(i).h.teff                      
          atm.h.logg = atm_int(i).h.logg                      
          atm.h.z    = atm_int(i).h.z                         

          atm.h.x    = atm_int(i).h.x                         
          atm.h.y    = atm_int(i).h.y                         
          atm.h.alpha= atm_int(i).h.alpha                     
          atm.h.eps_nr = atm_int(i).h.eps_nr
          DO k=1,atm.h.eps_nr 
              atm.h.eps(k) = atm_int(i).h.eps(k)
          ENDDO
          atm.h.ref_nr = atm_int(i).h.ref_nr                      
          atm.h.ref    = atm_int(i).h.ref                      
          atm.h.cmt    = atm_int(i).h.cmt

C --------- providing spectrum header ----------------------------
          head.atm     = atm_int(i).h

          DO n = 1, ndepth
            atm.logtau(n)   =        atm_int(i).tau(n)
            atm.tau(n)      = 10. ** atm_int(i).tau(n)
            atm.t(n)        = 10. ** atm_int(i).t(n)
            atm.kappa_ref(n)= 10. ** atm_int(i).kappa_ref(n)
            atm.pe(n)       = 10. ** atm_int(i).pe(n)
            atm.pk(n)       = 10. ** atm_int(i).pk(n)
            atm.pg(n)       = 10. ** atm_int(i).pg(n)
            atm.rho(n)      = 10. ** atm_int(i).rho(n)
          END DO

c------- MB 12.2010------------

          IF (nmol.GE.1) THEN        ! MAFAGS models with full data
            DO k = 1,nmol                        
             DO n = 1, ndepth
              atm.pmol(n,k)   = 10. ** atm_int(i).pmol(n,k)
              IF (atm.pmol(n,k).LT.EEEPS) atm.pmol(n,k) = 0.              
              atm.umol(n,k)   = 10. ** atm_int(i).umol(n,k)
             END DO
             atm.h.cmol(k) = lbl(NATOM+k)
            END DO
          ELSE IF (nmol.EQ.0) THEN   ! MARCS do not have pmol and umol --> set them to 0
            nmol = 14
            DO k = 1, nmol  ! set the max number of molecules - 14 (molec.ini)
              DO n = 1, ndepth
              atm.pmol(n,k) = 0.
              atm.umol(n,k) = 0.
              END DO
              atm.h.cmol(k) = lbl(NATOM+k)
             END DO
          ENDIF 

       PRINT*, 'Number of molecules', nmol
       PRINT*, 'PP for molecules? pmol(nd,nmol)',atm.pmol(ndepth,nmol)
       PRINT*, 'PF for molecules? umol(nd,nmol)',atm.umol(ndepth,nmol)
       PRINT*, 'NR_ION - number of ion stages', NR_ION

          DO m = 1,NR_ION             ! go from 3D array pion(x,y,z) to 2D
           DO n = 1,ndepth            ! new variable atm.pp - used in kappa.f for
              atm.pp(n,m) = 0.        !    free-free transitions  
           ENDDO
          ENDDO

          DO k = 1, nel                          ! index number
           atomic_nr = atm_int(1).idx_el(k)      ! atomic number (Fe 26)
           DO m = 1, atm_int(1).nion(k)          ! ionization stage, e.g. 1
            DO n = 1, ndepth
             atm.pion(n,m,atomic_nr) = 10.**atm_int(i).pion(n,m,k)
             IF (atm.pion(n,m,atomic_nr).LT.EEEPS)
     *           atm.pion(n,m,atomic_nr) = 0.              
             atm.uion(n,m,atomic_nr) = 10.**atm_int(i).uion(n,m,k)
             atm.pp(n,m) = atm.pp(n,m) + atm.pion(n,m,atomic_nr)  !for all el-s of each ion stage 
            END DO
           END DO
           atm.h.nion(k)   = atm_int(1).nion(k)
           atm.h.idx_el(k) = atm_int(1).idx_el(k)  ! atomic number
           atm.h.cel(k)    = lbl(atomic_nr)        ! label: 26 --> Fe
          END DO

         atm.h.nbreak      = nbreak
         DO k=1,atm.h.nbreak
            atm.h.nscale(k)   = atm_int(1).nscale(k)
         ENDDO
          atm.h.ndepth      = ndepth
          atm.h.nel         = nel
          atm.h.nmol        = nmol

C --- storing name of of atmosphere into atm reference list

          IF (BJTEST(lf_ctrl,22)) THEN  ! if atm specified by ATMOSPHERE variable
           atm.h.ref    = env(1:env_len) // ' // ' //
     *                    atm.h.ref(1:atm.h.ref_nr)
           atm.h.ref_nr = atm.h.ref_nr + 4 + env_len
          ELSE
             atm.h.ref    = atm_fn(1) // ' // ' //
     *                    atm.h.ref(1:atm.h.ref_nr)
             atm.h.ref_nr = atm.h.ref_nr + 4 + LEN(atm_fn(1))
          ENDIF                      

C --- add atmosphere to atm stack -------------------------

          ua_ptr = ua_ptr + 1
          IF (ua_ptr.GT.N_ATM_HOLD) ua_ptr = 1
          ua(ua_ptr) = atm

          GOTO 8000 ! jump to end - finish model input for an EXACT model
        ENDIF
      ENDDO   ! 8 models 

C --- add atmosphere gridpoints to atmgrid stack -------------------------

       DO i=1,8
        uag_ptr = uag_ptr + 1
        IF (uag_ptr.GT.N_ATMGRID_HOLD) uag_ptr = 1
        IF (idx_atmgrid(i).LE.0) uag(uag_ptr) = atm_int(i)
       ENDDO

C --------------------------------------------------------

      atm.h.ref_nr = atm_int(i).h.ref_nr  + 4
      atm.h.ref    = ' // ' // atm_int(i).h.ref
      DO i=1,8
        atm.h.ref    = atm_fn(i) // ',' // atm.h.ref(1:atm.h.ref_nr)
        atm.h.ref_nr = atm.h.ref_nr + 1 + LEN(atm_fn(i))
      ENDDO

 300  CONTINUE        

C-----------------------------------------------------------------------
C
C    Interpolation of data arrays (NDMAX) and (!) exponentation
C
C ----------------------------------------------------------------------
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                   'Atmosphere - Interpolation'

      DO  n = 1, ndepth
        atm.tau(n)        = 0.
        atm.t(n)          = 0.
        atm.kappa_ref(n)  = 0.
        atm.pe(n)         = 0.
        atm.pk(n)         = 0.
        atm.pg(n)         = 0.
        atm.rho(n)        = 0.
        DO i = 1,8
          atm.tau(n) = atm.tau(n) + co(i) * atm_int(i).tau(n)
          atm.t  (n) = atm.t  (n) + co(i) * atm_int(i).t  (n)
          atm.kappa_ref(n) = 
     *           atm.kappa_ref(n) + co(i) * atm_int(i).kappa_ref(n)
          atm.pe (n) = atm.pe (n) + co(i) * atm_int(i).pe (n)
          atm.pk (n) = atm.pk (n) + co(i) * atm_int(i).pk (n)
          atm.pg (n) = atm.pg (n) + co(i) * atm_int(i).pg (n)
          atm.rho(n) = atm.rho(n) + co(i) * atm_int(i).rho(n)
        ENDDO
        atm.logtau(n)     =     atm.tau(n)
        atm.tau(n)        = 10.**atm.tau(n)
        atm.t(n)          = 10.**atm.t  (n)
        atm.pe(n)         = 10.**atm.pe (n)
        atm.pk(n)         = 10.**atm.pk (n)
        atm.pg(n)         = 10.**atm.pg (n)
        atm.rho(n)        = 10.**atm.rho(n)
        atm.kappa_ref(n)  = 10.**atm.kappa_ref(n)
      END DO

      DO k = 1,nmol
        atm.h.cmol(k)  = lbl(NATOM+k)
        DO n = 1, ndepth
          sum1 = 0.
          sum2 = 0.
          DO ic = 1,8
            sum1 = sum1 + co(ic) * atm_int(ic).pmol(n,k)
            sum2 = sum2 + co(ic) * atm_int(ic).umol(n,k)
          ENDDO
          atm.pmol(n,k)   =    10 ** sum1
          IF (atm.pmol(n,k).LT.EEEPS) atm.pmol(n,k) = 0.              
          atm.umol(n,k)   =    10 ** sum2
        ENDDO
      ENDDO

      DO k = 1,nel
        atomic_nr =atm_int(1).idx_el(k)
        atm.h.cel(k)    = lbl(atomic_nr)
        atm.h.nion(k)   = atm_int(1).nion(k)
        atm.h.idx_el(k) = atm_int(1).idx_el(k)
        DO m = 1,atm.h.nion(k)
          DO n = 1, ndepth
            sum1 = 0.
            sum2 = 0.
            DO ic = 1,8
              sum1 = sum1 + co(ic) * atm_int(ic).pion(n,m,k)
              sum2 = sum2 + co(ic) * atm_int(ic).uion(n,m,k)
            ENDDO
            atm.pion(n,m,atomic_nr) = 10 ** sum1
             IF (atm.pion(n,m,atomic_nr).LT.EEEPS)
     *           atm.pion(n,m,atomic_nr) = 0.              
            atm.uion(n,m,atomic_nr) = 10 ** sum2
          ENDDO
        ENDDO
      ENDDO

      atm.h.nbreak      = nbreak
      DO i=1,nbreak
        atm.h.nscale(i) = atm_int(1).nscale(i)
      ENDDO
      atm.h.ndepth      = ndepth
      atm.h.nel         = nel
      atm.h.nmol        = nmol

      i = 1                                       ! only first header
      atm.h.date   = atm_int(i).h.date                      
      atm.h.user   = atm_int(i).h.user                      
      atm.h.x      = atm_int(i).h.x                         
      atm.h.y      = atm_int(i).h.y                         
      atm.h.alpha  = atm_int(i).h.alpha                     
      atm.h.eps_nr = atm_int(i).h.eps_nr
      atm.h.cmt_nr = atm_int(i).h.cmt_nr                      
      IF (i.EQ.1) THEN
            atm.h.cmt = atm_int(1).h.cmt(:atm.h.cmt_nr)// 
     *                        ' Grid interpolation'
            atm.h.cmt_nr = atm.h.cmt_nr + 19
            atm.h.status = IIBSET(atm.h.status,0)
      ENDIF        

C --- Interpolation of abundances
      DO k=1,atm.h.eps_nr
         atm.h.eps(k) = 0.
         DO n=1,8
          atm.h.eps(k) = atm.h.eps(k) + co(n) * atm_int(n).h.eps(k)
         ENDDO
      ENDDO

C --- add atmosphere to atm stack -------------------------
      ua_ptr = ua_ptr + 1
      IF (ua_ptr.GT.N_ATM_HOLD) ua_ptr = 1
      ua(ua_ptr) = atm
C --- -----------------------------------------------------
 8000 CONTINUE

C --- assign header information

      head.atm.teff   = atm.h.teff
      head.atm.logg   = atm.h.logg
      head.atm.z      = atm.h.z   

      head.atm.date   = atm.h.date
      head.atm.user   = atm.h.user
      head.atm.x      = atm.h.x   
      head.atm.y      = atm.h.y     
      head.atm.alpha  = atm.h.alpha
      head.atm.eps_nr = atm.h.eps_nr
      head.atm.cmt_nr = atm.h.cmt_nr
      head.atm.ref_nr = atm.h.ref_nr
      head.atm.cmt    = atm.h.cmt
      head.atm.ref    = atm.h.ref
      DO k=1,head.atm.eps_nr
         head.atm.eps(k) = atm.h.eps(k)
      ENDDO 

      RETURN
      END
