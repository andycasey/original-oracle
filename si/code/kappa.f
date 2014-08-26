      SUBROUTINE KAPPA(read_bf,wlam,wmin,wmax)

C --- Partial pressures  and T structure in 'atm'
C --- wmin > 0.: Result placed in atm.eta_cont_w,atm.eta_sig and atm.eta_cont
C ---            atm.eta_cont_w(1) = wstart and atm.eta_cont_w(2) = wend
C ---            wlam = minimum distance of two continuum points in v/c 
C --- wmin  .LE. 0. : Result placed in atm.eta_cont_w(1),atm.eta_sig(1) and atm.eta_cont(1)
C --- read_bf  : reread the bf data file if .TRUE.

C --- At the first call, eta_sig and eta_cont at 5000A is calculated additionally
C --- If kappa-5000 is inconsistent with kappa-ref(=5000!) from the atmosphere
C     by more than 1% you will be informed. Be aware that this is 
c --- may by a FATAL inconstistency.

C-----INPUT OF BOUND-FREE CROSS SECTION- AND TRANSITION-DATA FOR METALS
C     uphoto.asc contains photoionisation cross-sections and related links 
c
      INCLUDE 'lf_param.inc'
      INCLUDE 'atm_def.inc'
      INCLUDE 'mol_data.inc'
      INCLUDE 'lf_info.inc'

      CHARACTER*80  REC_FMT
      CHARACTER*160 INFO
      PARAMETER(REC_LEN=80,REC_FMT='(A80)')
      PARAMETER(METAL_MAX=300,METAL_MAX_VAL=3,MAX_WAVE=MAX_ETA_CONT)
      LOGICAL hminbf,h2plus,h2quasi,hebf,metbf,chbf,ohbf
      LOGICAL read_bf,gen_new_tau_scale,show_tau_one
      REAL*8 w,w2,wlam,wmin,wmax,wmax1
      REAL*4 xsig(NDMAX),xkap(NDMAX),pp(NDMAX,2),xkaptmp(NDMAX),ytmp
      REAL*4 xtau(NDMAX),fc_lam(NDMAX),dxtau(NDMAX),xta,xdiff
      INTEGER*4 ircnt,idx_tmp(MAX_WAVE)
 
      STRUCTURE /OPTION_TYP/
        LOGICAL h1,hm,hp,h2m,h2,hmff,h2mff,heff,
     *          h2p,he1,he2,metal,metff,chbf,ohbf,
     *          link,rayleigh,thompson,fac 
      END STRUCTURE

      STRUCTURE /METAL_TYP/
        REAL*4       w0
        INTEGER*2    id
        INTEGER*2    ion
        CHARACTER*10 desig
        REAL*4       elow
        REAL*4       g
        REAL*4       s
        INTEGER*2    maxval
        REAL*4       val(METAL_MAX_VAL)
      END STRUCTURE

      STRUCTURE /WAVE_TYP/    ! Grid with Edges to be definitly taken into account
        INTEGER*4   maxwave
        REAL*8      w(MAX_WAVE)
      ENDSTRUCTURE
       
      STRUCTURE /BF_TYP/
        INTEGER*2 metmax
        INTEGER*2 lnkmax
        REAL*4    zwmax
        REAL*4    bhmin(NDMAX)       ! departure factor
        REAL*4    wmin_fac
        REAL*4    wmax_fac
        REAL*4    fac(NDMAX)         ! additional factor between lammin,lammax 
        RECORD /OPTION_TYP/ with
        RECORD /METAL_TYP/  z(METAL_MAX)
        RECORD /WAVE_TYP/   edge
      ENDSTRUCTURE
       
      RECORD /BF_TYP/ bf

      INTEGER*4 lf_ctrl,lf_status
      COMMON/CONTROL/  lf_ctrl, lf_status
      COMMON /KAPPA_WORK/ bf,xsig,xkap,pp,xkaptmp,
     *        xtau,fc_lam,dxtau,xta

C --- FOR ff crossections

      REAL*4 XMICLAM,XMICLAM2,
     *  HFIT1(12),HFIT2(12),HFIT3(12),HFIT4(12),HFIT5(12),
     *  HFIT6(12),HEFIT1(5,4),HEFIT2(5,4),H2FIT1(5,4),
     *  H2FIT2(5,4)     

C-----FREE-FREE TRANSITION FOR NEGATIVE IONS HMINUS, HE AND H2
C-----(STIMULATED EMISSION ALREADY INCLUDED)

      DATA HFIT1/
     *     518.1021,473.2636,-482.2089,115.5291,0.,0.,
     *     0.,2483.346,-3449.889,2200.04,-696.271,88.283/
      DATA HFIT2/
     *     -734.8666,1443.4137,-737.1616,169.6374,0.,0.,
     *     0.,285.827,-1158.382,2427.719,-1841.4,444.517/
      DATA HFIT3/
     *     1021.1775,-1977.3395,1096.8827,-245.649,0.,0.,
     *     0.,-2054.291,8746.523,-13651.105,8624.97,-1863.864/
      DATA HFIT4/
     *     -479.0721,922.3575,-521.1341,114.243,0.,0.,
     *     0.,2827.776,-11485.632,16755.524,-10051.53,2095.288/
      DATA HFIT5/
     *     93.1373,-178.9275,101.7963,-21.9972,0.,0.,
     *     0.,-1341.537,5303.609,-7510.494,4400.067,-901.788/
      DATA HFIT6/
     *     -6.4285,12.36,-7.0571,1.5097,0.,0.,
     *     0.,208.952,-812.939,1132.738,-655.02,132.985/
      DATA HEFIT1/
     *     1.0933E-26,-2.5920E-21,4.9416E-16,6.8535E-13,3.4565E-08,
     *     -1.7515E-29,4.4136E-24,-9.1140E-19,-6.5192E-16,-1.9573E-11,
     *     7.3266E-33,-1.9652E-27,8.5717E-22,-2.5374E-18,3.0173E-14,
     *     -7.7148E-37,2.1065E-31,4.0995E-26,3.9213E-22,-2.7214E-18/
      DATA HEFIT2/
     *     -3.5558E-22,1.0612E-17,-1.1670E-13,5.6915E-10,-1.0012E-06, 
     *     3.9255E-25,-1.1970E-20,1.3401E-16,-6.6289E-13,1.1683E-09,
     *     -9.1739E-29,2.9436E-24,-3.4054E-20,1.7747E-16,-3.0801E-13, 
     *     1.1345E-32,-3.5124E-28,4.1100E-24,-2.0387E-20,3.6481E-17/
      DATA H2FIT1/
     *     -7.8332E-26,2.0112E-20,7.6805E-16,4.5822E-11,-3.5157E-07,
     *     8.1734E-29,-2.0809E-23,-2.0524E-18,-4.7964E-14,4.2256E-10, 
     *     -2.3488E-32,6.0786E-27,1.6335E-21,1.5977E-17,-9.1062E-14,
     *     1.5808E-36,-4.1194E-31,1.7169E-25,-1.1072E-21,6.3797E-18/
      DATA H2FIT2/
     *     -2.6829E-23,7.1315E-19,-2.3389E-15,-5.1760E-13,8.8265E-09, 
     *     2.7483E-26,-7.1468E-22,3.2818E-19,1.3923E-14,-3.0805E-11,
     *     7.9826E-30,-2.5980E-25,5.5306E-21,-1.6901E-17,3.6529E-14,
     *     2.8803E-33,-8.7588E-29,1.0659E-24,-4.2094E-21,6.0402E-18/

C --- for scattering coefficients

      REAL EDGE(3),SIGCO(3,3)
      DATA EDGE/ 1215.7E-8 , 584.5E-8 , 1108.6E-8/
  
      DATA SIGCO/                         ! for w*1.e8  [A]
     *     1.162E-12,5.490E-14,8.779E-13,
     *     2.452E+06,4.476E+05,1.507E+06,
     *     4.801E+12,1.688E+11,2.557E+12/

      CHARACTER*(REC_LEN) rec
      CHARACTER cdum

      INCLUDE 'physcnst.inc'                    ! physical constants

      wmax1 = wmax * 1.2D0
      
      IF (read_bf) THEN
       bf.zwmax         = 0.
       bf.metmax        = 0
       bf.lnkmax        = 0
       bf.with.h1       = .FALSE. 
       bf.with.hm       = .FALSE. 
       bf.with.hp       = .FALSE. 
       bf.with.h2m      = .FALSE. 
       bf.with.h2       = .FALSE. 
       bf.with.h2p      = .FALSE. 
       bf.with.he1      = .FALSE. 
       bf.with.he2      = .FALSE. 
       bf.with.hmff     = .FALSE. 
       bf.with.h2mff    = .FALSE. 
       bf.with.heff     = .FALSE. 
       bf.with.chbf     = .FALSE.
       bf.with.ohbf     = .FALSE.
       bf.with.metff    = .FALSE.
       bf.with.metal    = .FALSE.
       bf.with.link     = .FALSE.
       bf.with.rayleigh = .FALSE. 
       bf.with.thompson = .FALSE.
       bf.with.fac      = .FALSE. 
       bf.edge.maxwave  = 0

       DO n=1,NDMAX       
          bf.bhmin(n)   = 1.     ! H minus departure factor
          bf.fac(n)     = 1.     ! H minus departure factor
       ENDDO

       ircnt = 0
c       env_len = IGETENV('LF_INPATH',env)
c       OPEN(UNIT=24,FILE=env(1:env_len)//'/'//RBF_FILE,
c     *      STATUS = 'OLD',READONLY)
c       WRITE(*,*) 'Background opacity taken from: ', 
c     *             env(1:env_len)//'/'//RBF_FILE
c
       OPEN(UNIT=24,FILE='/data/mbergema/siu/ewgrid/code/f/'
     *      //RBF_FILE,
     *      STATUS = 'OLD',READONLY)
       WRITE(*,*) 'Background opacity taken from: ', 
     *             RBF_FILE
c
 1     CONTINUE
       READ(24,REC_FMT,END=999) rec
       ircnt = ircnt + 1
       IF (rec.LE.'     '.OR.rec(1:1).EQ.'#') GOTO 1
       IF (INDEX(REC,'\OPAC').GT.0) THEN
         IF (INDEX(rec,'H1') .GT.0)      bf.with.h1       = .TRUE. 
         IF (INDEX(rec,'H-bf').GT.0)     bf.with.hm       = .TRUE. 
         IF (INDEX(rec,'Hmbf').GT.0)     bf.with.hm       = .TRUE. 
         IF (INDEX(rec,'H+') .GT.0)      bf.with.hp       = .TRUE. 
         IF (INDEX(rec,'H2-bf').GT.0)    bf.with.h2m      = .TRUE. 
         IF (INDEX(rec,'H2mbf').GT.0)    bf.with.h2m      = .TRUE. 
         IF (INDEX(rec,'H2Q').GT.0)      bf.with.h2       = .TRUE. 
         IF (INDEX(rec,'H2+').GT.0)      bf.with.h2p      = .TRUE. 
         IF (INDEX(rec,'He1').GT.0)      bf.with.he1      = .TRUE. 
         IF (INDEX(rec,'He2').GT.0)      bf.with.he2      = .TRUE. 
         IF (INDEX(rec,'H-ff').GT.0)     bf.with.hmff     = .TRUE. 
         IF (INDEX(rec,'Hmff').GT.0)     bf.with.hmff     = .TRUE. 
         IF (INDEX(rec,'H2mff').GT.0)    bf.with.h2mff    = .TRUE. 
         IF (INDEX(rec,'H2-ff').GT.0)    bf.with.h2mff    = .TRUE. 
         IF (INDEX(rec,'Heff').GT.0)     bf.with.heff     = .TRUE. 
         IF (INDEX(rec,'METAL').GT.0)    bf.with.metal    = .TRUE.
         IF (INDEX(rec,'METFF').GT.0)    bf.with.metff    = .TRUE.
         IF (INDEX(rec,'CHBF').GT.0)     bf.with.chbf     = .TRUE.
         IF (INDEX(rec,'OHBF').GT.0)     bf.with.ohbf     = .TRUE.
         IF (INDEX(rec,'RAY').GT.0)      bf.with.rayleigh = .TRUE. 
         IF (INDEX(rec,'THOMP').GT.0)    bf.with.thompson = .TRUE.  
         IF (INDEX(rec,'LINK') .GT.0)    bf.with.link     = .TRUE.
         GOTO 1
       ELSE IF (INDEX(REC,'\HMINB').GT.0) THEN   ! departure coeff for Hminus
          READ(rec(7:REC_LEN),*) ndepth          ! tau(1)<tau(n)
          IF (ndepth.EQ. 1) THEN
            READ(24,*) bf.bhmin(1)               ! tau(1)<tau(n)
            DO n=2,NDMAX       
              bf.bhmin(n) = bf.bhmin(1)          ! H minus departure factor
            ENDDO
          ELSE            
            DO n = 1,ndepth
              READ(24,*) bf.bhmin(n)             ! tau(1)<tau(n)
            ENDDO
          ENDIF
          GOTO 1 
       ELSE IF (INDEX(REC,'\SCALE').GT.0) THEN   ! departure coeff for Hminus
          READ(rec(7:REC_LEN),*) ndepth,bf.wmin_fac,bf.wmax_fac ! tau(1)<tau(n)
          IF (ndepth.EQ. 1) THEN
            READ(24,*) bf.fac(1)                 ! tau(1)<tau(n)
            DO n=2,NDMAX       
              bf.fac(n) = bf.fac(1)              ! H minus departure factor
            ENDDO
          ELSE            
            DO n = 1,ndepth
              READ(24,*) bf.fac(n)               ! tau(1)<tau(n)
            ENDDO
          ENDIF
          bf.with.fac  = .TRUE.  
          GOTO 1 
       ELSE IF (INDEX(rec,'\METAL').GT.0) THEN   ! Photo. cross-sections of metals
c all hard-coded in uphoto.asc
          i = 0
 20       READ(24,REC_FMT,END=999) rec
          ircnt = ircnt + 1
          IF (rec.LE.'     '.OR.rec(1:1).EQ.'#') GOTO 20
          IF (INDEX(REC,'\END ').GT.0) GOTO 1
          i = i + 1
          READ(rec,2402) bf.z(i).w0,bf.z(i).desig,bf.z(i).id,
     *       bf.z(i).ion,bf.z(i).elow,bf.z(i).g,
     *       bf.z(i).s,(bf.z(i).val(k),k=1,METAL_MAX_VAL)
 2402     FORMAT(F9.2,A9,2I5,F7.3,F5.0,F7.2,E11.4,2F8.3)
          bf.z(i).maxval=1
          bf.z(i).w0    = bf.z(i).w0*1.E-8   ! A --> cm
          IF (bf.z(i).w0 .GT. bf.zwmax) bf.zwmax = bf.z(i).w0 
          DO k=2,METAL_MAX_VAL
           IF (bf.z(i).val(k).NE.0.) bf.z(i).maxval = bf.z(i).maxval+1
          ENDDO
          bf.metmax = i
          bf.z(i).elow = EV_K * bf.z(i).elow  ! * e/k
          a1 = 1.
          xs = 0.
          IF ( bf.z(i).maxval .GE. 2) THEN
            xs = xs + 1.
            a1 = bf.z(i).val(2) - 1.
            IF ( bf.z(i).maxval .GE. 3) THEN
               xs = xs + 1.
               a1 = a1 - bf.z(i).val(3)
               bf.z(i).val(2) = bf.z(i).val(2) - 2. * bf.z(i).val(3) 
               bf.z(i).val(3) = bf.z(i).val(3) * bf.z(i).w0 ** 2 / a1
            ENDIF
            bf.z(i).val(2) = bf.z(i).val(2) * bf.z(i).w0 / a1    
          ENDIF
          bf.z(i).val(1) = 
     *       bf.z(i).g * bf.z(i).val(1) * a1 / bf.z(i).w0 ** xs
         GOTO 20          
       ELSE IF (INDEX(REC,'\LINK').GT.0) THEN  ! links to OP Photoionisationscrossections (planned)
c        NOT available yet
 30      READ(24,REC_FMT,END=999) rec
         ircnt = ircnt + 1
         IF (rec.LE.'     '.OR.rec(1:1).EQ.'#') GOTO 30
         IF (INDEX(REC,'\END ').GT.0) GOTO 1
         bf.lnkmax = i
C        rec        ! links
         GOTO 30
       ELSE IF (INDEX(rec,'\WAVE').GT.0) THEN ! b-f edges or H, H-, etc...
         i = bf.edge.maxwave
 40      READ(24,REC_FMT,END=999) rec
         ircnt = ircnt + 1
         IF (rec.LE.'     '.OR.rec(1:1).EQ.'#') GOTO 40
         IF (INDEX(REC,'\END ').GT.0) THEN
           GOTO 1
         ELSE IF (INDEX(rec,'\ADD').GT.0) THEN 
           rec = rec((INDEX(rec,'\ADD')+4): REC_LEN)
 42        CONTINUE
             kk = INDEX(rec,',')
             IF (kk.LE.1) GOTO 40
             i = i + 1
             READ(rec(1:kk-1),*) xxx 
             bf.edge.w(i) = xxx  * 1.E-8
             rec = rec(kk+1:REC_LEN)
             bf.edge.maxwave=i
           GOTO 42
         ELSE IF (INDEX(rec,'\Z_EDGE').GT.0) THEN            
            READ(rec((INDEX(rec,'\Z_EDGE')+7): REC_LEN),*) dist
            dist = dist * 1.E-08
            DO k = 1,bf.metmax
              bf.edge.w(i+1) = bf.z(k).w0 - dist
              bf.edge.w(i+2) = bf.z(k).w0 + dist
              i = i + 2
            ENDDO
            bf.edge.maxwave=i
            GOTO 40
         ELSE
          READ(rec,*) xxx,dist
          bf.edge.w(i+1) = (xxx - dist) * 1.e-8          
          bf.edge.w(i+2) = (xxx + dist) * 1.e-8          
          i = i + 2
          bf.edge.maxwave=i
          GOTO 40
         ENDIF
       ELSE
         GOTO 1
       ENDIF
       GOTO 1 
 999   CLOSE(24)            
       IF (bf.metmax.LE.0) bf.with.metal = .FALSE.
       IF (bf.lnkmax.LE.0) bf.with.link  = .FALSE.
       CALL SORT(bf.edge.maxwave,bf.edge.w,idx_tmp)
      ENDIF   !    if read_bf .EQ. .TRUE. !

C---------- ATMOSPHERIC STRATIFICATION LOOP -----------

    
 500  IF (read_bf) THEN
        w = 0.5D-4        ! calculate kappa_5000A for reference
        GOTO 501
      ENDIF 

      nr_edges = 0
      iwcnt    = 1
      IF (wmin.GT.0.D0) THEN
        w = wmin
        DO i=1,bf.edge.maxwave
           IF (bf.edge.w(i) .GE. wmax1) GOTO 601
           IF (bf.edge.w(i) .GT. wmin) THEN
              nr_edges = nr_edges + 1    
           ELSE
              iwcnt = i
           ENDIF
        ENDDO
 601    CONTINUE
        iwstart  = iwcnt + 1
        iwend    = iwcnt + nr_edges
        iwcnt    = iwstart
        wlam     = DMAX1( wlam,
     *    LOG(wmax/wmin) /(MAX_ETA_CONT-nr_edges-3))
      ELSE IF (wmin.LE.0.D0) THEN
        w = wlam
      ELSE
        STOP 951
      ENDIF
      il = 0
 501  CONTINUE                                   ! wavelength loop
      il = il + 1                                ! wavelength counter
c usually 3 frequency points across a line

C --------------------------------------------------------------
      hminbf  = bf.with.hm .AND.
     *           ( w .GE. 1250.D-8) .AND. ( w .LE. 16419.D-8)
      h2plus  = bf.with.h2p .AND. 
     *           ( w .GT. 500.D-8 ) .AND. ( w .LT. 10000.D-8 )
      h2quasi = bf.with.h2  .AND.  ( w .GE. 1400.D-8) .AND. 
     *           ( w .LE. 10000.D-8 ) .AND. ( t .LE. 12000. )
      hebf    = bf.with.he1   .AND. W .LE. 3679.0D-8
      metbf   = bf.with.metal .AND. W. LE. bf.zwmax
      chbf    = bf.with.chbf  .AND. W. LE. 7000.D-08
      ohbf    = bf.with.ohbf  .AND. W. LE. 7000.D-08

C --- ----------------------------------------------------------

      DO n = 1,atm.h.ndepth                      ! begin stratification loop
       t      = atm.t(n)
       divt   = 1. / t
       div_kt = DIV_K * divt
       EMI = EXP(- HC_DIV_K/(t*w))
       EM = 1. - EMI 
       akap   = 0.
       aminff = 0.               
       sig    = 0.

C --- BOUND-FREE TRANSITIONS OF HMINUS
       IF ( hminbf )
     *    akap = akap + AHMIN(w) * atm.pion(n,3,1) * bf.bhmin(n)
c      IF (n.eq.5) PRINT*, 'b-f H-', akap     
C-----BOUND-FREE TRANSITIONS OF H I
       IF ( bf.with.h1 )
     *    akap = akap + AHYDBF(divt,w,1.,emi)*atm.pion(n,1,1) /
     *           atm.uion(n,1,1)
c      IF (n.eq.5) PRINT*, 'b-f H', akap     
C-----BOUND-FREE AND FREE-FREE ABSORPTION BY H2+ - MOLECULES

      IF ( h2plus )
     *  akap = akap+AH2PLUS(divt,w)*atm.pion(n,1,1)*atm.pion(n,2,1)
c      IF (n.eq.5) PRINT*, 'b-f H2+ mol', akap     
C-----ABSORPTION BY H2 - QUASI-MOLECULES
      IF ( h2quasi ) 
     *   akap = akap + QUASIH2(t,w)*divt*
     *                     (atm.pion(n,1,1)/atm.uion(n,1,1))**2
c      IF (n.eq.5) PRINT*, 'b-f H2-', akap           
C-----BOUND-FREE-TRANSITIONS OF HE I
      IF ( hebf )
     *  akap = akap + AHEIBF(divt,w) * atm.pion(n,1,2)/atm.uion(n,1,2)
c      IF (n.eq.5) PRINT*, 'b-f HeI', akap           
C-----BOUND-FREE-TRANSITIONS OF HE II
      IF (bf.with.he2)
     *    AKAP = AKAP + AHYDBF(divt,w,2.,emi)*atm.pion(n,2,2) /
     *                                        atm.uion(n,2,2)
c      IF (n.eq.5) PRINT*, 'b-f HeII', akap           
C-----FREE-FREE TRANSITION FOR NEGATIVE IONS HMINUS, HE AND H2
C-----(STIMULATED EMISSION ALREADY INCLUDED)
C
      IF (bf.with.hmff) THEN

C-----Free-free  absorption  coefficient  for  HMINUS
C-----according to T.L John: Astron. Astrophys. 193,189 (1988)

        IF (w .LT. 1823.D-8) GOTO 90
        ICOEFF = 6
        IF (w .GT. 3645.D-8) ICOEFF = 12
        XMICLAM = w * 10000.         ! FORMULAE HOLD FOR WLAM IN MICRONS
        XMICLAM2 = XMICLAM * XMICLAM
        FACT = ( 5040. * divt ) ** .5
        AHMINFF = 0.
        DO I=ICOEFF,ICOEFF-5,-1
          AHMINFF =
     *  (   HFIT1(I) * XMICLAM2          + HFIT2(I)
     *    + HFIT3(I) / XMICLAM           + HFIT4(I) / XMICLAM2
     *    + HFIT5(I) / XMICLAM2 / XMICLAM + HFIT6(I)/ XMICLAM2/XMICLAM2)
     *    + AHMINFF * FACT
        END DO
        AHMINFF = AHMINFF * 1.E-29 * FACT * FACT
        AHMINFF = AHMINFF * 2. * atm.pion(n,1,1)/atm.uion(n,1,1)
        AMINFF  = AMINFF + AHMINFF
c       IF (n.eq.5) PRINT*, 'f-f H-, aminff', aminff
C --- FREE-FREE ABSORPTION COEFFICIENT FOR HE
90      CONTINUE
        AHE = 0.
        DO I = 1, 4
         A = 0.
         IF ( w .GT. 9113.D-8 ) THEN
            DO J = 1, 5
               A = HEFIT1(J,I) + A * w * 1.e-8
            ENDDO
            A = A * 1.E-10
         ELSE
            DO J = 1, 5
               A = HEFIT2(J,I) + A * w * 1.e-8
            ENDDO
            A = A * 1.E-10
         ENDIF
         AHE = A + AHE * divt
        ENDDO
        AMINFF = AMINFF + AHE * atm.pion(n,1,2)/atm.uion(n,1,2)
c      IF (n.eq.5) PRINT*, 'f-f He aminff', aminff

C MB 12.2010
C-----FREE-FREE ABSORPTION COEFFICIENT FOR MOLECULAR HYDROGEN H2 
C --- (3/96 : mit John 1994 zu aktualisieren !)
        AH2 = 0.
        DO I = 1, 4
         A = 0.
         IF ( w .GT. 11391.D-8 ) THEN
           DO J = 1, 5
              A = H2FIT1(J,I) + A * w * 1.e-8
           ENDDO
           A = A * 1.E-10
         ELSE
           DO J = 1, 5
              A = H2FIT2(J,I) + A * w * 1.e-8
           ENDDO
           A = A * 1.E-10
         ENDIF
         AH2 = A + AH2 * divt
        ENDDO
        AMINFF = AMINFF + AH2 * atm.pmol(n,mol_idx(1,1))
c             IF (n.eq.5) PRINT*, 'f-f H2 mol aminff', aminff           
      ENDIF

      IF (bf.with.metff) THEN
C-----FREE-FREE-TRANSITIONS, FOR ALL IONS IN THE HYDROGENIC APPROXIMATION
        x = 0.
        DO k = 1, 2
          zz = FLOAT(k * k)
          x = x + zz * atm.pp(n,k) * GAHFF(t,w,zz)
        ENDDO
        AHFF = 9.93E-32 * ( W * SQRT(divt) ) ** 3 * X
        AKAP = AKAP + AHFF * atm.pe(n) 
c      IF (n.eq.5) PRINT*, 'ff met, hydro', akap           
      ENDIF

C MB 12.2010
C-----CH Photodissociation 
      IF (chbf) THEN
       IF (atm.pmol(1,1).NE.0.0) THEN
        nm = mol_idx(1,6)           ! 6 - for Carbon
        akap = akap +  atm.pmol(n,nm)/atm.umol(n,nm)*CHOP(t,w)*
     *           atm.h.eps_dev(1)*atm.h.eps_dev(6)
c             IF (n.eq.5) PRINT*, 'CH mol akap', akap
       ENDIF           
      ENDIF  

C-----OH Photodissociation 
      IF (ohbf) THEN
       IF (atm.pmol(1,1).NE.0.0) THEN
        nm = mol_idx(1,8)           ! 8 - for Oxygen
        akap = akap +  atm.pmol(n,nm)/atm.umol(n,nm)*OHOP(t,w) 
     *      * atm.h.eps_dev(1)*atm.h.eps_dev(8)
c             IF (n.eq.5) PRINT*, 'OH mol akap', akap
       ENDIF          
      ENDIF  

C-----ELECTRON AND RAYLEIGH-SCATTERING FOR H I, HE I AND H2 
      sig = 0.
      IF (bf.with.thompson) sig = SIGE * atm.pe(n)
c             IF (n.eq.5) PRINT*, ' thomson, sig', sig           

      IF (bf.with.rayleigh) THEN
        sray = 0.
        DIVLAM2 = 1.e-16  / ( w * w )         ! w in cm !
        DIVLAM4 = DIVLAM2 * DIVLAM2

C-----RAYLEIGH-SCATTERING FOR HYDROGEN (HI) AND HELIUM (HEI)

        DO I = 1, 2
         IF ( w .GE. EDGE(I) ) THEN
            SIGMA = SIGCO(I,1) *DIVLAM4 *( 1. + SIGCO(I,2) * DIVLAM2
     *                                        + SIGCO(I,3) * DIVLAM4 )
            sray  = sray + sigma * atm.pion(n,i,1)/atm.uion(n,i,1)
         ENDIF
        ENDDO
c             IF (n.eq.5) PRINT*, 'ray, sray', sray         

C-----RAYLEIGH-SCATTERING FOR MOLECULAR HYDROGEN (H2)
        IF ( w .GE. EDGE(3) ) THEN
         SIGMA = SIGCO(3,1) * DIVLAM4 * ( 1. + SIGCO(3,2) * DIVLAM2
     *                                       + SIGCO(3,3) * DIVLAM4 ) 
         sray = sray + sigma * atm.pmol(n,mol_idx(1,1)) 
         sig  = sig + sray
        ENDIF                              ! ----- Rayleigh
c            IF (n.eq.5) PRINT*, 'sig = sig +ray', sig
      ENDIF
C --- ---------------------------------------------------------------

C --- BOUND-FREE-TRANSITIONS OF METALS
      IF ( metbf ) THEN 
        akmetbf = 0.
        DO i = 1,bf.metmax
         mm = bf.z(i).id
         IF ( w .LE. bf.z(i).w0 ) THEN
           ion = bf.z(i).ion
           alpha = 1. 
           IF (bf.z(i).maxval.GT.1) THEN
              alpha = bf.z(i).val(2) - w
              IF (bf.z(i).maxval.GT.2) 
     *           alpha = bf.z(i).val(3) + alpha*w
           ENDIF
           alpha = ((w/bf.z(i).w0) ** bf.z(i).s) * 
     *      ( bf.z(i).val(1) * alpha * atm.h.eps_dev(mm) *
     *        atm.pion(n,ion,mm)/atm.uion(n,ion,mm)        )
           IF (bf.z(i).elow .GT. 0.) THEN
              ALPHA = ALPHA * EXP(-(bf.z(i).elow * divt))
           ENDIF
           akmetbf = akmetbf + alpha
         ENDIF
        ENDDO
        akap = akap + akmetbf 
c                  IF (n.eq.5) PRINT*, 'b-f metals', akmetbf
      ENDIF
C --- ---------------------------------------------------------------

      xkap (n) = (em * akap + aminff * atm.pe(n)) * div_kt
      xsig (n) = sig * div_kt

      ENDDO     !-- end depth loop

C --- ---------------------------------------------------------------
C check whether model atmosphere and SIU internal kappa_background 
C calculations are consistent. Use the old model tau scale for 
C discrepancy up to 1% MB 13.05.11
C --- ---------------------------------------------------------------

      IF (read_bf) THEN
        gen_new_tau_scale=.FALSE.
        xdiff=0.
        DO n=1,atm.h.ndepth        
          xkaptmp(n) = xkap(n) + xsig(n)
          ytmp = xkaptmp(n)/atm.kappa_ref(n) - 1.
          IF (abs(ytmp).GT.0.01.AND..NOT.BJTEST(lf_ctrl,26)) THEN
c           WRITE(CHN,*) 'inconsistent kappa_ref scale: ',
c     *                  n, atm.t(n), ytmp*100.,'%'
c           WRITE(CHN,*) 'kappa_5000: ',
c     *                  n, atm.t(n),xkaptmp(n),
c     *                  atm.kappa_ref(n)
            xdiff=MAX(xdiff,ABS(ytmp))
            gen_new_tau_scale=.TRUE.           
          ENDIF
        ENDDO
       IF (gen_new_tau_scale) THEN
             WRITE(*,*) 'New kappa/tau scale will be used'
             WRITE(*,*) gen_new_tau_scale
          ELSE 
            WRITE(*,*) 'Model atm. kappa/tau scale will be used'
            WRITE(*,*) gen_new_tau_scale 
       ENDIF
        IF (gen_new_tau_scale) THEN
        ! assume not equidistant tauscale
           atm.h.nbreak    = 0
           atm.h.nscale(1) = 0
           fc_lam(1) = xkaptmp(1)/atm.kappa_ref(1)  * atm.tau(1)
           dxtau(1)  = 1.
           DO n=2,atm.h.ndepth        
             fc_lam(n) = xkaptmp(n)/atm.kappa_ref(n)  * atm.tau(n)
             dxtau(n)  = LN10_DIV_2 * (atm.logtau(n)-atm.logtau(n-1))
           ENDDO
c call tau integrator with trapez (set nbreak, nscale to 0)
           CALL NEW_TAU_SCALE(xtau,fc_lam,dxtau,atm.h.ndepth,nn
     *          ,1.e35,2,atm.h.nscale,atm.h.nbreak)
c          IF (.NOT.BJTEST(lf_ctrl,26))
             write(CHN,*) 'New tau(5000) scale generated: ',
     *         'kappa_5000(atmosph.) is inconsistent with ',
     *         'kappa_5000(lineform) by 1 -',xdiff*100.,'%'
           show_tau_one = .NOT.BJTEST(lf_ctrl,26)
           DO n=1,atm.h.ndepth    
c              IF (show_tau_one.AND.xtau(n).gt.1) THEN
c              write(CHN,*) n, 'old log(tau5000): ',atm.logtau(n)   ! model atm tau
c               write(CHN,*) n, 'new log(tau5000): ',alog10(xtau(n)) ! new tau
c       write(CHN,*) n, 'tau5000,N,O: ',alog10(xtau(n)),atm.logtau(n)
               show_tau_one = .FALSE.
c              ENDIF
c
c ---- NEW optical depth scale -------------------------------- 
             atm.tau(n)       = xtau(n)
             atm.logtau(n)    = alog10(xtau(n))
             atm.kappa_ref(n) = xkaptmp(n)             ! new reference kappa(lambda)
           ENDDO
        ENDIF
        read_bf = .FALSE.
        GOTO 500                                       ! calculate kappa(w)
      ENDIF

      IF (wmin .LE. 0) THEN               ! only kappa(5000)
        atm.eta_cont_w(1) = w
         atm.h.eta_cont_nr = 0
         DO n=1,atm.h.ndepth        
           div_kapref        = 1./atm.kappa_ref(n)
           atm.eta_cont(n,1) = xkap(n)*div_kapref
           atm.etasig(n,1)   = xsig(n)*div_kapref
         ENDDO
         GOTO 8889                        ! kappa finished
      ELSE
c        PRINT*,'wmin GT 0, use new atm.kappa_ref'    ! used
         atm.eta_cont_w(il)   = w
         DO n=1,atm.h.ndepth
           div_kapref         = 1./atm.kappa_ref(n)
           atm.eta_cont(n,il) = xkap(n)*div_kapref
           atm.etasig(n,il)   = xsig(n)*div_kapref
         ENDDO
         IF (bf.with.fac) THEN
          IF (w.GE.bf.wmin_fac.AND.w.LE.bf.wmax_fac) THEN
           DO n=1,atm.h.ndepth
            atm.eta_cont(n,il) = atm.eta_cont(n,il)*bf.fac(n)
           ENDDO
          ENDIF
         ENDIF
         atm.h.eta_cont_nr=il
         IF (il .GE.MAX_ETA_CONT) THEN 
           WRITE(CHN,*) 'Too many continuum points. Increase dwmax !'
           GOTO 8889
         ENDIF
         IF (w .GE. wmax) GOTO 8889          ! kappa finished
         w = w*(1.D0 + wlam)
         IF (iwcnt.LT.iwend) THEN
           IF(w.GT.bf.edge.w(iwcnt)) THEN
             w  = bf.edge.w(iwcnt)
             iwcnt  = iwcnt + 1
             IF (iwcnt.GT.bf.edge.maxwave) iwcnt = bf.edge.maxwave
           ENDIF
         ENDIF          
         w = DMIN1(w,wmax)
      ENDIF
      GOTO 501                            ! calculate next wavelength
! usually 3 points across a line profile

 8889 CONTINUE

      info = ' ' 
      IF ( bf.with.h1       ) info = 'H bf,'  // info 
      IF ( bf.with.hm       ) info = 'H- bf,' // info  
      IF ( bf.with.hp       ) info = 'H+ bf,' // info  
      IF ( bf.with.h2m      ) info = 'H2- bf,'// info  
      IF ( bf.with.h2       ) info = 'H2 bf,' // info  
      IF ( bf.with.h2p      ) info = 'H2+ bf,'// info  
      IF ( bf.with.he1      ) info = 'HeI bf,'// info  
      IF ( bf.with.he2      ) info = 'He2 bf,'// info  
      IF ( bf.with.hmff     ) info = 'H- ff,' // info  
      IF ( bf.with.h2mff    ) info = 'H2- ff,'// info  
      IF ( bf.with.heff     ) info = 'He ff,' // info  
      IF ( bf.with.chbf     ) info = 'CH bf,'// info 
      IF ( bf.with.ohbf     ) info = 'OH bf,'// info 
      IF ( bf.with.metff    ) info = 'Metal ff,' // info 
      IF ( bf.with.metal    ) info = 'Metal bf,' // info 
      IF ( bf.with.link     ) info = 'Links,'    // info
      IF ( bf.with.rayleigh ) info = 'Rayleigh,' // info  
      IF ( bf.with.thompson ) info = 'Thompson,' // info 
      bf_info=info
      info = 'selected Background opacities: '//info

      IF (.NOT.BJTEST(lf_ctrl,26)) THEN
       WRITE(CHN,2405) atm.h.eta_cont_nr, nr_edges, RBF_FILE
 2405  FORMAT(1X,I4,' frequency points (background opacities), ',
     *       1X,I4,' points from ',A20) 
       WRITE(CHN,*) info
       IF ( bf.bhmin(1).NE.1.) WRITE(CHN,*) 
     *             'departures b(H-) applied !' 
       IF ( bf.with.fac) WRITE(CHN,*) 
     *             'scaling factor ',bf.fac(1),' applied between ',
     *  bf.wmin_fac*1.e8,' and ',bf.wmax_fac*1.e8,' A' 
        WRITE(CHN,*) ' '
      ENDIF
      RETURN
      END 
