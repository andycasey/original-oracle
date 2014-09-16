C --- Directives:
C --- SWCTRL = 1     stepwidth control checking kappa_line/kappa_cont
C --- SWCTRL = 2     stepwidth control checking curvature of profile
#define SWCTRL  2       
C --------------------------------------------------------------------------------
C <<<
C
C Lineform Documentation
C ======================
C NAME        : lineform
C FILENAME    : lineform.for
C DIRECTORY   : /usr/users/reetz/model/lineform/code
C
C PURPOSE     : This program calculates the emergent spectral energy distribution 
C               for a given planparallel, homogeneous LTE-atmosphere 
C               (temperature- and [partial-]pressure-structure and continuum
C                background opacities)
C               and a given list of atomic and molecular linedata.
C
C PROCEEDING  : 1) Data Input: atmosphere file
C                              linetable, 
C                              pretabulated Balmer-profiles 
C                              optional: departures,
C                                        depth-dependent microturbulence
C               2) Calculation of kappa-coefficients, dopplerbroadening and 
C                  dampingfactors for each line at each depth-point 
C               2b)NLTE correction of opacities and true absorption 
C                  Sourcefunction using departures (see data input)  
C               2c)Optional iteration of Sourcefunction Slam starting with
C                  true absorption Sourcefunction Strue
C               3) Integration of Slam(tau_lam)*E2(tau_lam)      (FLUX)
C                            or   Slam(tau_lam)*exp(-tau_lam/mu) (INTENSITY)
C                  (ref. B.Baschek, H.Holweger, G.Traving in 
C                   "Abhandlungen aus der Hamburger 
C                  Sternwarte", 1966) using weight factors 
C                  w_f and tau_f of 
C                        J.Reetz             (1992,unpublished) 
C                  w_i and tau_i of
C                        Laguerre coefficients for n=6
C                        (no difference to R.Cayrel, G.Traving 
C                            (Z. Astrophysik 50, 239, 1960) found)
C
C MODUL       : FORTRAN - SUBROUTINE
C
C INCLUDEFILES: physcnst.inc      physical constants
C               lf_decl0.inc      general declarations
C               lf_decl1.inc      'NLTE' declarations
C               lf_decl2.inc      'LINETABLE' declarations
C               lf_decl3.inc      coefficients for Lambda-Iteration (ATLAS9)
C               elemdata.inc      atomic data
C               mol_data.inc      molecule-id's and data
C
C OBJECTFILES : ~/reetz/obj/maths.o    (math utilities with routines of T.Gehren)
C               ~/reetz/obj/igetenv.o  (to query environment table)
C MAIN PGM    : ~/reetz/model/lineform/code/call_lineform.c    (interface to IDL)
C
C DOCUMENT    : see LINEFORM.DOC             (planned, not available)
C
C PHYS.UNITS  : ALL PHYSICAL QUANTITIES ARE REPRESENTED IN cgs-UNITS IF
C               NOTHING ELSE HAS BEEN REMARKED             ===
C
C REMARKS     : In opposite to usual FORTRAN Programming Practice, each 
C               variable has to be declared explicit before its use.
C               Declaration-Integrity should be checked at compilation-time
C               using the the qualifier /WARNINGS=(DECLARATIONS)
C
C DATAFILES   : 1) Depth-dependend microturbulence
C                  Optional entry of depth-dependent microturbulence
C                  if BIT 16 set.
C                  Microturbulence in [km/s]
C                  Filename  : LINFOR_MICROTURB  (environment variable)
C                  Fileformat: ASCII
C                  --------------------------------------------------------------
C                  1. row     <dim>                   1 or 2 dimensions
C                  2. row     <Xi_v_1>    <Xi_h_1>    outermost microturbulence
C                                                     vertical and horizontal
C                                                     (if dim=2)
C                  3. row     <Xi_v_2>    <Xi_h_2>    ..........
C                  .....
C                  NDEPTH+1.row   <Xi_v_n>    <Xi_h_n>    innermost microturb.
C                  --------------------------------------------------------------
C
C               2) Depth-dependend departures
C                  Optional entry of depth-dependent departures
C                  if only (!) BIT 20 set.
C                  Filename  : POP_FILE  (environment variable)
C                  Fileformat: ASCII, keyword processing
C                  --------------------------------------------------------------
C                  ATOM=<Element character>
C                  NDEPTH=<Nr of Atmosphere Depthpoints>
C                  LEVEL=<DESIGN.> <Order Nr> <Energie[eV] (dummy)> <Ionisationlevel>
C                  RBB=<Index of Lower Level> <Index of upper Level> <Wavelength[A]> <dummy>
C                  DEPARTURES=<Level Index> 
C                  .... NDEPTH departures and log_tau_5000
C                  --------------------------------------------------------------
C                  Example:
C
C                  %
C                  ATOM=O
C                  NDEPTH=80
C                  LEVEL=13P5PC      1        -1.00000   1
C                  LEVEL=13P5PE      2        -1.00000   1
C                  LEVEL=13P5PG      3        -1.00000   1
C                  LEVEL= 13S5S      4        -1.00000   1
C                  %RBB=<Lower Lev.> <Upper Lev.> <Wavelength> <gf>
C                  RBB=   4      3           7771.96          0.000
C                  RBB=   4      2           7774.18          0.000
C                  RBB=   4      1           7775.40          0.000
C                  DEPARTURES=1
C                       0.203270     -6.00000
C                       0.152480     -5.78571
C                         ........
C                       1.00000      1.43077
C                       1.00000      1.50000
C                  DEPARTURES=2
C                       0.116420     -6.00000
C                       0.118720     -5.78571
C                         ........
C                  DEPARTURES=3
C                         ........
C                  DEPARTURES=4
C                         ........
C
C               3) Depth-dependend departures
C                  Optional entry of depth-dependent departures
C                  if only BIT 20  and BIT 24 set.
C                  departure file : LINTAB_LOOKUP   
C                       Fileformat: BINARY  (as created by line_report.pro; mode=11)
C                  departure file : POP_FILE        Fileformat: BINARY
C                       Fileformat: BINARY  (as created by alinvel.al)
C
C CONVENTIONS : Following programming conventions are used:
C
C                    1) FORTRAN commands are written     uppercase
C                    2) Constants are written            uppercase
C                    3) Names of subroutines are written uppercase
C                    4) Variables are written            lowercase
C
C               Control-Bitmap LF_STATUS:
C               =========================
C               Bitnr                      meaning
C                             not set                    set
C               ------------------------------------------------------------
C                0            NO FATAL ERROR           FATAL ERROR
C                1            NO ERROR                 Teff Gridpoint missing
C                2            NO ERROR           Logg Gridpoint missing
C                3            NO ERROR           LogZ Gridpoint missing
C                4            NO ERROR           VCS: TEFF interp. out of bounds
C                5            NO ERROR           VCS: NE   interp. out of bounds
C
C
C               Control-Bitmap LF_CTRL:
C               =======================
C               Bitnr                      meaning
C                             not set                    set
C               ------------------------------------------------------------
C                0      stupid                     intelligent
C                1      irradiance                 flux
C                2      use all lines (see 13)     use selected (sun)lines
C                3                                 search exact atmosphere
C                4      parameter input            file input
C                5      absolute                   residual irradiance/flux 
C                6                                 file output
C                7      spectrum output            contribution function output
C                                                  (Bitfield 10,11:  > 0)
C                8,9    tau(lambda)-scale integration
C                       0: Simpson, 1: linear trapezinteg., 2: log. trapezinteg.
C                    
C               10,11   Integral(S(tau)*E(tau)dtau)
C                       0: 6 point Gaussquadrature
C                       1: Splineintegration (linear)
C                       2: Splineintegration (log.)
C
C               12                                 fit modus see <fit>.ini file
C               13      use all lines (see 2)      use only molecular lines
C               14      include Balmerlines        no Balmerlines
C               15      include C4-interaction     no C4-interaction
C               16                                 read microturbulence file
C               17      read linetable
C               18                                 calculate S/B 
C               19      dynamic step width control use wavelengths of wlam_out
C               20                                 NLTE-Lineformation
C               21                                 only lines within interval
C               22                                 use explicitly assigned atmosp.
C               23      (WMIN,WMAX) given          only line exactly at WMIN is
C                                                  used (if NOT marked as 
C                                                        "DELETED" in linetable)
C                                                  WMAX = (elem2)(elem1).(ion): e.g. 01.01
C                                                  Lineforming region is +/-5A
C               24                                 check NLTE lookup table
C               25                                 exact_line: -->
C                                                  wmin=rec_nr, wmax=rec_pos
C               26                                 no messages
C               27 not 'equidistant'               'equidist.' (res.dist.=dwmin at 5000A)
C               28 ignore deleted lines            include deleted lines
C               29        -                        ignore predicted lines
C               30 line symmetry                   line asymmetry (see below)
C               31 true absorption scatter (sig*S) Exact scattering term (sig*J)
C 
C AUTHOR      : Johannes K. Reetz
C DATE        : September 1992 (init)
C               
C HISTORY     : This Procedure was derived from the LIDA-Subroutine of ABUND.FOR
C               (T.Gehren, 1988) and from LINFOR.FOR (T.Gehren, C.Reile,
C                K.Fuhrmann, J.Reetz) which was derived from the Model-
C                Atmosphere-Programm MAFAGS.FOR
C               Lambda-Iteration (ITERATE_SOURCE_FCT) due to Kurucz' ATLAS9
C               
C               Init: 15.09.1992 (JKR)
C                     01.10.1992 (JKR) Diverse Tests:
C                                     a) zur Bestimmung ausreichender
C                                        Dichte des Tiefenpunkte-Rasters
C                                        mit dem Ergebnis 80 Tiefenpunkte
C                                        zu nehmen, bei feinerer Unterteilung
C                                        im Bereich tau_ref > 0.001
C                                     b) Die Gewichte und Stuetzstellen von
C                                        Cayrel zur Berechnung des Flusses 
C                                        mittels Gaussquadratur erwiesen sich
C                                        als mangelhaft (Vergleich mit den
C                                        Balmerlinien, die vom bisherigen
C                                        lineformation Pgm. (Feautrier Methode) 
C                                        geliefert wurden)
C
C                     22.06.1993 (JKR) Eingabe einer Mikroturbulenzschichtung
C                                     moeglich. Die Schichtung hat NDMAX
C                                     Tiefenpunkte, entsprechend der 
C                                     Atmosphaere.
C                     22.07.1993 (JKR) Linien werden in read_lines
C                                     nach aufsteigenden Wellenlaengen
C                                     sortiert.
C                     23.07.1993 (JKR) Verwertung von Departures ermoeglicht.
C                                     Departure-Datei befindet sich in
C                                     'linfor_populations'
C                     01.12.1993 (JKR) Interpolation auch bei kolinear angeordneten
C                                     Atmosphaeren Gitterpunkten moeglich
C                     19.12.1993 (JKR) NLTE-Lineformation ueberprueft, Meldungen
C                                     ueberarbeitet, Spezifizierung der
C                                     Position von Uebergaengen auf der Wellen-
C                                     laengenscala moeglich
C                     30.12.1993 (JKR) lf_status Ausgaben: Fataler Fehler, wenn
C                                     BitNr.0 gesetzt. Dieses Bit sollte immer
C                                     vom Clienten (dem aufrufenden Programm)
C                                     abgefragt werden.
C                     15.03.1994 (JKR) Atmosphaerenfile Format geaendert:
C                                     Nun koennen Files mit 
C                                     unterschiedlicher Anzahl von Tiefen-
C                                     punkten, Elementen, Stuetzstellen von
C                                     kontinuierlichen Opazitaeten (eta_cont)
C                                     eingelesen werden.
C                                     DENNOCH: Gitterpunkte muessen alle vom
C                                              SELBEN Typ sein !!!!!!
C                                     Die Flexibilisierung bezieht sich nur
C                                     auf die Verarbeitung einzelner "exakter"
C                                     Atmosphaeren.
C                                     Wichtige Angaben, wie z.B. die Anzahl
C                                     der Tiefenpunkte werden in der Struktur
C                                     ATM.H abgespeichert.
C                    
C                     22.03.1994      Atmosphaeren filename kann als
C                                     VMS LOGICAL  'atmosphere' explizit
C                                     spezifiziert werden.
C                     09.04.1994      Lineformation einer einzelnen Linie
C                                     moeglich. Bit 23 steuert die Bedeutung
C                                     von WMIN und WMAX
C                     14.04.1994      siehe Bit 25
C                     09.08.1994      siehe Bit 28 und 29
C                     12.08.1994      Lineprofile asymetries according to
C                                     Szudy & Baylis 1975, JQSRT 15,641
C                                     see BIT 30
C                                     (does not work currently)
C                     22.11.1994      Contributionfunction of the
C                                     emergent Linedepression
C                                       Intensity: Magain, 1986
C                                       Flux     : Achmad, de Jager 1991
C                     01.12.1994      Cayrels Gewichte zur Quadratur der
C                                     Intensitaet ersetzt durch Laguerre
C                                     gewichte. Tests haben allerdings
C                                     keine signifikanten Unterschiede
C                                     gezeigt.
C                     02.12.1994      Voigtprofile ueberprueft:
C                                     Fuer einige Linien (OI7774,Mgb5172
C                                     und OH3130) wurde die approximierte
C                                     Voigtfunktion durch die exakte
C                                     Funktion ersetzt. Fazit: keine
C                                     signifikanten Unterschiede in dem
C                                     diesbezueglichen (v,a) Gebiet.
C                     02.01.1995      Exact inclusion of the scattering term
C                                     sigma*J according 
C                                     to 'linfor' (Holweger et al.) and
C                                     ATLAS9 (Kurucz 1992)
C                                     see bit 31
C                                     Therefore: Format of Atmosphere had to 
C                                                be changed !
C                                                see array: etasig
C                     03.01.1995      Wenn Bit 24 und Bit 31 gesetzt sind:
C                                     Ausgabe des monochromatischen tiefen-
C                                     abhaengigen Strahlungsfeldes (u.ae.)
C                                     in lineform_rad.txt (current directory)
C                     07.01.1995      If linefile has extension '.linfor'
C                                     the 'qline.dat' ascii fileformat
C                                     is assumed (see ~/reetz/linedata/o.linfor)
C                     08.01.1995      Compiler Directive SCAT
C                                      SCAT= 0 (ATM input: eta_cont 
C                                             includes scattering coeff. etasig)
C                                      SCAT= 1 (ATM input includes etasig
C                                         
C                     24.06.1995      Wellenlaengen gitter kann fest vorgegeben werden
C                                     siehe Bit 19
C                     03.07.1995      Fit Modus (BIT 12): Die fuer den Fit relevanten
C                                     relevanten Definitionen werden in einem
C                                     ini-File abgelegt (siehe Environment variable
C                                     LF_FIT_DEF)
C                                     Gelesen wird das initializations file
C                                     durch READ_FIT_DEF(filename,fit_def)
C                     25.07.1995      Neue Behandlung der NLTE Quellfunktion:
C                                     Linienquellfunktionen, die zustarke
C                                     negative Beitrage bringen werden
C                                     ignoriert !
C                      1.10.1995      Mean Values of S/B, tau_5000,T and
C                                     dT/dtau are calculated and stored into
C                                     wlam_out if Bit 18 of lf_ctrl is set
C                     21.03.1996      Aenderung vom 25.07.1995 rueckgaengig
C                                     gemacht
C                                     Aenderung vom 08.01.1995 rueckgaengig
C                                     gemacht
C                           --->      Background-continuum opacities are
C                                     explicitly calculated (Routine KAPPA).
C                                     You have to assign a data file: uphoto.asc
c                   -------------------------------
c                   \OPAC H1,H-bf,H2-,H2Q,H2+,He1,H-ff,H2-ff,Heff,RAYLEIGH,THOMPSON,METFF
c                   \OPAC CHBF,OHBF,METAL
c                   \METAL
c                       2071.3  AL I 2P0   13    1  0.000   6.   4.40  6.500E-17   0.000   0.000
c                       2513.8  MG I 3P0   12    1  2.714   9.   2.70  2.000E-17   0.000   0.000
c                   \END METAL
c                   \WAVE 
c                    \Z_EDGE 0.05
c                      3647.04  0.1  H  I (Balmer)
c                      16419.00  0.5  H- bf
c                   \END WAVE
c                   -------------------------------
C                                     New format of atmospheres (see ATMOSIN1)
C                     05.12.1996      Interpolation of departures from *.pop files
C                     22.03.1997      waals(2,lc) * atm.pion(n,1,2) war
C                                     seit 21.3.96 falsch: 
C                                         waals(2,lc) * atm.pion(n,2,1) 
C                     03.02.1998      automatische Abschaetzung
C                                     maximaler Linienbreiten 
C                                   (siehe XKAPMIN1, XKAPMIN, ignore_cnt, ...)
C                     05.02.1998     Balmerlinien H11....H18 jetzt nach EWS
C                                    Balmerlinien H3 ....H10 wie bisher nach VCS                                       
C                     10.01.1999     NLTE lookup table (see BIT 24) to read NLTE-departures
C                                    filename of lookup table: $LINTAB_LOOKUP
C                                    IF BIT24 set $LINTAB_LOOKUP is inspected
C                                    to identify transitions which are actually
C                                    used for NLTE-Lineformation 
C                                    NOTE: BIT24 is only active when BIT20 is set
C
C                     02.12.2005              fug Frank GRUPP takes over this special version!
C                     02.12.2005     MOD01.01 BALMER-NLTE first coding (take over from Lyudmila)
C                     05.12.2005     MOD02.01 BALMER-NLTE parameters (AG-BPO, scaling, FineStructure on-off)
C                     12.12.2005     MOD03.01 NEW BALMER TABELS (Preconvolved)
C                     13.12.2005     MOD04.01 BPO
c                     16.12.2005     MOD05.01 Hydrogen-NLTE
c                     02.04.2008     MOD06.01 Lyudmila Mashonkina comments and corrections regarding BPO
C                     03.02.2011     MOD07.01 MB some modifications to make automatic EW calculations
C                                             MB note that this lineform.f can not be used with the 
C                                                IDL version of SIU, because it contains extra 
c                                                parameter: ewidth (eq. width)
C >>>
C     -------------------------------------------------------------------------
C

      SUBROUTINE lineform(wmin_in,wmax_in,teff_in,
     *                    glog_in,zlog_in,
     *                    ximicro_in,cos_theta_in,eps_dev_in,
     *                    lf_ctrl_in,sw_crit_in,dwmin_in,dwmax_in,          
     *                    linedata_file_in,lf_status_out,
     *                    wlam_cnt_out,wlam_out,flux_out,ewidth,
     *                    head_i_out,head_f_out,
     *                    head_atmdate,head_atmuser,head_atmref,
     *                    head_atmcmt,head_lindate,head_linuser,
     *                    head_linref,head_lincmt,EXCPOT)

C     ------------------------ Declarations -----------------------------------
      INCLUDE 'physcnst.inc'                    ! physical constants
      INCLUDE 'lf_decl0.inc'                    ! general declarations
      INCLUDE 'lf_decl1.inc'                    ! NLTE-Lineform,.....
      INCLUDE 'elemdata.inc'                    ! standard element abundances
      INCLUDE 'mol_data.inc'                    ! molecule-id's and data
      INCLUDE 'iondat.inc'                      ! partition functions and ionisation energies
      INCLUDE 'lf_info.inc'



C --- stuff for hydrogen NLTE                         ! MOD01.01 fug
      common/chieta/mlhf                              ! MOD01.01 fug
      dimension hyddep(21,ndmax)                      ! MOD01.01 fug

C --- Input- and Outputvariables ----------------
      INTEGER*4 lf_ctrl_in                      ! program control bitmask
      REAL*4    wmin_in,wmax_in                 ! min. and max. wavelength [Angstroem]
      REAL*4    teff_in                         ! Effective-temperature [K]
      REAL*4    glog_in                         ! gravity acceleration
      REAL*4    zlog_in                         ! metallicity
      REAL*4    ximicro_in                      ! microturbulence       [km/s]
      REAL*4    cos_theta_in                    ! cos(theta) - only for irradiance
      REAL*4    eps_dev_in(NATOM)               ! abundance deviations
      REAL*4    sw_crit_in                      ! stepwidth criterium      
      REAL*4    dwmin_in                        ! min. stepwidth [mA] 
      REAL*4    dwmax_in                        ! max. stepwidth [Angstroem]
      CHARACTER*(*) linedata_file_in            ! path of linetable

      INTEGER*4 head_i_out(NHEAD_I)             ! header  (integer)
      INTEGER*8 wlam_cnt_out                    ! Nr of flux points
      INTEGER*4 lf_status_out                   ! program status  bitmask
      REAL*4    head_f_out(NHEAD_F)             ! header  (real)
      REAL*8    flux_out(NWMAX)                 ! flux       points
      REAL*8    wlam_out(NWMAX)                 ! wavelength points     [Angstroem]

      CHARACTER*20  head_atmuser                ! header: atm user
      CHARACTER*20  head_linuser                ! header: lin user
      CHARACTER*30  head_atmdate                ! header: atm generation date
      CHARACTER*30  head_lindate                ! header: lin generation date
      CHARACTER*100 head_atmcmt                 ! header: atm comment
      CHARACTER*100 head_lincmt                 ! header: lin comment
      CHARACTER*255 head_atmref                 ! header: atm reference
      CHARACTER*255 head_linref                 ! header: lin reference
      CHARACTER*255 store_lup                   ! stores name of lup file MOD05.01 fug
      CHARACTER*255 store_dep                   ! stores name of dep file MOD05.01 fug
      REAL*4        dep_in(21,ndmax)            ! MOD05.01 fug
      REAL*8        nlte_hyd_fac                ! MOD05.01 fug 
      REAL*8        nlte_hyd_fac2               ! MOD05.01 fug 
      CHARACTER*10  atmbuilder                  ! MOD07.11 mb

    
C --- local variables (in alphabetic order) -----
      LOGICAL*4 no_atm_read
      LOGICAL*4 atm_exist,file_exist
      LOGICAL*4 interpol_departures             !
      LOGICAL*4 first_depart                    !
      LOGICAL*4 log_set                         ! = (CHN .NE. STD_OUTPUT_CHN) 
      LOGICAL*4 sort_lines                      ! flag for read_lines
      INTEGER*4 firstline,actline,lastline      ! Idx of first,current,last line
      INTEGER*4 iteration_step,max_iteration_step !
      LOGICAL*4 id_used(NATOM+NMOLMAX)          ! = TRUE, if element/molecule is concerned
      INTEGER*4 idx_w                           ! Index of current wavelengthpoint
      INTEGER*4 nwords                          ! for STRTRIM (number of non-Blank words)
      INTEGER*4 k_eta_cont                      ! current index of eta_cont  
      INTEGER*4 first_k_eta_cont                ! first index of eta_cont
      INTEGER*4 flux_integr_method              ! 0,1,2; for definition see LF_ORG.DOC
      INTEGER*4 tau_integr_method               ! 0,1,2; for definition see LF_ORG.DOC
      INTEGER*4 iximicro_dim                    ! 0:microturbulence is constant
                                                ! 1,2:1- or 2-dim. microturb.
      REAL*4    ximicro2                        ! = ximicro^2
      REAL*4    ximicro2_arr(NDMAX)             ! = xi(tau)^2, if depthdependent
      REAL*4    ximicro_arr(NDMAX,2)            ! xi(tau,1)=vertical turbulence 
                                                ! xi(tau,2)=horizontal turb.  
      REAL*4    qstark(NLINE)                   ! quadratic Stark dampingfactor
      REAL*4    waals (2,NLINE)                 ! van der Waals dampingfactor
      REAL*4    yint(NDMAX)                     ! Buffer     
      REAL*4    as(NDMAX),bs(NDMAX),cs(NDMAX)   ! Buffer for Spline Coefficients
      REAL*4    cont_flux_1,cont_flux_2         ! left and right continuumpoints of continuum
      REAL*4    div_cos_theta                   ! 1./cos(theta)  
      REAL*4    div_dlamdop(NDMAX,NLINE)        ! = 1./ (vdop(n,k)*line(lc).w)
      REAL*4    divt(NDMAX)                     ! = 1/T(n)
      REAL*4    divtn_FAC_B                     ! = FAC_B / T(n)
      REAL*4    xta                             ! Buffer
      REAL*4    dxtau(NDMAX)                    ! Buffer for delta log(tau) OR delta tau
      REAL*4    dtdtau(NDMAX)                   ! Buffer for dT/dlogtau
      REAL*4    f_lam (NDMAX)                   ! = eta_lam(n) * tau(n)
      REAL*4    eta_cf (NDMAX)                  ! = for contribution function
      REAL*4    S_B(NDMAX)                      ! = S/B  (if requested, see Bit 18)
      REAL*4    eta_line_vec(NDMAX)             ! = kappa_lam(line)/kappa_ref
      REAL*4    Sline(NDMAX)                    ! = line source function
      REAL*4    fc_lam(NDMAX)                   ! = (eta_cont_lam+etasig_lam) * tau(n)
      REAL*4    beta(NDMAX)                     ! = etasig / eta_lam
      REAL*4    strue(NDMAX)                    ! = S-true absorption
      REAL*4    gf(NLINE)                       ! gf    = 10**log(gf)
      REAL*4    elow_k(NLINE)                   ! lower energy in [K]
      REAL*4    hcc_2_div_lam5                  ! 2*h*c^2 / w^5  (for calc. Planckfunction)
      REAL*4    xx                              ! Buffer variable
      REAL*4    tim                             ! elapsed time (seconds)
      REAL*4    old_tim                         ! time buffer (seconds)
      REAL*4    tn_BETA                         ! = tn ** -0.7  (temporary)
      REAL*4    tn_R_2                          ! = 2 * T(n) * R
      REAL*4    eps_dev10                       ! peculiar abund. dev. (temp.)
      REAL*8    dx                              ! temporary buffer
      REAL*8    div_fac                         ! = 1./(wlam_cont_2 - wlam_cont_1)
      REAL*8    em                              ! 1 - exp(-E/kT) = factor for stim. emission
      REAL*8    emi                             ! exp(-E/kT)
      REAL*8    gamma_tmp                       ! temporary
      REAL*8    fsp(4)                          ! fluxpoint-buffer for splinecalc.
      REAL*8    rtmp                            ! temporary buffer
      REAL*8    xstart,yend  ,xend              !     "       "
      REAL*8    w_int,w1,w2                     ! for calculating eta_cont at w
      REAL*8    wfac                            ! for calculating eta_cont at w
      REAL*8    wlam_cont_1,wlam_cont_2         ! left and right wavelengthpoints of continuum
      REAL*8    wline_2                         ! = line(lc).w * line(lc).w
      REAL*8    wsp(4)                          ! wavelength-buffer for splinecalc.
      REAL*8    x,y                             ! temporary buffer
      REAL*8    dw_eps                          ! = MAX(dwmin*0.8,W_EPS)
      REAL*8    flux_tmp                        ! temporary flux
      REAL*8    a,v
      REAL*4    div_dlg(NDMAX,NLINE)            ! for asymetric line profiles
      REAL*4    gam2waals(NDMAX,NLINE)          ! for asymetric line profiles
      REAL*8    tmp_fac
      REAL*8    dlamx
      REAL*8    xdwmin                          ! actually used dwmin
      REAL*8    xdwmax                          ! actually used dwmax
      INTEGER   agbpo                           ! AG or BPO MOD02.01 fug
      INTEGER   fsyesno                         ! Hydrogen finestructure yes no MOD02.01 fug
      REAL*8    resscale                        ! Resonance profile scaling factor MOD02.01 fug
      INTEGER   agwrong                         ! Convolution or adding            MOD02.01.fug

      LOGICAL   read_bf           ! if TRUE reread RBF data from file (KAPPA)

      CHARACTER*10 buf_date                     ! buffer for date
      REAL*8    xxxx                            ! temporary buffer     MOD06.01 lm
      CHARACTER*20 lin_account                  ! user account (header item)
      CHARACTER*80 rec,substr
c      CHARACTER*80 substr
c      CHARACTER*56 rec
      CHARACTER*200 atm_name1,atm_name_old
 
      CHARACTER*1  str001
      CHARACTER*2  str002
      CHARACTER*3  str003
      CHARACTER*4  str004
      REAL*4    ewidth                          ! equivalent width - if negative - emission line

C --- Functions ---------------------------------
      LOGICAL*4 CHECK_STEPWIDTH_CRIT           ! stepwidth control (TRUE = one step back) 
      REAL*4    SECNDS                          ! 
      REAL*4    EXPINTEGRAL                     ! Integralexponentialfunction
      REAL*4    SPLINE_INTEGRAL                 ! integral of splinefunction 
      REAL*4    SPOL                            ! Splineinterpolation
      REAL*8    VOIGT                           ! 

      real*4    dum_ali           
 
      PARAMETER (max_iteration_step = 4)        ! for stepwidth control 

      COMMON /MAIN_WORK/ wlam_cont_1,
     *          wlam_cont_2,cont_flux_1,cont_flux_2,div_fac,
     *          k_eta_cont,wline_2,rtmp,x,w_int,w1,w2,wfac,
     *          gamma_tmp,emi,em,hcc_2_div_lam5,tn_BETA,ximicro2_arr,
     *          ximicro_arr,ximicro2,qstark,waals,elow_k,
     *          tn_R_2,divtn_FAC_B,f_lam,fc_lam,div_dlamdop,divt,gf,
     *          id_used,yint,dxtau,as,bs,cs,atm_name_old

C      ff(q)  = exp(-1.33*q)*log(1.+2.27/q) +
C     *              0.487*q / (0.153+q**(1.6667)) + q/(7.93+q**3)
      WRITE(*,*) 'LINEFORM.F'
      WRITE(*,*) wmin_in,wmax_in,teff_in,
     *                    glog_in,zlog_in,
     *                    ximicro_in

C     ------------------------ Data - Input -----------------------------------

      log_set = (CHN .NE. STD_OUTPUT_CHN)                              
      env_len = IGETENV('SIU_MAIN',env)
      IF (log_set)  OPEN(UNIT=CHN,FILE=env(1:env_len)//'/log',
     *                  STATUS='NEW',FORM='FORMATTED')

       tim =SECNDS(0.0)                                 ! Initialization of timebuffer
      old_tim = tim
      lf_ctrl = lf_ctrl_in

c --- hydrogen NLT INITS   MOD01.01 fug begin
      mlhf=0
      
      do i=1,20
         do ilm=1,ndmax
            hyddep(i,ilm)=1.0
         enddo
      enddo

c --- MOD01.01 fug end
c --- hydrogen NLTE parameter file read   MOD02.01 fug begin

      env_len = IGETENV('SIU_MAIN',env)
      OPEN(UNIT=99,FILE=env(1:env_len)//'/code/brd_control.in',
     *       FORM='FORMATTED',STATUS='OLD',READONLY)
      READ(99,'(I1)') agbpo
      READ(99,'(I1)') fsyesno
      READ(99,'(F5.3)') resscale
      READ(99,'(I1)') agwrong
      WRITE(*,*) '---------------------------'
      WRITE(*,*) agbpo,fsyesno,resscale
      WRITE(*,*) '---------------------------'
      CLOSE(99)

c --- MOD02.01 fug end

      IF (.NOT.BJTEST(lf_ctrl,26).AND.BJTEST(lf_ctrl,31)) 
     *        WRITE(CHN,*) 'with correct scattering term'

      IF (BJTEST(lf_ctrl,23)) THEN                   ! single lineformation
c MB 02.2011 - do not consider any other line, if LF of selected lines is requested
c       lf_ctrl = JIBCLR(lf_ctrl,2)                   ! consider all lines
       lf_ctrl = JIBSET(lf_ctrl,14)                  ! no balmerlines
       wline_center = DBLE(wmin_in)                  ! line wavelength [A]
       wline_code   = DBLE(wmax_in)                  ! line code
       wmin_in      = wline_center  - XMAX_LINE_RANGE*wline_center
       wmax_in      = wline_center  + XMAX_LINE_RANGE*wline_center
       IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *   'Forming only ',wline_code,' at ',wline_center 
      ENDIF
 
      IF (BJTEST(lf_ctrl,25)) THEN                   ! single lineformation
       lf_ctrl = JIBCLR(lf_ctrl,2)                   ! consider all lines
       lf_ctrl = JIBSET(lf_ctrl,14)                  ! no balmerlines
       lin_rec_nr = INT(wmin_in)+1
       lin_rec_pos= INT(wmax_in)+1
       wmin_in = 0.
       wmax_in = 0.
      ENDIF

C --- Define wavelength limits IF
      IF (BJTEST(lf_ctrl,19)) THEN   ! wlambda positions predefined
            IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                'using predefined lambda scale'
       wlam_cnt = 1                  ! see beginning of section LINEFORMATION
       i        = 1
       k        = wlam_cnt_out
        wmin_in = wlam_out(1)
        wmax_in = wlam_out(wlam_cnt_out)
      ENDIF
      wmin          = DBLE(wmin_in)    * 1.D-8
      wmax          = DBLE(wmax_in)    * 1.D-8

C ---
      IF (BJTEST(lf_ctrl,16)) THEN   ! microturbulence stratification ?
        head.lin.xi = -1.
      ELSE 
        head.lin.xi   = ximicro_in
        ximicro       = ximicro_in * 1.e5   ! [km/s] --> [cm/s]
      ENDIF

      cos_theta     = cos_theta_in
      linedata_file = linedata_file_in        
      CALL STRCOMPRESS(linedata_file,k)

      DO i = 1,NATOM
        atm.h.eps_dev(i) = 10.**eps_dev_in(i)
      END DO
        
C     RESET LF_STATUS
      lf_status = 0

C --- Check if atmosphere already exist --------------------------------------- 
C --- Get new atmosphere, if not       ---------------------------------------

      IF (BJTEST(lf_ctrl,22)) lf_ctrl = JIBCLR(lf_ctrl,3) ! ATM explicitely decl                                                                       
      no_atm_read = .FALSE.

      IF (BJTEST(lf_ctrl,22)) THEN         ! use explicitely assigned atmosphere                  
         env_len = IGETENV('ATMOSPHERE',env)    ! now hard-coded in lineform_cl.pl
         atm_name1 = env(1:env_len)
         INQUIRE(FILE = atm_name1,EXIST = atm_exist)
         IF (.NOT. atm_exist) THEN
           lf_status = JIBSET(lf_status,0)  ! Fatal Error
           IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *         'Atmosphere not found: ',env(1:env_len)
           GOTO 8999
         ENDIF
         IF (atm_name1.EQ.atm_name_old) THEN
           no_atm_read = .TRUE. 
         ELSE
           atm_name_old = atm_name1
         ENDIF   
      ELSE
         atm_name_old = ' '
      ENDIF
      IF ( (atm.h.teff.EQ.teff_in.AND.
     *      atm.h.logg.EQ.glog_in.AND.
     *      atm.h.z   .EQ.zlog_in.AND.
     *      BJTEST(lf_ctrl,3) )
     *      .OR. no_atm_read) THEN
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)  
     *   'Same atmosphere used: '//atm.h.ref(1:INDEX(atm.h.ref,'//'))
      ELSE
        atm.h.teff    = teff_in    
        atm.h.logg    = glog_in
        atm.h.z       = zlog_in

C --------------------------------------------------
         CALL ATMOSIN1                  !  atmosphere input

         atmbuilder = TRIM(head.atm.user)
         WRITE(*,*) atm.h.user
C         IF (atmbuilder.EQ.'Karin Lind') THEN
C           WRITE(CHN,*) 'Partition functions and partial pressures '
C     *             //'are computed by SIU with model atm abundances'
C           CALL PARTPRESS              !  partition function, part pressures
C         ENDIF
C --------------------------------------------------

        IF (BJTEST(lf_status,0)) THEN  !  Fatal Error ?
            atm.h.teff = 0    
            atm.h.logg = 0
            atm.h.z    = -9999
            GOTO 8999
        ENDIF
      ENDIF

      DO i=1,NATOM
        nelpos(i) = 0
      ENDDO

      DO i=1,atm.h.nel
        nelpos(atm.h.idx_el(i)) = i
      ENDDO

      pop.atom     = 0
      pop.ndepth   = 0
      pop.level_nr = 0
      pop.rbb_nr   = 0

C---------------------------------------------------------------------
C   NLTE lineformation without lookup table                  --- begin
C---------------------------------------------------------------------

      IF (BJTEST(lf_ctrl,20).AND..NOT.BJTEST(lf_ctrl,24)) THEN 
         WRITE(*,*) "DOING NLTE LINE FORMATION"
         env_len = IGETENV('POP_FILE',env)
         WRITE(*,*)  env(1:env_len)
         OPEN(UNIT=NLTE_UT,FILE=env(1:env_len),FORM='FORMATTED',
     *        STATUS='OLD',READONLY)
         
         interpol_departures  = .FALSE.
         first_depart = .TRUE.
        
 30      CONTINUE               ! loop
         READ(NLTE_UT,'(A80)',END=60) rec
         IF (rec(1:1) .NE. '%' .AND. rec(1:1) .NE.' ') THEN
            IF (rec(1:5) .EQ. 'ATOM=') THEN ! ATOM
               pop.atom=0
               substr = rec(6:LEN(rec))
               CALL STRCOMPRESS(substr,k)
               str002=substr(1:k)
               IF (str002(2:2).EQ.' ') THEN
                  str001 = str002(1:1)
                  str002 = ' '//str001
               ENDIF
               i = 1
               DO WHILE ((cel(i).NE.str002) .AND. i .LT. NATOM)
                  i = i + 1
               ENDDO
               IF (i.LE.NATOM) pop.atom = i
            ELSEIF (rec(1:7) .EQ. 'NDEPTH=') THEN ! NDEPTH
               READ(rec(8:LEN(rec)),*) pop.ndepth
            ELSEIF (rec(1:6) .EQ. 'LEVEL=') THEN ! LEVEL
               READ(rec(7:LEN(rec)),'(A6)')    desig
               READ(rec(7+NCHRLAB:LEN(rec)),*) k,x,i
               pop.level_nr = pop.level_nr + 1
               pop.level(k).desig = desig
               IF (desig(1:1).EQ.' ') pop.level(k).desig = 
     *              desig(2:NCHRLAB)//' '
               pop.level(k).ion   = i
               pop.level(k).ev    = x
            ELSEIF (rec(1:4) .EQ. 'RBB=') THEN
               pop.rbb_nr = pop.rbb_nr + 1
               i = pop.rbb_nr
               READ(rec(5:LEN(rec)),*) pop.rbb(i).ll,
     *              pop.rbb(i).ul,
     *              pop.rbb(i).w,dummy
               pop.rbb(i).w = pop.rbb(i).w * 1.D-8 ! [A] --> [cm] 
               pop.rbb(i).accept = .FALSE. ! init
            ELSEIF (rec(1:11) .EQ. 'DEPARTURES=') THEN
               READ(rec(12:LEN(rec)),*) i
               IF (first_depart) THEN 
                  READ(NLTE_UT,'(A80)') rec
                  CALL STRTRIM(rec,k,nwords)
                  IF (nwords.EQ.1) THEN ! test, whether tau-scales is present
                     READ(rec,*) pop.level(i).b(1)
                     DO k=2,pop.ndepth
                        READ(NLTE_UT,*,END=60) pop.level(i).b(k)
                     ENDDO
                     IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                'Departures: log(tau) scale not present ==> ',
     *               'atmosphere-tau-scale is presumed to be adequate'
                  ELSE
                     READ(rec,*) pop.level(i).b(1),dxtau(1)
                     DO k=2,pop.ndepth
                       READ(NLTE_UT,*,END=60) pop.level(i).b(k),dxtau(k)
                     ENDDO
                     interpol_departures = (atm.h.ndepth.NE.pop.ndepth)
                     IF (.NOT.interpol_departures) THEN
                        DO k=1,pop.ndepth ! test, whether tau-scales are different
                         IF (ABS(atm.logtau(k)-dxtau(k)).GT. 1.e-4) THEN
                              interpol_departures = .TRUE.
ccc    Mashonkina L.  August 9, 2005                         !MOD01.01 fug
ccc                         interpol_departures = .false.    !MOD01.01 fug
                              IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                             'DEPARTURES are spline-interpolated.'
                              GOTO 35
                           ENDIF
                        ENDDO
 35                     CONTINUE
                     ENDIF
                  ENDIF
                  first_depart = .FALSE.
               ELSE             ! other departures
                  IF (interpol_departures) THEN
                     DO k=1,pop.ndepth
                       READ(NLTE_UT,*,END=60) pop.level(i).b(k),dxtau(k)
                     ENDDO
                  ELSE
                     DO k=1,pop.ndepth
                        READ(NLTE_UT,*,END=60) pop.level(i).b(k)
                     ENDDO
                  ENDIF
               ENDIF
               
C     --- Interpolation of departures, if TAU-Scales are inconsistent
               IF (interpol_departures)  THEN
                  DO k=1,atm.h.ndepth
                   CALL SPLINE(dxtau,pop.level(i).b,as,bs,cs,pop.ndepth)
                  yint(k)=SPOL(atm.logtau(k),dxtau,pop.level(i).b,as,bs,
     *                    cs,pop.ndepth)
                  ENDDO
                  DO k=1,atm.h.ndepth
                  pop.level(i).b(k) = yint(k)
                  ENDDO
               ENDIF
            ENDIF
          ENDIF

      GOTO 30
 60   CONTINUE
      IF (interpol_departures) pop.ndepth = atm.h.ndepth
      CLOSE(NLTE_UT)        
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,4010) 
 4010 FORMAT(' NLTE-lineformation  -  selected populations')
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,4011) (cel(pop.atom),
     *     pop.level(i).ion,
     *                   pop.level(i).desig,
     *                   pop.level(i).ev,   i=1,pop.level_nr)
 4011   FORMAT((' Atom: ',A2,I1,'  Level: ',A6,
     *          '  Energy: ',F10.6,' eV'))

C --- Calculate possible transitions between given levels
C   
        IF (pop.rbb_nr .LE. 0) THEN  ! transitions not explicitly declared
          DO i = 1,pop.level_nr
            idx_arr(i) = i
            tmp_arr(i) = pop.level(i).ev
          ENDDO 

          CALL SORT(pop.level_nr,tmp_arr,idx_arr)
             
          DO i = 1,pop.level_nr-1
             DO k=i+1,pop.level_nr
                IF ((pop.level(idx_arr(k)).ev - 
     *               pop.level(idx_arr(i)).ev) .GT. MIN_EV) THEN                
                  n = pop.rbb_nr + 1
                  pop.rbb_nr = n
                  pop.rbb(n).ll=idx_arr(i)
                  pop.rbb(n).ul=idx_arr(k)
                ENDIF
             ENDDO                                                  
          ENDDO
        ENDIF

        DO n = 1,pop.rbb_nr
           IF (pop.rbb(n).w .LE. 0) THEN
             i = pop.rbb(n).ll
             k = pop.rbb(n).ul
             x = pop.level(i).ev
             y = pop.level(k).ev
             IF ( (i.NE.k) .AND. (x.LT.y)) pop.rbb(n).w = 
     *                                          HC_E_SI / (y - x) 
           ENDIF
           pop.rbb(n).accept = .FALSE.          ! init
        ENDDO
  
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,4012) 
     *                  (pop.level(pop.rbb(n).ll).desig,
     *                   pop.level(pop.rbb(n).ul).desig,
     *                   (pop.rbb(n).w*1.D8),   n=1,pop.rbb_nr)
 4012   FORMAT((' Transition: ',A6,' --> ',A6,'   Wavelength: ',
     *            G13.4,' A'))


      ENDIF   ! ---- end of NLTE LOOP ----------------------------

      sort_lines = .FALSE.
      EXCPOT = -1.
      IF (.NOT.BJTEST(lf_ctrl,17)) CALL read_lines(sort_lines, EXCPOT) !  linetable  input

C---------------------------------------------------------------------
C   NLTE lineformation using lookup table                  --- begin
C---------------------------------------------------------------------

      IF (BJTEST(lf_ctrl,20).AND.BJTEST(lf_ctrl,24)) THEN
c ---  open lookup tables
        env_len = IGETENV('LINTAB_LOOKUP',env)
        IF (env_len.LE.0) THEN
               WRITE(CHN,*) 'ERROR READING LINTAB_LOOKUP' 
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
        ENDIF
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)
     *     ' Inspecting lookup table: ',env(1:env_len)
        INQUIRE(FILE = env(1:env_len),EXIST = file_exist)
        IF (.NOT.file_exist) THEN
               WRITE(CHN,*) 'lookup table is missing (fatal error)' 
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
        ENDIF

        OPEN(UNIT = 33, FILE = env(1:env_len),STATUS = 'OLD', 
     *       FORM = 'UNFORMATTED',RECORDTYPE='STREAM', READONLY )
        store_lup = env(1:env_len)                                 !MOD05.01 fug
c --- open index file for departures
        env_len = IGETENV('POP_FILE',env)
        IF (env_len.LE.0) THEN
               WRITE(CHN,*) 'ERROR READING POP_FILE' 
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
        ENDIF
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*)
     *     ' inspecting departure database ',env(1:env_len)
        INQUIRE(FILE = env(1:env_len)//ext_depidx,EXIST = file_exist)
        IF (.NOT.file_exist) THEN
          WRITE(CHN,*) 'Index file of departure database missing: '
          WRITE(CHN,*) ' ',env(1:env_len)//ext_depidx
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
        ENDIF
        INQUIRE(FILE = env(1:env_len)//ext_dep,EXIST = file_exist)
        IF (.NOT.file_exist) THEN
               WRITE(CHN,*) 'departure main file missing: ' 
               WRITE(CHN,*) ' ',env(1:env_len)//ext_dep
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
        ENDIF
        WRITE(*,*) "DEP file: ",env(1:env_len)//ext_dep
c
c --- read idx file
c        
        OPEN(UNIT=34,FILE=env(1:env_len)//ext_depidx,STATUS='OLD', 
     *       FORM = 'UNFORMATTED',RECORDTYPE='STREAM', READONLY )
            READ(34) ndepth_ali,nlev_ali,natom_ali,niter_ali,
     *        teff_ali,logg_ali,feh_ali,xi_ali,
     *        (epsdev_ali(i),i=1,natom_ali),atm_ali,
     *        hbbfac_ali,hbffac_ali,escale_ali,rbfscl_ali,
     *        cbf_scl_ali,gmax_ali,gmin_ali,
     *        (term_ali(i),i=1,nlev_ali),
     *        (xdum(i),i=1,ndepth_ali),
     *        (xdum(i),i=1,ndepth_ali),
     *        (logtau5000_ali(i),i=1,ndepth_ali)
        CLOSE(34)
        pop.ndepth = ndepth_ali

        print*,'ALI/ATM depths: ', pop.ndepth, atm.h.ndepth
        interpol_departures = (atm.h.ndepth.NE.pop.ndepth)
        IF (.NOT.interpol_departures) THEN
          DO k=1,pop.ndepth          ! test, whether tau-scales are different
            IF (ABS(atm.logtau(k)-logtau5000_ali(k)).GT. 1.e-3) THEN
              interpol_departures = .TRUE.
              IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                 'DEPARTURES are spline-interpolated!'
              GOTO 9120
            ENDIF
          ENDDO
9120      CONTINUE
        ENDIF
c
c --- read lookup table
c
        READ(33) mult_code,wmin_lup,wmax_lup,n_nlte_lines,n_terms   

        WRITE(*,*) 'MOIN 1: ',mult_code
        dw_eps = w_eps/5000.D-8
        READ(33) (terms(i),i=1,n_terms)

c      in case of hydrogen nlte manipulate lookup table data
        IF (store_lup.EQ.                                          !MOD05.01 fug begin
     *      '/home/homam/fug/siu/linedata/master_line.dat.lup.01.01')
     *  THEN
           mult_code = 01.01
           wmin_lup = 4000.
           wmax_lup = 8000.
           n_nlte_lines = 20
           n_terms = 20
        ENDIF                                                      !MOD05.01 fug end
c
c --- read dep file (all supplementary info loaded from idx file)
c        
        store_dep=env(1:env_len)//ext_dep                          !MOD05.01 fug
        OPEN(UNIT=34,FILE=env(1:env_len)//ext_dep,
     *      ACCESS='DIRECT',RECL=NDMAX,STATUS='OLD',READONLY)
        w = 0.

        DO lc = 1,linemax         
         IF (line(lc).id.EQ.INT(mult_code)) THEN 
                    mult_ion=INT(100*(mult_code-line(lc).id))
         IF (mult_ion.LE.1.OR.line(lc).ion.EQ.mult_ion)THEN
           IF (w.GE. line(lc).w*(1.D0-dw_eps)) GOTO 4002 
 4001      CONTINUE
           READ(33,end=9115) w_lup,el_lup,  ! label 9115 --> no more lines need to be checked
     *                       xjl_lup,xju_lup,il_lup,iu_lup
           w = w_lup*1.D-8
           IF (w.lt. line(lc).w*(1.D0-dw_eps)) GOTO 4001
 4002      CONTINUE
           IF (w.gt. line(lc).w*(1.D0+dw_eps)) GOTO 9110
           IF (line(lc).xjl.GE.0.) THEN
              IF (.NOT.(xjl_lup.EQ.line(lc).xjl.AND.
     *             xju_lup.EQ.line(lc).xju)) GOTO 4001
           ENDIF
           IF (ABS(el_lup-line(lc).elow/(E_SI * 1.E7)).GT.EL_EPS) 
     *          GOTO 4001
         
C ----     transition found in lookup_table ---> processing
           IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(*,3045),
     *         terms(il_lup),terms(iu_lup),cel(line(lc).id),
     *         line(lc).ion,w_lup
 3045      FORMAT(' NLTE TRANSITION: ',A6,' - ',A6,
     *       ' (',A2,I1,' ',F9.3,' A)')
           pop.rbb_nr = pop.rbb_nr + 1
           pop.rbb(pop.rbb_nr).accept = .TRUE.
           line(lc).rbb_idx=pop.rbb_nr
           pop.rbb(pop.rbb_nr).w = line(lc).w
           DO idum = 1, 2
             IF (idum.EQ.1) THEN
               ikk=il_lup
               pop.rbb(pop.rbb_nr).ll = pop.level_nr+1
             ELSE
               ikk=iu_lup
               pop.rbb(pop.rbb_nr).ul = pop.level_nr+1
             ENDIF
C ---  check whether departures have been already read
             DO ip=1,pop.level_nr
               IF (idx_list(ip).eq.ikk) THEN
                 IF (idum.EQ.1) THEN
                   pop.rbb(pop.rbb_nr).ll = ip
                 ELSE
                   pop.rbb(pop.rbb_nr).ul = ip
                 ENDIF
                 GOTO 9060  ! try to read next departures
               ENDIF
             ENDDO
						 
             DO ip = 1,nlev_ali
               IF (term_ali(ip).EQ.terms(ikk)) THEN
                itmp=ikk
                ikk =ip
                GOTO 9050
               ENDIF
             ENDDO
             IF ( mult_code.EQ.01.01) GOTO 9050              !MOD05.01 fug
**             DO ip = 1,nlev_ali
**               write(*,*) term_ali(ip),'.',terms(ikk)
**             ENDDO
             pop.rbb(pop.rbb_nr).accept = .FALSE.
             WRITE(CHN,*) 'Designation ',terms(ikk),
     *          ' from $LINTAB_LOOKUP not found'
             WRITE(CHN,*) 'in departure database $POP_FILE.'
             lf_status = JIBSET(lf_status,0)         ! Fatal Error
             RETURN
             GOTO 9060
 9050        CONTINUE
             pop.level_nr = pop.level_nr + 1
             idx_list(pop.level_nr) = itmp
             pop.level(pop.level_nr).ev = line(lc).elow
             WRITE(*,*) 'level Nr and seq num: ', pop.level_nr, ikk
C
C
C
c             DO kmm = 1, nlev_ali              
c             READ(34,rec=kmm) (xdum(i),i=1,pop.ndepth)
c               print*, kmm, (xdum(i),  i=1,pop.ndepth)
c             ENDDO
C
C
C
              
             IF (interpol_departures)  THEN
C ---  Interpolation of departures, if TAU-Scales are inconsistent
                READ(34,rec=ikk) (xdum(i),i=1,pop.ndepth)
                DO k=1,atm.h.ndepth
                   CALL SPLINE(logtau5000_ali,xdum,as,bs,cs,pop.ndepth)
                   pop.level(pop.level_nr).b(k)=
     *                  SPOL(atm.logtau(k),logtau5000_ali,
     *                  xdum,as,bs,cs,pop.ndepth)
                ENDDO
             ELSE
                READ(34,rec=ikk) 
     *               (pop.level(pop.level_nr).b(i),i=1,pop.ndepth)
C                  print*, (pop.level(pop.level_nr).b(i),i=1,pop.ndepth)
             ENDIF
                
 9060        CONTINUE
           ENDDO
 9110      CONTINUE  ! transition not found in lookup table
         ENDIF
         ENDIF             ! IF (line(lc).id.EQ.INT(mult_code)) THEN 
        ENDDO              ! DO lc = 1,linemax
      ENDIF
 9115 CONTINUE
      CLOSE(33)
      CLOSE(34)


c ***  NLTE HI ---------------------------- MOD01.01 fug begin
      IF ( mult_code.EQ.01.01) pop.atom=1   
      WRITE(*,*) "POP.ATOM=",pop.atom   
      
      IF(pop.atom.eq.1) THEN
         OPEN(89,FILE=store_dep,STATUS='UNKNOWN',               
     &        RECORDTYPE='STREAM',FORM='UNFORMATTED')                               
         READ(89) ((dep_in(i,k),k=1,80),i=1,21) 
         WRITE(*,*) "READING IN HYDROGEN NLTE DEPARTURES"
         hyddep = 1.*dep_in
         CLOSE(89)    
      ENDIF
c END NLTE HI ----------------------------- MOD01.01 fug begin

C -----------------------------------------------------------------------------
C  NLTE lineformation using lookup tables --- end
C -----------------------------------------------------------------------------

      IF (.NOT.BJTEST(lf_ctrl,14)) THEN
       hyd = .NOT. (( wmin.GT. 1.500e-5).AND.(wmax.LT. 3.646e-5) .OR.
     *              ( wmin.GT. 5.100e-5).AND.(wmax.LT. 6.350e-5) .OR.
     *              ( wmin.GT. 6.800e-5).AND.(wmax.GT. 6.800e-5))

C ------ insert Balmer line data ----------------------------------------------

       IF (hyd) THEN                       
        DO i = 1,NBALMER_ALL                      
          IF (wmin.LE.wbalmer(i)+(2.0D-6).AND.
     *        wmax+(2.D-6).GE.wbalmer(i)) THEN    ! balmer line into array 'line.w'
            DO lc = linemax,1,-1
             IF (wbalmer(i).GT.line(lc).w) THEN 
               idx_w = lc + 1
               GOTO 250 
             ENDIF
            END DO   
            idx_w = 1

  250       CONTINUE
            IF (idx_w.LE.linemax) THEN
              DO lc = linemax,idx_w,-1
                line(lc+1) = line(lc)
              ENDDO
            ENDIF
            line(idx_w).w     = wbalmer(i)
            line(idx_w).id    = 1           ! H
            line(idx_w).ion   = 1
            line(idx_w).elow  = 0.              
            line(idx_w).gflog = balmer_gf(i)               
            line(idx_w).c4log = 0.        
            line(idx_w).c6log = 0.         
            line(idx_w).grad  = 0.         
            line(idx_w).width = 1.e-6          ! 100 Angstroem !
            line(idx_w).mode  = 0
            line(idx_w).use   = .true.
            line(idx_w).rbb_idx = 0
            linemax           = linemax + 1
          ENDIF
        ENDDO
       ENDIF

       IF (hyd) THEN                                ! MOD03.01 fug begin
          IF ((agwrong.eq.1).OR.(agbpo.EQ.1)) THEN
             CALL READ_VCS_HYD  !  STARK-Function data for Balmer lines
                                !  according to Vidal, Cooper and Smith
          ELSE
             WRITE(*,*) ' CALLING READ_HYDPROF'
             CALL READ_HYDPROF  ! Preconvolved tabels
             WRITE(*,*) ' Back FROM READ_HYDPROF'
          ENDIF
       ENDIF                                         ! MOD03.01 fug end
      ELSE
        hyd = hyd.AND.(.NOT.BJTEST(lf_ctrl,14))

      ENDIF

C ----------------- read microturbulence data (optional) ---------------------

      IF (BJTEST(lf_ctrl,16)) THEN

        env_len = IGETENV('LINFOR_MICROTURB',env)
        IF (env_len.LE.0) THEN 
               WRITE(CHN,*) 'ERROR READING LINFOR_MICROTURB' 
               lf_status = JIBSET(lf_status,0)  ! Fatal Error
               GOTO 8999
      ENDIF                        
      OPEN(UNIT = 30, FILE = env(1:env_len),                
     *       STATUS = 'UNKNOWN', FORM = 'FORMATTED', READONLY )
          READ(30,*) iximicro_dim
          IF (iximicro_dim .EQ. 1) THEN
            DO i = 1,atm.h.ndepth
              READ(30,*) ximicro_arr(i,1)       ! vertical direction
              ximicro_arr(i,1) = ximicro_arr(i,1) * 1.e5
              print*, ximicro_arr(i,1)
            ENDDO
          ELSE
            DO i = 1,atm.h.ndepth
              READ(30,*) ximicro_arr(i,1), ximicro_arr(i,2)
              ximicro_arr(i,1) = ximicro_arr(i,1) * 1.e5
              ximicro_arr(i,2) = ximicro_arr(i,2) * 1.e5
            ENDDO
          ENDIF
        CLOSE(30)
      ELSE
        iximicro_dim = 0
      ENDIF

       IF (.NOT.BJTEST(lf_ctrl,1)) THEN         ! irradiance calculation
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                  'IRRADIANCE CALCULATION AT MU= ',cos_theta
       ENDIF

C --------------------------- Preparations -----------------------------------

      sw_crit = DBLE(sw_crit_in)                     ! stepwidth criterium
      dwmin   = DBLE(dwmin_in)   * 1.D-11
      IF (dwmax_in.GT.0.) THEN
        dwmax   = DBLE(dwmax_in) * 1.D-8
      ELSE
        dwmax   = 1.D3*dwmin
      ENDIF
      IF (sw_crit.LE.0.D0) THEN sw_crit  = STD_SW_CRIT
      IF (dwmin.LE.0.D0) dwmin    = STD_DWMIN * (wmin+wmax) * .5
      IF (dwmax.LE.0.D0) dwmax    = STD_DWMAX * (wmin+wmax) * .5
      IF (BJTEST(lf_ctrl,27)) dwmax = dwmin
      IF (.NOT.BJTEST(lf_ctrl,26))
     *       WRITE(CHN,*) 'Max. point spacing',dwmax*1.D8,'/5000.' 

      dw_eps  = DMAX1(dwmin*0.8D0,W_EPS) ! real minimal stepwidth

C---------------------------------------------------------------------------
C --- Calculate Background opacities and new logtau scale

      read_bf = .TRUE.
      x = dwmax
      IF (BJTEST(lf_ctrl,27)) x = dwmin*1000.

      CALL KAPPA(read_bf,x*DIV_WLAMREF,wmin,wmax)     ! precalculation of bf opacities
                                                      ! and new log tau 5000 scale
C ---------------------------------------------------------------------------
C --- Define methods to calculate new opt.depth scales of continuum
C --- and lines, and method of integration over all depth points for intensity or flux

      tau_integr_method  = 0
      CALL MVBITS(lf_ctrl,8,2,tau_integr_method,0)   ! integration method for tau_lam 
      flux_integr_method = 0
      CALL MVBITS(lf_ctrl,10,2,flux_integr_method,0) ! integration method 
                                                     ! for E2(tau)*B(tau)

      IF (flux_integr_method.EQ.0) THEN              ! Gauss integration method
       IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 'Gauss quadrature' 
       IF (BJTEST(lf_ctrl,1)) THEN                   ! Flux calculation
        DO i = 1,6
          ta(i) = tf(i)
          wa(i) = wf(i)
        END DO
       ELSE                                          ! Irradiance calculation
        DO i = 1,6
          ta(i) = ti(i) * cos_theta
          wa(i) = wi(i)                             
        END DO
       ENDIF 
      ELSE                                           ! Spline integration method
       IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 'Spline-Integration' 
       ta(6)= 1.e35                                  ! see subroutine NEW_TAU_SCALE
       IF (.NOT.BJTEST(lf_ctrl,1)) div_cos_theta = 1./cos_theta   
      ENDIF                     
      
      IF (atm.h.nbreak.LE.1 .AND. atm.h.nscale(1).EQ.0) THEN
        WRITE(CHN,*) 'tau5000 is not equidistant.'
        tau_integr_method = 2
      ENDIF

c ---------------------------------------------------------------------------
C --- Preparation for new tau scale evaluations
c DEFAULT for MAFAGS models is Simpson (but set to Trapez for MARCS because
c it is unclear how to handle break points 
c
      IF (tau_integr_method.EQ.0) THEN               ! Simpson (log. scale)
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 'Simpson int. log tau'
        dxtau(1) = 1.                                ! see new_tau_scale
        k = 1
        DO j=1,atm.h.nbreak
          m = atm.h.nscale(j)                        ! breakpoint of tau scale
          dxtau(k+1) = (atm.logtau(k+1)-atm.logtau(k)) * LN10_DIV_2
          DO n = k+2,m
            dxtau(n) = LN10_DIV_3 * (atm.logtau(n)-atm.logtau(n-1)) ! (2/6-method)
            IF ( ABS(atm.logtau(n)+atm.logtau(n-2)  
     *               - 2.*atm.logtau(n-1)) .GT. 1.E-4) THEN 
                IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                'WARNING: Wrong Break found at ',
     *                 atm.logtau(n),' Do not use Simpson Integration !'
            ENDIF
          END DO
          k = m
        END DO               
      ELSEIF (tau_integr_method.EQ.1) THEN           ! Trapez (linear scale)
        dxtau(1) = atm.tau(1)                        ! see new_tau_scale
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *     'monochr. tau scale: lin. trapez integration'
        DO n = 2,atm.h.ndepth
          dxtau(n) = (atm.tau(n) - atm.tau(n-1)) * 0.5
        END DO
      ELSE                                           ! Trapez (log. scale)
        dxtau(1) = 1.                                ! see new_tau_scale
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *     'monochr. tau scale: log. trapez integration'
        DO n = 2,atm.h.ndepth
          dxtau(n) = LN10_DIV_2 * (atm.logtau(n)-atm.logtau(n-1))
        END DO
      ENDIF

C --- Linedata preparations -----------------------------------------------------

      IF (BJTEST(lf_ctrl,30).AND..NOT.BJTEST(lf_ctrl,26)) ! not used (see header)
     *      WRITE(CHN,*) 'ASYMMETRIC LINEFORMATION'

      IF (.NOT.BJTEST(lf_ctrl,17)) THEN              ! only if linetable has been read 
       DO i = 1,NATOM+atm.h.nmol
        id_used(i) = .FALSE.
       END DO

       DO lc = 1,linemax
        line(lc).use = .true.
        IF (line(lc).id.NE.1) THEN                   ! no calculations for H 
         id_used(line(lc).id) = .TRUE.               ! this element/molecule is used
         wline_2 = line(lc).w * line(lc).w

*** MB         print*, 'wav: ', line(lc).w, ' log gf: ', line(lc).gflog

         gf(lc)  = (10. ** line(lc).gflog) * wline_2

         elow_k(lc) = DIV_K  * line(lc).elow         ! energy in [K]

C ---    calculation of radiation-damping with classical formulae if not already done

c MB
         IF (line(lc).grad .LE. 0.0) line(lc).grad = .222342 / wline_2

C ---    van der Waals damping (preparation)
C ---    formula see Unsoeld (1955), p.332
       
         rtmp = 10. ** (0.4 * line(lc).c6log + 1.9583E+01)            ! dv=C6/r^6
         x    = 1./mass(line(lc).id)
         waals(1,lc) =          rtmp * (1.+ x * mass(1)) ** 0.3       ! H
         waals(2,lc) = 0.4133 * rtmp * (1.+ x * mass(2)) ** 0.3       ! He

         IF (BJTEST(lf_ctrl,30)) THEN               ! asymmetric lineformation
c           tmp_fac = 6.8685541D-08 * wline_2 * (1+x)**0.6
c           div_dlg(NDMAX,lc)=(10.**(0.2*line(lc).c6log))/tmp_fac
         ENDIF

C ---    quadratic Stark damping (preparation)
C ---    formula see Unsoeld (1955), p.327
         IF (.NOT.BJTEST(lf_ctrl,15)) THEN    ! C4-Interaction
          qstark(lc) = 10.**(0.666667 * line(lc).c4log + 1.9362E+01) ! dv=C4/r^4
         ENDIF
        ENDIF                                 ! IF (line(lc).id.NE.1) .....
       ENDDO                                                    
      ENDIF                                   ! if line table has been read


C --- microturbulence preparations --------------------------------------------
      IF (iximicro_dim.GT.0) THEN
       IF (iximicro_dim.EQ.2) THEN
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                  '2-dim microturbulence stratification'
         IF (.NOT.BJTEST(lf_ctrl,1)) THEN   ! irradiance
           IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                   'Calculation of Xi(mu,tau)'
           DO n = 1,atm.h.ndepth
             x = ximicro_arr(n,1) * ximicro_arr(n,1)    ! vertical
             y = ximicro_arr(n,2) * ximicro_arr(n,2)    ! horizontal
             ximicro2_arr(n) = (y + cos_theta*cos_theta*(x - y))
c             ximicro2_arr(n) = x*y/(y + cos_theta*cos_theta*(x - y))
           ENDDO
         ELSE            
           DO n = 1,atm.h.ndepth                    
             ximicro2_arr(n) = 0.5*(ximicro_arr(n,1)+ximicro_arr(n,2))  ! flux  
             ximicro2_arr(n) = ximicro2_arr(n) * ximicro2_arr(n)
           ENDDO
         ENDIF
       ELSE                ! v_turb = h_turb  ---> isotropic turbulence
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                'Isotropic microturbulence stratification'
         DO n = 1,atm.h.ndepth                    
           ximicro2_arr(n) = ximicro_arr(n,1) * ximicro_arr(n,1)
         ENDDO
       ENDIF
       DO n=1,atm.h.ndepth
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,1005) 
     *            n,atm.logtau(n),SQRT(ximicro2_arr(n))*1.e-5
1005     FORMAT(' TAU(',I3,')= ',F7.4,'    XI =  ',F8.3,' (km/s)') 
       ENDDO
      ELSE
        ximicro2 = ximicro*ximicro
        DO n = 1,atm.h.ndepth
          ximicro2_arr(n) = ximicro2
        ENDDO
      ENDIF
c
C ------------------- LINE ABSORPTION COEFFICIENT --------------------------

c      PRINT*, 'LTE number densities: ID, kap_l, N_tot (FeI), N_l'
c      PRINT*, '                      -- 6200 A (b3F2 - y3F3*)'
      DO n = 1, atm.h.ndepth                          
        tn      = atm.t(n)
        divtn   = 1. / tn
        divt(n) = divtn
        divtn_FAC_B = divtn * FAC_B
        tn_R_2  = tn * R_2
        tn_BETA = (tn ** -0.7)
        IF (hyd) THEN
            fh(n) = FAC_C * atm.pion(n,1,1)/atm.uion(n,1,1)*divtn
            f0(n) = FAC_D * (atm.pe(n) * divtn )**0.66667
        ENDIF

C ---   Doppler broadening of atomic lines for each element
        DO m = 1, atm.h.nel
          atomic_nr = atm.h.idx_el(m)
          IF (id_used(atomic_nr)) vdop(n,atomic_nr) = 
     *               SQRT(ximicro2_arr(n)+tn_R_2/mass(atomic_nr))/C
        END DO                     

C ---   Doppler broadening of molecular lines for each molecule
        DO m = 1, atm.h.nmol
          k = NATOM + m
          IF (id_used(k)) vdop(n,k) = 
     *                     SQRT(ximicro2_arr(n)+tn_R_2/mass(k))/C
        END DO                                 
      
C ---   Determination of the total damping-factor gamma for metallic lines,
        ! including radiative and van der-Waals damping (see Unsoeld 1966)

        DO lc = 1,linemax
         IF (id_used(line(lc).id)) THEN
          k = line(lc).id
          gamma_tmp  = 0.

C --- C6 --------------------------------------------
          gamma_tmp  = (waals(1,lc) * atm.pion(n,1,1) + 
     *                  waals(2,lc) * atm.pion(n,1,2)) * tn_BETA
          
          IF (BJTEST(lf_ctrl,30))   THEN    ! asymmetric lineformation
           WRITE(CHN,*) 'someone attempts asymmetric line formation?'
           WRITE(CHN,*) 'it doesnot work...'
c            tmp_fac = gamma_tmp*0.5*line(lc).w*line(lc).w/C
c            div_dlg(n,lc)   = div_dlg(NDMAX,lc)*(divtn**0.6)*tmp_fac
c            gam2waals(n,lc) = tmp_fac**2
          ENDIF

C --- C4 --------------------------------------------

          IF (.NOT.BJTEST(lf_ctrl,15)) THEN                  ! C4-Interaction
           E_up    = line(lc).elow + H*C/line(lc).w
           XIminusE= chi(line(lc).ion,line(lc).id)*E_SI*1.e7 - E_up
           xneff2  = line(lc).ion**2*2.178657e-11/XIminusE    ! (n*)^2
           IF (line(lc).c4log.LE.0.) THEN
C --- Unsoeld, p.327
             gamma_tmp1 = qstark(lc)*atm.pe(n)*divtn**0.83333
           ELSE
C           a0_bohr = 5.2917706e-09
           rrr     = 1.58114*xneff2/line(lc).ion              ! radius in a0_bohr
           xtmp    = 26299.702*divtn*line(lc).ion**2          ! eta*R

C --- Cowley 1971, Observatory
           gamma_tmp1  = 1.2e10*atm.pe(n)*divtn**1.5*xneff2**2

1410        FORMAT(' C4-Value missing for ',A4,I1,'(',F10.3,'). ',
     *              'Quadratic Starkeffect ignored.')
           ENDIF
           gamma_tmp = gamma_tmp + gamma_tmp1   
1119       FORMAT(1X,A4,I2,1X,F6.0,"K",3(1X,F7.4,"eV"),1X,F7.4,
     *            1X,E11.4,1X,F8.4)
          ENDIF
c
C --- GR --------------------------------------------
          gamma_tmp  = line(lc).grad + gamma_tmp 

c
C --- Gamma -----------------------------------------
          gamma(n,lc)= DIV_4PIC * gamma_tmp * line(lc).w / vdop(n,k)

c --- atomic line
C --- peculiar abundance deviations are taken into account (eps_dev)
c     atm.h.eps_dev are supplied by call to lineform.f 
 
         IF (k.LE.NATOM) THEN                
           div_dlamdop(n,lc) = 1./ (vdop(n,k)*line(lc).w)
           kappa_line(n,lc)  =  divtn_FAC_B
     *        * atm.pion(n,line(lc).ion,k)/atm.uion(n,line(lc).ion,k)
     *                         * exp(- elow_k(lc)*divtn)
     *                         * div_dlamdop(n,lc) * gf(lc)
     *                         * atm.h.eps_dev(k)               !n = 1, atm.h.ndepth 
c       IF (n.eq.5.OR.n.eq.55) PRINT*, 'LINE FORMATION 6200 Fe I'
c          PRINT*, 'LTE LINE FORMATION Q-S 6200 Fe I',                     ! for Fe I 6200 A line
c     *   kappa_line(n,lc), line(lc).ion, k, atm.pion(n,line(lc).ion,k),
c     *   atm.uion(n,line(lc).ion,k),elow_k(lc), gf(lc), atm.h.eps_dev(k)
c          WRITE(*,4609)
c     *   n, kappa_line(n,lc),
c     *   atm.pion(n,line(lc).ion,k)*divtn*DIV_K,
c     *   5*exp(- elow_k(lc)*divtn)*atm.pion(n,line(lc).ion,k)
c     *                *divtn*DIV_K/atm.uion(n,line(lc).ion,k)   ! 5 is g_l factor! for the line!!!
c 4609      FORMAT(5X,I2,3(E17.7))
c --- molecular line
C --- peculiar abundance deviations are taken into account (eps_dev)
c
          ELSE                                
           div_dlamdop(n,lc) = 1. / (vdop(n,k) * line(lc).w)
           kappa_line(n,lc)  = (divtn_FAC_B)
     *                       * atm.pmol(n,k-NATOM)/atm.umol(n,k-NATOM)
     *                       * exp(-elow_k(lc)*divtn)
     *                       * div_dlamdop(n,lc) * gf(lc)
     *                       * atm.h.eps_dev(mel1(k-NATOM))
     *                       * atm.h.eps_dev(mel2(k-NATOM))
          ENDIF
c
C --- Multiplication with departure coefficients (NLTE-Lineformation) ---------
          IF (line(lc).rbb_idx .GT. 0) THEN
             x = pop.level(pop.rbb(line(lc).rbb_idx).ll).b(n)
             y = pop.level(pop.rbb(line(lc).rbb_idx).ul).b(n)
             pop.rbb(line(lc).rbb_idx).q(n) = y / x  ! for calc. Sourcefunction
             kappa_line(n,lc) =  kappa_line(n,lc) * x     ! NLTE absorption coefficient 

c MB 2012.07.28
c          if (n.eq.30) then 
***              print*, x,y, pop.rbb(line(lc).rbb_idx).q(n)
c          endif
          ENDIF
         ENDIF                               ! IF (id_used) ..........
        END DO                               ! end wavelength     loop
      END DO                                 ! end id stratification loop

c --- END DEPTH LOOP
c------------------------------------------------------------------------------

      ignore_cnt = 0
      iki        = 1
      k_eta_cont=1
      DO lc = 1,linemax
        IF (id_used(line(lc).id)) THEN  
          iki=iki+1          
          IF (line(lc).id.EQ.1) GOTO 241
          w=line(lc).w
          a = 0.
          v = 0.
          itmp=atm.h.ndepth
          DO n = atm.h.ndepth,1,-1                           
           IF (n .EQ. atm.h.ndepth) THEN
             DO k = k_eta_cont, atm.h.eta_cont_nr - 1
              IF (k .EQ. 1) THEN
                IF (w.LT.atm.eta_cont_w(1)) THEN
                  k_eta_cont=1
                  w_int=0.
                  GOTO 240
                ENDIF
              ELSE IF (k .EQ. atm.h.eta_cont_nr-1) THEN
                IF (w.GT.atm.eta_cont_w(atm.h.eta_cont_nr)) THEN
                  k_eta_cont=atm.h.eta_cont_nr-1
                  w_int=1.
                  GOTO 240
                ENDIF
              ENDIF
              IF ((w.GE.atm.eta_cont_w(k)).AND.
     *            (w.LE.atm.eta_cont_w(k+1))) THEN
                w1 = atm.eta_cont_w(k)         
                w2 = atm.eta_cont_w(k+1)
                k_eta_cont = k
                GOTO 239
              ENDIF 
             END DO
 239         CONTINUE
             w_int = (w - w1) / (w2 - w1)
           ENDIF
 240       CONTINUE
           x = 
     *       (atm.eta_cont(n,k_eta_cont)+atm.etasig(n,k_eta_cont))
     *                                 * (1. -  w_int) + 
     *     (atm.eta_cont(n,k_eta_cont+1) + atm.etasig(n,k_eta_cont+1))
     *                                 *     w_int
           x = x*atm.kappa_ref(n)
           y = kappa_line(n,lc)
           IF (y .GT. x*XKAPMIN1) THEN
              itmp = n
              GOTO 246
           ENDIF
          END DO
          IF (BTEST(lf_ctrl,0)) THEN
           ignore_cnt   = ignore_cnt + 1
           line(lc).use = .false.
          ENDIF    
          GOTO 241
 246      CONTINUE
c
C --- Calculate dopplerwidth dependent line width at innermost border ---------
c
          IF (BTEST(lf_ctrl,0)) THEN
           a = y / x
           y = 3.*sqrt(a*gamma(n,lc)/(XKAPMIN*SQRT_PI))
     *              / div_dlamdop(n,lc)
           x = sqrt(log(a*div_dlamdop(n,lc)/(XKAPMIN*SQRT_PI)))/
     *             div_dlamdop(n,lc)
           a = DMIN1(DMAX1(x,y),0.1*w)
           IF (line(lc).width.LT.a) line(lc).width = 
     *                                DMAX1(a,STD_WIDTH*w)
          ELSE
           IF (line(lc).width.LE. 0.) line(lc).width = STD_WIDTH*w*4.
          ENDIF
C          WRITE(91,5555) itmp,w*1.e8,line(lc).width*1.e8,x*1.e8,y*1.e8
C 5555     FORMAT(' ',I2,2X,F9.3,3X,F7.3,3X,F7.3,3X,F7.3)
 241      CONTINUE
        ENDIF
      END DO
      IF (.NOT.BJTEST(lf_ctrl,26).AND.ignore_cnt.GT.0) 
     *     WRITE(CHN,*) ignore_cnt,' of ',iki,
     *                 ' lines are too faint --> ignored.'

C ------------ C O N T I N U U M - F O R M A T I O N --------------------------

      cont_nr     = 0  
      k_eta_cont  = 1
      IF (BJTEST(lf_ctrl,5)) THEN            ! if residual Flux/Irradiance
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 'calculate Continuum'
      DO kcnt = 1,atm.h.eta_cont_nr 
        w = atm.eta_cont_w(kcnt)
        hcc_2_div_lam5 = HCC_2 / (w**5)
        DO n = 1,atm.h.ndepth                       ! begin stratification loop
          tn   = atm.t(n) 
          divtn= divt(n)
          x    = HC_DIV_K * divtn / w               ! see physcnst.inc
          emi  = EXP(-x)
          em   = 1. - emi
          blam(n) = emi/em                          ! Planck function
          slam(n) = blam(n)                         
          jlam(n) = blam(n)

C --- Continuum opacity (kappa.f)

          eta_cont_lam(n) = 
     *              atm.eta_cont(n,kcnt)
          etasig_lam(n) = 
     *              atm.etasig(n,kcnt)
          eta_lam(n)=  eta_cont_lam(n) + etasig_lam(n)  ! => atm.kappa_ref (kappa.f)
          beta(n)   =  etasig_lam(n) / eta_lam(n)
          fc_lam(n) =  eta_lam(n) * atm.tau(n)       ! for tau_lam - generation
c        WRITE(*,*) atm.kappa_ref(n)
        END DO                                       ! end   stratification loop

C --- Evaluation of tau-continuum(lambda)-scale -----------------------------
  
        xta = ta(6)
        IF (BTEST(lf_ctrl,31)) xta = 1.e30

        IF (tau_integr_method.EQ.1) THEN                    ! LINEAR SCALE  
           CALL NEW_TAU_SCALE(tau_lam,eta_lam,dxtau,atm.h.ndepth,
     *           nn,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
        ELSE                                                ! LOGARITHMIC SCALE
           CALL NEW_TAU_SCALE(tau_lam,fc_lam,dxtau,atm.h.ndepth,nn
     *          ,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
        ENDIF

        IF (BTEST(lf_ctrl,31)) CALL            ! if exact scattering term
     *    ITERATE_SOURCE_FCT(atm.h.ndepth,tau_lam,beta,slam,jlam,ta(6))


C --- Integration over all depth points -------------------------------------

        IF (flux_integr_method.GT.0) THEN                   ! if spline integration
          IF (BJTEST(lf_ctrl,1)) THEN                       ! if Flux calculation
            IF     (flux_integr_method.EQ.2) THEN           ! LOGARITHMIC SCALE
             DO n = 1,atm.h.ndepth                                        
              yint(n) = slam(n) * fc_lam(n) *
     *                  EXPINTEGRAL(tau_lam(n),2) * LN10_2
             END DO
            ELSEIF (flux_integr_method.EQ.1) THEN           ! LINEAR      SCALE
             DO n = 1,atm.h.ndepth
              yint(n) = 2. * slam(n) * eta_lam(n) *
     *                            EXPINTEGRAL(tau_lam(n),2)
             END DO
            ENDIF
          ELSE                                              ! irradiance calc.
            IF     (flux_integr_method.EQ.2) THEN           ! LOGARITHMIC SCALE
             DO n = 1,atm.h.ndepth
              yint(n) = slam(n) * fc_lam(n) * LN10 *
     *             EXP(-tau_lam(n)*div_cos_theta) * div_cos_theta
             END DO
            ELSEIF (flux_integr_method.EQ.1) THEN           ! LINEAR      SCALE
             DO n = 1,atm.h.ndepth
              yint(n) = slam(n) * eta_lam(n) *
     *            EXP(-tau_lam(n)* div_cos_theta) * div_cos_theta
             END DO
            ENDIF
          ENDIF                                             ! flux-/or irradiance calc.
        ENDIF                                               ! if spline integration
        
        IF     (flux_integr_method.EQ.0) THEN               ! gauss quadrature (integration by weights)
          x = 0.
          DO k = 1,6
            x = x + wa(k)*QPOL(ta(k),tau_lam,slam,nn,m)
          END DO
          IF (ta(1).LT.tau_lam(1) ) 
     *      WRITE(CHN,99876) tau_lam(1),ta(1),w*1.e8
        ELSEIF (flux_integr_method.EQ.2) THEN
          x = SPLINE_INTEGRAL(atm.logtau,yint,as,bs,cs,atm.h.ndepth,1) ! spline integr. (log.scale)          
        ELSE      
          x = SPLINE_INTEGRAL(atm.tau,yint,as,bs,cs,atm.h.ndepth,1)    ! spline integr. (linear scale)          
        ENDIF      

        cont_nr          = cont_nr + 1
        cspc(cont_nr).f  = x * hcc_2_div_lam5 
        cspc(cont_nr).w  = w
 
C --- Determine next wavelength point w

        IF (w.GE.wmax) GOTO 1900               ! exit continuum loop
        IF (w*(1.+CRANGE).GT.w2.AND.w.LT.w2) THEN
           w = DMIN1(w2,wmax)                  ! exact continuum point assigned
        ELSE
           w = DMIN1(w*(1. + CRANGE),wmax+dwmin*w*DIV_WLAMREF)
        ENDIF
      ENDDO   ! wavelength loop
      ENDIF


C --- ------------------ L I N E - F O R M A T I O N --------------------------

 1900 CONTINUE     
      IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 'Lineformation'

C --- Startvalues ------

      ilstart          = 1
      ilend            = linemax
      k_eta_cont       = 1                    
      actline          = 0
      firstline        = 0                    
      lastline         = 0                    
      w                = wmin             
      IF (.NOT.BJTEST(lf_ctrl,19)) wlam_cnt = 1  ! wlambda positions not predefined    
      idx_w            = wlam_cnt                                                    
      spc(idx_w).w     = wmin
      iteration_step   = 1
      kcnt             = atm.h.eta_cont_nr ! counter for backgroundopacities
      first_k_eta_cont = 1
      
C --- Determination of nearest line ----------

      DO lc = 1,linemax
        IF (line(lc).w.GT.wmax) GOTO 150
        IF (firstline.EQ.0.AND.line(lc).w.GT.wmin) THEN
          firstline = lc
        ENDIF    
        IF (firstline.GT.0) THEN
          IF (line(lc).w.GT.wmax) GOTO 150
          lastline = lc                                    
        ENDIF  
      END DO
  150 CONTINUE 

      iki=0
      DO lc = firstline,lastline
       IF (line(lc).use) iki=iki+1   
      ENDDO

      IF (iki.GT.0) THEN
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,1151) 
     *        iki,wmin*1.e8,wmax*1.e8
 1151   FORMAT(' ',I5,
     *         ' line center within interval [',F10.2,',',F10.2,']') 
      ELSE
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,1152) 
     *            wmin*1.e8,wmax*1.e8
 1152   FORMAT(' No line center within interval [',F10.2,',',F10.2,']') 
      ENDIF

C --- Check if fluxpoint is already calculated
 2000 CONTINUE                                     ! begin wavelength loop

c        IF (MOD(wlam_cnt,1000).EQ.0 .AND. .NOT. BJTEST(lf_ctrl,26)) 
c     *       WRITE(CHN,*) wlam_cnt,'. wavelength point' 

        hcc_2_div_lam5 = HCC_2 / (w**5)

        DO n = 1,atm.h.ndepth                      ! begin stratification loop
          tn       = atm.t(n)                                
          divtn    = divt(n)
          x        = HC_DIV_K * divtn / w          ! see physcnst.inc
          emi      = EXP(-x)
          em       = 1.D0 - emi
          Blam(n)  = emi/em                        ! Planckfunction
          Jlam(n)  = Blam(n)
          kappa_lam= 0.D0                          ! kappa(lambda)

C --- continuum opacity ----------------------------
          IF (n .EQ. 1) THEN
            DO k = 1, atm.h.eta_cont_nr - 1
              IF ((w.GE.atm.eta_cont_w(k)).AND.
     *            (w.LE.atm.eta_cont_w(k+1))) THEN
                w1 = atm.eta_cont_w(k)         
                w2 = atm.eta_cont_w(k+1)
                k_eta_cont = k
                GOTO 200
              ENDIF 
            END DO
200         CONTINUE
            w_int = (w - w1) / (w2 - w1)
          ENDIF

c          eta_cont_lam(n) = exp(
c     *          alog(atm.eta_cont(n,k_eta_cont)) * (1. -  w_int)
c     *        + alog(atm.eta_cont(n,k_eta_cont+1)) *    w_int)
c          etasig_lam(n)   = exp(
c     *          alog(atm.etasig  (n,k_eta_cont))  * (1. -  w_int)
c     *       + alog(atm.etasig(n,k_eta_cont+1))   *    w_int)

          eta_cont_lam(n) = 
     *              atm.eta_cont(n,k_eta_cont)     * (1. -  w_int)
     *            + atm.eta_cont(n,k_eta_cont + 1) *    w_int
          etasig_lam(n)   = 
     *              atm.etasig  (n,k_eta_cont)     * (1. -  w_int)
     *            + atm.etasig  (n,k_eta_cont + 1) *    w_int

C --- final absorption coefficients
          x = 0.D0
          Slam(n)   = 0.
          kappa_lam = 0.D0
          IF (hyd.AND.w.GE.3.5D-5) THEN
             !CALL HYDLINE(divtn,w,n,x)
             !Slam(n)   = Slam(n)   + x * emi
             !kappa_lam = kappa_lam + x * em

c ---        NEW BALMER LINE TREATMENT                            MOD02.01 fug BEGIN

             !IF (pop.atom.eq.1) THEN

             ! AG65 LTE

!!!
             IF (agbpo.EQ.0.or.w.lt.4.24D-5) THEN       !Mashonkina  March 2008    MOD0601 lm begin             
                IF (fsyesno .EQ. 0) THEN !NO FineStructure
                   IF (agwrong .EQ. 0) THEN 
                   ! the data are available for Ha to 3682, Pa_a to 8392
                      CALL HYDL_AG_CON(divtn,w,n,x)       ! AG Preconvolved profiles
                      IF (w.lt.3.680D-5) THEN         !Mashonkina  March 2008                         
                         CALL HYD_EDGE(divtn,w,n,xxxx)
                        x=x+xxxx
                      endif                             
                   ELSE IF  (agwrong .EQ. 1) THEN                       
                      CALL HYDL_AG(divtn,w,n,x,resscale) ! AG simple adding
                   ENDIF
                ELSE IF (fsyesno .EQ. 1) THEN ! FineStructure
                   ! FS data are available for Ha, Hb, Hg and Hd
                   IF (w.gt.4.0D-5) THEN          !Mashonkina  March 2008                      
                      CALL HYDL_AG_CON_FS(divtn,w,n,x)
                   else
                     x = 0.D0
                     CALL HYD_ESW(divtn,w,n,x)
                   ENDIF
                ENDIF
c             ENDIF
              else
             ! BPO LTE                                            !MOD04.01 fug begin
c             IF (agbpo.EQ.1) THEN
                IF (fsyesno .EQ. 0) THEN                            !NO FineStructure
                   CALL HYDL_BPO_CON(divtn,w,n,x)
                ELSE IF (fsyesno .EQ. 1) THEN ! FineStructure
                   !  data are available for Ha, Hb, Hg                    
                   CALL HYDL_BPO_CON_FS(divtn,w,n,x)
                ENDIF
             ENDIF               ! NLTE-LTE
c             IF (pop.atom.EQ.1) THEN !MOD05.01 fug begin
                jup=2*dsqrt(w/(w-3646.0*1.0d-8))+0.5
             if(jup.le.20) then
                Slam(n)   = Slam(n)   + x * emi*hyddep(jup,n)
                kappa_lam = kappa_lam + x * (hyddep(2,n)-
     &               hyddep(jup,n)*emi)                  
             ELSE
                Slam(n)   = Slam(n)   + x * emi
                kappa_lam = kappa_lam + x * em
             ENDIF              !MOD05.01 fug end   !MOD0601 lm end      
          ENDIF   
          IF (hyd.AND.w.LE.4.0D-5) THEN            
             x = 0.D0
             CALL HYD_ESW(divtn,w,n,x)
             Slam(n)   = Slam(n)   + x * emi
             kappa_lam = kappa_lam + x * em
          ENDIF

C --- Calculation of line opacities for all metallic lines which are concerned
C --- at this wavelenght point with radiative and van-der-Waals-broadening
C --- (factors have been precalculated).        
         
          IF (n.EQ.1) THEN 
            ilstart = 1
            ilend   = linemax
            i1      = 0      ! at n=1: inspect always all available lines
            i2      = linemax
          ELSE IF (n.EQ.2) THEN
            ilstart = MAX(i1,1)
            ilend   = i2
          ENDIF
          DO lc = ilstart,ilend              
           IF (id_used(line(lc).id).AND.line(lc).use) THEN  
            dlam = DABS(w - line(lc).w)
            IF (dlam .LE. line(lc).width) THEN         ! check if line is concerned
             IF (n.EQ.1) THEN
               IF (i1.EQ.0) i1 = lc
               i2 = lc
             ENDIF 
             a = gamma(n,lc)
             v = dlam * div_dlamdop(n,lc)
             x = VOIGT(a,v)
             IF (BJTEST(lf_ctrl,30)) THEN            ! asymmetric lineformation
C               dlamx = line(lc).w - w 
C               x = x + div_dlg(n,lc)*dlamx/(dlamx**2+gam2waals(n,lc))
             ENDIF

             x = kappa_line(n,lc) * x 

             IF (line(lc).rbb_idx.GT.0) THEN            ! NLTE-Line
C ---           NLTE-Correction  B --> S
                emi_s     = emi * pop.rbb(line(lc).rbb_idx).q(n)
                em_s      = 1.D0 - emi_s
                y         = x * em_s
                Slam(n)   = Slam(n)   + x * emi_s 
                kappa_lam = kappa_lam + y
             ELSE
C ---           LTE-Case: S = B  <- Sourcefunction assumed to be the Planckfunction
                Slam(n)   = Slam(n)   + x * emi
                kappa_lam = kappa_lam + x * em
             ENDIF  
            ENDIF
           ENDIF                             ! IF (id_used) ......
          END DO                             ! endloop of calculation metallic absorption


          IF (BJTEST(lf_ctrl,7)) THEN        ! contribution function
           tau_lam(n) = 0.
           IF (BJTEST(lf_ctrl,1)) THEN       ! flux
             eta_line_vec(n) = kappa_lam / atm.kappa_ref(n)
             eta_cf(n)= eta_line_vec(n)
     *                   + eta_cont_lam(n) + etasig_lam(n)
           ELSE                                     ! intensity
             eta_line_vec(n) = kappa_lam / atm.kappa_ref(n)
             eta_cf(n)= eta_line_vec(n)     ! = k_r
     *                 + (eta_cont_lam(n)+etasig_lam(n)) 
     *                   * hcc_2_div_lam5 * Blam(n) / cspc(1).f
           ENDIF
           IF (ABS(eta_line_vec(n)).GT. 1.e-30) THEN
               Sline(n) = Slam(n) / kappa_lam       ! line source function
           ELSE
               Sline(n)        = 0.
               eta_line_vec(n) = 0.
           ENDIF
           ta(6)      = 1.E4                        ! big dummy value
          ENDIF

          IF (BJTEST(lf_ctrl,31)) THEN    
             eta_lam(n)= kappa_lam / atm.kappa_ref(n) + 
     *                     eta_cont_lam(n)
             Slam(n)   = Slam(n)   / atm.kappa_ref(n) +
     *                     eta_cont_lam(n) * Blam(n) ! add 
             Slam(n)   = Slam(n)   / eta_lam(n)     ! now true absorption S
             eta_lam(n)= eta_lam(n) + etasig_lam(n)
             beta(n)   =  etasig_lam(n) / eta_lam(n)
             strue(n)  = slam(n)
          ELSE                            
             eta_lam(n)= kappa_lam / atm.kappa_ref(n) + 
     *                   eta_cont_lam(n) + etasig_lam(n) 
             Slam(n)   = ( Slam(n)/atm.kappa_ref(n)
     *                 + (eta_cont_lam(n)+etasig_lam(n))*Blam(n))
             Slam(n)   = Slam(n) / eta_lam(n) 
          ENDIF    
          
          f_lam(n) = atm.tau(n) * eta_lam(n) ! for further use
        END DO                               ! end   stratification loop

C --- Evaluation of tau-line(lambda)-scale = monochromatic tau-scale      
        xta = ta(6)
        IF (BTEST(lf_ctrl,31)) xta = 1.e30
        IF (BJTEST(lf_ctrl,7)) THEN          ! contributionfunction
          DO n=1,atm.h.ndepth
            f_lam(n) = atm.tau(n) * eta_cf(n)
          ENDDO
          IF (tau_integr_method.EQ.1) THEN                   ! LINEAR SCALE
           CALL NEW_TAU_SCALE(tau_lam,eta_cf,dxtau,atm.h.ndepth,
     *         nn,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
          ELSE                                               ! LOGARITHMIC SCALE  
           CALL NEW_TAU_SCALE(tau_lam,f_lam,dxtau,atm.h.ndepth,
     *         nn,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
          ENDIF
          DO n=1,atm.h.ndepth
            IF (BJTEST(lf_ctrl,1)) THEN       ! flux
              Sline(n)   = Sline(n)  * atm.tau(n) * eta_line_vec(n)
            ELSE                                     ! intensity
              Sline(n)   = eta_line_vec(n) * atm.tau(n) * 
     *                    (1.-Sline(n)/cspc(1).f)    ! = S_line/I_c
            ENDIF
          ENDDO
        ELSE
          IF (tau_integr_method.EQ.1) THEN                   ! LINEAR SCALE
           CALL NEW_TAU_SCALE(tau_lam,eta_lam,dxtau,atm.h.ndepth,
     *         nn,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
          ELSE                                               ! LOGARITHMIC SCALE  
           CALL NEW_TAU_SCALE(tau_lam,f_lam,dxtau,atm.h.ndepth,
     *         nn,xta,tau_integr_method,atm.h.nscale,atm.h.nbreak)
          ENDIF
        ENDIF    

        IF (BTEST(lf_ctrl,31)) THEN                   ! exact scattering term ? 
C ---      correction of true absorption sourcefunction including exact scattering term
           CALL ITERATE_SOURCE_FCT(atm.h.ndepth,tau_lam,beta,slam,
     *                            jlam,ta(6))
        ENDIF
        

C --- Integration over all depth points ----------------------------------------
C --- --------------------------------------------------------------------------
                   
C --- Line(lambda) -------------------------------------------------------------

        IF (BJTEST(lf_ctrl,7)) THEN          ! output of contribution function
          wlam_cnt = NR_SPLINE
          xx        = atm.logtau(1)
          dx       = (atm.logtau(atm.h.ndepth)-xx)/(NR_SPLINE-1)
          IF (BJTEST(lf_ctrl,1)) THEN
            DO n = 1,atm.h.ndepth
              IF (tau_lam(n).GT.0.) THEN
                 yint(n)= Sline(n) * EXPINTEGRAL(tau_lam(n),2) * LN10
              ELSE
                 yint(n) = 0.
              ENDIF
            END DO
            IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(chn,*)
     *                               ' Flux Contribution Function'
          ELSE
            DO n = 1,atm.h.ndepth
              IF (tau_lam(n).GT.0.) THEN
                 yint(n) = Sline(n)  * EXP(-tau_lam(n)*div_cos_theta)
     *                      *  LN10 * div_cos_theta
              ELSE
                 yint(n) = 0.
              ENDIF
            END DO
            IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(chn,*)
     *                   ' Intensity Contribution Function'
          ENDIF
          y=SPLINE_INTEGRAL(atm.logtau,yint,as,bs,cs,atm.h.ndepth,1)
          IF (BTEST(lf_ctrl,18)) THEN 
            IF (y.GT.1.e-33) THEN
             CALL SPLINE(atm.logtau,atm.t,as,bs,cs,atm.h.ndepth)      
             IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(chn,*) '<S/B>,<tau>,<T>'
             DO n = 1,atm.h.ndepth
               S_B(n)   = Slam(n) / Blam(n)
               f_lam(n) = yint(n) * S_B(n)
               dtdtau(n)= as(n)
             END DO

c             n = atm.h.ndepth
c             dxtau(n) = atm.logtau(n)-atm.logtau(n-1)
c             dtdtau(n)= as(n-1) + dxtau(n)*(2.*bs(n-1) 
c     *                    + 3.*cs(3)*dxtau(n))

             spc(1).w = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                                  cs,atm.h.ndepth,1) / y
             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * S_B(n) * S_B(n)
             END DO
             spc(1).f = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                   cs,atm.h.ndepth,1) / y - spc(1).w**2 
             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * atm.logtau(n)
             END DO
             spc(2).w = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                                  cs,atm.h.ndepth,1) / y
             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n)*atm.logtau(n)*atm.logtau(n)
             END DO
             spc(2).f = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                   cs,atm.h.ndepth,1) / y - spc(2).w**2 

             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * atm.t(n)
             END DO
             spc(3).w = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                                  cs,atm.h.ndepth,1) / y
             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * atm.t(n) * atm.t(n)
             END DO
             spc(3).f = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                  cs,atm.h.ndepth,1) / y - spc(3).w**2

             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * dtdtau(n)
             END DO
             spc(4).w = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                                  cs,atm.h.ndepth,1) / y
             DO n = 1,atm.h.ndepth
               f_lam(n) = yint(n) * dtdtau(n) * dtdtau(n)
             END DO
             spc(4).f = SPLINE_INTEGRAL(atm.logtau,f_lam,as,bs,
     *                  cs,atm.h.ndepth,1) / y - spc(4).w**2

            ELSE
             spc(1).w = 0.   ! mean of S/B
             spc(2).w = 0.   ! mean of tau
             spc(3).w = 0.   ! mean of t
             spc(4).w = 0.   ! mean of dt/dlogtau
            ENDIF 
            IF (spc(1).f .LT. 0.) spc(1).f = 0.
            IF (spc(2).f .LT. 0.) spc(2).f = 0.
            IF (spc(3).f .LT. 0.) spc(3).f = 0.
            IF (spc(4).f .LT. 0.) spc(4).f = 0.

            call SPLINE(atm.logtau,S_B,as,bs,cs,atm.h.ndepth)
            wlam_cnt = NR_SPLINE+4 
            DO n=5,wlam_cnt
             spc(n).w = xx
             spc(n).f = SPOL(xx,atm.logtau,S_B,as,bs,cs,atm.h.ndepth)
             xx = xx + dx
            END DO
            wmin = spc(5).w
            wmax = spc(wlam_cnt).w
          ELSE
           DO n=1,NR_SPLINE
            spc(n).w = xx
            spc(n).f = SPOL(xx,atm.logtau,yint,as,bs,cs,atm.h.ndepth)
            xx = xx + dx
           END DO
           wmin = spc(1).w
           wmax = spc(wlam_cnt).w
          ENDIF
          GOTO 4000                                        ! finish program
        ENDIF

        IF (flux_integr_method.GT.0) THEN                  ! if spline integration
         IF (BJTEST(lf_ctrl,1)) THEN          ! flux calculation
          IF     (flux_integr_method.EQ.2) THEN            ! LOGARITHMIC SCALE
           DO n = 1,atm.h.ndepth
            yint(n)= Slam(n) * f_lam(n) * 
     *               EXPINTEGRAL(tau_lam(n),2) * LN10_2
           END DO                                                 
          ELSE                                             ! LINEAR      SCALE
           DO n = 1,atm.h.ndepth                 
            yint(n)= 2.*Slam(n)*eta_lam(n)*EXPINTEGRAL(tau_lam(n),2)
           END DO
          ENDIF
         ELSE                                 ! irradiance calculation
          IF     (flux_integr_method.EQ.2) THEN            ! LOGARITHMIC SCALE
           DO n = 1,atm.h.ndepth
            yint(n)=Slam(n) * f_lam(n) *  LN10 * div_cos_theta *
     *              EXP(-tau_lam(n)*div_cos_theta)
           END DO
          ELSE                                             ! LINEAR      SCALE
           DO n = 1,atm.h.ndepth      
            yint(n)=Slam(n)*eta_lam(n)*EXP(-tau_lam(n)*div_cos_theta)*
     *              div_cos_theta
           END DO
          ENDIF
         ENDIF   ! flux-/or irradiance calc.
        ENDIF ! if spline integration                                    
       
        IF     (flux_integr_method.EQ.0) THEN              ! gauss integration by weights
          x = 0.
          DO k = 1,6
            x = x + wa(k)*QPOL(ta(k),tau_lam,Slam,nn,m)
          END DO
*
*
          IF (ta(1).LT.tau_lam(1) ) 
     *      WRITE(CHN,99876) tau_lam(1),ta(1),w*1.e8
99876     FORMAT('WARNING: 1.depthpoint (',G12.5,') > (',
     *            G12.5,') at ',F11.3,' A')
*
*
        ELSE                                               ! spline integration
         IF (flux_integr_method.EQ.2) THEN                 ! spline integr. (log.scale)  
          x = SPLINE_INTEGRAL(atm.logtau,yint,as,bs,cs,atm.h.ndepth,1) 
         ELSE                                              ! spline integr. (linear scale)          
          x = SPLINE_INTEGRAL(atm.tau,yint,as,bs,cs,atm.h.ndepth,1)   
         ENDIF

        ENDIF
*
*
* MB 31.10.11
        
        if (ABS(w*1.e8-wline_center).lt.1.e-3) then
         print*, 'NLTE core: ',w*1.e8
         DO n = 1,atm.h.ndepth
          print*,'tau_cont(500): ', atm.logtau(n), 
     *           'tau_lam(n): ', tau_lam(n)
         ENDDO
        endif

c        WRITE(*,*) 'After spline integ of x:'
c        write(*,*) atm.logtau
c        write(*,*) x
        spc(idx_w).f = x * hcc_2_div_lam5

        IF (iteration_step.GT.max_iteration_step) GOTO 3000 ! finish wavel.-loop

C --- Determination of next wavelength point w --------------------------------
C --- spectrum is calculated in 4 steps
C ---
C --- 1) Calculation of spectrum at start-, end-point and line centers
C        and at continuum points
C --- 2) Calculation of spectrum at each distance dwmax (from right to left)
C ---    IF BIT 27 OF LF_CTRL SET ---> FINISHED HERE
C --- 3) Calculation "   "   "   from left to right, starting at each line
C ---    center from line(i) to 0.5(line(i)+line(i+1)
C --- 4) Calculation "   "   "   from right to left, starting at each line
C ---    center from line(i) to 0.5(line(i)+line(i-1)
C ---       
C --- (If bit 19 SET:----> get next w from wlam_out(i) --> only one iterationstep)
        
2500  CONTINUE               
      IF (wlam_cnt.GE.NWMAX) GOTO 4000                                        ! finish program
      IF (SECNDS(old_tim).GT.60.) THEN
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(chn,*)
     *        wlam_cnt,'. wavelength point (',w*1.D8,') Step: ',
     *        iteration_step 
          old_tim = SECNDS(0.0)
      ENDIF

      IF (BJTEST(lf_ctrl,19))  THEN   ! wlambda positions already predefined
        IF (wlam_cnt_out.LE.wlam_cnt) GOTO 3000 ! finish and output
        wlam_cnt     = wlam_cnt + 1
        idx_w        = idx_w    + 1
        w            = wlam_out(idx_w)*1.D-8    ! wlam_out in [A] !
        spc(idx_w).w = w
      ELSE
       IF (iteration_step.EQ.1) THEN
        IF (firstline.EQ.0) THEN
          w = wmax
          iteration_step = 2                    ! next iteration step
        ELSE
          IF (wlam_cnt.EQ.1) THEN
            actline = firstline      
            w = line(actline).w
          ELSE
            actline = actline + 1
            IF (actline.GT.lastline) THEN
              w = wmax
              iteration_step = 2               ! next iteration step
            ELSE
              w = DMIN1(line(actline).w,wmax)  ! line center
            ENDIF
          ENDIF
        ENDIF
       ELSEIF (iteration_step.EQ.2) THEN
        IF (BJTEST(lf_ctrl,27)) THEN
           w = w*(1.D0 - dwmin*DIV_WLAMREF)
        ELSE
           w = w*(1.D0 - dwmax*DIV_WLAMREF)
        ENDIF
        IF (kcnt.GT.0.AND.w.LE.atm.eta_cont_w(kcnt)) THEN
           w    = atm.eta_cont_w(kcnt)
           kcnt = kcnt - 1
        ENDIF
        IF (w.LE.wmin) THEN
          IF (BJTEST(lf_ctrl,27)) THEN
             IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *               'Fast lineformation finished'
             GOTO 3000            
          ENDIF
          IF (firstline.EQ.0) THEN                 ! if no line within interval
C ---       if no line in internal linelist exist, finish wavelength loop
            IF (linemax.LE.0) GOTO 3000            ! finished
            w = wmin + dwmin*w*DIV_WLAMREF
            IF (DABS(w-wmax).GE.dw_eps*w*DIV_WLAMREF) GOTO 3000   ! finished

C ---       check if line exist with linecenter < wmin
            i = 0
            DO lc = 1,linemax
              IF (line(lc).w .LT. wmin) THEN
                i = lc
              ELSE
                GOTO 2505
              ENDIF 
            ENDDO
            i = linemax
 2505       CONTINUE

            IF (i.EQ.0) THEN
              w      = wmax - dwmin*w*DIV_WLAMREF
              xstart = wmax
              xend   = wmin
              iteration_step = 4
            ELSE
C ---         check if line exist with linecenter > wmax
              xend   = line(i).w
              i    = 0
              DO lc = linemax,1,-1
                IF (line(lc).w .GT. wmax) THEN
                  i = lc
                ELSE
                  GOTO 2507
                ENDIF 
              ENDDO
              i = 1
 2507         CONTINUE
              IF (i.EQ.0) THEN
                xend   = wmax
              ELSE
                xend   = DMIN1(wmax,0.5*(xend  +line(i).w))                
              ENDIF
              iteration_step = 3
              xstart = wmin
            ENDIF
          ELSE                    
            actline = firstline 
            w      = line(firstline).w + dwmin*w*DIV_WLAMREF        ! starting at 1.line center
            xend   = wmax
            IF (firstline+1.LE. linemax) THEN
              yend   = (0.5*(line(firstline).w+line(firstline+1).w))
              IF (y.LT.wmax) xend   = yend  
            ENDIF
            
C ---       Search nearest line at left side of wmin
            IF (actline.GT.1) THEN
              IF((0.5*(line(actline-1).w+line(actline).w)).GE.wmin) THEN
                w    = wmin + dwmin*w*DIV_WLAMREF                 ! starting at wmin
                xend   = (0.5*(line(actline-1).w+line(actline).w))
              ENDIF
            ENDIF
            iteration_step = 3
            xstart = wmin
          ENDIF                
        ENDIF         
       ELSEIF (iteration_step.EQ.3) THEN
        IF (w.EQ.xend  ) THEN             ! goto next line
          IF (firstline.GT.0) THEN 
            IF (xend  .GE. line(firstline).w) THEN
              actline = actline + 1
            ENDIF
                 
            IF (actline.GT.lastline) THEN
              actline = lastline
              xstart  = wmax
              IF (xend  .LT.wmax) THEN
                w = wmax - dwmin*w*DIV_WLAMREF   ! xend   remains unchanged
              ELSE               
                w = line(actline).w - dwmin*w*DIV_WLAMREF
                IF (lastline.GT.1) THEN
                  xend   = (0.5*(line(actline-1).w+line(actline).w))
                ELSE
                  xend   = wmax
                ENDIF 
              ENDIF
              iteration_step = 4
            ELSE                            ! actline <= lastline
              xstart = line(actline).w
              IF (actline.GE.linemax) THEN
                xend   = wmax
              ELSE
                yend = (0.5*(xstart + line(actline+1).w))
                xend   = DMIN1(yend,wmax)
              ENDIF
              w      = DMIN1(xstart + dwmin*w*DIV_WLAMREF,xend)
            ENDIF
          ELSE                           ! firstline = 0, w > xend  
            IF (xend.GE.wmax) GOTO 3000       ! finish
            xstart = wmax
            w = xstart - dwmin*w*DIV_WLAMREF           ! xend remains unchanged
            iteration_step = 4
          ENDIF                          ! IF (firstline.GT.0)
        ELSE
          IF (DABS(w-xend).LT.dw_eps*w*DIV_WLAMREF) THEN
            w = xend
            GOTO 2550                    ! prepare next iteration step     
          ENDIF
          IF (spc(idx_w).w - xstart .LE. 
     *             (dwmin+dw_eps)*w*DIV_WLAMREF) THEN 
            w = DMIN1(spc(idx_w).w + dwmin*w*DIV_WLAMREF,xend)              
            GOTO 2550                    ! prepare next iteration step     
          ELSEIF ((spc(idx_w).w-spc(idx_w-1).w).GE.
     *             (dwmin+dw_eps)*w*DIV_WLAMREF) THEN
C ---       Check current stepwidth
#if SWCTRL==2
            wsp(1)  = spc(idx_w - 3).w
            wsp(2)  = spc(idx_w - 2).w
            wsp(3)  = spc(idx_w - 1).w
            wsp(4)  = spc(idx_w    ).w

            fsp(1)  = DBLE(spc(idx_w - 3).f)
            fsp(2)  = DBLE(spc(idx_w - 2).f)
            fsp(3)  = DBLE(spc(idx_w - 1).f)
            fsp(4)  = DBLE(spc(idx_w    ).f)
C ---       normalization
            rtmp = 0.25D0*(fsp(1)+fsp(2)+fsp(3)+fsp(4))
            IF (rtmp.GT.small) THEN
             DO i = 1,4
               fsp(i) = fsp(i) / rtmp
             ENDDO
            ENDIF
C ---       Check if next stepwidth can be in- or decremented
            IF (CHECK_STEPWIDTH_CRIT(wsp,fsp,sw_crit)) THEN
#else
C ---       Check if next stepwidth can be in- or decremented        
            IF (CHECK_STEPWIDTH_CRIT(atm.h.ndepth,eta_lam,
     *                           eta_cont_lam,tau_lam,sw_crit)) THEN
#endif
C ---         devide stepsize
              w = DMAX1(spc(idx_w - 1).w+dwmin*w*DIV_WLAMREF,
     *                  0.5D0*(spc(idx_w - 1).w + spc(idx_w).w))
            ELSE       ! double stepsize
              IF ((idx_w.LT.wlam_cnt) .AND. 
     *             (     (spc(idx_w+1).w - spc(idx_w).w) .LE.
     *              2.D0*(spc(idx_w).w - spc(idx_w-1).w)+
     *                    dw_eps*w*DIV_WLAMREF)) THEN
                  idx_w = idx_w + 1
                  w     = spc(idx_w).w 
                  GOTO 2500
              ELSE
C ---           double stepwidth
                w=DMIN1(spc(idx_w).w+2.D0*(spc(idx_w).w-spc(idx_w-1).w),
     *                  xend)
              ENDIF
            ENDIF
          ELSE         ! IF ((spc(idx_w).w-spc(idx_w-1).w) .LE. (dwmin+dw_eps)*w*DIV_WLAMREF)
            w = DMIN1(spc(idx_w).w + 2.D0 * dwmin*w*DIV_WLAMREF,xend)
          ENDIF  
        ENDIF

 2550   CONTINUE
        IF (w.GE.wmax) THEN     
          actline = lastline 
          xstart  = line(lastline).w 
          xend    = wmin
          w       = xstart - dwmin*w*DIV_WLAMREF      ! starting at last line center
          IF (lastline.GT.1) THEN
            yend = (0.5*(xstart + line(lastline-1).w))
            IF (y.GT.wmin) xend = yend
          ENDIF

C ---     Look for nearest line at right side of wmax
          IF (actline.LT.linemax) THEN
           IF((0.5*(xstart + line(actline+1).w)).LE.wmax) THEN
            xstart = wmax
            w      = xstart - dwmin*w*DIV_WLAMREF                 ! starting at wmax
            xend   = (0.5*(line(actline).w + line(actline+1).w))
           ENDIF
          ENDIF
          iteration_step = 4
        ENDIF
       ELSEIF (iteration_step.EQ.4) THEN   ! now iterate from right to left
        IF (w.EQ.xend  ) THEN              ! goto next line
          IF (firstline.GT.0) THEN 
            IF (xend  .LE. line(lastline).w) THEN
              actline = actline - 1
            ENDIF
                 
            IF (actline.LT.firstline) THEN
              GOTO 3000                     ! finish
            ELSE                            ! actline <= lastline
              xstart = line(actline).w
              IF (actline.LE.1) THEN
                xend   = wmin
              ELSE
                yend = (0.5*(xstart + line(actline-1).w))
                xend   = DMAX1(yend,wmin)
              ENDIF
              w      = DMAX1(xstart - dwmin*w*DIV_WLAMREF,xend)
            ENDIF
          ELSE                           ! firstline = 0, w > xend  
            GOTO 3000                    ! finish
          ENDIF                          ! IF (firstline.GT.0)
        ELSE
          IF (DABS(w-xend).LT.dw_eps*w*DIV_WLAMREF) THEN
            w = xend
            GOTO 2650                    ! prepare next iteration step     
          ENDIF
          IF (xstart - spc(idx_w).w .LE. 
     *               (dwmin+dw_eps)*w*DIV_WLAMREF) THEN 
            w = DMAX1(spc(idx_w).w - dwmin*w*DIV_WLAMREF,xend)              
            GOTO 2650                    ! prepare next iteration step     
          ELSEIF ((spc(idx_w+1).w-spc(idx_w).w).GE.
     *              (dwmin+dw_eps)*w*DIV_WLAMREF) THEN
C ---       Check current stepwidth
#if SWCTRL==2

            wsp(1)  = spc(idx_w + 3).w
            wsp(2)  = spc(idx_w + 2).w
            wsp(3)  = spc(idx_w + 1).w
            wsp(4)  = spc(idx_w    ).w

            fsp(1)  = DBLE(spc(idx_w + 3).f)
            fsp(2)  = DBLE(spc(idx_w + 2).f)
            fsp(3)  = DBLE(spc(idx_w + 1).f)
            fsp(4)  = DBLE(spc(idx_w    ).f)

C ---       normalization
            rtmp = 0.25D0*(fsp(1)+fsp(2)+fsp(3)+fsp(4))
            IF (rtmp.GT.small) THEN
             DO i = 1,4
               fsp(i) = fsp(i) / rtmp
             ENDDO
            ENDIF
C ---       Check if next stepwidth can be in- or decremented
            IF (CHECK_STEPWIDTH_CRIT(wsp,fsp,sw_crit)) THEN
#else
C ---       Check if next stepwidth can be in- or decremented        
            IF (CHECK_STEPWIDTH_CRIT(atm.h.ndepth,eta_lam,
     *                           eta_cont_lam,tau_lam,sw_crit)) THEN
#endif
C ---         devide stepsize
              w = DMAX1(spc(idx_w + 1).w - dwmin*w*DIV_WLAMREF,
     *                  0.5D0*(spc(idx_w + 1).w + spc(idx_w).w))
            ELSE       ! double stepsize
              IF ((idx_w.GT.1) .AND. 
     *             (     (spc(idx_w).w - spc(idx_w-1).w) .LE.
     *              2.D0*(spc(idx_w+1).w - spc(idx_w).w)+
     *                    dw_eps*w*DIV_WLAMREF)) THEN
                  idx_w = idx_w - 1
                  w     = spc(idx_w).w 
                  GOTO 2500             
              ELSE
C ---           double stepwidth
                w=DMAX1(spc(idx_w).w + 2.D0*(spc(idx_w).w
     *                                       - spc(idx_w+1).w),xend)
              ENDIF
            ENDIF
          ELSE         ! IF ((spc(idx_w+1).w-spc(idx_w).w) .LE. (dwmin+dw_eps)*w*DIV_WLAMREF)
            w = DMAX1(spc(idx_w).w - 2.D0 * dwmin*w*DIV_WLAMREF,xend)
          ENDIF  
        ENDIF            

 2650   CONTINUE
        IF (w.LE.wmin) GOTO 3000
       ELSE
        GOTO 3000 ! finish
       ENDIF

C ---  sort next wavelength point into array 'wlam' ---------------------------
C ---  search for insert index
       IF (iteration_step.EQ.1.OR.iteration_step.EQ.3) THEN   ! (left --> right)
        DO i = wlam_cnt,1,-1
          IF (wlam_cnt.GT.1.AND.DABS(w-spc(i).w).LE.
     *        dw_eps*w*DIV_WLAMREF) THEN
             idx_w = i
C            WRITE(CHN,*) 'Same wavelength ',w,spc(i).w
            GOTO 2500                         ! wavelength already exist
                                              ! assign new wavelength
          ENDIF
          IF (w.GT.spc(i).w) THEN 
            idx_w = i + 1
            GOTO 2800 
          ENDIF
        END DO
       ELSE                                                    ! (right --> left)
        DO i = 1,wlam_cnt
          IF ( iteration_step.EQ.4.AND.
     *         DABS(w-spc(i).w).LE.dw_eps*w*DIV_WLAMREF) THEN
            idx_w = i
C            WRITE(CHN,*) 'Same wavelength ',w,spc(i).w
            GOTO 2500                         ! wavelength already exist
                                              ! assign new wavelength
          ENDIF
          IF (w.LT.spc(i).w) THEN 
            idx_w = i 
            GOTO 2800 
          ENDIF
        END DO   
        idx_w = wlam_cnt + 1
        GOTO 2800 
       ENDIF
       WRITE(CHN,*) 'ERROR AT LABEL 2800'
       lf_status = JIBSET(lf_status,0)         ! Fatal Error
       GOTO 8999

C ---  revise spc array
 2800  CONTINUE
       IF (idx_w.LE.wlam_cnt) THEN
        DO i = wlam_cnt,idx_w,-1
           spc(i+1) = spc(i)

        ENDDO
       ENDIF

       spc(idx_w).w = w
       wlam_cnt     = wlam_cnt + 1

      ENDIF       ! bit 19 not set
C --- ------------------------------------------------------------------

      GOTO 2000                                ! next iteration

C     ------------------------ Output -----------------------------------------
3000  CONTINUE
      IF (wlam_out(wlam_cnt) .EQ. wlam_out(wlam_cnt-1)) THEN
        wlam_cnt = wlam_cnt-1
      ENDIF

      IF (BJTEST(lf_ctrl,5)) THEN              ! residual flux/irradiance ? 
        wlam_cont_1 = cspc(1).w
        cont_flux_1 = cspc(1).f
        wlam_cont_2 = cspc(2).w
        cont_flux_2 = cspc(2).f
        div_fac     = 1./(wlam_cont_2 - wlam_cont_1)

        n = 2               
        DO i = 1,wlam_cnt
          w = spc(i).w
          IF (w.GT.wlam_cont_2) THEN
            j = n                              ! last index of right freq.point
            DO k = j,cont_nr - 1
              IF (w .GT. cspc(k).w   .AND.
     *            w .LE. cspc(k+1).w) THEN 
                n = k+1
                wlam_cont_1 = wlam_cont_2
                cont_flux_1 = cont_flux_2
                wlam_cont_2 = cspc(n).w
                cont_flux_2 = cspc(n).f
                div_fac     = 1./(wlam_cont_2 - wlam_cont_1)
                GOTO 3100
              ENDIF
            END DO 
          WRITE(CHN,*) 'WARNING: w > wlam_cont_2 !  w= ',
     *                   w*1.D8,' wlam_cont_2= ',wlam_cont_2*1.e8
 3100       CONTINUE            
          ENDIF
          w_int   = (wlam_cont_2 - w) * div_fac
          spc(i).f = spc(i).f / (  w_int    * cont_flux_1 + 
     *                          (1. - w_int) * cont_flux_2  )

c          spc(i).f = spc(i).f / exp(  w_int * alog(cont_flux_1) + 
c     *                          (1. - w_int) * alog(cont_flux_2))


        END DO
      ENDIF                                

4000  CONTINUE
C --- providing header ---------------------------------------------------------
      IF (BJTEST(lf_ctrl,7)) GOTO 4400       ! output of contributionfunction

      head.lin.w_start = wmin*1.e8           ! minimum wavelength
      head.lin.w_end   = wmax*1.e8           ! max.    wavelength

      lf_output_filename(12:12) = '_'
      WRITE(lf_output_filename(13:18),'(I6.6)') NINT(head.lin.w_start)
      lf_output_filename(19:19) = '-'
      WRITE(lf_output_filename(20:25),'(I6.6)') NINT(head.lin.w_end)
      
      head.id          = lf_output_filename(1:25)! Spektrum-id
      head.lin.rdst    = -1.                 ! not equidistant
      head.lin.w_nr    = wlam_cnt            ! Nr of points
                   
      head.lin.eps_nr = head.atm.eps_nr      ! Nr of abundance deviations

      head.lin.eps(1) = 0.
      head.lin.eps(2) = 0.
      DO i = 3,head.lin.eps_nr
        IF (atm.h.eps(i).GT.0. .AND. (atm.h.eps_dev(i).NE.1..OR.
     *     (ABS(atm.h.eps(i)-epss(i)-atm.h.z).GT.(0.001)))) THEN 
          head.lin.eps(i) = log10(atm.h.eps_dev(i))
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,1987) 
     *  cel(i),epss(i),atm.h.eps(i) + 
     *     log10(atm.h.eps_dev(i))-epss(i)-atm.h.z
 1987     FORMAT(' ',A2,': ',F6.2,' dex solar --> ',F6.2,
     *    ' dex peculiar deviation')
        ELSE
          head.lin.eps(i) = 0.
        ENDIF
      END DO

      CALL date_and_time(buf_date)      

      head.lin.date = buf_date
      env_len = IGETENV('USER',env)
      IF (env_len.GT.0) head.lin.user = env(1:env_len)                      

      k = 1
      env_len = IGETENV('LINEFORM_SBR',env)
      IF (env_len.GT.0) THEN
         head.lin.ref = 'LF_SBR='//env(1:env_len)
         CALL STRTRIM(head.lin.ref,k,nwords)
      ENDIF

      IF (BJTEST(lf_ctrl,20)) THEN
        env_len = IGETENV('POP_FILE',env)
        IF (env_len.GT.0) THEN
           head.lin.ref = head.lin.ref(1:k)
     *                  //' DEPART='//env(1:env_len)
           CALL STRTRIM(head.lin.ref,k,nwords)
        ENDIF
      ENDIF

      IF (BJTEST(lf_ctrl,24)) THEN
        env_len = IGETENV('LINTAB_LOOKUP',env)
        IF (env_len.GT.0) THEN
           head.lin.ref = head.lin.ref(1:k)
     *                  //' LOOKUP='//env(1:env_len)
           CALL STRTRIM(head.lin.ref,k,nwords)
        ENDIF
      ENDIF

      IF (BJTEST(lf_ctrl,16)) THEN
        env_len = IGETENV('LINFOR_MICROTURB',env)
        IF (env_len.GT.0) THEN
           head.lin.ref = head.lin.ref(1:k)
     *                  //' VARTURB='//env(1:env_len)
           CALL STRTRIM(head.lin.ref,k,nwords)
        ENDIF
      ENDIF


      head.lin.ref = head.lin.ref(1:k)
     *                  //' LINES='//linedata_file
      CALL STRTRIM(head.lin.ref,k,nwords)

      IF (BJTEST(lf_ctrl,22)) THEN
         env_len = IGETENV('ATMOSPHERE',env)       
         IF (env_len.GT.0) THEN 
           head.lin.ref = head.lin.ref(1:k)
     *                   //' ATM='//env(1:env_len)
           CALL STRTRIM(head.lin.ref,k,nwords)
         ENDIF
      ENDIF

      env_len = IGETENV('LF_CODE',env)
      IF (env_len.GT.0) THEN
         head.lin.ref = head.lin.ref(1:k)
     *                  //' LF_CODE='//env(1:env_len)
         CALL STRTRIM(head.lin.ref,k,nwords)
      ENDIF

      CALL STRTRIM(head.lin.ref,k,nwords)
      head.lin.ref = head.lin.ref(1:k)
      head.lin.ref_nr = k
      head.lin.cmt = 'BF='//bf_info
      CALL STRTRIM(head.lin.cmt,k,nwords)
      head.lin.cmt_nr = k                    ! Nr of comment character
      head.lin.status = 0
C --- --------------------------------------------------------------------------
 4400 CONTINUE
      x = 1.D8
      IF (BJTEST(lf_ctrl,7)) x=1.            ! if yint-output 

      IF (BJTEST(lf_ctrl,6)) THEN ! file output ? 
        lf_output_filename(22:26) = '.tspc'
        IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *            'Spectrum ===> ',lf_output_filename(1:26)
        env_len = IGETENV('SIU_MAIN',env)
        OPEN(UNIT = 30, FILE = env(1:env_len)
     *       //'/out/'//lf_output_filename, 
     *       STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
        WRITE(30,*) head.id
        WRITE(30,4500) head.atm.teff,head.atm.logg,head.atm.z,
     *                 head.atm.alpha
 4500   FORMAT(' Teff:',I5,' log(g):',F5.2,' [Fe/H]:',F5.2,
     *         ' l/Hp:',F5.2)
        WRITE(30,*) head.atm.cmt(1:head.atm.cmt_nr)
        WRITE(30,*) wlam_cnt                                    
        DO i = 1,wlam_cnt                               
         WRITE(30,'(F14.5,3x,E15.7)') spc(i).w * x,spc(i).f
        END DO
        CLOSE(30)
      ENDIF

      head_i_out(1) = MIN(head.atm.ref_nr,LEN(head_atmref))
      head_i_out(2) = head.atm.eps_nr
      head_i_out(3) = head.atm.status
      head_i_out(4) = MIN(head.atm.cmt_nr,LEN(head_atmcmt))
      head_i_out(5) = MIN(head.lin.ref_nr,LEN(head_linref))
      head_i_out(6) = head.lin.eps_nr
      head_i_out(7) = head.lin.w_nr
      head_i_out(8) = head.lin.status
      head_i_out(9) = MIN(head.lin.cmt_nr,LEN(head_lincmt))
      
                    
      head_f_out(1) = head.atm.teff
      head_f_out(2) = head.atm.logg
      head_f_out(3) = head.atm.x
      head_f_out(4) = head.atm.y
      head_f_out(5) = head.atm.z
      head_f_out(6) = head.atm.alpha 
      head_f_out(7) = head.lin.xi

      head_f_out(8) = head.lin.rdst             ! (resampling distance)
      head_f_out(9) = head.lin.w_start 
      head_f_out(10)= head.lin.w_end

      DO i=1,MIN(head.atm.eps_nr,NATOM)
        head_f_out(10+i) = head.atm.eps(i)
      ENDDO

      DO i=1,MIN(head.lin.eps_nr,NATOM)
        head_f_out(10+NATOM+i) = head.lin.eps(i)
      ENDDO

      head_atmdate = head.atm.date         
      head_atmuser = head.atm.user         
      head_atmref  = head.atm.ref
      head_atmcmt  = head.atm.cmt
      head_lindate = head.lin.date
      head_linuser = head.lin.user
      head_linref  = head.lin.ref
      head_lincmt  = head.lin.cmt

      wlam_cnt_out = wlam_cnt               ! Nr of flux points
      DO i = 1,wlam_cnt      
        wlam_out(i) = spc(i).w * x          ! wavelength points     [Angstroem]
        flux_out(i) = spc(i).f              ! flux       points
      END DO
      IF (wlam_out(wlam_cnt) .EQ. wlam_out(wlam_cnt-1)) THEN
        WRITE(CHN,*) 'WARNING: Identical wavelength  ',
     *            wlam_out(wlam_cnt),wlam_cnt 
      ENDIF
c
c -------------------------- compute EW's --------------------------
c trapezoidal integration
c MB 15.08.14, excitation potential print

      print*, '', flux_out(1),flux_out(wlam_cnt)
      ewidth = EW(0,wlam_out,flux_out,
     *            wlam_cnt,wlam_out(1),wlam_out(wlam_cnt))
      ewidth = 1000.*ewidth                 ! in mA
         WRITE(CHN,*), ewidth
         WRITE(CHN,*), 'excitation potential: ', EXCPOT
c --------

 8999 CONTINUE
      lf_status_out = lf_status             ! program status  bitmask
      tim = SECNDS(tim)
      IF (.NOT.BJTEST(lf_ctrl,26)) THEN
         WRITE(CHN,*)  'elapsed time: ',tim,' [s]'
         WRITE(CHN,*)  ' '
      ENDIF
      IF (log_set) CLOSE(CHN)

      RETURN
      END                                   ! MAIN
      

C -------------------------------- NEW_TAU_SCALE ------------------------------

      SUBROUTINE NEW_TAU_SCALE(r,y,dx,ndepth,n,ta_last,mode,
     *                         nscale,n_break)

      REAL*4    r(ndepth)                 ! result (new tau-scale)
      REAL*4    y(ndepth)                 ! integrand array
      REAL*4    dx(ndepth)                ! delta x-array (dx is assumed to be
                                          ! prepared by factor 0.5, LN10_DIV_2 
                                          ! or LN10_DIV_3 and LN10_3DIV8
      INTEGER*4 ndepth                    ! nr of array elements
      INTEGER*4 n                         ! last index nr
      INTEGER*4 mode                      ! 0:Simpson-/ 1,2: Trapez-Integration
      REAL*4    ta_last                   ! last tf or ti
      INTEGER*4 i,j                       ! loop index
      INTEGER*4 k,m                       ! scale index (begin and end of tau scale section)
      INTEGER*4 n_break                   ! No of breakpoints
      INTEGER*4 nscale(n_break)           ! breakpoints of tau_scale

      r(1) = y(1)*dx(1)

      IF (mode .GT. 0) GOTO 100           ! if trapez integration
C --- SIMPSON INTEGRATION METHOD   (mode = 0) 
      k    = 1
      DO j = 1,n_break
        m  = nscale(j)                    ! breakpoint of tau scale
        r(k+1) = r(k) + (y(k) + y(k+1)) * dx(k+1) 
        DO i = k+2,m
          r(i) = r(i-2)+dx(i)*(y(i-2)+4.*y(i-1)+y(i)) 
          n = i
          IF (r(i-1).GT.ta_last) THEN
            r(1) = r(2)*0.5
            RETURN
          ENDIF
        END DO
        k = m
      END DO
      r(1) = r(2)*0.5
      RETURN

C --- TRAPEZOID INTEGRATION METHOD (mode = 1,2)
100   CONTINUE
      DO i = 2,ndepth                        
         n = i
         r(n) = r(n-1) + (y(n-1) + y(n)) * dx(n)
         IF (r(n-1).GT.ta_last) THEN
           r(1) = 0.       
           RETURN
         ENDIF
      END DO
      r(1) = r(2)*0.5
      RETURN
      END         

C --- ------------------------------- STRCOMPRESS ------------------------------
C --- Elimination of blanks
      SUBROUTINE STRCOMPRESS(a,k)
      ! a: String to compress
      ! k: Nr of non-blank character
      INTEGER*4 i,n,k
      CHARACTER*(*) a
      k = 0
      DO i = 1,LEN(a)
        IF (a(i:i).NE.' ') THEN
          k = k + 1
          a(k:k) = a(i:i)
        ENDIF
      END DO
      IF (k.LT.LEN(a)) a(k+1:LEN(a)) = ' '
      RETURN
      END

C --- ------------------------------- STRTRIM ------------------------------
C --- Elimination of blanks
      SUBROUTINE STRTRIM(a,k,nwords)
      ! a: String
      ! k: Nr of non-blank character
      LOGICAL*4 alpha_found
      INTEGER*4 i,n,k,m,nwords,nwords_old
      CHARACTER*(*) a
      k = 0
      m = 0
      nwords = 0
      DO i = LEN(a),1,-1
        IF (a(i:i).LE.' ') THEN
           m = m + 1
        ELSE
           GOTO 5
        ENDIF
      END DO
 5    alpha_found = .FALSE.
      nwords_old  = 0 
      DO i = 1,LEN(a) - m
        IF ( (a(i:i).GT.' ') .OR. (alpha_found) ) THEN
          alpha_found = .TRUE.
          k = k + 1
          a(k:k) = a(i:i)
          IF (a(i:i).GT.' ') THEN
             IF (nwords.EQ.nwords_old) nwords = nwords + 1
          ELSE
             nwords_old = nwords
          ENDIF
        ENDIF 
      END DO
      IF (k.LT.LEN(a)) a(k+1:LEN(a)) = ' '
      RETURN
      END
