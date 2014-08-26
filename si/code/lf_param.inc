* ----------------------------------------------------------------------
* >>>
* FILENAME  : LF_PARAM.INC
* DIRECTORY :
* PURPOSE   : Declaration of constants
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 28.12.95
* REMARK    : units are explizit declared
*
* <<<
* ----------------------------------------------------------------------

C---- PARAMETER settings refer to maximum number of array elements 

      INTEGER*4    NELMAX,NMOLMAX,NATOM,NR_ION,NDMAX,
     *             NHT,NHE,NHW,MAX_ETA_CONT,NBALMER,NBALMER_ALL,
     *             NWSAVMAX,NCMAX,CHN,NLINTAB,NMAXTERMS,
     *             STD_OUTPUT_CHN,NHEAD_I,NHEAD_F,NHEAD_C,
     *             NR_SPLINE,MAXBREAK,N_ATM_HOLD,N_ATMGRID_HOLD
      INTEGER*8    NWMAX,NLINE
      REAL*8       STD_WIDTH,STD_SW_CRIT,STD_DWMIN,STD_DWMAX,
     *             W_EPS,CRANGE,XMAX_LINE_RANGE,DIV_WLAMREF,XKAPMIN,
     *             DDEPS
      REAL*4       EEEPS
      CHARACTER*80 BALMER_FILE,LOG_MSG
      CHARACTER*20 RBF_FILE
      CHARACTER    dir_t*28
      CHARACTER*4  ext_dep
      CHARACTER*4  ext_depidx
      PARAMETER       (ext_dep   =  '.dep') 
      PARAMETER       (ext_depidx=  '.idx') 

      PARAMETER       (NELMAX   =    92)  ! Elements used in atmosphere calculation
      PARAMETER       (NMOLMAX  =    17)  ! Maximum-nr of Molecules
      PARAMETER       (NATOM    =    92)  ! Maximum-nr of Elements (s. ELEMDATA.INC)
      PARAMETER       (NR_ION   =     3)  ! Maximum-nr of Ionization stages in atm.structure
      PARAMETER       (NMAXTERMS=   700)  ! Maximum-nr of TERMS in each lookup table (only NLTE lineformation)
      PARAMETER       (NDMAX    =    80)  ! max.nr of atmospheric depth points
      PARAMETER       (MAXBREAK  =    2)  ! Nr of Breakpoints in tauref-scale
      PARAMETER       (DIV_WLAMREF = 1.D0/5000.D-8)  ! max. range to next continuum flux point [cm/5000A]
      PARAMETER       (CRANGE   = 5.D-8*DIV_WLAMREF) ! max. range to next continuum flux point [cm/5000A]
      PARAMETER       (NWMAX    =  1000000)  ! Frequency points
      PARAMETER       (NCMAX    =   500)  ! Frequency points of continuum
      PARAMETER       (NHMAX    =    18)  ! Max. number of VCS lines in series

      PARAMETER       (STD_WIDTH= 0.5D-8*DIV_WLAMREF)  ! Standard wavelength range [A/5000A] for line contributions
      PARAMETER       (NLINE    =  1000000)  ! Number of lines (see structure LINE)
      PARAMETER       (NHT      =     4)  ! Temperature grid points in VCS table
      PARAMETER       (NHE      =    22)  ! Density grid points in VCS table
      PARAMETER       (NHW      =    35)  ! Nr of wavelength points in VCS table
      PARAMETER       (NBALMER  =     8)  ! Nr of VCS-Balmerprofiles taken into account
      PARAMETER       (NBALMER_ALL=  16)  ! Nr of VCS+ESW-Balmerprofiles taken into account
      PARAMETER       (NNHT      =     4)  ! Temperature grid points in VCS table
      PARAMETER       (NNHE      =    13)  ! Electron density grid points in VCS table
      PARAMETER       (NNHH      =    13)  ! Max. number of hydrogen density grid points in VCS table
      PARAMETER       (NNHW      =    44)  ! Nr of wavelength points in VCS table
      PARAMETER       (NNHMAX    =    18)  ! Max. number of VCS lines in series

      PARAMETER       (NSTART_ESW=    9)  ! StartNr (m) of ESW-Balmerprofiles taken into account
      PARAMETER       (NLINTAB  =     1)  ! Dummy (see LF_DECL2.INC -> LINTAB_TYPE)
      PARAMETER       (XMAX_LINE_RANGE=25.D-8*DIV_WLAMREF) ! [A] wline +- XMAX_LINE_RANGE
C      PARAMETER       (XMAX_LINE_RANGE=10.D-8*DIV_WLAMREF)! [A] wline +- XMAX_LINE_RANGE
      PARAMETER       (MAX_ETA_CONT=NCMAX)  ! max number of continuum grid points
C                                              = lineforming Region in case of
C                                                single lineformation
C
      PARAMETER       (NR_SPLINE=   501)        ! Nr of Spline Interpolation points
      PARAMETER       (STD_SW_CRIT=0.001D0)     ! Standardkriterium for dlamda-dubbling [cm]
      PARAMETER       (STD_DWMIN=  5.D-11*DIV_WLAMREF) ! Minimum stepwidth of wavelengthscale  [cm]
      PARAMETER       (STD_DWMAX=100.D-8*DIV_WLAMREF) ! Maximum stepwidth of wavelengthscale  [cm]
      PARAMETER       (W_EPS    =  3.D-11)      ! Minimum resolution [cm] at 500 nm
      PARAMETER       (EL_EPS   =  5.E-03)      ! Max. deviation from lower energy level in [eV]
      PARAMETER       (XKAPMIN1 =  1.D-03)      ! minimum fraction of kappa_line/kappa_cont
      PARAMETER       (XKAPMIN  =  1.D-04)      ! minimum fraction of kappa_line/kappa_cont*Phi(delta lambda)
      PARAMETER       (EEEPS =  1.E-30)   ! numerical accuracy
      PARAMETER       (STD_OUTPUT_CHN= 6) ! Standard output channel (terminal)
      PARAMETER       (CHN           = 6) ! channel for log-messages
      PARAMETER       (LOG_MSG='lineform.log') ! Logfile (if CHN.NE.STD_OUTPUT_CHN)

      PARAMETER       (NHEAD_I =      9)   !  Nr of integer*2 values of header
      PARAMETER       (NHEAD_F =10+2*NATOM)! Nr of real*4 values of header

      PARAMETER       (N_ATM_HOLD     =  2)! MaxNr of stacked Atmospheres
      PARAMETER       (N_ATMGRID_HOLD =  2)! MaxNr of stacked Atm.gridpoints

      PARAMETER (BALMER_FILE='balmer_vcs.dat')              ! VCS-Profiles
      PARAMETER (RBF_FILE   ='uphoto.asc')                  ! VCS-Profiles
      PARAMETER (dir_t='/usr/users/arch3/atm_models/')      ! directory with grid
 
      INTEGER*4 env_len                         ! length of environment var.
      CHARACTER*255 env                         ! environment var. (s.IGETENV)
 
C ---- see MOL_DATA.INC ------------------------------------------------------
      BYTE        mol_idx(NATOM,NATOM)    ! 0=if not considered, otherwise the
                                          ! molecule-index ordered in MOL_DATA.INC
      CHARACTER*4 cmol   (NMOLMAX)        ! Molecule-Label
      REAL*4      molmass(NMOLMAX)        ! molecular mass     (see 'mass')
      INTEGER*2   mel1(NMOLMAX),mel2(NMOLMAX)
C ----------------------------------------------------------------------------

