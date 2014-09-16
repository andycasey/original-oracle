* ----------------------------------------------------------------------
* >>>
* FILENAME  : LF_DECL0.INC
* DIRECTORY :
* PURPOSE   : Declaration of variables, constants and global (COMMON)
*             variable
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR,T.Gehren
* DATE      : 03.09.92
* REMARK    : units are explizit declared
*
* <<<
* ----------------------------------------------------------------------

      INCLUDE 'lf_param.inc'                    ! physical constants
      INCLUDE 'atm_def.inc'                     ! physical constants
 
C --- Functions ------------------------------------------------------------
      REAL*4      QPOL                  ! Interpolation function
      INTEGER*4   IGETENV               ! get environment variables

      REAL*8 small
      PARAMETER  (small = 1.D-36) ! Minimum wavelength-spacing [cm] 

C --- Headerinfo -----------------------------------------------------------
   
      STRUCTURE /atm_head_type/
        CHARACTER*30  date              ! atmosphere generation date
        CHARACTER*20  user              ! atmosphere calc. by user ...
        REAL*4        teff              ! effective temp.
        REAL*4        logg              ! log(g) (gravity)
        REAL*4        x                 ! [H]
        REAL*4        y                 ! [He/H]
        REAL*4        z                 ! [Fe/H]
        REAL*4        alpha             ! l/Hp   Mixinglength
        INTEGER*2     status            ! Atm Status
        INTEGER*2     eps_nr            ! Nr of atm element abundances
        REAL*4        eps(NATOM)        ! atm element abundances
        INTEGER*2     ref_nr            ! Nr of Ref.-character (atm)
        CHARACTER*1000 ref              ! atm references   
        INTEGER*2      cmt_nr           ! Nr of atm comment character
        CHARACTER*1000 cmt              ! atm comment
      END STRUCTURE

      STRUCTURE /lin_head_type/
        CHARACTER*30  date              ! spectrum generation date (lin)
        CHARACTER*20  user              ! spectrum calc. by user ..(lin)
        REAL*4        xi                ! Microturbulence (km/s)   (lin)
        REAL*4        rdst              ! resampling distance      
        REAL*4        w_start           ! start wavelength
        REAL*4        w_end             ! end   wavelength
        INTEGER*2     w_nr              ! nr of wavelength     
        INTEGER*2     status            ! lin status
        INTEGER*2     eps_nr            ! Nr of element abundance differences
        REAL*4        eps(NATOM)        ! element abundance differences
        INTEGER*2     ref_nr            ! Nr of Ref.-character     (lin)
        CHARACTER*1000 ref              ! lin references           (lin)   
        INTEGER*2      cmt_nr           ! Nr of comment character  (lin)
        CHARACTER*1000 cmt              ! lin comment
      END STRUCTURE    

      STRUCTURE /head_type/             ! type declaration full header
        CHARACTER*20 id                 ! spectrum-id
        RECORD /atm_head_type/ atm        
        RECORD /lin_head_type/ lin
      END STRUCTURE 

      RECORD /head_type/ head           ! full header (ATM-& LIN-information)

      REAL*4    ximicro                   ! microturbulence             [cm/s]


C --- additional

      INTEGER*4   i,j,k,m,n,nn            ! loop-indices
      INTEGER*4   atomic_nr               ! = idx_el(k)
      INTEGER*8   wlam_cnt                ! number of currently calc. flux points
      INTEGER*8   cont_nr                 ! number of currently calc. continuumflux points
      REAL*4      tn                      ! = T(n)
      REAL*4      divtn                   ! = 1/T(n)
      REAL*4      f0(NDMAX)               ! electrical normal fieldstrength
      REAL*4      fh(NDMAX)               ! precalculated factor for H line opacity
      REAL*4      vdop(NDMAX,NATOM+NMOLMAX) ! doppler velocity of particle
      REAL*4      gamma(NDMAX,NLINE)      ! damping factor
      REAL*4      Blam   (NDMAX)          ! Planckfunktion:     1 / (exp(E/kT(n)) - 1)                 
      REAL*4      Slam   (NDMAX)          ! Sourcefunktion    
      REAL*4      Jlam   (NDMAX)          ! depthdependent mean intensity    
      REAL*4      tau_lam(NDMAX)          ! tau(lambda)   of line
      REAL*4      eta_lam(NDMAX)          ! kappa(lambda)/kappa_ref
      REAL*4      eta_cont_lam(NDMAX)     ! kappa_cont(lambda)/kappa_ref
      REAL*4      etasig_lam(NDMAX)       ! sigma(lambda)/kappa_ref
      REAL*4      kappa_line(NDMAX,NLINE) ! line opacity at n-th depthpoint (* 1/dopplerwidth)
                                          ! = (FAC_B * divtn) * pion_u(n,ion(lc),k)
                                          !   * exp(-(line(lc).elow(lc)+
                                          !          HC_DIV_K/line(lc).w)*divtn)
                                          !   * div_dlamdop(n,lc) * gf(lc)
      REAL*8      kappa_lam               ! sum over all kappa(lambda)*profile(lambda)
      REAL*8      w                       ! temporary wavelength
      REAL*8      dw                      ! current stepwidth
      REAL*8      dlam                    ! delta wavelength: ABS(w-line(lc).w)
      REAL*8      wlast_cont              ! last wavelength point of continuum
      CHARACTER*60 lf_output_filename     ! Outputfilename of LINEFORM defined in ATMOSIN

C ---- see ELEMDATA.INC ------------------------------------------------------


      INTEGER*4   nelpos(NATOM)           ! index-map
      INTEGER*4   nion  (NATOM)           ! Nr of Ionisationsteps
      REAL*4      epss  (NATOM)           ! solar abundance
      REAL*4      amass (NATOM)           ! atomic mass        (see 'mass')
      REAL*4      eps_dev(NATOM)          ! deviations from solar abundances
      REAL*4      eps   (NATOM)           ! rel. abundance (deviations from solar
      CHARACTER*4 cel   (NATOM)           ! Element-Label
                                          !                 abundances included)
 
      REAL*4      mass  (NATOM+NMOLMAX)   ! amass (1,NATOM) or molmass (1+NATOM,NATOM+NMOL)
      CHARACTER*4 lbl   (NATOM+NMOLMAX)   ! Label of Atom or Molecule (s. EQUIVALENCE)
 
C ---- LINE DATA -

      INTEGER*4     lc                    ! Line-Counter: Index of actual Line
      INTEGER*4     linemax               ! Nr. of Lines taken into accout
      CHARACTER*150 linedata_file         ! for filename of linedatafile in use

      STRUCTURE /line_type/
        LOGICAL*1   use                   ! Ionisation-step (0,1 or 2)
        REAL*8      w                     ! center wavelength of line
        INTEGER*4   id                    ! line-id: 1-NATOM - (atomic line)
                                          !          NATOM < - (molecular line)
                                          !      -> line_id = molecule-idx + NATOM
        INTEGER*2   ion                   ! Ionisation-step (0,1 or 2)
        REAL*4      elow                  ! lower exitation energy
        REAL*4      gflog                 ! log(gf)
        REAL*4      c4log                 ! log(C4)  (Unsoeld-Def.: dv = C4/r^4)
        REAL*4      c6log                 ! log(C6)  (Unsoeld-Def.: dv = C6/r^6)
        REAL*4      grad                  ! radiation damping
        REAL*8      width                 ! suggested maximum contribution area
        REAL*4      xjl                   ! J_lower
        REAL*4      xju                   ! J_upper
        INTEGER*2   mode                  ! for control (not used)
        INTEGER*4   rbb_idx               ! transition index for POP (NLTE!)
      END STRUCTURE

      RECORD /line_type/ line(NLINE)      ! linetable

      REAL*8      wmin,wmax               ! min. and max. wavelength of spectrum      
      REAL*8      wline_center,wline_code ! for single lineformation see BIT 23
      INTEGER*4   lin_rec_nr,lin_rec_pos  ! lineposition in linetable

C ---- Hydro  ----

      LOGICAL*4   hyd                     ! TRUE, if Balmerlines are concerned
      INTEGER*4   nhw_nr(NBALMER)         ! Nr of wavelengthpoints of Balmerline k
      INTEGER*4   balmer_min,balmer_max   ! min. and max. index of concerned balmerline (H-Alpha -> k=1)
      REAL*4      e_vcs(NHE,NBALMER)      ! electron density grid points
      REAL*4      t_vcs(NHT,NBALMER)      ! temperatur grid points
      REAL*4      w_vcs(NHW,NBALMER)      ! wavelength grid points      [cm]
      REAL*4      phi_vcs(NHW,NDMAX,NBALMER)! buffer for VCS-Profiles
      REAL*4      broad                   ! additional resonance broadening 
      REAL*8      wbalmer  (NBALMER_ALL)  ! Balmerlinepositions         [cm]
      REAL*8      balmer_gf(NBALMER_ALL)  ! gf-Values for Balmerlines

C ---- HydroNEw  ----
      INTEGER*4 nhyd                      ! Balmer (1) or Paschen (2) lines are concerned
      INTEGER*4 minh,maxh                 ! Min. and max. index of hydrogen lines
      INTEGER*4 nlow(NNHMAX),nup(NNHMAX),nhhnr(NNHMAX)
      REAL*4 nh_vcs(NNHH)                   ! Temperature grid points
      REAL*4 ne_vcs(NNHE)                   ! Electron density grid points
      REAL*4 nt_vcs(NNHT)                   ! Temperature grid points
      REAL*4 nw_vcs(NNHW)                   ! Wavelength grid points [cm]
      REAL*4 nphi_vcs(NNHW,NDMAX,NNHMAX)     ! Buffer for VCS-Profiles
      REAL*8 whyd(NNHMAX,2)                ! Hydrogen line wavelengths [cm]
      REAL*8 gfhyd(NNHMAX,2)               ! gf values for hydrogen lines
      REAL*4 gfine(7,4)                   ! Balmer line fine structure f-values
      REAL*8 wfine(7,4)                   ! Balmer line fine structure wavelengths [cm]

C ---- Energy ----

      REAL*4      cos_theta
      REAL*8      sw_crit                 ! crit. of max. deviation (see calc.of lambda-scale)
      REAL*8      dwmin                   ! min. stepwidth of lambda-scale
      REAL*8      dwmax                   ! max. stepwidth of lambda-scale

      STRUCTURE /spectrum/
       REAL*8      w                      ! wavelength
       REAL*4      f                      ! flux- or irradiance
      END STRUCTURE
      RECORD /spectrum/ spc(NWMAX)        ! Spectrum,  common linetable
      RECORD /spectrum/ cspc(NCMAX)       ! Continuum, common linetable

      REAL*8      wa(6),wi(6),wf(6)       ! integration weights
      REAL*4      ta(6),ti(6),tf(6)
C      	          wi,   ti                ! (gauss-) integration weights for irradiance
C      	          wf,   tf                ! integration weights for flux calculation
                  
C ---- Control ----
 
      INTEGER*4        lf_ctrl            ! Control-Bitmask ("... pgm has to do ...")
      INTEGER*4        lf_status          ! Status -Bitmask ("... pgm has already done ...")

      EQUIVALENCE      (lbl(1),cel),    (lbl(1 + NATOM),cmol)
      EQUIVALENCE      (mass(1),amass), (mass(1+NATOM) ,molmass)


      COMMON/ATMOS/    f0, fh, vdop,      
     *                 gamma, kappa_line, lf_output_filename

      COMMON/ELEMENT/  lbl, nelpos, epss, mass, eps, 
     *                 mol_idx, eps_dev, nion

      COMMON/LINES/    lc,linemax, wmin, wmax,line,linedata_file,
     *                 wline_center,wline_code,lin_rec_nr,lin_rec_pos

      COMMON/HYDRO/    hyd, balmer_min, balmer_max, nhw_nr,
     *                 w_vcs, phi_vcs, e_vcs, t_vcs

      COMMON/HYDRONEW/ nhyd,minh,maxh,nw_vcs,nh_vcs,ne_vcs,nt_vcs,
     *                 nphi_vcs

      COMMON/ENERGY/   cos_theta, sw_crit, dwmin, dwmax, wlam_cnt, 
     *                 cont_nr, spc, cspc, wa, ta

      COMMON/HEADER/   head                  ! see structure def. above
                  
      COMMON/CONTROL/  lf_ctrl, lf_status


C ---- DATA - SECTION                 

C --- Integration weights (see head-description of LINEFORM.FOR)
      ! -- irradiance -- (Cayrel & Traving)
C      DATA    wi/ .3699D0,.2623D0, .1360D0, .0965D0, .1155D0, .0198D0/
C      DATA    ti/ .1808,.7524,1.1808,1.7524,2.5858,5.4142/

C --- Laguerre Integration weights (see head-description of LINEFORM.FOR)
      ! -- irradiance -- (Abramowiz,Stegun)
      DATA    wi/ .4589647D0,.4170008D0,0.1133734D0,
     *            .0103991D0,.2610172D-3,.8985479D-6/
      DATA    ti/ .222847D0 ,1.188932D0, 2.992736D0, 
     *           5.775144D0 ,9.837467D0,15.982874D0 /

      ! -- flux       -- calculated by J.Reetz (1992)
      DATA    wf/ 0.543818D0  , 0.372530D0  , 0.777969D-01,
     *            0.573486D-02, 0.120155D-03, 0.353505D-06/
      DATA    tf/ 0.161342,0.94604, 2.524137, 5.066564,
     *            8.880943,10.4758785/

C --- wavelengths of Balmerlines          [cm]
      DATA wbalmer/ 6.56282D-5, 4.86133D-5, 4.34047D-5, 4.10174D-5
     *             ,3.97007D-5, 3.88905D-5, 3.83539D-5, 3.79790D-5
     *     ,3770.63D-8,3750.15D-8,3734.37D-8,3721.94D-8,3711.97D-8
     *     ,3703.86D-8,3697.15D-8,3691.56D-8/

C --- gf-values of Balmerlines
      DATA balmer_gf/ 5.12625, .95455, .35736, .17673,
     *                 .10165, .06428, .04343, .03080,8*0./

C-----Wavelengths of hydrogen lines [cm] NEW
      DATA whyd/
     &  6562.82D-8, 4861.33D-8, 4340.47D-8, 4101.74D-8, 3970.07D-8,
     &  3889.05D-8, 3835.39D-8, 3797.90D-8, 3770.63D-8, 3750.15D-8,
     &  3734.37D-8, 3721.94D-8, 3711.97D-8, 3703.86D-8, 3697.15D-8,
     &  3691.56D-8, 3686.83D-8, 3682.81D-8,
     & 18751.00D-8,12818.10D-8,10938.10D-8,10049.40D-8, 9545.98D-8,
     &  9229.02D-8, 9014.91D-8, 8862.79D-8, 8750.47D-8, 8665.02D-8,
     &  8598.39D-8, 8545.39D-8, 8502.49D-8, 8467.26D-8, 8437.96D-8,
     &  8413.32D-8, 8392.40D-8,     .0/
C-----gf-values of hydrogen lines
      DATA gfhyd/
     &  5.12625, .95455, .35736, .17673, .10165, .06428, .04343,
     &   .03080, .02268, .01721, .01388, .01061, .00856, .00701,
     &   .00582, .00488, .00413, .00353,
     & 15.15780,2.71080,1.00512, .49824, .28872, .18414, .12564,
     &   .08993, .06680, .05110, .04003, .03197, .02597, .02138,
     &   .01785, .01505, .01281, .0/
C-----Fine structure wavelengths for the first four Balmer lines
      DATA wfine/
     &  6562.772,6562.725,6562.752,6562.910,6562.710,6562.868,6562.852,
     &  4861.298,4861.287,4861.289,4861.375,4861.279,4861.365,4861.362,
     &  4340.438,4340.434,4340.431,4340.500,4340.427,4340.496,4340.494,
     &  4101.711,4101.709,4101.705,4101.766,4101.703,4101.764,4101.763/
C-----Fine structure oscillator strengths for the first four Balmer lines
      DATA gfine/      ! gf-values
     &  0.289936,0.579864,0.027176,0.054360,1.391600,0.278320,2.504880,
     &  0.068536,0.137064,0.006088,0.012176,0.243600,0.048720,0.438480,
     &  0.027952,0.055904,0.002424,0.004856,0.088736,0.017752,0.159728,
     &  0.014416,0.028840,0.001232,0.002472,0.043264,0.008656,0.077864/

