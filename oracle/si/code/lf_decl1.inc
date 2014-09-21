* ----------------------------------------------------------------------
* >>>
* FILENAME  : LF_DECL1.INC
* DIRECTORY :
* PURPOSE   : Declaration of variables, constants and global (COMMON)
*             variables concerning NLTE-Lineformation
*
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 22.07.93
* REMARK    : 
*
* <<<
* ----------------------------------------------------------------------

C---- PARAMETER
      INTEGER*4     NLTE_UT,MAXLEVEL,MAXRBB,NCHRLAB
      REAL*8        MIN_EV
      CHARACTER*120 pop_file
      PARAMETER (XDLAM    =  0.1)    ! [A], Min Difference between w(popfile/lookup)-w(linefile)
      PARAMETER (MAXLEVEL = 300)                  ! Max Nr of level
      PARAMETER (MAXRBB   =3000)                  ! Max Nr of Transitions
      PARAMETER (MIN_EV   = 0.0001D0)             ! Min diff. between Energylevels
      PARAMETER (NCHRLAB  =   6)                  ! Nr of Designation character
      PARAMETER (NLTE_UT  =  31)                  ! Input channel of POP-File

C---- TYPE DECLARATION

      STRUCTURE /rbb_typ/
        INTEGER*4     ll                      ! index of lower level
        INTEGER*4     ul                      ! index of upper level
        REAL*8        w                       ! [cm]
        LOGICAL*4     accept                  ! TRUE, if transition is taken
        REAL*4        q(NDMAX)                ! b_lower/b_upper
      ENDSTRUCTURE
        
      STRUCTURE /level_typ/
        CHARACTER*(NCHRLAB) desig             ! Designation / Label
        INTEGER*4     ion                     ! Ion Nr
        REAL*8        ev                      ! Energy [eV] 
        REAL*4        b(NDMAX)                ! departures 
      ENDSTRUCTURE

      STRUCTURE /pop_typ/
        INTEGER*4   atom                      ! Atomic Nr of Element
        INTEGER*4   ndepth                    ! Nr of atmospheric depthpoints
        INTEGER*4   level_nr                  ! Nr of Level         
        INTEGER*4   rbb_nr                    ! Nr of Transitions   
        RECORD /LEVEL_TYP/ level(MAXLEVEL)
        RECORD /RBB_TYP/   rbb(MAXRBB)
      ENDSTRUCTURE
      
C---- VARIABLE DECLARATION

      CHARACTER*(NCHRLAB) desig               ! Buffer for Label/Designation
      INTEGER*4  idx_arr(MAXLEVEL)
      INTEGER*4  rbb_idx
      REAL*8     tmp_arr(MAXLEVEL)
      REAL*8     em_s                         ! 1. - exp(-x)
      REAL*8     emi_s                        ! exp(-x)
      RECORD /POP_TYP/ pop

      COMMON /NLTE/ pop

      INTEGER*2 il_lup,iu_lup 
      INTEGER*4 n_nlte_lines,n_terms,idx_list(MAXRBB),
     *            ndepth_ali,nlev_ali,natom_ali,niter_ali      
      REAL*4    mult_code,w_lup,el_lup,xjl_lup,xju_lup,
     *            teff_ali,logg_ali,feh_ali,xi_ali,
     *            epsdev_ali(NATOM), hbbfac_ali,hbffac_ali,
     *            escale_ali,rbfscl_ali,
     *            cbf_scl_ali,gmax_ali,gmin_ali,xdum(NDMAX),
     *            logtau5000_ali(NDMAX)
      REAL*8    wmin_lup,wmax_lup
      CHARACTER*6 term_lup(NMAXTERMS),term_ali(NMAXTERMS),
     *            terms(1000)
      CHARACTER*200 atm_ali

