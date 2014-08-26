* ----------------------------------------------------------------------
* >>>
* FILENAME  : LF_DECL2.INC
* DIRECTORY :
* PURPOSE   : Declarations for READ_LINETABLE (lineform),........
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 22.07.93
* REMARK    : see also LF_DECL0.INC
* 
* <<<
* ----------------------------------------------------------------------


      INTEGER*4    LINDAT_SIZ,BLK_LEN,WREF_NR,REC_LEN,LNR,LS_MASK
      CHARACTER*7  LT_LINFOR      
      CHARACTER*10 MOMENTUM

      PARAMETER (LNR      =   200)              ! Lines per record in line list
      PARAMETER (LINDAT_SIZ = 76)               ! Size of data-segment of one line (byte-nr)
      PARAMETER (BLK_LEN    = LINDAT_SIZ * LNR) ! Blocklength (bytes)
      PARAMETER (WREF_NR    = BLK_LEN / 8)      ! Max.nr of reference wavelength in first block

      PARAMETER (REC_LEN    = BLK_LEN / 4)      ! C&VAX    Blocklength (4*byte)
C&HP  PARAMETER (REC_LEN    = BLK_LEN)          ! Blocklength 
      PARAMETER (LS_MASK    = '2C'X)            ! Line status reference (gf-& C6-Values must exist)
      PARAMETER (LT_LINFOR  = '.linfor')        ! extension of 'linfor' format
      PARAMETER (MOMENTUM='SPDFGHIJKL')

C --- Definitions for Linedata-Block
      REAL*8      wref(WREF_NR)       ! Wavelength-References (first block)
      BYTE        buf (BLK_LEN)       ! Buffer for reading each block of linedata

      REAL*8      w_ (LNR)            ! Wavelength of linecenter
      REAL*8      el_(LNR)            ! lower energylevel
      REAL*4      gf_(LNR)            ! log(gf) value (only used if xgf is NOT set)
      REAL*4      c4_(LNR)            ! log(C4) value
      REAL*4      c6_(LNR)            ! log(C6) value (only used if xc6 is NOT set)
      REAL*4      xjl_(LNR),xjh_(LNR) ! J-Nummer of angular momentum
      REAL*4      xgf_(LNR)           ! alternative log(gf)-Value
      REAL*4      xc6_(LNR)           ! alternative log(C6)-Value
      REAL*4      dw_(LNR)            ! wavelength area of significant linedepth
      REAL*4      gr_(LNR)            ! radiation damping
 
      CHARACTER*5 ip_(LNR)            ! multiplett
      CHARACTER*3 rgf_(LNR)           ! reference of gf-value
      CHARACTER*3 rc6_(LNR)           ! reference of c6-Value
      CHARACTER*3 ll_(LNR)            ! label of lower level
      CHARACTER*3 lh_(LNR)            ! label of upper level
      BYTE        iel1_(LNR)          ! atomic nr (of 1st element in diatomic molecules)
      BYTE        iel2_(LNR)          ! 2nd Element of diatomic molecule (otherwise = 0)
      BYTE        ion_(LNR)           ! Ionisationlevel
      INTEGER*4   istat_(LNR)         ! Line-Status
C ---
      INTEGER*4   maxrec              ! Max.Nr of records in linetable (<=WREF_NR)
      INTEGER*4   nrec                ! current recordnr
      INTEGER*2   mol_idx_tmp         ! buffer for molecule index
      CHARACTER*8 lbl_tmp             ! buffer for molecule label
      LOGICAL*4   nlte_line           ! =TRUE, NLTE-Lineformation for that line
      LOGICAL*4   sort_lines          ! =TRUE, if linelist must be sorted

      STRUCTURE /lintab_body_type/
        REAL*8      w                 ! Wavelength of linecenter
        BYTE        iel1              ! atomic nr (of 1st element in diatomic molecules)
        BYTE        iel2              ! 2nd Element of diatomic molecule (otherwise = 0)
        BYTE        ion               ! Ionisationlevel
        CHARACTER*5 ip                ! multiplett
        REAL*4      dw                ! wavelength area of significant linedepth
        REAL*8      el                ! lower energylevel [eV]
        REAL*4      gr                ! radiation damping
        REAL*4      gf                ! log(gf) value (only used if xgf is NOT set)
        CHARACTER*3 rgf               ! reference of gf-value
        REAL*4      xgf               ! log(gf) value (only used if xgf is NOT set)
        REAL*4      c4                ! log(C4) value
        REAL*4      c6                ! log(C6) value
        CHARACTER*3 rc6               ! reference of c6-Value
        REAL*4      xc6               ! log(C6) value
        CHARACTER*3 ll                ! label of lower level
        REAL*4      xjl               ! J-Nummer of angular momentum
        CHARACTER*3 lh                ! label of upper level
        REAL*4      xjh               ! J-Nummer of angular momentum
        INTEGER*4   istat             ! Line-Status
        INTEGER*4   mode              ! for further control (not used)
      END STRUCTURE

      STRUCTURE /lintab_head_type/
        REAL*8      wmin              ! min. Wavelength in [A]
        REAL*8      wmax              ! max. Wavelength in [A]
        BYTE        iel1              ! atomic nr (of 1st element in diatomic molecules)
        BYTE        iel2              ! 2nd Element of diatomic molecule (otherwise = 0)
        BYTE        ion               ! Ionisationlevel
        CHARACTER*5 ip                ! multiplett
        INTEGER*4   cnt               ! Nr of lines
        LOGICAL*4   unsorted          ! = 1; Linelist is unsorted (Wavelength!)
        INTEGER*4   select            ! selection criterium (0 = ALL)
        REAL*4      loggf             ! new offset or constant value
        REAL*4      logc6             ! new offset or constant value
      END STRUCTURE

      STRUCTURE /lintab_type/
        RECORD /lintab_head_type/ h            ! linetable
        RECORD /lintab_body_type/ tab(nlintab) ! linetable
      END STRUCTURE

      EQUIVALENCE (buf(1),w_)
      EQUIVALENCE (buf(1+ 8*LNR),el_)
      EQUIVALENCE (buf(1+16*LNR),xjl_)
      EQUIVALENCE (buf(1+20*LNR),xjh_)
      EQUIVALENCE (buf(1+24*LNR),dw_)
      EQUIVALENCE (buf(1+28*LNR),gf_)
      EQUIVALENCE (buf(1+32*LNR),gr_)
      EQUIVALENCE (buf(1+36*LNR),c4_)
      EQUIVALENCE (buf(1+40*LNR),c6_)
      EQUIVALENCE (buf(1+44*LNR),xgf_)
      EQUIVALENCE (buf(1+48*LNR),xc6_)
      EQUIVALENCE (buf(1+52*LNR),istat_)
      EQUIVALENCE (buf(1+56*LNR),ip_)
      EQUIVALENCE (buf(1+61*LNR),ll_)
      EQUIVALENCE (buf(1+64*LNR),lh_)
      EQUIVALENCE (buf(1+67*LNR),rgf_)
      EQUIVALENCE (buf(1+70*LNR),rc6_)
      EQUIVALENCE (buf(1+73*LNR),iel1_)
      EQUIVALENCE (buf(1+74*LNR),iel2_)
      EQUIVALENCE (buf(1+75*LNR),ion_)

