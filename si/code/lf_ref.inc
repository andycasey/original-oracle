* ----------------------------------------------------------------------
* >>>
* FILENAME  : LF_REF.INC
* DIRECTORY :
* PURPOSE   : Declaration of variables, constants concerning
*             the reference spectrum (see LINFIT.FOR,LINEFORM.FOR)
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 24.06.95
* REMARK    : 
* CONTEXT   : LINEFORM
*
* <<<
* ----------------------------------------------------------------------

      PARAMETER (MAX_REF_NR = 32768)    ! max.Nr of reference data points

      INTEGER*4  nref                   ! nr of ref points
      STRUCTURE /ref_type/
        REAL*8     w                    ! wavelength
        REAL*4     f                    ! flux, intensity
        REAL*4     s                    ! standard deviation
        REAL*4     g                    ! weights
      END STRUCTURE
