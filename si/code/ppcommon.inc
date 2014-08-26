* ----------------------------------------------------------------------
* >>>
* FILENAME  : ppcommon.inc
* DIRECTORY :
* PURPOSE   : Declarations for lineform (partpress section)
* MODULE    : FORTRAN - INCLUDE FILE
* AUTOR     : JKR
* DATE      : 28.12.95
* <<<
* ----------------------------------------------------------------------

      REAL*4 pion(NDMAX,NR_ION,NATOM)	   ! ions
      REAL*4 uion(NDMAX,NR_ION,NATOM)	   ! partition function of ions

*     REAL*4 eps(NATOM)	                   ! eps(n) = A(n)/(A(1) + A(2) + ... + A(n_eps))

      COMMON /pp/ pion

      COMMON /WORK/sion(2,NATOM),s(NATOM),zs(NATOM),ZZS(NATOM),
     *            ZZZS(NATOM)


