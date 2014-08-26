





C --- ------------------------------ READ_VCS_HYD ------------------------------
      SUBROUTINE READ_VCS_HYD
C --- Input of tabulated VCS Balmer profiles based on the UNIFIED THEORY as 
C --- developed by Vidal, Cooper and Smith (1970, JQSRT 10,1011), 
C --- This subroutine determines the STARK-Function PHI(NUE), by linear 
C --- interpolation of tables calculated by T.Schoening and K.Butler (1989).
C --- Subroutine adapted from Klaus Fuhrmann

      INCLUDE 'physcnst.inc'           ! Physical constants
      INCLUDE 'lf_decl0.inc'           ! common declarations

      INTEGER*4 itmp                   ! temporary integer
      INTEGER*4 nhe_nr(NBALMER)        ! Nr of electron-density points of line k
      INTEGER*4 nht_nr(NBALMER)        ! Nr of temperatur points of line k
      REAL*4    nelog                  ! log(e-density)
      REAL*4    tlog                   ! log(T(n))
      REAL*4    grid_vcs(NHW,NHE,NHT,NBALMER) ! buffer for VCS-profiles of all NBALMER lines
      REAL*4    grid_vcs_ijk(NHT)      ! temporary buffer for VCS-profile interpolation
      REAL*4    grid_vcs_ik (NHE)      ! temporary buffer for VCS-profile interpolation
      REAL*4    e_vcs_k     (NHE)      ! temporary buffer for e-density gridpoints of one line
      REAL*4    t_vcs_k     (NHT)      ! temporary buffer for t-gridpoints of one line

      COMMON/READ_VCS_HYD_WORK/itmp,nhe_nr,nht_nr,nelog,tlog,grid_vcs,
     *          grid_vcs_ijk,grid_vcs_ik,e_vcs_k,t_vcs_k
      env_len = IGETENV('SIU_MAIN',env)
      OPEN (UNIT = 28, FILE = env(1:env_len)//'/data/'//BALMER_FILE, 
     *      STATUS = 'UNKNOWN', FORM = 'UNFORMATTED',READONLY)
C&VAX OPEN (UNIT = 28, FILE = BALMER_FILE, STATUS = 'UNKNOWN',
C&VAX*      FORM = 'UNFORMATTED',READONLY,SHARED)
      READ (28) (nhw_nr(k),(w_vcs(j,k),j=1,nhw_nr(k)),
     *           nht_nr(k),(t_vcs(j,k),j=1,nht_nr(k)),
     *           nhe_nr(k),(e_vcs(j,k),j=1,nhe_nr(k)),
     *           (((grid_vcs(i,j,m,k),i=1,nhw_nr(k)),
     *             m=1,nht_nr(k)),j=1,nhe_nr(k)),k=1,NBALMER)
      CLOSE(28)

C --- Get relevant Balmer lines for quadratic interpolation
      balmer_max = 0      
      balmer_min = 1

      DO k = 1,NBALMER
        IF ( wmin .LE. wbalmer(k) + 1.e-6 )  balmer_max = k
      END DO

      DO k = NBALMER,1,-1
        IF ( wmax .GE. wbalmer(k) - 1.e-6 )  balmer_min = k
      END DO

      DO k = balmer_min,balmer_max           ! begin of balmerlines loop
        DO i = 1,nht_nr(k)
          t_vcs_k(i) = t_vcs(i,k)
        END DO

        DO i = 1,nhe_nr(k)
          e_vcs_k(i) = e_vcs(i,k)
        END DO

        DO n = 1, atm.h.ndepth                      ! begin of stratification loop
          tlog   = LOG10(atm.t(n))
          nelog  = LOG10(atm.pe(n)) - tlog - LOG_K

C --- Check temperature and density range

          IF (k .EQ. balmer_min) THEN
            IF ( (tlog .LT. t_vcs(1,k)) .OR.
     *           (tlog .GT. t_vcs(nht_nr(k),k)) ) THEN
              WRITE(CHN,*) 'temperature - interpolation out of bounds'
              lf_status = JIBSET(lf_status,0)         ! Fatal Error
              lf_status = JIBSET(lf_status,4)         ! Fatal Error
              RETURN
            ENDIF
            IF ((nelog .LT. e_vcs(1,k)) .OR.
     *                 (nelog .GT. e_vcs(nhe_nr(k),k))) THEN
              WRITE(6,*) 'electron density - interpolation out',
     *                   ' of bounds --> extrapolation !'
C              lf_status = JIBSET(lf_status,0)         ! Fatal Error
C              lf_status = JIBSET(lf_status,5)         ! Fatal Error
C              RETURN
            ENDIF
          ENDIF
                                             
C ---     Quadratic interpolation of PHI(NUE) as a function of (T,NE)
          DO i = 1, nhw_nr(k)
            DO j = 1,nhe_nr(k)
              DO m = 1,nht_nr(k)
                grid_vcs_ijk(m) = grid_vcs(i,j,m,k)
              END DO
              grid_vcs_ik(j) =
     *            QPOL(tlog,t_vcs_k,grid_vcs_ijk,nht_nr(k),itmp)
            END DO
            phi_vcs(i,n,k) =
     *          QPOL(nelog,e_vcs_k,grid_vcs_ik,nhe_nr(k),itmp)
          END DO
        END DO                               ! end of stratification loop
      ENDDO                                  ! end of balmerlines loop
      RETURN
      END
