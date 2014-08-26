C --- ---------------------------- read_lines ------------------------------
C

      SUBROUTINE read_lines(sort_lines, EXCPOT)

      INCLUDE 'physcnst.inc'                    ! Physical constants
      INCLUDE 'lf_decl0.inc'                    ! general declarations, Common-Blocks
      INCLUDE 'lf_decl1.inc'                    ! 'NLTE declarations'
      INCLUDE 'lf_decl2.inc'                    ! read_lines declarations
      RECORD /line_type/ line_tmp               ! linetable sort

      INTEGER*4 icode0,icode1,icode2            ! temp. var.
      INTEGER*4 mask
      INTEGER*4 lintab_fmt                      !
      REAL*8    xw                              ! buffer
      REAL*4    EXCPOT                          ! excitation potential

      COMMON /read_lines_WORK/ wref, buf

      IF (sort_lines) GOTO 2005                 ! only sorting lines
      lintab_fmt = 0                            ! standard AFR_LIN format
      sort_lines = .FALSE.
      mask       = LS_MASK
      IF (INDEX(linedata_file,LT_LINFOR).GT.1) lintab_fmt = 1

      IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*)
     *                 'Linefile= ',linedata_file

      wmin = wmin * 1.e8               ! [cm] -> [A]
      wmax = wmax * 1.e8               ! [cm] -> [A]

c MB 15.08.14
      EXCPOT = -1.

      IF (lintab_fmt.EQ.0) THEN        ! linetable format    - binary
c       PRINT*, 'binary masterline table'
        IF (BJTEST(lf_ctrl,28)) THEN
          IF (.NOT.BJTEST(lf_ctrl,26)) 
     *           write(CHN,*) 'all lines are considered'
          mask = 'C'X                  ! (gf-& C6-Values must exist)
        ENDIF
       IF (BJTEST(lf_ctrl,29).AND..NOT.BJTEST(lf_ctrl,26)) 
     *                  write(CHN,*) 'predicted lines ignored'
       IF (BJTEST(lf_ctrl,2).AND..NOT.BJTEST(lf_ctrl,26)) 
     *                  write(CHN,*) 'only selected lines used'

       OPEN(UNIT=25,FILE=linedata_file,ACCESS='DIRECT',RECL=REC_LEN,  ! C&HP
     *      STATUS='OLD',READONLY)

       maxrec = WREF_NR
       lc   = 0

       IF (BJTEST(lf_ctrl,25)) THEN            ! exact line formation
        READ(25,REC=lin_rec_nr) buf           
        n  = lin_rec_pos

        wmin = w_(n)*(1. - XMAX_LINE_RANGE)
        wmax = w_(n)*(1. + XMAX_LINE_RANGE)
        mol_idx_tmp = 0
        PRINT*, 'exact line formation for elements: ', iel2_(n)
        IF (iel2_(n) .LE. 0) THEN
          IF (nelpos(iel1_(n)).LE.0) THEN
                 PRINT*, 'nelpos(iel1_(n)).LE.0'
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *                lbl(iel1_(n)),'  Line: ',w_(n),' ignored'
            GOTO 999
          ELSEIF (iel1_(n).EQ.1) THEN     ! H-Lines are especially considered
            GOTO 999
          ENDIF
        ELSE
          mol_idx_tmp = mol_idx(iel1_(n), iel2_(n))
          IF (mol_idx_tmp .LE. 0) THEN
            lbl_tmp = cel(iel1_(n)) // cel(iel2_(n))
            CALL STRCOMPRESS(lbl_tmp,k)
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *                 lbl_tmp(1:k),'  Line: ',w_(n),' ignored'
            GOTO 999
          ENDIF
        ENDIF

        lc = 1
        line(lc).w = w_(n) * 1.D-8                         ! [cm]
        IF (lc .GT.1) THEN
          IF (line(lc).w .LT. line(lc-1).w) sort_lines = .TRUE.
        ENDIF
        IF (iel2_(n) .LE. 0) THEN
          line(lc).id = iel1_(n)                           ! Atom
        ELSE
          line(lc).id = NATOM + mol_idx(iel1_(n), iel2_(n))! Molecule
        ENDIF
        line(lc).ion   = ion_(n)
        line(lc).elow  = el_(n) * E_SI * 1.e7              ! [erg]
        line(lc).gflog = gf_(n)                  
        line(lc).c4log = c4_(n)    
        line(lc).c6log = c6_(n)
        line(lc).xjl = xjl_(n)
        line(lc).xju = xjh_(n)
        line(lc).mode  = 0
        IF (BJTEST(istat_(n),6)) line(lc).gflog = xgf_(n)  ! ALTERNATIVE VALUE
        IF (BJTEST(istat_(n),7)) line(lc).c6log = xc6_(n)  ! ALTERNATIVE VALUE
        line(lc).grad  = gr_(n)
        line(lc).width = DBLE(dw_(n)) * 1.D-8              ! [cm]
        line(lc).rbb_idx  = 0
        IF (.NOT.BJTEST(lf_ctrl,26)) THEN
                      WRITE(CHN,*) 
     *               'Exact line: ',line(lc).id,' at ',line(lc).w,
     *               'W-Range: ',wmin,' - ',wmax,' - ', el_(n)
                      EXCPOT = el_(n)
        ENDIF
        GOTO 999
       ENDIF                  ! if this section skipped than it was NOT exact line formation


       READ(25,REC=1) wref  ! REFERENCE-REC.(START-WAVELENGTH OF EACH REC.)
       maxrec = WREF_NR
       rbb_idx= pop.rbb_nr            ! Act. Idx of Transition in POP-Array (NLTE!)

C---   GET MAXIMUM RECORD-NUMBER (MUST BE < 1901)
       DO i = 1,WREF_NR
        IF (wref(i).EQ.0.) THEN
          maxrec = i
          GOTO 2
        ENDIF
       ENDDO
 
       write(CHN,*) 'Maximum Number of blocks detected !',
     *             ' Error in ',linedata_file

2      CONTINUE
       nrec = 2    !DEFAULT START-RECORD (START-WAVELENGTH = WREF(1))
       DO i = 1,maxrec-1
        IF (wref(i+1).GE.wmin-200.) THEN    ! CHECK W-START OF THE FOLLOWING RECORD
          nrec = i + 1                      ! WREF(I) IS W-START OF REC. (I+1)
          GOTO 3
        ENDIF
       ENDDO


3      CONTINUE
       lc = 0   !LINE-COUNTER

C---   RECORDS ARE SCANNED FOR SELECTED LINES (CHECK BITNR.5 OF LINE-STATUS)
C---   ALL LINES WITH MISSING GF- AND(!) C6- VALUES WILL BE IGNORED

10     CONTINUE
       IF (nrec.GT.maxrec) GOTO 999    !there exist only MAXREC Records
       READ(25,REC=nrec) buf           !READ RECORD WITH LNR LINES

       DO n = 1,LNR

c MB TESTING SECTION ********************************
c      IF (w_(n).lt.6500.0.and.w_(n).gt.6490.0) THEN
c       WRITE(*,*) nrec,n,w_(n),iel1_(n),ion_(n),gf_(n)
c     ENDIF

       IF (dw_(n) .LE. 0.) dw_(n) = STD_WIDTH*w_(n)
       IF (w_(n).EQ.0..OR.(w_(n)-200.).GT. wmax) GOTO 999
       IF ((w_(n) + dw_(n)) .LT. wmin) GOTO 899
       IF ((w_(n) - dw_(n)) .GT. wmax) GOTO 899
              
       IF (JIAND(istat_(n),mask) .EQ. 0) THEN

        IF (BJTEST(lf_ctrl,21)) THEN
         IF (w_(n) .LT. wmin .OR. w_(n) .GT. wmax) GOTO 899
         IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *      'Line ',w_(n),' within selected interval.'
        ENDIF

        nlte_line = .FALSE.

c          write(*,*)'pre', iel2_(n), iel1_(n), pop.atom

C---    CHECK if NLTE-Lineformation selected
        IF (BJTEST(lf_ctrl,20).AND..NOT.BJTEST(lf_ctrl,24)) THEN
         IF (rbb_idx .GT. 0) THEN
          ! Check whether Atom is correct
          IF (iel2_(n).LE.0.AND.iel1_(n).EQ.pop.atom) THEN
            DO i=1,pop.rbb_nr
              IF (ion_(n).EQ.pop.level(pop.rbb(i).ll).ion) THEN
               IF (DABS(w_(n) - 1.D8*pop.rbb(i).w) .LE. XDLAM) THEN
                 IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *                'NLTE-LINEFORMATION AT W= ',w_(n)
                 nlte_line = .TRUE.
                 rbb_idx   = i
                 pop.rbb(rbb_idx).accept = .TRUE.
                 GOTO 20
               ENDIF
              ENDIF
            ENDDO            
          ENDIF
         ENDIF
        ENDIF
  20    CONTINUE

C---    CHECK IF ONLY SELECTED LINES SHOULD BE CONSIDERED (ONLY IF BJTEST(lf_ctrl,2) EQ .TRUE.)
        IF (BJTEST(lf_ctrl,2) .AND. .NOT.BJTEST(ISTAT_(n),9))  GOTO 899

C---    CHECK IF PREDICTED LINES SHOULD BE IGNORED (ONLY IF BJTEST(lf_ctrl,2) EQ .TRUE.)
        IF (BJTEST(lf_ctrl,29) .AND. BJTEST(ISTAT_(n),17))  GOTO 899

C---    CHECK IF ONLY MOLECULAR LINES SHOUD BE CONSIDERED (ONLY IF BJTEST(lf_ctrl,13) EQ .TRUE.)
        IF (BJTEST(lf_ctrl,13).AND. .NOT.BJTEST(ISTAT_(n),15)) GOTO 899

        IF (BJTEST(lf_ctrl,23)) THEN                   ! single lineformation
C        PRINT*, 'Single line formation', w_(n), wline_center
          IF (ABS(w_(n) - wline_center).GT.1.D0) GOTO 899    ! condition of 1A on the wrange from line center
C---    wline_code/200 instead of /100 because high-Z ions were not treated properly (e.g. Ba - 56)
C       MB 2011
          icode2 = NINT(wline_code/200.)
          icode1 = NINT(wline_code - icode2 * 100.)
          icode0 = NINT((wline_code - icode1 - icode2 * 100.) * 100.)
C         PRINT*, icode2, icode1, nelpos(iel1_(n)), iel1_(n), wline_code
          IF (icode0.EQ.0) icode0 = ion_(n)
          IF ( icode2.NE.iel2_(n) .OR. icode1.NE.iel1_(n) .OR.
     *         icode0.NE.ion_(n) ) GOTO 899 
          IF (iel2_(n) .LE. 0) THEN
            IF (nelpos(iel1_(n)).LE.0 .OR. iel1_(n).EQ.1) GOTO 899
             IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,2411) 
     *            lbl(iel1_(n)),w_(n)
          ELSE
          IF (mol_idx(iel1_(n), iel2_(n)).LE.0) GOTO 899
             lbl_tmp = cel(iel1_(n)) // cel(iel2_(n))
             CALL STRCOMPRESS(lbl_tmp,k)
             IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,2411) 
     *            lbl_tmp(1:k),w_(n) 
          ENDIF
          IF (.NOT.BJTEST(lf_ctrl,26)) THEN
                      EXCPOT = el_(n)
                      WRITE(CHN,*) 
     *                'Exact line found at: ',w_(n),' ', 'Elow: ',
     *                 EXCPOT
          ENDIF
          dw_(n) = 0.5D0 * (wmax-wmin)  ! full width
        ENDIF
 2411   FORMAT(' Only Line ',A4,' (',F10.3,')  selected !')
 
        mol_idx_tmp = 0
        IF (iel2_(n) .LE. 0) THEN
ccc MB          PRINT*, 'selected line formation for: ',iel1_(n)
          IF (nelpos(iel1_(n)).LE.0) THEN
c          PRINT*,'nelpos' 
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *             lbl(iel1_(n)),'  Line: ',w_(n),' ignored'
            GOTO 899
          ELSEIF (iel1_(n).EQ.1.AND.w_(n).GT.7000.) THEN     ! Balmerlines are especially considered
            GOTO 899
          ENDIF
        ELSE
          mol_idx_tmp = mol_idx(iel1_(n), iel2_(n))
          IF (mol_idx_tmp .LE. 0) THEN
            lbl_tmp = cel(iel1_(n)) // cel(iel2_(n))
            CALL STRCOMPRESS(lbl_tmp,k)
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *              lbl_tmp(1:k),'  Line: ',w_(n),' ignored'
            GOTO 899
          ENDIF
        ENDIF

        lc = lc + 1
        IF (lc.GE.NLINE) THEN
         IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*)
     *       'Max. number (',NLINE,') of lines exceeded !'
           GOTO 999            ! finished
        ENDIF
        line(lc).w = w_(n) * 1.D-8                         ! [cm]
        IF (lc .GT.1) THEN
          IF (line(lc).w .LT. line(lc-1).w) sort_lines = .TRUE.
        ENDIF

        IF (mol_idx_tmp .LE. 0) THEN
          line(lc).id = iel1_(n)                           ! Atom
        ELSE
          line(lc).id = NATOM + mol_idx_tmp                ! Molecule
        ENDIF

        line(lc).ion   = ion_(n)
        line(lc).elow  = el_(n) * E_SI * 1.e7              ! [erg]
        line(lc).gflog = gf_(n)                  
        line(lc).c4log = c4_(n)    
        line(lc).c6log = c6_(n)
        line(lc).xjl   = xjl_(n)    
        line(lc).xju   = xjh_(n)
        line(lc).mode  = 0

        IF (BJTEST(istat_(n),6)) line(lc).gflog = xgf_(n)  ! ALTERNATIVE VALUE
        IF (BJTEST(istat_(n),7)) line(lc).c6log = xc6_(n)  ! ALTERNATIVE VALUE

        line(lc).grad  = gr_(n)
        line(lc).width = DBLE(dw_(n)) * 1.D-8              ! [cm]

        IF (nlte_line) THEN
          pop.rbb(rbb_idx).accept = .TRUE.
          pop.rbb(rbb_idx).w      = line(lc).w
C          pop.rbb(rbb_idx).dw     = line(lc).width
          line(lc).rbb_idx  = rbb_idx
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,4010) 
     *                    lbl(iel1_(n)),ion_(n),w_(n),
     *                    pop.level(pop.rbb(rbb_idx).ll).desig,
     *                    pop.level(pop.rbb(rbb_idx).ul).desig
 4010     FORMAT(' NLTE Line: ',A2,I1,3X,G13.3,' A   ',
     *             A6,' --> ',A6)
        ELSE
          line(lc).rbb_idx  = 0
        ENDIF
       ENDIF

899    CONTINUE

       ENDDO       ! handling one record

C---  read next record
       nrec = nrec + 1
       GOTO 10

999    CONTINUE                           
       close(25)
C --- ------------------------------------------------------------------------
      ELSE                                 ! ascii linetable format (qline...)
C --- ------------------------------------------------------------------------
       OPEN(UNIT=25,FILE=linedata_file,STATUS='OLD',READONLY)
       lc = 0
       k  = 0
 501   CONTINUE
         k = k + 1
         IF (k.GT.NLINE) GOTO 589            ! finished
         READ(25,*,END=589) line(1).w,
     *              line(1).id,
     *              line(1).ion,
     *              line(1).elow,
     *              line(1).gflog,
     *              line(1).c6log,
     *              line(1).grad,
     *              line(1).width, 
     *              line(1).c4log
        line(1).xjl = -1.
        line(1).xju = -1.

         IF (line(1).width .EQ. 0.D0) line(1).width = 
     *                                  STD_WIDTH*line(1).w
         xw = DABS(line(1).width)
         line(1).rbb_idx  = 0
         IF (line(1).w .GT. (wmax+xw) ) GOTO 589
       IF (line(1).w .LT. (wmin-xw) ) GOTO 501 
       lc = 1
 503   CONTINUE
         k = k + 1
         IF (k.GT.NLINE) GOTO 589            ! finished
         READ(25,*,END=589) line(lc+1).w,
     *              line(lc+1).id,
     *              line(lc+1).ion,
     *              line(lc+1).elow,
     *              line(lc+1).gflog,
     *              line(lc+1).c6log,
     *              line(lc+1).grad,
     *              line(lc+1).width, 
     *              line(lc+1).c4log
         line(lc+1).xjl = -1.
         line(lc+1).xju = -1.
         lc = lc + 1
         IF (line(lc).width .EQ. 0.D0) line(lc).width = 
     *                                        STD_WIDTH*line(lc).w
         xw = DABS(line(lc).width)
         line(lc).rbb_idx = 0
       IF (line(lc).w .LE. (wmax+xw)) GOTO 503
 589   CONTINUE
       CLOSE(25)

       rbb_idx = pop.rbb_nr
       DO n=1,lc 
        IF (line(n).width .LE. 0.) THEN
          line(n).id = 0
          GOTO 595
        ENDIF
        nlte_line = .FALSE.
C---    CHECK if NLTE-Lineformation selected
        IF (BJTEST(lf_ctrl,20).AND..NOT.BJTEST(lf_ctrl,24)) THEN
         IF (rbb_idx .GT. 0) THEN
          ! Check if Atom is right
          IF (line(n).id .EQ. pop.atom) THEN
cc            DO i=rbb_idx-1,1,-1
            DO i=1,pop.rbb_nr
cv             IF (.NOT. pop.rbb(i).accept) THEN
              IF (line(n).ion .EQ. pop.level(pop.rbb(i).ll).ion) THEN
               IF (DABS(line(n).w - 1.D8*pop.rbb(i).w).LE.XDLAM) THEN
                 IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,*) 
     *               'NLTE-LINEFORMATION AT W= ',line(n).w
                  nlte_line = .TRUE.
                  rbb_idx   = i
                  pop.rbb(rbb_idx).accept = .TRUE.
                  GOTO 591
               ENDIF
              ENDIF
cv             ENDIF
            ENDDO
          ENDIF
         ENDIF
        ENDIF
 591    CONTINUE
        IF (line(n).id .LE. NATOM) THEN
          IF (nelpos(line(n).id).LE.0) THEN
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *             lbl(line(n).id),'  Line: ',line(n).w,' ignored'
            line(n).id=0
            GOTO 595
          ELSEIF (line(n).id.EQ.1) THEN     ! H-Lines are especially considered
            line(n).id=0
            GOTO 595
          ENDIF
        ELSE
          mol_idx_tmp = mol_idx(MOD(line(n).id,100),line(n).id/100)
          IF (mol_idx_tmp .LE. 0) THEN
            lbl_tmp = cel(MOD(line(n).id,100))//cel(line(n).id/100)
            CALL STRCOMPRESS(lbl_tmp,k)
            IF (.NOT.BJTEST(lf_ctrl,26)) write(CHN,*) 
     *              lbl_tmp(1:k),'  Line: ',line(n).w,' ignored'
            line(n).id=0
            GOTO 595
          ENDIF
          line(n).id = NATOM + mol_idx_tmp                ! Molecule
        ENDIF
        line(n).w = line(n).w * 1.D-8                         ! [cm]
        IF (nlte_line) THEN
          pop.rbb(rbb_idx).accept = .TRUE.
          pop.rbb(rbb_idx).w      = line(n).w
C          pop.rbb(rbb_idx).dw     = line(n).width
          line(n).rbb_idx = rbb_idx
          IF (.NOT.BJTEST(lf_ctrl,26)) WRITE(CHN,4010) 
     *                    lbl(line(n).id),line(n).ion,line(n).w,
     *                    pop.level(pop.rbb(rbb_idx).ll).desig,
     *                    pop.level(pop.rbb(rbb_idx).ul).desig
        ELSE
          line(n).rbb_idx  = 0
        ENDIF
        line(n).elow = line(n).elow * E_SI * 1.e7        ! [erg]
 595    CONTINUE
        IF (line(n).id .LE. 0) sort_lines = .TRUE.
        line(n).width = line(n).width * 1.D-8              ! [cm]
        line(n).mode  = 0
       ENDDO

       IF (sort_lines) THEN
        sort_lines = .FALSE.
        linemax = lc
        lc = 0
        DO n=1,linemax
           IF (line(n).id .GT. 0) THEN
             lc = lc + 1
             line(lc) = line(n)
           ENDIF
        ENDDO
       ENDIF
      ENDIF                                ! linetable format

      wmin    = wmin * 1.e-8               ! [A] -> [cm]
      wmax    = wmax * 1.e-8               ! [A] -> [cm]
      linemax = lc

      IF (.NOT.sort_lines) GOTO 2050 
 2005 CONTINUE   ! jump only here for sorting lines

C --- SORT LINE                             
C --- see Heapsort: Numerical Recipes (1989,P.231)       
      n = linemax/2 + 1
      k = linemax
 2010 CONTINUE
        IF (n.GT.1) THEN
           n = n - 1
           line_tmp = line(n)
        ELSE
           line_tmp = line(k)
           line(k)  = line(1)
           k        = k - 1
           IF (k.EQ.1) THEN
              line(1) = line_tmp
              GOTO 2050                ! return
           ENDIF
        ENDIF
        i = n
        j = n + n
 2020   IF (j.LE.k) THEN
          IF (j.LT.k) THEN
            IF ( line(j).w .LT. line(j+1).w ) j = j + 1            
          ENDIF
          IF ( line_tmp.w .LT. line(j).w ) THEN
             line(i) = line(j)
             i = j
             j = j + j
          ELSE
             j = k + 1
          ENDIF
        GOTO 2020
        ENDIF
        line(i) = line_tmp
      GOTO 2010
 2050 CONTINUE

      IF (.NOT.nlte_line) THEN
       DO lc=1,linemax
         IF (line(lc).rbb_idx.GT.0) line(lc).rbb_idx = 0
       ENDDO
      ENDIF


      RETURN
      END
