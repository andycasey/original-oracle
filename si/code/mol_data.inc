C >>>
C     --------------- Lineformation - Includefile ---------------------------
C
C     Filename  :  MOL_DATA.INC
C     Purpose   :  Fixed Molecular Data for Lineformationprogramm Linfor
C     Modul     :  Fortran - Includefile
C     Data      :  cmol     Character-Label for molecule
C                  molmass  molecular mass     (atomic mass units)
C                  mol_idx  =0, if molecule is not considered, molecule-index
C                           otherwise.
C
C     Author    : JR
C     Date      : 20.09.92
C
C     History   : -
C
C <<<
C --------------------------------------------------------------------------

 
      INTEGER*4  NATOM_2
      PARAMETER (NATOM_2 = NATOM*NATOM)

c      DATA mol_idx         / NATOM_2 * 0 /          ! INIT
      DATA mol_idx( 1, 1)  /  1/                    ! H2
      DATA mol_idx( 6, 1)  /  2/                    ! CH
      DATA mol_idx( 1, 6)  /  2/                    ! CH
      DATA mol_idx( 7, 1)  /  3/                    ! NH
      DATA mol_idx( 1, 7)  /  3/                    ! NH

      DATA mol_idx( 8, 1)  /  4/                    ! OH
      DATA mol_idx( 1, 8)  /  4/                    ! OH

      DATA mol_idx( 6, 6)  /  5/                    ! C2

      DATA mol_idx( 6, 7)  /  6/                    ! CN
      DATA mol_idx( 7, 6)  /  6/                    ! CN

      DATA mol_idx( 6, 8)  /  7/                    ! CO
      DATA mol_idx( 8, 6)  /  7/                    ! CO
      DATA mol_idx( 7, 7)  /  8/                    ! N2
      DATA mol_idx( 7, 8)  /  9/                    ! NO
      DATA mol_idx( 8, 7)  /  9/                    ! NO

      DATA mol_idx( 8, 8)  / 10/                    ! O2

      DATA mol_idx(12, 1)  / 11/                    ! MgH
      DATA mol_idx( 1,12)  / 11/                    ! MgH

      DATA mol_idx(14, 1)  / 12/                    ! SiH
      DATA mol_idx( 1,14)  / 12/                    ! SiH

      DATA mol_idx(20, 1)  / 13/                    ! CaH
      DATA mol_idx( 1,20)  / 13/                    ! CaH

      DATA mol_idx(22, 8)  / 14/                    ! TiO
      DATA mol_idx( 8,22)  / 14/                    ! TiO

      DATA mol_idx(23, 8)  / 15/                    ! VO
      DATA mol_idx( 8,23)  / 15/                    ! VO

      DATA mol_idx(40, 8)  / 16/                    ! ZrO
      DATA mol_idx( 8,40)  / 16/                    ! ZrO

      DATA mol_idx(26, 1)  / 17/                    ! FeH
      DATA mol_idx( 1,26)  / 17/                    ! FeH

      DATA cmol( 1),molmass( 1),mel1( 1),mel2( 1) /'H2 ', 2.016, 1, 1/
      DATA cmol( 2),molmass( 2),mel1( 2),mel2( 2) /'CH ',13.019, 6, 1/
      DATA cmol( 3),molmass( 3),mel1( 3),mel2( 3) /'NH ',15.015, 7, 1/
      DATA cmol( 4),molmass( 4),mel1( 4),mel2( 4) /'OH ',17.007, 8, 1/
      DATA cmol( 5),molmass( 5),mel1( 5),mel2( 5) /'C2 ',24.022, 6, 6/
      DATA cmol( 6),molmass( 6),mel1( 6),mel2( 6) /'CN ',26.018, 6, 7/
      DATA cmol( 7),molmass( 7),mel1( 7),mel2( 7) /'CO ',28.010, 6, 8/
      DATA cmol( 8),molmass( 8),mel1( 8),mel2( 8) /'N2 ',28.014, 7, 7/
      DATA cmol( 9),molmass( 9),mel1( 9),mel2( 9) /'NO ',30.006, 7, 8/
      DATA cmol(10),molmass(10),mel1(10),mel2(10) /'O2 ',31.998, 8, 8/
      DATA cmol(11),molmass(11),mel1(11),mel2(11) /'MgH',25.313,12, 1/
      DATA cmol(12),molmass(12),mel1(12),mel2(12) /'SiH',29.094,14, 1/
      DATA cmol(13),molmass(13),mel1(13),mel2(13) /'CaH',41.088,20, 1/
      DATA cmol(14),molmass(14),mel1(14),mel2(14) /'TiO',63.899,22, 8/
      DATA cmol(15),molmass(15),mel1(15),mel2(15) /'VO ',66.941,23, 8/
      DATA cmol(16),molmass(16),mel1(16),mel2(16) /'ZrO',107.22,40, 8/
      DATA cmol(17),molmass(17),mel1(17),mel2(17) /'FeH',56.853,26, 1/