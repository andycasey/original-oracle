C >>>
C     --------------- Lineformation - Includefile ---------------------------
C
C     Filename  :  ELEMDATA.INC
C     Purpose   :  Fixed Data for Lineformationprogramm Linfor
C     Modul     :  Fortran - Includefile
C               :  LINFOR.FOR; Variable Declarations in LINFOR_DECL.INC
C
C     Data      : 
C                  cel      Character-Label for Element
C                  epss     Solar Abundance (log. scale !, see Ref.)
C                  amass    Atomic Mass     (atomic mass unit)
C                  nion     Number of ionisation-levels taken into account
C                           for each element (0 = not actually used)
C
C     Reference : Solar Abundances by:
C			Asplund et al. 2009, ARAA, 47, 481 (H,He,C,N,O+meteoritic)
C			Scott et al., 2014, A&A, in press (F-Ca)
C			Scott et al., 2014, A&A, in press (Sc-Ni)
C			Grevesse et al., 2014, A&A, in press (Cu-Th)
C
C     Author    : M. Asplund
C     Date      : 29.08.2014
C
C     History   : -
C
C <<<
C --------------------------------------------------------------------------

C     The atomic number is here mapped to the equivalent array-element-Index.

      DATA cel(1) ,epss(1) ,amass(1) ,nion(1)  /' H',12.00,  1.008,3/
      DATA cel(2) ,epss(2) ,amass(2) ,nion(2)  /'He',10.93,  4.003,2/
      DATA cel(3) ,epss(3) ,amass(3) ,nion(3)  /'Li', 1.05,  6.939,2/
      DATA cel(4) ,epss(4) ,amass(4) ,nion(4)  /'Be', 1.38,  9.012,3/
      DATA cel(5) ,epss(5) ,amass(5) ,nion(5)  /' B', 2.70, 10.811,3/
      DATA cel(6) ,epss(6) ,amass(6) ,nion(6)  /' C', 8.43, 12.011,2/
      DATA cel(7) ,epss(7) ,amass(7) ,nion(7)  /' N', 7.83, 14.007,2/
      DATA cel(8) ,epss(8) ,amass(8) ,nion(8)  /' O', 8.69, 15.999,2/
      DATA cel(9) ,epss(9) ,amass(9) ,nion(9)  /' F', 4.40, 18.998,2/
      DATA cel(10),epss(10),amass(10),nion(10) /'Ne', 7.93, 20.179,2/
      DATA cel(11),epss(11),amass(11),nion(11) /'Na', 6.21, 22.990,2/
      DATA cel(12),epss(12),amass(12),nion(12) /'Mg', 7.59, 24.305,3/
      DATA cel(13),epss(13),amass(13),nion(13) /'Al', 6.43, 26.982,3/
      DATA cel(14),epss(14),amass(14),nion(14) /'Si', 7.51, 28.086,3/
      DATA cel(15),epss(15),amass(15),nion(15) /' P', 5.41, 30.974,3/
      DATA cel(16),epss(16),amass(16),nion(16) /' S', 7.12, 32.060,3/
      DATA cel(17),epss(17),amass(17),nion(17) /'Cl', 5.50, 35.453,2/
      DATA cel(18),epss(18),amass(18),nion(18) /'Ar', 6.40, 39.948,2/
      DATA cel(19),epss(19),amass(19),nion(19) /' K', 5.04, 39.100,2/
      DATA cel(20),epss(20),amass(20),nion(20) /'Ca', 6.32, 40.080,3/
      DATA cel(21),epss(21),amass(21),nion(21) /'Sc', 3.16, 44.958,3/
      DATA cel(22),epss(22),amass(22),nion(22) /'Ti', 4.93, 47.900,3/
      DATA cel(23),epss(23),amass(23),nion(23) /' V', 3.89, 50.944,3/
      DATA cel(24),epss(24),amass(24),nion(24) /'Cr', 5.62, 51.996,3/
      DATA cel(25),epss(25),amass(25),nion(25) /'Mn', 5.42, 54.938,3/
      DATA cel(26),epss(26),amass(26),nion(26) /'Fe', 7.47, 55.847,3/
      DATA cel(27),epss(27),amass(27),nion(27) /'Co', 4.93, 58.940,3/
      DATA cel(28),epss(28),amass(28),nion(28) /'Ni', 6.20, 58.710,3/
      DATA cel(29),epss(29),amass(29),nion(29) /'Cu', 4.18, 63.550,3/
      DATA cel(30),epss(30),amass(30),nion(30) /'Zn', 4.56, 65.370,3/
      DATA cel(31),epss(31),amass(31),nion(31) /'Ga', 3.02, 69.723,3/
      DATA cel(32),epss(32),amass(32),nion(32) /'Ge', 3.63, 71.610,3/
      DATA cel(33),epss(33),amass(33),nion(33) /'As', 2.30, 74.922,3/
      DATA cel(34),epss(34),amass(34),nion(34) /'Se', 3.34, 78.960,3/
      DATA cel(35),epss(35),amass(35),nion(35) /'Br', 2.54, 79.904,3/
      DATA cel(36),epss(36),amass(36),nion(36) /'Kr', 3.25, 83.800,2/
      DATA cel(37),epss(37),amass(37),nion(37) /'Rb', 2.47, 85.468,3/
      DATA cel(38),epss(38),amass(38),nion(38) /'Sr', 2.83, 87.630,3/
      DATA cel(39),epss(39),amass(39),nion(39) /' Y', 2.21, 88.908,3/
      DATA cel(40),epss(40),amass(40),nion(40) /'Zr', 2.59, 91.220,3/
      DATA cel(41),epss(41),amass(41),nion(41) /'Nb', 1.47, 92.906,3/
      DATA cel(42),epss(42),amass(42),nion(42) /'Mo', 1.88, 95.940,3/
      DATA cel(43),epss(43),amass(43),nion(43) /'Tc’,-9.99, 98.906,3/
      DATA cel(44),epss(44),amass(44),nion(44) /'Ru', 1.75,101.070,3/
      DATA cel(45),epss(45),amass(45),nion(45) /'Rh', 0.89,102.906,3/
      DATA cel(46),epss(46),amass(46),nion(46) /'Pd', 1.55,106.420,3/
      DATA cel(47),epss(47),amass(47),nion(47) /'Ag', 0.96,107.868,3/
      DATA cel(48),epss(48),amass(48),nion(48) /'Cd', 1.77,112.411,3/
      DATA cel(49),epss(49),amass(49),nion(49) /'In', 0.80,114.818,3/
      DATA cel(50),epss(50),amass(50),nion(50) /'Sn', 2.02,118.710,3/
      DATA cel(51),epss(51),amass(51),nion(51) /'Sb', 1.01,121.757,3/
      DATA cel(52),epss(52),amass(52),nion(52) /'Te', 2.18,127.600,3/
      DATA cel(53),epss(53),amass(53),nion(53) /' I', 1.55,126.904,3/
      DATA cel(54),epss(54),amass(54),nion(54) /'Xe', 2.24,131.290,3/
      DATA cel(55),epss(55),amass(55),nion(55) /'Cs', 1.08,132.905,3/
      DATA cel(56),epss(56),amass(56),nion(56) /'Ba', 2.25,137.350,3/
      DATA cel(57),epss(57),amass(57),nion(57) /'La', 1.11,138.920,3/
      DATA cel(58),epss(58),amass(58),nion(58) /'Ce', 1.58,140.130,3/
      DATA cel(59),epss(59),amass(59),nion(59) /'Pr', 0.72,140.908,3/
      DATA cel(60),epss(60),amass(60),nion(60) /'Nd', 1.42,144.250,3/
      DATA cel(61),epss(61),amass(61),nion(61) /'Pm’,-9.99,146.915,3/
      DATA cel(62),epss(62),amass(62),nion(62) /'Sm', 0.95,150.360,3/
      DATA cel(63),epss(63),amass(63),nion(63) /'Eu', 0.52,151.960,3/
      DATA cel(64),epss(64),amass(64),nion(64) /'Gd', 1.08,157.250,3/
      DATA cel(65),epss(65),amass(65),nion(65) /'Tb', 0.31,158.925,3/
      DATA cel(66),epss(66),amass(66),nion(66) /'Dy', 1.10,162.500,3/
      DATA cel(67),epss(67),amass(67),nion(67) /'Ho', 0.48,164.930,3/
      DATA cel(68),epss(68),amass(68),nion(68) /'Er', 0.93,167.260,3/
      DATA cel(69),epss(69),amass(69),nion(69) /'Tm', 0.11,168.934,3/
      DATA cel(70),epss(70),amass(70),nion(70) /'Yb', 0.85,173.040,3/
      DATA cel(71),epss(71),amass(71),nion(71) /'Lu', 0.10,174.967,3/
      DATA cel(72),epss(72),amass(72),nion(72) /'Hf', 0.85,178.490,3/
      DATA cel(73),epss(73),amass(73),nion(73) /'Ta’,-0.12,180.948,3/
      DATA cel(74),epss(74),amass(74),nion(74) /' W', 0.83,183.840,3/
      DATA cel(75),epss(75),amass(75),nion(75) /'Re', 0.26,186.207,3/
      DATA cel(76),epss(76),amass(76),nion(76) /'Os', 1.40,190.230,3/
      DATA cel(77),epss(77),amass(77),nion(77) /'Ir', 1.42,192.220,3/
      DATA cel(78),epss(78),amass(78),nion(78) /'Pt', 1.62,195.080,3/
      DATA cel(79),epss(79),amass(79),nion(79) /'Au', 0.91,196.967,3/
      DATA cel(80),epss(80),amass(80),nion(80) /'Hg', 1.17,200.590,3/
      DATA cel(81),epss(81),amass(81),nion(81) /'Tl', 0.90,204.383,3/
      DATA cel(82),epss(82),amass(82),nion(82) /'Pb', 1.92,207.200,3/
      DATA cel(83),epss(83),amass(83),nion(83) /'Bi', 0.65,208.980,3/
      DATA cel(84),epss(84),amass(84),nion(84) /'Po',-9.99,208.982,3/
      DATA cel(85),epss(85),amass(85),nion(85) /'At',-9.99,209.987,3/
      DATA cel(86),epss(86),amass(86),nion(86) /'Rn',-9.99,222.018,3/
      DATA cel(87),epss(87),amass(87),nion(87) /'Fr',-9.99,223.020,3/
      DATA cel(88),epss(88),amass(88),nion(88) /'Ra',-9.99,226.025,3/
      DATA cel(89),epss(89),amass(89),nion(89) /'Ac',-9.99,227.028,3/
      DATA cel(90),epss(90),amass(90),nion(90) /'Th', 0.03,232.038,3/
      DATA cel(91),epss(91),amass(91),nion(91) /'Pa',-9.99,231.036,3/
      DATA cel(92),epss(92),amass(92),nion(92) /' U', 0.00,238.029,3/

