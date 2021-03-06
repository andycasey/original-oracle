C --- INCLUDE FILE ATMGRID.INC FOR LINEFORM.F ---------
C --- This Tables provides all accessible atmosphere grid points
C
      INTEGER*4 NTEFFS,NLOGGS,NZS
      PARAMETER (NTEFFS =  7) ! Number of Teff-gridpoints
      PARAMETER (NLOGGS =  8) ! Number of log(g)-gridpoints
      PARAMETER (NZS    =  7) ! Number of metallicity-gridpoints
      REAL*4 atm_teff(NTEFFS)
      REAL*4 atm_logg(NLOGGS)
      REAL*4 atm_z   (NZS)
      CHARACTER*25 atm_filename (NTEFFS,NLOGGS,NZS)
C
C --- AVAILABLE EFFECTIVE TEMPERATURES ----------------
      DATA atm_teff(  1)  / 4500.0/
      DATA atm_teff(  2)  / 5000.0/
      DATA atm_teff(  3)  / 5500.0/
      DATA atm_teff(  4)  / 5777.0/
      DATA atm_teff(  5)  / 6000.0/
      DATA atm_teff(  6)  / 6500.0/
      DATA atm_teff(  7)  / 7000.0/
C
C --- AVAILABLE GRAVITIES -----------------------------
      DATA atm_logg(  1)  / 2.000/
      DATA atm_logg(  2)  / 2.500/
      DATA atm_logg(  3)  / 3.000/
      DATA atm_logg(  4)  / 3.500/
      DATA atm_logg(  5)  / 4.000/
      DATA atm_logg(  6)  / 4.400/
      DATA atm_logg(  7)  / 4.500/
      DATA atm_logg(  8)  / 5.000/
C
C --- AVAILABLE METALLICITIES -------------------------
      DATA atm_z   (  1)  /-4.000/
      DATA atm_z   (  2)  /-3.000/
      DATA atm_z   (  3)  /-2.000/
      DATA atm_z   (  4)  /-1.000/
      DATA atm_z   (  5)  /-0.500/
      DATA atm_z   (  6)  /-0.000/
      DATA atm_z   (  7)  /+0.500/
C
C --- AVAILABLE ATMOSPHERES ---------------------------
      DATA atm_filename(  1,  1,  1)  /'                         '/
      DATA atm_filename(  1,  1,  2)  /'                         '/
      DATA atm_filename(  1,  1,  3)  /'t4500_2.0_-2.00_80.dat   '/
      DATA atm_filename(  1,  1,  4)  /'                         '/
      DATA atm_filename(  1,  1,  5)  /'                         '/
      DATA atm_filename(  1,  1,  6)  /'t4500_2.0_0.00_80.dat    '/
      DATA atm_filename(  1,  1,  7)  /'                         '/
      DATA atm_filename(  1,  2,  1)  /'                         '/
      DATA atm_filename(  1,  2,  2)  /'                         '/
      DATA atm_filename(  1,  2,  3)  /'                         '/
      DATA atm_filename(  1,  2,  4)  /'                         '/
      DATA atm_filename(  1,  2,  5)  /'                         '/
      DATA atm_filename(  1,  2,  6)  /'                         '/
      DATA atm_filename(  1,  2,  7)  /'                         '/
      DATA atm_filename(  1,  3,  1)  /'                         '/
      DATA atm_filename(  1,  3,  2)  /'                         '/
      DATA atm_filename(  1,  3,  3)  /'                         '/
      DATA atm_filename(  1,  3,  4)  /'                         '/
      DATA atm_filename(  1,  3,  5)  /'                         '/
      DATA atm_filename(  1,  3,  6)  /'                         '/
      DATA atm_filename(  1,  3,  7)  /'                         '/
      DATA atm_filename(  1,  4,  1)  /'                         '/
      DATA atm_filename(  1,  4,  2)  /'                         '/
      DATA atm_filename(  1,  4,  3)  /'                         '/
      DATA atm_filename(  1,  4,  4)  /'                         '/
      DATA atm_filename(  1,  4,  5)  /'                         '/
      DATA atm_filename(  1,  4,  6)  /'                         '/
      DATA atm_filename(  1,  4,  7)  /'                         '/
      DATA atm_filename(  1,  5,  1)  /'                         '/
      DATA atm_filename(  1,  5,  2)  /'                         '/
      DATA atm_filename(  1,  5,  3)  /'                         '/
      DATA atm_filename(  1,  5,  4)  /'                         '/
      DATA atm_filename(  1,  5,  5)  /'                         '/
      DATA atm_filename(  1,  5,  6)  /'                         '/
      DATA atm_filename(  1,  5,  7)  /'                         '/
      DATA atm_filename(  1,  6,  1)  /'                         '/
      DATA atm_filename(  1,  6,  2)  /'                         '/
      DATA atm_filename(  1,  6,  3)  /'                         '/
      DATA atm_filename(  1,  6,  4)  /'                         '/
      DATA atm_filename(  1,  6,  5)  /'                         '/
      DATA atm_filename(  1,  6,  6)  /'                         '/
      DATA atm_filename(  1,  6,  7)  /'                         '/
      DATA atm_filename(  1,  7,  1)  /'                         '/
      DATA atm_filename(  1,  7,  2)  /'                         '/
      DATA atm_filename(  1,  7,  3)  /'                         '/
      DATA atm_filename(  1,  7,  4)  /'                         '/
      DATA atm_filename(  1,  7,  5)  /'                         '/
      DATA atm_filename(  1,  7,  6)  /'                         '/
      DATA atm_filename(  1,  7,  7)  /'                         '/
      DATA atm_filename(  1,  8,  1)  /'                         '/
      DATA atm_filename(  1,  8,  2)  /'                         '/
      DATA atm_filename(  1,  8,  3)  /'t4500_5.0_-2.00_80.dat   '/
      DATA atm_filename(  1,  8,  4)  /'                         '/
      DATA atm_filename(  1,  8,  5)  /'                         '/
      DATA atm_filename(  1,  8,  6)  /'t4500_5.0_0.00_80.dat    '/
      DATA atm_filename(  1,  8,  7)  /'                         '/
      DATA atm_filename(  2,  1,  1)  /'                         '/
      DATA atm_filename(  2,  1,  2)  /'                         '/
      DATA atm_filename(  2,  1,  3)  /'                         '/
      DATA atm_filename(  2,  1,  4)  /'                         '/
      DATA atm_filename(  2,  1,  5)  /'                         '/
      DATA atm_filename(  2,  1,  6)  /'                         '/
      DATA atm_filename(  2,  1,  7)  /'                         '/
      DATA atm_filename(  2,  2,  1)  /'                         '/
      DATA atm_filename(  2,  2,  2)  /'                         '/
      DATA atm_filename(  2,  2,  3)  /'                         '/
      DATA atm_filename(  2,  2,  4)  /'                         '/
      DATA atm_filename(  2,  2,  5)  /'                         '/
      DATA atm_filename(  2,  2,  6)  /'                         '/
      DATA atm_filename(  2,  2,  7)  /'                         '/
      DATA atm_filename(  2,  3,  1)  /'t5000_3.0_-4.00_80.dat   '/
      DATA atm_filename(  2,  3,  2)  /'t5000_3.0_-3.00_80.dat   '/
      DATA atm_filename(  2,  3,  3)  /'t5000_3.0_-2.00_80.dat   '/
      DATA atm_filename(  2,  3,  4)  /'t5000_3.0_-1.00_80.dat   '/
      DATA atm_filename(  2,  3,  5)  /'t5000_3.0_-0.50_80.dat   '/
      DATA atm_filename(  2,  3,  6)  /'t5000_3.0_0.00_80.dat    '/
      DATA atm_filename(  2,  3,  7)  /'t5000_3.0_0.50_80.dat    '/
      DATA atm_filename(  2,  4,  1)  /'t5000_3.5_-4.00_80.dat   '/
      DATA atm_filename(  2,  4,  2)  /'t5000_3.5_-3.00_80.dat   '/
      DATA atm_filename(  2,  4,  3)  /'t5000_3.5_-2.00_80.dat   '/
      DATA atm_filename(  2,  4,  4)  /'t5000_3.5_-1.00_80.dat   '/
      DATA atm_filename(  2,  4,  5)  /'t5000_3.5_-0.50_80.dat   '/
      DATA atm_filename(  2,  4,  6)  /'t5000_3.5_0.00_80.dat    '/
      DATA atm_filename(  2,  4,  7)  /'t5000_3.5_0.50_80.dat    '/
      DATA atm_filename(  2,  5,  1)  /'t5000_4.0_-4.00_80.dat   '/
      DATA atm_filename(  2,  5,  2)  /'t5000_4.0_-3.00_80.dat   '/
      DATA atm_filename(  2,  5,  3)  /'t5000_4.0_-2.00_80.dat   '/
      DATA atm_filename(  2,  5,  4)  /'t5000_4.0_-1.00_80.dat   '/
      DATA atm_filename(  2,  5,  5)  /'t5000_4.0_-0.50_80.dat   '/
      DATA atm_filename(  2,  5,  6)  /'t5000_4.0_0.00_80.dat    '/
      DATA atm_filename(  2,  5,  7)  /'t5000_4.0_0.50_80.dat    '/
      DATA atm_filename(  2,  6,  1)  /'                         '/
      DATA atm_filename(  2,  6,  2)  /'                         '/
      DATA atm_filename(  2,  6,  3)  /'                         '/
      DATA atm_filename(  2,  6,  4)  /'                         '/
      DATA atm_filename(  2,  6,  5)  /'                         '/
      DATA atm_filename(  2,  6,  6)  /'                         '/
      DATA atm_filename(  2,  6,  7)  /'                         '/
      DATA atm_filename(  2,  7,  1)  /'t5000_4.5_-4.00_80.dat   '/
      DATA atm_filename(  2,  7,  2)  /'t5000_4.5_-3.00_80.dat   '/
      DATA atm_filename(  2,  7,  3)  /'t5000_4.5_-2.00_80.dat   '/
      DATA atm_filename(  2,  7,  4)  /'t5000_4.5_-1.00_80.dat   '/
      DATA atm_filename(  2,  7,  5)  /'t5000_4.5_-0.50_80.dat   '/
      DATA atm_filename(  2,  7,  6)  /'t5000_4.5_0.00_80.dat    '/
      DATA atm_filename(  2,  7,  7)  /'t5000_4.5_0.50_80.dat    '/
      DATA atm_filename(  2,  8,  1)  /'t5000_5.0_-4.00_80.dat   '/
      DATA atm_filename(  2,  8,  2)  /'t5000_5.0_-3.00_80.dat   '/
      DATA atm_filename(  2,  8,  3)  /'t5000_5.0_-2.00_80.dat   '/
      DATA atm_filename(  2,  8,  4)  /'t5000_5.0_-1.00_80.dat   '/
      DATA atm_filename(  2,  8,  5)  /'t5000_5.0_-0.50_80.dat   '/
      DATA atm_filename(  2,  8,  6)  /'t5000_5.0_0.00_80.dat    '/
      DATA atm_filename(  2,  8,  7)  /'t5000_5.0_0.50_80.dat    '/
      DATA atm_filename(  3,  1,  1)  /'                         '/
      DATA atm_filename(  3,  1,  2)  /'                         '/
      DATA atm_filename(  3,  1,  3)  /'                         '/
      DATA atm_filename(  3,  1,  4)  /'                         '/
      DATA atm_filename(  3,  1,  5)  /'                         '/
      DATA atm_filename(  3,  1,  6)  /'                         '/
      DATA atm_filename(  3,  1,  7)  /'                         '/
      DATA atm_filename(  3,  2,  1)  /'                         '/
      DATA atm_filename(  3,  2,  2)  /'                         '/
      DATA atm_filename(  3,  2,  3)  /'t5500_2.5_-2.00_80.dat   '/
      DATA atm_filename(  3,  2,  4)  /'                         '/
      DATA atm_filename(  3,  2,  5)  /'                         '/
      DATA atm_filename(  3,  2,  6)  /'                         '/
      DATA atm_filename(  3,  2,  7)  /'                         '/
      DATA atm_filename(  3,  3,  1)  /'t5500_3.0_-4.00_80.dat   '/
      DATA atm_filename(  3,  3,  2)  /'t5500_3.0_-3.00_80.dat   '/
      DATA atm_filename(  3,  3,  3)  /'t5500_3.0_-2.00_80.dat   '/
      DATA atm_filename(  3,  3,  4)  /'t5500_3.0_-1.00_80.dat   '/
      DATA atm_filename(  3,  3,  5)  /'t5500_3.0_-0.50_80.dat   '/
      DATA atm_filename(  3,  3,  6)  /'t5500_3.0_0.00_80.dat    '/
      DATA atm_filename(  3,  3,  7)  /'                         '/
      DATA atm_filename(  3,  4,  1)  /'t5500_3.5_-4.00_80.dat   '/
      DATA atm_filename(  3,  4,  2)  /'t5500_3.5_-3.00_80.dat   '/
      DATA atm_filename(  3,  4,  3)  /'t5500_3.5_-2.00_80.dat   '/
      DATA atm_filename(  3,  4,  4)  /'t5500_3.5_-1.00_80.dat   '/
      DATA atm_filename(  3,  4,  5)  /'t5500_3.5_-0.50_80.dat   '/
      DATA atm_filename(  3,  4,  6)  /'t5500_3.5_0.00_80.dat    '/
      DATA atm_filename(  3,  4,  7)  /'t5500_3.5_0.50_80.dat    '/
      DATA atm_filename(  3,  5,  1)  /'t5500_4.0_-4.00_80.dat   '/
      DATA atm_filename(  3,  5,  2)  /'t5500_4.0_-3.00_80.dat   '/
      DATA atm_filename(  3,  5,  3)  /'t5500_4.0_-2.00_80.dat   '/
      DATA atm_filename(  3,  5,  4)  /'t5500_4.0_-1.00_80.dat   '/
      DATA atm_filename(  3,  5,  5)  /'t5500_4.0_-0.50_80.dat   '/
      DATA atm_filename(  3,  5,  6)  /'t5500_4.0_0.00_80.dat    '/
      DATA atm_filename(  3,  5,  7)  /'t5500_4.0_0.50_80.dat    '/
      DATA atm_filename(  3,  6,  1)  /'                         '/
      DATA atm_filename(  3,  6,  2)  /'                         '/
      DATA atm_filename(  3,  6,  3)  /'                         '/
      DATA atm_filename(  3,  6,  4)  /'                         '/
      DATA atm_filename(  3,  6,  5)  /'                         '/
      DATA atm_filename(  3,  6,  6)  /'                         '/
      DATA atm_filename(  3,  6,  7)  /'                         '/
      DATA atm_filename(  3,  7,  1)  /'t5500_4.5_-4.00_80.dat   '/
      DATA atm_filename(  3,  7,  2)  /'t5500_4.5_-3.00_80.dat   '/
      DATA atm_filename(  3,  7,  3)  /'t5500_4.5_-2.00_80.dat   '/
      DATA atm_filename(  3,  7,  4)  /'t5500_4.5_-1.00_80.dat   '/
      DATA atm_filename(  3,  7,  5)  /'t5500_4.5_-0.50_80.dat   '/
      DATA atm_filename(  3,  7,  6)  /'t5500_4.5_0.00_80.dat    '/
      DATA atm_filename(  3,  7,  7)  /'t5500_4.5_0.50_80.dat    '/
      DATA atm_filename(  3,  8,  1)  /'t5500_5.0_-4.00_80.dat   '/
      DATA atm_filename(  3,  8,  2)  /'t5500_5.0_-3.00_80.dat   '/
      DATA atm_filename(  3,  8,  3)  /'t5500_5.0_-2.00_80.dat   '/
      DATA atm_filename(  3,  8,  4)  /'t5500_5.0_-1.00_80.dat   '/
      DATA atm_filename(  3,  8,  5)  /'t5500_5.0_-0.50_80.dat   '/
      DATA atm_filename(  3,  8,  6)  /'t5500_5.0_0.00_80.dat    '/
      DATA atm_filename(  3,  8,  7)  /'t5500_5.0_0.50_80.dat    '/
      DATA atm_filename(  4,  1,  1)  /'                         '/
      DATA atm_filename(  4,  1,  2)  /'                         '/
      DATA atm_filename(  4,  1,  3)  /'                         '/
      DATA atm_filename(  4,  1,  4)  /'                         '/
      DATA atm_filename(  4,  1,  5)  /'                         '/
      DATA atm_filename(  4,  1,  6)  /'                         '/
      DATA atm_filename(  4,  1,  7)  /'                         '/
      DATA atm_filename(  4,  2,  1)  /'                         '/
      DATA atm_filename(  4,  2,  2)  /'                         '/
      DATA atm_filename(  4,  2,  3)  /'                         '/
      DATA atm_filename(  4,  2,  4)  /'                         '/
      DATA atm_filename(  4,  2,  5)  /'                         '/
      DATA atm_filename(  4,  2,  6)  /'                         '/
      DATA atm_filename(  4,  2,  7)  /'                         '/
      DATA atm_filename(  4,  3,  1)  /'                         '/
      DATA atm_filename(  4,  3,  2)  /'                         '/
      DATA atm_filename(  4,  3,  3)  /'                         '/
      DATA atm_filename(  4,  3,  4)  /'                         '/
      DATA atm_filename(  4,  3,  5)  /'                         '/
      DATA atm_filename(  4,  3,  6)  /'                         '/
      DATA atm_filename(  4,  3,  7)  /'                         '/
      DATA atm_filename(  4,  4,  1)  /'                         '/
      DATA atm_filename(  4,  4,  2)  /'                         '/
      DATA atm_filename(  4,  4,  3)  /'                         '/
      DATA atm_filename(  4,  4,  4)  /'                         '/
      DATA atm_filename(  4,  4,  5)  /'                         '/
      DATA atm_filename(  4,  4,  6)  /'                         '/
      DATA atm_filename(  4,  4,  7)  /'                         '/
      DATA atm_filename(  4,  5,  1)  /'                         '/
      DATA atm_filename(  4,  5,  2)  /'                         '/
      DATA atm_filename(  4,  5,  3)  /'                         '/
      DATA atm_filename(  4,  5,  4)  /'                         '/
      DATA atm_filename(  4,  5,  5)  /'                         '/
      DATA atm_filename(  4,  5,  6)  /'                         '/
      DATA atm_filename(  4,  5,  7)  /'                         '/
      DATA atm_filename(  4,  6,  1)  /'t5777_4.4_-4.00_80.dat   '/
      DATA atm_filename(  4,  6,  2)  /'t5777_4.4_-3.00_80.dat   '/
      DATA atm_filename(  4,  6,  3)  /'t5777_4.4_-2.00_80.dat   '/
      DATA atm_filename(  4,  6,  4)  /'t5777_4.4_-1.00_80.dat   '/
      DATA atm_filename(  4,  6,  5)  /'t5777_4.4_-0.50_80.dat   '/
      DATA atm_filename(  4,  6,  6)  /'t5777_4.4_0.00_80.dat    '/
      DATA atm_filename(  4,  6,  7)  /'t5777_4.4_0.50_80.dat    '/
      DATA atm_filename(  4,  7,  1)  /'                         '/
      DATA atm_filename(  4,  7,  2)  /'                         '/
      DATA atm_filename(  4,  7,  3)  /'                         '/
      DATA atm_filename(  4,  7,  4)  /'                         '/
      DATA atm_filename(  4,  7,  5)  /'                         '/
      DATA atm_filename(  4,  7,  6)  /'                         '/
      DATA atm_filename(  4,  7,  7)  /'                         '/
      DATA atm_filename(  4,  8,  1)  /'                         '/
      DATA atm_filename(  4,  8,  2)  /'                         '/
      DATA atm_filename(  4,  8,  3)  /'                         '/
      DATA atm_filename(  4,  8,  4)  /'                         '/
      DATA atm_filename(  4,  8,  5)  /'                         '/
      DATA atm_filename(  4,  8,  6)  /'                         '/
      DATA atm_filename(  4,  8,  7)  /'                         '/
      DATA atm_filename(  5,  1,  1)  /'                         '/
      DATA atm_filename(  5,  1,  2)  /'                         '/
      DATA atm_filename(  5,  1,  3)  /'                         '/
      DATA atm_filename(  5,  1,  4)  /'                         '/
      DATA atm_filename(  5,  1,  5)  /'                         '/
      DATA atm_filename(  5,  1,  6)  /'                         '/
      DATA atm_filename(  5,  1,  7)  /'                         '/
      DATA atm_filename(  5,  2,  1)  /'                         '/
      DATA atm_filename(  5,  2,  2)  /'                         '/
      DATA atm_filename(  5,  2,  3)  /'                         '/
      DATA atm_filename(  5,  2,  4)  /'                         '/
      DATA atm_filename(  5,  2,  5)  /'                         '/
      DATA atm_filename(  5,  2,  6)  /'                         '/
      DATA atm_filename(  5,  2,  7)  /'                         '/
      DATA atm_filename(  5,  3,  1)  /'                         '/
      DATA atm_filename(  5,  3,  2)  /'                         '/
      DATA atm_filename(  5,  3,  3)  /'                         '/
      DATA atm_filename(  5,  3,  4)  /'                         '/
      DATA atm_filename(  5,  3,  5)  /'                         '/
      DATA atm_filename(  5,  3,  6)  /'                         '/
      DATA atm_filename(  5,  3,  7)  /'                         '/
      DATA atm_filename(  5,  4,  1)  /'t6000_3.5_-4.00_80.dat   '/
      DATA atm_filename(  5,  4,  2)  /'t6000_3.5_-3.00_80.dat   '/
      DATA atm_filename(  5,  4,  3)  /'t6000_3.5_-2.00_80.dat   '/
      DATA atm_filename(  5,  4,  4)  /'t6000_3.5_-1.00_80.dat   '/
      DATA atm_filename(  5,  4,  5)  /'t6000_3.5_-0.50_80.dat   '/
      DATA atm_filename(  5,  4,  6)  /'t6000_3.5_0.00_80.dat    '/
      DATA atm_filename(  5,  4,  7)  /'                         '/
      DATA atm_filename(  5,  5,  1)  /'t6000_4.0_-4.00_80.dat   '/
      DATA atm_filename(  5,  5,  2)  /'t6000_4.0_-3.00_80.dat   '/
      DATA atm_filename(  5,  5,  3)  /'t6000_4.0_-2.00_80.dat   '/
      DATA atm_filename(  5,  5,  4)  /'t6000_4.0_-1.00_80.dat   '/
      DATA atm_filename(  5,  5,  5)  /'t6000_4.0_-0.50_80.dat   '/
      DATA atm_filename(  5,  5,  6)  /'t6000_4.0_0.00_80.dat    '/
      DATA atm_filename(  5,  5,  7)  /'t6000_4.0_0.50_80.dat    '/
      DATA atm_filename(  5,  6,  1)  /'                         '/
      DATA atm_filename(  5,  6,  2)  /'                         '/
      DATA atm_filename(  5,  6,  3)  /'                         '/
      DATA atm_filename(  5,  6,  4)  /'                         '/
      DATA atm_filename(  5,  6,  5)  /'                         '/
      DATA atm_filename(  5,  6,  6)  /'                         '/
      DATA atm_filename(  5,  6,  7)  /'                         '/
      DATA atm_filename(  5,  7,  1)  /'t6000_4.5_-4.00_80.dat   '/
      DATA atm_filename(  5,  7,  2)  /'t6000_4.5_-3.00_80.dat   '/
      DATA atm_filename(  5,  7,  3)  /'t6000_4.5_-2.00_80.dat   '/
      DATA atm_filename(  5,  7,  4)  /'t6000_4.5_-1.00_80.dat   '/
      DATA atm_filename(  5,  7,  5)  /'t6000_4.5_-0.50_80.dat   '/
      DATA atm_filename(  5,  7,  6)  /'t6000_4.5_0.00_80.dat    '/
      DATA atm_filename(  5,  7,  7)  /'t6000_4.5_0.50_80.dat    '/
      DATA atm_filename(  5,  8,  1)  /'                         '/
      DATA atm_filename(  5,  8,  2)  /'                         '/
      DATA atm_filename(  5,  8,  3)  /'                         '/
      DATA atm_filename(  5,  8,  4)  /'                         '/
      DATA atm_filename(  5,  8,  5)  /'                         '/
      DATA atm_filename(  5,  8,  6)  /'                         '/
      DATA atm_filename(  5,  8,  7)  /'                         '/
      DATA atm_filename(  6,  1,  1)  /'                         '/
      DATA atm_filename(  6,  1,  2)  /'                         '/
      DATA atm_filename(  6,  1,  3)  /'                         '/
      DATA atm_filename(  6,  1,  4)  /'                         '/
      DATA atm_filename(  6,  1,  5)  /'                         '/
      DATA atm_filename(  6,  1,  6)  /'                         '/
      DATA atm_filename(  6,  1,  7)  /'                         '/
      DATA atm_filename(  6,  2,  1)  /'                         '/
      DATA atm_filename(  6,  2,  2)  /'                         '/
      DATA atm_filename(  6,  2,  3)  /'                         '/
      DATA atm_filename(  6,  2,  4)  /'                         '/
      DATA atm_filename(  6,  2,  5)  /'                         '/
      DATA atm_filename(  6,  2,  6)  /'                         '/
      DATA atm_filename(  6,  2,  7)  /'                         '/
      DATA atm_filename(  6,  3,  1)  /'                         '/
      DATA atm_filename(  6,  3,  2)  /'                         '/
      DATA atm_filename(  6,  3,  3)  /'                         '/
      DATA atm_filename(  6,  3,  4)  /'                         '/
      DATA atm_filename(  6,  3,  5)  /'                         '/
      DATA atm_filename(  6,  3,  6)  /'                         '/
      DATA atm_filename(  6,  3,  7)  /'                         '/
      DATA atm_filename(  6,  4,  1)  /'                         '/
      DATA atm_filename(  6,  4,  2)  /'                         '/
      DATA atm_filename(  6,  4,  3)  /'                         '/
      DATA atm_filename(  6,  4,  4)  /'                         '/
      DATA atm_filename(  6,  4,  5)  /'                         '/
      DATA atm_filename(  6,  4,  6)  /'                         '/
      DATA atm_filename(  6,  4,  7)  /'                         '/
      DATA atm_filename(  6,  5,  1)  /'t6500_4.0_-4.00_80.dat   '/
      DATA atm_filename(  6,  5,  2)  /'t6500_4.0_-3.00_80.dat   '/
      DATA atm_filename(  6,  5,  3)  /'t6500_4.0_-2.00_80.dat   '/
      DATA atm_filename(  6,  5,  4)  /'t6500_4.0_-1.00_80.dat   '/
      DATA atm_filename(  6,  5,  5)  /'t6500_4.0_-0.50_80.dat   '/
      DATA atm_filename(  6,  5,  6)  /'t6500_4.0_0.00_80.dat    '/
      DATA atm_filename(  6,  5,  7)  /'                         '/
      DATA atm_filename(  6,  6,  1)  /'                         '/
      DATA atm_filename(  6,  6,  2)  /'                         '/
      DATA atm_filename(  6,  6,  3)  /'                         '/
      DATA atm_filename(  6,  6,  4)  /'                         '/
      DATA atm_filename(  6,  6,  5)  /'                         '/
      DATA atm_filename(  6,  6,  6)  /'                         '/
      DATA atm_filename(  6,  6,  7)  /'                         '/
      DATA atm_filename(  6,  7,  1)  /'t6500_4.5_-4.00_80.dat   '/
      DATA atm_filename(  6,  7,  2)  /'t6500_4.5_-3.00_80.dat   '/
      DATA atm_filename(  6,  7,  3)  /'t6500_4.5_-2.00_80.dat   '/
      DATA atm_filename(  6,  7,  4)  /'t6500_4.5_-1.00_80.dat   '/
      DATA atm_filename(  6,  7,  5)  /'t6500_4.5_-0.50_80.dat   '/
      DATA atm_filename(  6,  7,  6)  /'t6500_4.5_0.00_80.dat    '/
      DATA atm_filename(  6,  7,  7)  /'t6500_4.5_0.50_80.dat    '/
      DATA atm_filename(  6,  8,  1)  /'                         '/
      DATA atm_filename(  6,  8,  2)  /'                         '/
      DATA atm_filename(  6,  8,  3)  /'                         '/
      DATA atm_filename(  6,  8,  4)  /'                         '/
      DATA atm_filename(  6,  8,  5)  /'                         '/
      DATA atm_filename(  6,  8,  6)  /'                         '/
      DATA atm_filename(  6,  8,  7)  /'                         '/
      DATA atm_filename(  7,  1,  1)  /'                         '/
      DATA atm_filename(  7,  1,  2)  /'                         '/
      DATA atm_filename(  7,  1,  3)  /'                         '/
      DATA atm_filename(  7,  1,  4)  /'                         '/
      DATA atm_filename(  7,  1,  5)  /'                         '/
      DATA atm_filename(  7,  1,  6)  /'                         '/
      DATA atm_filename(  7,  1,  7)  /'                         '/
      DATA atm_filename(  7,  2,  1)  /'                         '/
      DATA atm_filename(  7,  2,  2)  /'                         '/
      DATA atm_filename(  7,  2,  3)  /'                         '/
      DATA atm_filename(  7,  2,  4)  /'                         '/
      DATA atm_filename(  7,  2,  5)  /'                         '/
      DATA atm_filename(  7,  2,  6)  /'                         '/
      DATA atm_filename(  7,  2,  7)  /'                         '/
      DATA atm_filename(  7,  3,  1)  /'                         '/
      DATA atm_filename(  7,  3,  2)  /'                         '/
      DATA atm_filename(  7,  3,  3)  /'                         '/
      DATA atm_filename(  7,  3,  4)  /'                         '/
      DATA atm_filename(  7,  3,  5)  /'                         '/
      DATA atm_filename(  7,  3,  6)  /'                         '/
      DATA atm_filename(  7,  3,  7)  /'                         '/
      DATA atm_filename(  7,  4,  1)  /'                         '/
      DATA atm_filename(  7,  4,  2)  /'                         '/
      DATA atm_filename(  7,  4,  3)  /'                         '/
      DATA atm_filename(  7,  4,  4)  /'                         '/
      DATA atm_filename(  7,  4,  5)  /'                         '/
      DATA atm_filename(  7,  4,  6)  /'                         '/
      DATA atm_filename(  7,  4,  7)  /'                         '/
      DATA atm_filename(  7,  5,  1)  /'                         '/
      DATA atm_filename(  7,  5,  2)  /'                         '/
      DATA atm_filename(  7,  5,  3)  /'                         '/
      DATA atm_filename(  7,  5,  4)  /'                         '/
      DATA atm_filename(  7,  5,  5)  /'                         '/
      DATA atm_filename(  7,  5,  6)  /'                         '/
      DATA atm_filename(  7,  5,  7)  /'                         '/
      DATA atm_filename(  7,  6,  1)  /'                         '/
      DATA atm_filename(  7,  6,  2)  /'                         '/
      DATA atm_filename(  7,  6,  3)  /'                         '/
      DATA atm_filename(  7,  6,  4)  /'                         '/
      DATA atm_filename(  7,  6,  5)  /'                         '/
      DATA atm_filename(  7,  6,  6)  /'                         '/
      DATA atm_filename(  7,  6,  7)  /'                         '/
      DATA atm_filename(  7,  7,  1)  /'                         '/
      DATA atm_filename(  7,  7,  2)  /'                         '/
      DATA atm_filename(  7,  7,  3)  /'t7000_4.5_-2.00_80.dat   '/
      DATA atm_filename(  7,  7,  4)  /'t7000_4.5_-1.00_80.dat   '/
      DATA atm_filename(  7,  7,  5)  /'t7000_4.5_-0.50_80.dat   '/
      DATA atm_filename(  7,  7,  6)  /'t7000_4.5_0.00_80.dat    '/
      DATA atm_filename(  7,  7,  7)  /'                         '/
      DATA atm_filename(  7,  8,  1)  /'                         '/
      DATA atm_filename(  7,  8,  2)  /'                         '/
      DATA atm_filename(  7,  8,  3)  /'                         '/
      DATA atm_filename(  7,  8,  4)  /'                         '/
      DATA atm_filename(  7,  8,  5)  /'                         '/
      DATA atm_filename(  7,  8,  6)  /'                         '/
      DATA atm_filename(  7,  8,  7)  /'                         '/