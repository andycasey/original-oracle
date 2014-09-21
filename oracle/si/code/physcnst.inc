*
* ----------------------------------------------------------------------
* >>>
* FILENAME  : PHYSICAL_CONSTANTS.INC
* DIRECTORY :
* PURPOSE   : Declaration of physical and mathematical constants of interest
* MODULE    : FORTRAN-INCLUDE FILE
* ADAPTED TO: LINE_FORM.FOR,....
* AUTOR     : JR
* DATE      : 03.09.92
* REMARK    : units are explizit declared
* <<<
* ----------------------------------------------------------------------
                                              

      REAL*8          PI, E, H, H_EV, H_PI, M_E, KB, PI2, LOG_E,
     *                C, HC_DIV_K, SQRT_PI, HCC_2, E_SI, LN10,
     *                FAC_A, FAC_B, FAC_C,FAC_D, DIV_4PIC, R_2,
     *                LOG_K,LN10_DIV_3,H_E2,PI_2,PI_4,DIV_PIC,
     *                PI2_4_DIV_K,H_E2_DIV_K,RYD,KB23,
     *                LN10_DIV_2,LN10_3DIV8,LN10_2,U,M_H,A0,
     *                DIV_K, R,C3FAC,HC_E_SI,
     *                LOG2,
     *                EV_K,SIGE
                                
      PARAMETER       (C       =  2.99792456D 10)  ! vac. velocity of light [cm/s]
      PARAMETER       (E       =  4.80324235D-10)  ! electron charge   [esE]
      PARAMETER       (KB      =  1.380662  D-16)  ! k-Boltzmann       [erg/K]
      PARAMETER       (KB23    =  2.6713135 D-11)  ! k-Boltzmann^0.667 [erg/K]^(2/3)
      PARAMETER       (R       =  8.31441   D 07)  ! Gas-Constant      [erg/molK]
      PARAMETER       (H       =  6.626176  D-27)  ! planck's constant [erg s]
      PARAMETER       (PI      =  3.14159265359D0) ! Pi
      PARAMETER       (SQRT_PI =  1.77245385091D0) ! sqrt(pi)
      PARAMETER       (U       =  1.66057   D-24)  ! atomic mass unit  [g]
      PARAMETER       (M_E     =  9.109534  D-28)  ! electron mass     [g]
      PARAMETER       (LOG_E   =  0.434294482D0 )  ! log10(exp(1))
      PARAMETER       (LN10    =  2.30258509299D0) ! ln(10)
      PARAMETER       (LOG2    =  3.01029995664D -1)  ! log10(2.)
      PARAMETER       (LOG_K   =-15.8599126281D0)  ! log10(k)          [erg/K]
      PARAMETER       (SIGE    =  6.653     D-25)  !Thomson scattering coefficient
      PARAMETER       (C3FAC   =  4.9303    D-04)  ! Ali and Griem (1965,1966) [cgs]
                                                   ! + dipole-quadrupole interaction
      PARAMETER       (H_E2    = 1.6340556  D-11)  ! Energy of n=2 H-level     [erg]
      PARAMETER       (A0      = 5.2917715  D-09)  ! Bohr atom radius [cm]
      PARAMETER       (RYD     = 1.09737312 D 05)  ! Rydberg [1/cm]

C ---- derived quantities -----------

      PARAMETER       (M_H     = 1.008D0*U      )  ! mass of Hydrogen   [g]
      PARAMETER       (E_SI    = E*10./C        )  ! electron charge   [C] (SI-Unit)
      PARAMETER       (EV_K    = E_SI/KB*1.D7   )  ! eV --> K   [K] (SI-Unit)
      PARAMETER       (HC_E_SI = H*C*1.D-7/E_SI )  ! hc/e  in [eVcm]
      PARAMETER       (PI_2    = 2.D0*PI        )  !
      PARAMETER       (PI_4    = 4.D0*PI        )  !
      PARAMETER       (DIV_4PIC= 1.D0/(4.D0*PI*C)) !
      PARAMETER       (DIV_PIC = 1.D0/(PI*C)    )  !
      PARAMETER       (DIV_K   = 1.D0/KB        )  ! 1/k               [K/erg]
      PARAMETER       (PI2_4_DIV_K = 4.D0*PI*PI/KB)! 4pi^2/k           [K/erg]
      PARAMETER       (FAC_A   = PI*E*E/(M_E*C*C)) ! pi*e^2)/(m_e*c^2) [cgs-units]
      PARAMETER       (FAC_B   = FAC_A/(SQRT_PI*KB)) ! sqrt(pi)*e^2)/(k*m_e*c^2) [cgs-units]
      PARAMETER       (FAC_C   = FAC_A/ KB      )  ! pi*e^2)/(k*m_e*c^2) [cgs-units]
      PARAMETER       (FAC_D   = 2.5982699*E/KB23)  ! pi*e^2)/(k*m_e*c^2) [cgs-units]
      PARAMETER       (HCC_2   = 2.D0*H*C*C     )  !  2 * h * c * c    [cgs]
      PARAMETER       (HC_DIV_K= H*C/KB         )  ! (h*c)/k           [cm*K]
      PARAMETER       (H_E2_DIV_K = H_E2/KB     )  ! Excitation energy in [K]
      PARAMETER       (LN10_DIV_2  = LN10 / 2.D0)  ! ln(10) / 2
      PARAMETER       (LN10_2      = LN10 * 2.D0)  ! ln(10) * 2
      PARAMETER       (LN10_DIV_3  = LN10 / 3.D0)  ! ln(10) / 3     
      PARAMETER       (LN10_3DIV8 =LN10*3.D0/8.D0) ! ln(10) * 3/8  
      PARAMETER       (R_2        = 2.D0*R      )  ! 2* Gasconstant    [erg/molK]
