      INTEGER*4	    NMOLMAX		       ! max. nr of molecules
      INTEGER*4	    MAX_MOL_LEV		       ! max. nr of level per molecule
      CHARACTER*80  mol_data_fn		       ! filename for molecular data

      PARAMETER (MAX_MOL_LEV = 17)
      PARAMETER (mol_data_fn = 'molecule.asc')

      INTEGER*4	    elm_idx(2*NMOLMAX)	         ! atomic_nrs of elements involed in molecule formation
      INTEGER*4	    nmol		         ! Nr of molecules taken into account
      INTEGER*4	    nelem		         ! Nr of different elements taken into account for molecule formation

      STRUCTURE /mol_data_typ/
	 CHARACTER*4   nam			 ! molecule name
	 REAL*4	       tmass			 ! total mass
	 REAL*4	       rmass			 ! reduced mass
	 INTEGER*4     i1			 ! 1. atomic nr
	 INTEGER*4     i2			 ! 2. atomic nr
	 REAL*4	       d0			 ! dissociation constant
	 REAL*4	       sym			 ! symetry parameter
	 INTEGER*4     nlev			 ! nr of levels taken into account
	 CHARACTER*6   level(MAX_MOL_LEV)	 ! designation
	 INTEGER*4     case(MAX_MOL_LEV)	 ! coupling case
	 REAL*4	       te(MAX_MOL_LEV)		 ! electronic excitation energy in [1/cm]
	 REAL*4	       ge(MAX_MOL_LEV)		 ! electronic multiplicity
	 REAL*4	       we(MAX_MOL_LEV)		 ! vibrational constant in [1/cm]
	 REAL*4	       wexe(MAX_MOL_LEV)	 ! 1. vibrational anharmonicity constant in [1/cm]
	 REAL*4	       weye(MAX_MOL_LEV)	 ! 2. vibrational anharmonicity constant in [1/cm]
	 REAL*4	       weze(MAX_MOL_LEV)	 ! 3. vibrational anharmonicity constant in [1/cm]
	 REAL*4	       be(MAX_MOL_LEV)		 ! rotational constant in [1/cm]
	 REAL*4	       alpha_e(MAX_MOL_LEV)	 ! 1. correction of rotational constant in [1/cm]
	 REAL*4	       gamma_e(MAX_MOL_LEV)	 ! 2. correction of rotational constant in [1/cm]
	 REAL*4	       delta_e(MAX_MOL_LEV)	 ! 3. correction of rotational constant in [1/cm]
	 REAL*4	       epsilon_e(MAX_MOL_LEV)	 ! 4. correction of rotational constant in [1/cm]
	 REAL*4	       zeta_e(MAX_MOL_LEV)	 ! 5. correction of rotational constant in [1/cm]
	 REAL*4	       alpha(MAX_MOL_LEV)	 ! 1
	 REAL*4	       gamma(MAX_MOL_LEV)	 ! 2
	 REAL*4	       delta(MAX_MOL_LEV)	 ! 3
	 REAL*4	       epsilon(MAX_MOL_LEV)	 ! 4
	 REAL*4	       zeta(MAX_MOL_LEV)	 ! 5
      END STRUCTURE


      RECORD /mol_data_typ/ molec(NMOLMAX)
