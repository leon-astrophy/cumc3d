!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the collapsar problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit constants !

! Physical constants to be as one !
REAL*8, PARAMETER :: gconst = 6.67430D-8
REAL*8, PARAMETER :: clight = 2.99792458D10
REAL*8, PARAMETER :: solar = 1.98847D33

! Here, mu_0 is not in cgs !
REAL*8, PARAMETER :: mu_0 = 1.25663706212D-6

! Solar Radius !
REAL*8, PARAMETER :: rsolar = 6.96342D10

! 1 GeV !
REAL*8, PARAMETER :: GeV2gram = 1.78266191D-24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversions !

! Conversion between units !
REAL*8, PARAMETER :: lencgs2code = (clight**2)/(solar*gconst)
REAL*8, PARAMETER :: masscgs2code = (1.0D0/solar)
REAL*8, PARAMETER :: tcgs2code = (clight**3)/(solar*gconst)

! Derived conversion !
REAL*8, PARAMETER :: rhocgs2code = (masscgs2code/lencgs2code**3)
REAL*8, PARAMETER :: taucgs2code = (masscgs2code*lencgs2code**2/tcgs2code**2)
REAL*8, PARAMETER :: h_bar = (1.054571817D-27)*(lencgs2code**2*masscgs2code/tcgs2code)

! Current conversion !
REAL*8, PARAMETER :: amp2code = (mu_0*1.0D5*masscgs2code*lencgs2code)**(0.5D0)*tcgs2code

! Magnetic field !
REAL*8, PARAMETER :: gauss2code = 1.0D-1*masscgs2code/amp2code/tcgs2code**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the physical setup !
 
! which rotational profile? 1 = rigid, 2 = differential !
INTEGER, PARAMETER :: rotation_rule = 2

! Differential rotations, position of maximum angular velocity (in cgs) !
REAL*8, PARAMETER :: s_0 = (5.0D3*1.0D5)*lencgs2code

! Differential rotations, core angular velocity (in cgs) !
REAL*8, PARAMETER :: w_0 = 0.5D0/tcgs2code

! Maximum magnetic field strength (in gauss) !
REAL*8, PARAMETER :: b_0 = 0.0d0 !1.0D10*gauss2code

! Limiting rotation velocity for the progenitor
REAL*8, PARAMETER :: v_crit = 0.01D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 10

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 100000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
