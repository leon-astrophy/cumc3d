!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the collapsar problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit constants !

! Speed of light in cgs !
REAL*8, PARAMETER :: clight = 2.99792458D10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! floor values for alven speed (in speed of light) !
REAL*8 :: max_alv = 5.0d0*clight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the physical setup !
 
! which rotational profile? 1 = rigid, 2 = differential !
INTEGER, PARAMETER :: rotation_rule = 2

! Differential rotations, position of maximum angular velocity (in cgs) !
REAL*8 :: s_0 = (5.0D3)*(1.0D5)

! Differential rotations, core angular velocity (in cgs) !
REAL*8 :: w_0 = 10.0D0

! Maximum magnetic field strength (in gauss) !
REAL*8 :: b_0 = 1.0D10

! Limiting rotation velocity for the progenitor (in speed of light) !
REAL*8 :: v_crit = 0.01D0*clight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 10

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 100000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
