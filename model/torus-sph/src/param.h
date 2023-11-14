!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the accretion torus problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field properties !

! Magnetisation parameter !
REAL*8, PARAMETER :: p_beta = 350.0

! magnetic field normalization mode !
LOGICAL, PARAMETER :: normalize_by_vol = .false.
LOGICAL, PARAMETER :: normalize_by_minbeta = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Torus physical properties !

! schwarzschild radius !
REAL*8, PARAMETER :: r_sh = 1.0

! adiabatic index !
REAL*8, PARAMETER :: gamma = (5.0d0/3.0d0)

! angular velocity gradient !
REAL*8, PARAMETER :: q_grad = 2.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Torus size, and computational grid !

! inner (equatorial) radius of the torus !
REAL*8, PARAMETER :: s_in = 3.0d0

! (equatorial) radius where the density is at maximum !
REAL*8, PARAMETER :: s_max = 4.7d0

! Custom grid? !
LOGICAL, PARAMETER :: refine_grid = .true. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Corona !

! Want a corona? !
LOGICAL, PARAMETER :: corona = .false. 

! ideal gas m/kt !
REAL*8, PARAMETER :: a_corn = 1.0d0

! Density scale 
REAL*8, PARAMETER :: a_eta = 1.0d-4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Controlling density/magnetic field floor

! minimum density (fraction) to define the last contour of the vector potential !
REAL*8, PARAMETER :: rho_cut = 5.0d-1

! atmospheric density !
REAL*8, PARAMETER :: rho_fac = 1.0d-4

! maximum density !
REAL*8, PARAMETER :: rho_max = 1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
