!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the accretion torus problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field properties !

! Magnetisation parameter !
REAL*8, PARAMETER :: p_beta = 100.0

! magnetic field normalization mode !
LOGICAL, PARAMETER :: normalize_by_vol = .false.
LOGICAL, PARAMETER :: normalize_by_minbeta = .true.

! use refined gird? !
LOGICAL, PARAMETER :: refined_grid = .false.

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
REAL*8, PARAMETER :: s_0 = 3.0d0

! (equatorial) radius where the density is at maximum !
REAL*8, PARAMETER :: s_max = 4.7d0

! outer (equatorial) radius of the torus !
REAL*8, PARAMETER :: s_1 = 8.0d0

! Outer domain of the box in the r-direction to contain finest resolution !
REAL*8, PARAMETER :: x_fine = 11.5D0

! Inner domain of the box in the theta-direction to contain finest resolution !
REAL*8, PARAMETER :: y_fine = pi/4.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Controlling density/magnetic field floor

! minimum density (fraction) to define the last contour of the vector potential !
REAL*8, PARAMETER :: rho_cut = 5.0d-1

! atmospheric density !
REAL*8, PARAMETER :: rho_fac = 1.0d-4

! maximum density !
REAL*8, PARAMETER :: rho_max = 1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!