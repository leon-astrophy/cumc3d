!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the tidal disruption event problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! adiabatic index !
REAL*8, PARAMETER :: gamma = (5.0d0/3.0d0)

! power law for density profiles !
REAL*8, PARAMETER :: rho_grad = 2.5d0

! power law for angular velocity !
REAL*8, PARAMETER :: q_grad = 2.0d0

! schwarzschild radius !
REAL*8, PARAMETER :: r_sh = 1.0

! atmospheric density !
REAL*8, PARAMETER :: rho_fac = 1.0d-4

! reference radius !
REAL*8, PARAMETER :: r_0 = 1.5d0

! maximum density !
REAL*8, PARAMETER :: rho_max = 1.0d0

! maximum pressure !
REAL*8, PARAMETER :: p_max = 1.0d0

! angular velocity !
REAL*8, PARAMETER :: omega_max = 1.0d0

! minimum magnetisation !
REAL*8, PARAMETER :: p_beta = 100.0d0
