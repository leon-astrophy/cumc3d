!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This module contains all the arrays and variables that are neccessary to run the hydro
! simulations of either pure hydro, discs, or self gravitating fluid objects like stars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary flag for variables !

INTEGER :: bfac_x(100)
INTEGER :: bfac_y(100)
INTEGER :: bfac_z(100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Minimum/Maximum Eq. number for DM
INTEGER :: imin1
INTEGER :: imax1

! Identifiers for the DM variables !

! DM density !
INTEGER :: irho1

! DM x-velocity
INTEGER :: ivel1_x

! DM y-velocity
INTEGER :: ivel1_y

! DM z-velocity
INTEGER :: ivel1_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Minimum/Maximum Eq. number for NM
INTEGER :: imin2
INTEGER :: imax2

! Identifiers for the NM variables !

! NM density !
INTEGER :: irho2

! NM x-velocity
INTEGER :: ivel2_x

! NM y-velocity
INTEGER :: ivel2_y

! NM z-velocity
INTEGER :: ivel2_z

! NM total energy density
INTEGER :: itau2

! NM internal energy density
INTEGER :: ieps2

! NM electron fractions
INTEGER :: iye2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
! Grid variables !

! Grid coordinates for DM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: x1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: y1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: z1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: dr1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: xF1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: yF1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: zF1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: cos1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: sin1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: vol1

! R and Z coordinate of the grid for NM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: x2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: y2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: z2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: dr2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: xF2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: yF2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: zF2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: cos2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: sin2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: vol2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrodynamical variables

! DM atmospheric primitive variables !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: prim1_a

! DM primitive and conservative variables !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: prim1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: cons1

! DM pressure, speed of sound, and pressure derivatives !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: p1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: cs1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: dpdrho1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: dpdeps1

! DM multipole moments !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: qpole1

! DM total gravitational potentials !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi1

! Potentials from the DM contribution in the DM grid !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi1_dm

! Potentials from the NM contribution in the DM grid !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi1_nm

! DM potential derivatives !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: dphi_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are the hydro set for NM sector

! NM atmospheric primitive variables !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: prim2_a

! NM primitive and conservative variables !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: prim2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: cons2

! NM speed of sound, and pressure derivatives !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: cs2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: dpdrho2
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: dpdeps2

! NM temperature !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: temp2

! NM internal energy !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: epsilon2

! NM multipole moments !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: qpole2

! NM total gravitational potentials !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi2

! Potentials from the NM contribution in the NM grid !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi2_dm

! Potentials from the NM contribution in the NM grid !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: phi2_nm

! NM potential derivatives !
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: dphi_2

! Pressure gradients
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: dp_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for time-evolution !

! For RK-Time evolution 
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: u_old1, u2_dm, u3_dm, l3_dm, l1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:,:) :: u_old2, u2_nm, u3_nm, l3_nm, l2

! The auxillary array for the flux term, DM
REAL (DP), allocatable, DIMENSION (:,:,:,:) :: sa1, sb1
REAL (DP), allocatable, DIMENSION (:,:,:,:,:) :: flux_1
REAL (DP), allocatable, DIMENSION (:,:,:,:,:) :: dflux_1

! Flux arrays for NM !
REAL (DP), allocatable, DIMENSION (:,:,:,:) :: sa2, sb2
REAL (DP), allocatable, DIMENSION (:,:,:,:,:) :: flux_2
REAL (DP), allocatable, DIMENSION (:,:,:,:,:) :: dflux_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real/integer variable scalars !

! The number of highest grid along three axis !
INTEGER :: nx_part_1, ny_part_1, nz_part_1
INTEGER :: nx_part_2, ny_part_2, nz_part_2
							
! The number of lowest grid along three axis !
INTEGER :: nx_min_1, ny_min_1, nz_min_1
INTEGER :: nx_min_2, ny_min_2, nz_min_2

! Time step !
REAL (DP) :: dt

! Total grid size !
REAL (DP) :: delta1	
REAL (DP) :: delta2

! Polytropic index !
REAL (DP) :: kgas1, ggas1
REAL (DP) :: kgas2, ggas2

! Atmospheric temperature !
REAL (DP) :: temp2_a

! Atmospheric epsilon !
REAL (DP) :: eps1_a, eps2_a

! Masses ! 
REAL (DP) :: mass1, mass2

! Individual energyies !
REAL (DP) :: energy1_kin, energy1_int, energy1_pot
REAL (DP) :: energy2_kin, energy2_int, energy2_pot

! Some logical count flag !
INTEGER :: potential_flag

! Global simulation time !
REAL (DP) :: global_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DEFINITION