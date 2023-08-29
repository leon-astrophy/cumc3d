!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up everything you need for running simulations 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initial_model
USE DEFINITION
USE MHD_MODULE
USE PPM_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set minimum and maximum domain !
nx_min_2 = 1
nx_part_2 = nx_2
ny_min_2 = 1
ny_part_2 = ny_2
nz_min_2 = 1
nz_part_2 = nz_2

! Set minimum and maximum domain !
nx_min_1 = 1
nx_part_1 = nx_1
ny_min_1 = 1
ny_part_1 = ny_1
nz_min_1 = 1
nz_part_1 = nz_1

! First build up all database, EOS table and arrays !
WRITE(*,*) 'We shall now setup everything for the simulations'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for building arrays !

WRITE(*,*) 'Build hydro equations'
CALL SETUP_EQN
WRITE(*,*) 'End building hydro equations'
WRITE(*,*)

WRITE(*,*) 'Build hydro arrays'
CALL BUILD_HYDRO
WRITE(*,*) 'End building hydro arrays'
WRITE(*,*)

WRITE(*,*) 'Build riemann solver variables'
CALL BUILDRIEMANN
WRITE(*,*) 'End building riemann solver variables'
WRITE(*,*)

! MHD !
WRITE(*,*) 'Build MHD variables'
CALL buildMHD
WRITE(*,*) 'End building MHD variables'
WRITE(*,*)

WRITE(*,*) 'Build grid variables'
CALL GETGRID_NM
WRITE(*,*) 'Done building grid variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Setup PPM weight 
CALL PPM_WEIGHT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for setting initial conditions !

! Set initial conservative/primitive variables !
cons2 = 0.0D0
prim2 = 0.0D0

! Build initial models !
WRITE(*,*) 'Now build the initial model'
CALL GET_MODEL
WRITE(*,*) 'Finished building initial model'
WRITE(*,*)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do initial updates after setting up the inital model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initial_update
USE DEFINITION
USE MHD_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare everything ... !

! Check density !
IF (checkrho_flag) THEN
	CALL CHECKRHO
END IF

WRITE(*,*) 'Build conservative variables'
CALL FROMRVETOU
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

! set boundary conditions !
call BOUNDARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL
WRITE(*,*) 'Do Update'
CALL UPDATE (0)
WRITE(*,*) 'Done initial update'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish initial...'
WRITE(*,*)

END SUBROUTINE