!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up everything you need for running simulations 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initial
USE DEFINITION
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

WRITE(*,*) 'Build grid variables'
CALL GETGRID_NM
IF(DM_flag) THEN
	CALL GETGRID_DM
END IF
WRITE(*,*) 'Done building grid variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for setting initial conditions !

! Set initial grid size !
delta1 = dx1
delta2 = dx2

! Set initial conservative/primitive variables !
cons2 = 0.0D0
prim2 = 0.0D0
IF(RUNDM_flag) THEN
    cons1 = 0.0D0
    prim1 = 0.0D0
END IF

! Build initial models !
WRITE(*,*) 'Now build the initial model'

! Do it case by case !
IF(hydro_test) THEN 
   	IF(n_dim == 1) THEN
        CALL Riemann_1d
   	ELSEIF(n_dim == 2) THEN
        CALL Riemann_2d
    END IF
END IF

WRITE(*,*) 'Finished building initial model'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Build conservative variables'
CALL FROMRVETOU
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL
WRITE(*,*) 'Do Update'
CALL UPDATE (1) 
WRITE(*,*) 'Done initial update'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish initial...'
WRITE(*,*)

END SUBROUTINE