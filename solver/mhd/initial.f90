!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up everything you need for running simulations 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL_MODEL
USE DEFINITION
USE MHD_MODULE
USE TVD_MODULE
USE WENO_MODULE
USE PPMC_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First build up all database, EOS table and arrays !
WRITE(*,*) 'We shall now setup everything for the simulations'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for building arrays !

WRITE(*,*) 'Build hydro equations'
CALL SETUP_EQNS
WRITE(*,*) 'End building hydro equations'
WRITE(*,*)

WRITE(*,*) 'Build hydro arrays'
CALL BUILD_HYDRO
WRITE(*,*) 'End building hydro arrays'
WRITE(*,*)

WRITE(*,*) 'Build riemann solver variables'
CALL BUILD_RIEMANN
WRITE(*,*) 'End building riemann solver variables'
WRITE(*,*)

! MHD !
WRITE(*,*) 'Build MHD variables'
CALL BUILD_MHD
WRITE(*,*) 'End building MHD variables'
WRITE(*,*)

WRITE(*,*) 'Build grid variables'
CALL GETGRID
WRITE(*,*) 'Done building grid variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set reconstruction weight !

WRITE(*,*) 'Build reconstruction weight'
IF(weno_flag) THEN
  CALL WENO_WEIGHT
ELSEIF(ppmc_flag) THEN
  CALL PPM_WEIGHT
ELSE
  CALL TVD_WEIGHT
END IF
WRITE(*,*) 'Done building reconstruction weight'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for setting initial conditions !

! Set initial conservative/primitive variables !
cons = 0.0D0
prim = 0.0D0

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

SUBROUTINE INITIAL_UPDATE
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare everything ... !

! Check density !
CALL CUSTOM_CHECKRHO

WRITE(*,*) 'Build conservative variables'
CALL FROMRVETOU
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

! set boundary conditions !
CALL BOUNDARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL
WRITE(*,*) 'Do update'
CALL UPDATE (0)
WRITE(*,*) 'Done initial update'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish initial...'
WRITE(*,*)

END SUBROUTINE