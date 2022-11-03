!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up everything you need for running simulations 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initial
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First build up all database, EOS table and arrays !
WRITE(*,*) 'We shall now setup everything for the simulations'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for building arrays !

WRITE(*,*) 'Build hydro equations'
CALL SETUP_EQN
WRITE(*,*) 'End building hydro equations'
WRITE(*,*)

WRITE(*,*) 'Build arrays'
CALL BUILD_HYDRO
WRITE(*,*) 'End building arrays'
WRITE(*,*)

WRITE(*,*) 'Build grid variables'
CALL GETGRIDNM
IF(DM_flag == 1) THEN
	CALL GETGRIDDM
END IF
WRITE(*,*) 'Done building grid variables'
WRITE(*,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for setting initial conditions !

! Build initial models !
WRITE(*,*) 'Now build the initial model'

! Do it case by case !
!IF(hydro_test) THEN 
!    IF(dimension_flag == 1) THEN
!        CALL Riemann_1d
!    END IF
!END IF

WRITE(*,*) 'Finished building initial model'
WRITE(*,*)

WRITE(*,*) 'Build atmospheric primitive variables'
!CALL GET_ATM
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Build conservative variables'
!CALL FROMRVETOU
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL
WRITE(*,*) 'Do Update'
!CALL UPDATE (1) 
WRITE(*,*) 'Done initial update'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish initial...'
WRITE(*,*)

END SUBROUTINE