!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine prepare the data necessary for constructing
! the flux for the spatial discretization.
! It takes input/output of the U array and the 
! mode p which decides whether or not to update
! the gravitational potential
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE UPDATE (p_in)
USE DEFINITION
USE LEVELSET_MODULE
USE NUCEOS_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Signal for finding potential
INTEGER, INTENT (IN) :: p_in

! Dummy variable
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find temperature !
IF(helmeos_flag == 1) THEN
	CALL FINDHELMTEMP
ELSEIF(nuceos_flag == 1) THEN
	call findnuctemp
END IF

! Find pressure
CALL FINDPRESSURE

! For dual energy
IF(dual_energy == 1) THEN
	CALL FINDRHOE
END IF

! Speed of sound !
CALL SOUNDSPEED

! Update turbulence transport term
IF(turb_flag == 1) THEN
	CALL FINDTURBULENCE
END IF

! Update the potential
IF (p_in == 1 .OR. (testmodel_flag == 5 .AND. initmodel_flag == 0)) THEN
	CALL FINDPOTENTIAL
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE