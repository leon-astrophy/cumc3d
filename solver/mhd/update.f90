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
IMPLICIT NONE

! Signal for finding potential
INTEGER, INTENT (IN) :: p_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find pressure
CALL FINDPRESSURE

! Speed of sound !
CALL SOUNDSPEED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE