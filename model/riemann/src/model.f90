!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
IMPLICIT NONE

! Do it case by case !
IF(n_dim == 1) THEN
  CALL Riemann_1d
ELSEIF(n_dim == 2) THEN
  CALL Riemann_2d
ELSEIF(n_dim == 3) THEN
  CALL Riemann_3d
END IF

END SUBROUTINE