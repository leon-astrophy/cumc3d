!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the maximum effective speed along a constant x
! Written by Leung Shing Chi in 2016
! The effective speed is obtained by solving the determinant of the 
! Jacobian of the flux vector written in terms of the primitive variables
! Modification of this subroutine is necessary when we modify the 
! Euler equation to include other physics, such as B-field
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALPHASPLIT_X (alpha1, alpha2, k_in, l_in, type_in)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type_in

! Input: The row number
INTEGER, INTENT(IN) :: k_in, l_in

! Output: The effective speed
REAL (DP), INTENT (OUT) :: alpha1(imin1:imax1)
REAL (DP), INTENT (OUT) :: alpha2(imin2:imax2)

! The candidate effective speeds for non-MHD
real (DP), DIMENSION(-2:nx_1+3) :: lambda1
real (DP), DIMENSION(-2:nx_2+3) :: lambda2

! Include DM component
IF(type_in == 1) THEN
   lambda1(:) = ABS(prim1(:,k_in,l_in,ivel1_x)) + DSQRT (cs1(:,k_in,l_in))

   ! Send the DM output
   alpha1(imin1:imax1) = MAXVAL(lambda1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Include NM component
IF(type_in == 2) THEN         
   lambda2(:) = ABS(prim2(:,k_in,l_in,ivel2_x)) + DSQRT (cs2(:,k_in,l_in))

   ! Send the DM output
   alpha2(imin2:imax2) = MAXVAL(lambda2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the maximum effective speed along a constant y
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALPHASPLIT_Y (alpha1, alpha2, j_in, l_in, type_in)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type_in

! Input: The row number
INTEGER, INTENT(IN) :: j_in, l_in

! Output: The effective speed
REAL (DP), INTENT (OUT) :: alpha1(imin1:imax1)
REAL (DP), INTENT (OUT) :: alpha2(imin2:imax2)

! The candidate effective speeds for non-MHD
real (DP), DIMENSION(-2:ny_1+3) :: lambda1
real (DP), DIMENSION(-2:ny_2+3) :: lambda2

! Include DM component
IF(type_in == 1) THEN
   lambda1(:) = ABS(prim1(j_in,:,l_in,ivel1_y)) + DSQRT (cs1(j_in,:,l_in))

   ! Send the DM output
   alpha1(imin1:imax1) = MAXVAL(lambda1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Include NM component
IF(type_in == 2) THEN         
   lambda2(:) = ABS(prim2(j_in,:,l_in,ivel2_y)) + DSQRT (cs2(j_in,:,l_in))

   ! Send the DM output
   alpha2(imin2:imax2) = MAXVAL(lambda2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the maximum effective speed along a constant y
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALPHASPLIT_Z (alpha1, alpha2, j_in, k_in, type_in)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type_in

! Input: The row number
INTEGER, INTENT(IN) :: j_in, k_in

! Output: The effective speed
REAL (DP), INTENT (OUT) :: alpha1(imin1:imax1)
REAL (DP), INTENT (OUT) :: alpha2(imin2:imax2)

! The candidate effective speeds for non-MHD
real (DP), DIMENSION(-2:nz_1+3) :: lambda1
real (DP), DIMENSION(-2:nz_2+3) :: lambda2

! Include DM component
IF(type_in == 1) THEN
   lambda1(:) = ABS(prim1(j_in,k_in,:,ivel1_z)) + DSQRT (cs1(j_in,k_in,:))

   ! Send the DM output
   alpha1(imin1:imax1) = MAXVAL(lambda1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Include NM component
IF(type_in == 2) THEN         
   lambda2(:) = ABS(prim2(j_in,k_in,:,ivel2_z)) + DSQRT (cs2(j_in,k_in,:))

   ! Send the DM output
   alpha2(imin2:imax2) = MAXVAL(lambda2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE to find the maximum effective speed !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDALPHA
USE OMP_LIB
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Dummy variables !
REAL (DP), DIMENSION(imin1:imax1) :: dummy1
REAL (DP), DIMENSION(imin2:imax2) :: dummy2

! Do for DM !
!$OMP PARALLEL SHARED(alpha1_x, alpha1_y, alpha1_z, alpha2_x, alpha2_y, alpha2_z)
IF(DM_flag) THEN
   !$OMP DO SIMD COLLAPSE(2) SCHEDULE(STATIC)
   DO k = ny_min_1, ny_part_1
      DO l = nz_min_1, nz_part_1
		   CALL ALPHASPLIT_X (alpha1_x(:,k,l), dummy2(:), k, l, dm_f)
	   END DO
   END DO
   !$OMP END DO
   IF(n_dim > 1) THEN
      !$OMP DO SIMD COLLAPSE(2) SCHEDULE(STATIC)
      DO l = nz_min_1, nz_part_1
         DO j = nx_min_1, nx_part_1
		      CALL ALPHASPLIT_Y (alpha1_y(:,j,l), dummy2(:), j, l, dm_f)
	      END DO
      END DO
      !$OMP END DO
   END IF
   IF(n_dim > 2) THEN
      !$OMP DO SIMD COLLAPSE(2) SCHEDULE(STATIC)
      DO k = ny_min_1, ny_part_1
	      DO j = nx_min_1, nx_part_1
		      CALL ALPHASPLIT_Z (alpha1_z(:,j,k), dummy2(:), j, k, dm_f)
	      END DO
      END DO
      !$OMP END DO
   END IF
END IF

! Do for NM !
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
   DO k = ny_min_2, ny_part_2
		CALL ALPHASPLIT_X (dummy1(:), alpha2_x(:,k,l), k, l, nm_f)
   END DO
END DO
!$OMP END DO
IF(n_dim > 1) THEN
   !$OMP DO SIMD COLLAPSE(2) SCHEDULE(STATIC)
   DO l = nz_min_2, nz_part_2
      DO j = nx_min_2, nx_part_2
         CALL ALPHASPLIT_Y (dummy1(:), alpha2_y(:,j,l), j, l, nm_f)
      END DO
   END DO
   !$OMP END DO
END IF
IF(n_dim > 2) THEN
   !$OMP DO SIMD COLLAPSE(2) SCHEDULE(STATIC)
   DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
		   CALL ALPHASPLIT_Z (dummy1(:), alpha2_z(:,j,k), j, k, nm_f)
	   END DO
   END DO
   !$OMP END DO
END IF
!$OMP END PARALLEL

END SUBROUTINE