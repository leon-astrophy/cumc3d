!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for WENO (DM and NM seperately)
! Written by Leung Shing Chi in 2016
! The subroutine takes ARRAY as input/output and SIGN
! for doing odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, check Boundary1D_FULL.f90
! For hybrid boundaries, such as the quadrant star 
! Specific modifications are needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_DM (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_part_1+1-j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = array(j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = fac_r * array(length_step_r_part_1+1-j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = array(length_step_r_part_1,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_min_part_1-j) = array(:,length_step_z_part_1+1-j)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
      array(:,length_step_z_min_part_1-j) = fac_z * array(:,length_step_z_min_part_1-1+j)
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      array(:,length_step_z_min_part_1-j) = array(:,length_step_z_min_part_1)
   ENDDO
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = array(:,length_step_z_min_part_1-1+j)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = fac_z * array(:,length_step_z_part_1+1-j)
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = array(:,length_step_z_part_1)
   ENDDO
ENDIF

END SUBROUTINE boundary1D_DM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for dark matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_NM (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_part_2+1-j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = array(j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = fac_r * array(length_step_r_part_2+1-j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = array(length_step_r_part_2,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_min_part_2-j) = array(:,length_step_z_part_2+1-j)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
      array(:,length_step_z_min_part_2-j) = fac_z * array(:,length_step_z_min_part_2-1+j)
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      array(:,length_step_z_min_part_2-j) = array(:,length_step_z_min_part_2)
   ENDDO
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = array(:,length_step_z_min_part_2-1+j)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = fac_z * array(:,length_step_z_part_2+1-j)
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = array(:,length_step_z_part_2)
   ENDDO
ENDIF

END SUBROUTINE boundary1D_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for WENO
! Written by Leung Shing Chi in 2016
! The subroutine takes ARRAY as input/output and SIGN
! for doing odd/even parity extension
! Notice that this subroutines worked for full size arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_DMFULL (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_1+1-j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      array(length_step_r_1+j,:) = array(j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      array(length_step_r_1+j,:) = fac_r * array(length_step_r_1+1-j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      array(length_step_r_1+j,:) = array(length_step_r_1,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      array(:,1-j) = array(:,length_step_z_1+1-j)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
      array(:,1-j) = fac_z * array(:,j)
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      array(:,1-j) = array(:,1)
   ENDDO
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = array(:,j)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = fac_z * array(:,length_step_z_1+1-j)
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = array(:,length_step_z_1)
   ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for Normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_NMFULL (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_2+1-j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      array(length_step_r_2+j,:) = array(j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      array(length_step_r_2+j,:) = fac_r * array(length_step_r_2+1-j,:)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      array(length_step_r_2+j,:) = array(length_step_r_2,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      array(:,1-j) = array(:,length_step_z_2+1-j)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
      array(:,1-j) = fac_z * array(:,j)
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      array(:,1-j) = array(:,1)
   ENDDO
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = array(:,j)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = fac_z * array(:,length_step_z_2+1-j)
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = array(:,length_step_z_2)
   ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the primitive variables
! Written by Leung Shing Chi in 2016
! The subroutine takes the full U arrays as input
! and do the odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, switch the checkstep_flag = 0 
! For hybrid boundaries, such as the quadrant star
! Specific modifications are needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE boundary2D_DM
USE DEFINITION
USE helmeos_module
USE Levelset_module
USE Turb_module
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the first (r-inner) boundary
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      u1(1-j,:,:) = u1(length_step_r_part_1+1-j,:,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin1, imax1
         u1(1-j,:,i) = bfac_r(i) * u1(j,:,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      u1(1-j,:,:) = u1(1,:,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      u1(length_step_r_part_1+j,:,:) = u1(j,:,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin1, imax1
         u1(length_step_r_part_1+j,:,i) = bfac_r(i) * u1(length_step_r_part_1+1-j,:,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      u1(length_step_r_part_1+j,:,:) = u1(length_step_r_part_1,:,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      u1(:,length_step_z_min_part_1-j,:) = u1(:,length_step_z_part_1+1-j,:)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
    DO j = 1, 5, 1
      DO i = imin1, imax1
         u1(:,length_step_z_min_part_1-j,i) = bfac_z(i) * u1(:,length_step_z_min_part_1-1+j,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      u1(:,length_step_z_min_part_1-j,:) = u1(:,length_step_z_min_part_1,:)
   ENDDO             
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      u1(:,length_step_z_part_1+j,:) = u1(:,length_step_z_min_part_1-1+j,:)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin1, imax1
         u1(:,length_step_z_part_1+j,i) = bfac_z(i) * u1(:,length_step_z_part_1+1-j,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      u1(:,length_step_z_part_1+j,:) = u1(:,length_step_z_part_1,:)
   ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for Normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE boundary2D_NM
USE DEFINITION
USE helmeos_module
USE Levelset_module
USE Turb_module
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the first (r-inner) boundary
IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      u2(1-j,:,:) = u2(length_step_r_part_2+1-j,:,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin2, imax2
         u2(1-j,:,i) = bfac_r(i) * u2(j,:,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      u2(1-j,:,:) = u2(1,:,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      u2(length_step_r_part_2+j,:,:) = u2(j,:,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin2, imax2
         u2(length_step_r_part_2+j,:,i) = bfac_r(i) * u2(length_step_r_part_2+1-j,:,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      u2(length_step_r_part_2+j,:,:) = u2(length_step_r_part_2,:,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      u2(:,length_step_z_min_part_2-j,:) = u2(:,length_step_z_part_2+1-j,:)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
      DO i = imin2, imax2
         u2(:,length_step_z_min_part_2-j,i) = bfac_z(i) * u2(:,length_step_z_min_part_2-1+j,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      u2(:,length_step_z_min_part_2-j,:) = u2(:,length_step_z_min_part_2,:)
   ENDDO             
ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      u2(:,length_step_z_part_2+j,:) = u2(:,length_step_z_min_part_2-1+j,:)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      DO i = imin2, imax2
         u2(:,length_step_z_part_2+j,i) = bfac_z(i) * u2(:,length_step_z_part_2+1-j,i)
      END DO
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN
   DO j = 1, 5, 1
      u2(:,length_step_z_part_2+j,:) = u2(:,length_step_z_part_2,:)
   ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the chemical isotopes
!
! The subroutine automatically extends the isotopes
! assuming all isotopes behave like scalars
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, switch the checkstep_flag = 0
! For hybrid boundaries, such as the quadrant star
! Specific modifications are needed
!
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY2D_X()
USE definition
USE helmeos_module
IMPLICIT NONE

!Dummy variables
integer :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(boundary_flag(1) == 0) THEN
   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(length_step_r_part_2+1-j,:,:)
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(j,:,:)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(1,:,:)
   ENDDO
ENDIF

! Do the second (r-outer) boundary
IF(boundary_flag(2) == 0) THEN
   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(j,:,:)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(length_step_r_part_2+1-j,:,:)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(length_step_r_part_2,:,:)                  
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary
IF(boundary_flag(3) == 0) THEN
   DO j = 1, 5, 1
      xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_part_2+1-j,:)                     
   ENDDO
ELSEIF(boundary_flag(3) == 1) THEN                 
   DO j = 1, 5, 1
       xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_min_part_2-1+j,:)
   ENDDO
ELSEIF(boundary_flag(3) == 2) THEN
   DO j = 1, 5, 1              
      xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_min_part_2,:)
   ENDDO             

ENDIF

! Do the fourth (z-outer) boundary
IF(boundary_flag(4) == 0) THEN
   DO j = 1, 5, 1
      xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_min_part_2-1+j,:)
   ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
   DO j = 1, 5, 1
      xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_part_2+1-j,:)
   ENDDO
ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_part_2,:)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary2D_X