!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine ensures that the density does not go below rho_atm
! and replace the grids with atmospheric density (tempearture and 
! chemical composition and so on) if found
! Written by Leung Shing Chi in 2016
! The subroutine do all the check automatically
! Notice that this subroutines also check the size of 
! the hydro array, and reduced the simulation grid-number
! to boost the calculation,
! i.e. (1:length_step_r_part, 1:length_step_z_part).
! For full array extension, switch the checkstep_flag = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKRHO
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Threshold for atmosphere density
real (DP) :: rho_min1, rho_min2

! Temporal variables !
INTEGER :: temp_step_1, temp_step_2
INTEGER :: x_grid1, y_grid1, z_grid1
INTEGER :: x_grid2, y_grid2, z_grid2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First check the DM
if(DM_flag) then

   ! Set up the threshold
   rho_min1 = 1.1D0 * prim1_a(irho1)

	! IF low-density grid is found
	! Replace them with stmospheric
   DO j = nx_min_1, nx_part_1
      DO k = ny_min_1, ny_part_1
         DO l = nz_min_1, nz_part_1
            if(prim1(j,k,l,irho1) <= rho_min1) then
               prim1(j,k,l,:) = prim1_a(:)
            endif
         END DO
      ENDDO
   ENDDO

   ! Make sure the update also applies to ghost cells
   CALL BOUNDARYP_DM 

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do the normal matter
rho_min2 = 1.1D0 * prim2_a(irho2)

! Check the density of normal matter
DO j = nx_min_2, nx_part_2
   DO k = ny_min_2, ny_part_2
      DO l = nz_min_2, nz_part_2
         if(prim2(j,k,l,irho2) <= rho_min2) then   
            prim2(j,k,l,:) = prim2_a(:)
            epsilon2(j,k,l) = epsilon2_a
            temp2(j,k,l) = temp2_a
         endif
      END DO
   enddo
enddo

! Make sure the ghost cell knows the udpate
CALL BOUNDARYP_NM  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation

IF(checkstepdm_flag) THEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain nx_part_1, nx_min_1

   ! Initialized r_grid1
   x_grid1 = 0          

   ! Find the outermost grid which is the star surface    
   DO k = 1, ny_1
      DO l = 1, nz_1
         DO j = nx_1, 1, -1
            if(prim1(j,k,l,irho1) == rho1_a .and. prim1(j-1,k,l,irho1) >= rho_min1) then
	            if(j - 1 > x_grid1) then
                  x_grid1 = j - 1
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_1 = nx_part_1
   ! Set the effective nx_part_1
   IF(x_grid1 == 0) then
      nx_part_1 = nx_1
   ELSEIF(x_grid1 /= 0 .and. x_grid1 < nx_1 - 15) then
      nx_part_1 = x_grid1 + 15
   ELSE
      nx_part_1 = nx_1
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(full_x) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      x_grid1 = nx_1/2

      ! Find the outermost grid which is the star surface
      DO k = 1, ny_1
         DO l = 1, nz_1
            DO j = 1, nx_1
               if(prim1(j,k,l,irho1) >= rho_min1 .and. prim1(j-1,k,l,irho1) == rho1_a) then
	               if(j - 1 > x_grid1) then
                     x_grid1 = j - 1
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_1 = length_step_z_min_part_1
      ! Set the effective length_step_z
      IF(z_grid1 == length_step_z_1/2-1) then
         length_step_z_min_part_1 = 1
      ELSEIF(z_grid1 /= length_step_z_1/2-1 .and. z_grid1 > 15) THEN
         length_step_z_min_part_1 = z_grid1 - 15
      ELSE
         length_step_z_min_part_1 = 1
      ENDIF
      IF(global_time > 0.0D0) THEN
         length_step_z_min_part_1 = min(temp_step_1, length_step_z_min_part_1)
      END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_z_part_1

   IF(hemisphere_flag == 1) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      z_grid1 = length_step_z_1/2-1

      ! Find the outermost grid which is the star surface
      DO j = 1, length_step_r_1, 1
         DO k = length_step_z_1/2-1, 2, -1
            IF(rho1(j,k) >= rho_min1 .and. rho1(j,k-1) == rho1_a) THEN
	            IF(k < z_grid1) THEN
                  z_grid1 = k
   	         ENDIF
	            !exit
            ENDIF
         ENDDO
      ENDDO 

      temp_step_1 = length_step_z_min_part_1
      ! Set the effective length_step_z
      IF(z_grid1 == length_step_z_1/2-1) then
         length_step_z_min_part_1 = 1
      ELSEIF(z_grid1 /= length_step_z_1/2-1 .and. z_grid1 > 15) THEN
         length_step_z_min_part_1 = z_grid1 - 15
      ELSE
         length_step_z_min_part_1 = 1
      ENDIF
      IF(global_time > 0.0D0) THEN
         length_step_z_min_part_1 = min(temp_step_1, length_step_z_min_part_1)
      END IF

   ELSEIF(hemisphere_flag == 0) THEN

      length_step_z_min_part_1 = 1

   ELSE

      STOP 'Check the value of hemisphere_flag'

   ENDIF

   ! Now find the upper limit of length_step_z
   IF(coordinate_flag == 2) THEN 
   length_step_z_part_1 = length_step_z_1
   ELSE
   ! Initialized z_grid1
   z_grid1 = 0

   ! Find the outermost grid which is the star surface
   DO j = 1, length_step_r_1, 1
      DO k = 1, length_step_z_1 - 1, 1
         IF(rho1(j,k) >= rho_min1 .and. rho1(j,k+1) == rho1_a) THEN
            IF(k > z_grid1) then
               z_grid1 = k
            ENDIF
            !exit
         ENDIF
      ENDDO
   ENDDO

   temp_step_1 = length_step_z_part_1
   ! Set the effective length_step_z
   IF(z_grid1 == 0) then
      length_step_z_part_1 = length_step_z_1
   ELSEIF(z_grid1 /= 0 .and. z_grid1 < length_step_z_1 - 15) THEN
      length_step_z_part_1 = z_grid1 + 15
   ELSE
      length_step_z_part_1 = length_step_z_1
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_z_part_1 = max(temp_step_1, length_step_z_part_1)
   END IF
   END IF

ELSE

   ! Nothing changed if you do not need the check
   length_step_z_min_part_1 = 1
   length_step_z_part_1 = length_step_z_1
   length_step_r_part_1 = length_step_r_1

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation

IF(checkstepnm_flag == 1) THEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_r_part

   ! Initialized r_grid2
   r_grid2 = 0          

   ! Find the outermost grid which is the star surface    
   DO k = 1, length_step_z_2, 1     
      DO j = 1, length_step_r_2 - 1, 1
         if(rho2(j,k) >= rho_min2 .and. rho2(j+1,k) == rho2_a) then
	         if(j > r_grid2) then
               r_grid2 = j
	         endif
	         !exit
         endif                 
      ENDDO
   ENDDO

   temp_step_2 = length_step_r_part_2
   ! Set the effective length_step_r
   IF(r_grid2 == 0) then
      length_step_r_part_2 = length_step_r_2
   ELSEIF(r_grid2 /= 0 .and. r_grid2 < length_step_r_2 - 15) then
      length_step_r_part_2 = r_grid2 + 15
   ELSE
      length_step_r_part_2 = length_step_r_2
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_r_part_2 = max(temp_step_2, length_step_r_part_2)
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_z_part

   IF(hemisphere_flag == 1) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid2
      z_grid2 = length_step_z_2/2-1

      ! Find the outermost grid which is the star surface
      DO j = 1, length_step_r_2, 1
         DO k = length_step_z_2/2-1, 2, -1
            IF(rho2(j,k) >= rho_min2 .and. rho2(j,k-1) == rho2_a) THEN
	            IF(k < z_grid2) THEN
                  z_grid2 = k
   	         ENDIF
	            !exit
            ENDIF
         ENDDO
      ENDDO 

      temp_step_2 = length_step_z_min_part_2
      ! Set the effective length_step_z
      IF(z_grid2 == length_step_z_2/2-1) then
         length_step_z_min_part_2 = 1
      ELSEIF(z_grid2 /= length_step_z_2/2-1 .and. z_grid2 > 15) THEN
         length_step_z_min_part_2 = z_grid2 - 15
      ELSE
         length_step_z_min_part_2 = 1
      ENDIF
      IF(global_time > 0.0D0) THEN
         length_step_z_min_part_2 = min(temp_step_2, length_step_z_min_part_2)
      END IF

   ELSEIF(hemisphere_flag == 0) THEN

      length_step_z_min_part_2 = 1

   ELSE

      STOP 'Check the value of hemisphere_flag'

   ENDIF

   ! Now find the upper limit of length_step_z
   IF(coordinate_flag == 2) THEN 
   length_step_z_part_2 = length_step_z_2
   ELSE
   ! Initialized z_grid2
   z_grid2 = 0

   ! Find the outermost grid which is the star surface
   DO j = 1, length_step_r_2, 1
      DO k = 1, length_step_z_2 - 1, 1
         IF(rho2(j,k) >= rho_min2 .and. rho2(j,k+1) == rho2_a) THEN
            IF(k > z_grid2) then
               z_grid2 = k
            ENDIF
            !exit
         ENDIF
      ENDDO
   ENDDO

   temp_step_2 = length_step_z_part_2
   ! Set the effective length_step_z
   IF(z_grid2 == 0) then
      length_step_z_part_2 = length_step_z_2
   ELSEIF(z_grid2 /= 0 .and. z_grid2 < length_step_z_2 - 15) THEN
      length_step_z_part_2 = z_grid2 + 15
   ELSE
      length_step_z_part_2 = length_step_z_2
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_z_part_2 = max(temp_step_2, length_step_z_part_2)
   END IF
   END IF

ELSE

   ! Nothing changed if you do not need the check
   length_step_z_min_part_2 = 1
   length_step_z_part_2 = length_step_z_2
   length_step_r_part_2 = length_step_r_2

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CHECKRHO