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
USE OMP_LIB
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
   DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
      if(prim1(j,k,l,irho1) <= rho_min1) then
         prim1(j,k,l,:) = prim1_a(:)
      endif
   END DO

   ! Make sure the update also applies to ghost cells
   CALL BOUNDARYP_DM

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the normal matter

! Set up the threshold
rho_min2 = 1.1D0 * prim2_a(irho2)

! Check the density of normal matter
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
   if(prim2(j,k,l,irho2) <= rho_min2) then   
      prim2(j,k,l,:) = prim2_a(:)
      epsilon2(j,k,l) = eps2_a
      temp2(j,k,l) = temp2_a
   endif
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
            if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j-1,k,l,irho1) >= rho_min1) then
	            if(j > x_grid1) then
                  x_grid1 = j
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
   IF(fullx_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      x_grid1 = nx_1/2

      ! Find the outermost grid which is the star surface
      DO k = 1, ny_1
         DO l = 1, nz_1
            DO j = 1, nx_1
               if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j+1,k,l,irho1) >= rho_min1) then
	               if(j < x_grid1) then
                     x_grid1 = j
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_1 = nx_min_1
      ! Set the effective length_step_z
      IF(x_grid1 == nx_1/2) then
         nx_min_1 = 1
      ELSEIF(x_grid1 /= nx_1/2 .and. x_grid1 > 15) THEN
         nx_min_1 = x_grid1 - 15
      ELSE
         nx_min_1 = 1
      ENDIF

   ELSE
      nx_min_1 = 1
   END IF

   IF(n_dim > 1) THEN
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain ny_part_1, ny_min_1

   ! Initialized r_grid1
   y_grid1 = 0          

   ! Find the outermost grid which is the star surface    
   DO j = 1, nx_1
      DO l = 1, nz_1
         DO k = ny_1, 1, -1
            if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j,k-1,l,irho1) >= rho_min1) then
	            if(k > y_grid1) then
                  y_grid1 = k
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_1 = ny_part_1
   ! Set the effective nx_part_1
   IF(y_grid1 == 0) then
      ny_part_1 = ny_1
   ELSEIF(y_grid1 /= 0 .and. y_grid1 < ny_1 - 15) then
      ny_part_1 = y_grid1 + 15
   ELSE
      ny_part_1 = ny_1
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(fully_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      y_grid1 = ny_1/2

      ! Find the outermost grid which is the star surface
      DO j = 1, nx_1
         DO l = 1, nz_1
            DO k = 1, ny_1
               if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j,k+1,l,irho1) >= rho_min1) then
	               if(k < y_grid1) then
                     y_grid1 = k
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_1 = ny_min_1
      ! Set the effective length_step_z
      IF(y_grid1 == ny_1/2) then
         ny_min_1 = 1
      ELSEIF(y_grid1 /= ny_1/2 .and. y_grid1 > 15) THEN
         ny_min_1 = y_grid1 - 15
      ELSE
         ny_min_1 = 1
      ENDIF

   ELSE
      ny_min_1 = 1
   END IF
   END IF

   IF(n_dim > 2) THEN
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain nz_part_1, nz_min_1

   ! Initialized r_grid1
   z_grid1 = 0          

   ! Find the outermost grid which is the star surface    
   DO j = 1, nx_1
      DO k = 1, ny_1
         DO l = nz_1, 1, -1
            if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j,k,l-1,irho1) >= rho_min1) then
	            if(l > z_grid1) then
                  z_grid1 = l
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_1 = nz_part_1
   ! Set the effective nx_part_1
   IF(z_grid1 == 0) then
      nz_part_1 = nz_1
   ELSEIF(z_grid1 /= 0 .and. z_grid1 < nz_1 - 15) then
      nz_part_1 = z_grid1 + 15
   ELSE
      nz_part_1 = nz_1
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(fullz_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      z_grid1 = nz_1/2

      ! Find the outermost grid which is the star surface
      DO j = 1, nx_1
         DO k = 1, ny_1
            DO l = 1, nz_1
               if(prim1(j,k,l,irho1) == prim1_a(irho1) .and. prim1(j,k,l+1,irho1) >= rho_min1) then
	               if(l < z_grid1) then
                     z_grid1 = l
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_1 = nz_min_1
      ! Set the effective length_step_z
      IF(z_grid1 == nz_1/2) then
         nz_min_1 = 1
      ELSEIF(z_grid1 /= nz_1/2 .and. z_grid1 > 15) THEN
         nz_min_1 = z_grid1 - 15
      ELSE
         nz_min_1 = 1
      ENDIF

   ELSE
      nz_min_1 = 1
   END IF
   END IF

ELSE

   ! Nothing changed if you do not need the check
   nx_min_1 = 1
   ny_min_1 = 1
   nz_min_1 = 1
   nx_part_1 = nx_1
   ny_part_1 = ny_1
   nz_part_1 = nz_1

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation

IF(checkstepnm_flag) THEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain nx_part_2, nx_min_2

   ! Initialized r_grid1
   x_grid2 = 0          

   ! Find the outermost grid which is the star surface    
   DO k = 1, ny_2
      DO l = 1, nz_2
         DO j = nx_2, 1, -1
            if(prim2(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j-1,k,l,irho2) >= rho_min2) then
	            if(j > x_grid2) then
                  x_grid2 = j
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_2 = nx_part_2
   ! Set the effective nx_part_2
   IF(x_grid2 == 0) then
      nx_part_2 = nx_2
   ELSEIF(x_grid2 /= 0 .and. x_grid2 < nx_2 - 15) then
      nx_part_2 = x_grid2 + 15
   ELSE
      nx_part_2 = nx_2
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(fullx_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      x_grid2 = nx_2/2

      ! Find the outermost grid which is the star surface
      DO k = 1, ny_2
         DO l = 1, nz_2
            DO j = 1, nx_2
               if(prim2(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j+1,k,l,irho2) >= rho_min2) then
	               if(j < x_grid2) then
                     x_grid2 = j
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_2 = nx_min_2
      ! Set the effective length_step_z
      IF(x_grid2 == nx_2/2) then
         nx_min_2 = 1
      ELSEIF(x_grid2 /= nx_2/2 .and. x_grid2 > 15) THEN
         nx_min_2 = x_grid2 - 15
      ELSE
         nx_min_2 = 1
      ENDIF

   ELSE
      nx_min_2 = 1
   END IF

   IF(n_dim > 1) THEN
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain ny_part_2, ny_min_2

   ! Initialized r_grid1
   y_grid2 = 0          

   ! Find the outermost grid which is the star surface    
   DO j = 1, nx_2
      DO l = 1, nz_2
         DO k = ny_2, 1, -1
            if(prim2(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j,k-1,l,irho2) >= rho_min2) then
	            if(k > y_grid2) then
                  y_grid2 = k
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_2 = ny_part_2
   ! Set the effective nx_part_2
   IF(y_grid2 == 0) then
      ny_part_2 = ny_2
   ELSEIF(y_grid2 /= 0 .and. y_grid2 < ny_2 - 15) then
      ny_part_2 = y_grid2 + 15
   ELSE
      ny_part_2 = ny_2
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(fully_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      y_grid2 = ny_2/2

      ! Find the outermost grid which is the star surface
      DO j = 1, nx_2
         DO l = 1, nz_2
            DO k = 1, ny_2
               if(prim2(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j,k+1,l,irho2) >= rho_min2) then
	               if(k < y_grid2) then
                     y_grid2 = k
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_2 = ny_min_2
      ! Set the effective length_step_z
      IF(y_grid2 == ny_2/2) then
         ny_min_2 = 1
      ELSEIF(y_grid2 /= ny_2/2 .and. y_grid2 > 15) THEN
         ny_min_2 = y_grid2 - 15
      ELSE
         ny_min_2 = 1
      ENDIF

   ELSE
      ny_min_2 = 1
   END IF
   END IF

   IF(n_dim > 2) THEN
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain nz_part_2, nz_min_2

   ! Initialized r_grid1
   z_grid2 = 0          

   ! Find the outermost grid which is the star surface    
   DO j = 1, nx_2
      DO k = 1, ny_2
         DO l = nz_2, 1, -1
            if(prim1(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j,k,l-1,irho2) >= rho_min2) then
	            if(l > z_grid2) then
                  z_grid2 = l
                  CYCLE
	            endif
            endif
         END DO                 
      ENDDO
   ENDDO

   temp_step_2 = nz_part_2
   ! Set the effective nx_part_2
   IF(z_grid2 == 0) then
      nz_part_2 = nz_2
   ELSEIF(z_grid2 /= 0 .and. z_grid2 < nz_2 - 15) then
      nz_part_2 = z_grid2 + 15
   ELSE
      nz_part_2 = nz_2
   ENDIF

   ! Find the innermost grid which is the star surface  
   IF(fullz_flag) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid2
      z_grid2 = nz_2/2

      ! Find the outermost grid which is the star surface
      DO j = 1, nx_2
         DO k = 1, ny_2
            DO l = 1, nz_2
               if(prim2(j,k,l,irho2) == prim2_a(irho2) .and. prim2(j,k,l+1,irho2) >= rho_min2) then
	               if(l < z_grid2) then
                     z_grid2 = l
                     CYCLE
	               endif
               endif
            END DO                 
         ENDDO
      ENDDO

      temp_step_2 = nz_min_2
      ! Set the effective length_step_z
      IF(z_grid2 == nz_2/2) then
         nz_min_2 = 1
      ELSEIF(z_grid2 /= nz_2/2 .and. z_grid2 > 15) THEN
         nz_min_2 = z_grid2 - 15
      ELSE
         nz_min_2 = 1
      ENDIF
      END IF
   ELSE
      nz_min_2 = 1
   END IF

ELSE

   ! Nothing changed if you do not need the check
   nx_min_2 = 1
   ny_min_2 = 1
   nz_min_2 = 1
   nx_part_2 = nx_2
   ny_part_2 = ny_2
   nz_part_2 = nz_2

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CHECKRHO