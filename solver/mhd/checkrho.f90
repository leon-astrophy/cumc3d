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
REAL*8 :: rho_min1, rho_min2, factor, diff

! Temporal variables !
INTEGER :: x_grid1, y_grid1, z_grid1
INTEGER :: x_grid2, y_grid2, z_grid2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the normal matter

! Check the density of normal matter, assign floor for hydrodynamic variables !
IF(custom_floor) THEN

  ! custom variable floor !
  CALL CUSTOMFLOOR

ELSE

  !$OMP PARALLEL DO PRIVATE(diff, factor, rho_min2) COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        rho_min2 = 1.1D0 * prim2_a(irho2,j,k,l)
        diff = prim2(irho2,j,k,l) - rho_min2
        factor = MAX(SIGN(1.0D0, diff), 0.0D0)
        prim2(irho2:ivel2_z,j,k,l) = factor*prim2(irho2:ivel2_z,j,k,l) + (1.0D0 - factor)*prim2_a(irho2:ivel2_z,j,k,l)
        epsilon2(j,k,l) = factor*epsilon2(j,k,l) + (1.0D0 - factor)*eps2_a(j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

END IF

! Make sure the ghost cell knows the udpate
CALL BOUNDARY1D_NM (epsilon2, part, even, even, even, even, even, even)
CALL BOUNDARYP_NM  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation
IF(checkstepnm_flag) THEN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Obtain nx_part_2, nx_min_2

  ! Initialized r_grid1
  x_grid2 = nx_2

  ! Find the outermost grid which is the star surface    
  DO k = 1, ny_2
    DO l = 1, nz_2
      DO j = 1, nx_2
        if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
          x_grid2 = MIN(j, x_grid2)  
          CYCLE
        endif
      END DO               
    ENDDO
  ENDDO

  ! Set the effective length_step_z
  nx_min_2 = MAX(x_grid2 - 15, 1)

  ! Find the innermost grid which is the star surface  
  ! Now find the lower limit of length_step_z
  ! Initialized z_grid1
  x_grid2 = 1

  ! Find the outermost grid which is the star surface
  DO k = 1, ny_2
    DO l = 1, nz_2
      DO j = nx_2, 1, -1
        if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
          x_grid2 = MAX(j, x_grid2) 
          CYCLE
        endif
      END DO                 
    ENDDO
  ENDDO

  ! Set the effective nx_part_2
  nx_part_2 = MIN(x_grid2 + 15, nx_2)
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Obtain ny_part_2, ny_min_2
  IF(n_dim > 1) THEN

    ! Find the innermost grid which is the star surface  
    ! Now find the lower limit of length_step_z
    ! Initialized z_grid1
    y_grid2 = ny_2

    ! Find the outermost grid which is the star surface
    DO j = 1, nx_2
      DO l = 1, nz_2
        DO k = 1, ny_2
          if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
            y_grid2 = MIN(k, y_grid2)  
            CYCLE
          endif
        END DO                 
      ENDDO
    ENDDO

    ! Set the effective length_step_z
    ny_min_2 = MAX(y_grid2 - 15, 1)

    ! Initialized r_grid1
    y_grid2 = 1 

    ! Find the outermost grid which is the star surface    
    DO j = 1, nx_2
      DO l = 1, nz_2
        DO k = ny_2, 1, -1
          if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
            y_grid2 = MAX(k, y_grid2)  
            CYCLE
          endif
        END DO                 
      ENDDO
    ENDDO

    ! Set the effective nx_part_2
    ny_part_2 = MIN(y_grid2 + 15, ny_2)

  END IF
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Obtain nz_part_2, nz_min_2
  IF(n_dim > 2) THEN

    ! Now find the lower limit of length_step_z
    ! Initialized z_grid2
    z_grid2 = nz_2

    ! Find the outermost grid which is the star surface
    DO j = 1, nx_2
      DO k = 1, ny_2
        DO l = 1, nz_2
          if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
            z_grid2 = MIN(l, z_grid2)  
            CYCLE
          endif
        END DO                 
      ENDDO
    ENDDO

    ! Set the effective length_step_z
    nz_min_2 = MAX(z_grid2 - 15, 1)

    ! Initialized r_grid1
    z_grid2 = 1          

    ! Find the outermost grid which is the star surface    
    DO j = 1, nx_2
      DO k = 1, ny_2
        DO l = nz_2, 1, -1
          if(prim2(irho2,j,k,l) > prim2_a(irho2,j,k,l)) then
            z_grid2 = MAX(l, z_grid2)  
            CYCLE
          end if
        END DO                 
      ENDDO
    ENDDO

    ! Set the effective nx_part_2
    nz_part_2 = MIN(z_grid2 + 15, nz_2)

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