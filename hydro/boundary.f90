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

SUBROUTINE BOUNDARY1D_DM (array, sign, domain)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-2:nx_1+3,-2:ny_1+3,-2:nz_1+3) :: array

! Input parity
INTEGER, INTENT (IN) :: sign, domain

! Dummy variables
INTEGER :: i, j, k, l

! Integer for domain size !
INTEGER :: nx_min, nx_max
INTEGER :: ny_min, ny_max
INTEGER :: nz_min, nz_max

! Parity factor
INTEGER :: fac_x, fac_y, fac_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Setup domain size !
IF(domain == 0) THEN
   nx_min = nx_min_1
   ny_min = ny_min_1
   nz_min = nz_min_1
   nx_max = nx_part_1
   ny_max = ny_part_1
   nz_max = nz_part_1
ELSEIF(domain == 1) THEN
   nx_min = 1
   ny_min = 1
   nz_min = 1
   nx_max = nx_1
   ny_max = ny_1
   nz_max = nz_1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_x = 1
   fac_y = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_x = -1
   fac_y = 1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_x = 1
   fac_y = -1
   fac_z = 1
ELSEIF(sign == 3) THEN
   fac_x = 1
   fac_y = 1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = array(nx_max+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = fac_x * array(nx_min-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = array(nx_min,k,l)
   ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = array(nx_min-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = fac_x * array(nx_max+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = array(nx_max,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = array(j,ny_max+1-k,l)                     
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = fac_y * array(j,ny_min-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = array(j,ny_min,l)
      ENDDO
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = array(j,ny_min-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = fac_y * array(j,ny_max+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = array(j,ny_max,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = array(j,k,nz_max+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = fac_z * array(j,k,nz_min-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = array(j,k,nz_min)
      ENDDO
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = array(j,k,nz_min-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = fac_z * array(j,k,nz_max+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = array(j,k,nz_max)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary1D_DM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_NM (array, sign, domain)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3) :: array

! Input parity
INTEGER, INTENT (IN) :: sign, domain

! Dummy variables
INTEGER :: i, j, k, l

! Integer for domain size !
INTEGER :: nx_min, nx_max
INTEGER :: ny_min, ny_max
INTEGER :: nz_min, nz_max

! Parity factor
INTEGER :: fac_x, fac_y, fac_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Setup domain size !
IF(domain == 0) THEN
   nx_min = nx_min_2
   ny_min = ny_min_2
   nz_min = nz_min_2
   nx_max = nx_part_2
   ny_max = ny_part_2
   nz_max = nz_part_2
ELSEIF(domain == 1) THEN
   nx_min = 1
   ny_min = 1
   nz_min = 1
   nx_max = nx_2
   ny_max = ny_2
   nz_max = nz_2
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_x = 1
   fac_y = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_x = -1
   fac_y = 1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_x = 1
   fac_y = -1
   fac_z = 1
ELSEIF(sign == 3) THEN
   fac_x = 1
   fac_y = 1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = array(nx_max+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = fac_x * array(nx_min-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_min-j,k,l) = array(nx_min,k,l)
   ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = array(nx_min-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = fac_x * array(nx_max+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
      array(nx_max+j,k,l) = array(nx_max,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = array(j,ny_max+1-k,l)                     
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = fac_y * array(j,ny_min-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_min-k,l) = array(j,ny_min,l)
      ENDDO
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = array(j,ny_min-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = fac_y * array(j,ny_max+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
         array(j,ny_max+k,l) = array(j,ny_max,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = array(j,k,nz_max+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = fac_z * array(j,k,nz_min-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_min-l) = array(j,k,nz_min)
      ENDDO
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = array(j,k,nz_min-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = fac_z * array(j,k,nz_max+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
         array(j,k,nz_max+l) = array(j,k,nz_max)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary1D_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the conservative variables
! Written by Leung Shing Chi in 2016
! The subroutine takes the full U arrays as input
! and do the odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, switch the checkstep_flag = 0 
! For hybrid boundaries, such as the quadrant star
! Specific modifications are needed
! Note: DM hasnt parallelize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYU_DM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      cons1(i,nx_min_1-j,k,l) = cons1(i,nx_part_1+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      cons1(i,nx_min_1-j,k,l) = bfac_x(i) * cons1(i,nx_min_1-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      cons1(i,nx_min_1-j,k,l) = cons1(i,nx_min_1,k,l)
   ENDDO             
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      cons1(i,nx_part_1+j,k,l) = cons1(i,nx_min_1-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,nx_part_1+j,k,l) = bfac_x(i) * cons1(i,nx_part_1+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      cons1(i,nx_part_1+j,k,l) = cons1(i,nx_part_1,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_min_1-k,l) = cons1(i,j,ny_part_1+1-k,l)                    
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_min_1-k,l) = bfac_y(i) * cons1(i,j,ny_min_1-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_min_1-k,l) = cons1(i,j,ny_min_1,l)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_part_1+k,l) = cons1(i,j,ny_min_1-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
     DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_part_1+k,l) = bfac_y(i) * cons1(i,j,ny_part_1+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         cons1(i,j,ny_part_1+k,l) = cons1(i,j,ny_part_1,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_min_1-l) = cons1(i,j,k,nz_part_1+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_min_1-l) = bfac_z(i) * cons1(i,j,k,nz_min_1-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_min_1-l) = cons1(i,j,k,nz_min_1)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_part_1+l) = cons1(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_part_1+l) = bfac_z(i) * cons1(i,j,k,nz_part_1+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         cons1(i,j,k,nz_part_1+l) = cons1(i,j,k,nz_part_1)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the primitive variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYP_DM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      prim1(i,nx_min_1-j,k,l) = prim1(i,nx_part_1+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      prim1(i,nx_min_1-j,k,l) = bfac_x(i) * prim1(i,nx_min_1-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      prim1(i,nx_min_1-j,k,l) = prim1(i,nx_min_1,k,l)
   ENDDO             
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      prim1(i,nx_part_1+j,k,l) = prim1(i,nx_min_1-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,nx_part_1+j,k,l) = bfac_x(i) * prim1(i,nx_part_1+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
      prim1(i,nx_part_1+j,k,l) = prim1(i,nx_part_1,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_min_1-k,l) = prim1(i,j,ny_part_1+1-k,l)                    
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_min_1-k,l) = bfac_y(i) * prim1(i,j,ny_min_1-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_min_1-k,l) = prim1(i,j,ny_min_1,l)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_part_1+k,l) = prim1(i,j,ny_min_1-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
     DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_part_1+k,l) = bfac_y(i) * prim1(i,j,ny_part_1+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = 1:3, l = nz_min_1:nz_part_1, i = imin1:imax1)
         prim1(i,j,ny_part_1+k,l) = prim1(i,j,ny_part_1,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_min_1-l) = prim1(i,j,k,nz_part_1+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_min_1-l) = bfac_z(i) * prim1(i,j,k,nz_min_1-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_min_1-l) = prim1(i,j,k,nz_min_1)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_part_1+l) = prim1(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_part_1+l) = bfac_z(i) * prim1(i,j,k,nz_part_1+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = 1:3, i = imin1:imax1)
         prim1(i,j,k,nz_part_1+l) = prim1(i,j,k,nz_part_1)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for Normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYU_NM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      cons2(i,nx_min_2-j,k,l) = cons2(i,nx_part_2+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      cons2(i,nx_min_2-j,k,l) = bfac_x(i) * cons2(i,nx_min_2-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      cons2(i,nx_min_2-j,k,l) = cons2(i,nx_min_2,k,l)
   ENDDO             
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      cons2(i,nx_part_2+j,k,l) = cons2(i,nx_min_2-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,nx_part_2+j,k,l) = bfac_x(i) * cons2(i,nx_part_2+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      cons2(i,nx_part_2+j,k,l) = cons2(i,nx_part_2,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_min_2-k,l) = cons2(i,j,ny_part_2+1-k,l)                    
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_min_2-k,l) = bfac_y(i) * cons2(i,j,ny_min_2-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_min_2-k,l) = cons2(i,j,ny_min_2,l)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_part_2+k,l) = cons2(i,j,ny_min_2-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
     DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_part_2+k,l) = bfac_y(i) * cons2(i,j,ny_part_2+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         cons2(i,j,ny_part_2+k,l) = cons2(i,j,ny_part_2,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_min_2-l) = cons2(i,j,k,nz_part_2+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_min_2-l) = bfac_z(i) * cons2(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_min_2-l) = cons2(i,j,k,nz_min_2)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_part_2+l) = cons2(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_part_2+l) = bfac_z(i) * cons2(i,j,k,nz_part_2+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         cons2(i,j,k,nz_part_2+l) = cons2(i,j,k,nz_part_2)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for normal matter primitive variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYP_NM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,nx_min_2-j,k,l) = prim2(i,nx_part_2+1-j,k,l)                     
   ENDDO
ELSEIF(boundary_flag(1) == 1) THEN                 
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,nx_min_2-j,k,l) = bfac_x(i) * prim2(i,nx_min_2-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,nx_min_2-j,k,l) = prim2(i,nx_min_2,k,l)
   ENDDO             
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,nx_part_2+j,k,l) = prim2(i,nx_min_2-1+j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,nx_part_2+j,k,l) = bfac_x(i) * prim2(i,nx_part_2+1-j,k,l)
   ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
   DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,nx_part_2+j,k,l) = prim2(i,nx_part_2,k,l)
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

! Do the inner boundary
IF(n_dim > 1) THEN
   IF(boundary_flag(3) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_min_2-k,l) = prim2(i,j,ny_part_2+1-k,l)                    
      ENDDO
   ELSEIF(boundary_flag(3) == 1) THEN                 
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_min_2-k,l) = bfac_y(i) * prim2(i,j,ny_min_2-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(3) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_min_2-k,l) = prim2(i,j,ny_min_2,l)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(4) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_part_2+k,l) = prim2(i,j,ny_min_2-1+k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 1) THEN
     DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_part_2+k,l) = bfac_y(i) * prim2(i,j,ny_part_2+1-k,l)
      ENDDO
   ELSEIF(boundary_flag(4) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
         prim2(i,j,ny_part_2+k,l) = prim2(i,j,ny_part_2,l)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

! Do the inner boundary
IF(n_dim > 2) THEN
   IF(boundary_flag(5) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_min_2-l) = prim2(i,j,k,nz_part_2+1-l)                     
      ENDDO
   ELSEIF(boundary_flag(5) == 1) THEN                 
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_min_2-l) = bfac_z(i) * prim2(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(5) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_min_2-l) = prim2(i,j,k,nz_min_2)
      ENDDO             
   ENDIF

   ! Do the outer boundary
   IF(boundary_flag(6) == 0) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_part_2+l) = prim2(i,j,k,nz_min_2-1+l)
      ENDDO
   ELSEIF(boundary_flag(6) == 1) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_part_2+l) = bfac_z(i) * prim2(i,j,k,nz_part_2+1-l)
      ENDDO
   ELSEIF(boundary_flag(6) == 2) THEN
      DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
         prim2(i,j,k,nz_part_2+l) = prim2(i,j,k,nz_part_2)
      ENDDO
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE