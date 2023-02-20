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

SUBROUTINE BOUNDARY1D_NM (array, domain, signx_in, signx_out, signy_in, signy_out, signz_in, signz_out)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL*8, INTENT (INOUT), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3) :: array

! Input parity
INTEGER, INTENT (IN) :: domain
INTEGER, INTENT (IN) :: signx_in, signx_out
INTEGER, INTENT (IN) :: signy_in, signy_out
INTEGER, INTENT (IN) :: signz_in, signz_out

! Dummy variables
INTEGER :: i, j, k, l

! Integer for domain size !
INTEGER :: nx_min, nx_max
INTEGER :: ny_min, ny_max
INTEGER :: nz_min, nz_max

! Parity factor
INTEGER :: fac_xin, fac_yin, fac_zin
INTEGER :: fac_xout, fac_yout, fac_zout

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
IF(signx_in == 0) THEN
  fac_xin = 1
ELSEIF(signx_in == 1) THEN
  fac_xin = -1
END IF

IF(signx_out == 0) THEN
  fac_xout = 1
ELSEIF(signx_in == 1) THEN
  fac_xout = -1
END IF

! y-direction !
IF(signy_in == 0) THEN
  fac_yin = 1
ELSEIF(signy_in == 1) THEN
  fac_yin = -1
END IF

IF(signy_out == 0) THEN
  fac_yout = 1
ELSEIF(signy_in == 1) THEN
  fac_yout = -1
END IF

! z-direciton !
IF(signz_in == 0) THEN
  fac_zin = 1
ELSEIF(signz_in == 1) THEN
  fac_zin = -1
END IF

IF(signz_out == 0) THEN
  fac_zout = 1
ELSEIF(signz_in == 1) THEN
  fac_zout = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN     
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_min-j,k,l) = array(nx_max+1-j,k,l)
  ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_min-j,k,l) = array(nx_min,k,l)
  ENDDO
ELSEIF(boundary_flag(1) >= 2) THEN                 
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_min-j,k,l) = fac_xin * array(nx_min-1+j,k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_max+j,k,l) = array(nx_min-1+j,k,l)
  ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_max+j,k,l) = array(nx_max,k,l)
  ENDDO
ELSEIF(boundary_flag(2) >= 2) THEN 
  DO CONCURRENT(j = 1:3, k = ny_min:ny_max, l = nz_min:nz_max)
    array(nx_max+j,k,l) = fac_xout * array(nx_max+1-j,k,l)
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
      array(j,ny_min-k,l) = array(j,ny_min,l)
    ENDDO
  ELSEIF(boundary_flag(3) >= 2) THEN                 
    DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
      array(j,ny_min-k,l) = fac_yin * array(j,ny_min-1+k,l)
    ENDDO
  ENDIF

  ! Do the outer boundary
  IF(boundary_flag(3) == 0) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
      array(j,ny_max+k,l) = array(j,ny_min-1+k,l)
    ENDDO
  ELSEIF(boundary_flag(4) == 1) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
      array(j,ny_max+k,l) = array(j,ny_max,l)
    ENDDO
  ELSEIF(boundary_flag(4) >= 2) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = 1:3, l = nz_min:nz_max)
      array(j,ny_max+k,l) = fac_yout * array(j,ny_max+1-k,l)
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
      array(j,k,nz_min-l) = array(j,k,nz_min)
    ENDDO
  ELSEIF(boundary_flag(5) >= 2) THEN         
    DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
      array(j,k,nz_min-l) = fac_zin * array(j,k,nz_min-1+l)
    ENDDO
  ENDIF

  ! Do the outer boundary
  IF(boundary_flag(6) == 0) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
      array(j,k,nz_max+l) = array(j,k,nz_min-1+l)
    ENDDO
  ELSEIF(boundary_flag(6) == 1) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
      array(j,k,nz_max+l) = array(j,k,nz_max)
    ENDDO
  ELSEIF(boundary_flag(6) >= 2) THEN
    DO CONCURRENT(j = nx_min:nx_max, k = ny_min:ny_max, l = 1:3)
      array(j,k,nz_max+l) = fac_zout * array(j,k,nz_max+1-l)
    ENDDO
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary1D_NM

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
    prim2(i,nx_min_2-j,k,l) = prim2(i,nx_min_2,k,l)
  ENDDO     
ELSEIF(boundary_flag(1) >= 2) THEN               
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
    prim2(i,nx_min_2-j,k,l) = bfac_xin(i) * prim2(i,nx_min_2-1+j,k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
    prim2(i,nx_part_2+j,k,l) = prim2(i,nx_min_2-1+j,k,l)
  ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
    prim2(i,nx_part_2+j,k,l) = prim2(i,nx_part_2,k,l)
  ENDDO
ELSEIF(boundary_flag(2) >= 2) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
    prim2(i,nx_part_2+j,k,l) = bfac_xout(i) * prim2(i,nx_part_2+1-j,k,l)
  ENDDO
ENDIF

! Or, add your custom boundary conditions !
CALL CUSTOM_BOUNDARY_X

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
      prim2(i,j,ny_min_2-k,l) = prim2(i,j,ny_min_2,l)
    ENDDO     
  ELSEIF(boundary_flag(3) >= 2) THEN                 
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,j,ny_min_2-k,l) = bfac_yin(i) * prim2(i,j,ny_min_2-1+k,l)
    ENDDO          
  ENDIF

  ! Do the outer boundary
  IF(boundary_flag(4) == 0) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,j,ny_part_2+k,l) = prim2(i,j,ny_min_2-1+k,l)
    ENDDO
  ELSEIF(boundary_flag(4) == 1) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,j,ny_part_2+k,l) = prim2(i,j,ny_part_2,l)
    ENDDO
  ELSEIF(boundary_flag(4) >= 2) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2, i = imin2:imax2)
      prim2(i,j,ny_part_2+k,l) = bfac_yout(i) * prim2(i,j,ny_part_2+1-k,l)
    ENDDO
  ENDIF

  ! Or, add your custom boundary conditions !
  CALL CUSTOM_BOUNDARY_Y

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
      prim2(i,j,k,nz_min_2-l) = prim2(i,j,k,nz_min_2)
    ENDDO
  ELSEIF(boundary_flag(5) >= 2) THEN                 
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
      prim2(i,j,k,nz_min_2-l) = bfac_zin(i) * prim2(i,j,k,nz_min_2-1+l)
    ENDDO            
  ENDIF

  ! Do the outer boundary
  IF(boundary_flag(6) == 0) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
      prim2(i,j,k,nz_part_2+l) = prim2(i,j,k,nz_min_2-1+l)
    ENDDO
  ELSEIF(boundary_flag(6) == 1) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
      prim2(i,j,k,nz_part_2+l) = prim2(i,j,k,nz_part_2)
    ENDDO
  ELSEIF(boundary_flag(6) >= 2) THEN
    DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3, i = imin2:imax2)
      prim2(i,j,k,nz_part_2+l) = bfac_zout(i) * prim2(i,j,k,nz_part_2+1-l)
    ENDDO
  ENDIF

  ! Or, add your custom boundary conditions !
  CALL CUSTOM_BOUNDARY_Z

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE