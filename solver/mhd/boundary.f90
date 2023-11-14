!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign boundary conditions in a lump-sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY 
USE DEFINITION
IMPLICIT NONE

! Call boundary condition !
call BOUNDARYP_NM
call BOUNDARY1D_NM (epsilon, even, even, even, even, even, even)

END SUBROUTINE

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

SUBROUTINE BOUNDARY1D_NM (array, signx_in, signx_out, signy_in, signy_out, signz_in, signz_out)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL*8, INTENT (INOUT), DIMENSION (-2:nx+3,-2:ny+3,-2:nz+3) :: array

! Input parity
INTEGER, INTENT (IN) :: signx_in, signx_out
INTEGER, INTENT (IN) :: signy_in, signy_out
INTEGER, INTENT (IN) :: signz_in, signz_out

! Dummy variables
INTEGER :: i, j, k, l

! Parity factor
INTEGER :: fac_xin, fac_yin, fac_zin
INTEGER :: fac_xout, fac_yout, fac_zout

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

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

!$OMP PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(1) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(1-j,k,l) = array(nx+1-j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(1) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(1-j,k,l) = array(1,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(1) >= 2) THEN     

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)   
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)       
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(1-j,k,l) = fac_xin * array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(2) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(nx+j,k,l) = array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(2) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(nx+j,k,l) = array(nx,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(2) >= 2) THEN 

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k =  -2, ny + 3 
      DO j = 1, 3 
        array(nx+j,k,l) = fac_xout * array(nx+1-j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(3) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) 
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,1-k,l) = array(j,ny+1-k,l)                     
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(3) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)   
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,1-k,l) = array(j,1,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(3) >= 2) THEN   

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)                  
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,1-k,l) = fac_yin * array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(3) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,ny+k,l) = array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(4) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,ny+k,l) = array(j,ny,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(4) >= 2) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        array(j,ny+k,l) = fac_yout * array(j,ny+1-k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(5) == 0) THEN   

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)   
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)      
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,1-l) = array(j,k,nz+1-l)                     
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(5) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,1-l) = array(j,k,1)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(5) >= 2) THEN   

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)         
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,1-l) = fac_zin * array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(6) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)    
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,nz+l) = array(j,k,l)
      END DO
    END DO
  ENDDO
  !$OMP END DO 

ELSEIF(boundary_flag(6) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,nz+l) = array(j,k,nz)
      END DO
    END DO
  ENDDO
  !$OMP END DO

ELSEIF(boundary_flag(6) >= 2) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)  
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, 3
    DO k =  -2, ny + 3 
      DO j = -2, nx + 3
        array(j,k,nz+l) = fac_zout * array(j,k,nz+1-l)
      END DO
    END DO
  ENDDO
  !$OMP END DO

END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'boundary_nm = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE boundary1D_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for normal matter primitive variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYP_NM
USE DEFINITION
USE MHD_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

!$OMP PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(1) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,1-j,k,l) = prim(imin:ibx-1,nx+1-j,k,l)    
        bcell(ibx:ibz,1-j,k,l) = bcell(ibx:ibz,nx+1-j,k,l)
        prim(iby,1-j,k,l) = prim(iby,nx+1-j,k,l)
        prim(ibz,1-j,k,l) = prim(ibz,nx+1-j,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 

ELSEIF(boundary_flag(1) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,1-j,k,l) = prim(imin:ibx-1,1,k,l)
        bcell(ibx:ibz,1-j,k,l) = bcell(ibx:ibz,1,k,l)
        prim(iby,1-j,k,l) = prim(iby,1,k,l)
        prim(ibz,1-j,k,l) = prim(ibz,1,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 

ELSEIF(boundary_flag(1) >= 2) THEN    

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,1-j,k,l) = bfac_xin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,1-j,k,l) = bfac_xin(ibx:ibz) * bcell(ibx:ibz,j,k,l)
        prim(iby,1-j,k,l) = bfac_xin(iby) * prim(iby,j,k,l)
        prim(ibz,1-j,k,l) = bfac_xin(ibz) * prim(ibz,j,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(2) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,nx+j,k,l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,nx+j,k,l) = bcell(ibx:ibz,j,k,l)
        prim(iby,nx+j,k,l) = prim(iby,j,k,l)
        prim(ibz,nx+j,k,l) = prim(ibz,j,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 
    
ELSEIF(boundary_flag(2) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,nx+j,k,l) = prim(imin:ibx-1,nx,k,l)
        bcell(ibx:ibz,nx+j,k,l) = bcell(ibx:ibz,nx,k,l)
        prim(iby,nx+j,k,l) = prim(iby,nx,k,l)
        prim(ibz,nx+j,k,l) = prim(ibz,nx,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 

ELSEIF(boundary_flag(2) >= 2) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = -2, ny + 3
      DO j = 1, 3
        prim(imin:ibx-1,nx+j,k,l) = bfac_xout(imin:ibx-1) * prim(imin:ibx-1,nx+1-j,k,l)
        bcell(ibx:ibz,nx+j,k,l) = bfac_xout(ibx:ibz) * bcell(ibx:ibz,nx+1-j,k,l)
        prim(iby,nx+j,k,l) = bfac_xout(iby) * prim(iby,nx+1-j,k,l)
        prim(ibz,nx+j,k,l) = bfac_xout(ibz) * prim(ibz,nx+1-j,k,l)
      END DO
    END DO               
  ENDDO
  !$OMP END DO 

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(3) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,1-k,l) = prim(imin:ibx-1,j,ny+1-k,l) 
        bcell(ibx:ibz,j,1-k,l) = bcell(ibx:ibz,j,ny+1-k,l)            
        prim(ibx,j,1-k,l) = prim(ibx,j,ny+1-k,l)
        prim(ibz,j,1-k,l) = prim(ibz,j,ny+1-k,l)               
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(3) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,1-k,l) = prim(imin:ibx-1,j,1,l)
        bcell(ibx:ibz,j,1-k,l) = bcell(ibx:ibz,j,1,l)
        prim(ibx,j,1-k,l) = prim(ibx,j,1,l)
        prim(ibz,j,1-k,l) = prim(ibz,j,1,l)          
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(3) >= 2) THEN    

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,1-k,l) = bfac_yin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,1-k,l) = bfac_yin(ibx:ibz) * bcell(ibx:ibz,j,k,l)
        prim(ibx,j,1-k,l) = bfac_yin(ibx) * prim(ibx,j,k,l)
        prim(ibz,j,1-k,l) = bfac_yin(ibz) * prim(ibz,j,k,l)         
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(4) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,ny+k,l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,ny+k,l) = bcell(ibx:ibz,j,k,l)
        prim(ibx,j,ny+k,l) = prim(ibx,j,k,l)
        prim(ibz,j,ny+k,l) = prim(ibz,j,k,l)        
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(4) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,ny+k,l) = prim(imin:ibx-1,j,ny,l)
        bcell(ibx:ibz,j,ny+k,l) = bcell(ibx:ibz,j,ny,l)
        prim(ibx,j,ny+k,l) = prim(ibx,j,ny,l)
        prim(ibz,j,ny+k,l) = prim(ibz,j,ny,l)       
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(4) >= 2) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = -2, nz + 3
    DO k = 1, 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,ny+k,l) = bfac_yout(imin:ibx-1) * prim(imin:ibx-1,j,ny+1-k,l)
        bcell(ibx:ibz,j,ny+k,l) = bfac_yout(ibx:ibz) * bcell(ibx:ibz,j,ny+1-k,l)
        prim(ibx,j,ny+k,l) = bfac_yout(ibx) * prim(ibx,j,ny+1-k,l)
        prim(ibz,j,ny+k,l) = bfac_yout(ibz) * prim(ibz,j,ny+1-k,l)       
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(boundary_flag(5) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,1-l) = prim(imin:ibx-1,j,k,nz+1-l)           
        bcell(ibx:ibz,j,k,1-l) = bcell(ibx:ibz,j,k,nz+1-l)               
        prim(ibx,j,k,1-l) = prim(ibx,j,k,nz+1-l)
        prim(iby,j,k,1-l) = prim(iby,j,k,nz+1-l)                   
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(5) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,1-l) = prim(imin:ibx-1,j,k,1)
        bcell(ibx:ibz,j,k,1-l) = bcell(ibx:ibz,j,k,1)
        prim(ibx,j,k,1-l) = prim(ibx,j,k,1)
        prim(iby,j,k,1-l) = prim(iby,j,k,1)                  
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(5) >= 2) THEN  

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,1-l) = bfac_zin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,k,1-l) = bfac_zin(ibx:ibz) * bcell(ibx:ibz,j,k,l)
        prim(ibx,j,k,1-l) = bfac_zin(ibx) * prim(ibx,j,k,l)
        prim(iby,j,k,1-l) = bfac_zin(iby) * prim(iby,j,k,l)              
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the outer boundary
IF(boundary_flag(6) == 0) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,nz+l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,k,nz+l) = bcell(ibx:ibz,j,k,l)
        prim(ibx,j,k,nz+l) = prim(ibx,j,k,l)
        prim(iby,j,k,nz+l) = prim(iby,j,k,l)            
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(6) == 1) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,nz+l) = prim(imin:ibx-1,j,k,nz)   
        bcell(ibx:ibz,j,k,nz+l) = bcell(ibx:ibz,j,k,nz)
        prim(ibx,j,k,nz+l) = prim(ibx,j,k,nz)
        prim(iby,j,k,nz+l) = prim(iby,j,k,nz)         
      END DO
    END DO               
  ENDDO 
  !$OMP END DO

ELSEIF(boundary_flag(6) >= 2) THEN

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, 3
    DO k = -2, ny + 3
      DO j = -2, nx + 3
        prim(imin:ibx-1,j,k,nz+l) = bfac_zout(imin:ibx-1) * prim(imin:ibx-1,j,k,nz+1-l)   
        bcell(ibx:ibz,j,k,nz+l) = bfac_zout(ibx:ibz) * bcell(ibx:ibz,j,k,nz+1-l) 
        prim(ibx,j,k,nz+l) = bfac_zout(ibx) * prim(ibx,j,k,nz+1-l)
        prim(iby,j,k,nz+l) = bfac_zout(iby) * prim(iby,j,k,nz+1-l)       
      END DO
    END DO               
  ENDDO 
  !$OMP END DO  

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Or, add your custom boundary conditions !
CALL CUSTOM_BOUNDARY_X

! Or, add your custom boundary conditions !
CALL CUSTOM_BOUNDARY_Y

! Or, add your custom boundary conditions !
CALL CUSTOM_BOUNDARY_Z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'boundary_p = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE