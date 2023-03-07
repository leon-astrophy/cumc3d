!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! MHD module, written by H.S. Leon Chan
! Uses Flux-CT (Balsara 1999 b, Toth 2000) to ensure divergence-less 
! 2nd order accurate scheme
! Updated to the Upwind-CT Scheme (Gardiner and Stone 2015)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MHD_module
use definition   
implicit none

! Maximum divergence B
REAL*8 :: maxDivB

! indices for face-centered electric fields !
! the structure is i(e-field component)(at coordinate inferfaces !
! e.g. iyx means the yth component efield at the xth interface !
INTEGER :: iyx = 1
INTEGER :: izx = 2
INTEGER :: ixy = 3
INTEGER :: izy = 4
INTEGER :: ixz = 5
INTEGER :: iyz = 6

! Sign operator for electric fields !
INTEGER :: efac_xin(1:6)
INTEGER :: efac_xout(1:6)
INTEGER :: efac_yin(1:6)
INTEGER :: efac_yout(1:6)
INTEGER :: efac_zin(1:6)
INTEGER :: efac_zout(1:6)

! face-centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: eface

! electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: efield_x, efield_y, efield_z

! line integrals !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: lint_x, lint_y, lint_z

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine allocates the array necessary for
! the calculation of MHD
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildMHD
use definition
implicit none

! Allocate !
ALLOCATE(eface(1:6,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(efield_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(efield_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(efield_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(lint_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(lint_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(lint_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

end subroutine buildMHD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open file for output 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_MHD
USE DEFINITION
IMPLICIT NONE

! Open !
OPEN (UNIT = 999, FILE = './outfile/star_weno_divb.dat')

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find divergence of magnetic fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FIND_DIVB
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l

! Dummy !
REAL*8 :: dbdx, dbdy, dbdz, divb

! Initialize !
maxDivB = 0.0d0

!$OMP PARALLEL PRIVATE(dbdx, dbdy, dbdz, divb)
IF(coordinate_flag == 0) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2
        dbdx = (lint_x(j,k,l) - lint_x(j-1,k,l))/(dx2(j))
        dbdy = (lint_y(j,k,l) - lint_y(j,k-1,l))/(dy2(k))
        dbdz = (lint_z(j,k,l) - lint_z(j,k,l-1))/(dz2(l))
        divb = dbdx + dbdy + dbdz
        maxDivB = max(maxDivB, divb)
      END DO
    END DO
  END DO
  !$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2
        dbdx = (xF2(j)*lint_x(j,k,l) - xF2(j-1)*lint_x(j-1,k,l))/(x2(j)*dx2(j))
        dbdy = (lint_y(j,k,l) - lint_y(j,k-1,l))/(x2(j)*dy2(k))
        dbdz = (lint_z(j,k,l) - lint_z(j,k,l-1))/(dz2(l))
        divb = dbdx + dbdy + dbdz
        maxDivB = max(maxDivB, divb)
      END DO
    END DO
  END DO
  !$OMP END DO
  ELSEIF(coordinate_flag == 2) THEN
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
    DO l = 1, nz_2
      DO k = 1, ny_2
        DO j = 1, nx_2
          dbdx = (xF2(j)**2*lint_x(j,k,l) - xF2(j-1)**2*lint_x(j-1,k,l))/(x2(j)**2*dx2(j))
          dbdy = (SIN(yF2(k))*lint_y(j,k,l) - SIN(yF2(k-1))*lint_y(j,k-1,l))/(x2(j)*SIN(y2(k))*dy2(k))
          dbdz = (lint_z(j,k,l) - lint_z(j,k,l-1))/(x2(j)*SIN(y2(k))*dz2(l))
          divb = dbdx + dbdy + dbdz
          maxDivB = max(maxDivB, divb)
        END DO
      END DO
    END DO
    !$OMP END DO
END IF
!$OMP END PARALLEL

! output !
WRITE (999, *) global_time, ABS(maxDivB)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backup mhd fluxes after each directional sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MHD_FLUX(dir_in)
use definition
implicit none

! Integer !
INTEGER :: j, k, l

! Integer !
INTEGER, INTENT(IN) :: dir_in

! Do it by direction !
IF(dir_in == x_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  DO l = nz_min_2 - 1, nz_part_2 + 1
    DO k = ny_min_2 - 1, ny_part_2 + 1
      DO j = nx_min_2 - 1, nx_part_2
        eface (iyx,j,k,l) = flux_2 (ibz,j,k,l)
        eface (izx,j,k,l) = - flux_2 (iby,j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(dir_in == y_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  DO l = nz_min_2 - 1, nz_part_2 + 1
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2 + 1
        eface (ixy,j,k,l) = - flux_2 (ibz,j,k,l)
        eface (izy,j,k,l) = flux_2 (ibx,j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(dir_in == z_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2 + 1
      DO j = nx_min_2 - 1, nx_part_2 + 1
        eface (ixz,j,k,l) = flux_2 (iby,j,k,l)
        eface (iyz,j,k,l) = - flux_2 (ibx,j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constrained transport on the mangetic field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE flux_ct
use definition
implicit none

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constrained transport !

! Add emf !
IF(n_dim == 1) THEN

  ! Boundary conditions !
  CALL BOUNDARY_EFIELDS_X

  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2
        efield_y(j,k,l) = 0.5D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1))
        efield_z(j,k,l) = 0.5D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l))
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(n_dim == 2) THEN

  ! Boundary conditions !
  CALL BOUNDARY_EFIELDS_X
  CALL BOUNDARY_EFIELDS_Y

  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2
        efield_x(j,k,l) = 0.5D0*(eface(ixy,j,k,l) + eface(ixy,j,k,l+1))
        efield_y(j,k,l) = 0.5D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1))
        efield_z(j,k,l) = 0.25D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l) + eface(izy,j,k,l) + eface(izy,j+1,k,l))
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(n_dim == 3) THEN

  ! Boundary conditions !
  CALL BOUNDARY_EFIELDS_X
  CALL BOUNDARY_EFIELDS_Y
  CALL BOUNDARY_EFIELDS_Z

  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2
        efield_x(j,k,l) = 0.25D0*(eface(ixy,j,k,l) + eface(ixy,j,k,l+1) + eface(ixz,j,k,l) + eface(ixz,j,k+1,l))
        efield_y(j,k,l) = 0.25D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1) + eface(iyz,j,k,l) + eface(iyz,j+1,k,l))
        efield_z(j,k,l) = 0.25D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l) + eface(izy,j,k,l) + eface(izy,j+1,k,l))
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions to the fluxes (electric fields)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_EFIELDS_X
use definition
implicit none

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the y-electric fields at the x-interface !
! This is for the z-direction boundaries

! Do the inner boundary
IF(boundary_flag(5) == 0) THEN      
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = eface(iyx,j,k,nz_part_2+1-l)                     
  ENDDO
ELSEIF(boundary_flag(5) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = eface(iyx,j,k,nz_min_2)
  ENDDO
ELSEIF(boundary_flag(5) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = efac_zin(iyx) * eface(iyx,j,k,nz_min_2-1+l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(6) == 0) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = eface(iyx,j,k,nz_min_2-1+l)
  ENDDO
ELSEIF(boundary_flag(6) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = eface(iyx,j,k,nz_part_2)
  ENDDO
ELSEIF(boundary_flag(6) >= 2) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = efac_zout(iyx) * eface(iyx,j,k,nz_part_2+1-l)
  ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the z-electric fields at the x-interface !
! This is for the y-direction boundaries

! Do the inner boundary
IF(boundary_flag(3) == 0) THEN    
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = eface(izx,j,ny_part_2+1-k,l)                     
  ENDDO
ELSEIF(boundary_flag(3) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = eface(izx,j,ny_min_2,l)
  ENDDO
ELSEIF(boundary_flag(3) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = efac_yin(izx) * eface(izx,j,ny_min_2-1+k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(4) == 0) THEN    
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = eface(izx,j,ny_min_2-1+k,l)
  ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = eface(izx,j,ny_part_2,l)
  ENDDO
ELSEIF(boundary_flag(4) >= 2) THEN          
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = efac_yout(izx) * eface(izx,j,ny_part_2+1-k,l)
  ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions to the fluxes (electric fields)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_EFIELDS_Y
use definition
implicit none

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the x-electric fields at the y-interface !
! This is for the z-direction boundaries

! Do the inner boundary
IF(boundary_flag(5) == 0) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_min_2-l) = eface(ixy,j,k,nz_part_2+1-l)                     
  ENDDO
ELSEIF(boundary_flag(5) == 1) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_min_2-l) = eface(ixy,j,k,nz_min_2)
  ENDDO
ELSEIF(boundary_flag(5) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_min_2-l) = efac_zin(ixy) * eface(ixy,j,k,nz_min_2-1+l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(6) == 0) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_part_2+l) = eface(ixy,j,k,nz_min_2-1+l)
  ENDDO
ELSEIF(boundary_flag(6) == 1) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_part_2+l) = eface(ixy,j,k,nz_part_2)
  ENDDO
ELSEIF(boundary_flag(6) >= 2) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(ixy,j,k,nz_part_2+l) = efac_zout(ixy) * eface(ixy,j,k,nz_part_2+1-l)
  ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the z-electric fields at the y-interface !
! This is for the x-direction boundaries

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN  
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = eface(izy,nx_part_2+1-j,k,l)                     
  ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = eface(izy,nx_min_2,k,l)
  ENDDO
ELSEIF(boundary_flag(1) >= 2) THEN                 
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = efac_xin(izy) * eface(izy,nx_min_2-1+j,k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = eface(izy,nx_min_2-1+j,k,l)
  ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = eface(izy,nx_part_2,k,l)
  ENDDO
ELSEIF(boundary_flag(2) >= 2) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = efac_xout(izy) * eface(izy,nx_part_2+1-j,k,l)
  ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions to the fluxes (electric fields)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_EFIELDS_Z
use definition
implicit none

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the x-electric fields at the z-interface !
! This is for the y-direction boundaries

! Do the inner boundary
IF(boundary_flag(3) == 0) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_min_2-k,l) = eface(ixz,j,ny_part_2+1-k,l)                     
  ENDDO
ELSEIF(boundary_flag(3) == 1) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_min_2-k,l) = eface(ixz,j,ny_min_2,l)
  ENDDO
ELSEIF(boundary_flag(3) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_min_2-k,l) = efac_yin(ixz) * eface(ixz,j,ny_min_2-1+k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(4) == 0) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_part_2+k,l) = eface(ixz,j,ny_min_2-1+k,l)
  ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_part_2+k,l) = eface(ixz,j,ny_part_2,l)
  ENDDO
ELSEIF(boundary_flag(4) >= 2) THEN     
  DO CONCURRENT(j = nx_min_2:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(ixz,j,ny_part_2+k,l) = efac_yout(ixz) * eface(ixz,j,ny_part_2+1-k,l)
  ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the y-electric fields at the z-interface !
! This is for the x-direction boundaries

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_min_2-j,k,l) = eface(iyz,nx_part_2+1-j,k,l)                     
  ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_min_2-j,k,l) = eface(iyz,nx_min_2,k,l)
  ENDDO
ELSEIF(boundary_flag(1) >= 2) THEN                 
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_min_2-j,k,l) = efac_xin(iyz) * eface(iyz,nx_min_2-1+j,k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_part_2+j,k,l) = eface(iyz,nx_min_2-1+j,k,l)
  ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_part_2+j,k,l) = eface(iyz,nx_part_2,k,l)
  ENDDO
ELSEIF(boundary_flag(2) >= 2) THEN       
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(iyz,nx_part_2+j,k,l) = efac_xout(iyz) * eface(iyz,nx_part_2+1-j,k,l)
  ENDDO
ENDIF

END SUBROUTINE
   
end module MHD_module

