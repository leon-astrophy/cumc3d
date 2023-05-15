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
! the structure is i(e-field component)(at coordinate inferfaces) !
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

! indices for cell-centered electric fields !
INTEGER :: ix = 1
INTEGER :: iy = 2
INTEGER :: iz = 3

! Sign operator for cell-centered electric fields !
INTEGER :: ecell_xin(1:3)
INTEGER :: ecell_xout(1:3)
INTEGER :: ecell_yin(1:3)
INTEGER :: ecell_yout(1:3)
INTEGER :: ecell_zin(1:3)
INTEGER :: ecell_zout(1:3)

! cell and face-centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: eface

! cell centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecell_x, ecell_y, ecell_z

! electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: efield_x, efield_y, efield_z

! line integrals !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: lint_x, lint_y, lint_z

! mass flux !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: mflux_x, mflux_y, mflux_z

! electric field gradients at interface !
! the structure is i(e-field component)(differentiate axis) !
! e.g. iyx means the dey/dx !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: dedsface_down, dedsface_up
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: dedscorn_down, dedscorn_up

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
ALLOCATE(ecell_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(ecell_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(ecell_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

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

! Allocate !
ALLOCATE(mflux_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(mflux_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(mflux_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(dedscorn_up(1:6,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(dedsface_up(1:6,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(dedscorn_down(1:6,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(dedsface_down(1:6,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

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
          dbdy = (SIN(yF2(k))*lint_y(j,k,l) - SIN(yF2(k-1))*lint_y(j,k-1,l))/(x2(j)*dcos2(k))!*SIN(y2(k))*dy2(k))
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
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2
        eface (iyx,j,k,l) = flux_2 (ibz,j,k,l)
        eface (izx,j,k,l) = - flux_2 (iby,j,k,l)
        mflux_x (j,k,l) = flux_2 (irho2,j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(dir_in == y_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2, nx_part_2
        eface (ixy,j,k,l) = - flux_2 (ibz,j,k,l)
        eface (izy,j,k,l) = flux_2 (ibx,j,k,l)
        mflux_y (j,k,l) = flux_2 (irho2,j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
ELSEIF(dir_in == z_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        eface (ixz,j,k,l) = flux_2 (iby,j,k,l)
        eface (iyz,j,k,l) = - flux_2 (ibx,j,k,l)
        mflux_z (j,k,l) = flux_2 (irho2,j,k,l)
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

! Real !
REAL*8 :: signs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find cell-centered electric fields !
!$OMP PARALLEL 
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - 1, nz_part_2 
  DO k = ny_min_2 - 1, ny_part_2 
    DO j = nx_min_2 - 1, nx_part_2 
      ecell_x (j,k,l) = prim2(ivel2_z,j,k,l)*prim2(iby,j,k,l) - prim2(ivel2_y,j,k,l)*prim2(ibz,j,k,l)
      ecell_y (j,k,l) = prim2(ivel2_x,j,k,l)*prim2(ibz,j,k,l) - prim2(ivel2_z,j,k,l)*prim2(ibx,j,k,l)
      ecell_z (j,k,l) = prim2(ivel2_y,j,k,l)*prim2(ibx,j,k,l) - prim2(ivel2_x,j,k,l)*prim2(iby,j,k,l)
    END DO
  END DO
END DO
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constrained transport !

! Assign cell-interface emf for 1D/2D problem !
IF(n_dim == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2 - 1, nz_part_2 
    DO k = ny_min_2 - 1, ny_part_2 
      DO j = nx_min_2 - 1, nx_part_2 
        eface (ixy,j,k,l) = ecell_x (j,k,l)
        eface (ixz,j,k,l) = ecell_x (j,k,l)
        eface (izy,j,k,l) = ecell_z (j,k,l)
        eface (iyz,j,k,l) = ecell_y (j,k,l)
      END DO
    END DO
  END DO
  !$OMP END DO
ELSEIF(n_dim == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2 
      DO j = nx_min_2 - 1, nx_part_2
        eface (ixz,j,k,l) = ecell_x (j,k,l)
        eface (iyz,j,k,l) = ecell_y (j,k,l)
      END DO
    END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

! Boundary conditions !
CALL BOUNDARY_EFIELDS

! Get electric field gradients !
!$OMP PARALLEL PRIVATE(signs)
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
      dedsface_up(ixy,j,k,l) = (ecell_x(j,k+1,l) - eface(ixy,j,k,l))/(0.5D0*dy2(k+1))
      dedsface_down(ixy,j,k,l) = (eface(ixy,j,k,l) - ecell_x(j,k,l))/(0.5D0*dy2(k))
      dedsface_up(ixz,j,k,l) = (ecell_x(j,k,l+1) - eface(ixz,j,k,l))/(0.5D0*dz2(l+1))
      dedsface_down(ixz,j,k,l) = (eface(ixz,j,k,l) - ecell_x(j,k,l))/(0.5D0*dz2(l))
      dedsface_up(iyx,j,k,l) = (ecell_y(j+1,k,l) - eface(iyx,j,k,l))/(0.5D0*dx2(j+1))
      dedsface_down(iyx,j,k,l) = (eface(iyx,j,k,l) - ecell_y(j,k,l))/(0.5D0*dx2(j))
      dedsface_up(iyz,j,k,l) = (ecell_y(j,k,l+1) - eface(iyz,j,k,l))/(0.5D0*dz2(l+1))
      dedsface_down(iyz,j,k,l) = (eface(iyz,j,k,l) - ecell_y(j,k,l))/(0.5D0*dz2(l))
      dedsface_up(izx,j,k,l) = (ecell_z(j+1,k,l) - eface(izx,j,k,l))/(0.5D0*dx2(j+1))
      dedsface_down(izx,j,k,l) = (eface(izx,j,k,l) - ecell_z(j,k,l))/(0.5D0*dx2(j))
      dedsface_up(izy,j,k,l) = (ecell_z(j,k+1,l) - eface(izy,j,k,l))/(0.5D0*dy2(k+1))
      dedsface_down(izy,j,k,l) = (eface(izy,j,k,l) - ecell_z(j,k,l))/(0.5D0*dy2(k))
    END DO
  END DO
END DO
!$OMP END DO

! Grid corner electric fiedl gradients using mass_flux as upwinding !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
      signs = SIGN(1.0d0, mflux_x(j,k,l))
      dedscorn_up(iyz,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(iyz,j,k,l) + (1.0D0 - signs)*dedsface_up(iyz,j+1,k,l))
      dedscorn_up(izy,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(izy,j,k,l) + (1.0D0 - signs)*dedsface_up(izy,j+1,k,l))
      dedscorn_down(iyz,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(iyz,j,k,l) + (1.0D0 - signs)*dedsface_down(iyz,j+1,k,l))
      dedscorn_down(izy,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(izy,j,k,l) + (1.0D0 - signs)*dedsface_down(izy,j+1,k,l))
      signs = SIGN(1.0d0, mflux_y(j,k,l))
      dedscorn_up(ixz,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(ixz,j,k,l) + (1.0D0 - signs)*dedsface_up(ixz,j,k+1,l))
      dedscorn_up(izx,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(izx,j,k,l) + (1.0D0 - signs)*dedsface_up(izx,j,k+1,l))
      dedscorn_down(ixz,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(ixz,j,k,l) + (1.0D0 - signs)*dedsface_down(ixz,j,k+1,l))
      dedscorn_down(izx,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(izx,j,k,l) + (1.0D0 - signs)*dedsface_down(izx,j,k+1,l))
      signs = SIGN(1.0d0, mflux_z(j,k,l))
      dedscorn_up(ixy,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(ixy,j,k,l) + (1.0D0 - signs)*dedsface_up(ixy,j,k,l+1))
      dedscorn_up(iyx,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_up(iyx,j,k,l) + (1.0D0 - signs)*dedsface_up(iyx,j,k,l+1))
      dedscorn_down(ixy,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(ixy,j,k,l) + (1.0D0 - signs)*dedsface_down(ixy,j,k,l+1))
      dedscorn_down(iyx,j,k,l) = 0.5D0*((1.0D0 + signs)*dedsface_down(iyx,j,k,l) + (1.0D0 - signs)*dedsface_down(iyx,j,k,l+1))
    END DO
  END DO
END DO
!$OMP END DO

! Add emf !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
      efield_x(j,k,l) = 0.25D0*(eface(ixy,j,k,l) + eface(ixy,j,k,l+1) + eface(ixz,j,k,l) + eface(ixz,j,k+1,l)) &
                       + 0.125D0*((dy2(k)*dedscorn_down(ixy,j,k,l) - dy2(k+1)*dedscorn_up(ixy,j,k,l)) & 
                       + (dz2(l)*dedscorn_down(ixz,j,k,l) - dz2(l+1)*dedscorn_up(ixz,j,k,l)))
      efield_y(j,k,l) = 0.25D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1) + eface(iyz,j,k,l) + eface(iyz,j+1,k,l)) &
                       + 0.125D0*((dx2(j)*dedscorn_down(iyx,j,k,l) - dx2(j+1)*dedscorn_up(iyx,j,k,l)) &
                       + (dz2(l)*dedscorn_down(iyz,j,k,l) - dz2(l+1)*dedscorn_up(iyz,j,k,l)))
      efield_z(j,k,l) = 0.25D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l) + eface(izy,j,k,l) + eface(izy,j+1,k,l)) &
                       + 0.125D0*((dx2(j)*dedscorn_down(izx,j,k,l) - dx2(j+1)*dedscorn_up(izx,j,k,l)) &
                       + (dy2(k)*dedscorn_down(izy,j,k,l) - dy2(k+1)*dedscorn_up(izy,j,k,l)))
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions to the fluxes (electric fields)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_EFIELDS
use definition
implicit none

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the y-electric fields at the z-interface !
! Apply boundary conditions, for the z-electric fields at the y-interface !
! Apply boundary conditions, for the cell-centered x,y,z electric fields !
! This is for the x-direction boundaries

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN  
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = eface(izy,nx_part_2+1-j,k,l)     
    eface(iyz,nx_min_2-j,k,l) = eface(iyz,nx_part_2+1-j,k,l)
    ecell_x (nx_min_2-j,k,l) = ecell_x (nx_part_2+1-j,k,l)
    ecell_y (nx_min_2-j,k,l) = ecell_y (nx_part_2+1-j,k,l)
    ecell_z (nx_min_2-j,k,l) = ecell_z (nx_part_2+1-j,k,l)                  
  ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = eface(izy,nx_min_2,k,l)
    eface(iyz,nx_min_2-j,k,l) = eface(iyz,nx_min_2,k,l)
    ecell_x (nx_min_2-j,k,l) = ecell_x (nx_min_2,k,l)
    ecell_y (nx_min_2-j,k,l) = ecell_y (nx_min_2,k,l)
    ecell_z (nx_min_2-j,k,l) = ecell_z (nx_min_2,k,l)
  ENDDO
ELSEIF(boundary_flag(1) >= 2) THEN                 
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_min_2-j,k,l) = efac_xin(izy) * eface(izy,nx_min_2-1+j,k,l)
    eface(iyz,nx_min_2-j,k,l) = efac_xin(iyz) * eface(iyz,nx_min_2-1+j,k,l)
    ecell_x (nx_min_2-j,k,l) = ecell_xin(ix) * ecell_x (nx_min_2-1+j,k,l)
    ecell_y (nx_min_2-j,k,l) = ecell_xin(iy) * ecell_y (nx_min_2-1+j,k,l)
    ecell_z (nx_min_2-j,k,l) = ecell_xin(iz) * ecell_z (nx_min_2-1+j,k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = eface(izy,nx_min_2-1+j,k,l)
    eface(iyz,nx_part_2+j,k,l) = eface(iyz,nx_min_2-1+j,k,l)
    ecell_x (nx_part_2+j,k,l) = ecell_x (nx_min_2-1+j,k,l)
    ecell_y (nx_part_2+j,k,l) = ecell_y (nx_min_2-1+j,k,l)
    ecell_z (nx_part_2+j,k,l) = ecell_z (nx_min_2-1+j,k,l)
  ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = eface(izy,nx_part_2,k,l)
    eface(iyz,nx_part_2+j,k,l) = eface(iyz,nx_part_2,k,l)
    ecell_x (nx_part_2+j,k,l) = ecell_x (nx_part_2,k,l)
    ecell_y (nx_part_2+j,k,l) = ecell_y (nx_part_2,k,l)
    ecell_z (nx_part_2+j,k,l) = ecell_z (nx_part_2,k,l)
  ENDDO
ELSEIF(boundary_flag(2) >= 2) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
    eface(izy,nx_part_2+j,k,l) = efac_xout(izy) * eface(izy,nx_part_2+1-j,k,l)
    eface(iyz,nx_part_2+j,k,l) = efac_xout(iyz) * eface(iyz,nx_part_2+1-j,k,l)
    ecell_x (nx_part_2+j,k,l) = ecell_xout(ix) * ecell_x (nx_part_2+1-j,k,l)
    ecell_y (nx_part_2+j,k,l) = ecell_xout(iy) * ecell_y (nx_part_2+1-j,k,l)
    ecell_z (nx_part_2+j,k,l) = ecell_xout(iz) * ecell_z (nx_part_2+1-j,k,l)
  ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the x-electric fields at the z-interface !
! Apply boundary conditions, for the z-electric fields at the x-interface !
! Apply boundary conditions, for the cell-centered x,y,z electric fields !
! This is for the y-direction boundaries

! Do the inner boundary
IF(boundary_flag(3) == 0) THEN    
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = eface(izx,j,ny_part_2+1-k,l)  
    eface(ixz,j,ny_min_2-k,l) = eface(ixz,j,ny_part_2+1-k,l)   
    ecell_x (j,ny_min_2-k,l) = ecell_x (j,ny_part_2+1-k,l)
    ecell_y (j,ny_min_2-k,l) = ecell_y (j,ny_part_2+1-k,l)
    ecell_z (j,ny_min_2-k,l) = ecell_z (j,ny_part_2+1-k,l)                
  ENDDO
ELSEIF(boundary_flag(3) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = eface(izx,j,ny_min_2,l)
    eface(ixz,j,ny_min_2-k,l) = eface(ixz,j,ny_min_2,l)
    ecell_x (j,ny_min_2-k,l) = ecell_x (j,ny_min_2,l)
    ecell_y (j,ny_min_2-k,l) = ecell_y (j,ny_min_2,l)
    ecell_z (j,ny_min_2-k,l) = ecell_z (j,ny_min_2,l)
  ENDDO
ELSEIF(boundary_flag(3) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_min_2-k,l) = efac_yin(izx) * eface(izx,j,ny_min_2-1+k,l)
    eface(ixz,j,ny_min_2-k,l) = efac_yin(ixz) * eface(ixz,j,ny_min_2-1+k,l)
    ecell_x (j,ny_min_2-k,l) = ecell_yin(ix) * ecell_x (j,ny_min_2-1+k,l)
    ecell_y (j,ny_min_2-k,l) = ecell_yin(iy) * ecell_y (j,ny_min_2-1+k,l)
    ecell_z (j,ny_min_2-k,l) = ecell_yin(iz) * ecell_z (j,ny_min_2-1+k,l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(4) == 0) THEN    
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = eface(izx,j,ny_min_2-1+k,l)
    eface(ixz,j,ny_part_2+k,l) = eface(ixz,j,ny_min_2-1+k,l)
    ecell_x (j,ny_part_2+k,l) = ecell_x (j,ny_min_2-1+k,l)
    ecell_y (j,ny_part_2+k,l) = ecell_y (j,ny_min_2-1+k,l)
    ecell_z (j,ny_part_2+k,l) = ecell_z (j,ny_min_2-1+k,l)
  ENDDO
ELSEIF(boundary_flag(4) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = eface(izx,j,ny_part_2,l)
    eface(ixz,j,ny_part_2+k,l) = eface(ixz,j,ny_part_2,l)
    ecell_x (j,ny_part_2+k,l) = ecell_x (j,ny_part_2,l)
    ecell_y (j,ny_part_2+k,l) = ecell_y (j,ny_part_2,l)
    ecell_z (j,ny_part_2+k,l) = ecell_z (j,ny_part_2,l)
  ENDDO
ELSEIF(boundary_flag(4) >= 2) THEN          
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = 1:3, l = nz_min_2-1:nz_part_2)
    eface(izx,j,ny_part_2+k,l) = efac_yout(izx) * eface(izx,j,ny_part_2+1-k,l)
    eface(ixz,j,ny_part_2+k,l) = efac_yout(ixz) * eface(ixz,j,ny_part_2+1-k,l)
    ecell_x (j,ny_part_2+k,l) = ecell_yout(ix) * ecell_x (j,ny_part_2+1-k,l)
    ecell_y (j,ny_part_2+k,l) = ecell_yout(iy) * ecell_y (j,ny_part_2+1-k,l)
    ecell_z (j,ny_part_2+k,l) = ecell_yout(iz) * ecell_z (j,ny_part_2+1-k,l)
  ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply boundary conditions, for the x-electric fields at the y-interface !
! Apply boundary conditions, for the y-electric fields at the x-interface !
! Apply boundary conditions, for the cell-centered x,y,z electric fields !
! This is for the z-direction boundaries

! Do the inner boundary
IF(boundary_flag(5) == 0) THEN      
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = eface(iyx,j,k,nz_part_2+1-l)
    eface(ixy,j,k,nz_min_2-l) = eface(ixy,j,k,nz_part_2+1-l)       
    ecell_x (j,k,nz_min_2-l) = ecell_x (j,k,nz_part_2+1-l)
    ecell_y (j,k,nz_min_2-l) = ecell_y (j,k,nz_part_2+1-l)
    ecell_z (j,k,nz_min_2-l) = ecell_z (j,k,nz_part_2+1-l)           
  ENDDO
ELSEIF(boundary_flag(5) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = eface(iyx,j,k,nz_min_2)
    eface(ixy,j,k,nz_min_2-l) = eface(ixy,j,k,nz_min_2)
    ecell_x (j,k,nz_min_2-l) = ecell_x (j,k,nz_min_2)
    ecell_y (j,k,nz_min_2-l) = ecell_y (j,k,nz_min_2)
    ecell_z (j,k,nz_min_2-l) = ecell_z (j,k,nz_min_2)
  ENDDO
ELSEIF(boundary_flag(5) >= 2) THEN                 
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_min_2-l) = efac_zin(iyx) * eface(iyx,j,k,nz_min_2-1+l)
    eface(ixy,j,k,nz_min_2-l) = efac_zin(ixy) * eface(ixy,j,k,nz_min_2-1+l)
    ecell_x (j,k,nz_min_2-l) = ecell_zin(ix) * ecell_x (j,k,nz_min_2-1+l)
    ecell_y (j,k,nz_min_2-l) = ecell_zin(iy) * ecell_y (j,k,nz_min_2-1+l)
    ecell_z (j,k,nz_min_2-l) = ecell_zin(iz) * ecell_z (j,k,nz_min_2-1+l)
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(6) == 0) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = eface(iyx,j,k,nz_min_2-1+l)
    eface(ixy,j,k,nz_part_2+l) = eface(ixy,j,k,nz_min_2-1+l)
    ecell_x (j,k,nz_part_2+l) = ecell_x (j,k,nz_min_2-1+l)
    ecell_y (j,k,nz_part_2+l) = ecell_y (j,k,nz_min_2-1+l)
    ecell_z (j,k,nz_part_2+l) = ecell_z (j,k,nz_min_2-1+l)
  ENDDO
ELSEIF(boundary_flag(6) == 1) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = eface(iyx,j,k,nz_part_2)
    eface(ixy,j,k,nz_part_2+l) = eface(ixy,j,k,nz_part_2)
    ecell_x (j,k,nz_part_2+l) = ecell_x (j,k,nz_part_2)
    ecell_y (j,k,nz_part_2+l) = ecell_y (j,k,nz_part_2)
    ecell_z (j,k,nz_part_2+l) = ecell_z (j,k,nz_part_2)
  ENDDO
ELSEIF(boundary_flag(6) >= 2) THEN
  DO CONCURRENT(j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = 1:3)
    eface(iyx,j,k,nz_part_2+l) = efac_zout(iyx) * eface(iyx,j,k,nz_part_2+1-l)
    eface(ixy,j,k,nz_part_2+l) = efac_zout(ixy) * eface(ixy,j,k,nz_part_2+1-l)
    ecell_x (j,k,nz_part_2+l) = ecell_zout(ix) * ecell_x (j,k,nz_part_2+1-l)
    ecell_y (j,k,nz_part_2+l) = ecell_zout(iy) * ecell_y (j,k,nz_part_2+1-l)
    ecell_z (j,k,nz_part_2+l) = ecell_zout(iz) * ecell_z (j,k,nz_part_2+1-l)
  ENDDO
ENDIF

END SUBROUTINE
    
end module MHD_module