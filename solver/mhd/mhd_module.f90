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

! cell and face-centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: eface

! cell centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecell_x, ecell_y, ecell_z

! electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: efield_x, efield_y, efield_z

! mass flux !
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: mflux_x

! electric field gradients at interface !
! the structure is i(e-field component)(differentiate axis) !
! e.g. iyx means the dey/dx !
REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_f_d_p, deds_f_u_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_f_d_m, deds_f_u_m
REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_c_d, deds_c_u

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
ALLOCATE(eface(iyx:iyz,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(ecell_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(ecell_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(ecell_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(efield_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(efield_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(efield_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
ALLOCATE(mflux_x(-2:ny_2+3,-2:nz_2+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate !
!ALLOCATE(sign_x(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
!ALLOCATE(sign_y(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
!ALLOCATE(sign_z(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate !
!ALLOCATE(deds_c_u(iyx:iyz))
!ALLOCATE(deds_c_d(iyx:iyz))
!ALLOCATE(deds_f_d_p(iyx:iyz))
!ALLOCATE(deds_f_u_p(iyx:iyz))
!ALLOCATE(deds_f_d_m(iyx:iyz))
!ALLOCATE(deds_f_u_m(iyx:iyz))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do it by direction !
IF(dir_in == x_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = nz_min_2 - 1, nz_part_2 + 1
    DO k = ny_min_2 - 1, ny_part_2 + 1
      DO j = nx_min_2 - 1, nx_part_2
        eface (iyx,j,k,l) = flux_2 (ibz,j,k,l)
        eface (izx,j,k,l) = - flux_2 (iby,j,k,l)
        mflux_x (k,l) = flux_2 (irho2,0,k,l)
        !sign_x (j,k,l) = SIGN(1.0d0, flux_2 (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(dir_in == y_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
  DO l = nz_min_2 - 1, nz_part_2 + 1
    DO k = ny_min_2 - 1, ny_part_2
      DO j = nx_min_2 - 1, nx_part_2 + 1
        eface (ixy,j,k,l) = - flux_2 (ibz,j,k,l)
        eface (izy,j,k,l) = flux_2 (ibx,j,k,l)
        !sign_y (j,k,l) = SIGN(1.0d0, flux_2 (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(dir_in == z_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2 + 1
      DO j = nx_min_2 - 1, nx_part_2 + 1
        eface (ixz,j,k,l) = flux_2 (iby,j,k,l)
        eface (iyz,j,k,l) = - flux_2 (ibx,j,k,l)
        !sign_z (j,k,l) = SIGN(1.0d0, flux_2 (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'mhd_flux = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constrained transport on the mangetic field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE flux_ct
use definition
implicit none

! Integer !
INTEGER :: j, k, l

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find cell-centered electric fields !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = nz_min_2 - 1, nz_part_2 + 1
  DO k = ny_min_2 - 1, ny_part_2 + 1
    DO j = nx_min_2 - 1, nx_part_2 + 1
      ecell_x (j,k,l) = prim2(ivel2_z,j,k,l)*prim2(iby,j,k,l) - prim2(ivel2_y,j,k,l)*prim2(ibz,j,k,l)
      ecell_y (j,k,l) = prim2(ivel2_x,j,k,l)*prim2(ibz,j,k,l) - prim2(ivel2_z,j,k,l)*prim2(ibx,j,k,l)
      ecell_z (j,k,l) = prim2(ivel2_y,j,k,l)*prim2(ibx,j,k,l) - prim2(ivel2_x,j,k,l)*prim2(iby,j,k,l)
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign cell-interface emf for 1D/2D problem !
IF(n_dim == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
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
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(n_dim == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
  DO l = nz_min_2 - 1, nz_part_2
    DO k = ny_min_2 - 1, ny_part_2 
      DO j = nx_min_2 - 1, nx_part_2
        eface (ixz,j,k,l) = ecell_x (j,k,l)
        eface (iyz,j,k,l) = ecell_y (j,k,l)
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! constrained transport !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! First, get electric field gradient at cell interface !
      !deds_f_u_m(ixy) = (ecell_x(j,k+1,l) - eface(ixy,j,k,l))/(0.5D0*dy2(k+1))
      !deds_f_d_m(ixy) = (eface(ixy,j,k,l) - ecell_x(j,k,l))/(0.5D0*dy2(k))
      !deds_f_u_m(ixz) = (ecell_x(j,k,l+1) - eface(ixz,j,k,l))/(0.5D0*dz2(l+1))
      !deds_f_d_m(ixz) = (eface(ixz,j,k,l) - ecell_x(j,k,l))/(0.5D0*dz2(l))
      !deds_f_u_m(iyx) = (ecell_y(j+1,k,l) - eface(iyx,j,k,l))/(0.5D0*dx2(j+1))
      !deds_f_d_m(iyx) = (eface(iyx,j,k,l) - ecell_y(j,k,l))/(0.5D0*dx2(j))
      !deds_f_u_m(iyz) = (ecell_y(j,k,l+1) - eface(iyz,j,k,l))/(0.5D0*dz2(l+1))
      !deds_f_d_m(iyz) = (eface(iyz,j,k,l) - ecell_y(j,k,l))/(0.5D0*dz2(l))
      !deds_f_u_m(izx) = (ecell_z(j+1,k,l) - eface(izx,j,k,l))/(0.5D0*dx2(j+1))
      !deds_f_d_m(izx) = (eface(izx,j,k,l) - ecell_z(j,k,l))/(0.5D0*dx2(j))
      !deds_f_u_m(izy) = (ecell_z(j,k+1,l) - eface(izy,j,k,l))/(0.5D0*dy2(k+1))
      !deds_f_d_m(izy) = (eface(izy,j,k,l) - ecell_z(j,k,l))/(0.5D0*dy2(k))

      ! Also get the graident at an upper grid !
      !deds_f_u_p(ixy) = (ecell_x(j,k+1,l+1) - eface(ixy,j,k,l+1))/(0.5D0*dy2(k+1))
      !deds_f_d_p(ixy) = (eface(ixy,j,k,l+1) - ecell_x(j,k,l+1))/(0.5D0*dy2(k))
      !deds_f_u_p(ixz) = (ecell_x(j,k+1,l+1) - eface(ixz,j,k+1,l))/(0.5D0*dz2(l+1))
      !deds_f_d_p(ixz) = (eface(ixz,j,k+1,l) - ecell_x(j,k+1,l))/(0.5D0*dz2(l))
      !deds_f_u_p(iyx) = (ecell_y(j+1,k,l+1) - eface(iyx,j,k,l+1))/(0.5D0*dx2(j+1))
      !deds_f_d_p(iyx) = (eface(iyx,j,k,l+1) - ecell_y(j,k,l+1))/(0.5D0*dx2(j))
      !deds_f_u_p(iyz) = (ecell_y(j+1,k,l+1) - eface(iyz,j+1,k,l))/(0.5D0*dz2(l+1))
      !deds_f_d_p(iyz) = (eface(iyz,j+1,k,l) - ecell_y(j+1,k,l))/(0.5D0*dz2(l))
      !deds_f_u_p(izx) = (ecell_z(j+1,k+1,l) - eface(izx,j,k+1,l))/(0.5D0*dx2(j+1))
      !deds_f_d_p(izx) = (eface(izx,j,k+1,l) - ecell_z(j,k+1,l))/(0.5D0*dx2(j))
      !deds_f_u_p(izy) = (ecell_z(j+1,k+1,l) - eface(izy,j+1,k,l))/(0.5D0*dy2(k+1))
      !deds_f_d_p(izy) = (eface(izy,j+1,k,l) - ecell_z(j+1,k,l))/(0.5D0*dy2(k))

      ! Grid corner electric fiedl gradients using mass_flux as upwinding !
      !deds_c_u(iyz) = 0.5D0*((1.0D0 + sign_x(j,k,l))*deds_f_u_m(iyz) + (1.0D0 - sign_x(j,k,l))*deds_f_u_p(iyz))
      !deds_c_u(izy) = 0.5D0*((1.0D0 + sign_x(j,k,l))*deds_f_u_m(izy) + (1.0D0 - sign_x(j,k,l))*deds_f_u_p(izy))
      !deds_c_d(iyz) = 0.5D0*((1.0D0 + sign_x(j,k,l))*deds_f_d_m(iyz) + (1.0D0 - sign_x(j,k,l))*deds_f_d_p(iyz))
      !deds_c_d(izy) = 0.5D0*((1.0D0 + sign_x(j,k,l))*deds_f_d_m(izy) + (1.0D0 - sign_x(j,k,l))*deds_f_d_p(izy))
      !deds_c_u(ixz) = 0.5D0*((1.0D0 + sign_y(j,k,l))*deds_f_u_m(ixz) + (1.0D0 - sign_y(j,k,l))*deds_f_u_p(ixz))
      !deds_c_u(izx) = 0.5D0*((1.0D0 + sign_y(j,k,l))*deds_f_u_m(izx) + (1.0D0 - sign_y(j,k,l))*deds_f_u_p(izx))
      !deds_c_d(ixz) = 0.5D0*((1.0D0 + sign_y(j,k,l))*deds_f_d_m(ixz) + (1.0D0 - sign_y(j,k,l))*deds_f_d_p(ixz))
      !deds_c_d(izx) = 0.5D0*((1.0D0 + sign_y(j,k,l))*deds_f_d_m(izx) + (1.0D0 - sign_y(j,k,l))*deds_f_d_p(izx))
      !deds_c_u(ixy) = 0.5D0*((1.0D0 + sign_z(j,k,l))*deds_f_u_m(ixy) + (1.0D0 - sign_z(j,k,l))*deds_f_u_p(ixy))
      !deds_c_u(iyx) = 0.5D0*((1.0D0 + sign_z(j,k,l))*deds_f_u_m(iyx) + (1.0D0 - sign_z(j,k,l))*deds_f_u_p(iyx))
      !deds_c_d(ixy) = 0.5D0*((1.0D0 + sign_z(j,k,l))*deds_f_d_m(ixy) + (1.0D0 - sign_z(j,k,l))*deds_f_d_p(ixy))
      !deds_c_d(iyx) = 0.5D0*((1.0D0 + sign_z(j,k,l))*deds_f_d_m(iyx) + (1.0D0 - sign_z(j,k,l))*deds_f_d_p(iyx))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Add emf !
      efield_x(j,k,l) = 0.5D0*(eface(ixy,j,k,l) + eface(ixy,j,k,l+1) + eface(ixz,j,k,l) + eface(ixz,j,k+1,l)) &
                      - 0.25D0*(ecell_x(j,k,l) + ecell_x(j,k,l+1) + ecell_x(j,k+1,l) + ecell_x(j,k+1,l+1))
                      !+ 0.125D0*((dy2(k)*deds_c_d(ixy) - dy2(k+1)*deds_c_u(ixy)) & 
                      !+ (dz2(l)*deds_c_d(ixz) - dz2(l+1)*deds_c_u(ixz)))
      efield_y(j,k,l) = 0.5D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1) + eface(iyz,j,k,l) + eface(iyz,j+1,k,l)) &
                      - 0.25D0*(ecell_y(j,k,l) + ecell_y(j,k,l+1) + ecell_y(j+1,k,l) + ecell_y(j+1,k,l+1))
                      !+ 0.125D0*((dx2(j)*deds_c_d(iyx) - dx2(j+1)*deds_c_u(iyx)) &
                      !+ (dz2(l)*deds_c_d(iyz) - dz2(l+1)*deds_c_u(iyz)))
      efield_z(j,k,l) = 0.5D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l) + eface(izy,j,k,l) + eface(izy,j+1,k,l)) &
                      - 0.25D0*(ecell_z(j,k,l) + ecell_z(j,k+1,l) + ecell_z(j+1,k,l) + ecell_z(j+1,k+1,l))
                      !+ 0.125D0*((dx2(j)*deds_c_d(izx) - dx2(j+1)*deds_c_u(izx)) &
                      !+ (dy2(k)*deds_c_d(izy) - dy2(k+1)*deds_c_u(izy)))

    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END PARALLEL
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'flux_ct = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE
    
end module MHD_module
