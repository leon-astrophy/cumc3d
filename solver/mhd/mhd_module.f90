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

! Indices for edge and cell-center electric fields !
INTEGER :: iex, iey, iez

! cell and face-centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: eface

! cell centered magnetic fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: bcell

! cell centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecell

! electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: efield_x, efield_y, efield_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! signed flux !
!REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: sign_x, sign_y, sign_z
! electric field gradients at interface !
! the structure is i(e-field component)(differentiate axis) !
! e.g. iyx means the dey/dx !
!REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_f_d_p, deds_f_u_p
!REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_f_d_m, deds_f_u_m
!REAL*8, ALLOCATABLE, DIMENSION(:) :: deds_c_d, deds_c_u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine allocates the array necessary for
! the calculation of MHD
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BUILD_MHD
use definition
implicit none

! Allocate !
ALLOCATE(eface(iyx:iyz,-2:nx+3,-2:ny+3,-2:nz+3))

! Allocate !
ALLOCATE(bcell(ibx:ibz,-2:nx+3,-2:ny+3,-2:nz+3))

! Allocate !
ALLOCATE(ecell(iex:iez,-2:nx+3,-2:ny+3,-2:nz+3))

! Allocate !
ALLOCATE(efield_x(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(efield_y(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(efield_z(-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate !
!ALLOCATE(sign_x(-2:nx+3,-2:ny+3,-2:nz+3))
!ALLOCATE(sign_y(-2:nx+3,-2:ny+3,-2:nz+3))
!ALLOCATE(sign_z(-2:nx+3,-2:ny+3,-2:nz+3))
! Allocate !
!ALLOCATE(deds_c_u(iyx:iyz))
!ALLOCATE(deds_c_d(iyx:iyz))
!ALLOCATE(deds_f_d_p(iyx:iyz))
!ALLOCATE(deds_f_u_p(iyx:iyz))
!ALLOCATE(deds_f_d_m(iyx:iyz))
!ALLOCATE(deds_f_u_m(iyx:iyz))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine build_MHD

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

! local divergence !
REAL*8 :: divb

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
maxDivB = 0.0d0

! Loop over the domain !
IF(coordinate_flag == 0) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        divb = (prim(ibx,j,k,l) - prim(ibx,j-1,k,l))/(dx(j)) &

             + (prim(iby,j,k,l) - prim(iby,j,k-1,l))/(dy(k)) &

             + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))/(dz(l))
        maxDivB = MAX(maxDivB, ABS(divb))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(coordinate_flag == 1) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        divb = (xF(j)*prim(ibx,j,k,l) - xF(j-1)*prim(ibx,j-1,k,l))/(x(j)*dx(j)) &

             + (prim(iby,j,k,l) - prim(iby,j,k-1,l))/(x(j)*dy(k)) &

             + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))/(dz(l))

        maxDivB = MAX(maxDivB, ABS(divb))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(coordinate_flag == 2) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(max:maxDivB)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        divb = (xF(j)*xF(j)*prim(ibx,j,k,l) - xF(j-1)*xF(j-1)*prim(ibx,j-1,k,l))/(dx_cb(j)/3.0d0) &

             + (sinf(k)*prim(iby,j,k,l) - sinf(k-1)*prim(iby,j,k-1,l))*(x(j)*dx(j))/(dx_cb(j)*dcose(k)/3.0d0) &
             
             + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))*(x(j)*dx(j)*dy(k))/(dx_cb(j)*dcose(k)*dz(l)/3.0d0)
             
        maxDivB = MAX(maxDivB, ABS(divb))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
END IF 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'finddivb = ', REAL(time_end - time_start) / rate
#endif

! output !
WRITE (999, *) global_time, maxDivB

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
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do it by direction !
IF(dir_in == x_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 0, nx
        eface (iyx,j,k,l) = flux (ibz,j,k,l)
        eface (izx,j,k,l) = - flux (iby,j,k,l)
        !sign_x (j,k,l) = SIGN(1.0d0, flux (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(dir_in == y_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
  DO l = 1, nz
    DO k = 0, ny
      DO j = 1, nx
        eface (ixy,j,k,l) = - flux (ibz,j,k,l)
        eface (izy,j,k,l) = flux (ibx,j,k,l)
        !sign_y (j,k,l) = SIGN(1.0d0, flux (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
ELSEIF(dir_in == z_dir) THEN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) 
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
  DO l = 0, nz
    DO k = 1, ny
      DO j = 1, nx
        eface (ixz,j,k,l) = flux (iby,j,k,l)
        eface (iyz,j,k,l) = - flux (ibx,j,k,l)
        !sign_z (j,k,l) = SIGN(1.0d0, flux (irho2,j,k,l))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
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

! Real !
REAL :: sum1, sum2 

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find cell-centered electric fields !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      ecell (iex,j,k,l) = prim(ivz,j,k,l)*bcell(iby,j,k,l) - prim(ivy,j,k,l)*bcell(ibz,j,k,l)
      ecell (iey,j,k,l) = prim(ivx,j,k,l)*bcell(ibz,j,k,l) - prim(ivz,j,k,l)*bcell(ibx,j,k,l)
      ecell (iez,j,k,l) = prim(ivy,j,k,l)*bcell(ibx,j,k,l) - prim(ivx,j,k,l)*bcell(iby,j,k,l)
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign cell-interface emf for 1D/2D problem !
IF(n_dim == 1) THEN
  !$OMP DO COLLAPSE(1) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(1) DEFAULT(PRESENT) 
  DO j = 1, nx 
    eface (ixy,j,1,1) = ecell (iex,j,1,1)
    eface (izy,j,1,1) = ecell (iez,j,1,1)
    eface (ixy,j,0,1) = ecell (iex,j,1,1)
    eface (izy,j,0,1) = ecell (iez,j,1,1)
    eface (ixz,j,1,1) = ecell (iex,j,1,1)
    eface (iyz,j,1,1) = ecell (iey,j,1,1)
    eface (ixz,j,1,0) = ecell (iex,j,1,1)
    eface (iyz,j,1,0) = ecell (iey,j,1,1)
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(n_dim == 2) THEN
  !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT) 
  DO k = 1, ny 
    DO j = 1, nx 
      eface (ixz,j,k,1) = ecell (iex,j,k,1)
      eface (iyz,j,k,1) = ecell (iey,j,k,1)
      eface (ixz,j,k,0) = ecell (iex,j,k,1)
      eface (iyz,j,k,0) = ecell (iey,j,k,1)      
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions for electric field !
CALL BOUNDARY_EFIELD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! constrained transport !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! First, get electric field gradient at cell interface !
      !deds_f_u_m(ixy) = (ecell(iex,j,k+1,l) - eface(ixy,j,k,l))/(0.5D0*dy(k+1))
      !deds_f_d_m(ixy) = (eface(ixy,j,k,l) - ecell(iex,j,k,l))/(0.5D0*dy(k))
      !deds_f_u_m(ixz) = (ecell(iex,j,k,l+1) - eface(ixz,j,k,l))/(0.5D0*dz(l+1))
      !deds_f_d_m(ixz) = (eface(ixz,j,k,l) - ecell(iex,j,k,l))/(0.5D0*dz(l))
      !deds_f_u_m(iyx) = (ecell(iey,j+1,k,l) - eface(iyx,j,k,l))/(0.5D0*dx(j+1))
      !deds_f_d_m(iyx) = (eface(iyx,j,k,l) - ecell(iey,j,k,l))/(0.5D0*dx(j))
      !deds_f_u_m(iyz) = (ecell(iey,j,k,l+1) - eface(iyz,j,k,l))/(0.5D0*dz(l+1))
      !deds_f_d_m(iyz) = (eface(iyz,j,k,l) - ecell(iey,j,k,l))/(0.5D0*dz(l))
      !deds_f_u_m(izx) = (ecell(iez,j+1,k,l) - eface(izx,j,k,l))/(0.5D0*dx(j+1))
      !deds_f_d_m(izx) = (eface(izx,j,k,l) - ecell(iez,j,k,l))/(0.5D0*dx(j))
      !deds_f_u_m(izy) = (ecell(iez,j,k+1,l) - eface(izy,j,k,l))/(0.5D0*dy(k+1))
      !deds_f_d_m(izy) = (eface(izy,j,k,l) - ecell(iez,j,k,l))/(0.5D0*dy(k))
      ! Also get the graident at an upper grid !
      !deds_f_u_p(ixy) = (ecell(iex,j,k+1,l+1) - eface(ixy,j,k,l+1))/(0.5D0*dy(k+1))
      !deds_f_d_p(ixy) = (eface(ixy,j,k,l+1) - ecell(iex,j,k,l+1))/(0.5D0*dy(k))
      !deds_f_u_p(ixz) = (ecell(iex,j,k+1,l+1) - eface(ixz,j,k+1,l))/(0.5D0*dz(l+1))
      !deds_f_d_p(ixz) = (eface(ixz,j,k+1,l) - ecell(iex,j,k+1,l))/(0.5D0*dz(l))
      !deds_f_u_p(iyx) = (ecell(iey,j+1,k,l+1) - eface(iyx,j,k,l+1))/(0.5D0*dx(j+1))
      !deds_f_d_p(iyx) = (eface(iyx,j,k,l+1) - ecell(iey,j,k,l+1))/(0.5D0*dx(j))
      !deds_f_u_p(iyz) = (ecell(iey,j+1,k,l+1) - eface(iyz,j+1,k,l))/(0.5D0*dz(l+1))
      !deds_f_d_p(iyz) = (eface(iyz,j+1,k,l) - ecell(iey,j+1,k,l))/(0.5D0*dz(l))
      !deds_f_u_p(izx) = (ecell(iez,j+1,k+1,l) - eface(izx,j,k+1,l))/(0.5D0*dx(j+1))
      !deds_f_d_p(izx) = (eface(izx,j,k+1,l) - ecell(iez,j,k+1,l))/(0.5D0*dx(j))
      !deds_f_u_p(izy) = (ecell(iez,j+1,k+1,l) - eface(izy,j+1,k,l))/(0.5D0*dy(k+1))
      !deds_f_d_p(izy) = (eface(izy,j+1,k,l) - ecell(iez,j+1,k,l))/(0.5D0*dy(k))
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
      efield_x(j,k,l) = 0.50D0*(eface(ixy,j,k,l) + eface(ixy,j,k,l+1) + eface(ixz,j,k,l) + eface(ixz,j,k+1,l)) &
                      - 0.25D0*(ecell(iex,j,k,l) + ecell(iex,j,k,l+1) + ecell(iex,j,k+1,l) + ecell(iex,j,k+1,l+1))
                      !+ 0.125D0*((dy(k)*deds_c_d(ixy) - dy(k+1)*deds_c_u(ixy)) & 
                      !+ (dz(l)*deds_c_d(ixz) - dz(l+1)*deds_c_u(ixz)))
      efield_y(j,k,l) = 0.50D0*(eface(iyx,j,k,l) + eface(iyx,j,k,l+1) + eface(iyz,j,k,l) + eface(iyz,j+1,k,l)) &
                      - 0.25D0*(ecell(iey,j,k,l) + ecell(iey,j,k,l+1) + ecell(iey,j+1,k,l) + ecell(iey,j+1,k,l+1))
                      !+ 0.125D0*((dx(j)*deds_c_d(iyx) - dx(j+1)*deds_c_u(iyx)) &
                      !+ (dz(l)*deds_c_d(iyz) - dz(l+1)*deds_c_u(iyz)))
      efield_z(j,k,l) = 0.50D0*(eface(izx,j,k,l) + eface(izx,j,k+1,l) + eface(izy,j,k,l) + eface(izy,j+1,k,l)) &
                      - 0.25D0*(ecell(iez,j,k,l) + ecell(iez,j,k+1,l) + ecell(iez,j+1,k,l) + ecell(iez,j+1,k+1,l))
                      !+ 0.125D0*((dx(j)*deds_c_d(izx) - dx(j+1)*deds_c_u(izx)) &
                      !+ (dy(k)*deds_c_d(izy) - dy(k+1)*deds_c_u(izy)))

    
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fix the electric field at coordinate axis !
#ifdef FIXPOLE
IF(coordinate_flag == 2) THEN
  ! First, average across the phi direction !
  DO j = 1, nx
    sum1 = 0.0d0
    sum2 = 0.0d0
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:sum1, sum2)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR DEFAULT(PRESENT) REDUCTION(+:sum1, sum2)
    DO l = 1, nz
      sum1 = sum1 + efield_x(j,0,l)
      sum2 = sum2 + efield_x(j,ny,l)
    END DO 
    !$ACC END PARALLEL
    !$OMP END DO 
    !$OMP DO SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR DEFAULT(PRESENT)
    DO l = -2, nz + 3
      efield_x(j,0,l) = sum1/DBLE(nz)
      efield_x(j,ny,l) = sum2/DBLE(nz)
    END DO
  !$ACC END PARALLEL
  !$OMP END DO 
  END DO
END IF
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'flux_ct = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign boundary conditions for electric fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_EFIELD
use definition
implicit none

! Dummy variables
INTEGER :: i, j, k, l

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
        ecell(iex:iez,1-j,k,l) = ecell(iex:iez,nx+1-j,k,l)
        eface(ixz,1-j,k,l) = eface(ixz,nx+1-j,k,l)
        eface(ixy,1-j,k,l) = eface(ixy,nx+1-j,k,l)
        eface(iyz,1-j,k,l) = eface(iyz,nx+1-j,k,l)
        eface(izy,1-j,k,l) = eface(izy,nx+1-j,k,l)
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
        ecell(iex:iez,1-j,k,l) = ecell(iex:iez,1,k,l)
        eface(ixz,1-j,k,l) = eface(ixz,1,k,l)
        eface(ixy,1-j,k,l) = eface(ixy,1,k,l)
        eface(iyz,1-j,k,l) = eface(iyz,1,k,l)
        eface(izy,1-j,k,l) = eface(izy,1,k,l)
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
        ecell(iex:iez,1-j,k,l) = bfac_xin(iex:iez) * ecell(iex:iez,j,k,l)
        eface(ixz,1-j,k,l) = bfac_xin(iex) * eface(ixz,j,k,l)
        eface(ixy,1-j,k,l) = bfac_xin(iex) * eface(ixy,j,k,l)
        eface(iyz,1-j,k,l) = bfac_xin(iey) * eface(iyz,j,k,l)
        eface(izy,1-j,k,l) = bfac_xin(iez) * eface(izy,j,k,l)
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
        ecell(iex:iez,nx+j,k,l) = ecell(iex:iez,j,k,l)
        eface(ixz,nx+j,k,l) = eface(ixz,j,k,l)
        eface(ixy,nx+j,k,l) = eface(ixy,j,k,l)
        eface(iyz,nx+j,k,l) = eface(iyz,j,k,l)
        eface(izy,nx+j,k,l) = eface(izy,j,k,l)
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
        ecell(iex:iez,nx+j,k,l) = ecell(iex:iez,nx,k,l)
        eface(ixz,nx+j,k,l) = eface(ixz,nx,k,l)
        eface(ixy,nx+j,k,l) = eface(ixy,nx,k,l)
        eface(iyz,nx+j,k,l) = eface(iyz,nx,k,l)
        eface(izy,nx+j,k,l) = eface(izy,nx,k,l)
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
        ecell(iex:iez,nx+j,k,l) = bfac_xout(iex:iez) * ecell(iex:iez,nx+1-j,k,l)
        eface(ixz,nx+j,k,l) = bfac_xout(iex) * eface(ixz,nx+1-j,k,l)
        eface(ixy,nx+j,k,l) = bfac_xout(iex) * eface(ixy,nx+1-j,k,l)
        eface(iyz,nx+j,k,l) = bfac_xout(iey) * eface(iyz,nx+1-j,k,l)
        eface(izy,nx+j,k,l) = bfac_xout(iez) * eface(izy,nx+1-j,k,l)
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
        ecell(iex:iez,j,1-k,l) = ecell(iex:iez,j,ny+1-k,l)            
        eface(ixz,j,1-k,l) = eface(ixz,j,ny+1-k,l)
        eface(iyx,j,1-k,l) = eface(iyx,j,ny+1-k,l) 
        eface(iyz,j,1-k,l) = eface(iyz,j,ny+1-k,l) 
        eface(izx,j,1-k,l) = eface(izx,j,ny+1-k,l)               
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
        ecell(iex:iez,j,1-k,l) = ecell(iex:iez,j,1,l)
        eface(ixz,j,1-k,l) = eface(ixz,j,1,l)
        eface(iyx,j,1-k,l) = eface(iyx,j,1,l)       
        eface(iyz,j,1-k,l) = eface(iyz,j,1,l)   
        eface(izx,j,1-k,l) = eface(izx,j,1,l)      
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
        ecell(iex:iez,j,1-k,l) = bfac_yin(iex:iez) * ecell(iex:iez,j,k,l)
        eface(ixz,j,1-k,l) = bfac_yin(iex) * eface(ixz,j,k,l)
        eface(iyx,j,1-k,l) = bfac_yin(iey) * eface(iyx,j,k,l)   
        eface(iyz,j,1-k,l) = bfac_yin(iey) * eface(iyz,j,k,l) 
        eface(izx,j,1-k,l) = bfac_yin(iez) * eface(izx,j,k,l)       
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
        ecell(iex:iez,j,ny+k,l) = ecell(iex:iez,j,k,l)
        eface(ixz,j,ny+k,l) = eface(ixz,j,k,l)
        eface(iyx,j,ny+k,l) = eface(iyx,j,k,l)  
        eface(iyz,j,ny+k,l) = eface(iyz,j,k,l)
        eface(izx,j,ny+k,l) = eface(izx,j,k,l)      
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
        ecell(iex:iez,j,ny+k,l) = ecell(iex:iez,j,ny,l)
        eface(ixz,j,ny+k,l) = eface(ixz,j,ny,l)
        eface(iyx,j,ny+k,l) = eface(iyx,j,ny,l)
        eface(iyz,j,ny+k,l) = eface(iyz,j,ny,l)
        eface(izx,j,ny+k,l) = eface(izx,j,ny,l)       
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
        ecell(iex:iez,j,ny+k,l) = bfac_yout(iex:iez) * ecell(iex:iez,j,ny+1-k,l)
        eface(ixz,j,ny+k,l) = bfac_yout(iex) * eface(ixz,j,ny+1-k,l)
        eface(iyx,j,ny+k,l) = bfac_yout(iey) * eface(iyx,j,ny+1-k,l)
        eface(iyz,j,ny+k,l) = bfac_yout(iey) * eface(iyz,j,ny+1-k,l)
        eface(izx,j,ny+k,l) = bfac_yout(iez) * eface(izx,j,ny+1-k,l)       
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
        ecell(iex:iez,j,k,1-l) = ecell(iex:iez,j,k,nz+1-l)               
        eface(ixy,j,k,1-l) = eface(ixy,j,k,nz+1-l)
        eface(iyx,j,k,1-l) = eface(iyx,j,k,nz+1-l)    
        eface(izx,j,k,1-l) = eface(izx,j,k,nz+1-l)
        eface(izy,j,k,1-l) = eface(izy,j,k,nz+1-l)               
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
        ecell(iex:iez,j,k,1-l) = ecell(iex:iez,j,k,1)
        eface(ixy,j,k,1-l) = eface(ixy,j,k,1)
        eface(iyx,j,k,1-l) = eface(iyx,j,k,1) 
        eface(izx,j,k,1-l) = eface(izx,j,k,1)
        eface(izy,j,k,1-l) = eface(izy,j,k,1)                 
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
        ecell(iex:iez,j,k,1-l) = bfac_zin(iex:iez) * ecell(iex:iez,j,k,l)
        eface(ixy,j,k,1-l) = bfac_zin(iex) * eface(ixy,j,k,l)
        eface(iyx,j,k,1-l) = bfac_zin(iey) * eface(iyx,j,k,l)   
        eface(izx,j,k,1-l) = bfac_zin(iez) * eface(izx,j,k,l)   
        eface(izy,j,k,1-l) = bfac_zin(iez) * eface(izy,j,k,l)              
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
        ecell(iex:iez,j,k,nz+l) = ecell(iex:iez,j,k,l)
        eface(ixy,j,k,nz+l) = eface(ixy,j,k,l)
        eface(iyx,j,k,nz+l) = eface(iyx,j,k,l)       
        eface(izx,j,k,nz+l) = eface(izx,j,k,l)       
        eface(izy,j,k,nz+l) = eface(izy,j,k,l)            
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
        ecell(iex:iez,j,k,nz+l) = ecell(iex:iez,j,k,nz)
        eface(ixy,j,k,nz+l) = eface(ixy,j,k,nz)
        eface(iyx,j,k,nz+l) = eface(iyx,j,k,nz)         
        eface(izx,j,k,nz+l) = eface(izx,j,k,nz)         
        eface(izy,j,k,nz+l) = eface(izy,j,k,nz)         
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
        ecell(iex:iez,j,k,nz+l) = bfac_zout(iex:iez) * ecell(iex:iez,j,k,nz+1-l) 
        eface(ixy,j,k,nz+l) = bfac_zout(iex) * eface(ixy,j,k,nz+1-l)
        eface(iyx,j,k,nz+l) = bfac_zout(iey) * eface(iyx,j,k,nz+1-l)  
        eface(izx,j,k,nz+l) = bfac_zout(iez) * eface(izx,j,k,nz+1-l) 
        eface(izy,j,k,nz+l) = bfac_zout(iez) * eface(izy,j,k,nz+1-l)      
      END DO
    END DO               
  ENDDO 
  !$OMP END DO  

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end module MHD_module
