!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! atmospheric values !
ALLOCATE (prim_a(imin:imax))

! mass fluxes !
ALLOCATE (mflux_x(-2:ny+3,-2:nz+3))

! gravitational potential energy !
ALLOCATE (phi(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (phi_old(-2:nx+3,-2:ny+3,-2:nz+3))

! for poisson equation !
ALLOCATE (ajp1(1:nx))
ALLOCATE (ajm1(1:nx))
ALLOCATE (bkp1(1:nx,1:ny))
ALLOCATE (bkm1(1:nx,1:ny))
ALLOCATE (clp1(1:nx,1:ny,1:nz))
ALLOCATE (clm1(1:nx,1:ny,1:nz))
ALLOCATE (epsc(1:nx,1:ny,1:nz))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(prim_a, mflux_x, phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(prim_a, mflux_x, phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, nlines
INTEGER :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = './profile/r_grid.dat') 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Error message !
IF(nlines .ne. nx+7) THEN
  WRITE (*,*) 'number of r-grid faces from files', nlines-7
  WRITE (*,*) 'number of r-grid faces in the program', nx
  STOP 'inconsistent number of r-grid faces, exit'
END IF

! Read !
OPEN(UNIT=999, FILE = './profile/r_grid.dat', ACTION='READ')
DO i = -3, nx+3
	READ(999,*) xF(i)
ENDDO
CLOSE(999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = './profile/th_grid.dat') 
DO 
  READ (999,*, END=20) 
  nlines = nlines + 1 
END DO 
20 CLOSE (999) 

! Error message !
IF(nlines .ne. ny+7) THEN
  WRITE (*,*) 'number of theta-grid faces from files', nlines-7
  WRITE (*,*) 'number of theta-grid faces in the program', ny
  STOP 'inconsistent number of theta-grid faces, exit'
END IF

! Read !
OPEN(UNIT=999, FILE = './profile/th_grid.dat', ACTION='READ')
DO i = -3, ny+3
	READ(999,*) yF(i)
ENDDO
CLOSE(999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read the number of lines in the file !
IF(n_dim > 2) THEN
  nlines = 0 
  OPEN (999, file = './profile/phi_grid.dat') 
  DO 
    READ (999,*, END=30) 
    nlines = nlines + 1 
  END DO 
  30 CLOSE (999) 

  ! Error message !
  IF(nlines .ne. nz+7) THEN
    WRITE (*,*) 'number of phi-grid faces from files', nlines-7
    WRITE (*,*) 'number of phi-grid faces in the program', nz
    STOP 'inconsistent number of phi-grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=999, FILE = './profile/phi_grid.dat', ACTION='READ')
  DO i = -3, nz+3
    READ(999,*) zF(i)
  ENDDO
  CLOSE(999)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_X
USE DEFINITION
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)   
DO l = - 2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = MIN(prim(ivx,1-j,k,l), 0.0D0)
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO               
ENDDO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom boundary = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Y
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Z
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_X
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do it by direction !
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
DO l = 1, nz
  DO k = 1, ny
      mflux_x (k,l) = flux (irho,0,k,l)
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Y
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Z
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE CUSTOM_DEF
USE MHD_MODULE
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Dummy variables
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: factor, diff, bfield, alven, rho_old, m_local

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO PRIVATE(diff, factor, bfield, alven, rho_old, m_local) COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(+:m_inj)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff, factor, bfield, alven, rho_old, m_local) REDUCTION(+:m_inj)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Standard !
      diff = prim(irho,j,k,l) - prim_a(irho)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      prim(irho:ivz,j,k,l) = factor*prim(irho:ivz,j,k,l) + (1.0D0 - factor)*prim_a(irho:ivz)
      epsilon(j,k,l) = factor*epsilon(j,k,l) + (1.0D0 - factor)*eps_a
      IF(epsilon(j,k,l) < 0.0d0) epsilon(j,k,l) = eps_a

      ! Check alven speed !
      !bfield = SQRT(dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l)))
      !alven = bfield/SQRT(prim(irho,j,k,l))

      ! Check !
      !IF(alven >= max_alv) THEN
      !  rho_old = prim(irho,j,k,l)
      !  prim(irho,j,k,l) = (bfield/max_alv)**2
      !  m_local = (prim(irho,j,k,l) - rho_old)*vol2(j,k,l)
      !  m_inj = m_inj + m_local
      !END IF

    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: dphidx, dphidy, dphidz

! Threshold for atmosphere density
REAL*8 :: factor, diff

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add black hole gravity !

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Standard !
      diff = prim(irho,j,k,l) - prim_a(irho)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)

      ! Gravitational potential of the matter !
      dphidx = first_derivative (x(j-1), x(j), x(j+1), phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
      dphidy = first_derivative (y(k-1), y(k), y(k+1), phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))
      dphidz = first_derivative (z(l-1), z(l), z(l+1), phi(j,k,l-1), phi(j,k,l), phi(j,k,l+1))

      ! Add black hole force !
      dphidx = dphidx + m_bh/((x(j) - r_sh)*(x(j) - r_sh))

      ! Add them to the source term !
      sc(ivx,j,k,l) = sc(ivx,j,k,l) - factor*prim(irho,j,k,l)*dphidx
      sc(ivy,j,k,l) = sc(ivy,j,k,l) - factor*prim(irho,j,k,l)*dphidy/x(j)
      sc(ivz,j,k,l) = sc(ivz,j,k,l) - factor*prim(irho,j,k,l)*dphidz/x(j)/sine(k)
      sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)* &
                        (prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/x(j) + &
                         prim(ivz,j,k,l)*dphidz/x(j)/sine(k))

    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom source = ', REAL(time_end - time_start) / rate
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function first_derivative (xm1, xc, xp1, fm1, fc, fp1)
	!$acc routine seq
	implicit none
	REAL*8 :: xm1, xc, xp1, fm1, fc, fp1, h1, h2
  h2 = xp1 - xc
  h1 = xc - xm1
	first_derivative = ((fp1-fc)*h1*h1+(fc-fm1)*h2*h2)/(h1*h2*(h1+h2))
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
USE CUSTOM_DEF 
USE MHD_MODULE 
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER, INTENT (IN) :: p_in

! Integer !
INTEGER :: j, k, l, n

! Real, mdot !
REAL*8 :: mdot

! For poisson solver !
REAL*8 :: abserror, rhs

! Density threshold !
REAL*8 :: rho_in, factor, diff

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update black hole masses !
IF (p_in > 0) THEN

  ! Find the mass accretion rate !
  mdot = 0.0d0
  bflux = 0.0d0
  !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) REDUCTION(+:mdot)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT) REDUCTION(+:mdot)
  DO l = 1, nz 
    DO k = 1, ny 
      mdot = mdot + (-mflux_x(k,l))*xF(0)*xF(0)*dcose(k)*dz(l)
      bflux = bflux + prim(ibx,0,k,l)*xF(0)*xF(0)*dcose(k)*dz(l)
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
  mdot = MAX(mdot, 0.0d0)
  
  ! Update black hole masses !
  IF (p_in == 1) THEN
    mbh_old = m_bh
    m_bh = mbh_old + dt * mdot
  ELSEIF (p_in == 2) THEN
    m_bh = rk20 * mbh_old + rk21 * m_bh + rk22 * dt * mdot
  ELSEIF (p_in == 3) THEN
    m_bh = rk30 * mbh_old + rk31 * m_bh + rk32 * dt * mdot
  END IF
  
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stop condition !
r_sh = 2.0d0*m_bh/(clight*vel2code)**2
IF(r_sh >= xF(0)) THEN
  STOP 'Black hole event horizon exceeds computational grid'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update gravitational potentials !
IF (p_in == 0 .OR. (p_in == 3 .AND. MOD(n_step, n_pot) == 0)) THEN

  ! special treatment for initial model !
  IF(p_in == 0) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, give a guessing potential !
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
    DO l = 0, nz + 1
      DO k = 0, ny + 1
        DO j = 0, nx + 1
          phi(j,k,l) = 0.0d0
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END PARALLEL DO
    
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calucaltes potential by RBSOR
  DO n = 1, relax_max
    
    !$OMP PARALLEL PRIVATE(diff,rhs,rho_in,factor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Back up potential !
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          phi_old(j,k,l) = phi(j,k,l)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set error !
    !$OMP SINGLE
    !$ACC SERIAL
	  abserror = 1.0D-50
	  !$ACC END SERIAL
    !$OMP END SINGLE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Red chess !
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF ((-1)**(j+k+l)>0) THEN	
            diff = prim(irho,j,k,l) - prim_a(irho)
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim(irho,j,k,l)
            rhs = (4.0d0*pi*rho_in - & 
								(ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
								 bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
								 clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
					  phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
          ELSE 
            CYCLE
          END IF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Black chess !
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF ((-1)**(j+k+l)<0) THEN
            diff = prim(irho,j,k,l) - prim_a(irho)
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim(irho,j,k,l)
            rhs = (4.0d0*pi*rho_in - & 
								(ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
								 bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
								 clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
					  phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
          ELSE 
            CYCLE
          END IF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Look for maximum abserror !
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) Reduction(MAX:abserror)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) Reduction(MAX:abserror)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          abserror = max(abserror, abs((phi(j,k,l) - phi_old(j,k,l)) / phi_old(j,k,l)))
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Boundary conditions !
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO l = 1, nz
      DO k = 1, ny
        phi(0,k,l) = phi(1,k,l)
        phi(nx+1,k,l) = -1.0d0
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO l = 1, nz
      DO j = 1, nx
        phi(j,0,l) = phi(j,1,l)
        phi(j,ny+1,l) = phi(j,ny,l)
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO k = 1, ny
      DO j = 1, nx
        phi(j,k,0) = phi(j,k,nz)
        phi(j,k,nz+1) = phi(j,k,1)
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$OMP END PARALLEL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Debug and exit !
	  !WRITE (*,*) n, abserror
    IF(abserror <= tolerance) EXIT 

  END DO
  !write (*,*) n, abserror
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Stop condition !
  IF(n == relax_max) THEN
    WRITE (*,*) n, relax_max
    STOP 'Convergence error in poisson solver'
  END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom update = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
IMPLICIT NONE

! Open !
OPEN (UNIT = 689, FILE = './outfile/star_weno_mbh.dat')
OPEN (UNIT = 777, FILE = './outfile/star_weno_minj.dat')
OPEN (UNIT = 888, FILE = './outfile/star_weno_bflux.dat')

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_ANALYSIS
USE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! WRITE !
WRITE (689, *) global_time, m_bh
WRITE (777, *) global_time, m_inj
WRITE (888, *) global_time, bflux

END SUBROUTINE
