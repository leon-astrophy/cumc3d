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

! gravitational potential energy !
ALLOCATE (phi(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (phi_old(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! for poisson equation !
ALLOCATE (ajp1(1:nx_2))
ALLOCATE (ajm1(1:nx_2))
ALLOCATE (bkp1(1:nx_2,1:ny_2))
ALLOCATE (bkm1(1:nx_2,1:ny_2))
ALLOCATE (clp1(1:nx_2,1:ny_2,1:nz_2))
ALLOCATE (clm1(1:nx_2,1:ny_2,1:nz_2))
ALLOCATE (epsc(1:nx_2,1:ny_2,1:nz_2))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, nlines
REAL :: dr_temp 

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = './profile/grid.dat') 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Error message !
IF(nlines .ne. nx_2+7) THEN
  WRITE (*,*) 'number of grid faces from files', nlines
  WRITE (*,*) 'number of grid faces in the program', nx_2+6
  STOP 'inconsistent number of grid faces, exit'
END IF

! Read !
OPEN(UNIT=999, FILE = './profile/grid.dat', ACTION='READ')
DO i = -3, nx_2+3
	READ(999,*) xF2(i)
ENDDO
CLOSE(999)

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
DO l = nz_min_2 - 3, nz_part_2 + 3
  DO k = ny_min_2 - 3, ny_part_2 + 3
    DO j = 1, 3
      prim2(ivel2_x,1-j,k,l) = MIN(prim2(ivel2_x,1-j,k,l), 0.0D0)
      prim2(ivel2_x,nx_2+j,k,l) = MAX(prim2(ivel2_x,nx_2+j,k,l), 0.0D0)
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

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOMFLOOR
USE CUSTOM_DEF
USE MHD_MODULE
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Dummy variables
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: rho_min2, factor, diff, bfield, alven, rho_old, m_local

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

! assign threshold density !
rho_min2 = 1.1D0 * prim2_a(irho2)

!$OMP PARALLEL DO FIRSTPRIVATE(rho_min2) PRIVATE(diff, factor, bfield, alven, rho_old, m_local) COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(+:m_inj)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) FIRSTPRIVATE(rho_min2) PRIVATE(diff, factor, bfield, alven, rho_old, m_local) REDUCTION(+:m_inj)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2

      ! Standard !
      diff = prim2(irho2,j,k,l) - rho_min2
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      prim2(irho2:ivel2_z,j,k,l) = factor*prim2(irho2:ivel2_z,j,k,l) + (1.0D0 - factor)*prim2_a(irho2:ivel2_z)
      epsilon2(j,k,l) = factor*epsilon2(j,k,l) + (1.0D0 - factor)*eps2_a

      ! Check alven speed !
      bfield = SQRT(dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l)))
      alven = bfield/SQRT(prim2(irho2,j,k,l))

      ! Check !
      IF(alven >= max_alv) THEN
        rho_old = prim2(irho2,j,k,l)
        prim2(irho2,j,k,l) = (bfield/max_alv)**2
        m_local = (prim2(irho2,j,k,l) - rho_old)*vol2(j,k,l)
        m_inj = m_inj + m_local
      END IF

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
REAL*8 :: rho_min1, rho_min2, factor, diff 

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

! threshold !
rho_min2 = 1.1D0 * prim2_a(irho2)

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(factor, diff, dphidx, dphidy, dphidz) FIRSTPRIVATE(rho_min2)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(factor, diff, dphidx, dphidy, dphidz) FIRSTPRIVATE(rho_min2)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2

      ! Include only non-atmosphere !
			diff = prim2(irho2,j,k,l) - rho_min2
      factor = MAX(SIGN(1.0D0, diff), 0.0D0) 

      ! Gravitational potential of the matter !
      dphidx = first_derivative (x2(j-1), x2(j), x2(j+1), phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
      dphidy = first_derivative (y2(k-1), y2(k), y2(k+1), phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))
      dphidz = first_derivative (z2(l-1), z2(l), z2(l+1), phi(j,k,l-1), phi(j,k,l), phi(j,k,l+1))

      ! Add black hole force !
      dphidx = dphidx + m_bh/((x2(j) - 2.0d0*m_bh)*(x2(j) - 2.0d0*m_bh))

      ! Add them to the source term !
      sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) - factor*prim2(irho2,j,k,l)*dphidx
      sc2(ivel2_y,j,k,l) = sc2(ivel2_y,j,k,l) - factor*prim2(irho2,j,k,l)*dphidy/x2(j)
      sc2(ivel2_z,j,k,l) = sc2(ivel2_z,j,k,l) - factor*prim2(irho2,j,k,l)*dphidz/x2(j)/sin2(k)
      sc2(itau2,j,k,l) = sc2(itau2,j,k,l) - factor*prim2(irho2,j,k,l)* &
                        (prim2(ivel2_x,j,k,l)*dphidx + prim2(ivel2_y,j,k,l)*dphidy/x2(j) + &
                         prim2(ivel2_z,j,k,l)*dphidz/x2(j)/sin2(k))

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
  !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) REDUCTION(+:mdot)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT) REDUCTION(+:mdot)
  DO l = nz_min_2 , nz_part_2 
    DO k = ny_min_2 , ny_part_2 
      mdot = mdot + (-mflux_x(k,l))*x2(0)*x2(0)*dcos2(k)*dz2(l)
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
IF(2.0d0*m_bh >= xF2(0)) THEN
  STOP 'Black hole event horizon exceeds computational grid'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update gravitational potentials !
IF (p_in == 0 .OR. MOD(n_step, n_pot) == 0) THEN

  ! special treatment for initial model !
  IF(p_in == 0) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, give a guessing potential !
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
    DO l = nz_min_2-1, nz_part_2+1
      DO k = ny_min_2-1, ny_part_2+1
        DO j = nx_min_2-1, nx_part_2+1
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
    DO l = nz_min_2, nz_part_2
      DO k = ny_min_2, ny_part_2
        DO j = nx_min_2, nx_part_2
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
    DO l = nz_min_2, nz_part_2
      DO k = ny_min_2, ny_part_2
        DO j = nx_min_2, nx_part_2
          IF ((-1)**(j+k+l)>0) THEN	
            diff = prim2(irho2,j,k,l) - prim2_a(irho2)
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim2(irho2,j,k,l)
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
    DO l = nz_min_2, nz_part_2
      DO k = ny_min_2, ny_part_2
        DO j = nx_min_2, nx_part_2
          IF ((-1)**(j+k+l)<0) THEN
            diff = prim2(irho2,j,k,l) - prim2_a(irho2)
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim2(irho2,j,k,l)
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
    DO l = nz_min_2, nz_part_2
      DO k = ny_min_2, ny_part_2
        DO j = nx_min_2, nx_part_2
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
    DO l = nz_min_2, nz_part_2
      DO k = ny_min_2, ny_part_2
        phi(0,k,l) = phi(1,k,l)
        phi(nx_part_2+1,k,l) = -1.0d0
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO l = nz_min_2, nz_part_2
      DO j = nx_min_2, nx_part_2
        phi(j,0,l) = phi(j,1,l)
        phi(j,ny_part_2+1,l) = phi(j,ny_part_2,l)
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END DO
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        phi(j,k,0) = phi(j,k,1)
        phi(j,k,nz_part_2+1) = phi(j,k,nz_part_2)
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

END SUBROUTINE
