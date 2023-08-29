!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l, nlines

! Distances !
REAL*8 :: rmax
REAL*8 :: s_eq
REAL*8 :: s_core

! Angular velocity !
REAL*8 :: omega

! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b

! Random number !
REAL*8 :: rand

! Vector potential !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: a_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !

! Allocate
Allocate(a_phi(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Poisson interpolation coefficient !
call get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read parameter !
OPEN(UNIT=999, FILE = './profile/par.dat', ACTION='READ')
DO j = 1, 1
	READ(999,*) ggas2, m_bh
ENDDO
CLOSE(999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = './profile/hydro.dat') 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Error message !
IF(nlines .ne. nx_2+6) THEN
  WRITE (*,*) 'number of simulation grids from files', nlines
  WRITE (*,*) 'number of simulation grids in the program', nx_2+6
  STOP 'inconsistent number of simulation grids, exit'
END IF

! Read !
OPEN(UNIT=999, FILE = './profile/hydro.dat', ACTION='READ')
DO j = -2, nx_2 + 3
	READ(999,*) prim2(irho2,j,1,1), prim2(itau2,j,1,1)
ENDDO
CLOSE(999)

! Assign density profile !
DO j = -2, nx_2 + 3
  prim2(irho2,j,:,:) = prim2(irho2,j,1,1)
  prim2(itau2,j,:,:) = prim2(itau2,j,1,1)
END DO

! Assign internal energy #
epsilon2(:,:,:) = prim2(itau2,:,:,:)/prim2(irho2,:,:,:)/(ggas2 - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose rotation rules !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core radius !
DO j = 1, nx_2
  IF(prim2(irho2,j-1,1,1) <= prim2(irho2,nx_2,1,1) .AND. &
     prim2(irho2,j,1,1) > prim2(irho2,nx_2,1,1)) THEN
    s_core = x2(j)
    EXIT
  END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rigid rotations !
IF(rotation_rule == 1) THEN

  ! Find outer most radius !
  rmax = MAXVAL(x2)

  ! Find angular velocity !
  omega = v_crit/rmax

  ! Assign angular velocity !
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2

        ! polar radius !
        s_eq = x2(j)*sin2(k)

        ! Select based on conditions !
        IF(prim2(irho2,j,k,l) > prim2(irho2,nx_2,1,1)) THEN
          prim2(ivel2_z,j,k,l) = omega*s_eq
        ELSE
          prim2(ivel2_z,j,k,l) = 0.0d0                   
        END IF
        
      END DO
    END DO
  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Differntial !
ELSEIF(rotation_rule == 2) THEN

  ! Assign angular velocity !
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2

        ! polar radius !
        s_eq = x2(j)*sin2(k)

        ! get omega !
        omega = w_0*s_core**2/(s_core**2 + x2(j)**2)

        ! Select based on conditions !
        IF(prim2(irho2,j,k,l) > prim2(irho2,nx_2,1,1)) THEN
          prim2(ivel2_z,j,k,l) = omega*s_eq
        ELSE
          prim2(ivel2_z,j,k,l) = 0.0d0                   
        END IF
        
      END DO
    END DO
  END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field !

! vector potential !
DO l = 0, nz_2
  DO k = 0, ny_2
    DO j = 0, nx_2
      a_phi(j,k,l) = xF2(j)*sin2f(k)*0.5*b_0 !*s_core**3/(s_core**3 + xF2(j)**3)*xF2(j)*sin2f(k)
    END DO
  END DO
END DO

! magnetic field !
prim2(ibx:ibz,:,:,:) = 0.0d0
DO l = 0, nz_2
  DO k = 0, ny_2
    DO j = 0, nx_2
      prim2(ibx,j,k,l) = (sin2f(k)*a_phi(j,k,l) - sin2f(k-1)*a_phi(j,k-1,l))*xF2(j)/(xF2(j)**2*dcos2(k) + small_num)
      prim2(iby,j,k,l) = -(xF2(j)*a_phi(j,k,l) - xF2(j-1)*a_phi(j-1,k,l))*sin2f(k)/(0.5d0*dx2_sq(j)*sin2f(k) + small_num)
    END DO
  END DO
END DO

! Find DIVB !
div_b = 0.0D0
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
      div_b = (xF2(j)**2*prim2(ibx,j,k,l) - xF2(j-1)**2*prim2(ibx,j-1,k,l))*dcos2(k)*dz2(l) &
            + (sin2f(k)*prim2(iby,j,k,l) - sin2f(k-1)*prim2(iby,j,k-1,l))*0.5d0*dx2_sq(j)*dz2(l) &
            + (prim2(ibz,j,k,l) - prim2(ibz,j,k,l-1))*0.5d0*dx2_sq(j)*dy2(k)
      maxdb = MAX(maxdb, ABS(div_b))
    END DO
  END DO
END DO
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Add perturbation to the density !
DO l = 1, nz_2
  DO k = 1, ny_2 
    DO j = 1, nx_2
      !IF(prim2(irho2,j,k,l) > prim2(irho2,nx_2,1,1)) THEN
      !  CALL RANDOM_NUMBER(rand)
      !  prim2(irho2,j,k,l) = prim2(irho2,j,k,l)*(1.0D0 + 1.0d-4*(rand - 0.5D0)/(0.5D0))
      !END IF
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate
Deallocate(a_phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = prim2(irho2,nx_2,1,1)
prim2_a(itau2) = prim2(itau2,nx_2,1,1)
eps2_a = prim2_a(itau2)/prim2_a(irho2)/(ggas2 - 1.0d0)

! Initialize !
m_inj = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 15.0d0*tcgs2code
output_profiletime = 0.03D0*tcgs2code 

END SUBROUTINE

