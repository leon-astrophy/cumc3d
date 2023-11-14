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
Allocate(a_phi(-2:nx+3,-2:ny+3,-2:nz+3))

! Poisson interpolation coefficient !
call get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read parameter !
OPEN(UNIT=999, FILE = './profile/par.dat', ACTION='READ')
DO j = 1, 1
	READ(999,*) ggas, m_bh, masscgs2code, lengthcgs2code, tcgs2code
ENDDO
CLOSE(999)

! Velocity !
vel2code = (lengthcgs2code/tcgs2code)

! Magnetic field !
gauss2code = (masscgs2code/lengthcgs2code)**(0.5D0)/tcgs2code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversion

! They are all in params.h !
s_0 = s_0*lengthcgs2code
w_0 = w_0/tcgs2code
b_0 = b_0*gauss2code
v_crit = v_crit*vel2code
max_alv = max_alv*vel2code

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
IF(nlines .ne. nx+6) THEN
  WRITE (*,*) 'number of simulation grids from files', nlines
  WRITE (*,*) 'number of simulation grids in the program', nx+6
  STOP 'inconsistent number of simulation grids, exit'
END IF

! Read !
OPEN(UNIT=999, FILE = './profile/hydro.dat', ACTION='READ')
DO j = -2, nx + 3
	READ(999,*) prim(irho,j,1,1), prim(itau,j,1,1)
ENDDO
CLOSE(999)

! Assign density profile !
DO j = -2, nx + 3
  prim(irho,j,:,:) = prim(irho,j,1,1)
  prim(itau,j,:,:) = prim(itau,j,1,1)
END DO

! Assign internal energy #
epsilon(:,:,:) = prim(itau,:,:,:)/prim(irho,:,:,:)/(ggas - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose rotation rules !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core radius !
DO j = 1, nx
  IF(prim(irho,j-1,1,1) <= prim(irho,nx,1,1) .AND. &
     prim(irho,j,1,1) > prim(irho,nx,1,1)) THEN
    s_core = x(j)
    EXIT
  END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rigid rotations !
IF(rotation_rule == 1) THEN

  ! Find outer most radius !
  rmax = MAXVAL(x)

  ! Find angular velocity !
  omega = v_crit/rmax

  ! Assign angular velocity !
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx

        ! polar radius !
        s_eq = x(j)*sine(k)

        ! Select based on conditions !
        IF(prim(irho,j,k,l) > prim(irho,nx,1,1)) THEN
          prim(ivz,j,k,l) = omega*s_eq
        ELSE
          prim(ivz,j,k,l) = 0.0d0                   
        END IF
        
      END DO
    END DO
  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Differntial !
ELSEIF(rotation_rule == 2) THEN

  ! Assign angular velocity !
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx

        ! polar radius !
        s_eq = x(j)*sine(k)

        ! get omega !
        omega = w_0*s_core**2/(s_core**2 + s_eq**2)

        ! Select based on conditions !
        IF(prim(irho,j,k,l) > prim(irho,nx,1,1)) THEN
          prim(ivz,j,k,l) = omega*s_eq
        ELSE
          prim(ivz,j,k,l) = 0.0d0                   
        END IF

      END DO
    END DO
  END DO

END IF
s_core = 5000.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field !

! vector potential !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      a_phi(j,k,l) = xF(j)*sinf(k)*0.5d0 
      !a_phi(j,k,l) = s_core**3/(s_core**3 + xF(j)**3)*xF(j)*sinf(k)
    END DO
  END DO
END DO

! magnetic field !
prim(ibx:ibz,:,:,:) = 0.0d0
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,l) = (sinf(k)*a_phi(j,k,l) - sinf(k-1)*a_phi(j,k-1,l))/(xF(j)*dcose(k)+small_num)
      prim(iby,j,k,l) = - (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
    END DO
  END DO
END DO

! Find DIVB !
div_b = 0.0D0
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      div_b = (xF(j)*xF(j)*prim(ibx,j,k,l) - xF(j-1)*xF(j-1)*prim(ibx,j-1,k,l))/(dx_cb(j)/3.0d0) &

            + (sinf(k)*prim(iby,j,k,l) - sinf(k-1)*prim(iby,j,k-1,l))*(x(j)*dx(j))/(dx_cb(j)*dcose(k)/3.0d0) &
      
            + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))*(x(j)*dx(j)*dy(k))/(dx_cb(j)*dcose(k)*dz(l)/3.0d0)
            
      maxdb = MAX(maxdb, div_b)
    END DO
  END DO
END DO
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

prim(ibx:ibz,:,:,:) = prim(ibx:ibz,:,:,:)*b_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Add perturbation to the density !
DO l = 1, nz
  DO k = 1, ny 
    DO j = 1, nx
      IF(prim(irho,j,k,l) > prim(irho,nx,1,1)) THEN
        CALL RANDOM_NUMBER(rand)
        prim(irho,j,k,l) = prim(irho,j,k,l)*(1.0D0 + 1.0d-4*(rand - 0.5D0)/(0.5D0))
      END IF
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate
Deallocate(a_phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim_a(:) = 0.0D0
prim_a(irho) = prim(irho,nx,1,1)
prim_a(itau) = prim(itau,nx,1,1)
eps_a = prim_a(itau)/prim_a(irho)/(ggas - 1.0d0)

! Initialize !
m_inj = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 15.0d0*tcgs2code
output_profiletime = 0.03D0*tcgs2code 

END SUBROUTINE

