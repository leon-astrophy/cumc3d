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
! no magnetic field !
prim2(ibx:ibz,:,:,:) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate
Deallocate(a_phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = prim2(irho2,nx_2,1,1)
prim2_a(itau2) = prim2(itau2,nx_2,1,1)
eps2_a = prim2_a(itau2)/prim2_a(irho2)/(ggas2 - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 1000.0d0 
output_profiletime = total_time/50.0d0

END SUBROUTINE

