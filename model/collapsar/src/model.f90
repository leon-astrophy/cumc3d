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

! Angular velocity !
REAL*8 :: omega

! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !

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
        s_eq = x2(j)*DSIN(y2(k))

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
        s_eq = x2(j)*DSIN(y2(k))

        ! get omega !
        omega = w_0*s_0**2/(s_0**2 + x2(j)**2)!s_eq**2)

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

! Simlpe magnetic fied !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      
      ! polar radius ! 
      s_eq = x2(j)*DSIN(y2(k))

      ! Poloridal Magnetic field !
      prim2(ibz,j,k,l) = b_0*s_0**2/(s_0**2 + x2(j)**2)!s_eq**2)        

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = prim2(irho2,nx_2,1,1)
prim2_a(itau2) = prim2(itau2,nx_2,1,1)
eps2_a = prim2_a(itau2)/prim2_a(irho2)/(ggas2 - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 15.0d0*tcgs2code
output_profiletime = 0.01D0*tcgs2code !total_time/50.0d0

END SUBROUTINE

