!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Real !
REAL*8 :: s_eq
REAL*8 :: omega
REAL*8 :: B_0
REAL*8 :: maxdb
REAL*8 :: div_b

! Vector potential !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: A_phi, A_corner

! Magnetic fields (surface) !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: bx_face, by_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate ! 
ALLOCATE(A_phi(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(A_corner(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(bx_face(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(by_face(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! override adiabatic index !
ggas2 = gamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! solve for the density profile !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2

      ! equatorial radius !
      s_eq = x2(j)*SIN(y2(k))

      ! Assign density profile !
      prim2(irho2,j,k,l) = rho_max*(x2(j)/r_0)**(-rho_grad)

      ! Angular velocity
      omega = omega_max*(r_0**q_grad/(x2(j)**q_grad + r_0**q_grad))
      prim2(ivel2_z,j,k,l) = omega*s_eq
      
      ! Set density floor !
      IF(prim2(irho2,j,k,l) < rho_fac) THEN
        prim2(irho2,j,k,l) = rho_fac
        prim2(ivel2_z,j,k,l) = 0.0d0
      END IF

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find pressure and epsilon
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      prim2(itau2,j,k,l) = p_max*(prim2(irho2,j,k,l)/rho_max)**(ggas2)
      epsilon2(j,k,l) = prim2(itau2,j,k,l)/prim2(irho2,j,k,l)/(ggas2 - 1.0d0)
    END DO
  END DO
END DO

! Mangetic field !
B_0 = SQRT(2.0D0*p_max/p_beta)
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      s_eq = x2(j)*SIN(y2(k))
      A_phi(j,k,l) = 0.5D0*s_eq*B_0
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recompute magnetic field to avoid numerical errors !

A_phi(:,:,:) = A_phi(:,:,:)*B_0
! Get cell-corner vector potential !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      A_corner(j,k,l) = 0.25D0*(A_phi(j,k,l) + A_phi(j+1,k,l) + A_phi(j,k+1,l) + A_phi(j+1,k+1,l))
    END DO
  END DO
END DO

! First, initialize magnetic fields !
prim2(ibx:ibz,:,:,:) = 0.0d0

! Get face-magnetic fields by cross product !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      bx_face(j,k,l) = (SIN(yF2(k))*A_corner(j,k,l) - SIN(yF2(k-1))*A_corner(j,k-1,l))/dy2(k)/(xF2(j))/SIN(y2(k))
      by_face(j,k,l) = -((xF2(j))*A_corner(j,k,l) - (xF2(j-1))*A_corner(j-1,k,l))/dx2(j)/x2(j)
    END DO
  END DO
END DO

! Check divergence-B = 0 constraint !
div_b = 0.0d0
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      div_b = (xF2(j)**2*bx_face(j,k,l) - xF2(j-1)**2*bx_face(j-1,k,l))/(x2(j)**2*dx2(j)) + & 
              (SIN(yF2(k))*by_face(j,k,l) - SIN(yF2(k-1))*by_face(j,k-1,l))/(x2(j)*SIN(y2(k))*dy2(k))
      maxdb = MAX(maxdb, div_b)
    END DO
  END DO
END DO
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

! Get cell-centered magnetic fields by averaging !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      prim2(ibx,j,k,l) = 0.5D0*(bx_face(j,k,l) + bx_face(j-1,k,l))
      prim2(iby,j,k,l) = 0.5D0*(by_face(j,k,l) + by_face(j,k-1,l))
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = rho_fac
prim2_a(itau2) = p_max*(rho_fac/rho_max)**(ggas2)
eps2_a = prim2_a(itau2)/prim2_a(irho2)/(ggas2 - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate !
DEALLOCATE(A_phi)
DEALLOCATE(A_corner)
DEALLOCATE(bx_face)
DEALLOCATE(by_face)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 1000.0d0
output_profiletime = total_time/100.0d0

END SUBROUTINE
