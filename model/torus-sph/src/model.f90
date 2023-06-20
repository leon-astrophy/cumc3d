!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Integer !
INTEGER :: i, j, k, l

! Real variables !
REAL*8 :: rand
REAL*8 :: h_square
REAL*8 :: c_const
REAL*8 :: k_poly
REAL*8 :: brac
REAL*8 :: rho
REAL*8 :: rho_a
REAL*8 :: s_eq
REAL*8 :: omega
REAL*8 :: t_scale

! Normalization !
REAL*8 :: pmag
REAL*8 :: pgas
REAL*8 :: beta_old
REAL*8 :: beta_min
REAL*8 :: beta_loc
REAL*8 :: B_0
REAL*8 :: b_2
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

! Find h**2, for which w = hs**(-q) !
h_square = (2.0D0 - 2.0D0*q_grad)*(schwarzschild(s_0,r_sh) - schwarzschild(s_1,r_sh)) & 
          /(s_0**(2.0D0 - 2.0D0*q_grad) - s_1**(2.0D0 - 2.0D0*q_grad))
          
! Find integration constant !
c_const = (s_0**(2.0D0 - 2.0D0*q_grad)*schwarzschild(s_1,r_sh) - s_1**(2.0D0 - 2.0D0*q_grad)*schwarzschild(s_0,r_sh)) & 
        /(s_0**(2.0D0 - 2.0D0*q_grad) - s_1**(2.0D0 - 2.0D0*q_grad))

! Find polytropic constant, defined at the position of maximum density along equator, and we set the max density to be 1 !
k_poly = c_const + h_square*s_max**(2.0D0 - 2.0D0*q_grad)/(2.0D0 - 2.0D0*q_grad) - schwarzschild(s_max,r_sh)

! Exit condition !
IF(k_poly < 0.0d0) THEN
  STOP 'The polytropic constant is negative'
END IF

! Find polytropic constant !
k_poly = k_poly*(gamma - 1.0d0)/gamma/(rho_max**(gamma - 1.0d0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! solve for the density profile !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2

      ! polar radius !
      s_eq = x2(j)*SIN(y2(k))

      ! bracket term !
      brac = c_const + h_square*s_eq**(2.0D0 - 2.0D0*q_grad)/(2.0D0 - 2.0D0*q_grad) - schwarzschild(x2(j),r_sh)
      IF(brac < 0.0d0) THEN
        brac = 0.0d0
      END IF
      brac = brac*(gamma - 1.0d0)/gamma/k_poly

      ! torus density and pressure
      rho = brac**(1.0d0/(gamma - 1.0d0))
      rho_a = rho_max*rho_fac

      ! Select based on conditions !
      IF(rho >= rho_a .AND. s_eq >= s_0 .AND. s_eq <=s_1) THEN
        prim2(irho2,j,k,l) = rho
        omega = SQRT(h_square)*s_eq**(-q_grad)
        prim2(ivel2_z,j,k,l) = omega*s_eq
      ELSE
        prim2(irho2,j,k,l) = rho_a
        prim2(ivel2_z,j,k,l) = 0.0d0                   
      END IF
      
      ! vector potential !
      IF(prim2(irho2,j,k,l) >= rho_cut) THEN
        A_phi(j,k,l) = (prim2(irho2,j,k,l) - rho_cut)/rho_max
      END IF
      
    END DO
  END DO
END DO

! Find pressure and epsilon
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      prim2(itau2,j,k,l) = k_poly*prim2(irho2,j,k,l)**(gamma)
      epsilon2(j,k,l) = k_poly*prim2(irho2,j,k,l)**(gamma - 1.0D0)/(gamma - 1.0D0)
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      bx_face(j,k,l) = (SIN(yF2(k))*A_corner(j,k,l) - SIN(yF2(k-1))*A_corner(j,k-1,l))/(xF2(j))/dcos2(k)!dy2(k)/SIN(y2(k))
      by_face(j,k,l) = -((xF2(j))*A_corner(j,k,l) - (xF2(j-1))*A_corner(j-1,k,l))/dx2(j)/x2(j)
    END DO
  END DO
END DO

! Get cell-centered magnetic fields by averaging !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      prim2(ibx,j,k,l) = 1.5D0*dx2(j)/dx2_cb(j)*(xF2(j)**2*bx_face(j,k,l) + xF2(j-1)**2*bx_face(j-1,k,l))
      prim2(iby,j,k,l) = 0.5D0*dy2(k)/dcos2(k)*(SIN(yF2(k))*by_face(j,k,l) + SIN(yF2(k-1))*by_face(j,k-1,l))
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Normalize magnetic field so that minimum magnetisation is beta !
IF(normalize_by_minbeta) THEN
  beta_min = 1.0D30
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2
        b_2 = prim2(ibx,j,k,l)**2 + prim2(iby,j,k,l)**2 + prim2(ibz,j,k,l)**2
        IF(b_2 > 0.0d0) THEN
          beta_loc = 2.0d0*prim2(itau2,j,k,l)/b_2
          beta_min = MIN(beta_min, beta_loc)
        END IF
      END DO
    END DO
  END DO

  ! Normalization factor !
  B_0 = SQRT(beta_min/p_beta)

! Or by volume average plasma beta !
ELSE IF(normalize_by_vol) THEN
  pgas = 0.0d0
  pmag = 0.0d0
  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2
        IF(prim2(irho2,j,k,l) > rho_a) THEN
          pgas = pgas + prim2(itau2,j,k,l)*vol2(j,k,l)
          pmag = pmag + 0.5D0*(prim2(ibx,j,k,l)**2 + prim2(iby,j,k,l)**2 + prim2(ibz,j,k,l)**2)*vol2(j,k,l)
        END IF
      END DO
    END DO
  END DO

  ! Normalization factor !
  beta_old = pgas/pmag
  B_0 = SQRT(beta_old/p_beta)

END IF

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
      bx_face(j,k,l) = (SIN(yF2(k))*A_corner(j,k,l) - SIN(yF2(k-1))*A_corner(j,k-1,l))/(xF2(j))/dcos2(k)!dy2(k)/SIN(y2(k))
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
              (SIN(yF2(k))*by_face(j,k,l) - SIN(yF2(k-1))*by_face(j,k-1,l))/(x2(j)*dcos2(k))!*SIN(y2(k))*dy2(k))
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
      prim2(ibx,j,k,l) = 1.5D0*dx2(j)/dx2_cb(j)*(xF2(j)**2*bx_face(j,k,l) + xF2(j-1)**2*bx_face(j-1,k,l))
      prim2(iby,j,k,l) = 0.5D0*dy2(k)/dcos2(k)*(SIN(yF2(k))*by_face(j,k,l) + SIN(yF2(k-1))*by_face(j,k-1,l))
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = rho_max*rho_fac
prim2_a(itau2) = k_poly*prim2_a(irho2)**(gamma)
eps2_a = k_poly*prim2_a(irho2)**(gamma - 1.0D0)/(gamma - 1.0D0)

! orbital time scasle
t_scale = 2.0d0*pi/(SQRT(h_square)*s_max**(-q_grad))

! Set output profile interval !
total_time = 1000.0d0 !*t_scale
output_profiletime = 20.0d0 !total_time/20.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate ! 
DEALLOCATE(A_phi)
DEALLOCATE(A_corner)
DEALLOCATE(bx_face)
DEALLOCATE(by_face)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! choose your pseudo-newtonian potentials !

contains

  REAL*8 function newton(r_in)
  implicit none
  REAL*8 :: r_in
  newton = -1.0d0/r_in
  end function

  REAL*8 function schwarzschild(r_in,r_g)
  implicit none
  REAL*8 :: r_in, r_g
  schwarzschild = -1.0d0/(r_in - r_g)
  end function

END SUBROUTINE
