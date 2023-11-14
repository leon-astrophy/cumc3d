!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Real variables !
REAL*8 :: rand_num
REAL*8 :: c_const
REAL*8 :: k_poly
REAL*8 :: enthalpy
REAL*8 :: rho_a
REAL*8 :: rho_lim
REAL*8 :: rho_local
REAL*8 :: omega_vel
REAL*8 :: t_scale
REAL*8 :: lkep2
REAL*8 :: omega

! Normalization !
REAL*8 :: pgas
REAL*8 :: pmag
REAL*8 :: beta_min
REAL*8 :: beta_loc
REAL*8 :: beta_old
REAL*8 :: B_0
REAL*8 :: b_2
REAL*8 :: maxdb
REAL*8 :: div_b

! Vector potential !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: A_phi, A_corner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate ! 
ALLOCATE(A_phi(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(A_corner(-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! override adiabatic index !
ggas = gamma

! Square of kepler angular momentum !
lkep2 = s_max**(2.0d0*q_grad - 1.0d0)/(s_max - r_sh)**2

! Find integration constant !
c_const = schwarzschild(s_in,0.0d0,r_sh) - lkep2*s_in**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad)

! Find enthalpy !
enthalpy = - schwarzschild(s_max,0.0d0,r_sh) + lkep2*s_max**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad) + c_const

! Find polytropic constant, defined at the position of maximum density along equator, and we set the max density to be 1 !
k_poly = enthalpy*(ggas - 1.0d0)/ggas/rho_max**(ggas - 1.0d0)

! Exit condition !
IF(k_poly < 0.0d0) THEN
  STOP 'The polytropic constant is negative'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! solve for the density profile !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Find enthalpy !
      enthalpy = - schwarzschild(x(j),z(l),r_sh) + lkep2*x(j)**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad) + c_const

      ! Find density !
      IF(enthalpy < 0.0d0) THEN
        enthalpy = 0.0d0
      END IF

      ! torus density and pressure
      rho_local = (enthalpy*(ggas - 1.0d0)/ggas/k_poly)**(1.0d0/(gamma - 1.0d0))

      ! Choose !
      IF(corona) THEN
        rho_a = a_eta*rho_max*EXP(-a_corn*(schwarzschild(x(j),z(l),r_sh) - schwarzschild(s_in,0.0d0,r_sh)))
      ELSE
        rho_a = rho_fac*rho_max
      ENDIF
      rho_lim = rho_fac*rho_max

      ! Select based on conditions !
      IF (rho_local > rho_lim .AND. x(j) >= s_in) THEN
        prim(irho,j,k,l) = rho_local
        omega = DSQRT(lkep2)*x(j)**(-q_grad)
        prim(ivy,j,k,l) = omega*x(j)
        prim(itau,j,k,l) = k_poly*prim(irho,j,k,l)**(gamma)
        epsilon(j,k,l) = k_poly*prim(irho,j,k,l)**(gamma - 1.0D0)/(gamma - 1.0D0)
      ELSE
        prim(irho,j,k,l) = rho_a
        prim(ivy,j,k,l) = 0.0d0   
        IF(corona) THEN
          prim(itau,j,k,l) = prim(irho,j,k,l)/a_corn
          epsilon(j,k,l) = prim(itau,j,k,l)/prim(irho,j,k,l)/(ggas - 1.0d0)
        ELSE
          prim(itau,j,k,l) = k_poly*prim(irho,j,k,l)**(gamma)
          epsilon(j,k,l) = k_poly*prim(irho,j,k,l)**(gamma - 1.0D0)/(gamma - 1.0D0)
        END IF            
      END IF

      ! vector potential !
      IF(prim(irho,j,k,l) >= rho_cut) THEN
        A_phi(j,k,l) = (prim(irho,j,k,l) - rho_cut)/rho_max
      END IF

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get cell-corner vector potential !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      A_corner(j,k,l) = 0.25D0*(A_phi(j,k,l) + A_phi(j+1,k,l) + A_phi(j,k,l+1) + A_phi(j+1,k,l+1))
    END DO
  END DO
END DO

! First, initialize magnetic fields !
prim(ibx:ibz,:,:,:) = 0.0d0
bcell(ibx:ibz,:,:,:) = 0.0d0

! Get face-magnetic fields by cross product !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,l) = - (A_corner(j,k,l) - A_corner(j,k,l-1))/(dz(l))
      prim(ibz,j,k,l) = (xF(j)*A_corner(j,k,l) - xF(j-1)*A_corner(j-1,k,l))/(x(j)*dx(j))
    END DO
  END DO
END DO

! Get cell-centered magnetic fields by averaging !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = (xF(j)*prim(ibx,j,k,l) + xF(j-1)*prim(ibx,j-1,k,l))/(xF(j) + xF(j-1))
      bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Normalize magnetic field so that minimum magnetisation is beta !
IF(normalize_by_minbeta) THEN
  beta_min = 1.0D30
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        b_2 = bcell(ibx,j,k,l)**2 + bcell(iby,j,k,l)**2 + bcell(ibz,j,k,l)**2
        IF(b_2 > 0.0d0) THEN
          beta_loc = 2.0d0*prim(itau,j,k,l)/b_2
          beta_min = MIN(beta_min, beta_loc)
        END IF
      END DO
    END DO
  END DO

  ! Normalization factor !
  B_0 = DSQRT(beta_min/p_beta)

! Or by volume average plasma beta !
ELSE IF(normalize_by_vol) THEN
  pgas = 0.0d0
  pmag = 0.0d0
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        IF(prim(irho,j,k,l) > rho_a) THEN
          pgas = pgas + prim(itau,j,k,l)*vol(j,k,l)
          pmag = pmag + 0.5D0*(bcell(ibx,j,k,l)**2 + bcell(iby,j,k,l)**2 + bcell(ibz,j,k,l)**2)*vol(j,k,l)
        END IF
      END DO
    END DO
  END DO

  ! Normalization factor !
  beta_old = pgas/pmag
  B_0 = DSQRT(beta_old/p_beta)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recompute magnetic field to avoid numerical errors !

A_phi(:,:,:) = A_phi(:,:,:)*B_0
! Get cell-corner vector potential !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      A_corner(j,k,l) = 0.25D0*(A_phi(j,k,l) + A_phi(j+1,k,l) + A_phi(j,k,l+1) + A_phi(j+1,k,l+1))
    END DO
  END DO
END DO

! First, initialize magnetic fields !
prim(ibx:ibz,:,:,:) = 0.0d0
bcell(ibx:ibz,:,:,:) = 0.0d0

! Get face-magnetic fields by cross product !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,l) = -(A_corner(j,k,l) - A_corner(j,k,l-1))/(dz(l)) 
      prim(ibz,j,k,l) = (xF(j)*A_corner(j,k,l) - xF(j-1)*A_corner(j-1,k,l))/(x(j)*dx(j))
    END DO
  END DO
END DO

! Check divergence-B = 0 constraint !
maxdb = 0.0d0       
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      div_b = (xF(j)*prim(ibx,j,k,l) - xF(j-1)*prim(ibx,j-1,k,l))*dy(k)*dz(l) &
            + (prim(iby,j,k,l) - prim(iby,j,k-1,l))*dx(j)*dz(l) &
            + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))*(x(j)*dx(j))*dy(k)
      maxdb = MAX(maxdb, div_b)
    END DO
  END DO
END DO
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! orbital time scasle
t_scale = 2.0d0*pi/(DSQRT(lkep2)*s_max**(-q_grad))

! Set output profile interval !
total_time = 20.0D0*t_scale
output_profiletime = 10.0d0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! coronal floor !
IF(corona) THEN
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        rho_floor(j,k,l) = a_eta*rho_max*EXP(-a_corn*(schwarzschild(x(j),z(l),r_sh) - schwarzschild(s_in,0.0d0,r_sh)))
        p_floor(j,k,l) = rho_floor(j,k,l)/a_corn
        eps_floor(j,k,l) = p_floor(j,k,l)/rho_floor(j,k,l)/(ggas - 1.0d0)
      END DO
    END DO
  END DO
ELSE
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        rho_floor(j,k,l) = rho_max*rho_fac
        p_floor(j,k,l) = k_poly*rho_floor(j,k,l)**(ggas)
        eps_floor(j,k,l) = k_poly*rho_floor(j,k,l)**(ggas - 1.0D0)/(ggas - 1.0D0)
      END DO 
    END DO
  END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! perturb the torus !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      IF(prim(irho,j,k,l) > rho_floor(j,k,l)) THEN
        CALL RANDOM_NUMBER(rand_num)
        prim(irho,j,k,l) = prim(irho,j,k,l) + prim(irho,j,k,l)*(rand_num - 0.5d0)/(0.5d0)*1.0d-4
      END IF
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate !
DEALLOCATE(A_phi)
DEALLOCATE(A_corner)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! choose your pseudo-newtonian potentials !

contains

  REAL*8 function newton(r_in,z_in)
  implicit none
  REAL*8 :: r_in,z_in
  newton = -1.0d0/(DSQRT(r_in**2 + z_in**2))
  end function

  REAL*8 function schwarzschild(r_in,z_in,r_g)
  implicit none
  REAL*8 :: r_in, z_in, r_g
  schwarzschild = -1.0d0/(DSQRT(r_in**2 + z_in**2) - r_g)
  end function

END SUBROUTINE
