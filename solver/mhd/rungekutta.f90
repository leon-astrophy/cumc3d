!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine does one single Runge-Kutta full step
! It uses the opeator splitting and separate
! all non-gravitational source term to be done 
! after the hydro step.
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RUNGEKUTTA (n_in)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! The input step number
INTEGER, INTENT (IN) :: n_in

! Dummy variables
INTEGER :: i, j, k, l
! Dummy !
REAL*8 :: rhoaold, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

! Backup old arrays !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2 
				u_old2 (i,j,k,l) = cons2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

! Discretize !
CALL SPATIAL

! NM sector !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = u_old2 (i,j,k,l) + dt * l2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

! Convert from conservative to primitive
CALL FROMUTORVE

! Update 
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

! Discretize !
CALL SPATIAL

! NM sector !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = rk20 * u_old2(i,j,k,l) + rk21 * cons2 (i,j,k,l) + rk22 * dt * l2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

! Convert from conservative to primitive
CALL FROMUTORVE

! Update 
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

CALL SPATIAL

! NM sector !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = rk30 * u_old2(i,j,k,l) + rk31 * cons2 (i,j,k,l) + rk32 * dt * l2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

! Convert from conservative to primitive
CALL FROMUTORVE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for operator splitting

CALL OPERATOR_SPLIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for adjusting atmospheric density !

! Do for NM !
IF(fixrhonm_flag) THEN

	! look for minimum atmospheri density !
	rhoaold = minval(prim2(irho2,:,:,:))
	!prim2_a(irho2) = min(min(maxval(prim2(irho2,:,:,:)), rho2_c)*rhofac_2, prim2_a(irho2))

	! Adjust density !
	!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				IF(prim2(irho2,j,k,l) == rhoaold) THEN
					prim2(irho2,j,k,l) = prim2_a(irho2)
				END IF
			END DO
		END DO
	END DO
	!$OMP END PARALLEL DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check density !
IF (checkrho_flag) THEN
	CALL CHECKRHO
END IF

! Update physical quantities
IF (MOD (n_in, 2) == 0) THEN
	CALL UPDATE (1)
ELSE
	CALL UPDATE (0)
END IF

! Update again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! This subroutine calculates the maximum time step
! which satisfies the Courant condition 
! Written by Leung Shing Chi in 2016  
! If you modify the Euler equation, make sure you change this 
! part to include the new effective sound speed
! Limiters are posed based on output time and running time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE finddt
USE definition
USE MHD_module
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! For MHD speed !
REAL*8 :: a2_mhd, b2_mhd
REAL*8 :: bx_mhd, by_mhd, bz_mhd
REAL*8 :: cfx_mhd, cfy_mhd, cfz_mhd

! Local maximum effective speed
REAL*8 :: lambda, lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_temp1, dt_temp2

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_out1, dt_out2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set !
dt_out1 = 1.0D5
dt_out2 = 1.0D5

! Now we find the minimum time constrained by NM sector
!$OMP PARALLEL PRIVATE(a2_mhd, b2_mhd, bx_mhd, by_mhd, bz_mhd, cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, dt_temp1, dt_temp2)
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(MIN:dt_out1,dt_out2)
DO j = nx_min_2, nx_part_2
	DO k = ny_min_2, ny_part_2
		DO l = nz_min_2, nz_part_2

			! Only grid with density above threshold density is counted
			a2_mhd = cs2(j,k,l)**2
			bx_mhd = SQRT(prim2(ibx,j,k,l)**2/prim2(irho2,j,k,l))
			by_mhd = SQRT(prim2(iby,j,k,l)**2/prim2(irho2,j,k,l))
			bz_mhd = SQRT(prim2(ibz,j,k,l)**2/prim2(irho2,j,k,l))
			b2_mhd = bx_mhd**2 + by_mhd**2 + bz_mhd**2
			cfx_mhd = SQRT(0.5D0*(a2_mhd + b2_mhd + SQRT((a2_mhd + b2_mhd)**2 - 4.0D0*a2_mhd*bx_mhd**2)))
			cfy_mhd = SQRT(0.5D0*(a2_mhd + b2_mhd + SQRT((a2_mhd + b2_mhd)**2 - 4.0D0*a2_mhd*by_mhd**2)))
			cfz_mhd = SQRT(0.5D0*(a2_mhd + b2_mhd + SQRT((a2_mhd + b2_mhd)**2 - 4.0D0*a2_mhd*bz_mhd**2)))
			lambda1 = ABS(prim2(ivel2_x,j,k,l)) + cfx_mhd
			lambda2 = ABS(prim2(ivel2_y,j,k,l)) + cfy_mhd
			lambda3 = ABS(prim2(ivel2_z,j,k,l)) + cfz_mhd
			lambda = MAX(lambda1, lambda2, lambda3)

			! Look for minimum grid size !
			dt_temp2 = dx2(j)
			IF(coordinate_flag == 0) THEN
				IF(n_dim > 1) THEN
					dt_temp2 = MIN(dt_temp2, dy2(k))
				END IF
				IF(n_dim > 2) THEN
					dt_temp2 = MIN(dt_temp2, dz2(l))
				END IF
			ELSEIF(coordinate_flag == 1) THEN
				IF(n_dim > 1 .AND. ny_2 > 1) THEN
					dt_temp2 = MIN(dt_temp2, x2(j)*dy2(k))
				END IF
				IF(n_dim > 2) THEN
					dt_temp2 = MIN(dt_temp2, dz2(l))
				END IF
			ELSEIF(coordinate_flag == 2) THEN
				IF(n_dim > 1) THEN
					dt_temp2 = MIN(dt_temp2, x2(j)*dy2(k))
				END IF
				IF(n_dim > 2) THEN
					dt_temp2 = MIN(dt_temp2, x2(j)*SIN(y2(k))*dz2(l))
				END IF
			END IF
			dt_temp2 = dt_temp2*cfl/lambda
			dt_out2 = MIN(dt_out2, dt_temp2)
			
		END DO
	ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Only the minimum one is chosen
dt = MIN(dt_out1, dt_out2)

END SUBROUTINE FindDt
