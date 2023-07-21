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

SUBROUTINE RUNGEKUTTA
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l
! Dummy !
REAL*8 :: rhoaold, dummy

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end, time1, time0
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time0)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! Backup old arrays !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
			DO i = imin2, imax2 
				u_old2 (i,j,k,l) = cons2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'backup = ', REAL(time_end - time_start) / rate
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

! Discretize !
CALL SPATIAL

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = u_old2 (i,j,k,l) + dt * l2 (i,j,k,l) 
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk1 = ', REAL(time_end - time_start) / rate
#endif

! Convert from conservative to primitive
CALL FROMUTORVE

! Check density !
IF (checkrho_flag) THEN
	CALL CHECKRHO
END IF

! Update 
CALL UPDATE (1)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

! Discretize !
CALL SPATIAL

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = rk20 * u_old2(i,j,k,l) + rk21 * cons2 (i,j,k,l) + rk22 * dt * l2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk2 = ', REAL(time_end - time_start) / rate
#endif

! Convert from conservative to primitive
CALL FROMUTORVE

! Check density !
IF (checkrho_flag) THEN
	CALL CHECKRHO
END IF

! Update 
CALL UPDATE (2)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

CALL SPATIAL

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
			DO i = imin2, imax2 
				cons2 (i,j,k,l) = rk30 * u_old2(i,j,k,l) + rk31 * cons2 (i,j,k,l) + rk32 * dt * l2 (i,j,k,l)
			END DO
		END DO
	END DO
END DO 
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk3 = ', REAL(time_end - time_start) / rate
#endif

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

	! Adjust density !
	!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				IF(prim2(irho2,j,k,l) == rhoaold) THEN
					prim2(irho2,j,k,l) = prim2_a(irho2)
				END IF
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END PARALLEL DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check density !
IF (checkrho_flag) THEN
	CALL CHECKRHO
END IF

! Update physical quantities
CALL UPDATE (3)

! Update again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time1)
WRITE(*,*) 'rk total = ', REAL(time1 - time0) / rate
#endif

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
REAL*8 :: a4_mhd, b4_mhd
REAL*8 :: b2x_mhd, b2y_mhd, b2z_mhd
REAL*8 :: cfx_mhd, cfy_mhd, cfz_mhd

! Local maximum effective speed
REAL*8 :: lambda, lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_temp1, dt_temp2

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_out1, dt_out2

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set !
dt_out1 = 1.0D5
dt_out2 = 1.0D5

! Now we find the minimum time constrained by NM sector
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) & 
!$OMP PRIVATE(a2_mhd, b2_mhd, a4_mhd, b4_mhd, b2x_mhd, b2y_mhd, b2z_mhd, &
!$OMP cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, dt_temp2) REDUCTION(MIN:dt_out2)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
!$ACC PRIVATE(a2_mhd, b2_mhd, a4_mhd, b4_mhd, b2x_mhd, b2y_mhd, b2z_mhd, & 
!$ACC cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, dt_temp2) REDUCTION(MIN:dt_out2)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2

			! Only grid with density above threshold density is counted
			a2_mhd = cs2(j,k,l)*cs2(j,k,l)
			a4_mhd = a2_mhd*a2_mhd
			b2x_mhd = (bcell(ibx,j,k,l)*bcell(ibx,j,k,l)/prim2(irho2,j,k,l))
			b2y_mhd = (bcell(iby,j,k,l)*bcell(iby,j,k,l)/prim2(irho2,j,k,l))
			b2z_mhd = (bcell(ibz,j,k,l)*bcell(ibz,j,k,l)/prim2(irho2,j,k,l))
			b2_mhd = b2x_mhd + b2y_mhd + b2z_mhd
			b4_mhd = b2_mhd*b2_mhd
			cfx_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2x_mhd)))
			cfy_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2y_mhd)))
			cfz_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2z_mhd)))
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
					dt_temp2 = MIN(dt_temp2, x2(j)*sin2(k)*dz2(l))
				END IF
			END IF
			dt_temp2 = dt_temp2*cfl/lambda
			dt_out2 = MIN(dt_out2, dt_temp2)
			
		END DO
	ENDDO
ENDDO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Only the minimum one is chosen
dt = MIN(dt_out1, dt_out2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'finddt = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE FindDt
