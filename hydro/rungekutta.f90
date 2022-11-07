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
USE DEFINITION
IMPLICIT NONE

! The input step number
INTEGER, INTENT (IN) :: n_in

! Dummy variables
INTEGER :: i, j, k, l
! Dummy !
REAL (DP) :: rhoaold, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

! Backup old arrays !
IF(DM_flag)THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		u_old1 (i,j,k,l) = cons1 (i,j,k,l)
	END DO
END IF
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	u_old2 (i,j,k,l) = cons2 (i,j,k,l)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

! Discretize !
CALL SPATIAL

IF (DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
   		cons1 (i,j,k,l) = u_old1 (i,j,k,l) + 0.391752226571890D0 * dt * l1 (i,j,k,l)
	END DO
   	CALL BOUNDARYU_DM
END IF

! NM sector !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	cons2 (i,j,k,l) = u_old2 (i,j,k,l) + 0.391752226571890D0 * dt * l2 (i,j,k,l)
END DO

! Copy the data to ghost cells in U
CALL BOUNDARYU_NM

! Convert from conservative to primitive
CALL FROMUTORVE

! Update 
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

CALL SPATIAL

IF (DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		cons1 (i,j,k,l) = 0.444370493651235D0 * u_old1 (i,j,k,l) + 0.555629506348765D0 * cons1 (i,j,k,l) + 0.368410593050371D0 * dt * l1 (i,j,k,l)
		u2_dm (i,j,k,l) = cons1 (i,j,k,l)
	END DO
 	CALL BOUNDARYU_DM
END IF

! NM sector !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	cons2 (i,j,k,l) = 0.444370493651235D0 * u_old2 (i,j,k,l) + 0.555629506348765D0 * cons2 (i,j,k,l) + 0.368410593050371D0 * dt * l2 (i,j,k,l)
	u2_nm (i,j,k,l) = cons2 (i,j,k,l)
END DO

! Copy the data to the ghost cells in U
CALL BOUNDARYU_NM

! Convert from conservative to primitive
CALL FROMUTORVE

! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3rd iteration

CALL SPATIAL

IF (DM_flag) THEN   
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		cons1 (i,j,k,l) = 0.620101851488403D0 * u_old1 (i,j,k,l) + 0.379898148511597D0 * cons1 (i,j,k,l) + 0.251891774271694D0 * dt * l1 (i,j,k,l)
		u3_dm (i,j,k,l) = cons1 (i,j,k,l)
		l3_dm (i,j,k,l) = l1 (i,j,k,l)	
	END DO
 	CALL BOUNDARYU_DM
END IF

! NM sector ! 
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	cons2 (i,j,k,l) = 0.620101851488403D0 * u_old2 (i,j,k,l) + 0.379898148511597D0 * cons2 (i,j,k,l) + 0.251891774271694D0 * dt * l2 (i,j,k,l)
	u3_nm (i,j,k,l) = cons2 (i,j,k,l)
	l3_nm (i,j,k,l) = l2 (i,j,k,l)
ENDDO

! Copy the data to the ghost cells in U
CALL BOUNDARYU_NM

! Convert from conservative to primitive
CALL FROMUTORVE

! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourth iteration

CALL SPATIAL

IF (DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		cons1 (i,j,k,l) = 0.178079954393132D0 * u_old1 (i,j,k,l)  + 0.821920045606868D0 * cons1 (i,j,k,l)  + 0.544974750228521D0 * dt * l1 (i,j,k,l) 
	END DO
 	CALL BOUNDARYU_DM
END IF

! NM sector !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	cons2 (i,j,k,l) = 0.178079954393132D0 * u_old2 (i,j,k,l) + 0.821920045606868D0 * cons2 (i,j,k,l) + 0.544974750228521D0 * dt * l2 (i,j,k,l)
END DO

! Copy the data to ghost cells in U
CALL BOUNDARYU_NM

! Convert from conservative to primitive
CALL FROMUTORVE

! Update physical quantities
CALL UPDATE (0)

! Convert again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

CALL SPATIAL

IF (DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		cons1 (i,j,k,l) = 0.517231671970585D0 * u2_dm (i,j,k,l) + 0.096059710526147D0 * u3_dm (i,j,k,l) + 0.386708617503269D0 * cons1 (i,j,k,l) &
						+ 0.063692468666290D0 * dt * l3_dm (i,j,k,l) + 0.226007483236906D0 * dt * l1 (i,j,k,l)
	END DO
 	CALL BOUNDARYU_DM
END IF

! NM sector !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	cons2 (i,j,k,l) = 0.517231671970585D0 * u2_nm (i,j,k,l) + 0.096059710526147D0 * u3_nm (i,j,k,l) + 0.386708617503269D0 * cons2 (i,j,k,l) &
					+ 0.063692468666290D0 * dt * l3_nm (i,j,k,l) + 0.226007483236906D0 * dt * l2 (i,j,k,l)
END DO

! Copy the data to ghost cells in U
CALL BOUNDARYU_NM

! Convert from conservative to primitive
CALL FROMUTORVE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for operator splitting

!CALL OPERATOR_SPLIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for adjusting atmospheric density !

! Do for NM !
IF(fixrhonm_flag) THEN

	! look for minimum atmospheri density !
	rhoaold = minval(prim2(:,:,:,irho2))
	!prim2_a(irho2) = min(min(maxval(prim2(:,:,:,irho2)), rho2_c)*rhofac_2, prim2_a(irho2))

	! Adjust density !
	DO CONCURRENT (j = 1:nx_2, k = 1:ny_2, l = 1:nz_2)
		IF(prim2(j,k,l,irho2) == rhoaold) THEN
			prim2(j,k,l,irho2) = prim2_a(irho2)
		END IF
	END DO

END IF

! For DM !
IF(fixrhodm_flag) THEN

	! look for minimum atmospheri density !
	rhoaold = minval(prim1(:,:,:,irho1))
	!prim1_a(irho1) = min(min(maxval(prim1(:,:,:,irho1)), rho1_c)*rhofac_1, prim1_a(irho1))

	! Adjust density !
	DO CONCURRENT (j = 1:nx_1, k = 1:ny_1, l = 1:nz_1)
		IF(prim1(j,k,l,irho1) == rhoaold) THEN
			prim1(j,k,l,irho1) = prim1_a(irho1)
		END IF
	END DO

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
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Local maximum effective speed
REAL (DP) :: lambda, lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL (DP) :: dt_temp1, dt_temp2

! Local minimum dt for DM, NM and 1st overlayer
REAL (DP) :: dt_out1, dt_out2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set !
dt_out1 = 1.0D5
dt_out2 = 1.0D5

! Now we find the minimum time constrained by DM sector
IF(DM_flag) THEN
   	DO j = nx_min_1, nx_part_1
		DO k = ny_min_1, ny_part_1
			DO l = nz_min_1, nz_part_1

	   			! Only grid with density above threshold density is counted
	   			IF(prim1(j,k,l,irho1) > prim1_a(irho1)) THEN
	      			lambda1 = ABS(prim1(j,k,l,ivel1_x)) + cs1(j,k,l)
         			lambda2 = ABS(prim1(j,k,l,ivel1_y)) + cs1(j,k,l)
					lambda3 = ABS(prim1(j,k,l,ivel1_z)) + cs1(j,k,l)
         			lambda = MAX(lambda1, lambda2, lambda3)
					dt_temp1 = dr1(j)
	     			IF(coordinate_flag == 2) THEN
						IF(n_dim > 1) THEN
							dt_temp1 = MIN(dt_temp1, x1(j)*sin1(j,k,l)*dy1)
						END IF
						IF(n_dim > 2) THEN
							dt_temp1 = MIN(dt_temp1, x1(j)*dz1)
						END IF
	      			ELSEIF(coordinate_flag == 1) THEN
						IF(n_dim > 2) THEN
	    	   				dt_temp1 = MIN(dt_temp1, x1(j)*dz1)
						END IF
					END IF
					dt_temp1 = dt_temp1*cfl/lambda
					dt_out1 = MIN(dt_out1, dt_temp1)
	    		ENDIF

			END DO
      	ENDDO
   	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we find the minimum time constrained by NM sector
DO j = nx_min_2, nx_part_2
	DO k = ny_min_2, ny_part_2
		DO l = nz_min_2, nz_part_2

	   		! Only grid with density above threshold density is counted
	   		IF(prim2(j,k,l,irho2) > prim2_a(irho2)) THEN
	      		lambda1 = ABS(prim2(j,k,l,ivel2_x)) + cs2(j,k,l)
         		lambda2 = ABS(prim2(j,k,l,ivel2_y)) + cs2(j,k,l)
				lambda3 = ABS(prim2(j,k,l,ivel2_z)) + cs2(j,k,l)
         		lambda = MAX(lambda1, lambda2, lambda3)
				dt_temp2 = dr2(j)
	     		IF(coordinate_flag == 2) THEN
					IF(n_dim > 1) THEN
						dt_temp2 = MIN(dt_temp2, x2(j)*sin2(j,k,l)*dy2)
					END IF
					IF(n_dim > 2) THEN
						dt_temp2 = MIN(dt_temp2, x2(j)*dz2)
					END IF
	      		ELSEIF(coordinate_flag == 1) THEN
					IF(n_dim > 2) THEN
	    	   			dt_temp2 = MIN(dt_temp2, x2(j)*dz2)
					END IF
				END IF
				dt_temp2 = dt_temp2*cfl/lambda
				dt_out2 = MIN(dt_out2, dt_temp2)
	    	ENDIF

		END DO
    ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Only the minimum one is chosen
dt = MIN(dt_out1, dt_out2)

END SUBROUTINE FindDt