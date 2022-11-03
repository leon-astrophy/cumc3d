!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine does one single Runge-Kutta full step
! It uses the opeator splitting and separate
! all non-gravitational source term to be done 
! after the hydro step.
!
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RUNGEKUTTA (n_in)
USE DEFINITION
USE PPT_MODULE
use ecap_module
USE TURB_MODULE
use nuceos_module
USE HELMEOS_MODULE
USE LEVELSET_MODULE
IMPLICIT NONE

! The input step number
INTEGER, INTENT (IN) :: n_in

! Dummy variables
INTEGER :: i, j, k

! Dummy dx due to Runge Kutta scheme
REAL (DP) :: dx1_old, dx1_two, dx1_three, dxdt1_three
REAL (DP) :: dx2_old, dx2_two, dx2_three, dxdt2_three

! Dummy !
REAL (DP) :: rhoaold, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Necessary variables for NUCEOS
! dummies for EOS call
real*8 eosdummy, ent_tmp
integer keyerr,keytemp
! end dummies for EOS call

! Setting the NUCEOS reading mode
keytemp = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

! Update the first half of particle tracer
IF(tracer_flag == 1) THEN
	call evolve_p_1st
END IF

! Find the boundaries !
IF(movinggridnm_flag == 1) THEN
	CALL FINDRADIUS_NM
END IF
IF(movinggriddm_flag == 1) THEN
	CALL FINDRADIUS_DM
END IF

! Backup !
IF(RUNDM_flag == 1)THEN
	u_old1 = u1
END IF
u_old2 = u2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

! Movinggrid !
IF(movinggridnm_flag == 1) THEN
	dx2_old = dx2
END IF
IF(movinggriddm_flag == 1) THEN
	dx1_old = dx1
END IF

! Discretize !
CALL SPATIAL

IF (RUNDM_flag == 1) THEN
   u1 (:,:,:) = u_old1 (:,:,:) + 0.391752226571890D0 * dt * l1 (:,:,:)
 	CALL BOUNDARY2D_DM
END IF

! NM sector !
u2 (:,:,:) = u_old2 (:,:,:) + 0.391752226571890D0 * dt * l2 (:,:,:)

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = dx1_old + 0.391752226571890D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = dx2_old + 0.391752226571890D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
   CALL FIND_AZBAR()
END IF

! Update 
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

IF(tracer_flag == 1) THEN
	call evolve_p_2nd
END IF

CALL SPATIAL

IF (RUNDM_flag == 1) THEN
	u1 (:,:,:) = 0.444370493651235D0 * u_old1 (:,:,:) + 0.555629506348765D0 * u1 (:,:,:) + 0.368410593050371D0 * dt * l1 (:,:,:)
	u2_dm (:,:,:) = u1 (:,:,:)
 	CALL BOUNDARY2D_DM
END IF

! NM sector !
u2 (:,:,:) = 0.444370493651235D0 * u_old2 (:,:,:) + 0.555629506348765D0 * u2 (:,:,:) + 0.368410593050371D0 * dt * l2 (:,:,:)
u2_nm (:,:,:) = u2 (:,:,:)

! Copy the data to the ghost cells in U
CALL BOUNDARY2D_NM

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.444370493651235D0 * dx1_old + 0.555629506348765D0 * dx1 + 0.368410593050371D0 * dt * vel1_max * dx1 / radius1
   dx1_two = dx1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.444370493651235D0 * dx2_old + 0.555629506348765D0 * dx2 + 0.368410593050371D0 * dt * vel2_max * dx2 / radius2
   dx2_two = dx2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3rd iteration

IF(tracer_flag == 1) THEN
	call evolve_p_3rd
END IF

CALL SPATIAL

IF (RUNDM_flag == 1) THEN   
	u1 (:,:,:) = 0.620101851488403D0 * u_old1 (:,:,:) + 0.379898148511597D0 * u1 (:,:,:) + 0.251891774271694D0 * dt * l1 (:,:,:)
	u3_dm (:,:,:) = u1 (:,:,:)
	l3_dm (:,:,:) = l1 (:,:,:)	
 	CALL BOUNDARY2D_DM
END IF

! NM sector ! 
u2 (:,:,:) = 0.620101851488403D0 * u_old2 (:,:,:) + 0.379898148511597D0 * u2 (:,:,:) + 0.251891774271694D0 * dt * l2 (:,:,:)
u3_nm (:,:,:) = u2 (:,:,:)
l3_nm (:,:,:) = l2 (:,:,:)

! Copy the data to the ghost cells in U
CALL BOUNDARY2D_NM

! Update dx
IF(movinggriddm_flag == 1) then
   dxdt1_three = vel1_max * dx1 / radius1
   dx1 = 0.620101851488403D0 * dx1_old + 0.379898148511597D0 * dx1 + 0.251891774271694D0 * dt * vel1_max * dx1 / radius1
   dx1_three = dx1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dxdt2_three = vel2_max * dx2 / radius2
   dx2 = 0.620101851488403D0 * dx2_old + 0.379898148511597D0 * dx2 + 0.251891774271694D0 * dt * vel2_max * dx2 / radius2
   dx2_three = dx2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourth iteration

IF(tracer_flag == 1) THEN
	call evolve_p_4th
END IF

CALL SPATIAL

IF (RUNDM_flag == 1) THEN
	u1 (:,:,:) = 0.178079954393132D0 * u_old1 (:,:,:)  + 0.821920045606868D0 * u1 (:,:,:)  + 0.544974750228521D0 * dt * l1 (:,:,:) 
 	CALL BOUNDARY2D_DM
END IF

! NM sector !
u2 (:,:,:) = 0.178079954393132D0 * u_old2 (:,:,:) + 0.821920045606868D0 * u2 (:,:,:) + 0.544974750228521D0 * dt * l2 (:,:,:)

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.178079954393132D0 * dx1_old + 0.821920045606868D0 * dx1 + 0.544974750228521D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.178079954393132D0 * dx2_old + 0.821920045606868D0 * dx2 + 0.544974750228521D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

! Update physical quantities
CALL UPDATE (0)

! Convert again !
CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

IF(tracer_flag == 1) THEN
	call evolve_p_5th
END IF

CALL SPATIAL 

IF (RUNDM_flag == 1) THEN
	u1 (:,:,:) = 0.517231671970585D0 * u2_dm (:,:,:) + 0.096059710526147D0 * u3_dm (:,:,:) + 0.386708617503269D0 * u1 (:,:,:) &
		+ 0.063692468666290D0 * dt * l3_dm (:,:,:) + 0.226007483236906D0 * dt * l1 (:,:,:)
 	CALL BOUNDARY2D_DM
END IF

! NM sector !
u2 (:,:,:) = 0.517231671970585D0 * u2_nm (:,:,:) + 0.096059710526147D0 * u3_nm (:,:,:) + 0.386708617503269D0 * u2 (:,:,:) &
		+ 0.063692468666290D0 * dt * l3_nm (:,:,:) + 0.226007483236906D0 * dt * l2 (:,:,:)

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.517231671970585D0 * dx1_two + 0.096059710526147D0 * dx1_three + 0.386708617503269D0 * dx1 &
		+ 0.063692468666290D0 * dt * dxdt1_three + 0.226007483236906D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.517231671970585D0 * dx2_two + 0.096059710526147D0 * dx2_three + 0.386708617503269D0 * dx2 &
		+ 0.063692468666290D0 * dt * dxdt2_three + 0.226007483236906D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE 

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN 
	CALL FIND_AZBAR()
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update temperature
IF (helmeos_flag == 1) THEN
	CALL FindhelmTEMP
END IF

! NUCEOS !
if(nuceos_flag == 1) THEN
	call findnuctemp
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for burning

! For Type Ia supernovae 
IF(fusion_flag == 1) THEN

   ! If there is level-set, update it
   IF(flame_flag == 1) CALL UPDATE_FLAME_RADIUS 

   ! This trigger the burning package proposed by
   ! Reinecke 1999b
   IF(xisotran_flag == 1) THEN

      ! This does the Carbon burning
      IF(carburn_flag == 1) CALL BURN_PHASE1B

      ! This do the O- and Si- burning
      IF(advburn_flag == 1) CALL BURN_PHASE2B

      ! Update the AZbar and temperature accordingly
      CALL FIND_AZBAR
      CALL FINDHELMTEMP

      ! For completely burnt zone, check if NSE applies
      IF(convert_nse_flag == 1) CALL NSE2

      ! Copy the new Xiso and epsilon to ghost cells
      CALL BOUNDARY2D_X()
      CALL BOUNDARY1D_NM(epsilon2, even)

      ! Check if the change of isotope perserve the sum
      CALL CHECKXISOTOPE

      ! Update the burntimw
      last_burntime = global_time

      ! Update Abar and Zbar and temperature again
      CALL Find_AZBAR()
      CALL FindhelmTemp

   ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For AIC, the electron capture occurs as 
! a function of density
if(ecap_flag == 1) then 
   if(bounce_flag == 0) then
      call findEcap
   endif
endif

! Check if the bounce occured or not
if(bounce_flag == 0 .and. ecap_flag == 1) then
   ! entropy condition for bounce
   keytemp = 1
   k = length_step_z_part_2
   do j=1,length_step_r_part_2,1
      if(rho2(j,1)>1.62d-8) then
         call nuc_eos_short(rho2(j,k)*6.171e17, temp2(j,k),ye2(j,k),eosdummy,eosdummy,ent_tmp, &
         eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,keytemp,keyerr,eos_rf_prec)
	      if(ent_tmp>=3) then
	         write(*,*) 'Bounce at', r2(j)
            bounce_flag = 1
            output_file = .true.
            total_time = global_time + 2.0d4 !total_time
	         exit
	      endif
      endif
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Section for adjusting atmospheric density !
IF(fixrhonm_flag == 1) THEN
	rhoaold = minval(rho2)
	rho2_a = min(min(maxval(rho2), rho2_c)*rhofac_2, rho2_a)

	! Assign atmoshperic epsilon !
	IF(nm_epsilon == 1 .AND. polyeosnm_flag == 1) THEN
		epsilon2_a = k_2 * rho2_a ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
	ELSEIF(helmeos_flag == 1) THEN
		CALL HELMEOS_RtoE (rho2_a, temp_a, abar_ini, zbar_ini, ye2_ini, epsilon2_a, dummy)
	ELSE
		CALL EOSEPSILON (rho2_a, dummy, epsilon2_a, 2)
	END IF

	! Adjust density !
	DO j = -4, length_step_r_2 + 5
		DO k = -4, length_step_z_2 + 5
			IF(rho2(j,k) == rhoaold) THEN
				rho2(j,k) = rho2_a
			END IF
		END DO
	END DO

	! Pressure !
	IF(helmeos_flag == 1) THEN
		CALL HELMEOS_RtoP(rho2_a, temp_a, abar_ini, zbar_ini, ye2_ini, dummy, dummy, dummy)
	END IF
END IF

! For DM !
IF(fixrhodm_flag == 1) THEN
	rhoaold = minval(rho1)
	rho1_a = min(min(maxval(rho1), rho1_c)*rhofac_1, rho1_a)

	! Assign atmoshperic epsilon !
	CALL EOSEPSILON (rho1_a, dummy, epsilon1_a, 1)

	! Adjust density !
	DO j = -4, length_step_r_1 + 5
		DO k = -4, length_step_z_1 + 5
			IF(rho1(j,k) == rhoaold) THEN
				rho1(j,k) = rho1_a
			END IF
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check density !
IF (checkrho_flag == 1) THEN
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

! Prepare data for the output purpose
IF (tracer_flag == 1) THEN
	call evolve_p_final
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Format !
100 format (20ES15.7)

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
INTEGER :: i, j, k

! Dummy
REAL (DP) :: dummy

! Local effective speed
REAL (DP) :: lambda2a, lambda2b

! Local maximum effective speed
REAL (DP) :: lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL (DP) :: dt_temp1, dt_temp2

! Initialize by setting an arbitrarily large number
dt_temp1 = 1.0D5
dt_temp2 = 1.0D5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we find the minimum time constrained by DM sector
! For non-MHD case
! For dark matter
IF(runDM_flag == 1) THEN
   DO k = length_step_z_min_part_1, length_step_z_part_1, 1
      DO j = 1, length_step_r_part_1, 1

	   ! Only grid with density above threshold density is counted
	   IF(rho1(j,k) > rho1_a) THEN
	      lambda2a = ABS(vel1_r(j,k)) + cs1(j,k)
         lambda2b = ABS(vel1_z(j,k)) + cs1(j,k)
         lambda1 = MAX(lambda2a, lambda2b)
	      IF(coordinate_flag == 2) THEN
	 	      dt_temp1 = MIN(cfl * MIN(dr1(j) / lambda1, r1(j) * dx1_z / lambda1) ,dt_temp1)
	      ELSE
	    	   dt_temp1 = MIN(cfl * dx1 / lambda1, dt_temp1)
	      END IF
	    ENDIF

      ENDDO
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we find the minimum time constrained by NM sector
! For non-MHD case
! For normal matter
DO k = length_step_z_min_part_2, length_step_z_part_2, 1
   DO j = 1, length_step_r_part_2, 1

      ! Only grid with density above threshold density is counted
      IF(rho2(j,k) > rho2_a) THEN
         lambda2a = ABS(vel2_r(j,k)) + cs2(j,k)
         lambda2b = ABS(vel2_z(j,k)) + cs2(j,k)
         lambda2 = MAX(lambda2a, lambda2b)
	      IF(coordinate_flag == 2) THEN
	 	      dt_temp1 = MIN(cfl * MIN(dr2(j) / lambda2, r2(j) * dx2_z / lambda2) ,dt_temp1)
	      ELSE
            dt_temp2 = MIN(dt_temp2, cfl * dx2 / lambda2)
	      END IF
      ENDIF

   ENDDO
ENDDO

! Only the minimum one is chosen
dt = MIN(dt_temp1, dt_temp2)

END SUBROUTINE FindDt