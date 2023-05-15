!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the PPM (Piecewise Parabolic Method) Module that reconstruct   
! interface values of primitive variables. Three different appoarch 
! can be chosen, and they differ by how to handle monotone conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPM_MODULE
USE DEFINITION
IMPLICIT NONE

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the reconstructions weights for 3rd order PPM reconstruction  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_WEIGHT
USE DEFINITION
IMPLICIT NONE

! Matrix !
REAL*8, DIMENSION(1:5,1:5) :: a_matrix, a_inverse

! Row vector !
REAL*8, DIMENSION(1:5) :: b_vector

! work array for LAPACK
REAL*8, DIMENSION(1:5) :: work    

! pivot indices   
INTEGER, DIMENSION(1:5) :: ipiv

! Power for jacobian derivatives !
REAL*8 :: pow_x, pow_y, pow_z

! Angular differences
REAL*8 :: dmu, dmubar

! Integer !
INTEGER :: i, j, k, info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initailize all reconstruction weight !

! X-direction !
lxm2 = 0.0d0
lxm1 = 0.0d0
lxc = 0.0d0
lxp1 = 0.0d0
lxp2 = 0.0d0
rxm2 = 0.0d0
rxm1 = 0.0d0
rxc = 0.0d0
rxp1 = 0.0d0
rxp2 = 0.0d0
hpx = 0.0d0
hmx = 0.0d0

! Y-direction !
lym2 = 0.0d0
lym1 = 0.0d0
lyc = 0.0d0
lyp1 = 0.0d0
lyp2 = 0.0d0
rym2 = 0.0d0
rym1 = 0.0d0
ryc = 0.0d0
ryp1 = 0.0d0
ryp2 = 0.0d0
hpy = 0.0d0
hmy = 0.0d0

! Z-direction !
lzm2 = 0.0d0
lzm1 = 0.0d0
lzc = 0.0d0
lzp1 = 0.0d0
lzp2 = 0.0d0
rzm2 = 0.0d0
rzm1 = 0.0d0
rzc = 0.0d0
rzp1 = 0.0d0
rzp2 = 0.0d0
hpz = 0.0d0
hmz = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, assign the jacobian derivative power !

IF(coordinate_flag == 0) THEN
	pow_x = 0.0d0
	pow_y = 0.0d0
	pow_z = 0.0d0
ELSEIF(coordinate_flag == 1) THEN
	pow_x = 1.0d0
	pow_y = 0.0d0
	pow_z = 0.0d0
ELSEIF(coordinate_flag == 2) THEN
	pow_x = 2.0d0
	pow_y = 0.0d0
	pow_z = 0.0d0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! h-parameter !

! x-direction !
IF(coordinate_flag == 0) THEN
	DO j = 0, nx_2 + 1
		hpx(j) = 3.0d0
		hmx(j) = 3.0d0
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO j = 0, nx_2 + 1
		hpx(j) = 3.0d0 + dx2(j)/(2.0d0*x2(j))
		hmx(j) = 3.0d0 - dx2(j)/(2.0d0*x2(j))
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO j = 0, nx_2 + 1
		hpx(j) = 3.0d0 + 2.0d0*dx2(j)*(dx2(j) + 10.0d0*x2(j))/(20.0d0*x2(j)**2 + dx2(j)**2)
		hmx(j) = 3.0d0 + 2.0d0*dx2(j)*(dx2(j) - 10.0d0*x2(j))/(20.0d0*x2(j)**2 + dx2(j)**2)
	END DO
END IF

! y-direction !
IF(coordinate_flag == 2) THEN
	DO j = 0, ny_2 + 1
		dmu = COS(yF2(k-1)) - COS(yF2(k))
		dmubar = SIN(yF2(k-1)) - SIN(yF2(k))
		hpy(j) = dy2(k)*(dmubar + dy2(k)*COS(yF2(k)))/(dy2(k)*(SIN(yF2(k-1)) + SIN(yF2(k))) - 2.0d0*dmu)
		hmy(j) = -dy2(k)*(dmubar + dy2(k)*COS(yF2(k-1)))/(dy2(k)*(SIN(yF2(k-1)) + SIN(yF2(k))) - 2.0d0*dmu)
	END DO
ELSE
	DO j = 0, ny_2 + 1
		hpy(j) = 3.0d0
		hmy(j) = 3.0d0
	END DO
END IF

! Z-direction !
DO j = 0, nz_2 + 1
	hpz(j) = 3.0d0
	hmz(j) = 3.0d0
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do for the x-direction !

! Loop over the complete domain !
DO j = 0, nx_2 + 1

	! Get matrix !
	CALL PPMBETA_SIMPLE(x_dir, pow_x, j, a_matrix)

	! Backup !
	a_inverse(:,:) = a_matrix(:,:)		

	! Invert matrix !
	CALL DGETRF(5,5,a_inverse,5,ipiv,info)
	CALL DGETRI(5,a_inverse,5,ipiv,work,5,info)

	! Assign column vector !
	DO i = 1, 5
		b_vector(i) = xF2(j)**(i - 1)
	END DO

	! Get weight !
	DO i = 1, 5
		rxm2(j) = rxm2(j) + a_inverse(1, i)*b_vector(i)
		rxm1(j) = rxm1(j) + a_inverse(2, i)*b_vector(i)
		rxc(j) = rxc(j) + a_inverse(3, i)*b_vector(i)
		rxp1(j) = rxp1(j) + a_inverse(4, i)*b_vector(i)
		rxp2(j) = rxp2(j) + a_inverse(5, i)*b_vector(i)
	END DO

	! Assign column vector !
	DO i = 1, 5
		b_vector(i) = xF2(j-1)**(i - 1)
	END DO

	! Get weight !
	DO i = 1, 5
		lxm2(j) = lxm2(j) + a_inverse(1, i)*b_vector(i)
		lxm1(j) = lxm1(j) + a_inverse(2, i)*b_vector(i)
		lxc(j) = lxc(j) + a_inverse(3, i)*b_vector(i)
		lxp1(j) = lxp1(j) + a_inverse(4, i)*b_vector(i)
		lxp2(j) = lxp2(j) + a_inverse(5, i)*b_vector(i)
	END DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do for the z-direction !

! Loop over the complete domain !
DO j = 0, nz_2 + 1

	! Get matrix !
	CALL PPMBETA_SIMPLE(z_dir, pow_z, j, a_matrix)

	! Backup !
	a_inverse(:,:) = a_matrix(:,:)		

	! Invert matrix !
	CALL DGETRF(5,5,a_inverse,5,ipiv,info)
	CALL DGETRI(5,a_inverse,5,ipiv,work,5,info)

	! Assign column vector !
	DO i = 1, 5
		b_vector(i) = zF2(j)**(i - 1)
	END DO

	! Get weight !
	DO i = 1, 5
		rzm2(j) = rzm2(j) + a_inverse(1, i)*b_vector(i)
		rzm1(j) = rzm1(j) + a_inverse(2, i)*b_vector(i)
		rzc(j) = rzc(j) + a_inverse(3, i)*b_vector(i)
		rzp1(j) = rzp1(j) + a_inverse(4, i)*b_vector(i)
		rzp2(j) = rzp2(j) + a_inverse(5, i)*b_vector(i)
	END DO

	! Assign column vector !
	DO i = 1, 5
		b_vector(i) = zF2(j-1)**(i - 1)
	END DO

	! Get weight !
	DO i = 1, 5
		lzm2(j) = lzm2(j) + a_inverse(1, i)*b_vector(i)
		lzm1(j) = lzm1(j) + a_inverse(2, i)*b_vector(i)
		lzc(j) = lzc(j) + a_inverse(3, i)*b_vector(i)
		lzp1(j) = lzp1(j) + a_inverse(4, i)*b_vector(i)
		lzp2(j) = lzp2(j) + a_inverse(5, i)*b_vector(i)
	END DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do for the y-direction !

! Choose coordinate system !
IF(coordinate_flag == 2) THEN

	! Loop over the complete domain !
	DO j = 0, ny_2 + 1

		! Get matrix !
		CALL PPMBETA_POLAR(j, a_matrix)
			
		! Backup !
		a_inverse(:,:) = a_matrix(:,:)		

		! Invert matrix !
		CALL DGETRF(5,5,a_inverse,5,ipiv,info)
		CALL DGETRI(5,a_inverse,5,ipiv,work,5,info)

		! Assign column vector !
		DO i = 1, 5
			b_vector(i) = yF2(j)**(i - 1)
		END DO

		! Get weight !
		DO i = 1, 5
			rym2(j) = rym2(j) + a_inverse(1, i)*b_vector(i)
			rym1(j) = rym1(j) + a_inverse(2, i)*b_vector(i)
			ryc(j) = ryc(j) + a_inverse(3, i)*b_vector(i)
			ryp1(j) = ryp1(j) + a_inverse(4, i)*b_vector(i)
			ryp2(j) = ryp2(j) + a_inverse(5, i)*b_vector(i)
		END DO

		! Assign column vector !
		DO i = 1, 5
			b_vector(i) = yF2(j-1)**(i - 1)
		END DO

		! Get weight !
		DO i = 1, 5
			lym2(j) = lym2(j) + a_inverse(1, i)*b_vector(i)
			lym1(j) = lym1(j) + a_inverse(2, i)*b_vector(i)
			lyc(j) = lyc(j) + a_inverse(3, i)*b_vector(i)
			lyp1(j) = lyp1(j) + a_inverse(4, i)*b_vector(i)
			lyp2(j) = lyp2(j) + a_inverse(5, i)*b_vector(i)
		END DO

	END DO

ELSE

	! Loop over the complete domain !
	DO j = 0, ny_2 + 1

		! Get matrix !
		CALL PPMBETA_SIMPLE(y_dir, pow_y, j, a_matrix)
			
		! Backup !
		a_inverse(:,:) = a_matrix(:,:)		

		! Invert matrix !
		CALL DGETRF(5,5,a_inverse,5,ipiv,info)
		CALL DGETRI(5,a_inverse,5,ipiv,work,5,info)

		! Assign column vector !
		DO i = 1, 5
			b_vector(i) = yF2(j)**(i - 1)
		END DO

		! Get weight !
		DO i = 1, 5
			rym2(j) = rym2(j) + a_inverse(1, i)*b_vector(i)
			rym1(j) = rym1(j) + a_inverse(2, i)*b_vector(i)
			ryc(j) = ryc(j) + a_inverse(3, i)*b_vector(i)
			ryp1(j) = ryp1(j) + a_inverse(4, i)*b_vector(i)
			ryp2(j) = ryp2(j) + a_inverse(5, i)*b_vector(i)
		END DO

		! Assign column vector !
		DO i = 1, 5
			b_vector(i) = yF2(j-1)**(i - 1)
		END DO

		! Get weight !
		DO i = 1, 5
			lym2(j) = lym2(j) + a_inverse(1, i)*b_vector(i)
			lym1(j) = lym1(j) + a_inverse(2, i)*b_vector(i)
			lyc(j) = lyc(j) + a_inverse(3, i)*b_vector(i)
			lyp1(j) = lyp1(j) + a_inverse(4, i)*b_vector(i)
			lyp2(j) = lyp2(j) + a_inverse(5, i)*b_vector(i)
		END DO

	END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix element for simple geometry i.e. dv/dx = x**(power)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMBETA_SIMPLE(dir_in, pow_in, j_in, matrix_out)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: j_in

! integer !
INTEGER, INTENT(IN) :: dir_in

! Input power !
REAL*8, INTENT(IN) :: pow_in

! Ouput matrix !
REAL*8, INTENT(OUT), DIMENSION(1:5,1:5) :: matrix_out

! Integer !
INTEGER :: i, j, k
INTEGER :: n_temp, s_temp

! Left and right interface coordinate  !
REAL*8 :: left, right

! Temporal variables !
REAL*8 :: numerator, denominator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
matrix_out = 0.0d0

! Loop over the whole matrix !
DO i = 1, 5
	DO j = 1, 5
		n_temp = i - 1
		s_temp = -2 + j - 1
		numerator = DBLE(pow_in + 1)
		denominator = DBLE(pow_in + n_temp + 1)
		IF(dir_in == x_dir) THEN
			left = xF2(j_in + s_temp - 1)
			right = xF2(j_in + s_temp)
		ELSEIF(dir_in == y_dir) THEN
			left = yF2(j_in + s_temp - 1)
			right = yF2(j_in + s_temp)
		ELSEIF(dir_in == z_dir) THEN
			left = zF2(j_in + s_temp - 1)
			right = zF2(j_in + s_temp)
		END IF
		matrix_out(i, j) = (numerator/denominator) & 
						*(right**(denominator) - left**(denominator)) &
						/(right**(numerator) - left**(numerator))
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix element for polar geometry i.e. dv/dx = sin(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMBETA_POLAR(j_in, matrix_out)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: j_in

! Ouput matrix !
REAL*8, INTENT(OUT), DIMENSION(1:5,1:5) :: matrix_out

! Integer !
INTEGER :: i, j, k
INTEGER :: n_temp, s_temp

! Left and right interface coordinate  !
REAL*8 :: left, right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
matrix_out = 0.0d0

! Loop over the whole matrix !
DO i = 1, 5
	DO j = 1, 5
		n_temp = i - 1
		s_temp = -2 + j - 1
		left = yF2(j_in + s_temp - 1)
		right = yF2(j_in + s_temp)
		DO k = 0, n_temp
			matrix_out(i, j) = matrix_out(i, j) + factorial(n_temp)/factorial(k) &
			*(left**(DBLE(n_temp - k))*COS(left + DBLE(k)*0.5d0*pi) - right**(DBLE(n_temp - k))*COS(right + DBLE(k)*0.5d0*pi))
		END DO
		matrix_out(i, j) = matrix_out(i, j)/(COS(left) - COS(right))
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	! Factorial !
	REAL*8 Function factorial(n_in)
	INTEGER, INTENT(IN) :: n_in
	INTEGER :: i
	factorial = 1.0d0
	DO i = 2, n_in
		factorial = factorial * i
	END DO
	factorial = DBLE(factorial)
	END function factorial

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the x-direction, and get states at the x-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_reconx
USE MHD_MODULE
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2 
		DO j = nx_min_2 - 1, nx_part_2 + 1
			DO i = imin2, imax2
				CALL PPM (lxm2(j), lxm1(j), lxc(j), lxp1(j), lxp2(j), & 
									rxm2(j), rxm1(j), rxc(j), rxp1(j), rxp2(j), & 
									hmx(j), hpx(j), & 
									prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), & 
									primR2(i,j-1,k,l), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2 
			DO j = nx_min_2 - 1, nx_part_2 + 1
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL PPM (lxm2(j), lxm1(j), lxc(j), lxp1(j), lxp2(j), & 
									rxm2(j), rxm1(j), rxc(j), rxp1(j), rxp2(j), & 
									hmx(j), hpx(j), & 
									cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), & 
									cs2R(j-1,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2 
			DO j = nx_min_2 - 1, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the y-direction, and get states at the y-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_recony
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2, nz_part_2 
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2
				CALL PPM (lym2(k), lym1(k), lyc(k), lyp1(k), lyp2(k), & 
									rym2(k), rym1(k), ryc(k), ryp1(k), ryp2(k), & 
									hmy(k), hpy(k), & 	
									prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), & 
									primR2(i,j,k-1,l), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2, nz_part_2 
		DO k = ny_min_2 - 1, ny_part_2 + 1
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL PPM (lym2(k), lym1(k), lyc(k), lyp1(k), lyp2(k), & 
									rym2(k), rym1(k), ryc(k), ryp1(k), ryp2(k), & 
									hmy(k), hpy(k), &
									cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), &
									cs2R(j,k-1,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2, nz_part_2 
		DO k = ny_min_2 - 1, ny_part_2
			DO j = nx_min_2, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the z-direction, and get states at the z-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_reconz
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2
				CALL PPM (lzm2(l), lzm1(l), lzc(l), lzp1(l), lzp2(l), & 
									rzm2(l), rzm1(l), rzc(l), rzp1(l), rzp2(l), & 
									hmz(l), hpz(l), &
									prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), & 
									primR2(i,j,k,l-1), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2 - 1, nz_part_2 + 1
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL PPM (lzm2(l), lzm1(l), lzc(l), lzp1(l), lzp2(l), & 
									rzm2(l), rzm1(l), rzc(l), rzp1(l), rzp2(l), & 
									hmz(l), hpz(l), &
									cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), & 
									cs2R(j,k,l-1), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2 - 1, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Getting interface primitive variables using the PPM method
! with the variants proposed by Mignone 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM (lm2, lm1, lc, lp1, lp2, & 
							  rm2, rm1, rc, rp1, rp2, &
								hm1, hp1, &
								vm2, vm1, vc, vp1, vp2, &
								vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL*8, INTENT (IN) :: hm1, hp1
REAL*8, INTENT (IN) :: lm2, lm1, lc, lp1, lp2
REAL*8, INTENT (IN) :: rm2, rm1, rc, rp1, rp2
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL*8, INTENT (OUT) :: vp_out, vm_out

! Temporal variables !
REAL*8 :: dqp, dqm
REAL*8 :: condition, condition1, condition2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the interpolant !
vm_out = lm2*vm2 + lm1*vm1 + lc*vc + lp1*vp1 + lp2*vp2
vp_out = rm2*vm2 + rm1*vm1 + rc*vc + rp1*vp1 + rp2*vp2

! Strict monotonic contrains !
vp_out = MIN(vp_out, MAX(vc, vp1))
vp_out = MAX(vp_out, MIN(vc, vp1))
vm_out = MIN(vm_out, MAX(vc, vm1))
vm_out = MAX(vm_out, MIN(vc, vm1))

! Differences !
dqp = vp_out - vc
dqm = vm_out - vc
condition = dqp*dqm
condition1 = ABS(dqp) - (hm1 + 1.0d0)/(hp1 - 1.0d0)*ABS(dqm)
condition2 = ABS(dqm) - (hp1 + 1.0d0)/(hm1 - 1.0d0)*ABS(dqp)

! Parabolic limiters, for the right state !
IF(condition >= 0.0d0) THEN
	vp_out = vc
ELSEIF(condition < 0.0d0 .AND. condition1 >= 0.0D0) THEN
	vp_out = vc - (hm1 + 1.0d0)/(hp1 - 1.0d0)*dqm
ELSE
	vp_out = vc + dqp
END IF

! Parabolic limiters, for the left state !
IF(condition >= 0.0d0) THEN
	vm_out = vc
ELSEIF(condition < 0.0d0 .AND. condition2 >= 0.0D0) THEN
	vm_out = vc - (hp1 + 1.0d0)/(hm1 - 1.0d0)*dqp
ELSE
	vm_out = vc + dqm
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
