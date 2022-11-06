!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! This module contains all the tools for including 
! magnetic field into hydrodynamics, which give
! rises to MHD. It uses the correction scheme
! to reduce the divergence problem of B field
! Written by Leung Shing Chi in 2016
!
! This module contains the following subroutines:
! 1.  subroutine buildMHD
! 2.  subroutine destoryMHD
! 3.  subroutine findcH
! 4.  subroutine getMHD
! 5.  subroutine getMHD2
! 6.  subroutine RelaxB2
! 7.  subroutine outputMHD_log
! 8.  subroutine outputMHD_profile(n)
! 9.  subroutine findDivB
! 10. subroutine relaxB
! 11. subroutine FindMaxDivB
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MHD_module
use definition   
implicit none

! The call-code for the magnetic field
integer :: iBr
integer :: iBz
integer :: iBp
integer :: iPsi

!
real (DP), allocatable, dimension(:,:) :: Br
real (DP), allocatable, dimension(:,:) :: Bz
real (DP), allocatable, dimension(:,:) :: Bp
real (DP), allocatable, dimension(:,:) :: psi
real (DP), allocatable, dimension(:,:) :: divB

! The highest effective speed
real (DP) :: c_h = 0.0D0
real (DP) :: c_p = 0.0D0

! Maximum divergence B
real (DP) :: maxDivB

! The averaged divergance B
real (DP) :: normB

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine allocates the array necessary for
   ! the calculation of MHD
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine buildMHD
   use definition
   implicit none

   allocate(Br(-4:length_step_r+5, -4:length_Step_z+5))
   allocate(Bz(-4:length_step_r+5, -4:length_Step_z+5))
   allocate(Bp(-4:length_step_r+5, -4:length_Step_z+5))
   allocate(Psi(-4:length_step_r+5, -4:length_Step_z+5))
   allocate(divB(-4:length_step_r+5, -4:length_Step_z+5))

   end subroutine buildMHD
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine deallocate the array to save up 
   ! space and avoid unnecessary blocking of memory
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine destroyMHD
   use definition
   implicit none

   deallocate(Br)
   deallocate(Bz)
   deallocate(Bp)
   deallocate(Psi)
   deallocate(DivB)
   
   end subroutine destroyMHD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the effective scaled
   ! CFL number which are important for the correction
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine findcH
   use definition
   implicit none

   c_h = cfl * dx / dt
   c_p = DSQRT(0.18D0 * c_h)

   end subroutine findCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine assigns the initial B-field
   ! for the evolution
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getMHD
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Assign the value
   do k = 1, length_Step_z, 1
      do j = 1, length_step_r, 1
	 Br(j,k) = 0.0D0
	 Bz(j,k) = 1.0D-6
	 Bp(j,k) = 0.0D0
      enddo
   enddo

   ! Copy the cell to ghost cell
   call boundary1D(Br, oddR)
   call boundary1D(Bz, oddZ)
   call boundary1D(Bp, oddR)

   end subroutine getMHD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine assign the initial correction term
   ! according to Migone (2010)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine GetMHD2
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k
  
   ! Assign the values
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1

         psi(j,k) = -c_p**2 * (0.5D0 * (Br(j+1,k) - Br(j-1,k)) /dx + 0.5D0 * (Bz(j,k+1) - Bz(j,k-1)) / dx)

      enddo
   enddo

   ! Copy the cell to ghost cells
   call boundary1D(psi, even)

   end subroutine GetMHD2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine updates the local correction 
   ! factor.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine RelaxB2
   use definition
   implicit none

   ! dummy variables
   integer :: j, k

   ! Update the factor
   do k = -4, length_step_z_part + 5, 1
      do j = -4, length_step_r_part + 5, 1

 	 psi(j,k) = psi(j,k) * EXP(-dt * c_h ** 2 / c_p ** 2)

      enddo
   enddo

   ! Copy the new result to ghost cells
   call boundary1d(psi, even)

   end subroutine RelaxB2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find the global quantity and 
   ! output the log file
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine outputMHD_log
   use definition
   implicit none

   ! First find the maximum divergence term
   call findMaxDivB

   ! Then do all the output 
   open(unit=601, file='Star_WENO_CentralB_0.dat', action='write', position='append')
   write(601,100) global_time, Br(1,1), Bz(1,1), Bp(1,1)
   close(601)

   open(unit=601, file='Star_WENO_MaxDivB_0.dat', action='write', position='append')
   write(601,100) global_time, maxDivB, normB
   close(601)

   100 format(5ES18.8)

   end subroutine outputMHD_log

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine opens the necessary file and 
   ! output the MHD profile
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine outputMHD_profile(n)
   use definition
   implicit none
  
   ! dummy variables
   integer :: j, k, n
   
   ! Package for file name
   INTEGER :: fileno_len
   CHARACTER (len = 256) :: fileno   

   ! Assign the file name
   WRITE (fileno, *) n
   fileno = ADJUSTL (fileno)
   fileno_len = LEN_TRIM (fileno)

   ! Do the output
   open(unit=601, file='Star_WENO_Br_'//fileno (1 : fileno_len)//'.dat', action='write', position='append')
   write(601,*) n, global_time
   write(601,*) length_step_r_part, length_step_z_part
   write(601,*) dx
   do k = length_step_z_part, 1, -1
      write(601,100) (Br(j,k), j = 1, length_step_r_part)
   enddo
   write(601,*)
   close(601)

   open(unit=601, file='Star_WENO_Bz_'//fileno (1 : fileno_len)//'.dat', action='write', position='append')
   write(601,*) n, global_time
   write(601,*) length_step_r_part, length_step_z_part
   write(601,*) dx
   do k = length_step_z_part, 1, -1
      write(601,100) (Bz(j,k), j = 1, length_step_r_part)
   enddo
   write(601,*)
   close(601)

   open(unit=601, file='Star_WENO_Bp_'//fileno (1 : fileno_len)//'.dat', action='write', position='append')
   write(601,*) n, global_time
   write(601,*) length_step_r_part, length_step_z_part
   write(601,*) dx
   do k = length_step_z_part, 1, -1
      write(601,100) (Bp(j,k), j = 1, length_step_r_part)
   enddo
   write(601,*)
   close(601)

   open(unit=601, file='Star_WENO_Psi_'//fileno (1 : fileno_len)//'.dat', action='write', position='append')
   write(601,*) n, global_time
   write(601,*) length_step_r_part, length_step_z_part       
   write(601,*) dx
   do k = length_step_z_part, 1, -1
      write(601,100) (DSQRT(Br(j,k)**2 + Bz(j,k)**2 + Bp(j,k)**2), j = 1, length_step_r_part)
   enddo          
   write(601,*)
   close(601)

   open(unit=601, file='Star_WENO_DivB_'//fileno (1 : fileno_len)//'.dat', action='write', position='append')
   write(601,*) n, global_time
   write(601,*) length_step_r_part, length_step_z_part
   write(601,*) dx   
   do k = length_step_z_part, 1, -1
      write(601,100) (divB(j,k), j = 1, length_step_r_part)
   enddo          
   write(601,*)
   close(601)

   100 format(750ES18.8)

   end subroutine outputMHD_profile

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the divergence B
   ! which is suppose to be zero for ideal case
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine findDivB
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Calculate by simple discrete differentiation
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1
         divB(j,k) = 0.5D0 * (Br(j+1,k) - Br(j-1,k)) / dx + 0.5D0 * (Bz(j,k+1) - Bz(j,k-1)) / dx
      enddo
   enddo

   end subroutine findDivB

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine 
   ! 
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine RelaxB
   use definition
   implicit none

   integer :: j, k
   real (DP), dimension(:,:,:), allocatable :: Bmid_r, Bmid_z

   allocate(Bmid_r(-4:length_step_r+5, -4:length_step_z+5, 2))
   allocate(Bmid_z(-4:length_step_r+5, -4:length_step_z+5, 2))

   Bmid_r(:,:,:) = 0.0D0
   Bmid_z(:,:,:) = 0.0D0

   ! First compute the Br-flux using the right-eigenstate
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1
	 Bmid_r(j,k,1) = Br(j,k) + (0.5D0 * (Br(j+1,k) - Br(j,k)) - 0.5D0 * (psi(j+1,k) - psi(j,k)) / c_h)
	 Bmid_r(j,k,2) = psi(j,k) + (0.5D0 * (psi(j+1,k) - psi(j,k)) - 0.5D0 * c_h * (Br(j+1,k) - Br(j,k)))
      enddo
   enddo

   ! Second compute the Bz-flux using the right-eigenstate
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1
	 Bmid_z(j,k,1) = Bz(j,k) + (0.5D0 * (Bz(j,k+1) - Bz(j,k)) - 0.5D0 * (psi(j,k+1) - psi(j,k)) / c_h)
         Bmid_z(j,k,2) = psi(j,k) + (0.5D0 * (psi(j,k+1) - psi(j,k)) - 0.5D0 * c_h * (Bz(j,k+1) - Bz(j,k)))
      enddo
   enddo

   ! COmpute the corrected B-field and psi
   ! No correction in Bp due to symmetry
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1

	 ! Relax the B field according to the prescription
	 Br(j,k) = Br(j,k) - dt * (Bmid_r(j+1,k,2) - Bmid_r(j,k,2)) / dx
	 Bz(j,k) = Bz(j,k) - dt * (Bmid_z(j,k+1,2) - Bmid_z(j,k,2)) / dx
	 psi(j,k) = psi(j,k) - dt * c_h ** 2 * &
		    ((Bmid_r(j+1,k,1) - Bmid_r(j,k,1)) / dx + &
		     (Bmid_z(j,k+1,1) - Bmid_z(j,k,1)) / dx)

      enddo
   enddo

   ! Copy the results to ghost cells
   call boundary1D(Br, even)
   call boundary1D(Bz, even)
   call boundary1D(psi, even)

   ! Deallocate everywhere 
   deallocate(Bmid_r)
   deallocate(Bmid_z)

   end subroutine RelaxB

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find the maximum diV B
   ! which is used for output
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine FindMaxDivB
   use definition
   implicit none

   ! dummy variable
   integer :: j, k

   ! Local variable
   real (DP) :: absDivB

   ! First find divB everywhere
   call findDivB

   ! Initialization
   normB = 0.0D0
   maxDivB = 0.0D0

   ! Do the search
   do k = 1, length_step_z, 1
      do j = 1, length_step_r, 1

	 absDivB = ABS(divB(j,k))
	 normB = normB + absDivB
         maxDivB = MAX(maxDivB, absDivB)

      enddo
   enddo

   ! Do the average
   normB = normB / (DBLE(length_step_r) * DBLE(length_step_z))

   end subroutine FindMaxDivB




end module MHD_module

