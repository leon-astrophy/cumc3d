!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
! Welcome to CUMC3D_Ver2022-1.42 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This code is a three-dimensional hydrodynamic, multiphysics code 
! which can simulate a broad range of astrophysical problems
! Written by Leung Shing Chi in 2016, major update by H.S. Leon Chan
! Questions and comments towards the code are welcome.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM CUMC3D
USE MHD_module
USE DEFINITION
IMPLICIT NONE

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate

CALL system_clock(count_rate=cr)
rate = REAL(cr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set number of threads !
CALL OMP_SET_NUM_THREADS(16)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*)
WRITE(*,*) '--------------------------------'
WRITE(*,*) '-Welcome to CUMC3D_Ver2022-1.42-'
WRITE(*,*) '--------------------------------'
WRITE(*,*)
WRITE(*,*) '-----------------------------------------'
WRITE(*,*) '- First written by Ka Wing Wong in 2010 -'
WRITE(*,*) '- Updated by Shing Chi Leung in 2016    -'
WRITE(*,*) '- Rewritten by H.S. Leon Chan in 2022   -'
WRITE(*,*) '-----------------------------------------'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial settings !

! Initialise !
n_iter = 0
n_step = 0
global_time = 0.0D0

! setup initial conditions !
CALL initial_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for openacc !

#ifdef GPU
CALL POPULATE_DEVICE
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! do initial updates !
CALL initial_update

! print grid variables !
CALL print_grid

! print primitive profiles !
CALL print_hydroprofile

! Divergence B !
CALL OPENFILE_MHD
CALL FIND_DIVB

! Custom quantities !
CALL OPENFILE_CUSTOM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main rungekutta loop !

! Print out !
WRITE (*,*) 'We now start the time evolution'
WRITE (*,*)
WRITE (*,*) '------------------------------------'
WRITE (*,*) 'iteration, timestep, simulation time'
WRITE (*,*) '------------------------------------'

CALL system_clock(time_start)

! Loop !
n_iter = 1
DO while (global_time < total_time)

  ! Find time step !
  CALL finddt

  ! Rungekutta step !
  CALL Rungekutta(n_step)

  ! Update !
  n_step = n_step + 1
  global_time = global_time + dt

  ! Print out !
  WRITE (*,*) n_step, dt, global_time

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!s	
  ! Section for data I/O !
  if(ABS(global_time - output_profiletime_last) >= output_profiletime .or. output_file .eqv. .true.) then

    ! print out !
	  output_profiletime_last = global_time
    call print_hydroprofile 

    ! divergence B !
    CALL FIND_DIVB

    ! custom analysis !
    CALL CUSTOM_ANALYSIS

    ! Increase the outfile index !
    n_iter = n_iter + 1

  end if

END DO

CALL system_clock(time_end)
!#ifdef DEBUG
WRITE(*,*) 'total time = ', REAL(time_end - time_start) / rate
!#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output again

call print_hydroprofile
CALL FIND_DIVB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for open acc

#ifdef GPU
CALL CLEAR_DEVICE
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out message 

WRITE(*,*) '----------------'
WRITE(*,*) '-Simulation End-'
WRITE(*,*) '----------------'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM
