!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
! Welcome to CUMC3D-Ver1.65 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
USE MHD_MODULE
USE DEFINITION
IMPLICIT NONE

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set number of threads !
CALL OMP_SET_NUM_THREADS(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*)
WRITE(*,*) '---------------------------'
WRITE(*,*) '-Welcome to CUMC3D-Ver1.65-'
WRITE(*,*) '---------------------------'
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

! print primitive profiles !
CALL print_hydroprofile

! Divergence B !
#ifdef DIVB
CALL OPENFILE_MHD
CALL FIND_DIVB
#endif

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

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! Loop !
n_iter = 1
DO while (global_time < total_time)

  ! Find time step !
  CALL finddt

  ! Adjust !
  IF(global_time + dt > total_time) THEN
    dt = total_time - global_time
  END IF

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
    #ifdef DIVB
    CALL FIND_DIVB 
    #endif

    ! custom analysis !
    CALL CUSTOM_ANALYSIS

    ! Increase the outfile index !
    n_iter = n_iter + 1

  end if

END DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'total time = ', REAL(time_end - time_start) / rate
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output again

call print_hydroprofile
#ifdef DIVB
CALL FIND_DIVB
#endif

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
