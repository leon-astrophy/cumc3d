!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
! Welcome to CUMC3D_Ver2022-1.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Iteration step number
INTEGER :: n_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*)
WRITE(*,*) '-------------------------------'
WRITE(*,*) '-Welcome to CUMC3D_Ver2022-1.0-'
WRITE(*,*) '-------------------------------'
WRITE(*,*)
WRITE(*,*) '-----------------------------------------'
WRITE(*,*) '- First written by Ka Wing Wong in 2010 -'
WRITE(*,*) '- Updated by Shing Chi Leung in 2016    -'
WRITE(*,*) '- Rewritten by H.S. Leon Chan in 2022   -'
WRITE(*,*) '-----------------------------------------'
WRITE(*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial settings !

CALL Initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out message 

WRITE(*,*) '----------------'
WRITE(*,*) '-Simulation End-'
WRITE(*,*) '----------------'

END PROGRAM