!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the global stellar data
! Written by Leung Shing CHi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUT_LOG
USE DEFINITION
USE TURB_MODULE
USE NUSPEC_MODULE
USE LEVELSET_MODULE
USE GW_MODULE, ONLY : OUTPUTGRAVQ
USE HELMEOS_MODULE, ONLY : OUTPUT_XMASS
IMPLICIT NONE

! dummy variable
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the information first
CALL PREPARE_LOGDATA()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output mandatory data

! DM log files 
If (DM_flag == 1) THEN
	WRITE (301, 701) global_time, mass1
	WRITE (302, 702) global_time, energy1, energy1_kin, energy1_int, energy1_pot
	WRITE (303, 702) global_time, centralrho1, maxval(rho1)
END IF

! NM log files 
WRITE (401, 701) global_time, mass2
WRITE (402, 702) global_time, energy2, energy2_kin, energy2_int, energy2_pot	
WRITE (403, 701) global_time, centralrho2, MAXVAL(rho2)
WRITE (404, 701) global_time, centraltemperature, MAXVAL(temp2)
WRITE (405, 701) global_time, ye2(1,1), MINVAL(ye2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output other quantities when needed
if(turb_flag == 1) THEN
	call outputturb_log
END IF
if(nuspec_flag == 1) THEN
	call output_nuphi
END IF
if(helmeos_flag == 1) THEN
	CALL output_xmass()
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Format !
701 FORMAT (F23.15, 2ES18.10)
702 FORMAT (F23.15, 30ES18.10)
703 FORMAT (F23.15, 30ES18.10)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares all the information which 
! are presented in the log file.
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PREPARE_LOGDATA() 
USE DEFINITION
USE GW_MODULE, ONLY : FINDGRAVQ
USE NUSPEC_MODULE, ONLY : FINDNUSPEC
USE HELMEOS_MODULE, ONLY : FINDNEUTRINOLOSS, FINDTOTALNEUTRINOLOSS, FINDLUMINOSITY, FINDASH
USE TURB_MODULE
IMPLICIT NONE

! Find Epsilon !
CALL FINDEPSILON

! Find related to hydro !
CALL FINDMASS
CALL FINDENERGY
CALL FINDCENTRALDENSITY
CALL FINDCENTRALTEMPERATURE

! Related to nuclear reaction or finite temperature eos
if(helmeos_flag == 1) THEN
	CALL FINDNEUTRINOLOSS
	CALL FINDTOTALNEUTRINOLOSS
 	CALL FINDLUMINOSITY
END IF

! Other stuff !
if(turb_flag == 1) THEN
	call findturbenergy
END IF
if(flame_flag == 1) THEN
	CALL FINDASH
END IF
if(nuspec_flag == 1) THEN
	call findNuSpec
END IF

END SUBROUTINE