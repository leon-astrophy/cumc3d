!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUMC3D_Ver2022-1.0 (Last Modified: Leon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Foreword:
! Developed by Leung Shing Chi based on the WENO prototype developed by Wong Ka Wing in 2010
! c.f. Leung et al., MNRAS 454, 1238 (2015).
! Further developed by Leon H.S. Chan to extend to 3D, and include MHD/GRMHD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Constant and universal values 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Boundary flag notation for boundary conditions
INTEGER, PARAMETER :: even = 0, odd = 1

! Which direction of reconstruction
INTEGER, PARAMETER :: x_dir = 1, y_dir = 2, z_dir = 3

! Apply boundary conditions to the full/partial box
INTEGER, PARAMETER :: part = 0, full = 1	

! Which type of matter considering
INTEGER, PARAMETER :: dm_f = 1, nm_f = 2

! pi constant 
REAL*8, PARAMETER :: pi = 3.1415926535897932384626433832795D0

! define small number to avoid coordinate singularity !
REAL*8, PARAMETER :: small_num = 1.0d-50

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Core part of the simulation box
! This is where the system variables are controled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! 0 = Cartesian coordinate
! 1 = Cylindrical coordinate
! 2 = Spherical coordinate
INTEGER, PARAMETER :: coordinate_flag = 0

! Dimension, 1 = 1D, 2 = 2D, 3 = 3D
INTEGER, PARAMETER :: n_dim = 2

! Flag for boundary condition
! The boundary flag is defined by four scalar
! 1st one for x-inner boundary
! 2nd one for x-outer boundary
! 3rd one for y-inner boundary
! 4th one for y-outer boundary
! 5rd one for z-inner boundary
! 6th one for z-outer boundary
! 0 = periodic
! 1 = outgoing (1st derivative = 0)
! 2 = reflecting boundary (depend on scalar/vector)
! 3 = axis-symmetric
! 4 = equatorial-symmetric
INTEGER :: boundary_flag(6) = (/0,0,0,0,1,1/)

! Flag for simulating with the full box (extend to negative x,y,z)
LOGICAL, PARAMETER :: fullx_flag = .false.
LOGICAL, PARAMETER :: fully_flag = .false.
LOGICAL, PARAMETER :: fullz_flag = .false.

! Starting position of the grid !
REAL*8, PARAMETER :: x1_start = 0.0d0
REAL*8, PARAMETER :: y1_start = 0.0d0
REAL*8, PARAMETER :: z1_start = 0.0d0

! Ending position of the grid !
REAL*8, PARAMETER :: x1_end = 1.0d0
REAL*8, PARAMETER :: y1_end = 1.0d0
REAL*8, PARAMETER :: z1_end = 1.0d0

! The number of grid in the x,y,z direction for DM
INTEGER, PARAMETER :: nx_1 = 1
INTEGER, PARAMETER :: ny_1 = 1
INTEGER, PARAMETER :: nz_1 = 1

! Grid sizes for DM 
! Hint: The standard grid size (Note: 1 unit = 1.4774 km)
REAL*8, PARAMETER :: dx1_ini = (x1_end - x1_start)/DBLE(nx_1)	
REAL*8, PARAMETER :: dy1_ini = (y1_end - y1_start)/DBLE(ny_1)	
REAL*8, PARAMETER :: dz1_ini = (z1_end - z1_start)/DBLE(nz_1)	

! Starting position of the grid !
REAL*8, PARAMETER :: x2_start = 0.0d0
REAL*8, PARAMETER :: y2_start = 0.0d0
REAL*8, PARAMETER :: z2_start = 0.0d0

! Ending position of the grid !
REAL*8, PARAMETER :: x2_end = 1.0d0
REAL*8, PARAMETER :: y2_end = 1.0d0
REAL*8, PARAMETER :: z2_end = 1.0d0

! The number of grid in the x,y,z direction for NM
INTEGER, PARAMETER :: nx_2 = 128
INTEGER, PARAMETER :: ny_2 = 128
INTEGER, PARAMETER :: nz_2 = 1

! Grid sizes for NM
REAL*8, PARAMETER :: dx2_ini = (x2_end - x2_start)/DBLE(nx_2)	
REAL*8, PARAMETER :: dy2_ini = (y2_end - y2_start)/DBLE(ny_2)	
REAL*8, PARAMETER :: dz2_ini = (z2_end - z2_start)/DBLE(nz_2)	

! Cournat-Friedrich-Levy constant
! Defined as dt = cfl * dx / MAX(vel + cs)
REAL*8, PARAMETER :: cfl = 0.20D0			

! Maximum time to be simulated in the model
REAL*8 :: total_time = 120.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for core hydrodynamic solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the LF (two-shocks) Riemann solver !
LOGICAL, PARAMETER :: LF_flag = .false.

! Use the HLL Riemann solver !
LOGICAL, PARAMETER :: HLL_flag = .true.

! Use the HLLC Riemann solver !
LOGICAL, PARAMETER :: HLLC_flag = .false.

! Use the HLLD Riemann solver !
LOGICAL, PARAMETER :: HLLD_flag = .false.

! Use the TVD (Mignone 2014) reconstruction scheme with MC limiter !
LOGICAL, PARAMETER :: tvdmc_flag = .false.

! Use the TVD (Mignone 2014) reconstruction scheme with Van-Leer limiter !
LOGICAL, PARAMETER :: tvdvl_flag = .true.

! Use the PPM (Mignone 2014) reconstruction scheme !
LOGICAL, PARAMETER :: ppm_flag = .false.

! Use the PPM (Colella 1984) reconstruction scheme !
LOGICAL, PARAMETER :: ppmc_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Realistic treatment of normal matter 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Flag for solving the internal energy equation !
LOGICAL, PARAMETER :: dual_energy = .false.   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for grid/atmospheric settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check if the NM/DM density fall below the lower limit !
LOGICAL, PARAMETER :: checkrho_flag = .false.

! Using custom variables floor !
LOGICAL, PARAMETER :: custom_floor = .false.

! Fix the DM/NM atmospheric density !
LOGICAL, PARAMETER :: fixrhonm_flag = .false.

! Check the computational box for DM/NM !
LOGICAL, PARAMETER :: checkstepnm_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Output setting
! This sets how frequent each type of profile is output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output configurations
! Notice they are not parameter because
! there are cases you want more frequent
! output at some dynamical scenarios
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for a global output
! This can be manipulated if one needs an output
! in all profiles outside the regular output time
LOGICAL :: output_file = .false.	

! Physical time interval for all log file
REAL*8 :: output_logtime = 1.0D1                            
REAL*8 :: output_logtime_last = 0.0D0

! Physical time interval for each hydro profile
REAL*8 :: output_profiletime = 0.05d0
REAL*8 :: output_profiletime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
