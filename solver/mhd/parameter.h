!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUMC3D-Ver1.65 (Last Modified: Leon)
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

! pi constant 
REAL*8, PARAMETER :: pi = 4.D0*DATAN(1.D0)

! define small number to avoid coordinate singularity !
REAL*8, PARAMETER :: small_num = TINY(1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Core part of the simulation box
! This is where the system variables are controled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! 0 = Cartesian coordinate
! 1 = Cylindrical coordinate
! 2 = Spherical coordinate
INTEGER, PARAMETER :: coordinate_flag = 1

! Dimension, 1 = 1D, 2 = 2D, 3 = 3D
INTEGER, PARAMETER :: n_dim = 3

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
INTEGER :: boundary_flag(6) = (/1,1,0,0,1,1/)

! Flag for simulating with the full box (extend to negative x,y,z)
LOGICAL, PARAMETER :: fullx_flag = .false.
LOGICAL, PARAMETER :: fully_flag = .false.
LOGICAL, PARAMETER :: fullz_flag = .false.

! Starting position of the grid !
REAL*8, PARAMETER :: x_start = 1.5d0
REAL*8, PARAMETER :: y_start = 0.0d0
REAL*8, PARAMETER :: z_start = -10.0d0

! Ending position of the grid !
REAL*8, PARAMETER :: x_end = 21.5d0
REAL*8, PARAMETER :: y_end = 2.0d0*pi
REAL*8, PARAMETER :: z_end = 10.0d0

! The number of grid in the x,y,z direction for NM
INTEGER, PARAMETER :: nx = 128 
INTEGER, PARAMETER :: ny = 128
INTEGER, PARAMETER :: nz = 128

! Grid sizes for NM
REAL*8, PARAMETER :: dx_ini = (x_end - x_start)/DBLE(nx)	
REAL*8, PARAMETER :: dy_ini = (y_end - y_start)/DBLE(ny)	
REAL*8, PARAMETER :: dz_ini = (z_end - z_start)/DBLE(nz)	

! Cournat-Friedrich-Levy constant
! Defined as dt = cfl * dx / MAX(vel + cs)
REAL*8, PARAMETER :: cfl = 0.80D0			

! Maximum time to be simulated in the model
REAL*8 :: total_time = 120.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for core hydrodynamic solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the LF (two-shocks) Riemann solver !
LOGICAL, PARAMETER :: LF_flag = .false.

! Use the HLL Riemann solver !
LOGICAL, PARAMETER :: HLL_flag = .false.

! Use the HLLC Riemann solver !
LOGICAL, PARAMETER :: HLLC_flag = .false.

! Use the HLLD Riemann solver !
LOGICAL, PARAMETER :: HLLD_flag = .true.

! Use the TVD (Mignone 2014) reconstruction scheme with Min-Mod limiter !
LOGICAL, PARAMETER :: tvdmm_flag = .false.

! Use the TVD (Mignone 2014) reconstruction scheme with MC limiter !
LOGICAL, PARAMETER :: tvdmc_flag = .false.

! Use the TVD (Mignone 2014) reconstruction scheme with Van-Leer limiter !
LOGICAL, PARAMETER :: tvdvl_flag = .false.

! Use the PPM (Colella 1984) reconstruction scheme !
LOGICAL, PARAMETER :: ppmc_flag = .true.

! Use the WENO (Shu 1997) reconstruction scheme !
LOGICAL, PARAMETER :: weno_flag = .false.

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
