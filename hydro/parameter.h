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

! Double precision declarer, Define double precision
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)	

! Boundary flag notation for boundary conditions
! Range (0,1,2,3) stands for the scalar, X-type vector, Y-type vector and Z-type vector
INTEGER, PARAMETER :: even = 0, oddX = 1, oddY = 2, oddZ = 3	

! Which direction of reconstruction
INTEGER, PARAMETER :: x_dir = 1, y_dir = 2, z_dir = 3

! Apply boundary conditions to the full/partial box
INTEGER, PARAMETER :: part = 0, full = 1	

! Which type of matter considering
INTEGER, PARAMETER :: dm_f = 1, nm_f = 2

! pi constant 
REAL (DP), PARAMETER :: pi = 3.1415926535897932384626433832795E0_DP

! Planck constant
REAL (DP), PARAMETER :: hbar = 1.1965E-76_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Core part of the simulation box
! This is where the system variables are controled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! 0 = Cartesian coordinate
! 1 = Cylindrical coordinate
! 2 = Spherical coordinate
INTEGER, PARAMETER :: coordinate_flag = 0

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
! 1 = reflecting boundary (depend on scalar/vector)
! 2 = outgoing (1st derivative = 0)
INTEGER, PARAMETER :: boundary_flag(6) = (/0,0,0,0,0,0/)

! Flag for simulating with the full box (extend to negative x,y,z)
LOGICAL, PARAMETER :: fullx_flag = .false.
LOGICAL, PARAMETER :: fully_flag = .false.
LOGICAL, PARAMETER :: fullz_flag = .false.

! Want DM and NM to occupy different computational box?
LOGICAL, PARAMETER :: dmnm_diffbox = .false.

! The number of grid in the x,y,z direction for DM
INTEGER, PARAMETER :: nx_1 = 100
INTEGER, PARAMETER :: ny_1 = 100
INTEGER, PARAMETER :: nz_1 = 100

! The number of grid in the x,y,z direction for NM
INTEGER, PARAMETER :: nx_2 = 100
INTEGER, PARAMETER :: ny_2 = 100
INTEGER, PARAMETER :: nz_2 = 100

! Grid sizes for DM 
! Hint: The standard grid size (Note: 1 unit = 1.4774 km)
REAL (DP), PARAMETER :: dx1 = 0.01D0	
REAL (DP), PARAMETER :: dy1 = 0.01D0		
REAL (DP), PARAMETER :: dz1 = 0.01D0	

! Grid sizes for NM
REAL (DP), PARAMETER :: dx2 = 0.01D0	
REAL (DP), PARAMETER :: dy2 = 0.01D0	
REAL (DP), PARAMETER :: dz2 = 0.01D0	

! Section for fornax !
! Use non-constant grid? !
LOGICAL, PARAMETER :: fornax_flag = .false.

! For fornax grid
REAL (DP), PARAMETER :: A_fornax = 0.5d0
REAL (DP), PARAMETER :: xt = 1.5d2

! Cournat-Friedrich-Levy constant
! Definedas dt = cfl * dx / MAX(vel + cs)
REAL (DP), PARAMETER :: cfl = 0.20E0_DP			

! Maximum time to be simulated in the model
REAL (DP) :: total_time = 0.12D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for core hydrodynamic solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the LF Riemann solver !
LOGICAL, PARAMETER :: LF_flag = .true.

! Use the HLL Riemann solver !
LOGICAL, PARAMETER :: HLL_flag = .false.

! Use the WENO-Z (non-equidistant mesh) reconstruction scheme !
LOGICAL, PARAMETER :: weno_flag = .true.

! Use the PPM (Colella 1984) reconstruction scheme !
LOGICAL, PARAMETER :: ppm_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Realistic treatment of normal matter 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for including temperature effects !
LOGICAL, PARAMETER :: have_temp = .false.

! Flag for solving the internal energy equation !
LOGICAL, PARAMETER :: dual_energy = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Gravity solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for allowing gravity
LOGICAL, PARAMETER :: w_gravity = .false.

! Maximum number of multipole moments !
INTEGER, PARAMETER :: lmax = 16

! Maximum run time in relaxation of the potential !
INTEGER , PARAMETER :: relax_max = 300000    

! Tolerance = Maximum residue allowed for approval in relaxation !              
REAL (DP), PARAMETER :: tolerance = 3.0E-8_DP    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for The dark matter component 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Want to have a DM component ? !
LOGICAL, PARAMETER :: DM_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Temporaily
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! perform riemann problem hydro test? !
LOGICAL, PARAMETER :: hydro_test = .true.

! which test model ? !
INTEGER, PARAMETER :: test_model = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for grid/atmospheric settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check if the NM/DM density fall below the lower limit !
LOGICAL, PARAMETER :: checkrho_flag = .false.

! Fix the DM/NM atmospheric density !
LOGICAL, PARAMETER :: fixrhonm_flag = .false.
LOGICAL, PARAMETER :: fixrhodm_flag = .false.

! Check the computational box for DM/NM !
LOGICAL, PARAMETER :: checkstepnm_flag = .true.
LOGICAL, PARAMETER :: checkstepdm_flag = .false.

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
REAL (DP) :: output_logtime = 1.0D1                            
REAL (DP) :: output_logtime_last = 0.0D0

! Physical time interval for each hydro profile
REAL (DP) :: output_profiletime = 0.1D0
REAL (DP) :: output_profiletime_last = 0.0D0