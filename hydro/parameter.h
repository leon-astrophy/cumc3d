!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUMC3D_Ver2022-1.0 (Last Modified: Leon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Foreword:
! Developed by Leung Shing Chi based on the WENO prototype developed by Wong Ka Wing in 2010
! c.f. Leung et al., MNRAS 454, 1238 (2015).
! Further developed by Leon H.S. Chan to extend to 3D, and include MHD/GRMHD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Constant and universal values 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Double precision declarer, Define double precision
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)	

! Boundary flag notation for boundary conditions
! Range (0,1,2,3) stands for the scalar, X-type vector, Y-type vector and Z-type vector
INTEGER, PARAMETER :: even = 0, oddX = 1, oddY = 2, oddZ = 3	

! Apply boundary conditions to the full/partial box
INTEGER, PARAMETER :: part = 0, full = 1	

! pi constant 
REAL (DP), PARAMETER :: pi = 3.1415926535897932384626433832795E0_DP

! Planck constant
REAL (DP), PARAMETER :: hbar = 1.1965E-76_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Core part of the simulation box
! This is where the system variables are controled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! 0 = Cartesian coordinate
! 1 = Cylindrical coordinate
! 2 = Spherical coordinate
INTEGER, PARAMETER :: coordinate_flag = 0

! Dimension, 1 = 1D, 2 = 2D, 3 = 3D
INTEGER, PARAMETER :: dimension_flag = 1

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
INTEGER, PARAMETER :: boundary_flag(6) = (/1,2,1,2,1,2/)

! Flag for simulating with the full box (extend to negative x,y,z)
LOGICAL, PARAMETER :: fullx_flag = .true.
LOGICAL, PARAMETER :: fully_flag = .true.
LOGICAL, PARAMETER :: fullz_flag = .true.

! The number of grid in the x,y,z direction for DM
INTEGER, PARAMETER :: nx_1 = 257
INTEGER, PARAMETER :: ny_1 = 129
INTEGER, PARAMETER :: nz_1 = 63

! The number of grid in the x,y,z direction for NM
INTEGER, PARAMETER :: nx_2 = 127
INTEGER, PARAMETER :: ny_2 = 63
INTEGER, PARAMETER :: nz_2 = 31

! Grid sizes for DM 
! Hint: The standard grid size (Note: 1 unit = 1.4774 km)
REAL (DP), PARAMETER :: dx1 = 0.001D0	
REAL (DP), PARAMETER :: dy1 = 0.01D0		
REAL (DP), PARAMETER :: dz1 = 0.1D0	

! Grid sizes for NM
REAL (DP), PARAMETER :: dx2 = 0.001D0	
REAL (DP), PARAMETER :: dy2 = 0.01D0		
REAL (DP), PARAMETER :: dz2 = 0.1D0	

! Section for fornax !
! Use non-constant grid? !
LOGICAL, PARAMETER :: fornax_flag = .true.

! For fornax grid
REAL (DP), PARAMETER :: A_fornax = 0.5d0
REAL (DP), PARAMETER :: xt = 1.5d2

! Cournat-Friedrich-Levy constant
! Definedas dt = cfl * dx / MAX(vel + cs)
REAL (DP), PARAMETER :: cfl = 0.20E0_DP			

! Maximum time to be simulated in the model
REAL (DP) :: total_time = 5.0D5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 3: Realistic treatment of normal matter 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for including temperature effects !
LOGICAL, PARAMETER :: have_temp = .true.

! Flag for solving the internal energy equation !
LOGICAL, PARAMETER :: dual_energy = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 4: Gravity solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for allowing gravity
LOGICAL, PARAMETER :: w_gravity = .true.

! Maximum number of multipole moments !
INTEGER, PARAMETER :: lmax = 16

! Maximum run time in relaxation of the potential !
INTEGER , PARAMETER :: relax_max = 300000    

! Tolerance = Maximum residue allowed for approval in relaxation !              
REAL (DP), PARAMETER :: tolerance = 3.0E-8_DP    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 5: The dark matter component 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Want to have a DM component ? !
LOGICAL, PARAMETER :: DM_flag = .true.

! Want the DM to be movable ? !
LOGICAL, PARAMETER :: runDM_flag = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 6: Temporaily
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER, PARAMETER :: test_model = 1

LOGICAL, PARAMETER :: fixrhonm_flag = .true.
LOGICAL, PARAMETER :: fixrhodm_flag = .true.
LOGICAL, PARAMETER :: checkrho_flag = .true.

LOGICAL, PARAMETER :: hydro_test = .true.

LOGICAL, PARAMETER :: etran_flag = .true.
LOGICAL, PARAMETER :: lapse_flag = .true.
LOGICAL, PARAMETER :: movinggridnm_flag = .true.
LOGICAL, PARAMETER :: movinggriddm_flag = .true.
LOGICAL, PARAMETER :: found_movinggridnm_flag = .true.
LOGICAL, PARAMETER :: found_movinggriddm_flag = .true.