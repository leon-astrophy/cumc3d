!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE DEFINITION
USE MHD_MODULE
USE HDF5
IMPLICIT NONE

! for HDF5 !
character(len=99) :: globalt
character(len=99) :: filename
integer :: error, space_rank
integer(HSIZE_T) :: eps_dims(3), data_dims(4), dist_dims(1)
integer(HID_T) :: file_id, dspace_id, dset_id1

! integer !
INTEGER :: j  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for GPU !

#ifdef GPU
!$ACC UPDATE HOST(prim2(imin2:ibx-1,:,:,:), epsilon2, bcell)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! write to character !
write(globalt,'(I)') n_iter

! assign !
filename = './outfile/rkiter-'// trim(adjustl(globalt)) //'-nm.hdf5'

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"time",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,global_time,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibx - 1 - imin2) + 1
data_dims(2) = nx_2
data_dims(3) = ny_2 
data_dims(4) = nz_2 

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitive",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim2(imin2:ibx-1,1:nx_2,1:ny_2,1:nz_2),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
eps_dims(1) = nx_2
eps_dims(2) = ny_2 
eps_dims(3) = nz_2 

! open dataspace !
call h5screate_simple_f(space_rank,eps_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"epsilon",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,epsilon2(1:nx_2,1:ny_2,1:nz_2),eps_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibz - ibx) + 1
data_dims(2) = nx_2
data_dims(3) = ny_2 
data_dims(4) = nz_2 

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"bfield",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,bcell(ibx:ibz,1:nx_2,1:ny_2,1:nz_2),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine print the grid variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_grid
USE DEFINITION
USE HDF5
IMPLICIT NONE

! for HDF5 !
character(len=99) :: filename
integer :: error, space_rank
integer(HSIZE_T) :: dist_dims(1), vol_dims(3)
integer(HID_T) :: file_id, dspace_id, dset_id1

! integer !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign !
filename = './outfile/grid_param.hdf5'

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"dimension",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,DBLE(n_dim),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"coordinate",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,DBLE(coordinate_flag+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nx_2

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"x-direction",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,x2(1:nx_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nx_2 + 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"x-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,xF2(0:nx_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = ny_2

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"y-direction",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,y2(1:ny_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = ny_2 + 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"y-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,yF2(0:ny_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nz_2

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"z-direction",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,z2(1:nz_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nz_2 + 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"z-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,zF2(0:nz_2),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
vol_dims(1) = nx_2
vol_dims(2) = ny_2 
vol_dims(3) = nz_2 

! open dataspace !
call h5screate_simple_f(space_rank,vol_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"volume",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,vol2(1:nx_2,1:ny_2,1:nz_2),vol_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

END SUBROUTINE
