!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE DEFINITION
USE HDF5
IMPLICIT NONE

! for HDF5 !
character(len=99) :: globalt
character(len=99) :: filename
integer :: error, space_rank
integer(HSIZE_T) :: data_dims(4), dist_dims(1)
integer(HID_T) :: file_id, dspace_id, dset_id1

! integer !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! write to character !
write(globalt,'(f0.10)') global_time

! assign !
filename = '../outfile/'// trim(adjustl(globalt)) //'.hdf5'

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
space_rank = 4
data_dims(1) = nx_2
data_dims(2) = ny_2 
data_dims(3) = nz_2 
data_dims(4) = (imax2 - imin2) + 1

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitive",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim2(1:nx_2,1:ny_2,1:nz_2,imin2:imax2),data_dims,error)

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