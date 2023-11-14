!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE HDF5
USE DEFINITION
IMPLICIT NONE

! for HDF5 !
integer :: error, space_rank
character(len=99) :: globalt
character(len=99) :: filename
integer(HID_T) :: file_id, dspace_id, dset_id1
integer(HSIZE_T) :: eps_dims(3), data_dims(4), dist_dims(1)

! integer !
INTEGER :: j  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for GPU !

#ifdef GPU
!$ACC UPDATE HOST(prim(imin:imax,:,:,:), epsilon)
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
dist_dims(1) = nx + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"x-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,xF(-1:nx+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = ny + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"y-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,yF(-1:ny+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nz + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"z-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,zF(-1:nz+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibx - 1 - imin) + 1
data_dims(2) = nx + 2
data_dims(3) = ny + 2
data_dims(4) = nz + 2

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitive",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(imin:ibx-1,0:nx+1,0:ny+1,0:nz+1),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
eps_dims(1) = nx + 2
eps_dims(2) = ny + 2
eps_dims(3) = nz + 2

! open dataspace !
call h5screate_simple_f(space_rank,eps_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"epsilon",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,epsilon(0:nx+1,0:ny+1,0:nz+1),eps_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibz - ibx) + 1
data_dims(2) = nx + 3
data_dims(3) = ny + 3
data_dims(4) = nz + 3

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"bfield",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(ibx:ibz,-1:nx+1,-1:ny+1,-1:nz+1),data_dims,error)

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
