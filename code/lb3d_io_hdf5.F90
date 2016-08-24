#include "lb3d.h"

!=========================================================================
!
! Copyright 1999-2012, Owners retain copyrights to their respective works.
!
! This file is part of lb3d.
!
! lb3d is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! lb3d is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
! License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with lb3d. If not, see <http://www.gnu.org/licenses/>.
!
!=========================================================================

!> hdf5 IO
module lb3d_io_hdf5_module
#ifdef USEHDF

  use HDF5
  use lb3d_mpi_module
  use lb3d_config_module, only: arg_input_dfile,arg_input_dfile_p,inp_file&
       &,lbeversion,write_AVS_fld
  use lb3d_io_helper_module, only: nflags,dump_avs_fld
  implicit none

  !Elena:  Define type to get all bit and pieces for metadata
  type hdf5_metadata
    character(len=80),   pointer :: inputfile(:)
    character(len=80),   pointer :: inputfile_diff(:)
#ifdef MD
    character(len=80),   pointer :: inputfile_md(:)
#endif
#ifdef ELEC
    character(len=80),   pointer :: inputfile_elec(:)
#endif
    character(len=80),   pointer :: mpidata(:)
    character(len=80),   pointer :: infodata(:)
  end type hdf5_metadata

  !FIXME not sure where to put the dimensionality, took it out, because it
  !is creating difficulties with the dependencies and is not implemented 
  !consistently anyway
 ! integer, parameter :: nd = 3

  save

  type (hdf5_metadata)           :: metadata
  character(len=80), pointer     :: fixed(:)
  integer :: n_meta_input, n_meta_mpi, n_meta_info, n_meta_total
  integer :: n_meta_input_diff = 0
#ifdef MD
  integer :: n_meta_input_md
#endif
#ifdef ELEC
  integer :: n_meta_input_elec
#endif
  logical :: made_metadata_hdf5 = .false.

#ifdef DEBUG_HDF5
  integer :: n_calls_read_iscalar_phdf5 = 0
  integer :: n_calls_read_rock_phdf5 = 0
  integer :: n_calls_dump_scalar_phdf5 = 0
  integer :: n_calls_dump_iscalar_phdf5 = 0
  integer :: n_calls_dump_vector_phdf5 = 0
  integer :: n_calls_total_phdf5 = 0
#endif

  character(len=32) :: hdfversion

contains

!>This routine calls \c h5open_f() to initialise HDF5. It should
!>only be called once during a given program run.
subroutine lb3d_io_init_hdf5()
  integer :: err
  integer :: majnum, minnum, relnum
  character(len=256) :: msgstr
  
  CALL log_msg("Initializing HDF5...",.false.)
  CALL h5open_f(err)
  CALL h5get_libversion_f(majnum, minnum, relnum, err)
  write(hdfversion,"(I0,'.',I0,'.',I0)") majnum, minnum, relnum
  write(msgstr,"('  Using HDF5 version ',I0,'.',I0,'.',I0,' .')") majnum, minnum, relnum
  CALL log_msg(trim(msgstr),.false.)

  ! Get all of the input file for the metadata later on, root only &
  ! only do it once (this is handled in lbe_make_metadata_hdf5()
  !N is local array of size N(0:nx+1,0:ny+1,0:nz+1) and each bit will be
  !kept on its PE
  if (myrankc==0) call lbe_make_metadata_hdf5

end subroutine lb3d_io_init_hdf5

!>This routine calls \c h5close_f() to initialise HDF5. It should
!>only be called once during a given program run.
subroutine lb3d_io_shutdown_hdf5()
  integer :: err
  CALL log_msg("Shutting down HDF5...",.false.)
  CALL h5close_f(err)
end subroutine lb3d_io_shutdown_hdf5

!> Reads an HDF5 file containig an integer scalar field in parallel.
!>
!> All CPUs open the file, but each CPU only loads its own subsection
!> of the file. The whole size of the file is supposed to be \c
!> tnx*tny*tnz data points.
!>
!> \param[in,out] iscalar The local chunk of the integer scalar field
!> is stored here, so the array ranges must be \c (1:nx,1:ny,1:nz) .
!>
!> \param[in] filename complete path and name of the HDF file to read from
subroutine read_iscalar_phdf5(iscalar,filename)
  integer,intent(inout) :: iscalar(1:,1:,1:)
  character(len=*),intent(in) :: filename
  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier
  integer(hid_t) :: mdsp_id      ! Dataspace identifier
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(3) :: stride = (/1,1,1/)    ! Stride
  integer(hsize_t), dimension(3) :: count = (/1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(3) :: blocksize             ! Size of a block
  integer(hsize_t), dimension(3) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(3) :: offset_m = (/0,0,0/)  ! Offset of the block in memory

  integer,dimension(nx,ny,nz) :: tmp_iscalar ! temporary storage
  integer(hsize_t),dimension(3) :: dims      ! dims of tmp_iscalar
  character(len=256) :: msgstr

  info = MPI_INFO_NULL ! MPI
#ifdef USE_IBM_LARGEBLOCK_IO
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

  ! Since we just want to take a block of data from the file and use
  ! it straightaway, these are all equal:
  dims = (/nx,ny,nz/)
  blocksize = (/nx,ny,nz/)

  ! Calculate which block to take
  offset = ccoords*dims

  DEBUG_HDF5_MSG_WS('------( HDF SCALAR READ DEBUG STARTED )------')

#ifdef DEBUG_HDF5
  n_calls_read_iscalar_phdf5 = n_calls_read_iscalar_phdf5+1
  n_calls_total_phdf5 = n_calls_total_phdf5 + 1
  CALL log_msg("HDF function call count:",.false.)
  write(msgstr,"('  read_iscalar_phdf5  : ',I0)") n_calls_read_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  read_rock_phdf5  : ',I0)") n_calls_read_rock_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_scalar_phdf5: ',I0)") n_calls_dump_scalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_iscalar_phdf5: ',I0)") n_calls_dump_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_vector_phdf5: ',I0)") n_calls_dump_vector_phdf5
  CALL log_msg(trim(msgstr),.false.)
#endif

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  DEBUG_HDF5_MSG('HDF opening file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,err,access_prp = plist_id)
  CALL h5pclose_f(plist_id, err)

  DEBUG_HDF5_MSG('HDF accessing dataset')
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  CALL h5dopen_f(file_id, dsetname, dset_id, err)
  CALL h5dget_space_f(dset_id, dsp_id, err)

  DEBUG_HDF5_MSG('HDF creating hyperslab 1')
  ! Select the data we want from the dataspace
  call h5sselect_hyperslab_f(dsp_id,H5S_SELECT_SET_F,offset,count,err&
       &,block=blocksize)

  ! Create a new dataspace to dump the data to - same shape as 
  ! iscalar, hence the offset is zero
  DEBUG_HDF5_MSG('HDF creating hyperslab 2')
  CALL h5screate_simple_f(ndim, dims, mdsp_id, err)
  call h5sselect_hyperslab_f(mdsp_id,H5S_SELECT_SET_F,offset_m,count,err&
       &,block=blocksize)

  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

  ! read the data into tmp_iscalar
  DEBUG_HDF5_MSG('HDF reading dataset')
  call H5dread_f(dset_id,H5T_NATIVE_INTEGER,tmp_iscalar,dims,err&
       &,file_space_id=dsp_id,mem_space_id=mdsp_id,xfer_prp=plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  DEBUG_HDF5_MSG('HDF closing')
  CALL h5pclose_f(plist_id, err)
  CALL h5sclose_f(dsp_id, err)
  CALL h5sclose_f(mdsp_id, err)
  CALL h5dclose_f(dset_id, err)
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')
  ! All done with HDF for now - all that remains is to move the data
  ! from the scalar array to the rock_state in a particular fashion.

  DEBUG_HDF5_MSG_WS('-----( HDF SCALAR READ DEBUG FINISHED )------')

  iscalar(1:nx,1:ny,1:nz) = tmp_iscalar
end subroutine read_iscalar_phdf5

!> Reads an HDF5 file into the rock matrix.
!> All CPUs open the file, but each CPU only loads its own
!> subsection of the rock.
subroutine read_rock_phdf5(filename,threshold,N)
  character(len=*) :: filename
  real(kind=rk) :: threshold
  type(lbe_site),dimension(0:,0:,0:) :: N

  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier
  integer(hid_t) :: mdsp_id      ! Dataspace identifier
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(3) :: stride = (/1,1,1/)    ! Stride
  integer(hsize_t), dimension(3) :: count = (/1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(3) :: blocksize             ! Size of a block
  integer(hsize_t), dimension(3) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(3) :: offset_m = (/0,0,0/)  ! Offset of the block in memory

  real(kind=rk), dimension(nx,ny,nz)    :: scalar                ! Temporary storage, real(kind=rk) is equivalent to double in C (?)
  integer(hsize_t), dimension(3) :: dims                  ! Dims of scalar
  real(kind=rk)                         :: rockfloat             ! Temporary variable to store a value of scalar in
  integer                        :: i,j,k                 ! Dummy loop variables

  character(len=256)     :: msgstr

  info = MPI_INFO_NULL ! MPI
#ifdef USE_IBM_LARGEBLOCK_IO
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/nx,ny,nz/)
  blocksize = (/nx,ny,nz/)

  ! Calculate which block to take
  offset(1) = ccoords(1)*dims(1)
  offset(2) = ccoords(2)*dims(2)
  offset(3) = ccoords(3)*dims(3)

  ! Specify that we'll be using MPIO to open the datafile and open it

  DEBUG_HDF5_MSG_WS('------( HDF SCALAR READ DEBUG STARTED )------')

#ifdef DEBUG_HDF5
  n_calls_read_rock_phdf5 = n_calls_read_rock_phdf5 + 1
  n_calls_total_phdf5 = n_calls_total_phdf5 + 1
  CALL log_msg("HDF function call count:",.false.)
  write(msgstr,"('  read_iscalar_phdf5  : ',I0)") n_calls_read_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  read_rock_phdf5  : ',I0)") n_calls_read_rock_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_scalar_phdf5: ',I0)") n_calls_dump_scalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_iscalar_phdf5: ',I0)") n_calls_dump_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_vector_phdf5: ',I0)") n_calls_dump_vector_phdf5
  CALL log_msg(trim(msgstr),.false.)
!   ! Magic reset button
!   if ( n_calls_total_phdf5 > 1000 ) then
!     CALL log_msg("HDF soft reset test...",.false.)
!     CALL lb3d_io_shutdown_hdf5()
!     CALL lb3d_io_init_hdf5()
!     n_calls_total_phdf5 = 0
!   endif
#endif

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  DEBUG_HDF5_MSG('HDF opening file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err, access_prp = plist_id)
  CALL h5pclose_f(plist_id, err)

  DEBUG_HDF5_MSG('HDF accessing dataset')
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  CALL h5dopen_f(file_id, dsetname, dset_id, err)
  CALL h5dget_space_f(dset_id, dsp_id, err)

  DEBUG_HDF5_MSG('HDF creating hyperslab 1')
  ! Select the data we want from the dataspace
  CALL h5sselect_hyperslab_f(dsp_id, H5S_SELECT_SET_F, offset, count, err, block = blocksize)

  ! Create a new dataspace to dump the data to - same shape as the 'scalar' array, hence the offset is zero
  DEBUG_HDF5_MSG('HDF creating hyperslab 2')
  CALL h5screate_simple_f(ndim, dims, mdsp_id, err)
  CALL h5sselect_hyperslab_f(mdsp_id, H5S_SELECT_SET_F, offset_m, count, err, block = blocksize)

  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

  ! Finally read the data into 'scalar', yay!
  DEBUG_HDF5_MSG('HDF reading dataset')
  CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, dims, err, file_space_id = dsp_id, mem_space_id = mdsp_id, xfer_prp = plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  DEBUG_HDF5_MSG('HDF closing')
  CALL h5pclose_f(plist_id, err)
  CALL h5sclose_f(dsp_id, err)
  CALL h5sclose_f(mdsp_id, err)
  CALL h5dclose_f(dset_id, err)
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')
  ! All done with HDF for now - all that remains is to move the data from the scalar array to the rock_state in a particular fashion.

  DEBUG_HDF5_MSG_WS('-----( HDF SCALAR READ DEBUG FINISHED )------')

#ifdef XDRROCKWET
  CALL log_msg("  Generating geometry using XDRROCKWET",.false.)
#else
  write(msgstr,"('  Generating geometry without using XDRROCKWET, threshold = ',F16.10)") threshold
  CALL log_msg(trim(msgstr),.false.)
#endif

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rockfloat = scalar(i,j,k)
        ! Code below copied from read_rock_xdrf_par function

#ifdef XDRROCKWET
        ! Set rock state to 1 if value in xdr file 
        ! is larger than 0.
        ! Colour the rock according to value in xdr file.
        ! To be able to distinguish between rock/non-rock
        ! subtract an offset of 5.

        if (rockfloat > 0.0) then
          rockfloat = rockfloat - 5.
          N(i,j,k)%rock_state = 1.d0
          N(i,j,k)%n_r(:nnonrest)=0.d0
#ifndef SINGLEFLUID
          N(i,j,k)%n_b(:nnonrest)=0.d0
#endif
#ifndef NOSURFACTANT
          N(i,j,k)%n_s(:nnonrest)=0.d0
          N(i,j,k)%n_s(restvec)=0.d0
#endif
          if (rockfloat > 0) then
            N(i,j,k)%n_r(restvec)=dble(rockfloat)
#ifdef WALLCONST
            N(i,j,k)%n_r(restvec)=dble(rockfloat)/tau_r
#endif
#ifndef SINGLEFLUID
            N(i,j,k)%n_b(restvec)=0.d0
#endif
          elseif (rockfloat < 0) then
#ifndef SINGLEFLUID
            N(i,j,k)%n_b(restvec)=-dble(rockfloat)
#ifdef WALLCONST
            N(i,j,k)%n_b(restvec)=-dble(rockfloat)/tau_b
#endif
#endif
            N(i,j,k)%n_r(restvec)=0.d0
          else
            N(i,j,k)%n_r(restvec)=0.d0
#ifndef SINGLEFLUID
            N(i,j,k)%n_b(restvec)=0.d0
#endif
          endif
        else 
          N(i,j,k)%rock_state = 0.d0
        endif
#else
        ! Set rock state to 1 if value in xdr file 
        ! is higher than threshold
        if (rockfloat > threshold) then
          N(i,j,k)%rock_state = 1.d0
        else
          N(i,j,k)%rock_state = 0.d0
        endif
#endif
      enddo
    enddo
  enddo
end subroutine read_rock_phdf5

!> Reads an HDF5 grid vector file.
!> All CPUs open the files, but each CPU only loads its own
!> subsection of the grid.
subroutine read_vector_phdf5(vector,filename,vdim)
!  implicit none
  real(kind=rk), dimension(:,:,:,:) :: vector
  integer :: vdim
  character(len=*) :: filename

  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 4 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier in file
  integer(hid_t) :: mdsp_id      ! Dataspace identifier in memory
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(4) :: stride = (/1,1,1,1/)    ! Stride
  integer(hsize_t), dimension(4) :: count = (/1,1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(4) :: blocksize               ! Size of a block
  integer(hsize_t), dimension(4) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(4) :: offset_m = (/0,0,0,0/)  ! Offset of the block in memory

!  real(kind=rk), dimension(nx,ny,nz)    :: scalar                ! Temporary storage, real(kind=rk) is equivalent to double in C (?)
  integer(hsize_t), dimension(4) :: dims                  ! Dims of scalar
!  real(kind=rk)                         :: rockfloat             ! Temporary variable to store a value of scalar in
  integer                        :: i,j,k                 ! Dummy loop variables

  character(len=256)     :: msgstr

  info = MPI_INFO_NULL ! MPI
#ifdef USE_IBM_LARGEBLOCK_IO
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

#ifndef HDF5_FLIP
  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/vdim,nx,ny,nz/)
  blocksize = (/vdim,nx,ny,nz/)

  ! Calculate which block to take
  offset(1) = 0
  offset(2) = ccoords(1)*dims(2)
  offset(3) = ccoords(2)*dims(3)
  offset(4) = ccoords(3)*dims(4)
#else
  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/nx,ny,nz,vdim/)
  blocksize = (/nx,ny,nz,vdim/)

  ! Calculate which block to take
  offset(1) = ccoords(1)*dims(1)
  offset(2) = ccoords(2)*dims(2)
  offset(3) = ccoords(3)*dims(3)
  offset(4) = 0
#endif

  ! Specify that we'll be using MPIO to open the datafile and open it

  DEBUG_HDF5_MSG_WS('------( HDF VECTOR READ DEBUG STARTED )------')

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  DEBUG_HDF5_MSG('HDF opening file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, err, access_prp = plist_id)
  CALL h5pclose_f(plist_id, err)

  DEBUG_HDF5_MSG('HDF accessing dataset')
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  CALL h5dopen_f(file_id, dsetname, dset_id, err)
  CALL h5dget_space_f(dset_id, dsp_id, err)

  DEBUG_HDF5_MSG('HDF creating hyperslab 1')
  ! Select the data we want from the dataspace
  CALL h5sselect_hyperslab_f(dsp_id, H5S_SELECT_SET_F, offset, count, err, block = blocksize)

  ! Create a new dataspace to dump the data to - same shape as the
  ! 'vector' array, hence the offset is zero
  DEBUG_HDF5_MSG('HDF creating hyperslab 2')
  CALL h5screate_simple_f(ndim, dims, mdsp_id, err)
  CALL h5sselect_hyperslab_f(mdsp_id, H5S_SELECT_SET_F, offset_m, count, err, block = blocksize)

  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

  ! Finally read the data into 'vector', yay!
  DEBUG_HDF5_MSG('HDF reading dataset')
  CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, vector, dims, err, file_space_id = dsp_id, mem_space_id = mdsp_id, xfer_prp = plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  DEBUG_HDF5_MSG('HDF closing')
  CALL h5pclose_f(plist_id, err)
  CALL h5sclose_f(dsp_id, err)
  CALL h5sclose_f(mdsp_id, err)
  CALL h5dclose_f(dset_id, err)
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')

  DEBUG_HDF5_MSG_WS('-----( HDF VECTOR READ DEBUG FINISHED )------')
end subroutine read_vector_phdf5

!>  Called by subroutine dump_scalar.
!>  For HDF5 formatted IO. Only available for postprocessing (post=.true.)
subroutine dump_scalar_phdf5(scalar,name)
  real(kind=rk), dimension(1:, 1:, 1:  ) :: scalar

  character(len=*)   :: name
  character(len=256) :: filename
  character(len=256) :: attrname

  character(len=8), parameter :: dsetname = 'OutArray' !Dataset name

  character(len=256) :: msgstr

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id       ! Datatype id for array (real or double)

  integer(hsize_t),  dimension(3) :: dimsf
  integer(hsize_t),  dimension(3) :: dims

  integer(hsize_t),  dimension(3) :: count
  integer(hssize_t), dimension(3) :: offset
  ! integer(hsize_t),  dimension(3) :: chunk_dims

  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

#ifdef DEBUG_HDF5
  ! For debugging purposes, we have a lot of variables storing time, some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s, t_f, t_init_s, t_init_f, t_write_s, t_write_f, t_close_s, t_close_f, t_closef_s, t_closef_f, t_wait, t_init, t_write, t_close, t_closef, t_total_nowait, t_total, p_wait
  integer :: tag = 1
  integer :: i
  integer status(MPI_STATUS_SIZE) 

  CALL MPI_Barrier(Comm_cart,err)
  t_s = MPI_Wtime()
#endif

  CALL lbe_make_filename_output(filename, trim(name), '.h5', nt)
  attrname = filename ! attrname has to be unique, can't have it twice

  DEBUG_HDF5_MSG_WS('--------( HDF SCALAR DEBUG STARTED )---------')

#ifdef DEBUG_HDF5
  n_calls_dump_scalar_phdf5 = n_calls_dump_scalar_phdf5 + 1
  n_calls_total_phdf5 = n_calls_total_phdf5 + 1
  CALL log_msg("HDF function call count:",.false.)
  write(msgstr,"('  read_iscalar_phdf5  : ',I0)") n_calls_read_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  read_rock_phdf5  : ',I0)") n_calls_read_rock_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_scalar_phdf5: ',I0)") n_calls_dump_scalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_iscalar_phdf5: ',I0)") n_calls_dump_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_vector_phdf5: ',I0)") n_calls_dump_vector_phdf5
  CALL log_msg(trim(msgstr),.false.)
!   ! Magic reset button
!   if ( n_calls_total_phdf5 > 1000 ) then
!     CALL log_msg("HDF soft reset test...",.false.)
!     CALL lb3d_io_shutdown_hdf5()
!     CALL lb3d_io_init_hdf5()
!     n_calls_total_phdf5 = 0
!   endif
#endif

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  info = MPI_INFO_NULL
#ifdef USE_IBM_LARGEBLOCK_IO
  DEBUG_HDF5_MSG('HDF using IBM_largeblock_io')
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

  dimsf  = (/tnx, tny, tnz/)
  dims   = (/nx, ny, nz/)
  count  = (/nx, ny, nz/)
  ! chunk_dims = (/nx, ny, nz/)

  ! Note: the above is in fact equivalent to the count/stride/block formulation:
  ! count  = (/1, 1, 1/)
  ! stride = (/1, 1, 1/)
  ! block  = (/nx, ny, nz/)
  ! (with matching parameters in h5sselect_hyperslab_f

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz

  if(dump_double)then
    DEBUG_HDF5_MSG('HDF datatype is H5T_NATIVE_DOUBLE')
    type_id = H5T_NATIVE_DOUBLE
  else
    DEBUG_HDF5_MSG('HDF datatype is H5T_NATIVE_REAL')
    type_id = H5T_NATIVE_REAL
  endif

#ifdef DEBUG_HDF5
  CALL MPI_Barrier(Comm_cart,err)
  t_init_s = MPI_Wtime()
#endif

  ! Setup file access property list with parallel I/O access.
  DEBUG_HDF5_MSG('HDF creating file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_cart, info, err)
  ! Create the file collectively.
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp = plist_id)
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  DEBUG_HDF5_MSG('HDF creating filespace')
  CALL h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! DEBUG_HDF5_MSG('HDF creating chunked dataset')
  ! CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! CALL h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! CALL h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  DEBUG_HDF5_MSG('HDF creating continuous dataset')
  CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  CALL h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
  DEBUG_HDF5_MSG('HDF creating memspace')
  CALL h5screate_simple_f(ndim, dims, memspace, err)

  ! Select hyperslab in the file.
  DEBUG_HDF5_MSG('HDF selecting hyperslab')
  CALL h5dget_space_f(dset_id, filespace, err)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

#ifdef DEBUG_HDF5
  t_init_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_write_s = MPI_Wtime()
#endif

  ! Different write calls for double or single data (have to convert real(kind=rk) scalar to real*4)
  if (dump_double) then
    DEBUG_HDF5_MSG('HDF writing double data')
    CALL h5dwrite_f(dset_id, type_id, scalar, dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  else
    DEBUG_HDF5_MSG('HDF writing single data')
    CALL h5dwrite_f(dset_id, type_id, real(scalar,4), dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  endif

#ifdef DEBUG_HDF5
  t_write_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_close_s = MPI_Wtime()
#endif

  ! Close dataspaces.
  DEBUG_HDF5_MSG('HDF closing filespace handle')
  CALL h5sclose_f(filespace, err)
  DEBUG_HDF5_MSG('HDF closing memspace handle')
  CALL h5sclose_f(memspace, err)
  ! Close the dataset.
  DEBUG_HDF5_MSG('HDF closing dataset handle')
  CALL h5dclose_f(dset_id, err)
  ! Close the property list.
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

#ifdef DEBUG_HDF5
  t_close_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_closef_s = MPI_Wtime()
#endif

  ! Close the file.
  DEBUG_HDF5_MSG('HDF closing file handle')
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')

#ifdef DEBUG_HDF5
  t_closef_f = MPI_Wtime()
#endif

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if (myrankc .eq. 0 ) then
    DEBUG_HDF5_MSG('HDF writing metadata')
    CALL lbe_write_attr_phdf5(filename, dsetname, attrname)
    DEBUG_HDF5_MSG('HDF finished writing metadata')
  endif

  DEBUG_HDF5_MSG_WS('-------( HDF SCALAR DEBUG FINISHED )---------')

#ifdef DEBUG_HDF5_TIMING
  ! This is a lot of debugging/timing stuff
  ! All processors create a string with their timings, then send it to rank 0.
  ! Rank 0 can then display them all in correct order
  CALL MPI_Barrier(Comm_cart,err)
  t_f = MPI_Wtime()

  DEBUG_HDF5_MSG('------------( HDF DEBUG TIMERS )-------------')
  DEBUG_HDF5_MSG('  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT')

  t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
  t_init = t_init_f - t_init_s
  t_write = t_write_f - t_write_s
  t_close = t_close_f - t_close_s
  t_closef = t_closef_f - t_closef_s

  t_total_nowait = t_init + t_write + t_close + t_closef
  t_total = t_f - t_s
  p_wait = t_wait / t_total


  ! Magic number '256' corresponds to the length of msgstr of course
  if ( myrankc .gt. 0 ) then
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    CALL MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
  else
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    DEBUG_HDF5_MSG(trim(msgstr))
    do i=1, nprocs-1
      CALL MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
      DEBUG_HDF5_MSG(trim(msgstr))
    enddo
  endif

  DEBUG_HDF5_MSG_WS('----------( END HDF DEBUG TIMERS )-----------')
#endif

end subroutine dump_scalar_phdf5

!> Called by subroutine dump_iscalar. For HDF5 formatted IO. Dumps
!> integer scalar fields.
subroutine dump_iscalar_phdf5(iscalar,name)
  integer,intent(in),dimension(1:,1:,1:) :: iscalar
  character(len=*),intent(in) :: name

  character(len=256) :: filename
  character(len=256) :: attrname

  character(len=8), parameter :: dsetname = 'OutArray' !Dataset name

  character(len=256) :: msgstr

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id       ! Datatype id for array (real or double)

  integer(hsize_t),  dimension(3) :: dimsf
  integer(hsize_t),  dimension(3) :: dims

  integer(hsize_t),  dimension(3) :: count
  integer(hssize_t), dimension(3) :: offset
  ! integer(hsize_t),  dimension(3) :: chunk_dims

  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

#ifdef DEBUG_HDF5
  ! For debugging purposes, we have a lot of variables storing time,
  ! some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s,t_f,t_init_s,t_init_f,t_write_s,t_write_f,t_close_s&
       &,t_close_f,t_closef_s,t_closef_f,t_wait,t_init,t_write,t_close,t_closef&
       &,t_total_nowait,t_total,p_wait
  integer,parameter :: tag=1
  integer :: i
  integer status(MPI_STATUS_SIZE)

  CALL MPI_Barrier(Comm_cart,err)
  t_s = MPI_Wtime()
#endif

  CALL lbe_make_filename_output(filename, trim(name), '.h5', nt)
  attrname = filename ! attrname has to be unique, can't have it twice

  DEBUG_HDF5_MSG_WS('--------( HDF ISCALAR DEBUG STARTED )---------')

#ifdef DEBUG_HDF5
  n_calls_dump_iscalar_phdf5 = n_calls_dump_iscalar_phdf5 + 1
  n_calls_total_phdf5 = n_calls_total_phdf5 + 1
  CALL log_msg("HDF function call count:",.false.)
  write(msgstr,"('  read_iscalar_phdf5  : ',I0)") n_calls_read_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  read_rock_phdf5  : ',I0)") n_calls_read_rock_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_scalar_phdf5: ',I0)") n_calls_dump_scalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_iscalar_phdf5: ',I0)") n_calls_dump_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_vector_phdf5: ',I0)") n_calls_dump_vector_phdf5
  CALL log_msg(trim(msgstr),.false.)
!   ! Magic reset button
!   if ( n_calls_total_phdf5 > 1000 ) then
!     CALL log_msg("HDF soft reset test...",.false.)
!     CALL lb3d_io_shutdown_hdf5()
!     CALL lb3d_io_init_hdf5()
!     n_calls_total_phdf5 = 0
!   endif
#endif

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  info = MPI_INFO_NULL
#ifdef USE_IBM_LARGEBLOCK_IO
  DEBUG_HDF5_MSG('HDF using IBM_largeblock_io')
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

  dimsf  = (/tnx, tny, tnz/)
  dims   = (/nx, ny, nz/)
  count  = (/nx, ny, nz/)
  ! chunk_dims = (/nx, ny, nz/)

  ! Note: the above is in fact equivalent to the count/stride/block formulation:
  ! count  = (/1, 1, 1/)
  ! stride = (/1, 1, 1/)
  ! block  = (/nx, ny, nz/)
  ! (with matching parameters in h5sselect_hyperslab_f

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz

  DEBUG_HDF5_MSG('HDF datatype is H5T_NATIVE_INTEGER')
  type_id = H5T_NATIVE_INTEGER

#ifdef DEBUG_HDF5
  CALL MPI_Barrier(Comm_cart,err)
  t_init_s = MPI_Wtime()
#endif

  ! Setup file access property list with parallel I/O access.
  DEBUG_HDF5_MSG('HDF creating file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_cart, info, err)
  ! Create the file collectively.
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,err,access_prp=plist_id)
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  DEBUG_HDF5_MSG('HDF creating filespace')
  CALL h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! DEBUG_HDF5_MSG('HDF creating chunked dataset')
  ! CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! CALL h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! CALL h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  DEBUG_HDF5_MSG('HDF creating continuous dataset')
  CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  CALL h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
  DEBUG_HDF5_MSG('HDF creating memspace')
  CALL h5screate_simple_f(ndim, dims, memspace, err)

  ! Select hyperslab in the file.
  DEBUG_HDF5_MSG('HDF selecting hyperslab')
  CALL h5dget_space_f(dset_id, filespace, err)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

#ifdef DEBUG_HDF5
  t_init_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_write_s = MPI_Wtime()
#endif

  DEBUG_HDF5_MSG('HDF writing integer data')
  call h5dwrite_f(dset_id,type_id,iscalar,dims,err,file_space_id=filespace&
       &,mem_space_id=memspace,xfer_prp=plist_id)

#ifdef DEBUG_HDF5
  t_write_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_close_s = MPI_Wtime()
#endif

  ! Close dataspaces.
  DEBUG_HDF5_MSG('HDF closing filespace handle')
  CALL h5sclose_f(filespace, err)
  DEBUG_HDF5_MSG('HDF closing memspace handle')
  CALL h5sclose_f(memspace, err)
  ! Close the dataset.
  DEBUG_HDF5_MSG('HDF closing dataset handle')
  CALL h5dclose_f(dset_id, err)
  ! Close the property list.
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

#ifdef DEBUG_HDF5
  t_close_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_closef_s = MPI_Wtime()
#endif

  ! Close the file.
  DEBUG_HDF5_MSG('HDF closing file handle')
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')

#ifdef DEBUG_HDF5
  t_closef_f = MPI_Wtime()
#endif

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if (myrankc .eq. 0 ) then
    DEBUG_HDF5_MSG('HDF writing metadata')
    CALL lbe_write_attr_phdf5(filename, dsetname, attrname)
    DEBUG_HDF5_MSG('HDF finished writing metadata')
  endif

  DEBUG_HDF5_MSG_WS('-------( HDF SCALAR DEBUG FINISHED )---------')

#ifdef DEBUG_HDF5_TIMING
  ! This is a lot of debugging/timing stuff
  ! All processors create a string with their timings, then send it to rank 0.
  ! Rank 0 can then display them all in correct order
  CALL MPI_Barrier(Comm_cart,err)
  t_f = MPI_Wtime()

  DEBUG_HDF5_MSG('------------( HDF DEBUG TIMERS )-------------')
  DEBUG_HDF5_MSG('  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT')

  t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
  t_init = t_init_f - t_init_s
  t_write = t_write_f - t_write_s
  t_close = t_close_f - t_close_s
  t_closef = t_closef_f - t_closef_s

  t_total_nowait = t_init + t_write + t_close + t_closef
  t_total = t_f - t_s
  p_wait = t_wait / t_total


  ! Magic number '256' corresponds to the length of msgstr of course
  if ( myrankc .gt. 0 ) then
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    CALL MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
  else
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    DEBUG_HDF5_MSG(trim(msgstr))
    do i=1, nprocs-1
      CALL MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
      DEBUG_HDF5_MSG(trim(msgstr))
    enddo
  endif

  DEBUG_HDF5_MSG_WS('----------( END HDF DEBUG TIMERS )-----------')
#endif
end subroutine dump_iscalar_phdf5

!> Postprocessed IO done with PHDF5.
!> Dump an array of vectors in hdf5 format.
subroutine dump_vector_phdf5(vector,filename)
  real(kind=rk), intent(inout),dimension(1:,1:,1:,1:) :: vector
  character(len=*) :: filename
  character(len=256) :: attrname
  integer :: vecdim
  integer :: nxi,nyi,nzi
  !HDF5 stuff
  character(len=8), parameter :: dsetname='OutArray' !Dataset name

  character(len=256) :: msgstr

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id

  integer(hsize_t), dimension(4) :: dims
  integer(hsize_t), dimension(4) :: dimsf

  integer(hsize_t),  dimension(4) :: count
  integer(hssize_t), dimension(4) :: offset
  ! integer(hsize_t),  dimension(4) :: chunk_dims

  integer , parameter :: ndim = 4 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

#ifdef DEBUG_HDF5
  ! For debugging purposes, we have a lot of variables storing time, some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s, t_f, t_init_s, t_init_f, t_write_s, t_write_f, t_close_s, t_close_f, t_closef_s, t_closef_f, t_wait, t_init, t_write, t_close, t_closef, t_total_nowait, t_total, p_wait
  integer :: tag = 1
  integer :: i
  integer status(MPI_STATUS_SIZE) 

  CALL MPI_Barrier(Comm_cart,err)
  t_s = MPI_Wtime()
#endif

#ifndef HDF5_FLIP
  vecdim = size(vector,1)
  nxi    = size(vector,2)
  nyi    = size(vector,3)
  nzi    = size(vector,4)
#else
  nxi    = size(vector,1)
  nyi    = size(vector,2)
  nzi    = size(vector,3)
  vecdim = size(vector,4)
#endif

  attrname = filename

  DEBUG_HDF5_MSG_WS('--------( HDF VECTOR DEBUG STARTED )---------')

#ifdef DEBUG_HDF5
  n_calls_dump_vector_phdf5 = n_calls_dump_vector_phdf5 + 1
  n_calls_total_phdf5 = n_calls_total_phdf5 + 1
  CALL log_msg("HDF function call count:",.false.)
  write(msgstr,"('  read_iscalar_phdf5  : ',I0)") n_calls_read_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  read_rock_phdf5  : ',I0)") n_calls_read_rock_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_scalar_phdf5: ',I0)") n_calls_dump_scalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_iscalar_phdf5: ',I0)") n_calls_dump_iscalar_phdf5
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('  dump_vector_phdf5: ',I0)") n_calls_dump_vector_phdf5
  CALL log_msg(trim(msgstr),.false.)
!   ! Magic reset button
!   if ( n_calls_total_phdf5 > 1000 ) then
!     CALL log_msg("HDF soft reset test...",.false.)
!     CALL lb3d_io_shutdown_hdf5()
!     CALL lb3d_io_init_hdf5()
!     n_calls_total_phdf5 = 0
!   endif
#endif

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))

  write(msgstr,"('nxi = ',I0,', nyi = ',I0,', nzi = ',I0,', vecdim = ',I0)") nxi, nyi, nzi, vecdim
  DEBUG_HDF5_MSG(trim(msgstr))

  info = MPI_INFO_NULL
#ifdef USE_IBM_LARGEBLOCK_IO
  DEBUG_HDF5_MSG('HDF using IBM_largeblock_io')
  CALL MPI_Info_create(info, err)
  CALL MPI_Info_set(info, "IBM_largeblock_io", "true", err)
#endif

#ifndef HDF5_FLIP
  dimsf  = (/vecdim, tnx, tny, tnz/)
  dims   = (/vecdim, nx, ny, nz/)
  count  = (/vecdim, nx, ny, nz/)
  ! chunk_dims = (/vecdim, nx, ny, nz/)

  offset(1) = 0
  offset(2) = ccoords(1)*nx
  offset(3) = ccoords(2)*ny
  offset(4) = ccoords(3)*nz
#else
  dimsf  = (/tnx, tny, tnz, vecdim/)
  dims   = (/nx, ny, nz, vecdim/)
  count  = (/nx, ny, nz, vecdim/)
  ! chunk_dims = (/nx, ny, nz, vecdim/)

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz
  offset(4) = 0

#endif


  if(dump_double)then
    DEBUG_HDF5_MSG('HDF datatype is H5T_NATIVE_DOUBLE')
    type_id = H5T_NATIVE_DOUBLE
  else
    DEBUG_HDF5_MSG('HDF datatype is H5T_NATIVE_REAL')
    type_id = H5T_NATIVE_REAL
  endif

#ifdef DEBUG_HDF5
  CALL MPI_Barrier(Comm_cart,err)
  t_init_s = MPI_Wtime()
#endif

  ! Setup file access property list with parallel I/O access.
  DEBUG_HDF5_MSG('HDF creating file')
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  CALL h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  ! Create the file collectively.
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp = plist_id)
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  DEBUG_HDF5_MSG('HDF creating filespace')
  CALL h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! DEBUG_HDF5_MSG('HDF creating chunked dataset')
  ! CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! CALL h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! CALL h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  DEBUG_HDF5_MSG('HDF creating continuous dataset')
  CALL h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  CALL h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
  DEBUG_HDF5_MSG('HDF creating memspace')
  CALL h5screate_simple_f(ndim, dims, memspace, err)  

  ! Select hyperslab in the file.
  DEBUG_HDF5_MSG('HDF selecting hyperslab')
  CALL h5dget_space_f(dset_id, filespace, err)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
#ifdef USE_HDF5_INDEPENDENT_IO
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_INDEPENDENT_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
#else
  DEBUG_HDF5_MSG('HDF using H5FD_MPIO_COLLECTIVE_F')
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
#endif

#ifdef DEBUG_HDF5
  t_init_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_write_s = MPI_Wtime()
#endif

  ! Different write calls for double or single data (have to convert real(kind=rk) scalar to real*4
  if (dump_double) then
    DEBUG_HDF5_MSG('HDF writing double data')
    CALL h5dwrite_f(dset_id, type_id, vector, dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  else
    DEBUG_HDF5_MSG('HDF writing single data')
    CALL h5dwrite_f(dset_id, type_id, real(vector,4), dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  endif

#ifdef DEBUG_HDF5
  t_write_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_close_s = MPI_Wtime()
#endif

  ! Close dataspaces.
  DEBUG_HDF5_MSG('HDF closing filespace handle')
  CALL h5sclose_f(filespace, err)
  DEBUG_HDF5_MSG('HDF closing memspace handle')
  CALL h5sclose_f(memspace, err)
  ! Close the dataset.
  DEBUG_HDF5_MSG('HDF closing dataset handle')
  CALL h5dclose_f(dset_id, err)
  ! Close the property list.
  DEBUG_HDF5_MSG('HDF closing property list handle')
  CALL h5pclose_f(plist_id, err)

#ifdef DEBUG_HDF5
  t_close_f = MPI_Wtime()
  CALL MPI_Barrier(Comm_cart,err)
  t_closef_s = MPI_Wtime()
#endif

  ! Close the file.
  DEBUG_HDF5_MSG('HDF closing file handle')
  CALL h5fclose_f(file_id, err)
  DEBUG_HDF5_MSG('HDF finished closing handles')

#ifdef DEBUG_HDF5
  t_closef_f = MPI_Wtime()
#endif

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if(myrankc.eq.0)then
    DEBUG_HDF5_MSG('HDF writing metadata')
    CALL lbe_write_attr_phdf5(filename, dsetname, attrname)
    DEBUG_HDF5_MSG('HDF finished writing metadata')
  endif

  DEBUG_HDF5_MSG_WS('-------( HDF VECTOR DEBUG FINISHED )---------')

#ifdef DEBUG_HDF5_TIMING
  ! This is a lot of debugging/timing stuff
  ! All processors create a string with their timings, then send it to rank 0.
  ! Rank 0 can then display them all in correct order
  CALL MPI_Barrier(Comm_cart,err)
  t_f = MPI_Wtime()

  DEBUG_HDF5_MSG('------------( HDF DEBUG TIMERS )-------------')
  DEBUG_HDF5_MSG('  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT')

  t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
  t_init = t_init_f - t_init_s
  t_write = t_write_f - t_write_s
  t_close = t_close_f - t_close_s
  t_closef = t_closef_f - t_closef_s

  t_total_nowait = t_init + t_write + t_close + t_closef
  t_total = t_f - t_s
  p_wait = t_wait / t_total


  ! Magic number '256' corresponds to the length of msgstr of course
  if ( myrankc .gt. 0 ) then
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    CALL MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
  else
    write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
    DEBUG_HDF5_MSG(trim(msgstr))
    do i=1, nprocs-1
      CALL MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
      DEBUG_HDF5_MSG(trim(msgstr))
    enddo
  endif

  DEBUG_HDF5_MSG_WS('----------( END HDF DEBUG TIMERS )-----------')
#endif

  if (write_AVS_fld) CALL dump_avs_fld(trim(filename),nxi,nyi,nzi,int(4))
end subroutine dump_vector_phdf5

!>This subroutine adds the checkpointing parameters to a file.
!>
!>It is called from root (rank 0) of the checkpointing subroutines.
subroutine lbe_write_cp_param_phdf5(data,filename,aname)
  implicit none

  character(len=*)   :: filename  ! File name
  character(len=*)   :: aname     ! Attribute / parameter name

  character(len=8), parameter :: dsetname='OutArray' ! Dataset name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer            :: error ! Error flag

  ! HDF wants an array, so we make a 1d array of length 1
  real(kind=rk) :: data
  real(kind=rk), dimension(1) :: data_array
  character(len=256) :: msgstr

  adims(1) = 1
  data_array(1) = data


  DEBUG_HDF5_MSG_WS('--------( HDF PARAMS DEBUG STARTED )---------')

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))


  DEBUG_HDF5_MSG('HDF opening file')
  CALL h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  DEBUG_HDF5_MSG('HDF opening dataset')
  CALL h5dopen_f(file_id,dsetname,dset_id,error)
  ! Create scalar data space for the attribute.
  DEBUG_HDF5_MSG('HDF creating attribute dataspace')

  write(msgstr,"('  Param length: ',I0)") adims(1)
  DEBUG_HDF5_MSG(trim(msgstr))

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  ! CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  ! CALL h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  DEBUG_HDF5_MSG('HDF creating attribute')
  CALL h5acreate_f(dset_id, trim(aname), H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)

  ! Write the attribute data.
  write(msgstr,"('HDF writing attribute <',A,'> = ',F16.10)") trim(aname), data_array(1)
  DEBUG_HDF5_MSG(trim(msgstr))
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data_array, adims, error)

  ! Close the attribute. 
  DEBUG_HDF5_MSG('HDF closing attribute handle')
  CALL h5aclose_f(attr_id, error)
  DEBUG_HDF5_MSG('HDF closing attribute space handle')
  CALL h5sclose_f(aspace_id, error)
  DEBUG_HDF5_MSG('HDF closing dataset handle')
  CALL h5dclose_f(dset_id,error)
  DEBUG_HDF5_MSG('HDF closing file handle')
  CALL h5fclose_f(file_id,error)

  DEBUG_HDF5_MSG_WS('-------( HDF PARAMS DEBUG FINISHED )---------')

end subroutine lbe_write_cp_param_phdf5

!>This subroutine reads the checkpointing parameters from a file.
!>
!>It is called from root (rank 0) of the checkpointing subroutines.
subroutine lbe_read_cp_param_phdf5(data,filename,aname)
  implicit none

  character(len=*) :: filename  ! File name
  character(len=*) aname  ! Attribute name

  character(len=8), parameter :: dsetname='OutArray' !Dataset name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer(hid_t)     :: type_id      ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer            :: error ! Error flag

  real(kind=rk) :: data
  real(kind=rk), dimension(1) :: data_array
  character(len=256) :: msgstr


  ! Because HDF5 wants arrays
  adims(1) = 1


  DEBUG_HDF5_MSG_WS('--------( HDF PARAMS DEBUG STARTED )---------')


  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  DEBUG_HDF5_MSG(trim(msgstr))


  DEBUG_HDF5_MSG('HDF opening file')
  CALL h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  DEBUG_HDF5_MSG('HDF opening dataset')
  CALL h5dopen_f(file_id,dsetname,dset_id,error)

  ! Create datatype for the attribute.
  ! CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  ! CALL h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  write(msgstr,"('HDF opening attribute with name <',A,'>')") trim(aname)
  CALL log_msg(trim(msgstr),.false.)
  CALL h5aopen_by_name_f(dset_id, ".", trim(aname), attr_id, error)
  DEBUG_HDF5_MSG('HDF opening attribute dataspace')
  CALL h5aget_space_f(attr_id, aspace_id, error)

  ! Read the attribute data.
  DEBUG_HDF5_MSG('HDF reading attribute')
  CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data_array, adims, error)
  data = data_array(1)

  write(msgstr,"('HDF read attribute <',A,'> = ',F16.10)") trim(aname), data
  DEBUG_HDF5_MSG(trim(msgstr))

  ! Close the attribute. 
  DEBUG_HDF5_MSG('HDF closing attribute handle')
  CALL h5aclose_f(attr_id, error)
  DEBUG_HDF5_MSG('HDF closing attribute space handle')
  CALL h5sclose_f(aspace_id, error)
  DEBUG_HDF5_MSG('HDF closing dataset handle')
  CALL h5dclose_f(dset_id,error)
  DEBUG_HDF5_MSG('HDF closing file handle')
  CALL h5fclose_f(file_id,error)

  DEBUG_HDF5_MSG_WS('-------( HDF PARAMS DEBUG FINISHED )---------')

end subroutine lbe_read_cp_param_phdf5


!>This subroutine adds the HDF metadata to each postprocessed file.
!>
!>It is called from root (rank 0) of the subroutines dump_..._phdf5.
!>It calls the subroutine lbe_get_input_phdf5 to get all of the input file.
subroutine lbe_write_attr_phdf5(filename,dsetname,aname)
  implicit none

  character(len=256) :: filename  ! File name
  character(len=8)   :: dsetname  ! Dataset name
  character(len=256) :: aname  ! Attribute name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer(hid_t)     :: atype_id      ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer(size_t)    :: attrlen = 80   ! Length of the attribute string
  integer(hsize_t), dimension(7) :: data_dims
  integer            :: error ! Error flag

  integer :: count = 1    !counts how often this subroutine is called
  integer :: blockstart, blockend

  adims = n_meta_total
  data_dims(1) = n_meta_total

  if (count .eq. 1 ) then     !Otherwise data should already exist

    blockstart = 1
    blockend   = n_meta_info
    fixed(blockstart:blockend) = metadata%infodata

    blockstart = blockstart + n_meta_info
    blockend   = blockend   + n_meta_mpi
    fixed(blockstart:blockend) = metadata%mpidata

    blockstart = blockstart + n_meta_mpi
    blockend   = blockend   + n_meta_input
    fixed(blockstart:blockend) = metadata%inputfile

    blockstart = blockstart + n_meta_input
    if (arg_input_dfile_p > 0) then
      blockend   = blockend   + n_meta_input_diff
      fixed(blockstart:blockend) = metadata%inputfile_diff
      blockstart = blockstart + n_meta_input_diff
    endif

#ifdef MD
    blockend = blockend + n_meta_input_md
    fixed(blockstart:blockend) = metadata%inputfile_md
#ifdef ELEC
    blockstart = blockstart + n_meta_input_md
    blockend = blockend + n_meta_input_elec
    fixed(blockstart:blockend) = metadata%inputfile_elec
#endif
#endif

    ! Cleanup
    deallocate(metadata%infodata, metadata%mpidata, metadata%inputfile)
    if (arg_input_dfile_p > 0) then
      deallocate(metadata%inputfile_diff)
    endif
#ifdef MD
    deallocate(metadata%inputfile_md)
#ifdef ELEC
    deallocate(metadata%inputfile_elec)
#endif
#endif

  endif

  count = count + 1

  CALL h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  CALL h5dopen_f(file_id,dsetname,dset_id,error)
  ! Create scalar data space for the attribute.
  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  CALL h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)
  ! Write the attribute data.
  CALL h5awrite_f(attr_id, atype_id, fixed, data_dims, error)
  ! Close the attribute. 
  CALL h5aclose_f(attr_id, error)
  ! Terminate access to the data space.
  CALL h5sclose_f(aspace_id, error)
  CALL h5dclose_f(dset_id,error)
  CALL h5fclose_f(file_id,error)
end subroutine lbe_write_attr_phdf5

subroutine lbe_make_metadata_hdf5()
  if (.not. made_metadata_hdf5) then
    CALL log_msg("  Creating HDF5 metadata (run only once)...",.false.)
    CALL lbe_more_metadata_phdf5()
    CALL lbe_mpidata_phdf5()
    CALL lbe_get_input_phdf5()
    if (arg_input_dfile_p > 0) then
      CALL lbe_get_input_diff_phdf5()
    endif
#ifdef MD
    CALL lbe_get_input_md_phdf5()
#ifdef ELEC
    CALL lbe_get_input_elec_phdf5()
#endif
#endif

    n_meta_total = n_meta_input + n_meta_input_diff + n_meta_mpi + n_meta_info
#ifdef MD
    n_meta_total = n_meta_total + n_meta_input_md
#ifdef ELEC
    n_meta_total = n_meta_total + n_meta_input_elec
#endif
#endif

    allocate(fixed(n_meta_total))   !Never deallocated
    made_metadata_hdf5 = .true.
  endif
end subroutine lbe_make_metadata_hdf5

!>Called by sb postprocessing once
!>Subroutine to get the username and the compiler flags
!>HDF5 metadata of each postprocessed output file
subroutine lbe_more_metadata_phdf5()
implicit none
  character(len=8) :: username
  integer :: pos
  integer :: nfixed = 6 ! Magic number corresponding to the number of lines hardcoded below

  ! Get the username
  CALL getenv('USER',username)

  if (arg_input_dfile_p > 0) then
    nfixed = nfixed + 1
  endif

  n_meta_info = nfixed + nflags
  allocate(metadata%infodata(n_meta_info))

  metadata%infodata(1) = '--------------------( Metadata )----------------------'
  metadata%infodata(2) = trim(lbeversion)
  metadata%infodata(3) = '  HDF version: '//trim(hdfversion)
  metadata%infodata(4) = 'User: '//username
  metadata%infodata(5) = 'Input-file: '//trim(inp_file)
  if (arg_input_dfile_p > 0) then
    metadata%infodata(6) = 'Differential input-file: '//trim(arg_input_dfile)
  endif
  metadata%infodata(nfixed) = '-----------------( Compiler flags used )--------------'

  pos = nfixed + 1

  ! Use the 'find_flags.sh' script to easily check if this is still up to date!
  ! Also, FORTRAN string handling still sucks
#ifdef BOUNCEBACK
  metadata%infodata(pos) = ' BOUNCEBACK'
  pos = pos + 1
#endif
#ifdef BUGGYIFORT11
  metadata%infodata(pos) = ' BUGGYIFORT11'
  pos = pos + 1
#endif
#ifdef BUGGYSENDINCOLLECT
  metadata%infodata(pos) = ' BUGGYSENDINCOLLECT'
  pos = pos + 1
#endif
#ifdef DEBUG_CHECKPOINT
  metadata%infodata(pos) = ' DEBUG_CHECKPOINT'
  pos = pos + 1
#endif
#ifdef DEBUG_HDF5
  metadata%infodata(pos) = ' DEBUG_HDF5'
  pos = pos + 1
#endif
#ifdef DEBUG_HDF5_TIMING
  metadata%infodata(pos) = ' DEBUG_HDF5_TIMING'
  pos = pos + 1
#endif
#ifdef DEBUG_LE
  metadata%infodata(pos) = ' DEBUG_LE'
  pos = pos + 1
#endif
#ifdef DEBUG_MPI
  metadata%infodata(pos) = ' DEBUG_MPI'
  pos = pos + 1
#endif
#ifdef DEBUG_REPORTMDCOMM
  metadata%infodata(pos) = ' DEBUG_REPORTMDCOMM'
  pos = pos + 1
#endif
#ifdef DIST
  metadata%infodata(pos) = ' DIST'
  pos = pos + 1
#endif
#ifdef ELEC
  metadata%infodata(pos) = ' ELEC'
  pos = pos + 1
#endif
#ifdef FASTBDIST2
  metadata%infodata(pos) = ' FASTBDIST2'
  pos = pos + 1
#endif
#ifdef HDF5_FLIP
  metadata%infodata(pos) = ' HDF5_FLIP'
  pos = pos + 1
#endif
#ifdef LB3D
  metadata%infodata(pos) = ' LB3D'
  pos = pos + 1
#endif
#ifdef MCMP
  metadata%infodata(pos) = ' MCMP'
  pos = pos + 1
#endif
#ifdef MD
  metadata%infodata(pos) = ' MD'
  pos = pos + 1
#endif
#ifdef MPI_ALLGV_FASTER_THAN_GV
  metadata%infodata(pos) = ' MPI_ALLGV_FASTER_THAN_GV'
  pos = pos + 1
#endif
#ifdef NOCALLSYSTEM
  metadata%infodata(pos) = ' NOCALLSYSTEM'
  pos = pos + 1
#endif
#ifdef NOEDGESTEP
  metadata%infodata(pos) = ' NOEDGESTEP'
  pos = pos + 1
#endif
#ifdef NOIEEEARITHMETIC
  metadata%infodata(pos) = ' NOIEEEARITHMETIC'
  pos = pos + 1
#endif
#ifdef NOISNAN
  metadata%infodata(pos) = ' NOISNAN'
  pos = pos + 1
#endif
#ifdef NOSURFACTANT
  metadata%infodata(pos) = ' NOSURFACTANT'
  pos = pos + 1
#endif
#ifdef OLDRRFORCE
  metadata%infodata(pos) = ' OLDRRFORCE'
  pos = pos + 1
#endif
#ifdef PARTICLESTRESS
  metadata%infodata(pos) = ' PARTICLESTRESS'
  pos = pos + 1
#endif
#ifdef RELTIME
  metadata%infodata(pos) = ' RELTIME'
  pos = pos + 1
#endif
#ifdef RWALK
  metadata%infodata(pos) = ' RWALK'
  pos = pos + 1
#endif
#ifdef SINGLEFLUID
  metadata%infodata(pos) = ' SINGLEFLUID'
  pos = pos + 1
#endif
#ifdef USEHDF
  metadata%infodata(pos) = ' USEHDF'
  pos = pos + 1
#endif
#ifdef USE_HDF5_INDEPENDENT_IO
  metadata%infodata(pos) = ' USE_HDF5_INDEPENDENT_IO'
  pos = pos + 1
#endif
#ifdef USE_IBM_LARGEBLOCK_IO
  metadata%infodata(pos) = ' USE_IBM_LARGEBLOCK_IO'
  pos = pos + 1
#endif
#ifdef USENEWMASSC
  metadata%infodata(pos) = ' USENEWMASSC'
  pos = pos + 1
#endif
#ifdef USEOLDROCK
  metadata%infodata(pos) = ' USEOLDROCK'
  pos = pos + 1
#endif
#ifdef USEXDRF
  metadata%infodata(pos) = ' USEXDRF'
  pos = pos + 1
#endif
#ifdef VECTORIZE
  metadata%infodata(pos) = ' VECTORIZE'
  pos = pos + 1
#endif
#ifdef WALLCONST
  metadata%infodata(pos) = ' WALLCONST'
  pos = pos + 1
#endif
#ifdef XDRROCKWET
  metadata%infodata(pos) = ' XDRROCKWET'
  pos = pos + 1
#endif
end subroutine lbe_more_metadata_phdf5

subroutine lbe_mpidata_phdf5()
  implicit none
  character(len=16) :: ccdims
  character(len=6)  :: cnprocs
  character(len=80) :: decomposition

  n_meta_mpi = 6
  allocate(metadata%mpidata(n_meta_mpi))

  write(decomposition,"('Decomposition: ',I0,' x ',I0,' x ', I0)") nx, ny, nz
  write(ccdims,"(I0, ' x ',I0, ' x ', I0)") cdims(1), cdims(2), cdims(3)
  write(cnprocs,'(I0)') nprocs

  metadata%mpidata(1) = '-----------------( MPI based data )-------------------'
  metadata%mpidata(2) = 'Hostname: '//hname
  metadata%mpidata(3) = 'Simulation started: '//startsimul
  metadata%mpidata(4) = 'Total number of processors: '//cnprocs
  metadata%mpidata(5) = 'Processors using a '//trim(ccdims)//' grid.'
  metadata%mpidata(6) = decomposition

end subroutine lbe_mpidata_phdf5

!>Called by sb postprocessing once
!>Subroutine to read in all of the input file in order to add it to the
!>HDF5 metadata of each postprocessed output file
!>
!>It will get all of the input-file till end of file is reached.
!>A comment of any length can be added at the very end of the input-file.
!>It won't be read by any of the namelists, just added to the metadata.
!>
!>IMPORTANT: The way the data are read in and the array is defined, it cannot
!>hold more than 78 characters. That is why you HAVE TO TYPE 'RETURN'
!>after each line of the comment. Otherwise you will loose all information
!>longer than 78 characters.
!> 2010-05-10 - Added a copy of the function for the MD input file
#ifdef MD
subroutine lbe_get_input_md_phdf5()
  implicit none
  character(len=78) :: counterchar
  integer :: iunit, unit
  integer :: i, eof
  logical :: op

  n_meta_input_md = 4 ! Starts at four to make room for the fixed lines
  eof = 0

  ! Assume that inp_file_md has been defined.
  inquire(FILE = trim(inp_file)//'.md', OPENED = op, NUMBER=iunit)

  !
  ! Usually the input file is already open; since namelists are used in order to
  ! read in the various input data, no crucial position is lost here by
  ! closing the file and re-opening it later to continue reading the namelists.
  ! But here we rewind the file instead!
  !
  if(op)then
  rewind(iunit)
  else
  iunit = 10
  open(iunit,file=trim(inp_file)//'.md',status='unknown',form='formatted')
  endif

  do while(eof .eq. 0)
  read(iunit,'(a78)',iostat = eof) counterchar
    if(eof.lt.0) then
      if(.not.op)then
        close(iunit)
      endif
      exit
    endif
  n_meta_input_md = n_meta_input_md + 1
  enddo

  if(.not.op)then
    open(iunit,file=trim(inp_file)//'.md',status='unknown',form='formatted')
  else
    rewind(iunit)
  endif

  allocate(metadata%inputfile_md(n_meta_input_md))

  metadata%inputfile_md(1) = ''
  metadata%inputfile_md(2) = '---------------( Start input file MD )----------------'

  do i = 3, n_meta_input_md - 2
    read(iunit,'(a78)')metadata%inputfile_md(i)
  enddo

  metadata%inputfile_md(n_meta_input_md-1) = '--------------( End input file MD )---------------'
  metadata%inputfile_md(n_meta_input_md) = ''

  if(.not.op)then
    close(iunit)
  endif

end subroutine lbe_get_input_md_phdf5
! endif MD
#endif

!>Called by sb postprocessing once
!>Subroutine to read in all of the input file in order to add it to the
!>HDF5 metadata of each postprocessed output file
!>
!>It will get all of the input-file till end of file is reached.
!>A comment of any length can be added at the very end of the input-file.
!>It won't be read by any of the namelists, just added to the metadata.
!>
!>IMPORTANT: The way the data are read in and the array is defined, it cannot
!>hold more than 78 characters. That is why you HAVE TO TYPE 'RETURN'
!>after each line of the comment. Otherwise you will loose all information
!>longer than 78 characters.
!> 2010-05-10 - Added a copy of the function for the MD input file
#ifdef ELEC
subroutine lbe_get_input_elec_phdf5()
  implicit none
  character(len=78) :: counterchar
  integer :: iunit, unit
  integer :: i, eof
  logical :: op

  n_meta_input_elec = 4 ! Starts at four to make room for the fixed lines
  eof = 0

  ! Assume that inp_file_elec has been defined.
  inquire(FILE = trim(inp_file)//'.elec', OPENED = op, NUMBER=iunit)

  !
  ! Usually the input file is already open; since namelists are used in order to
  ! read in the various input data, no crucial position is lost here by
  ! closing the file and re-opening it later to continue reading the namelists.
  ! But here we rewind the file instead!
  !
  if(op)then
  rewind(iunit)
  else
  iunit = 10
  open(iunit,file=trim(inp_file)//'.elec',status='unknown',form='formatted')
  endif

  do while(eof .eq. 0)
  read(iunit,'(a78)',iostat = eof) counterchar
    if(eof.lt.0) then
      if(.not.op)then
        close(iunit)
      endif
      exit
    endif
  n_meta_input_elec = n_meta_input_elec + 1
  enddo

  if(.not.op)then
    open(iunit,file=trim(inp_file)//'.elec',status='unknown',form='formatted')
  else
    rewind(iunit)
  endif

  allocate(metadata%inputfile_elec(n_meta_input_elec))

  metadata%inputfile_elec(1) = ''
  metadata%inputfile_elec(2) = '--------------( Start input file ELEC )---------------'

  do i = 3, n_meta_input_elec - 2
    read(iunit,'(a78)')metadata%inputfile_elec(i)
  enddo

  metadata%inputfile_elec(n_meta_input_elec-1) = '-------------( End input file ELEC )--------------'
  metadata%inputfile_elec(n_meta_input_elec) = ''

  if(.not.op)then
    close(iunit)
  endif

end subroutine lbe_get_input_elec_phdf5
! end ELEC
#endif

subroutine lbe_get_input_diff_phdf5()
  implicit none
  character(len=78) :: counterchar
  integer :: iunit, unit
  integer :: i, eof
  logical :: op

  n_meta_input_diff = 4 ! Starts at four to make room for the fixed lines
  eof = 0

  ! Assume that inp_file has been defined.
  inquire(FILE = arg_input_dfile, OPENED = op, NUMBER=iunit)

  !
  ! Usually the input file is already open; since namelists are used in order to
  ! read in the various input data, no crucial position is lost here by
  ! closing the file and re-opening it later to continue reading the namelists.
  ! But here we rewind the file instead!
  !
  if(op)then
  rewind(iunit)
  else
  iunit = 10
  open(iunit,file=arg_input_dfile,status='unknown',form='formatted')
  endif

  do while(eof .eq. 0)
  read(iunit,'(a78)',iostat = eof) counterchar
    if(eof.lt.0) then
      if(.not.op)then
        close(iunit)
      endif
      exit
    endif
  n_meta_input_diff = n_meta_input_diff + 1
  enddo

  if(.not.op)then
    open(iunit,file=arg_input_dfile,status='unknown',form='formatted')
  else
    rewind(iunit)
  endif

  allocate(metadata%inputfile_diff(n_meta_input_diff))

  metadata%inputfile_diff(1) = ''
  metadata%inputfile_diff(2) = '----------( Start differential input file )-----------'

  do i = 3, n_meta_input_diff - 2
    read(iunit,'(a78)')metadata%inputfile_diff(i)
  enddo 

  metadata%inputfile_diff(n_meta_input_diff-1) = '-----------( End differential input file )------------'
  metadata%inputfile_diff(n_meta_input_diff) = ''

  if(.not.op)then
    close(iunit)
  endif

end subroutine lbe_get_input_diff_phdf5

subroutine lbe_get_input_phdf5()
    integer,parameter :: iunit=10
    character(len=78) :: counterchar
    integer :: i, eof

    n_meta_input = 4 ! Starts at four to make room for the fixed lines
    eof = 0
    open (iunit,file=inp_file,status='unknown',form='formatted')
    do while (eof .eq. 0)
       read (iunit,'(a78)',iostat = eof) counterchar
       if (eof<0) exit
       n_meta_input = n_meta_input + 1
    enddo
    close(iunit)

    allocate(metadata%inputfile(n_meta_input))

    metadata%inputfile(1) = ''
    metadata%inputfile(2) = &
         &'-----------------( Start input file )-----------------'

    open(iunit,file=inp_file,status='unknown',form='formatted')
    do i = 3,n_meta_input-2
       read (iunit,'(a78)') metadata%inputfile(i)
    enddo
    close(iunit)

    metadata%inputfile(n_meta_input-1) = &
         &'------------------( End input file )------------------'
    metadata%inputfile(n_meta_input) = ''
end subroutine lbe_get_input_phdf5

#endif
end module lb3d_io_hdf5_module
