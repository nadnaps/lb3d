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

!> Contains MPI related parameters and provides the mpi headers
!> \todo More consistent naming for the mpi part
module lb3d_mpi_parameters_module


! use lb3d_lattice_module, only:nd


  implicit none
  include "mpif.h"

  !FIXME not sure where to put the dimensionality, took it out, because it
  !is creating difficulties with the dependencies and is not implemented 
  !consistently anyway
  integer, parameter :: nd = 3

  save ! Static variables.
 ! integer :: start_time_t, end_time_t
  character(len=16) :: input_source = '.input-file' !< Name of the .input-file
  !> \{
  !> \name Variables holding info about the processor topology
  !>
  !> Set to the number of CPUs taking part in the simulation.
  integer :: nprocs
  integer :: myrankw            !< rank of the CPU (MPI_COMM_WORLD)
  integer :: myrankc            !< rank of the CPU (Comm_Cart)
  integer :: myfarm             !< Task farming?
  integer :: cdims(nd) = 0      !< dimensions of the CPU topology.
  integer :: nnprocs(nd, 2)     !< ranks of the nearest-neighbour CPUs
  integer :: Comm_Cart      !< The Cartesian CPU topology communicator
  integer :: ccoords(nd) = 0    !< My coords in the Cartesian grid.
  integer :: start(nd) !< global coordinates of the start of the CPU's subdomain
  ! integer :: nx,ny,nz         !< Dimensions of my chunk of the lattice.
  integer :: tnx,tny,tnz        !< Total size of lattice.
  !> \}

  integer :: lbe_site_mpitype !< mpi datatype representing \c lbe_site

  !> \{
  !> \name MPI data types
  !> mpi datatypes representing the lower (l)/upper (u) part of the halo that
  !> is sent (s)/received (r) during the exchange for each direction (indices
  !> 1/2/3 represent x/y/z). Datatype definitions are relative to \c N as a
  !> whole.
  integer :: ls_mpitype(3),us_mpitype(3),lr_mpitype(3),ur_mpitype(3)
  !> \}


  !> Determines whether or not the Cartesian topology is periodic in
  !> each given dimension.
  !>
  !> \note This may change for different boundary conditions.
  logical, dimension(3) :: periodic_p = .true.


  !> \{
  !> \name MPI Derived Datatypes
  !>
  !> These datatypes describe the surface of interface between
  !> two processors' domains.
  integer :: interf_px,interf_mx
  integer :: interf_py,interf_my
  integer :: interf_pz,interf_mz
  !> \}

  !> MPI tag for messages related to initialisation
  integer, parameter :: tag_init = 1
  !> \{
  !> \name MPI tags used in the halo exchange
  integer, parameter :: tag_px = 2,tag_mx = 3
  integer, parameter :: tag_py = 4,tag_my = 5
  integer, parameter :: tag_pz = 6,tag_mz = 7
  !> \}
  integer, parameter :: tag_rock = 8
  !> MPI tag for messages related to postprocessing
  integer, parameter :: tag_post = 15


  !> \{
  !> \name EYB HDF5 added
  character(len=MPI_MAX_PROCESSOR_NAME) :: hname !< For hostname
  character(len=24) :: startsimul                !< For starttime
  !> \}

end module lb3d_mpi_parameters_module
