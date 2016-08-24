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

!> Contains subroutines for periodic boundary conditions
module lb3d_bc_periodic_module

  use lb3d_lattice_module
  use lb3d_mpi_parameters_module
  use lb3d_log_module

  implicit none

contains

  !> Perform halo exchange using the mpi datatypes defined in \c lb3d_mpi_module
  !> \todo Optimisation: Fine-grain the function to enable calls per lattice direction.
  subroutine lbe_halo_exchange(N)
    type(lbe_site),intent(inout) :: &
         &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer k,ierror
    integer status(MPI_STATUS_SIZE)

#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lbe_halo_exchange',.false.) 
#endif

    do k = 1,3
       call mpi_sendrecv&   ! send "downward"
            &(N,1,ls_mpitype(k),nnprocs(k,1),0&
            &,N,1,ur_mpitype(k),nnprocs(k,2),0&
            &,comm_cart,status,ierror)
       call mpi_sendrecv&   ! send "upward"
            &(N,1,us_mpitype(k),nnprocs(k,2),0&
            &,N,1,lr_mpitype(k),nnprocs(k,1),0&
            &,comm_cart,status,ierror)
    end do
  end subroutine lbe_halo_exchange


end module lb3d_bc_periodic_module
