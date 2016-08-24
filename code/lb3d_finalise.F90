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

module lb3d_finalise_module
!========================================================================    

use lb3d_config_module, only:boundary_cond, inv_fluid
use lb3d_bc_leesedwards_module,only: le_cleanup
use lb3d_timer_module, only: sum_timer, calculate_performance
#ifdef USEHDF
use lb3d_io_hdf5_module,only:lb3d_io_shutdown_hdf5
#endif
use lb3d_mpi_module, only:tnx,tny,tnz,finalizempicomms
use lb3d_log_module

implicit none

contains

  !> Calls subroutines to exit gracefully (MPI and HDF5)
  !> Triggers evaluation of timers.
  subroutine lb3d_finalise
!========================================================================    
    call sum_timer((/6/)) ! dump timer evaluation to standard output (unit 6)
    
    if ((boundary_cond.ne.0).and.(inv_fluid.eq.5.or.inv_fluid.eq.6)) then
       CALL le_cleanup()
    endif
    
#ifdef USEHDF
    CALL lb3d_io_shutdown_hdf5()
#endif

    CALL FinalizeMPIcomms()

    if (myrankc == 0) then
       CALL calculate_performance
    end if
    
    CALL log_msg("Exiting...",.false.)
    
  end subroutine lb3d_finalise

end module lb3d_finalise_module
