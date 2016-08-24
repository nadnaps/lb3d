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

!> Contains the main time-loop calling the different parts implementing
!> the lattice Boltzmann algorithm
!> \todo Optimisation: Move spatial loops here as well?
module lb3d_time_loop_module
!========================================================================
  use lb3d_global_module
  use lb3d_io_module
  use lb3d_timer_module
  use lb3d_force_module
  use lb3d_bc_module, only:lb3d_bc_apply
  use lb3d_advection_module
  use lb3d_collision_module
  use lb3d_helper_module, only:unixtime
  use lb3d_sanity_module
  implicit none

contains

  subroutine lb3d_time_loop
!========================================================================    
    implicit none
    
    integer             :: stage
   

    call lb3d_io_write_data
    call lb3d_io_write_checkpoint
    call lb3d_sanity_check

    timesteps_count = 0

    start_time_t = unixtime()
    call start_timer(ti_total)

    time_loop: do nt = tstart+1, n_iteration
       write(msgstr,"('Starting timestep ',i0)") nt
       call log_msg(trim(msgstr),.false.)

       timesteps_count = timesteps_count + 1
       stage = 6 
 

       call lb3d_force_apply(stage)
       call lb3d_bc_apply(stage)

       stage = 7 ! Advection
       !=================================================================
       
       call lb3d_advection

       stage = 8
       !=================================================================
       
       call lb3d_bc_apply(stage)
       call lb3d_force_apply(stage)

       stage = 9 ! Collision
       !=================================================================

       call lb3d_collision
       

       stage = 10 
       !=================================================================
       call lb3d_io_write_data
       call lb3d_io_write_checkpoint
       call lb3d_sanity_check

    end do time_loop
    
    end_time_t = unixtime()  
    call stop_timer(ti_total)

  end subroutine lb3d_time_loop

end module lb3d_time_loop_module
