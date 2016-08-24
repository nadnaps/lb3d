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

!> Calls the various initilisation routines
module lb3d_initialise_module
!========================================================================

  use lb3d_io_module     ,only: lb3d_io_init
  use lb3d_config_module ,only: lb3d_config_init
  use lb3d_mpi_module  !  ,only: lb3d_mpi_init
  use lb3d_random_module ,only: lb3d_random_init
  use lb3d_lattice_module,only: lb3d_lattice_init
  use lb3d_bc_module     ,only: lb3d_bc_init
  use lb3d_ic_module     ,only: lb3d_ic_init
  use lb3d_force_module  ,only: lb3d_force_init
  use lb3d_log_module    ,only: lb3d_log_init
  use lb3d_timer_module  ,only: lb3d_timer_init
  use lb3d_compatibility_module

  implicit none
  
contains 

  subroutine lb3d_init
    !====================================================================    
    implicit none
    
    integer             :: stage

    stage = 0
    !====================================================================
    call lb3d_mpi_init(stage)

    call lb3d_log_init(stage)
    call lb3d_timer_init(stage)
    call lb3d_io_init(stage)
    call lb3d_config_init(stage)
    call lb3d_compatibility_check

    call lb3d_random_init(stage)

    call lb3d_bc_init(stage)
    call lb3d_ic_init(stage)
    call lb3d_force_init(stage)

    stage = 1
    !====================================================================

    call lb3d_mpi_init(stage)
    
    stage = 2
    !====================================================================    

    call lb3d_io_init(stage)
    call lb3d_config_init(stage)

    call lb3d_random_init(stage)

    call lb3d_bc_init(stage)
    call lb3d_ic_init(stage)
    call lb3d_force_init(stage)

    stage = 3
    !====================================================================

    call lb3d_lattice_init(stage)

    stage = 4
    !====================================================================

    call lb3d_mpi_init(stage)

    call lb3d_ic_init(stage)
    call lb3d_bc_init(stage)
    call lb3d_force_init(stage)

  end subroutine lb3d_init

end module lb3d_initialise_module

