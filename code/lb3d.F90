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

!> The main program is merely calling the appropriate subroutines for
!> initialisation, the main time-loop and finalisation of the algorithm 
!> \todo Add Doxygen mainpage here?
!> \todo Add complete Authorlist here
!> \todo LGPL headers
program lb3d

  use lb3d_initialise_module
  use lb3d_time_loop_module
  use lb3d_finalise_module

  implicit none

  integer               :: err
  
  call lb3d_init
  call lb3d_time_loop
  call lb3d_finalise

end program

