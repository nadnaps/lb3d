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

!> Contains subroutines implementing bounce back boundary conditions
module lb3d_bc_bounce_module

use lb3d_config_module, only: acccoef,bcsel
use lb3d_lattice_module

implicit none

contains

  !> Wrapper function to provide the possibility of setting 
  !> the bounce back type via the input-file  
  subroutine lb3d_bc_bounce_select (N,x,y,z,Nbuf)
    
    implicit none
    type(lbe_site), dimension(0:,0:,-1:) :: Nbuf
    integer :: x,y,z
    type(lbe_site), dimension(0:,0:,0:) :: N

    select case (bcsel)
    case (0)
       call lb3d_bc_bounce_back_hw (N,x,y,z,Nbuf)

    case default
       !FIXME error
    end select

  end subroutine lb3d_bc_bounce_select

  !> Implements a half-way bounce back boundary condition
  subroutine lb3d_bc_bounce_back_hw (N,x,y,z,Nbuf)

    implicit none
    type(lbe_site), dimension(0:,0:,-1:) :: Nbuf


    integer :: x,y,z,xa,ya,za
    integer :: s
    type(lbe_site), dimension(0:,0:,0:) :: N


    do s=1,nnonrest
       xa = x + cx(s)
       ya = y + cy(s)
       za = z + cz(s)

       if (N(xa,ya,za)%rock_state .ne. 0.) then
          if (Nbuf(x,y,0)%rock_state .eq. 0.) then

             N(x,y,z)%n_r(bounce(s)) = Nbuf(x,y,0)%n_r(s)
#ifndef SINGLEFLUID
             N(x,y,z)%n_b(bounce(s))=Nbuf(x,y,0)%n_b(s)
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%n_s(bounce(s))=Nbuf(x,y,0)%n_s(s)

             N(x,y,z)%da = N(x,y,z)%da &
                  + Nbuf(x,y,0)%n_s(s)*g(s) &
                  * Nbuf(xa,ya,cz(s))%db

#endif

          endif
       endif

    end do

  end subroutine lb3d_bc_bounce_back_hw



end module lb3d_bc_bounce_module
