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

!> Contains subroutines implementing the advection step of the lattice Boltzmann algorithm
!> \todo Move the Buffer allocation to lb3d_lattice
!> \todo Optimisation: Integrate with collision step 
module lb3d_advection_module
!========================================================================

  use lb3d_lattice_module
  use lb3d_bc_module
  use lb3d_timer_module

  use lb3d_analysis_module

  implicit none

contains

  !> Wrapper to remove dependency on lattice information from the time-loop
  subroutine lb3d_advection
    
#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lb3d_advection',.false.) 
#endif
    
    call start_timer(ti_adv)
    call lb3d_advect(N)
    call stop_timer(ti_adv)

  end subroutine lb3d_advection

  !> Contains the saptial advection-step loop for the array \c N
  !> Calling both advection and bounce-back boundary condition
  subroutine lb3d_advect(N)
    implicit none
    integer :: x, y, z
    integer :: xa, ya, za, xb, yb, zb
    integer :: bcz
    integer :: s,i
    type(lbe_site), dimension(0:,0:,0:) :: N

    type(lbe_site), dimension(:,:,:),allocatable :: Nbuf

    real(kind=rk) :: init, finl, diff
    integer :: xx, yy, zz

    real(kind=rk) :: pred, pblue, psurf, qred, qblue, qsurf
    integer :: xp, yp, zp, xq, yq, zq, sp, sq
    integer :: fooflag=0

    if (.not. allocated(Nbuf)) allocate(Nbuf(0:nx+1,0:ny+1,-1:1))

    Nbuf(0:nx+1,0:ny+1,-1) = N(0:nx+1,0:ny+1,0)
    Nbuf(0:nx+1,0:ny+1,0) = N(0:nx+1,0:ny+1,1)

    do z=1,nz
       Nbuf(0:nx+1,0:ny+1,1) = N(0:nx+1,0:ny+1,z+1)
       do y=1,ny
          do x=1,nx
#ifndef NOSURFACTANT
             N(x,y,z)%da = 0.0_rk
#endif


#ifdef COMPAT_BCSEL
             call lb3d_bc_bounce_select (N,x,y,z,Nbuf)
#else
             ! Simple half way bounce back is the default for now
             call lb3d_bc_bounce_back_hw (N,x,y,z,Nbuf)
#endif

             call lb3d_advection_step (N,x,y,z,Nbuf)

             call lb3d_update_density (N,x,y,z)
             call lb3d_update_velocity (N,x,y,z)
             !N(x,y,z)%da = N(x,y,z)%da / N(x,y,z)%rho_s

          end do     !x
       end do        !y
       Nbuf(0:nx+1,0:ny+1,-1)=Nbuf(0:nx+1,0:ny+1,0)
       Nbuf(0:nx+1,0:ny+1,0)=Nbuf(0:nx+1,0:ny+1,1)
    end do           !z

  end subroutine lb3d_advect

  !> Implementing a pull-advection for one lattice site
  !> Calculating the advected dipole moments up to a norm for each lattice vector
  !> writing to \c N%da.
  subroutine lb3d_advection_step (N,x,y,z,Nbuf)

    implicit none 

    integer :: x, y, z
    integer :: xa, ya, za, xb, yb, zb
    integer :: bcz
    integer :: s
    type(lbe_site), dimension(0:,0:,-1:) :: Nbuf
    type(lbe_site), dimension(0:,0:,0:) :: N
    
    do s=1,nnonrest
       xb = x - cx(s)
       yb = y - cy(s)
       zb = z - cz(s)
       bcz = - cz(s)
       
       if (N(x,y,z)%rock_state .eq. 0.) then
          if (Nbuf(xb,yb,bcz)%rock_state .eq. 0.) then

             N(x,y,z)%n_r(s) = Nbuf(xb,yb,bcz)%n_r(s)
             
#ifndef SINGLEFLUID
             N(x,y,z)%n_b(s) = Nbuf(xb,yb,bcz)%n_b(s)
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%n_s(s) = Nbuf(xb,yb,bcz)%n_s(s)
             
             N(x,y,z)%da = N(x,y,z)%da &
                  + (Nbuf(xb,yb,bcz)%n_s(s)*g(s) &
                  * Nbuf(xb,yb,bcz)%db)

#endif
          endif
       end if
    end do  !s

#ifndef NOSURFACTANT
    N(x,y,z)%da = N(x,y,z)%da + Nbuf(x,y,0)%n_s(restvec)*g(restvec) &
         * Nbuf(x,y,0)%db
#endif 

end subroutine lb3d_advection_step



end module lb3d_advection_module

