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

!> Contains subroutines wrapping the initialisation and application of 
!> boundary conditions (peridiodic and influx)
!> \todo Make consistent for the bounce back boundaries
module lb3d_bc_module

  use lb3d_config_module, only: acccoef,bcsel, boundary_cond, shear_u, shear_tmp,init_cond, restore, inv_fluid
  use lb3d_timer_module
  use lb3d_log_module
 
  use lb3d_lattice_module
  use lb3d_analysis_module

  use lb3d_bc_invasion_module,only:lbe_invade
  use lb3d_force_module, only:lbe_add_force_halo
  
  use lb3d_bc_periodic_module
  use lb3d_bc_leesedwards_module
  use lb3d_bc_bounce_module

  implicit none

contains

  !> Calling subroutines needed to initialise boundary conditions
  !> Executing initial halo-exchange
  subroutine lb3d_bc_init (stage)

    integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_bc_init stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif    

    select case (stage)

    case (4) ! post mem alloc
       
       ! If we're shearing with Lees Edwards bc's
        if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
           ! Correct shear_u so that velocity on side planes in
           ! simple fluid tends to value in input file.
           shear_u = (tnx/(tnx-1.d0))*shear_u
           shear_tmp = shear_u
           if (init_cond==7.or.restore) then
              ! If you restore from checkpoint, subtract one shear_sum
              ! because it will be added a second time in
              ! le_halo_exchange.
              shear_sum = shear_sum - shear_u
           endif

           CALL le_init()
           ! Work out how processors lie
           CALL le_neighbours()
           ! Do special halo exchange
           CALL le_halo_exchange(N)
           call lb3d_update_density_loop (N)
           call lb3d_update_velocity_loop (N)
        else
           CALL lbe_halo_exchange(whole_N)
        endif
 
    case default


    end select    

  end subroutine lb3d_bc_init


  !> Calling boundary condition routines for influx and periodic 
  !> boundary conditions
  subroutine lb3d_bc_apply (stage)
  
    integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_bc_apply stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif

    select case (stage)

    case (6) ! pre advection
       
       call start_timer(ti_halo)
       if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
          shear_u = shear_tmp*cos(shear_omega*nt)
          if (myrankc==0.and.shear_omega>0.0_rk) then
             write(msgstr,"('SHEAR_U is set to ',ES15.8)") shear_u
             call log_msg(trim(msgstr),.false.)
          endif
          call le_neighbours()
          call le_halo_exchange(N)
          call lb3d_update_density_loop(N)
          call lb3d_update_velocity_loop(N)
       else
          call lbe_halo_exchange(whole_N)
       endif
       call stop_timer(ti_halo)
       
 case (8) ! post advection / pre collision

    ! Modified by Maddalena, 19-Oct-2004
    ! do not invade if inv_fluid<0
    ! Useful to have boundary_cond.ne.0  and disabled invasion
    if (boundary_cond/=0.and.inv_fluid>=0) then
       call start_timer(ti_inv)
       call lbe_invade(N)
       call stop_timer(ti_inv)
    end if
    
    call start_timer(ti_halo)
    if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
       call le_halo_exchange(N)
       call lb3d_update_density_loop (N)
       call lb3d_update_velocity_loop (N)
    else
       call lbe_halo_exchange(whole_N)
    endif

    call stop_timer(ti_halo)

    call lbe_add_force_halo()

 case (10) ! post collision / pre advection / default dump

 case default
    !FIXME Error

 end select

end subroutine lb3d_bc_apply

!!$
!!$subroutine lb3d_halo_exchange_wrapper(lbe_N,whole_N,d_adv) 
!!$  type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
!!$  type(lbe_site),intent(inout) :: &
!!$             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
!!$#ifndef NOSURFACTANT
!!$  real(kind=rk),dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
!!$#else
!!$  integer,intent(inout) :: d_adv
!!$#endif
!!$
!!$  call start_timer(ti_halo)
!!$  if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
!!$     shear_u = shear_tmp*cos(shear_omega*nt)
!!$     if (myrankc==0.and.shear_omega>0.0_rk) then
!!$        write(msgstr,"('SHEAR_U is set to ',ES15.8)") shear_u
!!$        call log_msg(trim(msgstr),.false.)
!!$     endif
!!$     call le_neighbours()
!!$     call le_halo_exchange(lbe_N)
!!$#ifndef NOSURFACTANT
!!$     call le_adv_dipole_exchange(d_adv)
!!$#endif
!!$  else
!!$     call lbe_halo_exchange(lbe_N)
!!$     
!!$#ifndef NOSURFACTANT
!!$   !  call lbe_adv_dipole_exchange(d_adv)
!!$     
!!$#endif
!!$        endif
!!$        call stop_timer(ti_halo)
!!$
!!$end subroutine lb3d_halo_exchange_wrapper
!!$

end module lb3d_bc_module
