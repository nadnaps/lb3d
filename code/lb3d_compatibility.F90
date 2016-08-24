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

!> Contains subroutines to check compatibility of selected parameters
module lb3d_compatibility_module

  use lb3d_config_module
  use lb3d_mpi_module

  implicit none

  contains

    !> Rudimentary routine to amass known compatibility problems
    !> \todo This can be handled more elegantly in a fine-grained
    !> Config read-in
    subroutine lb3d_compatibility_check

#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lb3d_compatibility_check',.false.) 
#endif
    
    if (myrankw.eq.0) then

       CALL log_msg_ws("-----------( Starting compatibility check )--",.false.)



       ! General warnings for untested functionality

       if (boundary_cond > 0) then
          call log_msg_ws('WARNING: Automated boundary patterns are untested.',.false.)
       end if
       if (sci_fluxz.eqv..true.) then 
          call log_msg_ws('WARNING: Flux output is untested.',.false.)
       end if
       if (sci_profile.eqv..true.) then 
          call log_msg_ws('WARNING: Profile output is untested.',.false.)
       end if
       if (sci_arrstats.eqv..true.) then
          call log_msg_ws('WARNING: Velocity array statistics output is untested.',.false.)
       end if
       if (sci_pressure.eqv..true.) then
          call log_msg_ws('WARNING: Pressure output is untested.',.false.)
       end if
       if (write_AVS_fld.eqv..true.) then
          call log_msg_ws('WARNING: AVS field output is untested.',.false.)
       end if
       if (use_lbe_force.eqv..true.) then
          call log_msg_ws('WARNING: lbe_force-style forcing is untested.',.false.)
       end if

       ! General warnings for parameter ranges

       if (g_br < 0.) then 
          call log_msg_ws('WARNING: Attractive binary interactions bear no meaning.',.false.)
       end if

       if (g_br .ge. 0.15) then 
          write(msgstr,"('g_br = ',F16.10,' is very high. Simulations are stable in the range 0.0 .. 0.14')") g_br
          CALL log_msg_ws(trim(msgstr),.false.)
       end if

       if (g_rr < 0.) then 
          call log_msg_ws('WARNING: repulsive self-interactions bear no meaning.',.false.)
       end if

       if (g_rr .ge. 0.15) then 
          write(msgstr,"('g_br = ',F16.10,' is very high. Simulations are stable in the range 0.0 .. 0.15')") g_br
          CALL log_msg_ws(trim(msgstr),.false.)
       end if

#ifdef SINGLEFLUID
       ! Warnings for settings incompatible with a single fluid
       

       !call log_msg('In lb3d_compatibility_check',.false.) 
       !endif #ifdef SINGLEFLUID 
#endif

#ifndef SINGLEFLUID
       ! Warnings for settings incompatible with a binary mixture


#ifndef NOSURFACTANT
       ! Warnings for settings incompatible with a ternary mixture
       if (boundary_cond.ne.0.and.(inv_fluid.eq.5.or.inv_fluid.eq.6)) then
          call log_msg_ws('FATAL: Lees Edwards BC are not supported for ternary mixtures',.false.)
          call abend
       end if
       if (obs_file.ne.'empty.dat') then 
          call log_msg_ws('FATAL: Obstacles are not supported for ternary mixtures',.false.)
          call abend
       end if

       !endif #ifndef NOSURFACTANT
#endif
       !endif #ifndef SINGLEFLUID
#endif


    end if
    end subroutine lb3d_compatibility_check

end module lb3d_compatibility_module
