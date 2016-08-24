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

!> Provides log functionality
module lb3d_log_module

    use lb3d_mpi_parameters_module, only:myrankc

    implicit none

  contains

    !      ==================================================================
    subroutine lb3d_log_init (stage)
      !      ==================================================================

      integer :: stage
      print*,'log_init',stage

      select case (stage)
      case (0)

      case default

      end select

    end subroutine lb3d_log_init

    !> Writes a timestamp to stdout. ws stands for whitespace
    !> If forall is true, every processor does this, otherwise only rank 0 will
    subroutine log_ws(forall)
      logical :: forall
      CALL log_msg("",forall)
    end subroutine log_ws

    !> Writes a message msg with a timestamp to stdout. ws stands for whitespace
    !> If forall is true, every processor does this, otherwise only rank 0 will
    subroutine log_msg_ws(msg,forall)
      character(len=*)     :: msg
      logical              :: forall
      CALL log_msg("",forall)
      CALL log_msg(msg,forall)
      CALL log_msg("",forall)
    end subroutine log_msg_ws

    !> Writes a message msg with a timestamp to stdout.
    !> If forall is true, every processor does this, otherwise only rank 0 will
    subroutine log_msg(msg,forall)
      character(len=*)     :: msg
      logical              :: forall
      integer,dimension(8) :: datevalues
      character(len=24)    :: datestr

      CALL date_and_time(values=datevalues)

      write(datestr,"(I0.2,':',I0.2,':',I0.2)") datevalues(5), datevalues(6), datevalues(7)

      if (forall) then
         write(*, "(A, ' - (', I6.6, ') ', A)" )  trim(datestr), myrankc, msg
      else if (myrankc == 0) then
         write(*, "(A, ' - ', A)") trim(datestr), msg
      endif

    end subroutine log_msg

  end module lb3d_log_module
