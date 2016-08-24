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

!> Contains subroutines setting up and providing random numbers
!> \todo Introduce a wrapper here to make RNG fully abstract.
module lb3d_random_module

  use lb3d_global_module
  use lb3d_log_module
  use lb3d_config_module, only: seed
  use lb3d_mpi_parameters_module, only: myrankc
  implicit none
  
contains

!      ==================================================================
  subroutine lb3d_random_init (stage)
!      ==================================================================
  
    integer :: stage

#ifdef LB3D_DEBUG_INFO
    write(msgstr,"('In lb3d_random_init stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif
    
    select case (stage)
    case (2)
       call lbe_unique_seeds()
    case default
#ifdef LB3D_DEBUG_INFO
       call log_msg('nothing to do.',.false.)
#endif
    end select           


  end subroutine lb3d_random_init

 
  subroutine lbe_unique_seeds()
    integer, dimension(10) :: seedarray
    integer                :: isize
    ! Now make sure that the initial random seed given to each processor
    ! is unique. BEWARE this is only a temporary measure - I know not
    ! if this will produce any hidden correlations between the random
    ! numbers. I leave this to the experts - it should not be difficult
    ! to modify the routine below.
    
    ! seed = seed * (myrankc + 1) * 10000
    
    ! To ensure the random number seed remains in range
    seed = MOD(seed * (myrankc + 1) * 10101, 30000) 
    
    ! Seed the PRNG.
    !! FIXME not sure what this is doing, exactly - taken
    !! from ME3D. Will isize ever be >10?
    seedarray = seed
    CALL random_seed(size=isize)
    CALL random_seed(put=seedarray(1:isize))
    CALL random_seed(get=seedarray(1:isize))
  end subroutine lbe_unique_seeds

 
end module lb3d_random_module

