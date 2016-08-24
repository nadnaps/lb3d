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

!> Contains global variables
module lb3d_global_module
!========================================================================
  
  implicit none
  !FIXME Excerpt from lbe.F90 declaration, to be cleaned up
  integer,parameter  :: stdout = 6

  integer            :: ierror
  integer            :: i,k,x,y
  integer            :: tstart
  character(len=24)  :: datestr
  character(len=256) :: message,msgstr
  logical            :: insanity
  logical            :: minz,maxz
 

  !FIXME move this to lb3d_io??
  !> unit used for input-file handle
  integer, parameter :: input_file_unit = 17
  !> unit used for differential input-file handle
  integer, parameter :: input_dfile_unit = 42

  !FIXME can be defined via pre-compiler constant?
  integer, parameter :: rk=8 !< type parameter (kind) of most of the reals

  !FIXME rather register and count? is this still used at all?
  !> number of lbe particle species
#ifdef SINGLEFLUID
  integer,parameter :: n_spec = 1
#elif NOSURFACTANT
  integer,parameter :: n_spec = 2
#else
  integer,parameter :: n_spec = 3
#endif

  !FIXME Now this is what I call a global (even universal?)
  real(kind=rk),  parameter :: pi = 3.1415926535897932384626433_rk

  !FIXME This has to be made accessible by invade again in an elegant way
  integer,save :: start_time_t, end_time_t
  integer, save  :: timesteps_count = 0 !< Global count of timesteps performed

  integer, save  :: nt = 0             ! No of timesteps since last output - not in namelist

  !FIXME Within the new logic this belongs with the force module
  !> \name common external forcing
  !> \{
  !> has to be set to \c .true. if \c lbe_force should be taken into
  !> account as additional forcing during collision
  logical,save :: use_lbe_force=.false.

  !> additional force on lbe particles, indices are: comp.,spec.,x,y,z
  real(kind=rk),save,allocatable,dimension(:,:,:,:,:) :: lbe_force

end module lb3d_global_module
