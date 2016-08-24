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

!> Contains subroutines to setup the boundary geometries, writing to
!> \N%rock_state
module lb3d_init_rock_module

  use lb3d_lattice_module!, only: nnonrest, restvec, nvecs

#ifdef USEHDF
  use lb3d_io_hdf5_module, only: read_rock_phdf5
#endif
#ifdef USEXDRF
  use lb3d_io_xdrf_module, only: read_bit_rock_xdrf_par, read_rock_xdrf_par
#endif

  use lb3d_mpi_module, only: abend, myrankc, recv_rock_par, tnx, tny, tnz, ccoords
  use lb3d_config_module, only: bcsel, boundary_cond, boundary_width, obs_file, obs_folder, nx, ny, nz, rock_colour_r, rock_colour_b
  use lb3d_log_module, only: log_msg

  implicit none
  private

  public lbe_init_rock

contains

!> This routine will handle initialization of the rock geometry.
!> It is the only routine that is exposed to other modules.
subroutine lbe_init_rock(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N

  call log_msg("Setting up rock geometry.",.false.)

  ! First, read rock from file, if any...
  if ( trim(obs_file) .ne. 'empty.dat') then
    call lbe_read_rock(N)
  else
    call log_msg("  Special value obs_file = <empty.dat> : no rock file will be read.",.false.)
  endif

  ! Then, modify, if boundary_cond != 0
  if ( boundary_cond .gt. 0 ) then
    call log_msg("  boundary_cond != 0: Adding additional rock geometry.",.false.)
    call lbe_modify_rock(N)
  else
    call log_msg("  boundary_cond == 0: No additional rock geometry will be added.",.false.)
  endif

#ifndef XDRROCKWET
    ! We only need this if we do not get wettability information from rock
    ! files and do not use read_rock_xdrf_par.
    call lbe_init_rock_colour(N)
#endif

end subroutine lbe_init_rock

!> Wraps calls to various rock reading functions.
subroutine lbe_read_rock(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  character(len=128) :: full_obs_fname
  character(len=128) :: msgstr
  real(kind=rk) :: rockthresh = 0.0

  full_obs_fname = trim(obs_folder)//'/'//trim(obs_file)

  write(msgstr,"('  Reading rock data file <',A,'>')") trim(full_obs_fname)
  CALL log_msg(trim(msgstr),.false.)

  ! Now adjust the obstacle matrix according to the specified
  ! boundary conditions.

  if (index(full_obs_fname,'.xdr') .gt. 0) then
#ifdef USEXDRF
    ! Call XDR-read code
    CALL read_rock_xdrf_par(full_obs_fname,rockthresh,N)

#ifdef DIST
    CALL read_dist_xdrf_par(full_obs_fname,N)
#endif

#ifdef RELTIME
    CALL read_rel_xdrf_par(full_obs_fname,N)
#endif

#else
    ! If USEXDRF is not set:
    CALL log_msg("FATAL ERROR: XDRF disabled, but xdr rock file supplied - can't initialize rock. Aborting...",.false.)
    CALL Abend
#endif

  else if (index(full_obs_fname,'.h5') .gt. 0) then
#ifdef USEHDF
    CALL read_rock_phdf5(full_obs_fname,rockthresh,N)
#else
    CALL log_msg("FATAL ERROR: HDF5 disabled, but hdf rock file supplied - can't initialize rock. Aborting...",.false.)
    CALL Abend
#endif

  else if (index(full_obs_fname,'.bit') .gt. 0) then
#ifdef USEXDRF
    ! Call XDR-read code
    CALL read_bit_rock_xdrf_par(full_obs_fname,N)
#else
    ! If USEXDRF is not set:
    CALL log_msg("FATAL ERROR: XDRF disabled, but xdr rock file supplied - can't initialize rock. Aborting...",.false.)
    CALL Abend
#endif

  else
    if (myrankc == 0) then
      ! Call parallel ASCII read code
      !CALL read_rock_all_par(full_obs_fname,N)
    else
      ! Receive rock from rank zero
      CALL recv_rock_par(N)
    endif
  endif

end subroutine lbe_read_rock

#ifndef XDRROCKWET
!> Sets rock colours if we don't have a separate rock wetting file.
subroutine lbe_init_rock_colour(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x,y,z

  ! Clear the nonrest vectors of each rock site.
  ! FIXME - this could be cleaned up, although it only
  ! gets called once, so it's not that important.

  do z = 1, nz
    do y = 1, ny
      do x = 1, nx
        if (N(x,y,z)%rock_state .ne. 0) then
 
            N(x,y,z)%n_r(:)=0.d0
#ifndef SINGLEFLUID
            N(x,y,z)%n_b(:)=0.d0
#endif

#ifndef NOSURFACTANT
            N(x,y,z)%n_s(:)=0.d0
            N(x,y,z)%da=(/0.d0,0.d0,0.d0/)            
            N(x,y,z)%db=(/0.d0,0.d0,0.d0/)
#endif

          ! If rock_colour is +ve, colour the rock red,
          ! else blue. XXX FIXME Note that this will change - 
          ! a much better way would be to specify occupation
          ! numbers in the input file. You could then have dipole
          ! rocks, although I'm not quite sure what the effect
          ! would be.
          if (rock_colour_r > 0) then
             N(x,y,z)%rock_colour_r = rock_colour_r
          end if
#ifndef SINGLEFLUID
          if (rock_colour_b > 0) then
             N(x,y,z)%rock_colour_b = rock_colour_b
          end if
#endif
       endif
      end do
    end do
  end do
end subroutine lbe_init_rock_colour
#endif

!> Modify geometries to add simple walls and tunnels.
subroutine lbe_modify_rock(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N

  ! Values to set the 3 sides to
  integer :: xval, yval, zval, mval

  ! Dummy indices
  integer :: i,j,k
  integer :: ti, tj, tk

  ! Values to set on the edges:
  !  1: rock
  !  0: no rock
  ! -1: do not modify
  ! The mapping to boundary_cond values are chosen to keep backwards compatibility
  select case(boundary_cond)
  case(1)
    call log_msg("  Adding rock cage.",.false.)
    xval = 1
    yval = 1
    zval = 1
  case(2)
    call log_msg("  Adding square duct in z-direction with cleared sites.",.false.)
    xval = 1
    yval = 1
    zval = 0
  case(3)
    call log_msg("  Clearing x and y sides.",.false.)
    xval = 0
    yval = 0
    zval = -1
  case(4)
    call log_msg("  Adding square duct in z-direction.",.false.)
    xval = 1
    yval = 1
    zval = -1
  case(5)
    call log_msg("  Clearing z sides.",.false.)
    xval = -1
    yval = -1
    zval = 0
  case(6)
    call log_msg("  Adding walls on x sides.",.false.)
    xval = 1
    yval = -1
    zval = -1
  case default
    return
  end select

  ! Loop over the local domain, and set relevant rock states
  ! Rock wins over no rock wins over NOOP
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ! Positions on the global lattice
        ti = i + ccoords(1) * nx
        tj = j + ccoords(2) * ny
        tk = k + ccoords(3) * nz
        
        if ( ( ti .le. boundary_width ) .or. ( ti .gt. tnx - boundary_width ) ) then

          if ( ( tj .le. boundary_width ) .or. ( tj .gt. tny - boundary_width ) ) then
            
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = max(xval, max(yval, zval) )
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            else
              mval = max(xval, yval)
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            endif

          else
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = max(xval, zval)
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            else
              mval = xval
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            endif
          endif
          
        else
          if ( ( tj .le. boundary_width ) .or. ( tj .gt. tny - boundary_width ) ) then

            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval= max(yval, zval)
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            else
              mval = yval
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            endif

          else
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = zval
              if (mval .ge. 0) N(i,j,k)%rock_state = mval
            else
              ! NOOP
            endif
          endif

        endif

      enddo
    enddo
  enddo

end subroutine lbe_modify_rock

end module lb3d_init_rock_module

