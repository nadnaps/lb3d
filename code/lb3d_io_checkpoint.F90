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

!> checkpointing for the lb part
!>
!> md checkpointing happens in \c lbe_md_input_module and \c
!> lbe_md_output_module.
module lb3d_io_checkpoint_module
  use lb3d_global_module
  use lb3d_io_arrstats_module, only: dump_arrstats_checkpoint&
       &,restore_arrstats_checkpoint
  use lb3d_config_module
  use lb3d_helper_module, only:unixtime
  use lb3d_io_helper_module
  use lb3d_bc_leesedwards_module!,only: shear_sum ! Need this to write the shear_sum
  
#ifdef USEXDRF
 use lb3d_io_xdrf_module
#endif

#ifdef USEHDF
 use lb3d_io_hdf5_module
#endif

  use lb3d_lattice_module
  use lb3d_mpi_parameters_module
  use lb3d_mpi_module!, only:calculate_displacements
  use lb3d_log_module!, only: log_msg

!use lb3d_mpi_parameters_module,only:myrankc
  implicit none

contains

!> \{
!> \name     Wrapper functions to dump, restore, delete checkpoints

!>  This is a wrapper for the different checkpointing routines
subroutine dump_checkpoint(N,chk_handle)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: chk_handle
  character(len=40)          :: tag
  integer :: err

  CALL log_msg("Checkpointing...",.false.)

  if ( .not. checkpoint_safe ) then
    if ( num_chkp_files .gt. 0 ) CALL delete_checkpoint()
  endif

  CALL log_msg("  Writing malleable checkpoint files...",.false.)
  CALL dump_lattice_MP(N)

  call dump_arrstats_checkpoint

  ! Put an MPI barrier here so the checkparams file only gets written
  ! after all the checkpoint files are done.
  CALL MPI_Barrier(Comm_cart,err)

  if ( checkpoint_format .eq. 'xdr' ) then
#ifdef USEXDRF
    if(myrankc .eq. 0) then
       CALL log_msg("  Writing checkparams file...",.false.)
       CALL checkpoint_params_xdr()
    endif
#else
    CALL log_msg("Cannot write checkparams file without USEXDRF set.",.false.)
#endif
  endif

  if ( checkpoint_safe ) then
    if ( num_chkp_files .gt. 0 ) CALL delete_checkpoint()
  endif

  CALL log_msg("Done checkpointing.",.false.)

end subroutine dump_checkpoint

!>  This is a wrapper for the different restore routines
subroutine restore_checkpoint(rnt, chkid, N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer                :: rnt, chkid, new_chk_uid, ierror
  character(len=128)     :: msgstr

  if (myrankc == 0) then
    n_restore = rnt
    chk_uid = chkid
    ! Will need a new UID once we've restarted from this 
    ! checkpoint so generate it and bcast it with other data
    CALL system_clock(new_chk_uid)

    write(msgstr,"('Restoring simulation from timestep ',I0,' with format ',A)") rnt, trim(checkpoint_format)
    CALL log_msg(trim(msgstr),.false.)

  endif

  CALL MPI_Bcast(n_restore,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(chk_uid,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(new_chk_uid,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL restore_lattice_MP(N)
  call restore_arrstats_checkpoint

  CALL log_msg("  Succesfully restored lattice.",.false.)
  ! Now we've restarted, change to new UID
  chk_uid = new_chk_uid
  write(msgstr,"('  New UID set to <',I10.10,'>.')") chk_uid
  CALL log_msg(trim(msgstr),.false.)
end subroutine restore_checkpoint

!>  This routine deletes old checkpoint files. 
subroutine delete_checkpoint()
  implicit none
  character(len=256) :: filename
  character(len=512) :: msgstr
  integer            :: last

  if (num_chkp_files == 0) return

  last = nt - num_chkp_files * n_checkpoint
  if (last .lt. 0) then
    return
  else
    write(msgstr,"('  Deleting checkpoint at t = ',I10.10,'...')") last
    CALL log_msg(trim(msgstr),.false.)

    if (checkpoint_format .eq. 'hdf') then
    ! For CHECKPOINT_HDF5 we don't have per-core files, so rank 0 can take care of everything.
      if (myrankc == 0) then
        CALL lbe_make_filename_cp(filename,'cp_n_r','.h5',last)
        CALL lbe_delete_file(filename)
        CALL lbe_make_filename_cp(filename,'cp_n_b','.h5',last)
        CALL lbe_delete_file(filename)
        CALL lbe_make_filename_cp(filename,'cp_n_s','.h5',last)
        CALL lbe_delete_file(filename)
        CALL lbe_make_filename_cp(filename,'cp_d__','.h5',last)
        CALL lbe_delete_file(filename)
        CALL lbe_make_filename_cp(filename,'cp_r_s','.h5',last)
        CALL lbe_delete_file(filename)
      endif
    else if (checkpoint_format .eq. 'xdr') then
      CALL lbe_make_filename_cp_rank(filename,'checkpoint' ,'.xdr', last, myrankc)
      CALL lbe_delete_file(filename)
      if (myrankc == 0) then
        CALL lbe_make_filename_cp(filename,'checkparams','.xdr',last)
        CALL lbe_delete_file(filename)
        CALL lbe_make_filename_cp(filename,'checktopo','.xdr',last)
        CALL lbe_delete_file(filename)
      endif
    else
      CALL lbe_make_filename_cp_rank(filename,'checkpoint','.bin',last,myrankc)
      CALL lbe_delete_file(filename)
    endif

  endif

end subroutine delete_checkpoint
!> \}

!> This dumps out the entire state of the system in binary format,
!> not including the halo region. Intended for checkpoint/restarts.
subroutine dump_lattice_bin(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi
  integer :: file_unit = 10

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  CALL lbe_make_filename_cp_rank(filename,'checkpoint','.bin',nt,myrankc)

  open(file_unit,file=filename,form='unformatted')

  do x=1,nxi
    do y=1,nyi
      do z=1,nzi
#ifdef SINGLEFLUID
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%rock_state
#else
#ifdef NOSURFACTANT
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%rock_state
#else
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%n_s(:), N(x,y,z)%da(:), N(x,y,z)%db(:), N(x,y,z)%rock_state
#endif
#endif
      end do
    end do
  end do
  close(file_unit)
end subroutine dump_lattice_bin

!>  Restores a binary lattice state, saved from a previous run.
subroutine restore_lattice_bin(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi
  integer :: file_unit = 10

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  nt = n_restore ! Let's do time warp, yeah...
  CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrankc)

  open(file_unit,file=filename,form='unformatted')

  do x=1,nxi
    do y=1,nyi
      do z=1,nzi
#ifdef SINGLEFLUID
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%rock_state
#else
#ifdef NOSURFACTANT
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%rock_state
#else
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%n_s(:), N(x,y,z)%da(:), N(x,y,z)%db(:), N(x,y,z)%rock_state
#endif
#endif
      end do
    end do
  end do
  close(file_unit)
end subroutine restore_lattice_bin

!> This dumps out the entire state of the system in binary, XDR, or
!> HDF5 format, not including the halo region. Intended for variable
!> number of processors checkpoint/restarts.
!>
!> \param[in] N local chunk of the system
subroutine dump_lattice_MP(N)
  implicit none
  type(lbe_site),intent(in),dimension(0:,0:,0:) :: N

  select case(checkpoint_format)
    case( 'ori' )
      ! old case
      CALL dump_lattice_bin(N)
    case( 'bin' )
      ! malleable numprocs, bin output
      CALL dump_lattice_MPbin(N)
    case( 'xdr' )
      ! malleable numprocs, XDR output
#ifdef USEXDRF
      call dump_lattice_MPrxdr(N)
#else
      call log_msg("Cannot write rXDRF checkpoint without USEXDRF set.",.false.)
#endif
    case( 'hdf' )
      ! malleable numprocs, HDF5 output
#ifdef USEHDF
      CALL dump_lattice_MPhdf5(N)
#else
      CALL log_msg("Cannot write HDF5 checkpoint without USEHDF set.",.false.)
#endif
  end select
end subroutine dump_lattice_MP

!>  Restores a binary, XDR or HDF5 lattice state, saved from a previous run.
!>  Intended for variable processors checkpoint/restarts.
subroutine restore_lattice_MP(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  ! For timing purposes
  ! real(kind=rk):: restore_start, restore_stop
  ! restore_start = MPI_Wtime()
  select case(checkpoint_format)
    case( 'ori' )
      ! old case
      CALL restore_lattice_bin(N)
    case( 'bin' )
      nt = n_restore
      ! malleable numprocs, bin input
      CALL restore_lattice_MPbin(N)
    case( 'xdr' )
#ifdef USEXDRF
      nt = n_restore
      ! malleable numprocs, XDR input
      call restore_lattice_MPrxdr(N)
#else
      call log_msg("Cannot restore rXDR checkpoint without USEXDRF set."&
           &,.false.)
#endif
    case( 'hdf' )
#ifdef USEHDF
      nt = n_restore
      ! malleable numprocs, HDF5 input
      CALL restore_lattice_MPhdf5(N)
#else
      CALL log_msg("Cannot restore HDF5 checkpoint without USEHDF set.",.false.)
#endif

  end select
  ! For timing purposes
  ! restore_stop =  MPI_Wtime()
  ! print*,'RESTORE++','myrank is ',myrankc,' checkpoint time = ', &
  ! restore_stop - restore_start
end subroutine restore_lattice_MP

!> This dumps out the entire state of the system in binary format,
!> not including the halo region. Intended for checkpoint/restarts.
!> Processor grid information is also dumped at the top the file for
!> processor 0.
subroutine dump_lattice_MPbin(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi
  integer :: file_unit = 10

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  CALL lbe_make_filename_cp_rank(filename,'checkpoint','.bin',nt,myrankc)

  open(file_unit,file=filename,form='unformatted')

  if (myrankc == 0) then
    write(file_unit) cdims(1), cdims(2), cdims(3)
  end if

  do x=1,nxi
    do y=1,nyi
      do z=1,nzi
#ifdef SINGLEFLUID
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%rock_state
#else
#ifdef NOSURFACTANT
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%rock_state
#else
        write(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%n_s(:), N(x,y,z)%da(:), N(x,y,z)%db(:), N(x,y,z)%rock_state
#endif
#endif
      end do
    end do
  end do
  close(file_unit)
end subroutine dump_lattice_MPbin

!> Restores a binary lattice state, saved from a previous run.
!> This includes the case where checkpoint was performed on fewer processors,
!> more processors and same number of processors.
subroutine restore_lattice_MPbin(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: cdims_old(3)
  integer :: Px, Py, Pz, ierror
  logical :: decrease_procs, equal_procs
  integer :: file_unit = 10

  if (myrankc == 0) then
    CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrankc)
    open(file_unit,file=filename,form='unformatted')
    read(file_unit) Px, Py, Pz  ! Checkpointed processor grid
    close(file_unit)

    call check_changed_decomposition((/Px,Py,Pz/))
  end if

  CALL MPI_Bcast(Px,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(Py,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(Pz,1,MPI_INTEGER,0,comm_cart,ierror)

  cdims_old(1) = Px
  cdims_old(2) = Py
  cdims_old(3) = Pz

  decrease_procs = ( Px > cdims(1) .or. Py > cdims(2) .or. Pz > cdims(3) )
  if (decrease_procs) then
    CALL decrease_restore_lattice_MPbin(N, cdims_old)
  else
    equal_procs = ( Px == cdims(1) .and. Py == cdims(2) .and. Pz == cdims(3) )
    if (equal_procs) then
      CALL equal_restore_lattice_MPbin(N)
    else
      CALL increase_restore_lattice_MPbin(N, cdims_old)
    end if
  end if
end subroutine restore_lattice_MPbin

!> Restores a binary lattice state, saved from a previous run.
!> This covers the case where the dump was performed on fewer processors.
subroutine increase_restore_lattice_MPbin(N, cdims_old)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: cdims_old(:)
  CALL par_restore_lattice_MPbin(N, cdims_old)
end subroutine increase_restore_lattice_MPbin

!> Restores a binary lattice state, saved from a previous run.
!> This covers the case where the dump was performed on fewer processors.
subroutine par_restore_lattice_MPbin(N, cdims_old)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: cdims_old(:)

  real(kind=rk), allocatable, dimension(:,:) :: tempBuffer ! to store inputs for MPI
  type(lbe_site) :: tempN ! to store inputs from file
  character(len=256) :: filename
  real(kind=rk) :: nred,nblue,nsurf
  integer :: x,y,z,s,nxi,nyi,nzi
  integer :: nx_old, ny_old, nz_old, myrank_old
  integer :: Px, Py, Pz, Rx, Ry, Rz,Tx, Ty, Tz, myrank_current
  integer :: i, j, k, ii, jj, kk
  integer :: target, source, ierror, stat(MPI_STATUS_SIZE)
  integer :: num_sent, num_own, num_received
  integer, dimension(3) :: garbage
  logical :: controlling_processor
  integer :: ibuf
  integer :: file_unit = 10

#ifdef SINGLEFLUID
  ibuf = 20
#else
#ifdef NOSURFACTANT
  ibuf = 39
#else
  ibuf = 64
#endif
#endif

  nt = n_restore
  num_sent = 0
  num_own = 0
  num_received = 0

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  Px = cdims_old(1)
  Py = cdims_old(2)
  Pz = cdims_old(3)

  Rx = cdims(1) / Px
  Ry = cdims(2) / Py
  Rz = cdims(3) / Pz

  ! Determine what processors will open checkpoint files
  if (myrankc == 0) then
    myrank_old = 0
    CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrank_old)
    open(file_unit,file=filename,form='unformatted')
    read(file_unit) garbage(1),garbage(2),garbage(3)  ! Checkpointed processor grid

    ii = 0
    jj = 0
    kk = 0

    i  = 0
    j  = 0
    k  = 0

    controlling_processor = .true.
    nx_old = nxi * Rx
    ny_old = nyi * Ry
    nz_old = nzi * Rz
  else
    myrank_current = myrankc

    ii = myrankc / (cdims(3)*cdims(2))
    jj = (myrankc - cdims(3)*cdims(2)*ii)/cdims(3)
    kk = mod(myrankc,cdims(3))

    controlling_processor = ( mod(ii,Rx) == 0 .and. mod(jj,Ry) == 0 .and. mod(kk,Rz) == 0 )
    if (controlling_processor) then
      i = ii / Rx
      j = jj / Ry
      k = kk / Rz
      myrank_old = k + (Pz)*j + Pz*Py*i

      CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrank_old)
      open(file_unit,file=filename,form='unformatted')
      nx_old = nxi * Rx
      ny_old = nyi * Ry
      nz_old = nzi * Rz
    end if
  end if

  !   Preliminaries are over with, now get down to reading files and
  !   distributing the data!!

  allocate(tempBuffer(nzi,ibuf))
  if (controlling_processor) then
    !   i, j and k were defined for controlling_processors earlier
    do Tx=Rx*i,Rx*(i+1)-1
      do x = 1,nxi
        do Ty=Ry*j,Ry*(j+1)-1
          do y = 1,nyi
            do Tz=(Rz*k),Rz*(k+1)-1
              do z= 1,nzi
#ifdef SINGLEFLUID
                read(file_unit) tempN%n_r(:), tempN%rock_state
#else
#ifdef NOSURFACTANT
                read(file_unit) tempN%n_r(:), tempN%n_b(:), tempN%rock_state
#else
                read(file_unit) tempN%n_r(:), tempN%n_b(:), tempN%n_s(:), tempN%da(:), tempN%db(:), tempN%rock_state
#endif
#endif

                tempBuffer(z,1:19)  = tempN%n_r(:)
#ifndef SINGLEFLUID
                tempBuffer(z,20:38) = tempN%n_b(:)
#endif
#ifndef NOSURFACTANT
                tempBuffer(z,39:57) = tempN%n_s(:)
                tempBuffer(z,58:60) = tempN%da(:)
                tempBuffer(z,61:63) = tempN%db(:)
#endif
                tempBuffer(z,ibuf)  = tempN%rock_state
              end do
              target = Tz+cdims(3)*Ty+cdims(3)*cdims(2)*Tx
              if (target == myrankc) then
                do z = 1, nzi
                  N(x,y,z)%n_r(:) = tempBuffer(z,1:19)
#ifndef SINGLEFLUID
                  N(x,y,z)%n_b(:) = tempBuffer(z,20:38)
#ifndef NOSURFACTANT
                  N(x,y,z)%n_s(:) = tempBuffer(z,39:57)
                  N(x,y,z)%da(:)   = tempBuffer(z,58:60)
                  N(x,y,z)%db(:)   = tempBuffer(z,61:63)
#endif
#endif
                  N(x,y,z)%rock_state = tempBuffer(z,ibuf)
                end do
              else
                CALL MPI_Send(  tempBuffer,         & ! buf
                                nzi*ibuf,           & ! length
                                LBE_REAL,           & ! datatype
                                target,             & ! dest
                                x+nxi*y,            & ! tag
                                Comm_Cart,          & ! communicator
                                ierror)
              end if
            end do ! for Tz
          end do   ! for y
        end do     ! for Ty
      end do       ! for x
    end do         ! for Tx
    close(file_unit)
  else ! just a receiving processor
    source = kk - mod(kk,Rz) + cdims(3)*(jj-mod(jj,Ry)) + cdims(3)*cdims(2)*(ii-mod(ii,Rx))
    ! write(*,*) myrankc,' receives from ', source, nxi*nyi
    do x = 1, nxi
      do y = 1, nyi
        CALL MPI_Recv(  tempBuffer,                  & ! buf
                        nzi*ibuf,                    & ! length
                        LBE_REAL,                    & ! datatype
                        source,                      & ! source
                        x+nxi*y,                     & ! tag
                        Comm_Cart,                   & ! communicator
                        stat,                        & ! status
                        ierror)
        do z = 1, nzi
          N(x,y,z)%n_r(:) = tempBuffer(z, 1:19)
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = tempBuffer(z,20:38)
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = tempBuffer(z,39:57)
          N(x,y,z)%da(:)  = tempBuffer(z,58:60)
          N(x,y,z)%db(:)  = tempBuffer(z,61:63)
#endif
#endif
          N(x,y,z)%rock_state = tempBuffer(z,ibuf)
        end do
      end do
    end do
  end if
  deallocate(tempBuffer)
end subroutine par_restore_lattice_MPbin

!>  Restores a binary lattice state, saved from a previous run.
subroutine equal_restore_lattice_MPbin(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi
  integer :: cdims_old(3)
  integer :: file_unit = 10

  nt = n_restore
  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrankc)

  open(file_unit,file=filename,form='unformatted')
  if (myrankc == 0) then
    ! read the three integers giving the old dimensions again
    ! They're not needed anymore, but they have to be gotten out of the way
    read(file_unit) cdims_old(1), cdims_old(2), cdims_old(3)
  end if

  do x=1,nxi
    do y=1,nyi
      do z=1,nzi
#ifdef SINGLEFLUID
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%rock_state
#else
#ifdef NOSURFACTANT
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%rock_state
#else
        read(file_unit) N(x,y,z)%n_r(:), N(x,y,z)%n_b(:), N(x,y,z)%n_s(:), N(x,y,z)%da(:), N(x,y,z)%db(:), N(x,y,z)%rock_state
#endif
#endif
      end do
    end do
  end do
  close(file_unit)
end subroutine equal_restore_lattice_MPbin

!> Restores a binary lattice state, saved from a previous run.
!> This covers the case where restore is performed on fewer processors.
subroutine decrease_restore_lattice_MPbin(N, cdims_old)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: cdims_old(:)

  character(len=256) :: filename
  integer :: file_unit = 10

  integer :: x,y,z,s,nxi,nyi,nzi
  integer :: nx_old, ny_old, nz_old, myrank_old
  integer :: Px, Py, Pz, Rx, Ry, Rz,Tx, Ty, Tz, myrank_current
  integer :: i, j, k, ii, jj, kk

  nt = n_restore
  ! nxi,nyi,nzi hold size of our own chunk of lattice without halo
  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  !Px,Py,Pz hold size of old CPU grid
  Px = cdims_old(1)
  Py = cdims_old(2)
  Pz = cdims_old(3)

  !Ratio x,y,z = ratio of old size to new.
  Rx = Px / cdims(1)
  Ry = Py / cdims(2)
  Rz = Pz / cdims(3)

   !nx_old contains X-dimension of single-CPU chunk from previous run.
  nx_old = nxi / Rx
  ny_old = nyi / Ry
  nz_old = nzi / Rz

  ! cdims contains current CPU grid dimensions; calculate the position
  ! of the CPU in the grid. This makes assumptions about how
  ! MPI orders the CPUs.
  ii = myrankc / (cdims(3)*cdims(2))
  jj = (myrankc - cdims(3)*cdims(2)*ii)/cdims(3)
  kk = mod(myrankc,cdims(3))
  i = Rx*ii
  j = Ry*jj
  k = Rz*kk

  do Tx = 0, Rx-1
    do Ty = 0, Ry-1
      do Tz = 0, Rz-1
        myrank_old = k+Tz + Pz*(j+Ty) + Pz*Py*(i+Tx)
        CALL lbe_make_filename_restore_rank(filename,'checkpoint','.bin',myrank_old)
        open(file_unit,file=filename,form='unformatted')
        if (myrankc == 0 .and.  Tx == 0 .and.  Ty ==0 .and.  Tz == 0) then
          ! read the three integers giving the old dimensions again, so
          ! all files now point to the first element of N_old to read.
          read(file_unit) (cdims_old(s), s=1,3)
        end if
        do x = 1, nx_old
          do y = 1, ny_old
            do z = 1, nz_old
#ifdef SINGLEFLUID
              read(file_unit) &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_r(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%rock_state
#else
#ifdef NOSURFACTANT
              read(file_unit) &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_r(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_b(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%rock_state
#else
              read(file_unit) &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_r(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_b(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%n_s(:),    &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%da(:),      &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%db(:),      &
                  N(x+Tx*nx_old,y+Ty*ny_old,z+Tz*nz_old)%rock_state
#endif
#endif
            end do
          end do
        end do
        close(file_unit)
      end do
    end do
  end do
end subroutine decrease_restore_lattice_MPbin
! end ML

#ifdef USEXDRF

!> This dumps out the entire state of the system in rXDR format,
!> not including the halo region.
!>
!> \param[in] N local chunk of the system
!>
!> Besides the lattice data, the root process creates an additional
!> file containing topology information.
subroutine dump_lattice_MPrxdr(N)
  type(lbe_site),intent(in),dimension(0:,0:,0:) :: N
  character(len=256) :: filename
  integer :: file_id,ierror,stat
  integer,allocatable :: all_ccoords(:)

  ! These would be function pointers in C. In fact, xdrf(double|float)
  ! ares external functions, but fortran does not seem to care about
  ! the type anyway...
  integer,external :: xdrfint


  prepare_receive_coords: if (myrankc==0) then
     allocate (all_ccoords(1:3*nprocs),stat=stat)
     call check_allocate(stat,'dump_lattice_MPrxdr(): all_ccoords')
  end if prepare_receive_coords

  ! not sure whether allgather would be faster on BG/P...
  call MPI_Gather(ccoords,3,MPI_INTEGER,all_ccoords,3,MPI_INTEGER,0,comm_cart&
       &,ierror)
  call checkmpi(ierror,'dump_lattice_MPrxdr(): MPI_Gather() failed')

  write_topo: if (myrankc==0) then
     call lbe_make_filename_cp(filename,'checktopo','.xdr',nt)
     call xdrfopen(file_id,filename,"w",ierror)
     call check_xdrfopen(ierror,filename)

     call xdrfvector(file_id,cdims,3,xdrfint,ierror)
     call xdrfvector(file_id,all_ccoords,3*nprocs,xdrfint,ierror)
     if (checkpoint_shearsum) then
       CALL xdrfdouble(file_id,shear_sum,ierror)
     end if
     call xdrfclose(file_id,ierror)
  end if write_topo

  if (myrankc==0) deallocate (all_ccoords)

  call dump_lattice_MPxdr(N,.false.)
end subroutine dump_lattice_MPrxdr

!> This dumps out the entire state of the system in XDR format, not
!> including the halo region. Intended for variable processor
!> checkpoint/restarts.
!>
!> \param[in] N local lattice chunk
!>
!> \param[in] pt write topology-related information in front of the
!> data (optional; defaults to \c .true. on \c myrankc==0 and \c
!> .false. elsewhere)
subroutine dump_lattice_MPxdr(N,pt)
  type(lbe_site),intent(in),dimension(0:,0:,0:) :: N
  logical,intent(in),optional :: pt
  logical prepend_topo
  character(len=256) :: filename
  real(kind=rk) :: nred,nblue,nsurf
  integer :: x,y,z,s,nxi,nyi,nzi
  integer :: size_n_r, size_n_b, size_n_s, size_da, size_db
  integer :: file_id, ierr

  if (present(pt)) then
     prepend_topo = pt
  else
     if (myrankc==0) then
        prepend_topo = .true.
     else
        prepend_topo = .false.
     end if
  end if

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  CALL lbe_make_filename_cp_rank(filename,'checkpoint','.xdr',nt,myrankc)
  CALL xdrfopen(file_id,filename,"w",ierr)

  if( ierr .eq. 0 ) then
    CALL log_msg("dump_lattice_MPxdr: xdrfopen returned error",.false.)
    CALL check_xdrfsloppy()
  endif

  if (prepend_topo) then
    do x = 1, size(cdims)
      CALL xdrfint(file_id,cdims(x),ierr)
    end do
    if (checkpoint_shearsum) then
      CALL xdrfdouble(file_id,shear_sum,ierr)
    end if
  end if
  !
  ! n_r, n_b, n_s and d are all  real 64 bit arrays
  ! Their output could be handled by a single call if we could assume
  ! they were contiguous in storage (which they probably are...)
  ! rock_state is real(kind=rk) now!
  !
  size_n_r = size(N(1,1,1)%n_r)
#ifndef SINGLEFLUID
  size_n_b = size(N(1,1,1)%n_b)
#endif
#ifndef NOSURFACTANT
  size_n_s = size(N(1,1,1)%n_s)
  size_da  = size(N(1,1,1)%da  )
  size_db  = size(N(1,1,1)%db  )
#endif

  do x=1,nxi
    do y=1,nyi
      do z=1,nzi
        do s = 1, size_n_r
          CALL xdrfdouble(file_id,N(x,y,z)%n_r(s)  ,ierr)
        end do
#ifndef SINGLEFLUID
        do s = 1, size_n_b
          CALL xdrfdouble(file_id,N(x,y,z)%n_b(s)  ,ierr)
        end do
#ifndef NOSURFACTANT
        do s = 1, size_n_s
          CALL xdrfdouble(file_id,N(x,y,z)%n_s(s)  ,ierr)
        end do
        do s = 1, size_da
          CALL xdrfdouble(file_id,N(x,y,z)%da(s)    ,ierr)
        end do
        do s = 1, size_db
          CALL xdrfdouble(file_id,N(x,y,z)%db(s)    ,ierr)
        end do
#endif
#endif
        CALL xdrfdouble(file_id,N(x,y,z)%rock_state,ierr)
      end do
    end do
  end do
  CALL xdrfclose(file_id,ierr)
end subroutine dump_lattice_MPxdr

!> Restores an rXDR lattice state, saved from a previous run.
!>
!> \param[in,out] N local lattice chunk
subroutine restore_lattice_MPrxdr(N)
  type(lbe_site),intent(inout),dimension(0:,0:,0:) :: N
  integer,allocatable :: all_files(:),all_files_start(:,:),counts(:),displs(:)&
       &,my_files(:),my_old_start(:,:),old_ccoords(:,:),old_ranks(:,:,:)&
       &,all_start(:,:)
  character(len=256) :: filename
  integer fbounds(2,3),file_id,i,ierror,lpos(3),max_oc(3),min_oc(3)&
       &,n_files_total,n_read,old_cdims(3),old_cs(3),old_nprocs,p,pp,stat,x,y,z
  character(len=512) :: msgstr

  ! These would be function pointers in C. In fact, xdrf(double|float)
  ! ares external functions, but fortran does not seem to care about
  ! the type anyway...
  integer,external :: xdrfint

  if (myrankc==0) then
     call lbe_make_filename_restore(filename,'checktopo','.xdr')
     call xdrfopen(file_id,filename,"r",ierror)
     call check_xdrfopen(ierror,filename)

     call xdrfvector(file_id,old_cdims,3,xdrfint,ierror)
     call check_xdrfint(ierror)

     old_nprocs = product(old_cdims)

     write(msgstr&
          &,"('  Previous decomposition: ',2(I0,' x '),I0,' ( ',I0,' files)')")&
          & old_cdims,old_nprocs
     call log_msg(trim(msgstr),.false.)

     allocate (old_ccoords(3,0:old_nprocs-1)&
          &,old_ranks(0:old_cdims(1)-1,0:old_cdims(2)-1,0:old_cdims(3)-1)&
          &,stat=stat)
     call check_allocate(stat,'restore_lattice_MPrxdr(): old_ccoords,old_ranks')

     do p=0,old_nprocs-1
        call xdrfvector(file_id,old_ccoords(1:3,p),3,xdrfint,ierror)
     end do

     if (checkpoint_shearsum) then
       call xdrfdouble(file_id,shear_sum,ierror)
     endif 
     
     call xdrfclose(file_id,ierror)

     ! old_ranks maps old process coordinates to old myrankc numbers
     do p=0,old_nprocs-1
        old_ranks(old_ccoords(1,p),old_ccoords(2,p),old_ccoords(3,p)) = p
     end do

     ! fortran's column-major order make each all_start(:,p)
     ! contiguous in memory (as needed for MPI_Gather())
     allocate (all_start(3,0:nprocs-1),counts(0:nprocs-1),displs(0:nprocs-1)&
          &,stat=stat)
     call check_allocate(stat&
          &,'restore_lattice_MPrxdr(): all_start,counts,displs')
  end if

  ! not sure whether allgather would be faster on BG/P...
  call MPI_Gather(start,3,MPI_INTEGER,all_start,3,MPI_INTEGER,0&
       &,comm_cart,ierror)
  call checkmpi(ierror,'restore_lattice_MPrxdr(): MPI_Gather() failed')

  if (myrankc==0) then
     old_cs = (/tnx,tny,tnz/)/old_cdims ! old chunksize (total size stays same)

     ! 1st pass: count number of files each process has to read
     do p=0,nprocs-1
        ! min and max process coordinates within old coordinates
        min_oc = (all_start(:,p)-1)/old_cs
        max_oc = (all_start(:,p)+(/nx,ny,nz/)-2)/old_cs

        counts(p) = product(1+max_oc-min_oc)
     end do
  end if

  call MPI_Scatter(counts,1,MPI_INTEGER,n_read,1,MPI_INTEGER,0,comm_cart&
       &,ierror)
  call checkmpi(ierror,'restore_lattice_MPrxdr(): MPI_Scatter() failed')

  allocate (my_files(n_read),my_old_start(3,n_read),stat=stat)
  call check_allocate(stat,'restore_lattice_MPrxdr(): my_files,my_old_start')

  if (myrankc==0) then
     n_files_total = sum(counts)

     write(msgstr,"('  Total number of file accesses: ',I0)") n_files_total
     call log_msg(trim(msgstr),.false.)

     allocate (all_files(n_files_total),all_files_start(3,n_files_total)&
          &,stat=stat)
     call check_allocate(stat&
          &,'restore_lattice_MPrxdr(): all_files,all_files_start')

     ! 2nd pass: collect and scatter list of files to read for each process
     pp = 0
     do p=0,nprocs-1
        ! min and max process coordinates within old coordinates
        min_oc = (all_start(:,p)-1)/old_cs
        max_oc = (all_start(:,p)+(/nx,ny,nz/)-2)/old_cs

        do x=min_oc(1),max_oc(1)
           do y=min_oc(2),max_oc(2)
              do z=min_oc(3),max_oc(3)
                 pp = pp+1
                 all_files(pp) = old_ranks(x,y,z)
                 all_files_start(:,pp) = (/x,y,z/)*old_cs+1
              end do
           end do
        end do
     end do

     call calculate_displacements(counts,displs)
  end if

  call MPI_Scatterv(all_files,counts,displs,MPI_INTEGER&
       &,my_files,n_read,MPI_INTEGER,0,comm_cart,ierror)
  call checkmpi(ierror,'restore_lattice_MPrxdr(): MPI_Scatterv() failed')

  if (myrankc==0) then
     ! scattering all_files_start means exactly 3 times more integers
     counts = counts*3
     displs = displs*3
  end if

  call MPI_Scatterv(all_files_start,counts,displs,MPI_INTEGER&
       &,my_old_start,n_read*3,MPI_INTEGER,0,comm_cart,ierror)
  call checkmpi(ierror,'restore_lattice_MPrxdr(): MPI_Scatterv() failed')

  call MPI_Bcast(old_cs,3,MPI_INTEGER,0,comm_cart,ierror)

  files: do i=1,n_read
     call lbe_make_filename_restore_rank(filename,'checkpoint','.xdr',my_files(i))

     ! which part of the file to read?
     fbounds(1,:) = max(1,1+start-my_old_start(:,i))
     fbounds(2,:) = min(old_cs,start+(/nx,ny,nz/)-my_old_start(:,i))

     ! and where to put it within N?
     lpos = 1+my_old_start(:,i)-start

     call restore_lattice_chunk_from_file_rxdr(N,filename,old_cs,fbounds,lpos)
  end do files

  if (myrankc==0) &
       &deallocate (old_ccoords,old_ranks,all_start,counts,all_files,displs)
  deallocate (my_files,my_old_start)
end subroutine restore_lattice_MPrxdr

!> restores a part of the local lattice chunk from a part of an XDR
!> checkpoint file
!>
!> \param[in,out] N local lattice chunk (including halo of depth 1)
!>
!> \param[in] filename XDR file to read from
!>
!> \param[in] fdim dimensions of the whole chunk in the file in all
!> three directions
!>
!> \param[in] fbounds first and last site to read in each direction
!> within \c fdim (1st index: first/last; 2nd index: x,y,z; possible
!> values: \c (1:fdim))
!>
!> \param[in] lpos local lattice position corresponding to position \c
!> (/1,1,1/) in the file. In general this is not the first position
!> read. \c lpos might well be negative.
subroutine restore_lattice_chunk_from_file_rxdr(N,filename,fdim,fbounds,lpos)
    type(lbe_site),intent(inout),dimension(0:,0:,0:) :: N
    character(len=256),intent(in) :: filename
    integer,intent(in) :: fdim(3),fbounds(2,3),lpos(3)
    integer file_id,ierror,l(3),x,y,z
    type(lbe_site) :: dummy

    call xdrfopen(file_id,filename,"r",ierror)
    call check_xdrfopen(ierror,filename)

    read_loops: do x=1,fdim(1)
       do y=1,fdim(2)
          do z=1,fdim(3)
             if (all((/x,y,z/)>=fbounds(1,:)).and.&
                  &all((/x,y,z/)<=fbounds(2,:))) then

                l = (/x,y,z/)+lpos-1 ! local position in N
                call lbe_xdrflbe_site(file_id,N(l(1),l(2),l(3)))

                ! skip the whole rest if the last site was just read
                if (all((/x,y,z/)==fbounds(2,:))) exit read_loops
             else
                call lbe_xdrflbe_site(file_id,dummy)
             end if
          end do
       end do
    end do read_loops

    call xdrfclose(file_id,ierror)
end subroutine restore_lattice_chunk_from_file_rxdr

!> reads or writes a lattice site from/to an XDR file
!>
!> \param[inout] file_id XDR file identifier as obtained from \c xdrfopen()
!>
!> \param[inout] site lattice site to read or write
!>
!> Whether data is read or written is determined by how \c file_id was
!> opened. In this sense, \c lbe_xdrflbe_site() works the same way as
!> the xdrf routines for builtin types (for instance \c xdrfdouble()).
subroutine lbe_xdrflbe_site(file_id,site)
    integer,intent(inout) :: file_id
    type(lbe_site),intent(inout) :: site
    integer ierror

    ! This would be a function pointer in C. In fact, xdrfdouble is an
    ! external functions, but fortran does not seem to care about the
    ! type anyway...
    integer,external :: xdrfdouble

    call xdrfvector(file_id,site%n_r,nvecs,xdrfdouble,ierror)
#ifndef SINGLEFLUID
    call xdrfvector(file_id,site%n_b,nvecs,xdrfdouble,ierror)
#ifndef NOSURFACTANT
    call xdrfvector(file_id,site%n_s,nvecs,xdrfdouble,ierror)
    call xdrfvector(file_id,site%da,3,xdrfdouble,ierror)
    call xdrfvector(file_id,site%db,3,xdrfdouble,ierror)
#endif
#endif
    ! this is stupid but having defined xdrfdouble above we cannot
    ! call it directly anymore
    call xdrfvector(file_id,site%rock_state,1,xdrfdouble,ierror)
end subroutine lbe_xdrflbe_site

! END USEXDRF
#endif

!> \{
!> \name CHECKPOINTING FACILITY Support Subroutines

!> Malleable checkpoint restoration does not work yet for equal number
!> of processes but different decomposition. This subroutine checks
!> both conditions and aborts if necessary.  ccd(1:3)  should contain
!> the decomposition read from the checkpoint.
subroutine check_changed_decomposition(ccd)
    integer,intent(in) :: ccd(3)
    character(len=256) :: msgstr

    if (product(ccd)==nprocs.and.any(ccd/=cdims)) then
       write(msgstr,"('Current decomposition does not match checkpoint. "&
            &//"Restart simulation with cdx = ',I0,', cdy = ',I0,', "&
            &//"and cdz = ',I0,'!')") ccd
       call log_msg(msgstr,.false.)
       call abend
    end if
end subroutine check_changed_decomposition
!> \}

#ifdef USEHDF

!> This dumps out the entire state of the system in HDF5 format,
!> not including the halo region. Intended for variable processor 
!> checkpoint/restarts.
subroutine dump_lattice_MPhdf5(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N

  real(kind=rk), dimension(:,:,:,:), allocatable :: v

  integer nxi, nyi, nzi, ierror, x, y, z

  character(len=256) :: msgstr
  character(len=256) :: filename

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

#ifndef HDF5_FLIP
  allocate(v(19,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,19),stat=ierror)
#endif
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate vector buffer for population checkpoints",.true.)
    return
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(:,x,y,z) = N(x,y,z)%n_r(:)
#else
        v(x,y,z,:) = N(x,y,z)%n_r(:)
#endif
      enddo
    enddo
  enddo
  CALL lbe_make_filename_cp(filename, 'cp_n_r', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)

#ifndef SINGLEFLUID
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(:,x,y,z) = N(x,y,z)%n_b(:)
#else
        v(x,y,z,:) = N(x,y,z)%n_b(:)
#endif
      enddo
    enddo
  enddo

  CALL lbe_make_filename_cp(filename, 'cp_n_b', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)
#ifndef NOSURFACTANT

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(:,x,y,z) = N(x,y,z)%n_s(:)
#else
        v(x,y,z,:) = N(x,y,z)%n_s(:)
#endif
      enddo
    enddo
  enddo
  CALL lbe_make_filename_cp(filename, 'cp_n_s', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)
#endif

#endif

  deallocate(v)


#ifndef NOSURFACTANT

#ifndef HDF5_FLIP
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,3),stat=ierror)
#endif
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate vector buffer for dipole checkpoint",.true.)
    return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(:,x,y,z) = N(x,y,z)%da(:)
#else
        v(x,y,z,:) = N(x,y,z)%da(:)
#endif
      enddo
    enddo
  enddo
  CALL lbe_make_filename_cp(filename, 'cp_da__', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)
  deallocate(v)


#ifndef HDF5_FLIP
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,3),stat=ierror)
#endif
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate vector buffer for dipole checkpoint",.true.)
    return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(:,x,y,z) = N(x,y,z)%db(:)
#else
        v(x,y,z,:) = N(x,y,z)%db(:)
#endif
      enddo
    enddo
  enddo
  CALL lbe_make_filename_cp(filename, 'cp_db__', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)
  deallocate(v)

#endif

#ifndef HDF5_FLIP
  allocate(v(1,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,1),stat=ierror)
#endif
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate vector buffer for rock checkpoint",.true.)
    return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        v(1,x,y,z) = N(x,y,z)%rock_state
#else
        v(x,y,z,1) = N(x,y,z)%rock_state
#endif
      enddo
    enddo
  enddo
  CALL lbe_make_filename_cp(filename, 'cp_r_s', '.h5', nt)
  CALL dump_vector_phdf5(v,filename)

  deallocate(v)

  ! Rank 0 writes additional parameters to the cp_n_r... file (which should always be there)

  if ( myrankc == 0) then
    write(msgstr,"('  Wrote parameter shear_sum = ',F16.10)") shear_sum
    CALL log_msg(trim(msgstr),.false.)

    CALL lbe_make_filename_cp(filename, 'cp_n_r', '.h5', nt)
    CALL lbe_write_cp_param_phdf5(shear_sum,filename,'shear_sum')
  endif

  !> \todo Add shear_sum parameter
end subroutine dump_lattice_MPhdf5

!> Restores an HDF5 lattice state, saved from a previous run.
!> This includes the case where checkpoint was performed on fewer processors,
!> more processors and same number of processors.
subroutine restore_lattice_MPhdf5(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N

  real(kind=rk), dimension(:,:,:,:), allocatable :: v

  integer :: nxi, nyi, nzi, ierror, x, y, z

  integer :: start_time, end_time

  character(len=256) :: msgstr
  character(len=256) :: filename

  nt = n_restore ! Let's do time warp, yeah...

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  start_time = unixtime()

#ifndef HDF5_FLIP
  allocate(v(19,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,19),stat=ierror)
#endif

  CALL lbe_make_filename_restore(filename, 'cp_n_r', '.h5')
  CALL read_vector_phdf5(v,filename,19)

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%n_r(:) = v(:,x,y,z) 
#else
        N(x,y,z)%n_r(:) = v(x,y,z,:) 
#endif
      enddo
    enddo
  enddo

#ifndef SINGLEFLUID
  CALL lbe_make_filename_restore(filename, 'cp_n_b', '.h5')
  CALL read_vector_phdf5(v,filename,19)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%n_b(:) = v(:,x,y,z) 
#else
        N(x,y,z)%n_b(:) = v(x,y,z,:) 
#endif
      enddo
    enddo
  enddo

#ifndef NOSURFACTANT
  CALL lbe_make_filename_restore(filename, 'cp_n_s', '.h5')
  CALL read_vector_phdf5(v,filename,19)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%n_s(:) = v(:,x,y,z) 
#else
        N(x,y,z)%n_s(:) = v(x,y,z,:) 
#endif
      enddo
    enddo
  enddo
#endif

#endif

  deallocate(v)

#ifndef NOSURFACTANT

#ifndef HDF5_FLIP
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,3),stat=ierror)
#endif

  CALL lbe_make_filename_restore(filename, 'cp_da__', '.h5')
  CALL read_vector_phdf5(v,filename,3)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%da(:) = v(:,x,y,z) 
#else
        N(x,y,z)%da(:) = v(x,y,z,:) 
#endif
      enddo
    enddo
  enddo

  deallocate(v)

#ifndef HDF5_FLIP
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,3),stat=ierror)
#endif

  CALL lbe_make_filename_restore(filename, 'cp_db__', '.h5')
  CALL read_vector_phdf5(v,filename,3)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%db(:) = v(:,x,y,z) 
#else
        N(x,y,z)%db(:) = v(x,y,z,:) 
#endif
      enddo
    enddo
  enddo

  deallocate(v)

#endif

#ifndef HDF5_FLIP
  allocate(v(1,1:nxi,1:nyi,1:nzi),stat=ierror)
#else
  allocate(v(1:nxi,1:nyi,1:nzi,1),stat=ierror)
#endif

  CALL lbe_make_filename_restore(filename, 'cp_r_s', '.h5')
  CALL read_vector_phdf5(v,filename,1)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
#ifndef HDF5_FLIP
        N(x,y,z)%rock_state = v(1,x,y,z) 
#else
        N(x,y,z)%rock_state = v(x,y,z,1) 
#endif
      enddo
    enddo
  enddo

  deallocate(v)

  ! Rank 0 writes additional parameters to the cp_n_r... file (which should always be there)
  if ( myrankc == 0) then

    CALL lbe_make_filename_restore(filename, 'cp_n_r', '.h5')
    CALL lbe_read_cp_param_phdf5(shear_sum,filename,'shear_sum')
    write(msgstr,"('  Read parameter shear_sum = ',F16.10)") shear_sum
    CALL log_msg(trim(msgstr),.false.)


  endif

  ! Now broadcast all parameters
  CALL MPI_Bcast(shear_sum,1,MPI_REAL8,0,comm_cart,ierror)

  if (myrankc == 0) then
    CALL lbe_make_metadata_hdf5()
  endif !myrankc == 0

  CALL MPI_Barrier(comm_cart,ierror)

  end_time = unixtime()

  write(msgstr,"('Checkpoint restore took ',I0,' seconds.')") (end_time-start_time)
  call log_msg(trim(msgstr),.false.)

end subroutine restore_lattice_MPhdf5

#endif

end module lb3d_io_checkpoint_module
