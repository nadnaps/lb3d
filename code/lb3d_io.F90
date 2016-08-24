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

!> Contains input/output related functions.
!> \todo This should be a mere interface.
module lb3d_io_module

  use lb3d_global_module

  use lb3d_lattice_module

  use lb3d_io_helper_module,only:nflags,check_dump_now,lbe_setup_fluxz
  use lb3d_io_arrstats_module,only:lb3d_io_setup_arrstats
  use lb3d_io_stress_module,only:setup_stress 

#ifdef USEHDF
  use lb3d_io_hdf5_module,only:lb3d_io_init_hdf5
#endif

  use lb3d_io_checkpoint_module,only: dump_checkpoint
  use lb3d_io_dump_data_module, only: dump_data

  use lb3d_config_module

  use lb3d_mpi_module,only:abend,recv_final_lattice,send_final_lattice
  use lb3d_mpi_parameters_module,only:input_source,tnx,tny,tnz,myrankc

  use lb3d_timer_module, only: start_timer, stop_timer,ti_dump
  use lb3d_log_module, only: log_msg, log_msg_ws, log_ws

  implicit none

contains

!> Setup of serial i/o
!      ==================================================================
  subroutine lb3d_io_init (stage)
!      ==================================================================

    integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_io_init stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif

    select case (stage)

    case (0) ! pre mpi init

    ! Get version and revision and log
    !FIXME move these constructs to lb3d_log?
    write(msgstr,"('Starting ',A,'.')") trim(lbeversion)
    CALL log_msg(trim(msgstr),.false.)
    CALL log_msg(trim(lbeflags),.false.)

    ! Parse command-line arguments
    CALL lbe_parse_arguments()

    ! Set up the input-file handles
    CALL lbe_define_inp_file()

    ! initialize fluxz region names before they are read from input file
    CALL lbe_setup_fluxz()

    case (2) ! post mpi init / pre memory alloc

    ! Do not call this before MPI has been initialized
    !FIXME Really depends on InitializeMPIcomms? 
    ! Maybe call it twice for verbosity as of what is set by MPI
    CALL lbe_detect_flags()

    call lb3d_io_halo_extent()

#ifdef USEHDF

    ! Set up HDF5 if required.
    CALL lb3d_io_init_hdf5()

#endif
    call lb3d_io_setup_arrstats()

    call setup_stress()

    case default

    end select    
    
  end subroutine lb3d_io_init

!> Reads any arguments passed on the command-line, and prints them out.
subroutine lbe_parse_arguments()
  implicit none
  integer, external  :: iargc
  integer            :: i, argc
  character(len=128) :: argbuf, msgstr

  CALL log_msg_ws("----------( Parsing CL arguments )-----------",.false.)
  argc = iargc() ! Number of arguments on command line.
  i = 1
  if (argc == 0) then
    CALL log_msg("No command line arguments supplied",.false.)
  endif
  do
    if (i .gt. argc) exit
    CALL getarg(i, argbuf)

    ! -f <input-file name>
    if ("-f" == argbuf) then
      if (i == argc) then
        CALL log_msg("Need filename after -f, ignoring",.false.)
        call Abend
      end if
      i = i+1
      CALL getarg(i,argbuf)
      arg_input_file = trim(argbuf)
      arg_input_file_p = 1
      write(msgstr,"('Found -f ',A)") arg_input_file
      CALL log_msg(trim(msgstr),.false.)
    end if

    ! -d <diff input-file name>
    if ("-d" == argbuf) then
      if (i == argc) then
        CALL log_msg("Need filename after -d, ignoring",.false.)
        call Abend
      end if
      i = i+1
      CALL getarg(i,argbuf)
      arg_input_dfile = trim(argbuf)
      arg_input_dfile_p = 1
      write(msgstr,"('Found -d ',A)") arg_input_dfile
      CALL log_msg(trim(msgstr),.false.)
    end if

    ! -r <restore-string>
    if ("-r" == argbuf) then
      if (i == argc) then
        CALL log_msg("Need restore-string after -r, ignoring",.false.)
        call Abend
      end if
      i = i+1
      CALL getarg(i,argbuf)
      arg_restore_string = trim(argbuf)
      arg_restore_string_p = 1
      write(msgstr,"('Found -r ',A)") arg_restore_string
      CALL log_msg(trim(msgstr),.false.)
    end if
    i = i+1
  end do
end subroutine lbe_parse_arguments


!> Sets the \c inp_file variable to contain the name of the input file.
!> If the -f option is passed on the command line, then
!> this is used. Otherwise, .input-file is checked -- if it contains the
!> word "INTERACTIVE", then the input file name is read from stdin.
!> Otherwise, the file named in .input-file is read.
subroutine lbe_define_inp_file()
  logical :: inp_files_p
  ! assume only called by rank 0
  if (0 .ne. arg_input_file_p) then
    ! Take it from command line
    inp_file = arg_input_file
  else
    inquire(FILE = input_source, EXIST = inp_files_p)
    if(inp_files_p) then
      open(UNIT = 10, FILE = input_source, STATUS = 'UNKNOWN')
      read(UNIT = 10, FMT = '(A)') inp_file
      close(UNIT = 10)
      if(index(inp_file,'INTERACTIVE').gt.0) then
        print *,'Input file ?'
        read(*,*) inp_file
      end if
    else
      CALL log_msg("FATAL ERROR: Could not open .input-file. Aborting...",.false.)
      CALL Abend
    end if
  endif
  inp_file = trim(inp_file)
  CALL log_msg("Variable inp_file = <"//trim(inp_file)//">.",.false.)
end subroutine lbe_define_inp_file



  subroutine lb3d_io_write_data
    call start_timer(ti_dump)
    if (nt >= n_sci_start) then
       if (post) then
          CALL postprocess(N)
       else
          CALL dump_data(N)
       endif
      ! CALL dump_parallel(N,whole_N)
    endif
    call stop_timer(ti_dump)
  end subroutine lb3d_io_write_data

  subroutine lb3d_io_write_checkpoint
    if (mod(nt,n_checkpoint)==0) then
       call start_timer(ti_dump)
       CALL dump_checkpoint(N,0)
       call stop_timer(ti_dump)
    end if
  end subroutine lb3d_io_write_checkpoint


!> write info on which compiler flags were used to stdout
!>
!> 2010-05-03: Added by Stefan, use the 'find_flags.sh' script to easily
!> check if this is still up to date!
subroutine lbe_detect_flags()
  nflags = 0 ! Declared in lb3d_io_helper so HDF5 metadata output can also use this value
  CALL log_msg_ws("--------( Reporting compiler flags )---------",.false.)

#ifdef BOUNCEBACK
  CALL log_msg("  BOUNCEBACK",.false.)
  nflags = nflags + 1
#endif
#ifdef BUGGYIFORT11
  CALL log_msg("  BUGGYIFORT11",.false.)
  nflags = nflags + 1
#endif
#ifdef BUGGYSENDINCOLLECT
  CALL log_msg("  BUGGYSENDINCOLLECT",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_CHECKPOINT
  CALL log_msg("  DEBUG_CHECKPOINT",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_HDF5
  CALL log_msg("  DEBUG_HDF5",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_HDF5_TIMING
  CALL log_msg("  DEBUG_HDF5_TIMING",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_LE
  CALL log_msg("  DEBUG_LE",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_MPI
  CALL log_msg("  DEBUG_MPI",.false.)
  nflags = nflags + 1
#endif
#ifdef DEBUG_REPORTMDCOMM
  CALL log_msg("  DEBUG_REPORTMDCOMM",.false.)
  nflags = nflags + 1
#endif
#ifdef DIST
  CALL log_msg("  DIST",.false.)
  nflags = nflags + 1
#endif
#ifdef ELEC
  CALL log_msg("  ELEC",.false.)
  nflags = nflags + 1
#endif
#ifdef FASTBDIST2
  CALL log_msg("  FASTBDIST2",.false.)
  nflags = nflags + 1
#endif
#ifdef HDF5_FLIP
  CALL log_msg("  HDF5_FLIP",.false.)
  nflags = nflags + 1
#endif
#ifdef LB3D
  CALL log_msg("  LB3D",.false.)
  nflags = nflags + 1
#endif
#ifdef MCMP
  CALL log_msg("  MCMP",.false.)
  nflags = nflags + 1
#endif
#ifdef MD
  CALL log_msg("  MD",.false.)
  nflags = nflags + 1
#endif
#ifdef MPI_ALLGV_FASTER_THAN_GV
  CALL log_msg("  MPI_ALLGV_FASTER_THAN_GV",.false.)
  nflags = nflags + 1
#endif
#ifdef NOCALLSYSTEM
  CALL log_msg("  NOCALLSYSTEM",.false.)
  nflags = nflags + 1
#endif
#ifdef NOEDGESTEP
  CALL log_msg("  NOEDGESTEP",.false.)
  nflags = nflags + 1
#endif
#ifdef NOIEEEARITHMETIC
  CALL log_msg("  NOIEEEARITHMETIC",.false.)
  nflags = nflags + 1
#endif
#ifdef NOISNAN
  CALL log_msg("  NOISNAN",.false.)
  nflags = nflags + 1
#endif
#ifdef NOSURFACTANT
  CALL log_msg("  NOSURFACTANT",.false.)
  nflags = nflags + 1
#endif
#ifdef OLDRRFORCE
  CALL log_msg("  OLDRRFORCE",.false.)
  nflags = nflags + 1
#endif
#ifdef PARTICLESTRESS
  CALL log_msg("  PARTICLESTRESS",.false.)
  nflags = nflags + 1
#endif
#ifdef RELTIME
  CALL log_msg("  RELTIME",.false.)
  nflags = nflags + 1
#endif
#ifdef RWALK
  CALL log_msg("  RWALK",.false.)
  nflags = nflags + 1
#endif
#ifdef SINGLEFLUID
  CALL log_msg("  SINGLEFLUID",.false.)
  nflags = nflags + 1
#endif
#ifdef USEHDF
  CALL log_msg("  USEHDF",.false.)
  nflags = nflags + 1
#endif
#ifdef USE_HDF5_INDEPENDENT_IO
  CALL log_msg("  USE_HDF5_INDEPENDENT_IO",.false.)
  nflags = nflags + 1
#endif
#ifdef USE_IBM_LARGEBLOCK_IO
  CALL log_msg("  USE_IBM_LARGEBLOCK_IO",.false.)
  nflags = nflags + 1
#endif
#ifdef USENEWMASSC
  CALL log_msg("  USENEWMASSC",.false.)
  nflags = nflags + 1
#endif
#ifdef USEOLDROCK
  CALL log_msg("  USEOLDROCK",.false.)
  nflags = nflags + 1
#endif
#ifdef USEXDRF
  CALL log_msg("  USEXDRF",.false.)
  nflags = nflags + 1
#endif
#ifdef VECTORIZE
  CALL log_msg("  VECTORIZE",.false.)
  nflags = nflags + 1
#endif
#ifdef WALLCONST
  CALL log_msg("  WALLCONST",.false.)
  nflags = nflags + 1
#endif
#ifdef XDRROCKWET
  CALL log_msg("  XDRROCKWET",.false.)
  nflags = nflags + 1
#endif
  CALL log_ws(.false.)
end subroutine lbe_detect_flags

!>This function collects the data from all processors and dumps the output
!>into a single file.
!>
!>Added 14.06.02 by Jens
subroutine postprocess(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  type(lbe_site), dimension(:,:,:),allocatable :: Nm
  integer :: ierror,x,y,z

   if (    ( check_dump_now(sci_int, n_sci_int) ) &
       .or.( check_dump_now(sci_sur, n_sci_sur) ) &
       .or.( check_dump_now(sci_od , n_sci_od ) ) &
       .or.( check_dump_now(sci_wd , n_sci_wd ) ) &
       .or.( check_dump_now(sci_dir, n_sci_dir) ) &
       .or.( check_dump_now(sci_vel, n_sci_vel) ) &
       .or.( check_dump_now(sci_flo, n_sci_flo) ) &
       .or.( check_dump_now(sci_arrows, n_sci_arrows) ) &
       .or.( check_dump_now(sci_velocities, n_sci_velocities) ) &
       .or.( check_dump_now(sci_rock, n_sci_rock) ) &
       .or.( check_dump_now(sci_pressure, n_sci_pressure) ) &
       .or.( check_dump_now(sci_fluxz, n_sci_fluxz) ) &
       .or.( check_dump_now(sci_massfluxz, n_sci_massfluxz) ) &
       .or.( check_dump_now(sci_perm, n_sci_perm) ) &
       ) then

    CALL log_msg("Started postprocessing:",.false.)

    if(index(dump_format,'hdf').gt.0)then
      CALL log_msg("Dumping HDF5 data:",.false.)
#ifdef USEHDF
      CALL dump_data(N)
#else
      CALL log_msg("FATAL ERROR: HDF5 support is switched off but HDF5 output is requested. Aborting...",.false.)
      CALL Abend
#endif
      CALL log_msg('Finished dumping HDF5 data.',.false.)
    else ! not phdf5
      CALL log_msg("Postprocessing, no HDF5",.false.)
      if (myrankc == 0) then
        ! Allocate the whole arena in the master CPU
        allocate(Nm(0:tnx+1,0:tny+1,0:tnz+1),stat=ierror)
        if (ierror .ne. 0) then
          CALL log_msg("FATAL ERROR: Unable to allocate memory for postprocessing data buffer. Aborting...",.true.)
          CALL Abend
        end if
        ! Collect the data
        CALL recv_final_lattice(N,Nm)
        ! Finally, dump data to disk.
        CALL log_msg("Dumping data.",.false.)
        CALL dump_data(Nm)

        deallocate(Nm)
      else    ! myrankc != 0
        ! Do nothing than sending my chunk to the master
        CALL send_final_lattice(N)
      end if  !myrankc
    end if !phdf5
  end if
end subroutine postprocess

    !> increase \c halo_extent if required by some of the IO routines
    subroutine lb3d_io_halo_extent
        if (sci_profile) then
           ! this is necessary because of
           ! fluid_velocity_and_density_and_site_occupation() which is
           ! called by dump_profile() and in turn calls
           ! md_calculate_forces()
           halo_extent = max(halo_extent,2)
        end if
        
    end subroutine lb3d_io_halo_extent

end module lb3d_io_module
