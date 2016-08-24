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

!> Contains function dependend on mpi
!> \todo More consistent naming for the mpi part
module lb3d_mpi_module

    use lb3d_config_module, only: post, dbg_report_topology, nt
    use lb3d_global_module,only: rk, msgstr!, tsize, border, chunksize,
    use lb3d_lattice_module

    use lb3d_mpi_parameters_module
    use lb3d_log_module
    use lb3d_helper_module, only:makedatestr

    use lb3d_config_module, only: chk_uid,dump_double,dump_format,folder,cpfolder&
       &,gr_out_file,nt,restore_string,srccpfolder

    implicit none

    !> Number of flags set - counted in lbe_detect_flags() and used for
    !> HDF metadata
    !integer :: nflags

    !> \{
    !> \name not in namelist
    character(len=32), save :: restore_fmt_t = '("t",i8.8,"-",i10.10)'
    character(len=32), save :: restore_fmt_p = '("p",i8.8,"-",i10.10)'
    !> \}
    
    type halo
       sequence
       integer, dimension(3) :: ls ! Lower Send
       integer, dimension(3) :: us ! Upper Send
       integer, dimension(3) :: lr ! Lower Recv
       integer, dimension(3) :: ur ! Upper Recv
       integer               :: mpitype ! Base MPI type
       integer               :: extent  ! Halo extent
    end type halo

    
    !> sizes of total system
    real(kind=rk),save :: tsize(3)
    !> sizes of my part
    real(kind=rk),save :: chunksize(3)
    !> lo/hi boundaries of my box in each dimension
    real(kind=rk),save :: border(2,3)  

  contains

    !> Calls subroutines to setup MPI parallelisation
    subroutine lb3d_mpi_init (stage)


      integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_mpi_init stage',I0)") stage
    CALL log_msg(trim(msgstr),.false.) 
#endif

      select case (stage)

      case (0)
         ! Initialise MPI and set up the topology.
         CALL InitializeMPIcomms()

      case (1) ! mpi init

         ! Change from MPI_COMM_WORLD to a Cartesian grid
         CALL ReorderMPIcomms()

         ! Divide the system into blocks
         CALL lbe_divide_grid()

         ! Set up the datatypes.
         CALL lbe_parallel_init()
           
         ! Write out the processor topology
         if ((myrankc == 0).and.(.not.(post))) then
            CALL lbe_write_topology()
          !  CALL lb3d_log_msg("Wrote topology",.false.)
         endif
         
      case (4)
         call setup_fluid(whole_N)
      case default
#ifdef LB3D_DEBUG_INFO
         !FIXME log_msg
         call log_msg('nothing to do.',.false.)
#endif
      end select

    end subroutine lb3d_mpi_init

!> This routine fills in an array, \c pcoords, with the coordinates
!> of each processor in the Cartesian topology.
!> Note that the array is zero-based, so that pcoords(:,i) holds
!> the coords of the ith processor.
subroutine find_topology(pcoords)

	implicit none
	integer :: np,pcoords(1:,0:) ! Processor rank, and coords in topology
	integer	:: sx,sy,sz	! Coords of start of a proc subdomain
	integer :: ierror
        character(LEN=80) bla

	do np=0,nprocs-1
		call MPI_cart_coords(	Comm_cart,		&
					np,			&
					nd,			&
					pcoords(:,np),		&
					ierror)
                bla = 'failed to find processor coords'
		call checkmpi(ierror,bla)
	end do

end subroutine find_topology

!> logs MPI cartesian coordinates and neighbors
!>
!> \warning This is likely to fail for extreme core counts.
subroutine debug_report_ccoords
  integer i,ierror,recv_ccoords(nd),recv_nnprocs(2*nd),status(MPI_STATUS_SIZE)
  integer :: tag = 1 ! MPI tag
  character(len=256) :: msgstr

  ! Report on the coordinates.
  if ( myrankc .gt. 0 ) then
    CALL MPI_Send(ccoords, nd , MPI_INTEGER, 0, tag, Comm_cart, ierror)
  else
    call log_msg("debug_report_ccoords(): reporting ccoords for each rank."&
         &,.false.)
    write (msgstr&
         &,"('  Rank ',I6.6,' has coordinates (',i0, ',', i0, ',', i0, ').')") &
         &myrankc, ccoords(1), ccoords(2), ccoords(3)
    call log_msg(trim(msgstr),.false.)
    do i=1,nprocs-1
      call mpi_recv(recv_ccoords,nd,MPI_INTEGER,i,tag,Comm_cart,status,ierror)
      write (msgstr&
           &,"('  Rank ',I6.6,' has coordinates (',i0, ',', i0, ',', i0, ').')"&
           &) i, recv_ccoords(1), recv_ccoords(2), recv_ccoords(3)
      call log_msg(trim(msgstr),.false.)
    enddo
  endif

  ! Report on nearest neighbours
  if ( myrankc .gt. 0 ) then
    CALL MPI_Send(nnprocs, 2*nd , MPI_INTEGER, 0, tag, Comm_cart, ierror)
  else
    call log_msg('debug_report_ccoords(): reporting nearest neighbours for '&
         &//'each rank (cx-1 cx+1 cy-1 cy+1 cz-1 cz+1).',.false.)
    write (msgstr,"('  Rank ',I6.6,' has neighbours ',6(I6.6,:,' '),'.')") &
         &myrankc,nnprocs(1,1:2),nnprocs(2,1:2),nnprocs(3,1:2)
    CALL log_msg(trim(msgstr),.false.)
    do i=1, nprocs - 1
      call mpi_recv(recv_nnprocs,2*nd,MPI_INTEGER,i,tag,Comm_cart,status,ierror)
      write (msgstr,"('  Rank ',I6.6,' has neighbours ',6(I6.6,:,' '),'.')") &
           &i,recv_nnprocs(1),recv_nnprocs(4),recv_nnprocs(2),recv_nnprocs(5)&
           &,recv_nnprocs(3),recv_nnprocs(6)
      CALL log_msg(trim(msgstr),.false.)
    enddo
  endif
end subroutine debug_report_ccoords

!> logs hostnames
!>
!> \warning This is likely to fail for extreme core counts.
subroutine debug_report_hostnames
  integer i,ierror,status(MPI_STATUS_SIZE)
  integer :: tag = 1 ! MPI tag
  character(len=MPI_MAX_PROCESSOR_NAME) :: recv_hname
  character(len=256) :: msgstr ! Hopefully this is large
                               ! enough. Should maybe calculate this
                               ! from MPI_MAX_PROCESSOR_NAME, but
                               ! Fortran is annoying with declarations

  if ( myrankc .gt. 0 ) then
    call mpi_send(hname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,tag,Comm_cart&
         &,ierror)
  else
    call log_msg("debug_report_hostnames(): reporting hostnames for each rank"&
         &,.false.)
    write (msgstr,"('  Rank ',I6.6,' has been assigned to host <',A,'>.')") &
         &myrankc,trim(hname)
    CALL log_msg(trim(msgstr),.false.)
    do i=1, nprocs-1
      call mpi_recv(recv_hname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,i,tag&
           &,Comm_cart,status,ierror)
      write (msgstr,"('  Rank ',I6.6,' has been assigned to host <',A,'>.')") &
           &i,trim(recv_hname)
      CALL log_msg(trim(msgstr),.false.)
    enddo
  endif
end subroutine debug_report_hostnames

!> fills the halo (extent 1) of a local chunk according to periodic
!> boundaries within itself
!>
!> \param[in,out] N local chunk of arbitrary size but expected to have
!> starting indices 0 and a halo of extent 1
subroutine lbe_halo_exchange_local(N)
    type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
    integer lx,ly,lz

    lx = size(N,1)-2
    ly = size(N,2)-2
    lz = size(N,3)-2

    N(     0,  1:ly,1:lz) = N(    lx,  1:ly,1:lz)
    N(  lx+1,  1:ly,1:lz) = N(     1,  1:ly,1:lz)
    N(0:lx+1,     0,1:lz) = N(0:lx+1,    ly,1:lz)
    N(0:lx+1,  ly+1,1:lz) = N(0:lx+1,     1,1:lz)
    N(0:lx+1,0:ly+1,   0) = N(0:lx+1,0:ly+1,  lz)
    N(0:lx+1,0:ly+1,lz+1) = N(0:lx+1,0:ly+1,   1)
end subroutine lbe_halo_exchange_local

    !> setup fluid coupling
    subroutine setup_fluid(N)
      type(lbe_site),intent(inout) :: &
           &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
      integer :: he,h1

      ! just abbreviations
      he = halo_extent
      h1 = halo_extent - 1

      ! setup custom mpi datatypes that simplify the halo exchange vastly

      call build_lbe_site_mpitype(lbe_site_mpitype)

      !! This type setup causes a full halo exchange, so  lbe_halo_exchange()
      !! may be replaced by  md_halo_exchange() .

      ! ! x swaps (still restricted y/z-range: data outside is out-dated)
      call build_lattice_chunk_mpitype&
           &(N,(/1,      he/),(/1,      ny/),(/1,      nz/),ls_mpitype(1))
      call build_lattice_chunk_mpitype&
           &(N,(/nx-h1,  nx/),(/1,      ny/),(/1,      nz/),us_mpitype(1))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1,     0/),(/1,      ny/),(/1,      nz/),lr_mpitype(1))
      call build_lattice_chunk_mpitype&
           &(N,(/nx+1,nx+he/),(/1,      ny/),(/1,      nz/),ur_mpitype(1))
      ! ! y swaps (full x-range - up-to-date data was received in x-swaps)
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/1,      he/),(/1,      nz/),ls_mpitype(2))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/ny-h1,  ny/),(/1,      nz/),us_mpitype(2))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/-h1,     0/),(/1,      nz/),lr_mpitype(2))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/ny+1,ny+he/),(/1,      nz/),ur_mpitype(2))

      ! ! z swaps (even full y-range - after z swap all data is complete)
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/1,      he/),ls_mpitype(3))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz-h1,  nz/),us_mpitype(3))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/-h1,     0/),lr_mpitype(3))
      call build_lattice_chunk_mpitype&
           &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz+1,nz+he/),ur_mpitype(3))

    end subroutine setup_fluid

    !> Build a custom mpi datatype that represents \c lbe_site, stores it in
    !> \c lsmt.
    subroutine build_lbe_site_mpitype(lsmt)
      ! include rock_colour_r and rock_colour_b in the hallo exchange may be innecessary,
      ! beacuse it only necessary to be done once at the beginning.
      integer,intent(out) :: lsmt ! type to be built
#ifdef SINGLEFLUID
      integer,parameter :: n_blocks = 6
#else
#ifdef NOSURFACTANT
      integer,parameter :: n_blocks = 10
#else
      integer,parameter :: n_blocks = 16
#endif
#endif
      integer lengths(n_blocks),types(n_blocks),nb
      integer(kind=MPI_ADDRESS_KIND) :: displs(n_blocks)
      integer(kind=MPI_ADDRESS_KIND) :: addrs(n_blocks)
      integer ierror
      type(lbe_site) :: sample(2)

      nb = 1

      lengths(nb) = 19         ! n_r
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%n_r(1),addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rho_r
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rho_r,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 3 ! u_r
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%u_r,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rock_colour_r
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rock_colour_r,addrs(nb),ierror)
      nb = nb + 1

#ifndef SINGLEFLUID

      lengths(nb) = 19         ! n_b
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%n_b(1),addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rho_b
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rho_b,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 3 ! u_b
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%u_b,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rock_colour_b
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rock_colour_b,addrs(nb),ierror)
      nb = nb + 1


#ifndef NOSURFACTANT
      lengths(nb) = 19         ! n_s
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%n_s(1),addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rho_s
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rho_s,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 3 ! u_s
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%u_s,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1 ! rock_colour_s
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rock_colour_s,addrs(nb),ierror)
      nb = nb + 1

!      lengths(nb) = 3          ! d
!      types(nb) = MPI_REAL8
!      call mpi_get_address(sample(1)%d(1),addrs(nb),ierror)
!      nb = nb + 1

      lengths(nb) = 3          ! da
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%da(1),addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 3          ! db
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%db(1),addrs(nb),ierror)
      nb = nb + 1


#endif
#endif

      lengths(nb) = 1 ! rock_state
      types(nb) = MPI_REAL8
      call mpi_get_address(sample(1)%rock_state,addrs(nb),ierror)
      nb = nb + 1

      lengths(nb) = 1   ! next lattice site
      types(nb) = MPI_UB
      call mpi_get_address(sample(2)%n_r(1),addrs(nb),ierror)

      displs(2:) = addrs(2:) - addrs(1)
      displs(1) = 0

      call mpi_type_create_struct(nb,lengths,displs,types,lsmt,ierror)
      call mpi_type_commit(lsmt,ierror)
    end subroutine build_lbe_site_mpitype

    !> Builds a custom mpi datatype that represents the part of  N  specified by
    !> the coordinate intervals  xr(2),yr(2),zr(2)  and stores it in  lcmt .
    subroutine build_lattice_chunk_mpitype(N,xr,yr,zr,lcmt)
      type(lbe_site),intent(in) :: &
           &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
      integer,intent(in) :: xr(2),yr(2),zr(2) ! chunk ranges for each dim.
      integer,intent(out) :: lcmt ! type to be built
      integer xrow_mt,xyplane_mt,xyzchunk_mt ! temporary mpi data types
      integer blocks,stride,ierror
      integer(kind=MPI_ADDRESS_KIND) adr1,adr2,base,offset
      integer lengths(1),displs(1) ! for  mpi_type_hindexed()

      ! mpi datatype for slices of N like  N(xr(1):xr(2),y,z)
      blocks = 1 + xr(2) - xr(1)
      call mpi_get_address(N(1,1,1),adr1,ierror)
      call mpi_get_address(N(2,1,1),adr2,ierror)
      stride = adr2 - adr1
      call mpi_type_hvector(blocks,1,stride,lbe_site_mpitype,xrow_mt,ierror)

      ! mpi datatype for slices of N like  N(xr(1):xr(2),yr(1):yr(2),z)
      blocks = 1 + yr(2) - yr(1)
      call mpi_get_address(N(1,1,1),adr1,ierror)
      call mpi_get_address(N(1,2,1),adr2,ierror)
      stride = adr2 - adr1
      call mpi_type_hvector(blocks,1,stride,xrow_mt,xyplane_mt,ierror)

      ! mpi datatype for whole chunk  N(xr(1):xr(2),yr(1):yr(2),zr(1):zr(2))
      blocks = 1 + zr(2) - zr(1)
      call mpi_get_address(N(1,1,1),adr1,ierror)
      call mpi_get_address(N(1,1,2),adr2,ierror)
      stride = adr2 - adr1
      call mpi_type_hvector(blocks,1,stride,xyplane_mt,xyzchunk_mt,ierror)

      ! position of the beginning of the chunk relative to the beginning of N
      call mpi_get_address(N,base,ierror)
      call mpi_get_address(N(xr(1),yr(1),zr(1)),offset,ierror)
      offset = offset - base

      ! lcmt becomes a datatype like  xyzchunk_mt  but relative to  base
      lengths = (/1/)
      displs = (/int(offset)/)
      call mpi_type_hindexed(1,lengths,displs,xyzchunk_mt,lcmt,ierror)
      call mpi_type_commit(lcmt,ierror)
    end subroutine build_lattice_chunk_mpitype


!===================================================================================
!>This routine was taken from ME3D and modified.
!>It initialises the MPI libraries, determines the rank
!>of each CPU, and sets up the CPU topology.
subroutine InitializeMPIcomms
  implicit none
  integer :: ierror, reslen ! reslen is just a dummy variable to be supplied as an argument
  character(len=128) :: msgstr

  ! Initialize MPI
  CALL MPI_Init(ierror)
  if (ierror .ne. MPI_SUCCESS) then
    CALL log_msg('MPI_Init() failed. Aborting...',.false.)
    call Abend
  endif

  ! Find out my rank
  CALL MPI_Comm_rank(MPI_COMM_WORLD, myrankw, ierror)
  if (myrankw == 0) then
    CALL log_msg_ws("-------( Initialized MPI_COMM_WORLD )--------",.false.)
  endif

  ! Find out the total number of processors.
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierror)
  if (myrankw == 0) then
    write(msgstr,"('Requested ',i0,' processors.')") nprocs
    CALL log_msg(trim(msgstr),.false.) 
  endif
  ! Added by Elena to be included in metadata
  CALL MPI_Get_processor_name(hname, reslen, ierror)

  ! Temporarily set Comm_Cart to world so we don't have to change the broadcast interface
  ! to exchange data before ReorderComms is called. Also myrankc = myrankw for now.
  Comm_Cart = MPI_COMM_WORLD
  myrankc = myrankw

end subroutine InitializeMPIcomms

subroutine ReorderMPIcomms
  integer :: ierror,i
  character(len=256) :: msgstr

  CALL log_msg_ws("-------( Creating MPI Cartesian grid )-------",.false.)
  ! Create the virtual topology.
  ! This call returns the dimensions of an appropriate
  ! Cartesian grid into which the processors may be divided.
  ! Would be nice to get a 'good' grid automatically
  CALL MPI_Dims_create(nprocs, nd, cdims, ierror)
  if (ierror .ne. MPI_SUCCESS) then
    CALL log_msg("MPI_Dims_create() failed. Aborting...",.false.)
    call Abend
  endif

  if (myrankw == 0) then
    CALL makedatestr(startsimul) ! For HDF5 output
    write (msgstr,"('Processors using a ',i0,' x ',i0,' x ',i0,' grid.')") cdims
    CALL log_msg(trim(msgstr),.false.)
  endif

  ! Now create a communicator to make it easy for each processor
  ! to talk only to its immediate neighbours.
  call MPI_Cart_create(MPI_COMM_WORLD& ! Make from all processors.
                     &,nd&             ! Dimension of lattice
                     &,cdims&          ! Dimensions of lattice
                     &,periodic_p&     ! Whether or not periodic
                     &,.true.&         ! Reorder ranks?
                     &,comm_cart&      ! New Cartesian comunicator.
                     &,ierror)

  if (ierror .ne. MPI_SUCCESS) then
    CALL log_msg("MPI_Cart_create() failed. Aborting...",.false.)
    call Abend
  endif

  ! Find out my new rank:- it may change if MPI_Cart_create()
  ! reorders the rankings so that the virtual topology corresponds
  ! more accurately to the physical topology.
  CALL MPI_Comm_rank(MPI_COMM_WORLD, myrankw, ierror)
  CALL MPI_Comm_rank(Comm_Cart, myrankc, ierror)

  ! Now fill in the ccoords array with my position in the Cartesian lattice.
  CALL MPI_Cart_get( Comm_Cart, nd, cdims, periodic_p, ccoords, ierror)

  ! Determine who my nearest neighbours are.
  do i = 1, nd
    CALL MPI_Cart_shift( Comm_Cart, i-1, 1, nnprocs(i,1), nnprocs(i,2), ierror )
  end do

  if ( dbg_report_topology ) then
    write(msgstr,"('DEBUG: comm_cart = ',I0)") comm_cart
    call log_msg(trim(msgstr),.false.)
    call debug_report_ccoords
    call debug_report_hostnames
  endif
end subroutine ReorderMPIcomms

!>This routine is called just before program termination - it
!>tells the MPI libraries to clean up and finish.
subroutine FinalizeMPIcomms
  implicit none
  integer :: ierror
  character(len=24) :: datestr

  CALL log_msg_ws("-------( Calling MPI_Finalize() )--------",.false.)
  CALL MPI_Finalize(ierror)
  
  end subroutine FinalizeMPIcomms


!>  Divides the tnx x tny x tnz grid into nx x ny x nz blocks assigned
!>  to the processors
subroutine lbe_divide_grid()
  character(len=128)     :: msgstr
  ! Must have the nx, ny, nz integer divisible by the
  ! number of processors assigned in that cartesian direction.

  if ( (MOD(nx,cdims(1))/= 0) .or. &
       (MOD(ny,cdims(2))/= 0) .or. &
       (MOD(nz,cdims(3))/= 0) ) then
    CALL log_msg("FATAL ERROR: Grid does not divide evenly by the number of processors. Aborting...",.false.)
    CALL Abend
  endif

  ! Remember the size of the original domain.
  tnx = nx; tny = ny; tnz = nz

  ! Now divide these digits so that a processor will only ever
  ! assign the amount of memory required for its own data.
  nx = nx/cdims(1)
  ny = ny/cdims(2)
  nz = nz/cdims(3)

  if(myrankc == 0 .and. myfarm == 0) then
    write(msgstr,"('Grid decomposition in blocks of size: ',i0,' x ',i0,' x ',i0,'.')") nx, ny, nz
    CALL log_msg(trim(msgstr),.false.)
  endif
  ! Find out where the local data lies in relation to the global data
  ! set. The array `start' will identify the start of the local data
  ! relative to the global data set.
  start(1) = ccoords(1)*nx + 1
  start(2) = ccoords(2)*ny + 1
  start(3) = ccoords(3)*nz + 1
end subroutine lbe_divide_grid


!> initialize stuff related to communication
subroutine lbe_parallel_init()
  tsize = real((/ tnx,tny,tnz /),kind=rk)
  chunksize = real((/ nx,ny,nz /),kind=rk)

  border(1,:) = start - 0.5_rk
  border(2,:) = start + chunksize - 0.5_rk

#ifdef MD
  ! Profile calculation requires the velocities of neighboring
  ! particles to be available. Of course, one could optimize
  ! this by introducing separate communication calls in
  !  dump_profile() .
  if (sci_profile) communicate_velocities = .true.
#endif
end subroutine lbe_parallel_init

!>Writes out an ASCII file detailing the position of each of
!>the processors in the virtual topology.
subroutine lbe_write_topology
  integer, dimension(:,:),allocatable :: pcoords
  character(len=256) :: filename
  integer :: i
  character(len=128)     :: msgstr

  call lbe_make_filename_output(filename,'coords','.txt',nt)

  write(msgstr,"('wtfile=<',A,'>')") trim(filename)
  CALL log_msg(trim(msgstr),.false.)

  open(unit=10,file=filename)
  allocate(pcoords(3,0:nprocs-1))
  CALL find_topology(pcoords)
  do i=nprocs-1,0,-1
    write(10,'(I6.6,a,I6.6,a,I6.6,a,I6.6)')	&
      i,' ',pcoords(1,i),' ',pcoords(2,i),' ',pcoords(3,i)
  end do
  deallocate(pcoords)
  close(unit=10)
end subroutine lbe_write_topology


!> Call to check the MPI return variable \c ierror.
!> Should probably be inlined.
!> If ierror is not MPI_SUCCESS, prints \c string and abends.
subroutine checkmpi(ierror,string)
  integer,intent(in)     :: ierror
  character*(*)         :: string
  character(len=128)     :: msgstr

  if (ierror .ne. MPI_SUCCESS) then
    write(msgstr,"('FATAL ERROR: MPI returned error:',A,'. Aborting...')") trim(string)
    CALL log_msg(trim(msgstr),.false.)
    CALL Abend
  end if
end subroutine checkmpi

!> print error message and exit
subroutine error(str)
    character*(*) str
    CALL log_msg("ERROR: "//trim(str),.true.)
    CALL Abend
end subroutine error


!> Abnormal ending
!>
!> This routine is called if something goes so disastrously wrong that
!> the simulation must be halted. At present, it's just a wrapper for
!> \c MPI_Abort(), but any required cleaning-up or checkpoint-on-error
!> code could be put here.
subroutine abend
    integer ierror

    call log_msg("Abnormal ending...",.true.)
    call MPI_Abort(MPI_COMM_WORLD,-1,ierror)
    stop
end subroutine abend


!> check if  stat , which should have been given to an  allocate -command
!> as the stat-argument, indicates an allocation failure; if yes, print
!>  msg  and stop the program.
subroutine check_allocate(stat,msg)
  integer,intent(in) :: stat
  character*(*),intent(in) :: msg

  if (stat/=0) then
     call log_msg('FATAL ERROR: allocate failed: '//msg,.true.)
     call Abend
  endif
end subroutine check_allocate


!> Find ave/min/max and histogram of one or more real(kind=rk) values
!> from all processors
subroutine stats_rk(data,ave,xmax,xmin,histo,histotmp,nhisto)
    real(kind=rk),intent(in) :: data
    real(kind=rk),intent(out) :: ave,xmax,xmin
    integer,intent(inout) :: histo(*),histotmp(*)
    integer,intent(in) :: nhisto
    integer j,ierror
    real(kind=rk) aave,del

    call mpi_allreduce(data,aave,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    ave = aave/nprocs
    call mpi_allreduce(data,xmax,1,MPI_REAL8,MPI_MAX,comm_cart,ierror)
    call mpi_allreduce(data,xmin,1,MPI_REAL8,MPI_MIN,comm_cart,ierror)

    histotmp(1:nhisto) = 0

    del = xmax-xmin
    if (del.eq.0.0) then
       j = 1
    else
       j = (data-xmin)/del * nhisto + 1
       if (j.gt.nhisto) j = nhisto
    endif
    histotmp(j) = histotmp(j) + 1

    call mpi_allreduce(histotmp,histo,nhisto,MPI_INTEGER,MPI_SUM,&
         &comm_cart,ierror)
end subroutine stats_rk

!> Take the state array N, and send it to rank 0.
!> Added 19.06.02 by Jens
subroutine send_final_lattice(N)
        type(lbe_site),dimension(0:,0:,0:) :: N
        integer :: ierror,stat,i,j,k
        character(LEN=80) :: bla
	real(kind=rk),dimension(:,:,:,:),allocatable :: buffer
        integer :: ibuf

#ifndef NOSURFACTANT
        ibuf = 64
#else
        ibuf = 39
#endif
#ifdef SINGLEFLUID
	ibuf = 20
#endif

	allocate(buffer(ibuf,nx,ny,nz))

                do k=1,nz
                 do j=1,ny
                  do i=1,nx
                   buffer(1:19,i,j,k) = N(i,j,k)%n_r(:)
#ifndef SINGLEFLUID
                   buffer(20:38,i,j,k) = N(i,j,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                   buffer(39:57,i,j,k) = N(i,j,k)%n_s(:)
                   buffer(58:60,i,j,k) = N(i,j,k)%da(:)
                   buffer(61:63,i,j,k) = N(i,j,k)%db(:)
#endif
                   buffer(ibuf,i,j,k)  = N(i,j,k)%rock_state
                  end do
                 end do
                end do

        call MPI_Send(  buffer,		      & ! buf
                        nx*ny*nz*ibuf,	      & ! length
                        LBE_REAL,             & ! datatype
                        0,                    & ! dest
                        tag_post,             & ! tag
                        Comm_Cart,       & ! communicator
                        ierror)
        bla = 'failed to send final data block'
        call checkmpi(ierror,bla)
	deallocate(buffer)

end subroutine send_final_lattice

!> Collects the whole lattice on \c myrankc==0, only to be called by rank 0
!>
!> \param[in] N local lattice (of rank 0)
!>
!> \param[out] Nm global lattice
!>
!> This should be called only by rankw 0, all other ranks should call
!> \c send_final_lattice(). Also rank 0's subdomain \c N will be
!> copied into \c Nm. Also the halo (extent 1) of Nm will be filled
!> accordingly.
subroutine recv_final_lattice(N,Nm)
    type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
    type(lbe_site),dimension(0:,0:,0:),intent(out) :: Nm
    real(kind=rk),dimension(:,:,:,:),allocatable :: buffer
    integer :: np,pcoords(nd) ! Processor rank, and coords in topology
    integer	:: sx,sy,sz	! Coords of start of a proc subdomain
    integer :: ierror,stat(MPI_STATUS_SIZE),i,j,k
    ! ibuf must have parameter attribute or its value changes after
    ! mpi_recv() . Obviously, ifort performs some optimizations and
    ! overlooks the side effects of mpi_recv() .
#ifdef SINGLEFLUID
    integer,parameter :: ibuf = 20
#else
#ifdef NOSURFACTANT
    integer,parameter :: ibuf = 39
#else
    integer,parameter :: ibuf = 64
#endif
#endif
    ! Same problem here: mpi_recv() / mpi_cart_coords() change their
    ! source/rank argument, so the loop counter must not be passed
    ! directly to these routines (Florian; ifort 9.1, mpich 1.0.4)
    integer tmpnp

    allocate(buffer(ibuf,nx,ny,nz))

    ! Fill big array with my own bits.
    sx=ccoords(1)*nx
    sy=ccoords(2)*ny
    sz=ccoords(3)*nz
    Nm(sx+1:sx+nx,sy+1:sy+ny,sz+1:sz+nz) = N(1:nx,1:ny,1:nz)

    ! collect other bits
    do np=nprocs-1,1,-1
       ! This works with ifort. Maybe outwitting other compilers needs
       ! something more tricky, like  tmpnp = nprocs-np .
       tmpnp = np
       call MPI_cart_coords(Comm_cart,tmpnp,nd,pcoords,ierror)
       call checkmpi(ierror,'failed to find processor coords')

       ! Get data from the appropriate processor.
       call MPI_Recv&
            &(buffer&                    ! buf
            &,nx*ny*nz*ibuf&             ! length
            &,LBE_REAL&                  ! datatype
            &,tmpnp&                     ! source
            &,tag_post&                  ! tag
            &,comm_cart&                 ! communicator
            &,stat&                      ! status
            &,ierror)
       call checkmpi(ierror,'failed to receive final data block')

       ! Fill big array with appropriate bits.
       sx=pcoords(1)*nx
       sy=pcoords(2)*ny
       sz=pcoords(3)*nz
       do k=1,nz
          do j=1,ny
             do i=1,nx
                Nm(sx+i,sy+j,sz+k)%n_r = buffer(1:19,i,j,k)
#ifndef SINGLEFLUID
                Nm(sx+i,sy+j,sz+k)%n_b = buffer(20:38,i,j,k)
#endif
#ifndef NOSURFACTANT
                Nm(sx+i,sy+j,sz+k)%n_s = buffer(39:57,i,j,k)
                Nm(sx+i,sy+j,sz+k)%da = buffer(58:60,i,j,k)
                Nm(sx+i,sy+j,sz+k)%db = buffer(61:63,i,j,k)
#endif
                Nm(sx+i,sy+j,sz+k)%rock_state = buffer(ibuf,i,j,k)
             end do
          end do
       end do
    end do
    deallocate(buffer)

    call lbe_halo_exchange_local(Nm)
end subroutine recv_final_lattice

!> Makes a filename by appending the given suffix and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_output(buffer, prefix, suffix, time)
  implicit none

  character(len=*)   :: buffer, prefix, suffix
  integer            :: time
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) time, chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'_',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine lbe_make_filename_output

!> generate displacement vector as needed by \c mpi_gatherv() or
!> \c mpi_scatterv() within \c comm_cart
!>
!> \param[in] counts number of elements for each rank
!>
!> \param[out] displs resulting displacement vector
subroutine calculate_displacements(counts,displs)
    integer,intent(in) :: counts(0:nprocs-1)
    integer,intent(out) :: displs(0:nprocs-1)
    integer i

    displs(0) = 0
    do i=0,nprocs-2
       displs(i+1) = displs(i) + counts(i)
    end do
end subroutine calculate_displacements


!> convert global to local coordinates
!>
!> \param[in] global global position vector
!>
!> \param[out] local local position vector
!>
!> \warning Thereby check for pbc, nevertheless result will be
!> worthless if \c global has no local representation but is in the
!> domain of another process (even if it is still in the local halo!).
subroutine local_coordinates(global,local)
  real(kind=rk),intent(in) :: global(3)
  real(kind=rk),intent(out) :: local(3)
  integer k

  ! If  global  is wrapped around due to pbc, wrap it back to get
  ! into the valid index range of  N(:,:,:) . Notice that this does not
  ! help if  global  is out of range for some other reason/bug!
  do k=1,3
    ! This is very tricky. Why was there halo_extent and -1.5 ? These
    ! values are surely wrong since they lead to a wrap-around already
    ! for global(1)==tnx if halo_extent==1 . This problem does not
    ! occur for the new code but I am not completely sure that this
    ! change will not break something else, maybe in connection with
    ! the stuff implemented by David. The reason for the 0.5 is that
    ! it is likely that particles are that much outside border because
    ! exchange() is called only on demand (see list_still_valid() !).
    ! FJ, 2009-09-03
!!$           if (global(k)>=border(2,k)+(real(halo_extent,kind=rk)-1.5_rk)) &
    if (global(k)>=border(2,k)+0.5_rk) then
      local(k) = global(k) - tsize(k)
!!$           else if (global(k)<border(1,k)-(real(halo_extent,kind=rk)-1.5_rk)) &
    else if (global(k)<border(1,k)-0.5_rk) then
      local(k) = global(k) + tsize(k)
    else
      local(k) = global(k)
    end if
  end do
  ! convert  local(:)  to local lbe coordinates
  local = local - (border(1,:)-0.5_rk)
end subroutine local_coordinates

!FIXME rock-related.. move out?


!> Called by CPU of zero rank. Sends the 'release' command. 
subroutine send_rock_release()
  integer :: buffer(4),ierror,np
  character(len=80) :: mpimsg

  ! SEND 99999 TO RELEASE WAITING PROCESSORS
  buffer = 99999
  do np=nprocs-1,1,-1
    call MPI_send(  &
      buffer,       & ! buf
      4,            & ! length
      MPI_INTEGER,  & ! datatype
      np,           & ! dest
      tag_rock,     & ! tag
      Comm_cart,    & ! communicator
      ierror)
    mpimsg = 'Failed to send release processor command'
    call checkmpi(ierror,mpimsg)
  enddo
end subroutine send_rock_release

!> Called by CPUs of nonzero rank. Waits for rank 0 to send the 
!> appropriate bit of rock data over, then fills in \c N%rock_state.
subroutine recv_rock_par(N)
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: stat(MPI_STATUS_SIZE),ierror,buffer(4)
  character(len=80) :: mpimsg

  buffer = 0
  ! MAGIC NUMBER 99999, see send_rock_release()
  do while (buffer(4).ne.99999)
    call MPI_Recv(  &	
      buffer,       & ! buf
      4,            & ! length
      MPI_INTEGER,  & ! datatype
      0,            & ! source
      tag_rock,     & ! tag
      Comm_Cart,    & ! communicator
      stat,         & ! status
      ierror)
    mpimsg = 'MPI_Recv() in recv_rock_par failed'
    call checkmpi(ierror,mpimsg)
    if (buffer(4).ne.99999) then
      N(buffer(1),buffer(2),buffer(3))%rock_state=dble(buffer(4))
    endif
  enddo
end subroutine recv_rock_par

end module lb3d_mpi_module
