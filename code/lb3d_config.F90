#include "lb3d.h"
#include "lb3d_revision.h"

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


!> Contain all configuration parameters as well as functions to read and distribute
!> \todo Clean up a lot
!> \todo Separate the different functionalities combined here. 
!> \todo This is a good candidate for automatic creation.
module lb3d_config_module

  use lb3d_log_module
  use lb3d_global_module!, only: rk
  use lb3d_mpi_parameters_module

  implicit none

  character(len=64)   :: lbeversion = 'LB3D version: v7.0' 
  character(len=1024) :: lbeflags   = 'LB3D flags: '//trim(LBEFLAGS)

  ! ======================================================================
  !                                FIXED INPUT
  ! ======================================================================

  integer :: nx = 16
  integer :: ny = 16
  integer :: nz = 16
  integer :: seed = 1111
  character(len=64) :: obs_file = 'empty.dat'
  character(len=64) :: obs_folder = '../rocks'
  character(len=3)  :: obs_rotation = 'xyz'
  integer :: boundary_cond = 0
  integer :: boundary_width = 2
  character(len=64) :: boundary = 'periodic'
  character(len=64) :: force = 'none'
  integer :: cdx = 0
  integer :: cdy = 0
  integer :: cdz = 0
  logical :: dbg_report_topology = .false.
  integer :: dbg_n_start_mpi_debug = 0

  namelist /FIXED_INPUT/ nx,ny,nz,seed,obs_file,obs_folder,obs_rotation,boundary_cond,boundary_width,boundary,force,cdx,cdy,cdz,dbg_report_topology,dbg_n_start_mpi_debug

  ! ======================================================================
  !                                VARIABLE INPUT
  ! ======================================================================

  logical :: sci_int = .false.
  logical :: sci_sur = .false.
  logical :: sci_od  = .false.
  logical :: sci_wd  = .false.
  logical :: sci_dir = .false.
  logical :: sci_vel = .false.
  logical :: sci_flo = .false.
  logical :: sci_arrows = .false.
  logical :: sci_velocities = .false.
  logical :: sci_rock = .false.
  logical :: sci_pressure = .false.
  logical :: sci_fluxz = .false.
  logical :: sci_massfluxz = .false.
  logical :: sci_profile = .false.
  logical :: sci_arrstats = .false.
  logical :: sci_perm = .false.
  logical :: sci_stress = .false.
  logical :: post = .false.

  character(len=64) :: gr_out_file = 'default'
  character(len=64) :: folder = 'default'
  character(len=64) :: cpfolder = '.'
  character(len=256) :: srccpfolder = '.'

  logical :: restore = .false.
  integer :: init_cond = 0
  integer :: sci_pressure_init = 1

  integer :: n_iteration = 100
  integer :: n_sci_start = 0
  integer :: n_sci_int   = 100
  integer :: n_sci_sur   = 100
  integer :: n_sci_od    = 100
  integer :: n_sci_wd    = 100
  integer :: n_sci_dir   = 100
  integer :: n_sci_vel   = 100
  integer :: n_sci_flo   = 100
  integer :: n_sci_arrows = 100
  integer :: n_sci_velocities = 100
  integer :: n_sci_rock  = 100
  integer :: n_sci_pressure = 100
  integer :: n_sci_fluxz = 100
  integer :: n_sci_massfluxz = 100
  integer :: n_sci_profile = 100
  !> \c n_sci_profile is the sampling interval for profile data. The
  !> data sampled since the last dump is averaged and dumped directly
  !> after sampling but only if \c mod(nt,n_sci_profile_dump)==0
  !> . This has two consequences: 1) As long as \c n_sci_profile_dump
  !> stays at its default value of \c 1, \c n_sci_profile controls the
  !> sampling and the dumping interval as well. 2) If averaged output
  !> is desired \c n_sci_profile_dump should be set to an integer
  !> multiple of \c n_sci_profile .
  integer :: n_sci_profile_dump = 1
  integer :: n_sci_arrstats = 100
  integer :: n_sci_perm = 100
  integer :: n_sci_stress = 100

  integer :: stress_minx = -1
  integer :: stress_maxx = -1

  !> maximum number of regions to calculate and dump fluxz data for
  integer,parameter,public :: fluxz_regions = 10

  !> Human-readable identifier for each fluxz region, default values
  !> are assigned later in case they are not overwritten in the input
  !> file.
  character(len=32),save,public :: fluxz_name(fluxz_regions) = ""

  !> \name lower boundaries for each fluxz region
  !> \{
  integer,save,public :: fluxz_xlo(fluxz_regions) = 1
  integer,save,public :: fluxz_ylo(fluxz_regions) = 1
  integer,save,public :: fluxz_zlo(fluxz_regions) = 1
  !> \}

  !> \name upper boundaries for each fluxz region
  !> \{
  !> value -1 for x-boundary indicates an unused slot, value 0 the
  !> maximum lattice position
  integer fxhidx                !< dummy index for in-line initialization
  integer,save,public :: fluxz_xhi(fluxz_regions)&
       & = (/0,(-1,fxhidx=2,fluxz_regions)/)
  integer,save,public :: fluxz_yhi(fluxz_regions) = 0
  integer,save,public :: fluxz_zhi(fluxz_regions) = 0
  !> \}

  !> maximum number of different intervals over which arrstats data is
  !> accumulated separately
  integer,parameter,public :: arrstats_intervals = 10
  integer :: asidx              !< dummy index for in-line initialization
  !> dump accumulated arrstats data each time step arrstats data is
  !> sampled and \c mod(nt,n_sci_arrstats_dump)==0. A value of \c 0
  !> disables the respective slot.
  integer,save,public :: n_sci_arrstats_dump(10)&
       & = (/1,(0,asidx=2,arrstats_intervals)/)
  !> Human-readable identifier for each arrstats region, default values
  !> are assigned later in case they are not overwritten in the input
  !> file.
  character(len=32),save,public :: arrstats_name(arrstats_intervals) = ""

  integer :: nr = 16
  real(kind=rk) :: fr = 1.0, fb = 0.5, fg = 0.1
  integer :: fd = 0
  real(kind=rk) :: fr1 = 0.2, fr2 = 0.3
  real(kind=rk) :: pr = 0.5, pb = 1.0, pg = 0.1
  integer :: pd = 0
  real(kind=rk) :: qr = 0.5, qb = 0.5, qg = 0.1
  integer :: qd = 0
  real(kind=rk) :: r1 = 10.0
  real(kind=rk) :: m_evp = 0.001
  real(kind=rk) :: mf_r = 0.001, mf_b = 0.001
  logical,save :: in_evp(3)=.false.
  logical,save :: out_evp(3)=.false.
  real(kind=rk) :: rock_colour = 0.0
  real(kind=rk) :: rock_colour_r = 0.0
  real(kind=rk) :: rock_colour_b = 0.0
  integer :: inv_fluid = 0
  integer :: inv_type = 0
  real(kind=rk), save :: beta = 1.0    !< Dipole temperature

  !added for shifting/cutting droplets
  integer :: drop_xshift = 0
  integer :: drop_yshift = 0
  integer :: drop_zshift = 0
  integer :: drop_xcut = 0
  integer :: drop_ycut = 0
  integer :: drop_zcut = 0

  namelist /VARIABLE_INPUT/ &
       &arrstats_name,sci_arrstats,sci_int,sci_sur,sci_od,sci_wd,sci_dir&
       &,sci_vel,sci_flo,sci_arrows,sci_velocities,sci_rock,sci_pressure&
       &,sci_fluxz,sci_massfluxz,sci_profile,sci_perm,sci_stress,post,gr_out_file&
       &,folder,cpfolder,srccpfolder,restore,init_cond,sci_pressure_init&
       &,n_iteration,n_sci_arrstats,n_sci_arrstats_dump,n_sci_start,n_sci_int&
       &,n_sci_sur,n_sci_od,n_sci_wd,n_sci_dir,n_sci_vel,n_sci_flo,n_sci_arrows&
       &,n_sci_velocities,n_sci_rock,n_sci_pressure,n_sci_fluxz,n_sci_massfluxz&
       &,n_sci_profile,n_sci_profile_dump,n_sci_perm,n_sci_stress,fluxz_name&
       &,fluxz_xlo,fluxz_ylo,fluxz_zlo,fluxz_xhi,fluxz_yhi,fluxz_zhi,fr,fb&
       &,fg,fd,fr1,fr2,pr,pb,pg,pd,qr,qb,qg,qd,rock_colour,rock_colour_r,rock_colour_b,inv_fluid,inv_type&
       &,beta,r1,m_evp,in_evp,out_evp,drop_xshift,drop_yshift,drop_zshift&
       &,drop_xcut,drop_ycut,drop_zcut,stress_minx,stress_maxx

  ! ======================================================================
  !                                LBE INPUT
  ! ======================================================================

  ! These variables are intrinsic to a given simulation run.

  ! Interaction types
  ! Single Component Multi Phase - pos. intra-comp. SC force
  logical        :: SCMP = .false. 
  ! Multi Component Multi Phase - pos. intra-comp.-, neg. inter-comp.-SC-force
  logical        :: MCMP = .false.

  ! Declare different collision model types:
  ! BGK - Single relaxation time model (=BGK)
  ! MRT - multiple relaxation time scheme (from Sebastian)
  ! FLUCTUATING_MRT - fluctuating mrt scheme from Duenweg et al (implemented by Philipp)
  integer, parameter :: BGK=1,MRT=2,FLUCTUATING_MRT=3
  integer        :: COLLISIONTYPE_ID=BGK
  character(len=32), save :: COLLISIONTYPE = 'BGK'
  ! For fluctuating MRT, we also need to know the value k_B*T
  ! where k_B is Boltzmann's constant and T is the temperature of the system. Default
  ! value is 1e-7
  real(kind=rk)         :: kbT=0.0000001d0
  logical        :: OXFORD = .false.
  logical        :: INDEWALL = .false.
  logical        :: ZEROFORCEOFFSET = .true.
  logical        :: ZFOSdiag = .false.
 
 ! Not in namelist
  integer, save  :: n_out_max     ! total number of outputs.
  integer, save  :: nc_max = 10000
  integer, save  :: nt_max        ! No of timesteps between outputs

  ! Molecular masses
  real(kind=rk), save   :: amass_r = 1.0, amass_b = 1.0 , amass_s = 1.0
  ! Relaxation times
  real(kind=rk), save   :: tau_r = 1.0, tau_b = 1.0, tau_s = 1.0, tau_d = 2.0
  real(kind=rk), save   :: taubulk_r = 0.84, taubulk_b = 0.84, taubulk_s = 0.84, taubulk_d = 1.68
  real(kind=rk), save   :: s03_r = 1.0, s05_r = 1.0, s11_r = 1.0, s14_r = 1.0, s17_r = 1.0
  real(kind=rk), save   :: s03_b = 1.0, s05_b = 1.0, s11_b = 1.0, s14_b = 1.0, s17_b = 1.0
  real(kind=rk), save   :: s03_s = 1.0, s05_s = 1.0, s11_s = 1.0, s14_s = 1.0, s17_s = 1.0
  integer, save  :: bcsel = 0
  real(kind=rk), save   :: acccoef = 2
  ! Inverse relaxation times
  real(kind=rk), save   :: omega_b, omega_r, omega_s, omega_d ! Not in namelist
  real(kind=rk), save   :: omegabulk_b, omegabulk_r, omegabulk_s, omegabulk_d  ! Not in namelist

  ! Interaction strengths
  real(kind=rk), save   :: g_br = 0.0 ,g_bs = 0.0 , g_ss = 0.0
  real(kind=rk), save   :: g_rr = 0.0, g_bb = 0.0
  real(kind=rk), save   :: g_wr = 0.0, g_wb = 0.0
  real(kind=rk), save   :: tau_wr = 1.0, tau_wb = 1.0

  real(kind=rk), save   :: g_accn = 0.  ! Acceleration due to gravity
  real(kind=rk), save   :: g_accn_x = 0.  ! Acceleration due to gravity
  real(kind=rk), save   :: g_accn_y = 0.  ! Acceleration due to gravity
  integer, save  :: g_accn_min = 0    ! Acceleration due to gravity
  integer, save  :: g_accn_max = 0    ! Acceleration due to gravity
  integer, save  :: g_accn_min_x = 0    ! Acceleration due to gravity in x
  integer, save  :: g_accn_max_x = 0    ! Acceleration due to gravity in x 
  integer, save  :: g_accn_min_y = 0    ! Acceleration due to gravity in y
  integer, save  :: g_accn_max_y = 0    ! Acceleration due to gravity in y

  integer, save  :: n_out            ! No of outputs written +1 - not in namelist
  real(kind=rk), save   :: perturbation = 0.00 ! Perturbation to initial state.
  integer, save  :: perm_iolet = 0
  integer, save  :: perm_wall = 0

  ! Edge fluid's speed in shear system and shear rate frequency.
  real(kind=rk), save   :: shear_u = 0.
  real(kind=rk), save   :: shear_omega = 0.
  real(kind=rk), save   :: shear_tmp = 0.

  integer, save  :: n_checkpoint = 300 ! Number of timesteps b/t checkpoints
  character(len=3)  :: checkpoint_format = 'xdr'
  logical, save  :: checkpoint_safe = .true.
  logical, save  :: checkpoint_shearsum = .true.
  integer, save  :: n_restore = 0 ! Timestep to restore from. - not in namelist
  ! Unique identifier for chkpoint files - not in namelist
  integer, save  :: chk_uid
  ! chk uid & timestep to restore from.
  character(len=32), save :: restore_string = '0000000000-t00000000'

  integer, save  :: num_chkp_files = 0

  integer, save  :: psifunc = 2    ! Which form of \psi to use in force terms
  integer, save  :: bdist = 0
  real(kind=rk), save   :: d_0 = 1        ! Maximum modulus of dipoles

  character(len=3)  :: dump_format = 'bin'
  logical :: xdrfsloppy = .false.
  logical :: write_AVS_fld = .false.
  logical :: dump_double = .true.

  character(len=80), save :: inp_file ! Name of the input file - not in namelist

  ! These variables are set through command-line options
  ! arg_foo_p is set to 1 if arg_foo is defined.
  ! Not in namelist!
  character(len=32), save :: arg_restore_string
  integer, save           :: arg_restore_string_p = 0
  character(len=80),save  :: arg_input_file
  integer, save           :: arg_input_file_p = 0
  character(len=80),save  :: arg_input_dfile
  integer, save           :: arg_input_dfile_p = 0

  ! interval after which sanity check is performed (0 means never)
  integer,save :: n_sanity_check = 1000

  namelist /LBE_INPUT/ &
    SCMP, MCMP, COLLISIONTYPE, kbT, OXFORD, INDEWALL, &
    ZEROFORCEOFFSET, ZFOSdiag, &
    amass_r, amass_b, amass_s, &
    tau_r, tau_b, tau_s, tau_d, &
    taubulk_r, taubulk_b, taubulk_s, taubulk_d, &
    s03_r, s05_r, s11_r, s14_r, s17_r, &
    s03_b, s05_b, s11_b, s14_b, s17_b, &
    s03_s, s05_s, s11_s, s14_s, s17_s, &
    bcsel, acccoef, &
    g_br, g_bs, g_ss, g_rr, g_bb, g_wr, g_wb, tau_wr, tau_wb, bdist, &
    g_accn, g_accn_x, g_accn_y, g_accn_min, g_accn_max, &
    g_accn_min_x, g_accn_max_x, g_accn_min_y, g_accn_max_y, &
    perturbation, perm_iolet, perm_wall,shear_omega,shear_u, &
    n_checkpoint, checkpoint_format, checkpoint_safe, checkpoint_shearsum, &
    restore_string, num_chkp_files, psifunc, d_0, &
    dump_format, xdrfsloppy, write_AVS_fld, dump_double, n_sanity_check

contains

  !> Reads all configuration information from the input-file
  subroutine lb3d_config_init (stage)

    integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_config_init stage',I0)") stage
    CALL log_msg(trim(msgstr),.false.) 
#endif    
    select case (stage)
       
    case (0) ! pre mpi init
       
       ! Read the fixed input.

       CALL lbe_get_fixed_input()
       call lbe_get_lbe_input()
       
       ! Read the variable input.
       call lbe_get_variable_input()
          
      
    case (2) ! post mpi init / pre memory alloc
       !FIXME detach Communicate config

     case default

#ifdef LB3D_DEBUG_INFO
        call log_msg('nothing to do.',.false.)
#endif    
    end select

  end subroutine lb3d_config_init

!> Reads the input file, and parses the \c /FIXED_INPUT/ namelist - see
!> the User's Guide for a description of the variables.
!>
!> Rank zero does the reading, and broadcasts the values to all other CPUs.
subroutine lbe_get_fixed_input()
  implicit none
  character(len=128)     :: msgstr
  integer                :: ierror

  ! Get fixed input:
  ! Read this data only on processor 0 and then broadcast it to the other processors.
  ! Processor reads in .input-file from where it expects to read the main file

  ! Use myrankw instead of myrankc here because the cartesian grid
  ! hasn't been set up yet
  if (myrankw == 0) then
    CALL log_msg_ws("-----------( Reading fixed input )-----------",.false.)
 !   CALL lbe_define_inp_file()
    open (unit=input_file_unit,file=inp_file,status='UNKNOWN')
    read (unit=input_file_unit,nml=FIXED_INPUT)
    close (unit=input_file_unit)

    if (arg_input_dfile_p > 0) then
      CALL log_msg("  Getting differential input...",.false.)
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = FIXED_INPUT, IOSTAT = ierror)
      if (ierror .ne. 0) then
        CALL log_msg("    WARNING: Differential namelist not found or errors encountered.",.false.)
      endif
      close(UNIT = input_dfile_unit)
      CALL log_ws(.false.)
    end if

    write(msgstr,"('nx = ',i0,', ny = ',i0,', nz = ', i0)") nx, ny, nz
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('seed = ',i0,', boundary_cond = ',i0)") seed, boundary_cond
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('boundary_width = ',i0)") boundary_width
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('boundary = <',A,'>')") trim(boundary)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('obs_file = <',A,'>')") trim(obs_file)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('obs_folder = <',A,'>')") trim(obs_folder)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('obs_rotation = <',A,'>')") trim(obs_rotation)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('force    = <',A,'>')") trim(force)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('cdx = ',i0,', cdy = ',i0,', cdz = ', i0)") cdx, cdy, cdz
    CALL log_msg(trim(msgstr),.false.)
    if ( dbg_report_topology ) then
      write(msgstr,"('dbg_report_topology = ',L1)") dbg_report_topology
      CALL log_msg(trim(msgstr),.false.)
    endif
#ifdef DEBUG_MPI
    write(msgstr,"('dbg_n_start_mpi_debug = ',I0)") dbg_n_start_mpi_debug
    CALL log_msg(trim(msgstr),.false.)
#endif
  endif

  CALL MPI_Bcast(nx,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(ny,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(nz,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(seed,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(obs_file,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(obs_folder,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(obs_rotation,3,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(boundary_cond,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(boundary_width,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(boundary,64,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(force,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(cdx,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(cdy,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(cdz,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(dbg_report_topology,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(dbg_n_start_mpi_debug,1,MPI_INTEGER,0,comm_cart,ierror)

  ! Set the suggested MPI Cartesian dimensions
  cdims(1) = cdx
  cdims(2) = cdy
  cdims(3) = cdz
end subroutine lbe_get_fixed_input

!> Reads the input file, and parses the \c /VARIABLE_INPUT/ namelist -
!> see the User's Guide for a description of the variables.
!>
!> Rank zero does the reading, and broadcasts the values to all other
!> CPUs.
subroutine lbe_get_variable_input()
  implicit none
  integer :: i,ierror,eof_err
  character(len=128)     :: msgstr

  if (myrankc == 0) then ! Let only processor 0 read from file.
    CALL log_msg_ws("----------( Reading variable input )---------",.false.)
    open (unit=input_file_unit,file=inp_file,status='UNKNOWN')
    read (unit=input_file_unit,nml=VARIABLE_INPUT,iostat=eof_err)
    if(eof_err .lt. 0) then
      call log_msg("FATAL ERROR: End of input file encountered, aborting..."&
           &,.false.)
      close (unit=input_file_unit)
    !FIXME This should be handled outside of this module
    !  call Abend
      call MPI_Abort(MPI_COMM_WORLD,-1,ierror)
    else if(eof_err .gt. 0) then
      call log_msg("FATAL ERROR: Error in input file encountered, aborting..."&
           &,.false.)
      call log_msg("NOTE: for consistency, the old sci_start is now "&
           &//"n_sci_start. Similarly, n_sci is now n_sci_int",.false.)
      close (unit=input_file_unit)
    !FIXME This should be handled outside of this module
    !  call Abend
      call MPI_Abort(MPI_COMM_WORLD,-1,ierror)
    end if
    close(UNIT = input_file_unit)

    if (arg_input_dfile_p > 0) then
      CALL log_msg("  Getting differential input...",.false.)
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = VARIABLE_INPUT, IOSTAT = ierror)
      if (ierror .ne. 0) then
        CALL log_msg("    WARNING: Differential namelist not found or errors encountered.",.false.)
      endif
      close(UNIT = input_dfile_unit)
      CALL log_ws(.false.)
    end if

    if (0 .ne. arg_restore_string_p) then
      restore = .true.
    endif

    write(msgstr,"('n_iteration    = ',I0)") n_iteration
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('n_sci_start    = ',I0)") n_sci_start
    CALL log_msg(trim(msgstr),.false.)

    CALL log_ws(.false.)
    write(msgstr,"('sci_od         = ',L1,', n_sci_od           = ',I0)") sci_od, n_sci_od
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('sci_wd         = ',L1,', n_sci_wd           = ',I0)") sci_wd, n_sci_wd
    CALL log_msg(trim(msgstr),.false.)
#endif
    write(msgstr,"('sci_int        = ',L1,', n_sci_int          = ',I0)") sci_int, n_sci_int
    CALL log_msg(trim(msgstr),.false.)
#ifndef NOSURFACTANT
    write(msgstr,"('sci_sur        = ',L1,', n_sci_sur          = ',I0)") sci_sur, n_sci_sur
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_dir        = ',L1,', n_sci_dir          = ',I0)") sci_dir, n_sci_dir
    CALL log_msg(trim(msgstr),.false.)
#endif
    write(msgstr,"('sci_vel        = ',L1,', n_sci_vel          = ',I0)") sci_vel, n_sci_vel
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_flo        = ',L1,', n_sci_flo          = ',I0)") sci_flo, n_sci_flo
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_arrows     = ',L1,', n_sci_arrows       = ',I0)") sci_arrows, n_sci_arrows
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_velocities = ',L1,', n_sci_velocities   = ',I0)") sci_velocities, n_sci_velocities
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_rock       = ',L1,', n_sci_rock         = ',I0)") sci_rock, n_sci_rock
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_pressure   = ',L1,', n_sci_pressure     = ',I0)") sci_pressure, n_sci_pressure
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('                    sci_pressure_init  = ',I0)") sci_pressure_init
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_profile    = ',L1,', n_sci_profile      = ',I0)") sci_profile, n_sci_profile
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('                    n_sci_profile_dump = ',I0)") n_sci_profile_dump
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_arrstats   = ',L1,', n_sci_arrstats     = ',I0)") sci_arrstats, n_sci_arrstats
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_fluxz      = ',L1,', n_sci_fluxz        = ',I0)") sci_fluxz, n_sci_fluxz
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_massfluxz  = ',L1,', n_sci_massfluxz    = ',I0)") sci_massfluxz, n_sci_massfluxz
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_perm       = ',L1,', n_sci_perm         = ',I0)") sci_perm, n_sci_perm
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('sci_stress     = ',L1,', n_sci_stress       = ',I0)") sci_stress, n_sci_stress
    CALL log_msg(trim(msgstr),.false.)
    CALL log_ws(.false.)

    if (any(fluxz_xhi /= -1)) then
      CALL log_msg("Regions for fluxz / massfluxz:",.false.)
      do i = 1, fluxz_regions
        if (fluxz_xhi(i) .ne. -1) then
          write(msgstr,"('  <',A,'>: x = [',I0,':',I0,'], y = [',I0,':',I0,'], z = [',I0,':',I0,']')") &
            trim(fluxz_name(i)), fluxz_xlo(i), fluxz_xhi(i), fluxz_ylo(i), fluxz_yhi(i), fluxz_zlo(i), fluxz_zhi(i)
          CALL log_msg(trim(msgstr),.false.)
        end if
      end do
    else
      CALL log_msg("No active fluxz regions defined",.false.)
    end if

    if (any(n_sci_arrstats_dump/=0)) then
      CALL log_msg("Intervals for arrstats:",.false.)
      do i = 1,arrstats_intervals
        if (n_sci_arrstats_dump(i)/=0) then
          write(msgstr,"('  <',A,'>: n_sci_arrstats_dump = ',I0)") &
               &trim(arrstats_name(i)),n_sci_arrstats_dump(i)
          CALL log_msg(trim(msgstr),.false.)
        end if
      end do
    else
      CALL log_msg("No active arrstats intervals defined",.false.)
    end if

    CALL log_ws(.false.)
    write(msgstr,"('post        = ',L1)") post
    CALL log_msg(trim(msgstr),.false.)

    write(msgstr,"('folder      = <',A,'>')") trim(folder)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('cpfolder    = <',A,'>')") trim(cpfolder)
    CALL log_msg(trim(msgstr),.false.)
    if (trim(srccpfolder) .ne. '') then
      write(msgstr,"('srccpfolder = <',A,'>')") trim(srccpfolder)
      CALL log_msg(trim(msgstr),.false.)
    endif
    write(msgstr,"('gr_out_file = <',A,'>')") trim(gr_out_file)
    CALL log_msg(trim(msgstr),.false.)

    write(msgstr,"('restore     = ',L1)") restore
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('init_cond   = ',I0)") init_cond
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('inv_fluid   = ',I0,', inv_type = ',I0)") inv_fluid, inv_type
    CALL log_msg(trim(msgstr),.false.)

    CALL log_ws(.false.)

    write(msgstr,"('fr          = ',F16.10)") fr
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('fb          = ',F16.10)") fb
    CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('fg          = ',F16.10)") fg
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('fd          = ',I0)") fd
    CALL log_msg(trim(msgstr),.false.)
#endif
    write(msgstr,"('fr1         = ',F16.10)") fr1
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('fr2         = ',F16.10)") fr2
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('pr          = ',F16.10)") pr
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('pb          = ',F16.10)") pb
    CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('pg          = ',F16.10)") pg
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('pd          = ',I0)") pd
    CALL log_msg(trim(msgstr),.false.)
#endif
    write(msgstr,"('qr          = ',F16.10)") qr
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('qb          = ',F16.10)") qb
    CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('qg          = ',F16.10)") qg
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('qd          = ',I0)") qd
    CALL log_msg(trim(msgstr),.false.)
#endif

    write(msgstr,"('r1          = ',F16.10)") r1
    CALL log_msg(trim(msgstr),.false.)

    write(msgstr,"('m_evp          = ',F16.10)") m_evp
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('in_evp  = (',L1,',',L1,',',L1,')')") in_evp(1), in_evp(2), in_evp(3)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('out_evp  = (',L1,',',L1,',',L1,')')") out_evp(1), out_evp(2), out_evp(3)
    CALL log_msg(trim(msgstr),.false.)

    !FIXME
    write(msgstr,"('rock_colour = ',F16.10)") rock_colour
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('WARNING: rock_colour is deprecated')")
    write(msgstr,"('Please use rock_colour_r and rock_colour_b')")
    if ( rock_colour > 0 ) then 
       rock_colour_r = rock_colour
    elseif ( rock_colour < 0 ) then
       rock_colour_b = rock_colour
    end if

    CALL log_msg(trim(msgstr),.false.)

    write(msgstr,"('rock_colour_r = ',F16.10)") rock_colour_r
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('rock_colour_r = ',F16.10)") rock_colour_b
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('beta        = ',F16.10)") beta
    CALL log_msg(trim(msgstr),.false.)

    if ( init_cond == 11 ) then
      CALL log_ws(.false.)
      write(msgstr,"('drop_xshift = ',I0,', drop_xcut = ',I0)") drop_xshift, drop_xcut
      CALL log_msg(trim(msgstr),.false.)
      write(msgstr,"('drop_yshift = ',I0,', drop_ycut = ',I0)") drop_yshift, drop_ycut
      CALL log_msg(trim(msgstr),.false.)
      write(msgstr,"('drop_zshift = ',I0,', drop_zcut = ',I0)") drop_zshift, drop_zcut
      CALL log_msg(trim(msgstr),.false.)
    endif

    if ( sci_stress ) then
      write(msgstr,"('stress_minx = ',I0,', stress_maxx = ',I0)") stress_minx,stress_maxx
      CALL log_msg(trim(msgstr),.false.)
    endif

  endif

  ! These variables match the order in lbe_params.F90, for clarity

  CALL MPI_Bcast(sci_int,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_sur,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_od,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_wd,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_dir,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_vel,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_flo,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_arrows,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_velocities,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_rock,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_pressure,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_fluxz,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_massfluxz,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_profile,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_arrstats,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_perm,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_stress,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(post,1,MPI_LOGICAL,0,comm_cart,ierror)

  CALL MPI_Bcast(gr_out_file,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(folder,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(cpfolder,64,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(srccpfolder,256,MPI_CHARACTER,0,comm_cart,ierror)

  CALL MPI_Bcast(restore,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(init_cond,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(sci_pressure_init,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(n_iteration,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_start,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_int,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_sur,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_od,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_wd,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_dir,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_vel,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_flo,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_arrows,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_velocities,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_rock,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_pressure,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_fluxz,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_massfluxz,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_profile,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_profile_dump,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_arrstats,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_perm,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_sci_stress,1,MPI_INTEGER,0,comm_cart,ierror)

  do i = 1, fluxz_regions
      CALL MPI_Bcast(fluxz_name(i),32,MPI_CHARACTER,0,comm_cart,ierror)
  end do

  CALL MPI_Bcast(fluxz_xlo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fluxz_ylo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fluxz_zlo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fluxz_xhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fluxz_yhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fluxz_zhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)

  call MPI_Bcast(n_sci_arrstats_dump,arrstats_intervals,MPI_INTEGER,0,comm_cart&
       &,ierror)
  do i = 1,arrstats_intervals
      call MPI_Bcast(arrstats_name(i),32,MPI_CHARACTER,0,comm_cart,ierror)
  end do

  CALL MPI_Bcast(nr,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(fb,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(fg,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(fd,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(fr1,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(fr2,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(pr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(pb,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(pg,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(pd,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(qr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(qb,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(qg,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(qd,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(r1,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(m_evp,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(in_evp,3,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(out_evp,3,MPI_LOGICAL,0,comm_cart,ierror)

  CALL MPI_Bcast(rock_colour,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(rock_colour_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(rock_colour_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(inv_fluid,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(inv_type,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(beta,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(stress_minx,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(stress_maxx,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(drop_xshift,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(drop_yshift,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(drop_zshift,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(drop_xcut,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(drop_ycut,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(drop_zcut,1,MPI_INTEGER,0,comm_cart,ierror)
end subroutine lbe_get_variable_input

!> This routine is called for all values of init_cond apart from
!> \c INIT_HUDONG; it reads the LBE-specific namelist from the input file and
!> broadcasts the data to all processors.
subroutine lbe_get_lbe_input()
  implicit none
  integer :: ierror,itmp
  character(len=128)     :: msgstr

  if (myrankc == 0) then
    CALL log_msg_ws("------------( Reading LBE input )------------",.false.)
    open(UNIT = input_file_unit, FILE = inp_file, STATUS = 'UNKNOWN')
    read(UNIT = input_file_unit, NML = lbe_input)
    close(unit=input_file_unit)

    if (arg_input_dfile_p > 0) then
      CALL log_msg("  Getting differential input...",.false.)
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = lbe_input, IOSTAT = ierror)
      if (ierror .ne. 0) then
        CALL log_msg("    WARNING: Differential namelist not found or errors encountered.",.false.)
      endif
      close(UNIT = input_dfile_unit)
      CALL log_ws(.false.)
    end if

    if (0 .ne. arg_restore_string_p) then
      restore_string = arg_restore_string
      write(msgstr,"('Setting restore_string from command-line to <',A,'>')") trim(restore_string)
      CALL log_msg(trim(msgstr),.false.)
    endif
  end if

  CALL MPI_Bcast(amass_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(tau_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(taubulk_r,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(amass_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(tau_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(taubulk_b,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(amass_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(tau_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(taubulk_s,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(tau_d,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(taubulk_d,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(g_rr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_br,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_bb,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_bs,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_ss,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(g_wr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_wb,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(tau_wr,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(tau_wb,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(d_0,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(perturbation,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(acccoef,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(g_accn,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_x,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_y,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(shear_u,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(shear_omega,1,MPI_REAL8,0,comm_cart,ierror)

  CALL MPI_Bcast(s03_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s05_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s11_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s14_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s17_r,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s03_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s05_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s11_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s14_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s17_b,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s03_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s05_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s11_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s14_s,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(s17_s,1,MPI_REAL8,0,comm_cart,ierror)

  omega_b = 1.0_rk/tau_b
  omega_r = 1.0_rk/tau_r
  omega_s = 1.0_rk/tau_s
  omega_d = 1.0_rk/tau_d
  omegabulk_b = 1.0_rk/taubulk_b
  omegabulk_r = 1.0_rk/taubulk_r
  omegabulk_s = 1.0_rk/taubulk_s
  omegabulk_d = 1.0_rk/taubulk_d

  if (myrankc == 0) then
    ! ARP - convert from string containing UID & timestep to two 
    ! separate integers
    itmp = INDEX(STRING=restore_string, SUBSTRING="-")

    if(itmp .eq. 0)then
      CALL log_msg("FATAL ERROR: restore_string does not contain '-'. Aborting...",.false.)
      !FIXME This should be handled outside of this module
    !  CALL Abend
      call MPI_Abort(MPI_COMM_WORLD,-1,ierror)
    end if

    ! First character should be 't' so skip it
    READ(restore_string(2:itmp-1), FMT='(I8)') n_restore

    if( ( init_cond /= 7 ) .and. ( .not. restore ) )then
      ! Generate new UID because we don't want the one
      ! from the input file since we're not restarting
      ! from a checkpoint
      CALL system_clock(chk_uid)
      write(msgstr,"('Generated UID <',I10.10,'>.')") chk_uid
      CALL log_msg(trim(msgstr),.false.)
    else
      ! Read value obtained from input file
      READ(restore_string(itmp+1:), FMT='(I10)') chk_uid
      write(msgstr,"('Retrieved UID <',I10.10,'>.')") chk_uid
      CALL log_msg(trim(msgstr),.false.)
    end if

    ! On vector machines we restrict ourselves to psifunc=2. See
    ! lbe_collision for details.
#ifdef VECTORIZE
    psifunc = 2
    CALL log_msg("NOTE: psifunc is set to 2 since we compiled with -DVECTORIZE.",.false.)
#endif

    if ( g_accn_max == 0 ) g_accn_max = nz
    if ( g_accn_max_x == 0 ) g_accn_max_x = nx
    if ( g_accn_max_y == 0 ) g_accn_max_y = ny
  end if

  ! Broadcast integer, logical and character data as well.

  CALL MPI_Bcast(n_sanity_check,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_checkpoint,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(num_chkp_files,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(n_restore,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(chk_uid,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(bcsel,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(bdist,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(psifunc,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(g_accn_min,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_max,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_min_x,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_max_x,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_min_y,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(g_accn_max_y,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(perm_wall,1,MPI_INTEGER,0,comm_cart,ierror)
  CALL MPI_Bcast(perm_iolet,1,MPI_INTEGER,0,comm_cart,ierror)

  CALL MPI_Bcast(COLLISIONTYPE,32,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(kbT,1,MPI_REAL8,0,comm_cart,ierror)
  CALL MPI_Bcast(INDEWALL,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(ZEROFORCEOFFSET,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(ZFOSdiag,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(write_AVS_fld,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(dump_double,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(SCMP,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(MCMP,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(OXFORD,1,MPI_LOGICAL,0,comm_cart,ierror)

  CALL MPI_Bcast(dump_format,3,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(xdrfsloppy,1,MPI_LOGICAL,0,comm_cart,ierror)

  CALL MPI_Bcast(checkpoint_format,3,MPI_CHARACTER,0,comm_cart,ierror)
  CALL MPI_Bcast(checkpoint_safe,1,MPI_LOGICAL,0,comm_cart,ierror)
  CALL MPI_Bcast(checkpoint_shearsum,1,MPI_LOGICAL,0,comm_cart,ierror)

  CALL MPI_Bcast(restore_string,32,MPI_CHARACTER,0,comm_cart,ierror)

  CALL log_ws(.false.)

  write(msgstr,"('COLLISIONTYPE   = <',A,'>')") trim(COLLISIONTYPE)
  CALL log_msg(trim(msgstr),.false.)
  if (COLLISIONTYPE=="BGK") then
    CALL log_msg("  -> Using BGK collision model.",.false.)
    COLLISIONTYPE_ID = BGK
  else if (COLLISIONTYPE=="MRT") then
    CALL log_msg("  -> Using MRT collision model.",.false.)
    COLLISIONTYPE_ID = MRT
  else if (COLLISIONTYPE=="FLUCTUATING_MRT") then
    CALL log_msg("  -> Using fluctuating MRT collision model.",.false.)
    COLLISIONTYPE_ID = FLUCTUATING_MRT
  endif
  CALL log_ws(.false.)

  ! parse kb*T parameter
  write(msgstr,"('kbT             = ',F16.10)") kbT
  CALL log_msg(trim(msgstr),.false.)
  
  write(msgstr,"('INDEWALL        = ',L1)") INDEWALL
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('ZEROFORCEOFFSET = ',L1)") ZEROFORCEOFFSET
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('write_AVS_fld   = ',L1)") write_AVS_fld
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('dump_double     = ',L1)") dump_double
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('SCMP            = ',L1)") SCMP
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('MCMP            = ',L1)") MCMP
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('OXFORD          = ',L1)") OXFORD
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('dump_format     = <',A,'>')") trim(dump_format)
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('xdrfsloppy      = ',L1)") xdrfsloppy
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  write(msgstr,"('n_sanity_check      = ',I0)") n_sanity_check
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('n_checkpoint        = ',I0)") n_checkpoint
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('checkpoint_format   = <',A,'>')") trim(checkpoint_format)
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('checkpoint_safe     = ',L1)") checkpoint_safe
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('checkpoint_shearsum = ',L1)") checkpoint_shearsum
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('num_chkp_files      = ',I0)") num_chkp_files
  CALL log_msg(trim(msgstr),.false.)
  if ( ( init_cond == 7 ) .or. ( restore ) ) then
    write(msgstr,"('restore_string      = <',A,'>')") trim(restore_string)
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('  -> n_restore      = ',I0)") n_restore
    CALL log_msg(trim(msgstr),.false.)
  endif
  CALL log_ws(.false.)

  write(msgstr,"('bcsel        = ',I0)") bcsel
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('bdist        = ',I0)") bdist
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('psifunc      = ',I0)") psifunc
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  write(msgstr,"('amass_r      = ',F16.10)") amass_r
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('tau_r        = ',F16.10)") tau_r
  CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
  write(msgstr,"('amass_b      = ',F16.10)") amass_b
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('tau_b        = ',F16.10)") tau_b
  CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
  write(msgstr,"('amass_s      = ',F16.10)") amass_s
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('tau_s        = ',F16.10)") tau_s
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('tau_d        = ',F16.10)") tau_d
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)
#endif
  write(msgstr,"('g_rr         = ',F16.10)") g_rr
  CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
  write(msgstr,"('g_br         = ',F16.10)") g_br
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_bb         = ',F16.10)") g_bb
  CALL log_msg(trim(msgstr),.false.)
#ifndef NOSURFACTANT
  write(msgstr,"('g_bs         = ',F16.10)") g_bs
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_ss         = ',F16.10)") g_ss
  CALL log_msg(trim(msgstr),.false.)
#endif
#endif

  if (INDEWALL) then
    write(msgstr,"('g_wr         = ',F16.10)") g_wr
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('g_wb         = ',F16.10)") g_wb
    CALL log_msg(trim(msgstr),.false.)
    CALL log_ws(.false.)
    write(msgstr,"('tau_wr       = ',F16.10)") tau_wr
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('tau_wb       = ',F16.10)") tau_wb
    CALL log_msg(trim(msgstr),.false.)
    CALL log_ws(.false.)
  endif

#ifndef NOSURFACTANT
  write(msgstr,"('d_0          = ',F16.10)") d_0
  CALL log_msg(trim(msgstr),.false.)
#endif
  write(msgstr,"('perturbation = ',F16.10)") perturbation
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('acccoef      = ',F16.10)") acccoef
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  write(msgstr,"('shear_u      = ',F16.10)") shear_u
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('shear_omega  = ',F16.10)") shear_omega
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  write(msgstr,"('g_accn       = ',F16.10)") g_accn
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_x     = ',F16.10)") g_accn_x
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_y     = ',F16.10)") g_accn_y
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_min   = ',I0)") g_accn_min
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_min_x = ',I0)") g_accn_min_x
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_min_y = ',I0)") g_accn_min_y
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_max   = ',I0)") g_accn_max
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_max_x = ',I0)") g_accn_max_x
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('g_accn_max_y = ',I0)") g_accn_max_y
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  write(msgstr,"('perm_wall    = ',I0)") perm_wall
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('perm_iolet   = ',I0)") perm_iolet
  CALL log_msg(trim(msgstr),.false.)
  CALL log_ws(.false.)

  if ( COLLISIONTYPE_ID .eq. MRT .or. COLLISIONTYPE_ID .eq. FLUCTUATING_MRT ) then
    write(msgstr,"('taubulk_r    = ',F16.10)") taubulk_r
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('taubulk_b    = ',F16.10)") taubulk_b
    CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('taubulk_s    = ',F16.10)") taubulk_s
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('taubulk_d    = ',F16.10)") taubulk_d
    CALL log_msg(trim(msgstr),.false.)
#endif
    write(msgstr,"('s03_r        = ',F16.10)") s03_r
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s05_r        = ',F16.10)") s05_r
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s11_r        = ',F16.10)") s11_r
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s14_r        = ',F16.10)") s14_r
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s17_r        = ',F16.10)") s17_r
    CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
    write(msgstr,"('s03_b        = ',F16.10)") s03_b
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s05_b        = ',F16.10)") s05_b
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s11_b        = ',F16.10)") s11_b
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s14_b        = ',F16.10)") s14_b
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s17_b        = ',F16.10)") s17_b
    CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('s03_s        = ',F16.10)") s03_s
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s05_s        = ',F16.10)") s05_s
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s11_s        = ',F16.10)") s11_s
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s14_s        = ',F16.10)") s14_s
    CALL log_msg(trim(msgstr),.false.)
    write(msgstr,"('s17_s        = ',F16.10)") s17_s
    CALL log_msg(trim(msgstr),.false.)
#endif

  endif

  CALL log_msg_ws("Successfully read LBE input.",.false.)

end subroutine lbe_get_lbe_input

end module lb3d_config_module

