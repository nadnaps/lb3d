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

!> Contains subroutines calculating data and writing files
!> \todo Calculations of data should be moved to the \c lb3d_analysis_module
module lb3d_io_dump_data_module

#ifdef USEHDF
  use lb3d_io_hdf5_module
#endif

#ifdef USEXDRF
  use lb3d_io_xdrf_module
#endif

  use lb3d_global_module

  use lb3d_config_module

  use lb3d_io_helper_module, only:nflags,check_dump_now

  use lb3d_force_shanchen_module, only: lb3d_calc_sc_forces
  use lb3d_analysis_module
  use lb3d_lattice_module
  use lb3d_mpi_module
  use lb3d_log_module, only: log_msg

#ifdef MD
  use map_module, only: Mii_map
#endif
  implicit none

  !> \name declarations for PROFILE output
  !> \{

  !> type representing one point of a profile containing data accumulated
  !> within the layer perpendicular to the profile direction
  type layer
     !> total fluid velocity
     real(kind=rk) :: vf(3)
     !> total mass flow
     real(kind=rk) ::  mf(3)
     !> total fluid density
     real(kind=rk) :: rhof
#ifndef SINGLEFLUID
     !> velocity of red fluid
     real(kind=rk) :: v_r(3)
     !> mass flow of red fluid
     real(kind=rk) :: m_r(3)
     !> density of red fluid
     real(kind=rk) :: rho_r
     !> velocity of blue fluid
     real(kind=rk) :: v_b(3)
     !> mass flow of blue fluid
     real(kind=rk) :: m_b(3)
     !> density of blue fluid
     real(kind=rk) :: rho_b
#endif
     !> number density of fluid sites
     real(kind=rk) :: nfs
     !> number density of solid rock sites
     real(kind=rk) :: nrs
#ifdef MD
     !> particle velocity
     real(kind=rk) :: vp(3)
     !> number densities of particle centers
     real(kind=rk) :: npc
     !> number density of particle sites
     real(kind=rk) :: nps
#endif
  end type layer

  !> type representing a complete profile along a specific direction
  !>
  !> \note maybe \c l should be pointer instead of allocatable, but it
  !> seems to work
  type profile
     character(len=5) :: name   !< distinct part of file names
     !> unit vectors in profile and averaging directions
     real(kind=rk) :: d(3,3)
     real(kind=rk) :: vd(3,3)    !< directions to project velocities at
     integer :: s(3)          !< starting point of profile
     integer :: size(3) !< dimensions in profile and averaging directions
     type(layer),allocatable :: l(:) !< data for each point of the profile
  end type profile

  !> profiles for PROFILE output
  type(profile),allocatable,save,private :: prf(:)

  integer,save,private :: layer_mpitype !< custom mpi data type for layer
  integer,save,private :: sum_layer_mpiop !< custom mpi reduction operation

  !> summation buffer, used only by root process
  type(profile),allocatable,save,private :: prfsum(:)
  !> number of samples taken (also included in output; used only by
  !> root process)
  integer,save,private :: n_samples
  !> \}

  !> For dump_pressure()
  integer, parameter :: droplet=1,nondroplet=0

#ifdef BUGGYIFORT11
  !> bugfix for ifort 11
  !>
  !> There is a bug in Intel Fortran Compiler 11.0 and 11.1 (but
  !> probably not in versions 9.1 and 10.1) which leads to erroneous
  !> optimization of
  !> fluid_velocity_and_density_and_site_occupation(). Doing
  !> meaningless operations with the variable below seems to prevent at
  !> least ifort 11.1 from showing this behavior. Essential for this
  !> bug is that
  !> \li the code is compiled with the options -O3 and -ipo and that
  !> \li somewhere else in the code there is a call to a not
  !>     user-provided subroutine with another subroutine as
  !>     argument. It is not necessary that this subroutine call will
  !>     be executed. Here mpi_op_create() in lb3d_io_module is the
  !>     call that causes the trouble.
  integer,save,public :: ifort_11_bug_dummy = 0
#endif

contains

!> write info on which compiler flags were used to stdout
!>
!> 2010-05-03: Added by Stefan, use the 'find_flags.sh' script to easily
!> check if this is still up to date!
subroutine lbe_detect_flags()
  nflags = 0 ! Declared in lb3d_io_helper so HDF5 metadata output can also use this value
  CALL log_msg_ws("--------( Reporting compiler flags )---------",.false.)

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
#ifdef MPI_ALLGV_FASTER_THAN_GV
  CALL log_msg("  MPI_ALLGV_FASTER_THAN_GV",.false.)
  nflags = nflags + 1
#endif
#ifdef NOCALLSYSTEM
  CALL log_msg("  NOCALLSYSTEM",.false.)
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
#ifdef USEXDRF
  CALL log_msg("  USEXDRF",.false.)
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


!> This subroutine just calls the specific routines to dump the data to disk
subroutine dump_data(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N

  if (check_dump_now(sci_od, n_sci_od) ) then
    CALL log_msg("  Dumping OD...",.false.)
    CALL dump_od(N)
  endif
#ifndef SINGLEFLUID
  if (check_dump_now(sci_int, n_sci_int) ) then
    CALL log_msg("  Dumping COLOUR...",.false.)
    CALL dump_colour(N)
  endif
  if (check_dump_now(sci_wd, n_sci_wd) ) then
    CALL log_msg("  Dumping WD...",.false.)
    CALL dump_wd(N)
  endif
#endif
#ifndef NOSURFACTANT
  if (check_dump_now(sci_sur, n_sci_sur) ) then
    CALL log_msg("  Dumping SUR...",.false.)
    CALL dump_sur(N)
  endif
  if (check_dump_now(sci_dir, n_sci_dir) ) then
    CALL log_msg("  Dumping DIR...",.false.)
    CALL dump_dir(N)
  endif
#endif
  if (check_dump_now(sci_vel, n_sci_vel) ) then
    CALL log_msg("  Dumping VEL...",.false.)
    CALL dump_vel(N)
  endif
  if (check_dump_now(sci_flo, n_sci_flo) ) then
    CALL log_msg("  Dumping FLO...",.false.)
    CALL dump_flo(N)
  endif
  if (check_dump_now(sci_arrows, n_sci_arrows) .or. &
      check_dump_now(sci_velocities, n_sci_velocities) ) then
    CALL log_msg("  Dumping ARROWS / VELOCITIES...",.false.)
    CALL dump_arrows(N)
  endif
  if (check_dump_now(sci_rock, n_sci_rock) ) then
    CALL log_msg("  Dumping ROCK...",.false.)
    CALL dump_rock(N)
  endif
  ! Note: dump_pressure is not called if the flag NOSURFACTANT is defined
  ! This is automatically excluded by the init conditions (see lbe_init).
  if (check_dump_now(sci_pressure, n_sci_pressure) ) then
    CALL log_msg("  Dumping PRESSURE...",.false.)
    if(sci_pressure_init == 1) then
      CALL dump_pressure(N,droplet)
      CALL dump_pressure(N,nondroplet)
    else
      ! if((sci_pressure_init == 8) .or.                     &
      !   (sci_pressure_init == 3) .or.                     &
      !   (sci_pressure_init == -2)) then
      CALL dump_pressure(N,nondroplet)
      ! CALL dump_popul_zx(N)
    endif
  endif
end subroutine dump_data

  !> dumps output of  fluxz()  (see there) for every region
  !> specified by  fluxz_[xyz](lo|hi)  to a file named according
  !> to  fluxz_name
subroutine dump_fluxz(N, mass_flux)
  type(lbe_site),intent(in) :: N(0:,0:,0:)
  logical, intent(in) :: mass_flux
  integer :: i,dz,nf,hipos(3)
  real(kind=rk) :: sf
#ifdef MD
  integer :: np
  real(kind=rk) :: sp
#endif
  integer,parameter :: fz_file_unit=12
  character(len=256) fz_file_name

  do i=1,fluxz_regions
    if (fluxz_xhi(i)==-1) cycle ! value -1 indicates unused slot

    ! value 0 means maximum position in respective direction
    hipos = (/fluxz_xhi(i),fluxz_yhi(i),fluxz_zhi(i)/)
    where (hipos==0) hipos = (/tnx,tny,tnz/)

    CALL fluxz(N, mass_flux, fluxz_xlo(i), fluxz_ylo(i), fluxz_zlo(i), &
      hipos(1), hipos(2), hipos(3), dz, nf, sf&
#ifdef MD
      &, np, sp&
#endif
      &)

    if (myrankc==0) then
      if (mass_flux) then
        CALL lbe_make_filename_output(fz_file_name,'massfluxz_'//trim(fluxz_name(i)),'.asc',nt)
      else
        CALL lbe_make_filename_output(fz_file_name,'fluxz_'//trim(fluxz_name(i)),'.asc',nt)
      end if
      open (unit=fz_file_unit,file=fz_file_name,status='REPLACE',action='WRITE',recl=80)
#ifdef MD
      write (unit=fz_file_unit,fmt='(SS,I9,X,SS,I4,X,I7,X,SP,ES15.8,X,SS,I7,X,SP,ES15.8)') nt, dz, nf, sf, np, sp
#else
      write (unit=fz_file_unit,fmt='(SS,I9,X,SS,I4,X,I7,X,SP,ES15.8)') nt, dz, nf, sf
#endif
      close (fz_file_unit)
    end if
  end do
end subroutine dump_fluxz

!> Returns in  dz  the dimension in z-direction, in  n[fp]  the number of
!> lattice nodes occupied by fluid/ladd particles, and in s[fp] the sum of
!> the z-velocities at all fluid/ladd particle sites for the cuboidal chunk
!> parallel to the coordinate axis spanned between the points  (xl,yl,zl)
!> and  (xh,yh,zh) .
!>
!> See documentation on how to use these values. If compiled
!> without  MD  set, only the fluid related data is calculated and returned.
!>
!> \warning possible particle rotations are not taken into account yet
subroutine fluxz(N,mass_flux, xl,yl,zl,xh,yh,zh,dz,nf,sf&
#ifdef MD
  &,np,sp&
#endif
  &)
  type(lbe_site),intent(in) :: N(0:,0:,0:)
  logical,intent(in) :: mass_flux
  integer,intent(in) :: xl,yl,zl,xh,yh,zh
  integer,intent(out) :: dz,nf
  real(kind=rk),intent(out) :: sf
#ifdef MD
  integer,intent(out) :: np
  real(kind=rk),intent(out) :: sp

  integer,parameter :: fzc = 2 ! fluid and ladd particles
#else
  integer,parameter :: fzc = 1 ! fluid only
#endif
  integer :: m(fzc),nsum(fzc)
  real(kind=rk) :: s(fzc),ssum(fzc)
  integer :: x,y,z,ierror
  integer :: minx(3),maxx(3)
  real(kind=rk), dimension(nx,ny,nz,3) :: f_r
#ifndef SINGLEFLUID
  real(kind=rk), dimension(nx,ny,nz,3) :: f_b
#ifndef NOSURFACTANT
  real(kind=rk), dimension(nx,ny,nz,3) :: f_s
#endif
#endif

  m = 0
  s = 0.0_rk

  minx = max(1,(/xl,yl,zl/)+1-start)
  maxx = min((/nx,ny,nz/),(/xh,yh,zh/)+1-start)

  do x=minx(1),maxx(1)
    do y=minx(2),maxx(2)
      do z=minx(3),maxx(3)
#ifdef MD
        no_solid_rock: if (N(x,y,z)%rock_state>=0.0_rk) then
#endif
          fluid_or_particle: if (N(x,y,z)%rock_state==0.0_rk) then
#ifdef SINGLEFLUID
            CALL lb3d_calc_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
            CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r)
#else
            CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
            m(1) = m(1) + 1

            if (mass_flux) then
              s(1) = s(1) + (&
                &amass_r*(sum(g(:)*N(x,y,z)%n_r(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_r(:))*&
                &((f_r(x,y,z,3)+g_accn)/2.0d0))&
#ifndef SINGLEFLUID
                &+amass_b*(sum(g(:)*N(x,y,z)%n_b(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_b(:))*&
                &((f_b(x,y,z,3)+g_accn)/2.0d0))&
#ifndef NOSURFACTANT
                &+amass_s*(sum(g(:)*N(x,y,z)%n_s(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_s(:))*&
                &((f_s(x,y,z,3)+g_accn)/2.0d0))&
#endif
#endif
                &)
            else
              s(1) = s(1) + ( ( &
                &amass_r*(sum(g(:)*N(x,y,z)%n_r(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_r(:))*&
                &((f_r(x,y,z,3)+g_accn)/2.0d0))&
#ifndef SINGLEFLUID
                &+amass_b*(sum(g(:)*N(x,y,z)%n_b(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_b(:))*&
                &((f_b(x,y,z,3)+g_accn)/2.0d0))&
#ifndef NOSURFACTANT
                &+amass_s*(sum(g(:)*N(x,y,z)%n_s(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_s(:))*&
                &((f_s(x,y,z,3)+g_accn)/2.0d0))&
#endif
#endif
                & )/( &
                &+sum(g(:)*N(x,y,z)%n_r(:))*amass_r&
#ifndef SINGLEFLUID
                &+sum(g(:)*N(x,y,z)%n_b(:))*amass_b&
#ifndef NOSURFACTANT
                &+sum(g(:)*N(x,y,z)%n_s(:))*amass_s&
#endif
#endif
                &) )
            end if
#ifdef MD
          else fluid_or_particle
            m(2) = m(2)+1
            s(2) = s(2)+P(Mii_map(uid2i,nint(N(x,y,z)%rock_state)))%v(3)
#endif
          end if fluid_or_particle
#ifdef MD
        end if no_solid_rock
#endif
      end do
    end do
  end do

  CALL MPI_Reduce(m,nsum,fzc,MPI_INTEGER,MPI_SUM,0,comm_cart,ierror)
  CALL MPI_Reduce(s,ssum,fzc,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)
  dz = 1+zh-zl
  nf = nsum(1)
  sf = ssum(1)
#ifdef MD
  np = nsum(2)
  sp = ssum(2)
#endif
end subroutine fluxz

!> \{
!> \name Debugging output routines
!>
!> These dump some part of the system state to a file, for debugging
!>purposes. Not intended for use in production code, but very useful to
!>see what's going on inside.

!> This routine is unique, and intended for debugging. It writes a file in
!> "all" format, ie containing one line per lattice site. Each line
!> begins with three integers giving the coordinates of that site in the
!> local lattice, and then six floating-point numbers for the red, blue
!> and surfactant densities and dipoles, followed by an integer containing
!> the rock state at that site.
! subroutine dump_sites_all(N)
! 	implicit none
! 	type(lbe_site), dimension(0:,0:,0:) :: N
! 	character(LEN=256) :: filename
! 	real(kind=rk) :: nred,nblue,nsurf
! 	integer :: x,y,z,s,nxi,nyi,nzi
! 	real(kind=rk) :: red,surf,blue
! 
!         nxi = size(N,1)-2
!         nyi = size(N,2)-2
!         nzi = size(N,3)-2
! 
! 	call lbe_make_filename(filename,'sites','.all',nt)
! 	
! 	open(10,file=filename)
! 
! 	do z=1,nzi
! 	 do y=1,nyi
! 	  do x=1,nxi
! 		red=sum(N(x,y,z)%n_r(:)*g)
! #ifndef SINGLEFLUID
! 		blue=sum(N(x,y,z)%n_b(:)*g)
! #endif
! #ifndef NOSURFACTANT
! 		surf=sum(N(x,y,z)%n_s(:)*g)
! 		write(10,'(3i3,6f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,blue,surf,					&
! 			N(x,y,z)%d(:),					&
! 			int(N(x,y,z)%rock_state)
! #else
! #ifndef SINGLEFLUID
! 		write(10,'(3i3,2f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,blue,					&
! 			int(N(x,y,z)%rock_state)
! #else
! 		write(10,'(3i3,1f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,					&
! 			int(N(x,y,z)%rock_state)
! #endif
! #endif
! 	  end do
! 	 end do
! 	end do
! 
! 	close(10)
! 
! end subroutine dump_sites_all

!> This routine dumps a subdomain-plus-haloes worth of reals.
!> It's intended for debugging purposes - useful to examine the state of a
!> single CPU's subdomain.
!> 
!> It also writes an AVS field file for easy visualization.
! subroutine dump_debug_real_all(prefix,dump)
! 	implicit none
! 	real(kind=rk), dimension(0:,0:,0:) :: dump
! 	character(LEN=256) :: filename
! 	character(LEN=*) :: prefix
! 	integer :: x,y,z,nxi,nyi,nzi
!         nxi = size(dump,1)-2
!         nyi = size(dump,2)-2
!         nzi = size(dump,3)-2
! 
! 	call lbe_make_filename(filename,prefix,'.all',nt)
! 	
! 	open(10,file=filename)
! 
! 	do z=0,nzi+1
! 	 do y=0,nyi+1
! 	  do x=0,nxi+1
! 		write(10,'(3i3,f24.12)')	x,y,z, dump(x,y,z)
! 	  end do
! 	 end do
! 	end do
! 
! 	close(10)
! 	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
! 
! 	! Write a field file, too.
!         call dump_avs_fld(prefix,nxi,nyi,nzi,int(1))
! 
! end subroutine dump_debug_real_all
!> \}

!> \{
!> \name Output routines - Science-specific routines
!>
!>These routines are responsible for writing the science outputs when
!>required, such as the colour at each site.
!>
!>If a given science output, \c foo, is required, then \p lbe.f90 makes
!>a call to the subroutine \c dump_foo.
!>
!>\c dump_foo then compiles all the data required into one or more arrays,
!>and passes these to \c dump_scalar, \c dump_nscalar, \c dump_vector, etc,
!>depending on how many arrays are to be dumped and of what type.
!>
!>These routines in turn call \c dump_scalar_all, \c dump_scalar_bin,
!>etc depending on whether binary or ASCII output has been specified.
!>
!>Each of these routines is responsible for a given sort of science
!>output. It calls the appropriate array routine, which in turn calls the
!>appropriate file-format routine.

!> Wrapper function for calc_perm(N)
!>
!> Added by Stefan 2010-04-28
subroutine dump_perm(N)
  type(lbe_site),intent(in) :: N(0:,0:,0:)
#ifdef SINGLEFLUID
  CALL calc_perm(N)
#else
  CALL log_msg("    WARNING: Permeability calculations are implemented for SINGLEFLUID only. Bypassing calculation.",.false.)
#endif
end subroutine dump_perm

!> Calculates permeabilities 'on the fly'
!> Writes profile data and time evolution of the permeability
!>
!> Added by Stefan 2010-04-28
#ifdef SINGLEFLUID
subroutine calc_perm(N)
  type(lbe_site),intent(in) :: N(0:,0:,0:)
  character(len=128) :: msgstr
  integer, dimension(3) :: offset, dims
  ! Small z-arrays
  integer, dimension(nz) :: rockz, porez
  real(kind=rk),  dimension(nz) :: odz, velz, mfz
  ! Small z-buffers
  integer, dimension(nz) :: rockz_buf, porez_buf
  real(kind=rk),  dimension(nz) :: odz_buf, velz_buf, mfz_buf
  ! Large z-arrays
  integer, dimension(tnz) :: trockz, tporez
  real(kind=rk),  dimension(tnz) :: todz, tvelz, tmfz
  ! Force correction term for velocities
  real(kind=rk),  dimension(nx,ny,nz,3) :: f_r
  ! Sums
  integer :: totalrock, totalpore, totalsites, zsites, Ls
  real(kind=rk)  :: totalod, localod, totalvel, localvel, totalmf, localmf

  ! Physical quantities - there's a lot of them for clarity of the code
  real(kind=rk) :: perm_au, perm_mu, delta_x, delta_t, tau, avg_vel_S, avg_vel_PS, avg_rho_PS, avg_rho_bottom_P, avg_rho_top_P
  !real(kind=rk)  :: nu, eta, perm, avg_grad_p_S, cs2, avg_p_top_P, avg_p_bottom_P, avg_rho_bottom_PS, avg_rho_top_PS,kappa, avg_rho_top, avg_rho_bottom, avg_mf_top, avg_mf_bottom

  ! Wall offsets
  integer :: offsetil, offsetir, offsetjl, offsetjr
  ! Dummies
  integer :: i, j, k, ci, cj, ck, p, tx,ty,tz
  ! MPI ranks
  integer :: targetrank, sourcerank
  ! MPI coordinates
  integer, dimension(3) :: targetcc, sourcecc
  integer :: ierr
  integer :: mpitag = 1

  ! Z-coordinates
  integer :: targetz, z_bottom, z_top

  real(kind=rk) :: dynvis, lattice, deltap_p

  integer status(MPI_STATUS_SIZE)

  real(kind=rk) :: porosity

  integer,parameter :: fileunit=12
  character(len=256) filename

  logical :: fexist

  dims = (/nx,ny,nz/)

  trockz = 0
  tporez = 0
  todz = 0.0
  tvelz = 0.0
  tmfz = 0.0

  offset(1) = ccoords(1)*dims(1)
  offset(2) = ccoords(2)*dims(2)
  offset(3) = ccoords(3)*dims(3)

  offsetil = 0
  offsetir = 0
  offsetjl = 0
  offsetjr = 0

  if ( ccoords(1) .eq. 0 ) then
    offsetil = perm_wall
  endif
  if ( ccoords(1) .eq. cdims(1) - 1 ) then
    offsetir = perm_wall
  endif
  if ( ccoords(2) .eq. 0 ) then
    offsetjl = perm_wall
  endif
  if ( ccoords(2) .eq. cdims(2) - 1 ) then
    offsetjr = perm_wall
  endif

  ! ------------------------------------------------------------------------
  !                    CALCULATE SUMS
  ! ------------------------------------------------------------------------

  do k = 1, nz
    rockz(k) = 0
    porez(k) = 0
    odz(k) = 0.0
    velz(k) = 0.0
    mfz(k) = 0.0
    ! Don't sum over the walls
    do j = offsetjl + 1, ny - offsetjr
      do i = offsetil + 1, nx - offsetir
        if (N(i,j,k)%rock_state==0.0_rk) then
          porez(k) = porez(k) + 1

          ! Calculate local density
          localod = amass_r*sum(N(i,j,k)%n_r*g)
          odz(k) = odz(k) + localod

          ! Calculate local velocity
          CALL lb3d_calc_sc_forces(N,i,j,k,f_r)
          localvel = sum(N(i,j,k)%n_r(:)*g*cz)*amass_r-((f_r(i,j,k,3)/2.0d0)*sum(N(i,j,k)%n_r(:)*g)*amass_r)
          localvel = localvel / max(10.D-9,dble( sum(N(i,j,k)%n_r(:)*g)*amass_r))

          ! force extra-term
          tx=i+ccoords(1)*nx
          ty=j+ccoords(2)*ny
          tz=k+ccoords(3)*nz
          if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
              (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
              (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
            localvel = localvel-g_accn/2.0d0
          endif

          velz(k) = velz(k) + localvel

          localmf = localod * localvel
          mfz(k) = mfz(k) + localmf

        else ! rock site
          rockz(k) = rockz(k) + 1
        endif
      enddo
    enddo
  enddo

  ! ------------------------------------------------------------------------
  !                    SEND SUMS TO RANK 0
  ! ------------------------------------------------------------------------  

  if (myrankc .eq. 0) then
    do k = 1, nz
      trockz(k) = rockz(k)
      tporez(k) = porez(k)
      todz(k)   = odz(k)
      tvelz(k)  = velz(k)
      tmfz(k)   = mfz(k)
    enddo
    do p = 1, nprocs - 1
      CALL MPI_Cart_coords( Comm_cart, p, nd, sourcecc, ierr)
      CALL MPI_Recv(rockz_buf,nz,MPI_INTEGER,p,mpitag,Comm_cart,status,ierr)
      CALL MPI_Recv(porez_buf,nz,MPI_INTEGER,p,mpitag,Comm_cart,status,ierr)
      CALL MPI_Recv(odz_buf,  nz,MPI_REAL8  ,p,mpitag,Comm_cart,status,ierr)
      CALL MPI_Recv(velz_buf, nz,MPI_REAL8  ,p,mpitag,Comm_cart,status,ierr)
      CALL MPI_Recv(mfz_buf,  nz,MPI_REAL8  ,p,mpitag,Comm_cart,status,ierr)
      do k = 1, nz
        targetz = k+nz*sourcecc(3)
        trockz(targetz) = trockz(targetz) + rockz_buf(k)
        tporez(targetz) = tporez(targetz) + porez_buf(k)
        todz(targetz)   = todz(targetz)   + odz_buf(k)
        tvelz(targetz)  = tvelz(targetz)  + velz_buf(k)
        tmfz(targetz)   = tmfz(targetz)   + mfz_buf(k)
      enddo
    enddo
  else
    CALL MPI_Send(rockz,nz,MPI_INTEGER,0,mpitag,Comm_cart,ierr)
    CALL MPI_Send(porez,nz,MPI_INTEGER,0,mpitag,Comm_cart,ierr)
    CALL MPI_Send(odz,  nz,MPI_REAL8  ,0,mpitag,Comm_cart,ierr)
    CALL MPI_Send(velz, nz,MPI_REAL8  ,0,mpitag,Comm_cart,ierr)
    CALL MPI_Send(mfz,  nz,MPI_REAL8  ,0,mpitag,Comm_cart,ierr)
  endif

  ! ------------------------------------------------------------------------
  !                    DO CALCULATIONS
  ! ------------------------------------------------------------------------

  if (myrankc .eq. 0) then

    CALL lbe_make_filename_output(filename,'perm_profile','.asc',nt)

    open (unit=fileunit,file=filename,status='REPLACE',action='WRITE',recl=80)

    totalrock = 0
    totalpore = 0
    totalod   = 0.0
    totalvel  = 0.0
    totalmf   = 0.0

    ! Sum only over the sample, which runs from perm_iolet + 1 to tnz - perm_iolet in Fortran array notation
    z_bottom = perm_iolet + 1
    z_top = tnz - perm_iolet
    Ls = tnz - 2*perm_iolet

    write (unit=fileunit,fmt="('#      z     rock     pore               od              vel               mf')")

    do k = 1, tnz
      if (k .ge. z_bottom .and. k .le. z_top) then
        totalrock = totalrock + trockz(k)
        totalpore = totalpore + tporez(k)
        totalod = totalod + todz(k)
        totalvel = totalvel + tvelz(k)
        totalmf = totalmf + tmfz(k)
      endif
      write (unit=fileunit,fmt='(I8,X,I8,X,I8,X,F16.10,X,F16.10,X,F16.10)') k, trockz(k), tporez(k), todz(k), tvelz(k), tmfz(k)
    enddo

    zsites = (tny - 2*perm_wall) * (tnx - 2*perm_wall)
    totalsites = zsites * Ls

    ! if ( totalsites .ne. totalpore + totalrock ) then
    !   write(msgstr,"('Totalsites = ',I0,' totalpore = ',I0,' totalrock = ',I0)") totalsites, totalpore, totalrock
    !   CALL log_msg(trim(msgstr),.false.)
    !   CALL Abend
    ! endif

    delta_x = 1.0
    delta_t = 1.0
    ! delta_t = 1.0/sqrt(3.0)

    ! tau_r is actually tau-hat_
    tau = tau_r * delta_t

    ! cs2 = (delta_x * delta_x) / (3.0 * delta_t * delta_t) ! Eqn (7)

    ! nu = cs2 * delta_t * ( tau / delta_t - 0.5 ) ! Eqn (15)

    avg_vel_S = totalvel / (1.0 * totalsites)
    avg_vel_PS = totalvel / (1.0 * totalpore)
    avg_rho_PS = totalod / (1.0 * totalpore)

    ! avg_rho_bottom = todz(z_bottom)/ (1.0 * zsites)
    ! avg_rho_top = todz(z_top)/ ( 1.0 * zsites )

    avg_rho_bottom_P = todz(z_bottom)/ ( 1.0 * tporez(z_bottom) )
    avg_rho_top_P = todz(z_top)/ ( 1.0 * tporez(z_top) )

    ! avg_p_bottom_P = cs2 * avg_rho_bottom_P
    ! avg_p_top_P = cs2 * avg_rho_top_P

    ! avg_mf_bottom = tmfz(z_bottom) / (1.0 * tporez(z_bottom) )
    ! avg_mf_top = tmfz(z_top) / (1.0 * tporez(z_top) )

    ! avg_grad_p_S = (avg_p_top_P - avg_p_bottom_P) / ( (1.0*Ls - 1.0 ) * delta_x ) ! Option c) in Ariel's paper

    ! eta = avg_rho_PS * nu

    ! perm = -eta * avg_vel_S / avg_grad_p_S


    ! kappa = ( -nu*(1.0*Ls - 1.0) * delta_x * avg_vel_S * (avg_rho_bottom_P + avg_rho_top_P) ) / &
    ! ( 2.0 * ( avg_p_top_P - avg_p_bottom_P ) )

    perm_au = (avg_vel_PS * avg_rho_PS * Ls) / ( 2.0 * (avg_rho_bottom_P - avg_rho_top_P ) ) ! As calculated in zvelocity.F90
    perm_mu = ((avg_vel_S * avg_rho_PS * Ls) / (avg_rho_bottom_P - avg_rho_top_P ))* ( delta_x * delta_x * (2.0*tau - 1.0) / 2.0 )

    ! write(msgstr,"('    perm_au: ',F16.10,', perm_mu: ',F16.10)") perm_au, perm_mu
    ! CALL log_msg(trim(msgstr),.false.)
    close (fileunit)

    ! Write the permeability with timestamp separately
    CALL lbe_make_filename_output(filename,'perm_evol','.asc',0)
    inquire(file=filename,exist=fexist)
    if (fexist) then
      open (unit=fileunit,file=filename,status='OLD',position='APPEND',recl=80)
    else
      open (unit=fileunit,file=filename,status='NEW',position='APPEND',recl=80)
    endif
    write (unit=fileunit,fmt='(I0,X,F16.10)') nt, perm_mu
    close (fileunit)
  endif
end subroutine calc_perm
! endif SINGLEFLUID
#endif

#ifndef SINGLEFLUID
!>Calculates the colour field at each site, and then calls \c dump_scalar.
subroutine dump_colour(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for colour output",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (N(x,y,z)%rock_state == 0.) then
          scalar(x,y,z)=sum(N(x,y,z)%n_r(:)*g)-sum(N(x,y,z)%n_b(:)*g)
        else
          scalar(x,y,z)=0.0
        endif
      end do
    end do
  end do
  CALL dump_scalar(scalar,'colour')
  deallocate(scalar)
end subroutine dump_colour
! endif nSINGLEFLUID
#endif

#ifndef NOSURFACTANT
!> Calculates the surfactant density at each site, then calls \c dump_scalar.
subroutine dump_sur(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for surfactant output",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        scalar(x,y,z)=amass_s*sum(N(x,y,z)%n_s(:)*g)
      end do
    end do
  end do
  CALL dump_scalar(scalar,'sur')
  deallocate(scalar)
end subroutine dump_sur
! endif NOSURFACTANT
#endif

#ifndef NOSURFACTANT
!> 
!> Calculates the dipole moment at each site, then calls \c dump_vector.
subroutine dump_dir(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:,:), allocatable :: d
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(d(3,1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for dir output",.true.)
    return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        d(:,x,y,z)=N(x,y,z)%db(:)
      end do
    end do
  end do
  CALL dump_vector(d,'dir')
  deallocate(d)
end subroutine dump_dir

#endif

!> Calculates the averaged Z velocity at each site, then calls
!> \c dump_scalar.
subroutine dump_vel(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi,tz,tx,ty
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_s
#endif

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for vel output"&
          &,.true.)
     return	! Inability to write output is not quite fatal.
  end if
  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for vel output"&
          &,.true.)
     return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        ! Note this is u not u_tilde (no funny averaging with
        ! \tau If all your \tau's are the same then u = u_tilde,
        ! if not and you want u_tilde, strip the code out of
        ! lbe_collision.F90

        ! Get Shan Chen forces
        ! only calculate forces on fluid at non-rock sites
        if (N(x,y,z)%rock_state == 0.d0) then
#ifdef SINGLEFLUID
          CALL lb3d_calc_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if

        if (N(x,y,z)%rock_state == 0.d0) then

           scalar(x,y,z) = sum(N(x,y,z)%n_r(:)*g*cz)*amass_r*omega_r 
           if (SCMP) then
             scalar(x,y,z) = sum(N(x,y,z)%n_r(:)*g*cz)*amass_r*omega_r - tau_r*f_r(x,y,z,3)/2.d0 
           end if
#ifndef SINGLEFLUID
           scalar(x,y,z) = scalar(x,y,z) + sum(N(x,y,z)%n_b(:)*g*cz)*amass_b*omega_b !- tau_b*f_b(x,y,z,3)/2.d0
#endif
#ifndef NOSURFACTANT
           scalar(x,y,z) = scalar(x,y,z) + sum(N(x,y,z)%n_s(:)*g*cz)*amass_s*omega_s !- tau_s*f_s(x,y,z,3)/2.d0
           scalar(x,y,z) = scalar(x,y,z) / &
                & max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*omega_r &
                + sum(N(x,y,z)%n_b(:)*g)*amass_b*omega_b + sum(N(x,y,z)%n_s(:)*g)*amass_s*omega_s) )
#else 
#ifndef SINGLEFLUID
           scalar(x,y,z) = scalar(x,y,z) / &
                & max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*omega_r  + sum(N(x,y,z)%n_b(:)*g)*amass_b*omega_b ))
#else
           scalar(x,y,z) = scalar(x,y,z) / max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*omega_r ))
#endif
#endif
           
           ! force extra-term
           tx = x + ccoords(1)*nx
           ty = y + ccoords(2)*ny
           tz = z + ccoords(3)*nz
           if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
                (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
                (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
              if (N(x,y,z)%rock_state==0.0d0) then
                 ! external force term, in the zone where the
                 ! external force works
                 scalar(x,y,z) = scalar(x,y,z) - g_accn/2.0d0
              end if
           endif
        else
           scalar(x,y,z) = 0.0d0
        end if

      end do
    end do
  end do
  CALL dump_scalar(scalar,'vel')
  deallocate(scalar)
end subroutine dump_vel

!> Calculates the oil density at each site, then calls
!> \c dump_scalar.
subroutine dump_od(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: oil
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(oil(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for od output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (N(x,y,z)%rock_state == 0.) then
          oil(x,y,z)=amass_r*sum(N(x,y,z)%n_r(:)*g)
        else
          oil(x,y,z)=0.0
        endif
      end do
    end do
  end do
  CALL dump_scalar(oil,'od')
  deallocate(oil)
end subroutine dump_od

#ifndef SINGLEFLUID
!> 
!> Calculates water density at each site, then calls
!> \c dump_scalar.
subroutine dump_wd(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: water
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(water(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for wd output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (N(x,y,z)%rock_state == 0.) then
          water(x,y,z)=amass_b*sum(N(x,y,z)%n_b(:)*g)
        else
          water(x,y,z)=0.0
        endif
      end do
    end do
  end do
  CALL dump_scalar(water,'wd')
  deallocate(water)
end subroutine dump_wd
#endif

!> Calculates the Z velocity of each species at each site, then calls
!> \c dump_3scalar.
subroutine dump_flo(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: oil,water,surf
  integer :: ierror,x,y,z,nxi,nyi,nzi,tz,tx,ty
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_s
#endif

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for flo output"&
          &,.true.)
     return
  end if

  allocate(oil(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for flo oil output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        ! Get Shan Chen forces
        ! only calculate forces on fluid at non-rock sites
        if (N(x,y,z)%rock_state == 0.d0) then
#ifdef SINGLEFLUID
          CALL lb3d_calc_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if
        oil(x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_r(:)*g)))-f_r(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if (N(x,y,z)%rock_state==0.0d0) then
            ! external force term, in the zone where the external force works
            oil(x,y,z) = oil(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  CALL dump_scalar(oil,'flooil')
  deallocate(oil)

#ifndef SINGLEFLUID
  allocate(water(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for flo water output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        water(x,y,z)=sum(N(x,y,z)%n_b(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_b(:)*g)))-f_b(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if (N(x,y,z)%rock_state==0.0d0) then
            ! external force term, in the zone where the external force works
            water(x,y,z) = water(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  CALL dump_scalar(water,'flowater')
  deallocate(water)
#endif
#ifndef NOSURFACTANT
  allocate(surf(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for flo surf output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        surf(x,y,z)=sum(N(x,y,z)%n_s(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_s(:)*g)))-f_s(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if (N(x,y,z)%rock_state==0.0d0) then
            ! external force term, in the zone where the external force works
            surf(x,y,z) = surf(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  CALL dump_scalar(surf,'flosurf')
  deallocate(surf)
#endif

end subroutine dump_flo

!>Calculates the nett flow direction at each site, then calls \c dump_vector.
!>
!>Added by Nelido, edited by Jens 26.06.02
!>Edited by FrankR 19.07.07 to calculate the averaged velocity, not the momentum
subroutine dump_arrows(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real(kind=rk), dimension(:,:,:,:), allocatable :: v
  integer :: ierror,x,y,z,nxi,nyi,nzi,tx,ty,tz
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real(kind=rk),dimension(:,:,:,:),allocatable :: f_s
#endif

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for arrows output"&
          &,.true.)
     return
  end if
  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for arrows output"&
          &,.true.)
     return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
      ! Get Shan Chen forces
      ! only calculate forces on fluid at non-rock sites
        if (N(x,y,z)%rock_state == 0.d0) then
#ifdef SINGLEFLUID
          CALL lb3d_calc_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lb3d_calc_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if

#ifndef SINGLEFLUID
#ifndef NOSURFACTANT
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r&
          &-(f_r(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cx)*amass_b&
          &-(f_b(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cx)*amass_s&
          &-(f_s(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r&
          &-(f_r(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cy)*amass_b&
          &-(f_b(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cy)*amass_s&
          &-(f_s(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r&
          &-(f_r(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cz)*amass_b&
          &-(f_b(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cz)*amass_s&
          &-(f_s(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)

        v(1,x,y,z)=v(1,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
        v(2,x,y,z)=v(2,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
        v(3,x,y,z)=v(3,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
#else
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r&
          &-(f_r(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cx)*amass_b&
          &-(f_b(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r&
          &-(f_r(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cy)*amass_b&
          &-(f_b(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r&
          &-(f_r(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cz)*amass_b&
          &-(f_b(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)

        v(1,x,y,z)=v(1,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
        v(2,x,y,z)=v(2,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
        v(3,x,y,z)=v(3,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
#endif
#else
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r

        v(1,x,y,z)=(v(1,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,1)/2.0d0
        v(2,x,y,z)=(v(2,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,2)/2.0d0
        v(3,x,y,z)=(v(3,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,3)/2.0d0
#endif

        ! force extra-term

        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        tz=z+ccoords(3)*nz
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if (N(x,y,z)%rock_state==0.0d0) then
            ! external force term, in the zone where the
            ! external force works
            v(1,x,y,z) = v(1,x,y,z)-g_accn_x/2.0d0
            v(2,x,y,z) = v(2,x,y,z)-g_accn_y/2.0d0
            v(3,x,y,z) = v(3,x,y,z)-g_accn/2.0d0
          end if
        endif

      end do
    end do
  end do
  if ( check_dump_now(sci_arrows, n_sci_arrows) ) then
    CALL dump_vector(v,'arr')
  endif

  if ( check_dump_now(sci_velocities, n_sci_velocities) ) then
    CALL dump_scalar(v(1,:,:,:),'velx')
    CALL dump_scalar(v(2,:,:,:),'vely')
    CALL dump_scalar(v(3,:,:,:),'velz')
  endif

  deallocate(v)
end subroutine dump_arrows

!>  dumps rock_state of whole system
subroutine dump_rock(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi
  real(kind=rk) :: rockfloat

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer for rock output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
       ! For XDRROCKWET this should be the inverse of read_rock
#ifdef XDRROCKWET
        if (N(x,y,z)%rock_state > 0.0) then
          if (N(x,y,z)%n_r(restvec) > 0.0) then
            rockfloat = N(x,y,z)%n_r(restvec)
#ifdef WALLCONST
            rockfloat = rockfloat*tau_r
#endif
#ifndef SINGLEFLUID
          elseif (N(x,y,z)%n_b(restvec) > 0.0) then
            rockfloat = -N(x,y,z)%n_b(restvec)
#ifdef WALLCONST
            rockfloat = rockfloat*tau_b
#endif
#endif
          else
            rockfloat = 0.0 ! Does this ever happen? I hope not...
          endif
          rockfloat = rockfloat + 5.0
          scalar(x,y,z) = rockfloat
        else ! rockstate !> 0.0
          scalar(x,y,z) = N(x,y,z)%rock_state
        endif
#else
        ! No XDRROCKWET
        scalar(x,y,z) = N(x,y,z)%rock_state
#endif
      end do
    end do
  end do

  CALL dump_scalar(scalar,'rock')

  deallocate(scalar)
end subroutine dump_rock

!>Calculates total fluid velocity at each site, then calls \c dump_vector.
!>
!>Based on the kinetic theory definition for gas mixtures
!>(cf. Chapman & Cowling, ``The mathematical theory of non-uniform
!>gases." (CUP: Cambridge, 1970, 3rd ed.)), which for a discrete
!>velocity set numbered by \f$k\f$ and species masses \f$m^{\alpha}\f$, reads:
!>\f[\mathbf{u}=\frac{\rho\mathbf{u}}{\rho}
!>\equiv
!>\frac{\sum_{\alpha}m^{\alpha}\sum_k\mathbf{c}_k n_k^{\alpha}}
!>{\sum_{\alpha}m^{\alpha}\sum_k n_k^{\alpha}}\f]
!>where \f$n_k^{alpha}\f$ is the number density, of the real var
!>\c N(x,y,z)%n_\f$\alpha\f$\c (k), \f$u\f$=total fluid velocity (the unknown),
!>\f$f\f$=single-particle d.f.
!>
!>Added by N. Gonzalez-Segredo 26.09.03.
subroutine tot_vel(N,u,switch,verb)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:,:) :: u
  integer :: ierror,x,y,z,nxi,nyi,nzi
  real(kind=rk)  :: dens
  real(kind=rk), parameter :: eps = 1.e-20
  logical :: switch, verb
  logical :: rep = .false.
  character(len=256) :: msgstr

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (N(x,y,z)%rock_state/=0.0_rk) then
          u(:,x,y,z) = 0.0_rk
        else
#ifndef NOSURFACTANT
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g) + &
                 amass_b*sum(N(x,y,z)%n_b(:)*g) + &
                 amass_s*sum(N(x,y,z)%n_s(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lb3d_io.F90: tot_vel(): dens < ',F16.10)") eps
            CALL log_msg(trim(msgstr),.false.)
            rep = .true.
          endif
          u(1,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cx) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cx) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cx)
          u(1,x,y,z)= u(1,x,y,z)/dens
          u(2,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cy) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cy) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cy)
          u(2,x,y,z)= u(2,x,y,z)/dens
          u(3,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cz) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cz) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cz)
          u(3,x,y,z)=	u(3,x,y,z)/dens
#else
#ifndef SINGLEFLUID
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g) + &
                 amass_b*sum(N(x,y,z)%n_b(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lb3d_io.F90: tot_vel(): dens < ',F16.10)") eps
            CALL log_msg(trim(msgstr),.false.)
            rep = .true.
          endif

          u(1,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cx) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cx)
          u(2,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cy) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cy)
          u(3,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cz) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cz)
#else
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lb3d_io.F90: tot_vel(): dens < ',F16.10)") eps
            CALL log_msg(trim(msgstr),.false.)
            rep = .true.
          endif
          u(1,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cx)
          u(2,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cy)
          u(3,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cz)
#endif
          if (dens.lt.eps) then
            u(1,x,y,z) = 0.
            u(2,x,y,z) = 0.
            u(3,x,y,z) = 0.
          else
            u(1,x,y,z) = u(1,x,y,z)/dens
            u(2,x,y,z) = u(2,x,y,z)/dens
            u(3,x,y,z) = u(3,x,y,z)/dens
          endif
#endif
        end if
      end do
    end do
  end do
  if(switch .eqv. .true.) then
    CALL dump_vector(u,'tvel')
  endif
end subroutine tot_vel
!> \}

!> \{
!> \name Array-specific stubroutines - Format-specific routines
!>
!> For an array of type \c foo, there will be a routine called \c
!> dump_foo(), which will in turn call \c dump_foo_all(), \c
!> dump_foo_bin(), or some other routine \c dump_foo_bar(), for file
!> format \c bar.

!> Dumps a floating point scalar field in a format depending on \c
!> dump_format.
!>
!> \param[in,out] scalar local portion of the scalar field
!>
!> \param[in] name string to be used to generate the filename from
subroutine dump_scalar(scalar,name)
    real(kind=rk), dimension(1:,1:,1:) :: scalar
    character(len=*)            :: name

    select case(trim(dump_format))
    case ('hdf')
#ifdef USEHDF
       CALL dump_scalar_phdf5(scalar, name)
#else
       CALL log_msg("HDF5 support switched off",.false.)
#endif
    case ('mpi')
#ifdef USEXDRF
       CALL dump_scalar_xdr_parallel(scalar, name)
#else
       CALL log_msg("XDRF support switched off",.false.)
#endif
    case ('all')
       CALL dump_scalar_all(scalar, name)
    case ('xdr')
#ifdef USEXDRF
       CALL dump_scalar_xdr(scalar, name)
#else
       CALL log_msg("XDRF support switched off",.false.)
#endif
    case ('vtk')
#ifdef USEXDRF
       CALL dump_scalar_vtk(scalar, name)
#else
       CALL log_msg("XDRF/VTK support switched off",.false.)
#endif
    case ('bin')
       CALL dump_scalar_bin(scalar, name)
    case default
       call error('unknown value: dump_format="'//trim(dump_format)//'"')
    end select
end subroutine dump_scalar

!> Dumps an integer scalar field in a format depending on \c
!> dump_format.
!>
!> \param[in,out] iscalar local portion of the scalar field
!>
!> \param[in] name string to be used to generate the filename from
subroutine dump_iscalar(iscalar,name)
    integer,dimension(1:,1:,1:) :: iscalar
    character(len=*)            :: name

    select case (trim(dump_format))
    case ('hdf')
#ifdef USEHDF
       call dump_iscalar_phdf5(iscalar,name)
#else
       call log_msg("HDF5 support switched off",.false.)
#endif
    case ('mpi')
#ifdef USEXDRF
       call error('dump_iscalar_xdr_parallel() not implemented yet')
!!$          call dump_iscalar_xdr_parallel(iscalar, name)
#else
       CALL log_msg("XDRF support switched off",.false.)
#endif
    case ('all')
       call error('dump_iscalar_all() not implemented yet')
!!$          call dump_iscalar_all(iscalar, name)
    case ('xdr')
#ifdef USEXDRF
       call dump_iscalar_xdr(iscalar,name)
#else
       call log_msg("XDRF support switched off",.false.)
#endif
    case ('vtk')
#ifdef USEXDRF
       call error('dump_iscalar_vtk() not implemented yet')
!!$          call dump_iscalar_vtk(iscalar, name)
#else
       CALL log_msg("XDRF/VTK support switched off",.false.)
#endif
    case ('bin')
       call error('dump_iscalar_bin() not implemented yet')
!!$          call dump_iscalar_bin(iscalar, name)
    case default
       call error('unknown value: dump_format="'//trim(dump_format)//'"')
    end select
end subroutine dump_iscalar

!> Dumps a scalar value in ASCII "all" form.
subroutine dump_scalar_all(scalar,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:) :: scalar
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi

        nxi = size(scalar,1)
        nyi = size(scalar,2)
        nzi = size(scalar,3)

  call lbe_make_filename_output(filename,trim(name),'.all',nt)
  
  open(10,file=filename)

  do z=1,nzi
   do y=1,nyi
    do x=1,nxi
             if (dump_double) then
    write(10,'(3i3,f12.8)')         &
      ccoords(1)*nxi+x,       &
      ccoords(2)*nyi+y,       &
      ccoords(3)*nzi+z,       &
        scalar(x,y,z)
             else
    write(10,'(3i3,f8.4)')          &
      ccoords(1)*nxi+x,       &
      ccoords(2)*nyi+y,       &
      ccoords(3)*nzi+z,       &
        scalar(x,y,z)
             end if
    end do
   end do
  end do

  close(10)
  !write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
#ifdef USEXDRF
  if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,1)
#else
    CALL log_msg("XDRF is disabled, can't dump avs field",.false.)
#endif
    
end subroutine dump_scalar_all

!> Dumps a scalar in binary format.
subroutine dump_scalar_bin(scalar,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:) :: scalar
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: x,y,z,nxi,nyi,nzi,ierror
        real*4, dimension(:,:,:),allocatable :: scalar2

        nxi = size(scalar,1)
        nyi = size(scalar,2)
        nzi = size(scalar,3)
 
  call lbe_make_filename_output(filename,trim(name),'.bin',nt)
  
  open(10,file=filename,form="unformatted")
          if (dump_double) then
          write(10) scalar(:,:,:)
          else

             ! FIXME: Very unelegegant way of doing this, 
             ! but conversion functions did not work
                allocate(scalar2(nxi,nyi,nzi),stat=ierror)
                 scalar2 = scalar
          write(10) scalar2(:,:,:)
                deallocate(scalar2)
          end if
  close(10)
        
  !write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'

#ifdef USEXDRF
  if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,1)
#else
  CALL log_msg("XDRF is disabled, can't dump avs field",.false.)
#endif
end subroutine dump_scalar_bin

!> Calls \c dump_vector_all or \c dump_vector_bin depending on the value
!> of the variable \c dump_format.
subroutine dump_vector(vector,name)
    real(kind=rk),intent(inout), dimension(1:,1:,1:,1:) :: vector
    character(len=*) :: name
    character(len=256) :: filename

    if (index(dump_format,'all').gt.0) then
       CALL dump_vector_all(vector,name)

    elseif(index(dump_format,'hdf').gt.0)then
#ifdef USEHDF
       CALL lbe_make_filename_output(filename, name, '.h5', nt)
       CALL dump_vector_phdf5(vector,filename)
#else
       CALL log_msg("HDF5 support switched off",.false.)
#endif

    elseif (index(dump_format,'mpi').gt.0) then
#ifdef USEXDRF
       CALL dump_vector_xdr_parallel(vector,name)
#else
       CALL log_msg("XDRF support switched off",.false.)
#endif

    elseif (index(dump_format,'xdr').gt.0) then
#ifdef USEXDRF
       CALL dump_vector_xdr(vector,name)
#else
       CALL log_msg("XDRF support switched off",.false.)
#endif

    else
       CALL dump_vector_bin(vector,name)
    endif
end subroutine dump_vector

!> Dumps a vector field in ASCII "all" form.
subroutine dump_vector_all(vector,name)
	implicit none
	real(kind=rk), dimension(1:,1:,1:,1:) :: vector
	character(LEN=*) :: name
	character(LEN=256) :: filename
	integer :: x,y,z,nxi,nyi,nzi

        nxi = size(vector,2)
        nyi = size(vector,3)
        nzi = size(vector,4)

	call lbe_make_filename_output(filename,trim(name),'.all',nt)
	open(10,file=filename)

	do z=1,nzi
	 do y=1,nyi
	  do x=1,nxi
             if (dump_double) then
		write(10,'(3i3,3f12.8)')				&
			ccoords(1)*nxi+x,				&
			ccoords(2)*nyi+y,				&
			ccoords(3)*nzi+z,				&
	  		vector(:,x,y,z)
             else
		write(10,'(3i3,3f8.4)')				&
			ccoords(1)*nxi+x,				&
			ccoords(2)*nyi+y,				&
			ccoords(3)*nzi+z,				&
	  		vector(:,x,y,z)
             end if
	  end do
	 end do
	end do

	close(10)
	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
#ifdef USEXDRF
        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
#else
        CALL log_msg("XDRF is disabled, can't dump avs field",.false.)
#endif

end subroutine dump_vector_all

!> Dump an array of vectors in binary format.
subroutine dump_vector_bin(vector,name)
	implicit none
	real(kind=rk), dimension(1:,1:,1:,1:) :: vector
        real*4, dimension(:,:,:,:),allocatable :: vector2 
	character(LEN=*) :: name
	character(LEN=256) :: filename
	integer :: x,y,z,nxi,nyi,nzi,ierror

        nxi = size(vector,2)
        nyi = size(vector,3)
        nzi = size(vector,4)

	call lbe_make_filename_output(filename,trim(name),'.bin',nt)
	
	open(10,file=filename,form="unformatted")

          if (dump_double) then
	     write(10) vector(:,:,:,:)
          else
             allocate(vector2(size(vector,1),nxi,nyi,nzi),stat=ierror)
             vector2 = vector
	     write(10) vector2(:,:,:,:)
             deallocate(vector2)
          end if

	close(10)
	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
      
#ifdef USEXDRF
        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
#else
        CALL log_msg("XDRF is disabled, can't dump avs field",.false.)
#endif
end subroutine dump_vector_bin
!> \}

!> \{
!> \name lbe_pressure routines
!>
!>T hese routines:
!>     [1] compute the components of the pressure tensor
!>     [2] compute Pzz-Pxx along a line perpendicular to
!>     a planar interface. (Not used.)
!>
!>Added by NGS. Modified by NGS on 26.09.03.

!> Calculates all the components of the pressure tensor.
!> 's' decides whether or not to include scalar pressure calculation.
subroutine dump_pressure(N,s)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: pxx, pyy, pzz,    &
                                           pxy, pyz, pxz,    &
                                           scalarpressure
  real(kind=rk), dimension(3,3) :: sum_psixc_r
  real(kind=rk), parameter :: omega = 0.25
  integer :: ierror,x,y,z,s,xc,yc,zc,i,j,nxi,nyi,nzi
  real(kind=rk) :: psixc_r
  real(kind=rk) :: psix_r
  real(kind=rk), dimension(:,:,:,:), allocatable :: u !Modif. NGS
#ifndef SINGLEFLUID
  real(kind=rk), dimension(3,3) :: sum_psixc_b, p_br
  real(kind=rk) :: psixc_b, psix_b
#ifndef NOSURFACTANT
  real(kind=rk), dimension(3,3) :: p_bs, p_ss, p_rs, sum_psixc_s
  real(kind=rk) :: psixc_s, psix_s
#endif
#endif

  ! !DEBUG
  !   integer, dimension(5) :: sample_point
  !   sample_point(1)=size(N,1)/8
  !   sample_point(2)=size(N,1)/4
  !   sample_point(3)=size(N,1)/2
  !   sample_point(4)=size(N,1)*(1.0-1.0/4.0)
  !   sample_point(5)=size(N,1)*(1.0-1.0/8.0)
  ! !END_DEBUG

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(u(3,1:nxi,1:nyi,1:nzi),stat=ierror) !Modif. NGS
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate total velocity real array u(:,:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxx(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_xx(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pyy(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_yy(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pzz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_zz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxy(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_xy(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pyz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_yz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer p_xz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  if(s == droplet) then
    allocate(scalarpressure(1:nxi,1:nyi,1:nzi),stat=ierror)
    if (ierror .ne. 0) then
      CALL log_msg("WARNING: unable to allocate scalar buffer scalarpressure(:,:,:)",.true.)
      return  ! Inability to write output is not quite fatal.
    end if
  endif

  ! Health measure:
  pxx = 0.
  pyy = 0.
  pzz = 0.
  pxy = 0.
  pyz = 0.
  pxz = 0.
  if (s == droplet) then
    scalarpressure = 0.
  end if

  CALL tot_vel(N,u,.false.,.false.)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi

        ! INTER-COMPONENT-INDEPENDENT PART
        ! OF THE PRESSURE TENSOR.
        ! DEFINED IN TERMS OF THE BARICENTRIC
        ! MICROSCOPIC VELOCITY    --NGS
        pxx(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#ifndef SINGLEFLUID
        pxx(x,y,z) = pxx(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = pyy(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = pzz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = pxy(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = pxz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = pyz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#ifndef NOSURFACTANT
        pxx(x,y,z) = pxx(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = pyy(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = pzz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = pxy(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = pxz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = pyz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#endif
#endif

        !DEBUG
        !do i=1,5
        !   if(x==sample_point(i) .and.				  &
        !   y==nyi/2 .and. z==nzi/2) then
        !      print*,'At (',x,',',y,',',z,'):'
        !      print*,'  pxz_kinetic = ', pxz(x,y,z)
        !      print*,'  u_x =', u(1,x,y,z)
        !      print*,'  u_z =', u(3,x,y,z)
        !      print*,'  amass_r * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_r * sum(N(x,y,z)%n_r(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'  amass_b * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_b * sum(N(x,y,z)%n_b(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'  amass_s * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_s * sum(N(x,y,z)%n_s(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'0 =? sum_alpha(mass * sum((c_k - u) * n)) = (',&
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cx-u(1,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cx-u(1,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cx-u(1,x,y,z))), &
        !      ',',						    &
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cy-u(2,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cy-u(2,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cy-u(2,x,y,z))), &
        !      ',',						    &
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cz-u(3,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cz-u(3,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cz-u(3,x,y,z))), &
        !      ')'
        !   endif
        !enddo
        !END

        ! INTER-COMPONENT-DEPENDENT PART
        ! OF THE PRESSURE TENSOR
        ! DOES NOT INCLUDE VELOCITY SINCE IT'S NOT A 
        ! KINETIC CONTRIBUTION    --NGS

        ! Set the contributions to zero initially
        !
        ! They are array(3x3):
        sum_psixc_r = 0.
        psix_r = sum( N(x,y,z)%n_r(:)*g(:) )
#ifndef SINGLEFLUID
        sum_psixc_b = 0.
        psix_b = sum( N(x,y,z)%n_b(:)*g(:) )
#ifndef NOSURFACTANT
        sum_psixc_s = 0.
        psix_s = sum( N(x,y,z)%n_s(:)*g(:) )
#endif
#endif

        select case(psifunc)
          case (0)
            ! Original code: psi=number density
            ! plus funny clipping routine
            psix_r = min(1., real(psix_r))
#ifndef SINGLEFLUID
            psix_b = min(1., real(psix_b))
#ifndef NOSURFACTANT
            psix_s = min(1., real(psix_s))
#endif
#endif
          case (1)
            ! psi = n, with no clipping
            ! nothing further to do.
          case (2)
            ! psi = 1 - exp(-n)
            ! no clipping
            psix_r = 1 - exp(-psix_r)
#ifndef SINGLEFLUID
            psix_b = 1 - exp(-psix_b)
#ifndef NOSURFACTANT
            psix_s = 1 - exp(-psix_s)
#endif
#endif
          case default
            CALL log_msg("ERROR: Unknown psi functional, aborting...",.false.)
            CALL Abend
        end select

        do i=1,nnonrest
          xc = x + cx(i)
          yc = y + cy(i)
          zc = z + cz(i)
          ! psi_b = no density of blue particles at the site (ixa,iya,iza),
          ! calculated by summing all those advecting from that site.

          psixc_r = sum( N(xc,yc,zc)%n_r(:)*g(:) )
#ifndef SINGLEFLUID
          psixc_b = sum( N(xc,yc,zc)%n_b(:)*g(:) )
#ifndef NOSURFACTANT
          psixc_s = sum( N(xc,yc,zc)%n_s(:)*g(:) )
#endif
#endif

          select case(psifunc)
            case (0)
              ! Original code: psi=number density
              ! plus funny clipping routine
              psixc_r = min(1., real(psixc_r))
#ifndef SINGLEFLUID
              psixc_b = min(1., real(psixc_b))
#ifndef NOSURFACTANT
              psixc_s = min(1., real(psixc_s))
#endif
#endif
            case (1)
              ! psi = n, with no clipping
              ! nothing further to do.
            case (2)
              ! psi = 1 - exp(-n)
              ! no clipping
              psixc_r = 1 - exp(-psixc_r)
#ifndef SINGLEFLUID
              psixc_b = 1 - exp(-psixc_b)
#ifndef NOSURFACTANT
              psixc_s = 1 - exp(-psixc_s)
#endif
#endif
            case default
              CALL log_msg("ERROR: Unknown psi functional, aborting...",.false.)
              CALL Abend
          end select

          sum_psixc_r(1,1) = sum_psixc_r(1,1) + psixc_r * g(i) * cx(i) * cx(i)
          sum_psixc_r(2,2) = sum_psixc_r(2,2) + psixc_r * g(i) * cy(i) * cy(i)
          sum_psixc_r(3,3) = sum_psixc_r(3,3) + psixc_r * g(i) * cz(i) * cz(i)
          sum_psixc_r(1,2) = sum_psixc_r(1,2) + psixc_r * g(i) * cx(i) * cy(i)
          sum_psixc_r(1,3) = sum_psixc_r(1,3) + psixc_r * g(i) * cx(i) * cz(i)
          sum_psixc_r(2,3) = sum_psixc_r(2,3) + psixc_r * g(i) * cy(i) * cz(i)
#ifndef SINGLEFLUID
          sum_psixc_b(1,1) = sum_psixc_b(1,1) + psixc_b * g(i) * cx(i) * cx(i)
          sum_psixc_b(2,2) = sum_psixc_b(2,2) + psixc_b * g(i) * cy(i) * cy(i)
          sum_psixc_b(3,3) = sum_psixc_b(3,3) + psixc_b * g(i) * cz(i) * cz(i)
          sum_psixc_b(1,2) = sum_psixc_b(1,2) + psixc_b * g(i) * cx(i) * cy(i)
          sum_psixc_b(1,3) = sum_psixc_b(1,3) + psixc_b * g(i) * cx(i) * cz(i)
          sum_psixc_b(2,3) = sum_psixc_b(2,3) + psixc_b * g(i) * cy(i) * cz(i)
#ifndef NOSURFACTANT
          sum_psixc_s(1,1) = sum_psixc_s(1,1) + psixc_s * g(i) * cx(i) * cx(i)
          sum_psixc_s(2,2) = sum_psixc_s(2,2) + psixc_s * g(i) * cy(i) * cy(i)
          sum_psixc_s(3,3) = sum_psixc_s(3,3) + psixc_s * g(i) * cz(i) * cz(i)
          sum_psixc_s(1,2) = sum_psixc_s(1,2) + psixc_s * g(i) * cx(i) * cy(i)
          sum_psixc_s(1,3) = sum_psixc_s(1,3) + psixc_s * g(i) * cx(i) * cz(i)
          sum_psixc_s(2,3) = sum_psixc_s(2,3) + psixc_s * g(i) * cy(i) * cz(i)
#endif
#endif
        end do !Next velocity i

        ! PRESSURE TENSOR IS DECOMPOSED INTO INTERACTIONS:
#ifndef SINGLEFLUID
        p_br = g_br * (psix_r * sum_psixc_b + psix_b * sum_psixc_r)
#ifndef NOSURFACTANT
        p_bs = g_bs * (psix_b * sum_psixc_s + psix_s * sum_psixc_b)
        p_ss = g_ss * (psix_s * sum_psixc_s) * 2.0
        p_rs = g_bs * (psix_r * sum_psixc_s + psix_s * sum_psixc_r)
#endif
#endif
        ! g_rs = g_bs

        ! NOW, FULL PRESSURE TENSOR
#ifndef SINGLEFLUID
        pxx(x,y,z) = pxx(x,y,z) + omega * (p_br(1,1))
        pyy(x,y,z) = pyy(x,y,z) + omega * (p_br(2,2))
        pzz(x,y,z) = pzz(x,y,z) + omega * (p_br(3,3))
        pxy(x,y,z) = pxy(x,y,z) + omega * (p_br(1,2))
        pxz(x,y,z) = pxz(x,y,z) + omega * (p_br(1,3))
        pyz(x,y,z) = pyz(x,y,z) + omega * (p_br(2,3))
#ifndef NOSURFACTANT
        pxx(x,y,z) = pxx(x,y,z) +                                   &
                    omega * (p_bs(1,1) + p_ss(1,1) + p_rs(1,1))
        pyy(x,y,z) = pyy(x,y,z) +                                   &
                    omega * (p_bs(2,2) + p_ss(2,2) + p_rs(2,2))
        pzz(x,y,z) = pzz(x,y,z) +                                   &
                    omega * (p_bs(3,3) + p_ss(3,3) + p_rs(3,3))
        pxy(x,y,z) = pxy(x,y,z) +                                   &
                    omega * (p_bs(1,2) + p_ss(1,2) + p_rs(1,2))
        pxz(x,y,z) = pxz(x,y,z) +                                   &
                    omega * (p_bs(1,3) + p_ss(1,3) + p_rs(1,3))
        pyz(x,y,z) = pyz(x,y,z) +                                   &
                    omega * (p_bs(2,3) + p_ss(2,3) + p_rs(2,3))
#endif
#endif
        if(s == droplet) then
          scalarpressure(x,y,z) = (pxx(x,y,z) + pyy(x,y,z) + pzz(x,y,z))/3.
        endif

      end do
    end do
  end do

  if(s == droplet) then
    CALL dump_scalar(scalarpressure,'scp')
  else if(s == nondroplet) then
    CALL dump_scalar(pxx,'pxx')
    CALL dump_scalar(pyy,'pyy')
    CALL dump_scalar(pzz,'pzz')
    CALL dump_scalar(pxy,'pxy')
    CALL dump_scalar(pyz,'pyz')
    CALL dump_scalar(pxz,'pxz')
  end if

  !DEBUG
  !do i=1,5
  !   print*,'At (',sample_point(i),',',nyi/2,',',nzi/2, &
  !   '):'
  !   print*,'  pxz_virial = ', &
  !   omega * (p_br(1,3) + p_bs(1,3) + p_ss(1,3) + p_rs(1,3))
  !   print*,'  pxz_kinetic+virial =', &
  !   pxz(sample_point(i),nyi/2,nzi/2)
  !enddo
  !END

  deallocate(pxx)
  deallocate(pyy)
  deallocate(pzz)
  deallocate(pxy)
  deallocate(pyz)
  deallocate(pxz)
  if(s == droplet) then
      deallocate(scalarpressure)
  endif

  deallocate(u)
end subroutine dump_pressure

!> Calculates \f$ p(tz)=(P_zz - P_xx)(tz) \f$ and then calls
!> \c dump_scalar to write it onto file in order to compute:
!> surface tension = \f$ \int pressure(z) dz \f$
!> at postprocessing time.
!>
!> s decides whether including scalar pressure calculation
subroutine dump_popul_zx(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(:,:,:), allocatable :: populz, populx
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(populz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer populz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if

  allocate(populx(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    CALL log_msg("WARNING: unable to allocate scalar buffer populx(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if


  ! Following is unnecessary if not using accumulators.
  do z=1,nzi
    do y=1,nyi
    do x=1,nxi
        populz(x,y,z) = 0.
        populx(x,y,z) = 0.
    enddo
    enddo
  enddo

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi

        populz(x,y,z) =   N(x,y,z)%n_r(13) + N(x,y,z)%n_r(15) + &
                          N(x,y,z)%n_r(9)  + N(x,y,z)%n_r(17) + &
                          N(x,y,z)%n_r(5)  + &
                          N(x,y,z)%n_r(14) + N(x,y,z)%n_r(16) + &
                          N(x,y,z)%n_r(10) + N(x,y,z)%n_r(18) + &
                          N(x,y,z)%n_r(6)
#ifndef SINGLEFLUID
        populz(x,y,z) =   populz(x,y,z)    + &
                          N(x,y,z)%n_b(13) + N(x,y,z)%n_b(15) + &
                          N(x,y,z)%n_b(9)  + N(x,y,z)%n_b(17) + &
                          N(x,y,z)%n_b(5)  + &
                          N(x,y,z)%n_b(14) + N(x,y,z)%n_b(16) + &
                          N(x,y,z)%n_b(10) + N(x,y,z)%n_b(18) + &
                          N(x,y,z)%n_b(6)
#endif
        populx(x,y,z) =   N(x,y,z)%n_r(7)  + N(x,y,z)%n_r(10) + &
                          N(x,y,z)%n_r(8)  + N(x,y,z)%n_r(9)  + &
                          N(x,y,z)%n_r(1)  + &
                          N(x,y,z)%n_r(11) + N(x,y,z)%n_r(14) + &
                          N(x,y,z)%n_r(12) + N(x,y,z)%n_r(13) + &
                          N(x,y,z)%n_r(2)
#ifndef SINGLEFLUID
        populx(x,y,z) = populx(x,y,z)      + &
                          N(x,y,z)%n_b(7)  + N(x,y,z)%n_b(10) + &
                          N(x,y,z)%n_b(8)  + N(x,y,z)%n_b(9)  + &
                          N(x,y,z)%n_b(1)  + &
                          N(x,y,z)%n_b(11) + N(x,y,z)%n_b(14) + &
                          N(x,y,z)%n_b(12) + N(x,y,z)%n_b(13) + &
                          N(x,y,z)%n_b(2)
#endif
      end do
    end do
  end do

  CALL dump_scalar(populz,'popz')
  CALL dump_scalar(populx,'popx')

  deallocate(populz)
  deallocate(populx)
end subroutine dump_popul_zx
!> \}

!> \{
!> \name PROFILE output

!!$!> dump profiles
!!$!>
!!$!> \param[in] N local lattice chunk with halo of depth \c halo_extent
!!$!>
!!$!> Dumps profiles along the coordinate axis and along and
!!$!> perpendicular to the body force vector (only if this itself is
!!$!> not parallel to one of the coordinate axis). Beside the
!!$!> self-explaining information in the profile header each file
!!$!> consist of one line per profile position containing profile
!!$!> position, averaged fluid velocity vector, averaged fluid
!!$!> density, averaged fluid site concentration, and rock site
!!$!> concentration. If compiled with MD additional columns containing
!!$!> averaged particle velocity, averaged particle center
!!$!> concentration, and averaged particle site concentration are
!!$!> printed. Depending on n_sci_profile and n_sci_profile_dump the
!!$!> dumped result is the time-average of several samples.
!!$subroutine dump_profile(N)
!!$  type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
!!$  integer,parameter :: profile_file_unit=12
!!$  character(len=256) profile_file_name
!!$  integer i, ii, k, l, lp(3), lx(3), opp_lp(3), u, v, w, ierror
!!$  real(kind=rk) :: rp(3), total_sites,weight
!!$  real(kind=rk) :: ifs, irs, nfs, nrs
!!$  real(kind=rk) :: ivf(3), vf(3), imf(3), mf(3), irf, rhof
!!$#ifndef SINGLEFLUID
!!$  real(kind=rk) :: iv_r(3), v_r(3), iv_b(3), v_b(3)
!!$  real(kind=rk) :: im_r(3), m_r(3), im_b(3), m_b(3)
!!$  real(kind=rk) :: ir_r, ir_b, rho_r, rho_b
!!$#endif
!!$#ifdef MD
!!$  real(kind=rk) :: ipc,ips,ivp(3),vp(3),npc,nps
!!$#endif
!!$
!!$  ! first call of this subroutine: perform initialization
!!$  if (.not.allocated(prf)) CALL setup_dump_profile()
!!$
!!$#ifdef MD
!!$  ! calculate on-site particle center density and velocity from
!!$  ! off-site particle positions
!!$  Nvp = 0.0_rk
!!$  Nnpc = 0.0_rk
!!$  i = atompnt
!!$  particles: do ii = 1,nlocal+nother
!!$    ! Iterate through all possible periodic images; this is
!!$    ! surely not efficient but most safe.
!!$    pbc_x: do u=-1,1
!!$      pbc_y: do v=-1,1
!!$        pbc_z: do w=-1,1
!!$          lx = floor(P(i)%x)+(/u*tnx,v*tny,w*tnz/)
!!$          lattice_points: do l=1,n_lp
!!$            ! lattice point to look at now, these are still
!!$            ! global coordinates
!!$            lp(:) = lx(:) + lp_sur(:,l)
!!$
!!$            ! transform to local coordinates, clipping to
!!$            ! local chunk follows
!!$            lp = lp-(start-1)
!!$
!!$            ! only calculate data for real nodes and halo
!!$            ! nodes (halo extent 1)
!!$            if (any(lp<0).or.any(lp>(/nx,ny,nz/)+1)) cycle
!!$
!!$            ! weight for each lattice point is the volume of the
!!$            ! hypercube spanned by the particle position and the
!!$            ! opposite lattice point
!!$            opp_lp = lx + opp_lp_sur(:,l)
!!$            weight = abs(product(real(opp_lp-(/u*tnx,v*tny,w*tnz/),kind=rk)-P(i)%x))
!!$            Nvp(lp(1),lp(2),lp(3),:) = Nvp(lp(1),lp(2),lp(3),:)+weight*P(i)%v
!!$            Nnpc(lp(1),lp(2),lp(3)) = Nnpc(lp(1),lp(2),lp(3))+weight
!!$          end do lattice_points
!!$        end do pbc_z
!!$      end do pbc_y
!!$    end do pbc_x
!!$
!!$    if (ii<=nlocal) then
!!$      i = list(i)
!!$    else
!!$      i = i+1
!!$    endif
!!$  enddo particles
!!$#endif
!!$
!!$  profiles: do k=1,size(prf)
!!$    prf_direction: do u=1,prf(k)%size(1)
!!$      avg1_direction: do v=1,prf(k)%size(2)
!!$        avg2_direction: do w=1,prf(k)%size(3)
!!$          ! position along profile directions
!!$          rp = prf(k)%s&
!!$            &+real(u-1,kind=rk)*prf(k)%d(1,:)&
!!$            &+real(v-1,kind=rk)*prf(k)%d(2,:)&
!!$            &+real(w-1,kind=rk)*prf(k)%d(3,:)
!!$
!!$          ! normalize according to periodic boundaries
!!$          where (rp<0.5_rk) rp = rp+real((/tnx,tny,tnz/),kind=rk)
!!$          where (rp>=real((/tnx,tny,tnz/),kind=rk)+0.5_rk) rp = rp-real((/tnx,tny,tnz/),kind=rk)
!!$
!!$          ! accumulate data only on that process that owns the
!!$          ! node which is closest to  rp
!!$          if (  any(rp<real(start,kind=rk)-0.5_rk).or.&
!!$              & any(rp>=real(start+(/nx,ny,nz/),kind=rk)-0.5_rk)) cycle
!!$
!!$          ! interpolate mass density, velocity, and site occupation
!!$          
!!$          !FIXME This is just commented out to make it compile for now
!!$          !profiles not supported in pre-release
!!$ 
!!$          !CALL fluid_velocity_and_density_and_site_occupation(N,rp,ivf,imf,irf&
!!$!#ifndef SINGLEFLUID
!!$!            &,iv_r,im_r,ir_r,iv_b,im_b,ir_b&
!!$!#endif
!!$!            &,ifs,irs&
!!$!#ifdef MD
!!$!            &,ips&
!!$!#endif
!!$!            &)
!!$
!!$          ! accumulate interpolated quantities
!!$          prf(k)%l(u)%rhof = prf(k)%l(u)%rhof + irf
!!$          prf(k)%l(u)%vf = prf(k)%l(u)%vf + matmul(prf(k)%vd,ivf)
!!$          prf(k)%l(u)%mf = prf(k)%l(u)%mf + matmul(prf(k)%vd,imf)
!!$#ifndef SINGLEFLUID
!!$          prf(k)%l(u)%rho_r = prf(k)%l(u)%rho_r + ir_r
!!$          prf(k)%l(u)%v_r = prf(k)%l(u)%v_r + matmul(prf(k)%vd,iv_r)
!!$          prf(k)%l(u)%m_r = prf(k)%l(u)%m_r + matmul(prf(k)%vd,im_r)
!!$          prf(k)%l(u)%rho_b = prf(k)%l(u)%rho_b + ir_b
!!$          prf(k)%l(u)%v_b = prf(k)%l(u)%v_b + matmul(prf(k)%vd,iv_b)
!!$          prf(k)%l(u)%m_b = prf(k)%l(u)%m_b + matmul(prf(k)%vd,im_b)
!!$#endif
!!$          prf(k)%l(u)%nfs = prf(k)%l(u)%nfs + ifs
!!$          prf(k)%l(u)%nrs = prf(k)%l(u)%nrs + irs
!!$#ifdef MD
!!$          prf(k)%l(u)%nps = prf(k)%l(u)%nps + ips
!!$
!!$          ! interpolate particle center density and velocity
!!$          ! at the continuous profile position  rp.
!!$          CALL interpolate(Nnpc,rp,ipc)
!!$          CALL interpolate(Nvp,rp,ivp)
!!$
!!$          ! accumulate additional interpolated quantities
!!$          prf(k)%l(u)%npc = prf(k)%l(u)%npc + ipc
!!$          prf(k)%l(u)%vp = prf(k)%l(u)%vp + matmul(prf(k)%vd,ivp)
!!$#endif
!!$        end do avg2_direction
!!$      end do avg1_direction
!!$    end do prf_direction
!!$  end do profiles
!!$
!!$  if (myrankc==0) n_samples = n_samples + 1
!!$
!!$  dump: if (mod(nt,n_sci_profile_dump)==0) then
!!$    do k=1,size(prf)
!!$      CALL MPI_Reduce(prf(k)%l,prfsum(k)%l,size(prf(k)%l),layer_mpitype,sum_layer_mpiop,0,comm_cart,ierror)
!!$
!!$      rank0: if (myrankc==0) then
!!$        CALL lbe_make_filename_output(profile_file_name,'profile-'//trim(prf(k)%name),'.asc',nt)
!!$        open (unit=profile_file_unit,file=profile_file_name,status='REPLACE',action='WRITE',recl=500)
!!$        write (unit=profile_file_unit,fmt='("# direction=(/",2(ES15.8,","),ES15.8,"/)")') prf(k)%d(1,:)
!!$        write (unit=profile_file_unit,fmt='("# position 1 corresponds to (/",2(I5,","),I5,"/)")') prf(k)%s
!!$        write (unit=profile_file_unit,fmt='("# profile length: ",I5)') prf(k)%size(1)
!!$        write (unit=profile_file_unit,fmt='("# averaging lengths: ",I5,"x",I5)') prf(k)%size(2:3)
!!$        write (unit=profile_file_unit,fmt='("# base for velocity components: '//'(/",2(ES15.8,","),ES15.8,"/), '&
!!$            &//'(/",2(ES15.8,","),ES15.8,"/), '//'(/",2(ES15.8,","),ES15.8,"/)")') prf(k)%vd(1,:),prf(k)%vd(2,:),prf(k)%vd(3,:)
!!$        write (unit=profile_file_unit,fmt='("# average of ",I9," samples")') n_samples
!!$        write (unit=profile_file_unit,fmt='("# ")')
!!$#ifndef SINGLEFLUID
!!$#ifdef MD
!!$        write (unit=profile_file_unit,fmt='("# pos'&
!!$             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
!!$             &//'            vr1             vr2             vr3             mr1             mr2             mr3             rho_r'&
!!$             &//'           vb1             vb2             vb3             mb1             mb2             mb3             rho_b'&
!!$             &//'           vp1             vp2             vp3             npc             nps'&
!!$             &//'              nfs             nrs'&
!!$             &//'")')
!!$#else
!!$        write (unit=profile_file_unit,fmt='("# pos'&
!!$             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
!!$             &//'            vr1             vr2             vr3             mr1             mr2             mr3             rho_r'&
!!$             &//'           vb1             vb2             vb3             mb1             mb2             mb3             rho_b'&
!!$             &//'           nfs             nrs'&
!!$             &//'")')
!!$#endif
!!$#else
!!$#ifdef MD
!!$        write (unit=profile_file_unit,fmt='("# pos'&
!!$             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
!!$             &//'            vp1             vp2             vp3             npc             nps'&
!!$             &//'              nfs             nrs'&
!!$             &//'")')
!!$#else
!!$        write (unit=profile_file_unit,fmt='("# pos'&
!!$             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
!!$             &//'            nfs             nrs'&
!!$             &//'")')
!!$#endif
!!$#endif
!!$
!!$        do i=lbound(prfsum(k)%l,1),ubound(prfsum(k)%l,1)
!!$          ! normalize and write cumulated numbers
!!$          total_sites = prfsum(k)%l(i)%nfs + prfsum(k)%l(i)%nrs
!!$#ifdef MD
!!$          total_sites = total_sites + prfsum(k)%l(i)%nps
!!$#endif
!!$
!!$          ! avoid output of NaN if site or particle counts
!!$          ! are zero, instead output zero, this is still
!!$          ! understandable from the output since also the
!!$          ! counts are dumped (as fractions, however).
!!$          if (prfsum(k)%l(i)%nfs==0.0_rk) then
!!$            vf = 0.0_rk
!!$            mf = 0.0_rk
!!$            rhof = 0.0_rk
!!$#ifndef SINGLEFLUID
!!$            v_r = 0.0_rk
!!$            m_r = 0.0_rk
!!$            rho_r = 0.0_rk
!!$            v_b = 0.0_rk
!!$            m_b = 0.0_rk
!!$            rho_b = 0.0_rk
!!$#endif
!!$          else
!!$            vf = prfsum(k)%l(i)%vf/prfsum(k)%l(i)%nfs
!!$            mf = prfsum(k)%l(i)%mf
!!$            rhof = prfsum(k)%l(i)%rhof/prfsum(k)%l(i)%nfs
!!$#ifndef SINGLEFLUID
!!$            v_r = prfsum(k)%l(i)%v_r/prfsum(k)%l(i)%nfs
!!$            m_r = prfsum(k)%l(i)%m_r
!!$            rho_r = prfsum(k)%l(i)%rho_r/prfsum(k)%l(i)%nfs
!!$            v_b = prfsum(k)%l(i)%v_b/prfsum(k)%l(i)%nfs
!!$            m_b = prfsum(k)%l(i)%m_b
!!$            rho_b = prfsum(k)%l(i)%rho_b/prfsum(k)%l(i)%nfs
!!$#endif
!!$          end if
!!$
!!$          if (total_sites==0.0_rk) then
!!$            nfs = 0.0_rk
!!$            nrs = 0.0_rk
!!$#ifdef MD
!!$            npc = 0.0_rk
!!$            nps = 0.0_rk
!!$#endif
!!$          else
!!$            nfs = prfsum(k)%l(i)%nfs/total_sites
!!$            nrs = prfsum(k)%l(i)%nrs/total_sites
!!$#ifdef MD
!!$            npc = prfsum(k)%l(i)%npc/total_sites
!!$            nps = prfsum(k)%l(i)%nps/total_sites
!!$#endif
!!$          end if
!!$#ifdef MD
!!$          if (prfsum(k)%l(i)%npc==0.0_rk) then
!!$            vp = 0.0_rk
!!$          else
!!$            vp = prfsum(k)%l(i)%vp/prfsum(k)%l(i)%npc
!!$          end if
!!$#endif
!!$
!!$#ifndef MD
!!$#ifndef SINGLEFLUID
!!$          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,23(ES15.8,:,X))') i,vf,mf,rhof,v_r,m_r,rho_r,v_b,m_b,rho_b,nfs,nrs
!!$#else
!!$          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,9(ES15.8,:,X))') i,vf,mf,rhof,nfs,nrs
!!$#endif
!!$#else
!!$#ifndef SINGLEFLUID
!!$          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,28(ES15.8,:,X))') i,vf,mf,rhof,v_r,m_r,rho_r,v_b,m_b,rho_b,vp,npc,nps,nfs,nrs
!!$#else
!!$          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,14(ES15.8,:,X))') i,vf,mf,rhof,vp,npc,nps,nfs,nrs
!!$#endif
!!$#endif
!!$        end do
!!$
!!$        close (profile_file_unit)
!!$      end if rank0
!!$    end do
!!$    CALL reset_profile(prf)
!!$
!!$    if (myrankc==0) n_samples = 0
!!$  end if dump
!!$end subroutine dump_profile
!!$
!!$!> reset accumulation buffers to zero
!!$elemental subroutine reset_profile(p)
!!$  ! Be careful:  intent(out)  would mean to assign a new object and thus
!!$  ! loose everything not explicitly defined here, namely the allocation
!!$  ! status of the array  p%l .
!!$  type(profile),intent(inout) :: p
!!$
!!$  p%l(:)%vf(1) = 0.0_rk
!!$  p%l(:)%vf(2) = 0.0_rk
!!$  p%l(:)%vf(3) = 0.0_rk
!!$  p%l(:)%mf(1) = 0.0_rk
!!$  p%l(:)%mf(2) = 0.0_rk
!!$  p%l(:)%mf(3) = 0.0_rk
!!$  p%l(:)%rhof = 0.0_rk
!!$#ifndef SINGLEFLUID
!!$  p%l(:)%v_r(1) = 0.0_rk
!!$  p%l(:)%v_r(2) = 0.0_rk
!!$  p%l(:)%v_r(3) = 0.0_rk
!!$  p%l(:)%m_r(1) = 0.0_rk
!!$  p%l(:)%m_r(2) = 0.0_rk
!!$  p%l(:)%m_r(3) = 0.0_rk
!!$  p%l(:)%rho_r = 0.0_rk
!!$  p%l(:)%v_b(1) = 0.0_rk
!!$  p%l(:)%v_b(2) = 0.0_rk
!!$  p%l(:)%v_b(3) = 0.0_rk
!!$  p%l(:)%m_b(1) = 0.0_rk
!!$  p%l(:)%m_b(2) = 0.0_rk
!!$  p%l(:)%m_b(3) = 0.0_rk
!!$  p%l(:)%rho_b = 0.0_rk
!!$#endif
!!$  p%l(:)%nfs = 0.0_rk
!!$  p%l(:)%nrs = 0.0_rk
!!$#ifdef MD
!!$  p%l(:)%vp(1) = 0.0_rk
!!$  p%l(:)%vp(2) = 0.0_rk
!!$  p%l(:)%vp(3) = 0.0_rk
!!$  p%l(:)%npc = 0.0_rk
!!$  p%l(:)%nps = 0.0_rk
!!$#endif
!!$end subroutine reset_profile
!!$
!!$!> initialization of profiles and associated mpi stuff
!!$subroutine setup_dump_profile
!!$  integer ierror,k,stat
!!$  character,parameter :: axischar(3)=(/'x','y','z'/) ! for file names
!!$  ! unit vectors parallel to coordinate axes
!!$  real(kind=rk),parameter :: e(3,3) = reshape( (/1.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk/),(/3,3/) )
!!$  real(kind=rk) :: norms(3) ! for choice of base containing body force
!!$  integer :: maxnorm(1)    ! index of direction to build 5th vector
!!$
!!$  ! creation of custom mpi data type for  type(layer)
!!$#ifndef SINGLEFLUID
!!$#ifndef MD
!!$  integer,parameter :: lm_n=13
!!$#else
!!$  integer,parameter :: lm_n=16
!!$#endif
!!$#else
!!$#ifndef MD
!!$  integer,parameter :: lm_n=7
!!$#else
!!$  integer,parameter :: lm_n=10
!!$#endif
!!$#endif
!!$  integer lm_lengths(lm_n),lm_types(lm_n),lm_displs(lm_n)
!!$  integer(kind=MPI_ADDRESS_KIND) :: lm_base,lm_addrs(lm_n)
!!$  ! external sum_layer
!!$
!!$  ! calculate additional profiles for each vector of the
!!$  ! orthogonal base containing the body force vector which is
!!$  ! not parallel to one of the coordinate axes
!!$  select case (count((/g_accn_x,g_accn_y,g_accn/)==0.0_rk))
!!$  case (3,2) ! body force switched off or parallel to one axis
!!$    allocate(prf(3),stat=stat)
!!$    CALL check_allocate(stat,'setup_dump_profile(): prf(3)')
!!$  case (1)             ! body force in plane with two axes
!!$    allocate(prf(5),stat=stat)
!!$    CALL check_allocate(stat,'setup_dump_profile(): prf(5)')
!!$  case default         ! general case
!!$    allocate(prf(6),stat=stat)
!!$    CALL check_allocate(stat,'setup_dump_profile(): prf(6)')
!!$  end select
!!$
!!$  ! Setup profile names, directions, starting points, and
!!$  ! dimensions. There are 3, 5, or 6 profiles. Choose the same
!!$  ! starting point  prf%s  for all profiles even if that means
!!$  ! that the associated direction might point out of the real
!!$  ! system, but we assume pbc.
!!$  do k=1,3
!!$      prf(k)%name = axischar(k)
!!$      prf(k)%d = cshift(e,k-1)
!!$      prf(k)%size = cshift((/tnx,tny,tnz/),k-1)
!!$      prf(k)%vd(:,:) = e(:,:)
!!$      prf(k)%s = (/1,1,1/)
!!$  end do
!!$
!!$  if (size(prf)>=4) then
!!$    ! 4th profile
!!$    prf(4)%name = 'g'
!!$    prf(4)%d(1,:) = unit_vector((/g_accn_x,g_accn_y,g_accn/))
!!$    prf(4)%s = (/1,1,1/)
!!$
!!$    ! choose system diagonal as off-axis profile or averaging length
!!$    prf(4)%size(1) = nint(norm(real((/tnx,tny,tnz/),kind=rk)))
!!$
!!$    ! fifth profile: direction is the cross product of the body
!!$    ! force vector with that unit vector which results in the
!!$    ! largest absolute value
!!$    do k=1,3
!!$      norms(k) = norm(cross_product((/g_accn_x,g_accn_y,g_accn/),e(k,:)))
!!$    end do
!!$
!!$    maxnorm = maxloc(norms)
!!$    prf(5)%name = 'gx'//axischar(maxnorm(1))
!!$    prf(5)%d(1,:) = unit_vector(cross_product((/g_accn_x,g_accn_y,g_accn/),e(maxnorm(1),:)))
!!$    prf(5)%s = (/1,1,1/)
!!$
!!$    ! maybe sixth profile
!!$    if (size(prf)==6) then
!!$      prf(6)%name = 'gx'//trim(prf(5)%name)
!!$      prf(6)%d(1,:) = unit_vector(cross_product((/g_accn_x,g_accn_y,g_accn/),prf(5)%d(1,:)))
!!$      prf(6)%s = (/1,1,1/)
!!$
!!$      ! profile and averaging directions
!!$      prf(4)%d(2,:) = prf(5)%d(1,:)
!!$      prf(4)%d(3,:) = prf(6)%d(1,:)
!!$      prf(5)%d = cshift(prf(4)%d,1)
!!$      prf(6)%d = cshift(prf(4)%d,2)
!!$
!!$      ! profile and averaging dimensions
!!$      prf(4)%size(2:3) = prf(4)%size(1)
!!$      prf(5)%size = prf(4)%size
!!$      prf(6)%size = prf(4)%size
!!$
!!$      ! directions for velocity projection
!!$      prf(4)%vd = prf(4)%d
!!$      prf(5)%vd = prf(4)%vd
!!$      prf(6)%vd = prf(4)%vd
!!$    else
!!$      ! make sure  d  and  vd  contain a right-hand-rule-oriented
!!$      ! system also for the case of 5 profiles
!!$      prf(4)%d(2,:) = e(maxnorm(1),:)
!!$      prf(4)%d(3,:) = prf(5)%d(1,:)
!!$      prf(5)%d = cshift(prf(4)%d,2)
!!$
!!$      prf(4)%vd = prf(4)%d
!!$      prf(5)%vd = prf(4)%vd
!!$
!!$      ! profile and averaging dimensions
!!$      prf(4)%size(2) = prf(1)%size(maxnorm(1))
!!$      prf(4)%size(3) = prf(4)%size(1)
!!$      prf(5)%size = cshift(prf(4)%size,2)
!!$    end if
!!$  end if
!!$
!!$  ! allocate accumulation buffers
!!$  do k=1,size(prf)
!!$    allocate(prf(k)%l(prf(k)%size(1)),stat=stat)
!!$    CALL check_allocate(stat,'setup_dump_profile(): prf(k)%l(prf(k)%size(1))')
!!$  end do
!!$
!!$  ! reset buffers
!!$  CALL reset_profile(prf)
!!$
!!$  ! build custom mpi data type for  layer
!!$  lm_lengths(1) = 1    ! start of layer in memory
!!$  lm_types(1) = MPI_LB
!!$  CALL MPI_Get_address(prf(1)%l(1),lm_addrs(1),ierror)
!!$  lm_lengths(2) = 3    ! fluid velocity
!!$  lm_types(2) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%vf(1),lm_addrs(2),ierror)
!!$  lm_lengths(3) = 3    ! mass flow
!!$  lm_types(3) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%mf(1),lm_addrs(3),ierror)
!!$  lm_lengths(4) = 1    ! fluid mass density
!!$  lm_types(4) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%rhof,lm_addrs(4),ierror)
!!$  lm_lengths(5) = 1    ! number of fluid sites
!!$  lm_types(5) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%nfs,lm_addrs(5),ierror)
!!$  lm_lengths(6) = 1    ! number of solid rock sites
!!$  lm_types(6) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%nrs,lm_addrs(6),ierror)
!!$#ifndef SINGLEFLUID
!!$  lm_lengths(7) = 3    ! fluid velocity r
!!$  lm_types(7) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%v_r(1),lm_addrs(7),ierror)
!!$  lm_lengths(8) = 3    ! mass flow r
!!$  lm_types(8) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%m_r(1),lm_addrs(8),ierror)
!!$  lm_lengths(9) = 1    ! fluid mass density r
!!$  lm_types(9) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%rho_r,lm_addrs(9),ierror)
!!$  lm_lengths(10) = 3    ! fluid velocity b
!!$  lm_types(10) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%v_b(1),lm_addrs(10),ierror)
!!$  lm_lengths(11) = 3    ! mass flow b
!!$  lm_types(11) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%m_b(1),lm_addrs(11),ierror)
!!$  lm_lengths(12) = 1    ! fluid mass density b
!!$  lm_types(12) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%rho_b,lm_addrs(12),ierror)
!!$#ifdef MD
!!$  lm_lengths(13) = 3    ! particle velocity
!!$  lm_types(13) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%vp(1),lm_addrs(13),ierror)
!!$  lm_lengths(14) = 1    ! number of particle centers
!!$  lm_types(14) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%npc,lm_addrs(14),ierror)
!!$  lm_lengths(15) = 1    ! number of moving rock sites
!!$  lm_types(15) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%nps,lm_addrs(15),ierror)
!!$#endif
!!$#else
!!$#ifdef MD
!!$  lm_lengths(7) = 3    ! particle velocity
!!$  lm_types(7) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%vp(1),lm_addrs(7),ierror)
!!$  lm_lengths(8) = 1    ! number of particle centers
!!$  lm_types(8) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%npc,lm_addrs(8),ierror)
!!$  lm_lengths(9) = 1    ! number of moving rock sites
!!$  lm_types(9) = MPI_REAL8
!!$  CALL MPI_Get_address(prf(1)%l(1)%nps,lm_addrs(9),ierror)
!!$#endif
!!$#endif
!!$  lm_lengths(lm_n) = 1
!!$  lm_types(lm_n) = MPI_UB
!!$  CALL MPI_Get_address(prf(1)%l(2),lm_addrs(lm_n),ierror)
!!$
!!$  CALL MPI_Get_address(prf(1)%l(1),lm_base,ierror) ! base address
!!$  lm_displs(1:lm_n) = lm_addrs(1:lm_n) - lm_base
!!$
!!$  CALL MPI_Type_struct(lm_n,lm_lengths,lm_displs,lm_types,layer_mpitype,ierror)
!!$  CALL MPI_Type_commit(layer_mpitype,ierror)
!!$
!!$  ! register custom mpi summation operation for  layer
!!$  CALL MPI_Op_create(sum_layer,.true.,sum_layer_mpiop,ierror)
!!$
!!$  ! allocate more profiles just to have the summation buffers
!!$  allocate (prfsum(size(prf)),stat=stat)
!!$  CALL check_allocate(stat,'setup_dump_profile(): prfsum(size(prf))')
!!$  do k=1,size(prf)
!!$    allocate (prfsum(k)%l(lbound(prf(k)%l,1):ubound(prf(k)%l,1)),stat=stat)
!!$    CALL check_allocate(stat,'setup_dump_profile(): '//'prfsum(k)%l(lbound(prf(k)%l,1):ubound(prf(k)%l,1))')
!!$  end do
!!$
!!$#ifdef MD
!!$  ! allocate arrays to store particle data on-site
!!$  allocate (Nvp(0:nx+1,0:ny+1,0:nz+1,1:3),Nnpc(0:nx+1,0:ny+1,0:nz+1),stat=stat)
!!$  CALL check_allocate(stat,'Nvp(0:nx+1,0:ny+1,0:nz+1,1:3),Nnpc(0:nx+1,0:ny+1,0:nz+1)')
!!$#endif
!!$
!!$  if (myrankc==0) then
!!$    n_samples = 0
!!$
!!$    if (n_sci_profile_dump==0) then 
!!$      CALL log_msg("FATAL ERROR: n_sci_profile_dump == 0 but sci_profile == .true. . Aborting...",.false.)
!!$      CALL Abend
!!$    end if
!!$  end if
!!$end subroutine setup_dump_profile
!!$
!!$!> custom mpi reduction operation to sum objects of  type(layer)
!!$subroutine sum_layer(invec,inoutvec,len,type)
!!$
!!$  implicit none
!!$  type(layer),intent(in) :: invec(len)
!!$  type(layer),intent(inout) :: inoutvec(len)
!!$  integer,intent(in) :: len,type
!!$  integer i
!!$
!!$  do i=1,len
!!$    inoutvec(i)%vf(:) = invec(i)%vf(:) + inoutvec(i)%vf(:)
!!$    inoutvec(i)%mf(:) = invec(i)%mf(:) + inoutvec(i)%mf(:)
!!$    inoutvec(i)%rhof = invec(i)%rhof + inoutvec(i)%rhof
!!$#ifndef SINGLEFLUID
!!$    inoutvec(i)%v_r(:) = invec(i)%v_r(:) + inoutvec(i)%v_r(:)
!!$    inoutvec(i)%m_r(:) = invec(i)%m_r(:) + inoutvec(i)%m_r(:)
!!$    inoutvec(i)%rho_r = invec(i)%rho_r + inoutvec(i)%rho_r
!!$    inoutvec(i)%v_b(:) = invec(i)%v_b(:) + inoutvec(i)%v_b(:)
!!$    inoutvec(i)%m_b(:) = invec(i)%m_b(:) + inoutvec(i)%m_b(:)
!!$    inoutvec(i)%rho_b = invec(i)%rho_b + inoutvec(i)%rho_b
!!$#endif
!!$    inoutvec(i)%nfs = invec(i)%nfs + inoutvec(i)%nfs
!!$    inoutvec(i)%nrs = invec(i)%nrs + inoutvec(i)%nrs
!!$#ifdef MD
!!$    inoutvec(i)%vp(:) = invec(i)%vp(:) + inoutvec(i)%vp(:)
!!$    inoutvec(i)%npc = invec(i)%npc + inoutvec(i)%npc
!!$    inoutvec(i)%nps = invec(i)%nps + inoutvec(i)%nps
!!$#endif
!!$  end do
!!$end subroutine sum_layer
!!$!> \}

end module lb3d_io_dump_data_module
