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

!> XDR output routines.
!>
!>Note that these routines use the xdrf library which has to
!>be compiled independently.
!>See http://md.chem.rug.nl/~hoesel/xdrf.html
!>USEXDRF has to be set in the source.
module lb3d_io_xdrf_module
#ifdef USEXDRF

  use lb3d_mpi_module
  use lb3d_log_module
  use lb3d_lattice_module

  use lb3d_config_module, only: amass_b,amass_r,amass_s,beta,fb,fg,fr,g_accn&
       &,g_accn_max,g_accn_max_x,g_accn_max_y,g_accn_min,g_accn_min_x&
       &,g_accn_min_y,g_accn_x,g_accn_y,g_br,g_bs,g_ss,g_wb,g_wr,n_checkpoint&
       &,n_sci_arrows,n_sci_dir,n_sci_flo,n_sci_fluxz,n_sci_int,n_sci_massfluxz&
       &,n_sci_od,n_sci_pressure,n_sci_profile,n_sci_profile_dump,n_sci_rock&
       &,n_sci_sur,n_sci_vel,n_sci_wd,post,s03_b,s03_r,s03_s,s05_b,s05_r,s05_s&
       &,s11_b,s11_r,s11_s,s14_b,s14_r,s14_s,s17_b,s17_r,s17_s,shear_omega&
       &,shear_u,tau_b,tau_d,tau_r,tau_s,tau_wb,tau_wr,taubulk_b,taubulk_d&
       &,taubulk_r,taubulk_s,write_AVS_fld,xdrfsloppy,nt,dump_double
  
  use lb3d_io_helper_module,only:nflags,dump_avs_fld,die_unless_exists,die_unless_nonzero,lbe_make_filename_cp
  implicit none

contains

subroutine check_xdrfopen(ierr,msg)
  integer,intent(in) :: ierr
  character(len=*),intent(in) :: msg

  call die_unless_nonzero(ierr,'xdrfopen returned an error ("'//msg//'")')
end subroutine check_xdrfopen

subroutine check_xdrfint(ierr)
  integer :: ierr
  character(len=64) :: string = 'xdrfint returned an error'
  CALL die_unless_nonzero(ierr,string)
end subroutine check_xdrfint

subroutine check_xdrfdouble(ierr)
  integer :: ierr
  character(len=64) :: string = 'xdrfdouble returned an error'
  CALL die_unless_nonzero(ierr,string)
end subroutine check_xdrfdouble

subroutine check_xdrfsloppy()
  if ( xdrfsloppy ) then
    CALL log_msg("Sloppy xdrf enabled. Going on.",.false.)
    return
  else
    CALL Abend
  endif
end subroutine check_xdrfsloppy

!> dumps all important simulation parameters and is to be called on the
!> root processor only
subroutine checkpoint_params_xdr
  character(len=256)     :: filename
  integer :: ierr
  integer :: ixdrs

  CALL lbe_make_filename_cp(filename,'checkparams','.xdr',nt)
  CALL xdrfopen(ixdrs,filename,"w",ierr)
  if(ierr .eq. 0) then
    CALL log_msg("checkpoint_params_xdr: xdrfopen failed.",.false.)
    CALL check_xdrfsloppy()
  end if

  CALL xdrfdouble(ixdrs,fr, ierr)
  CALL xdrfdouble(ixdrs,fg, ierr)
  CALL xdrfdouble(ixdrs,fb, ierr)
  CALL xdrfdouble(ixdrs,beta, ierr)
  CALL xdrfdouble(ixdrs,amass_b, ierr)
  CALL xdrfdouble(ixdrs,amass_r, ierr)
  CALL xdrfdouble(ixdrs,amass_s, ierr)
  CALL xdrfdouble(ixdrs,tau_b, ierr)
  CALL xdrfdouble(ixdrs,tau_r, ierr)
  CALL xdrfdouble(ixdrs,tau_s, ierr)
  CALL xdrfdouble(ixdrs,tau_d, ierr)
  CALL xdrfdouble(ixdrs,taubulk_b, ierr)
  CALL xdrfdouble(ixdrs,taubulk_r, ierr)
  CALL xdrfdouble(ixdrs,taubulk_s, ierr)
  CALL xdrfdouble(ixdrs,taubulk_d, ierr)
  CALL xdrfdouble(ixdrs,g_br, ierr)
  CALL xdrfdouble(ixdrs,g_bs, ierr)
  CALL xdrfdouble(ixdrs,g_ss, ierr)
  CALL xdrfdouble(ixdrs,g_wr, ierr)
  CALL xdrfdouble(ixdrs,g_wb, ierr)
  CALL xdrfdouble(ixdrs,tau_wr, ierr)
  CALL xdrfdouble(ixdrs,tau_wb, ierr)
  CALL xdrfdouble(ixdrs,shear_u, ierr)
  CALL xdrfdouble(ixdrs,shear_omega, ierr)
  CALL xdrfdouble(ixdrs,g_accn, ierr)
  CALL xdrfdouble(ixdrs,g_accn_x, ierr)
  CALL xdrfdouble(ixdrs,g_accn_y, ierr)
  CALL xdrfdouble(ixdrs,s03_r, ierr)
  CALL xdrfdouble(ixdrs,s05_r, ierr)
  CALL xdrfdouble(ixdrs,s11_r, ierr)
  CALL xdrfdouble(ixdrs,s14_r, ierr)
  CALL xdrfdouble(ixdrs,s17_r, ierr)
  CALL xdrfdouble(ixdrs,s03_b, ierr)
  CALL xdrfdouble(ixdrs,s05_b, ierr)
  CALL xdrfdouble(ixdrs,s11_b, ierr)
  CALL xdrfdouble(ixdrs,s14_b, ierr)
  CALL xdrfdouble(ixdrs,s17_b, ierr)
  CALL xdrfdouble(ixdrs,s03_s, ierr)
  CALL xdrfdouble(ixdrs,s05_s, ierr)
  CALL xdrfdouble(ixdrs,s11_s, ierr)
  CALL xdrfdouble(ixdrs,s14_s, ierr)
  CALL xdrfdouble(ixdrs,s17_s, ierr)
  CALL xdrfint(ixdrs,n_sci_int,ierr)
  CALL xdrfint(ixdrs,n_sci_sur,ierr)
  CALL xdrfint(ixdrs,n_sci_od,ierr)
  CALL xdrfint(ixdrs,n_sci_wd,ierr)
  CALL xdrfint(ixdrs,n_sci_dir,ierr)
  CALL xdrfint(ixdrs,n_sci_vel,ierr)
  CALL xdrfint(ixdrs,n_sci_flo,ierr)
  CALL xdrfint(ixdrs,n_sci_arrows,ierr)
  CALL xdrfint(ixdrs,n_sci_rock,ierr)
  CALL xdrfint(ixdrs,n_sci_pressure,ierr)
  CALL xdrfint(ixdrs,n_sci_fluxz,ierr)
  CALL xdrfint(ixdrs,n_sci_massfluxz,ierr)
  CALL xdrfint(ixdrs,n_sci_profile,ierr)
  CALL xdrfint(ixdrs,n_sci_profile_dump,ierr)
  CALL xdrfint(ixdrs,n_checkpoint,ierr)
  CALL xdrfint(ixdrs,g_accn_min, ierr)
  CALL xdrfint(ixdrs,g_accn_max, ierr)
  CALL xdrfint(ixdrs,g_accn_min_x, ierr)
  CALL xdrfint(ixdrs,g_accn_max_x, ierr)
  CALL xdrfint(ixdrs,g_accn_min_y, ierr)
  CALL xdrfint(ixdrs,g_accn_max_y, ierr)

  CALL xdrfclose(ixdrs,ierr)
end subroutine checkpoint_params_xdr

!> Restores a file written by checkpoint_params_xdr and distributes
!> values. To be called on all CPUs
subroutine restore_params_xdr
  character(len=256)     :: filename
  integer                :: ierr
  integer                :: ixdrs

  rank0: if (myrankc == 0) then
    CALL lbe_make_filename_cp(filename,'checkparams','.xdr',nt)
    CALL xdrfopen(ixdrs,filename,"r",ierr)
    if(ierr .eq. 0)then
      CALL log_msg("checkpoint_params_xdr: xdrfopen failed.",.false.)
      CALL check_xdrfsloppy()
    endif

    CALL xdrfdouble(ixdrs,fr, ierr)
    CALL xdrfdouble(ixdrs,fg, ierr)
    CALL xdrfdouble(ixdrs,fb, ierr)
    CALL xdrfdouble(ixdrs,beta, ierr)
    CALL xdrfdouble(ixdrs,amass_b, ierr)
    CALL xdrfdouble(ixdrs,amass_r, ierr)
    CALL xdrfdouble(ixdrs,amass_s, ierr)
    CALL xdrfdouble(ixdrs,tau_b, ierr)
    CALL xdrfdouble(ixdrs,tau_r, ierr)
    CALL xdrfdouble(ixdrs,tau_s, ierr)
    CALL xdrfdouble(ixdrs,tau_d, ierr)
    CALL xdrfdouble(ixdrs,taubulk_b, ierr)
    CALL xdrfdouble(ixdrs,taubulk_r, ierr)
    CALL xdrfdouble(ixdrs,taubulk_s, ierr)
    CALL xdrfdouble(ixdrs,taubulk_d, ierr)
    CALL xdrfdouble(ixdrs,g_br, ierr)
    CALL xdrfdouble(ixdrs,g_bs, ierr)
    CALL xdrfdouble(ixdrs,g_ss, ierr)
    CALL xdrfdouble(ixdrs,g_wr, ierr)
    CALL xdrfdouble(ixdrs,g_wb, ierr)
    CALL xdrfdouble(ixdrs,tau_wr, ierr)
    CALL xdrfdouble(ixdrs,tau_wb, ierr)
    CALL xdrfdouble(ixdrs,shear_u, ierr)
    CALL xdrfdouble(ixdrs,shear_omega, ierr)
    CALL xdrfdouble(ixdrs,g_accn, ierr)
    CALL xdrfdouble(ixdrs,g_accn_x, ierr)
    CALL xdrfdouble(ixdrs,g_accn_y, ierr)

    CALL xdrfdouble(ixdrs,s03_r, ierr)
    CALL xdrfdouble(ixdrs,s05_r, ierr)
    CALL xdrfdouble(ixdrs,s11_r, ierr)
    CALL xdrfdouble(ixdrs,s14_r, ierr)
    CALL xdrfdouble(ixdrs,s17_r, ierr)
    CALL xdrfdouble(ixdrs,s03_b, ierr)
    CALL xdrfdouble(ixdrs,s05_b, ierr)
    CALL xdrfdouble(ixdrs,s11_b, ierr)
    CALL xdrfdouble(ixdrs,s14_b, ierr)
    CALL xdrfdouble(ixdrs,s17_b, ierr)
    CALL xdrfdouble(ixdrs,s03_s, ierr)
    CALL xdrfdouble(ixdrs,s05_s, ierr)
    CALL xdrfdouble(ixdrs,s11_s, ierr)
    CALL xdrfdouble(ixdrs,s14_s, ierr)
    CALL xdrfdouble(ixdrs,s17_s, ierr)

    CALL xdrfint(ixdrs,n_sci_int,ierr)
    CALL xdrfint(ixdrs,n_sci_sur,ierr)
    CALL xdrfint(ixdrs,n_sci_od,ierr)
    CALL xdrfint(ixdrs,n_sci_wd,ierr)
    CALL xdrfint(ixdrs,n_sci_dir,ierr)
    CALL xdrfint(ixdrs,n_sci_vel,ierr)
    CALL xdrfint(ixdrs,n_sci_flo,ierr)
    CALL xdrfint(ixdrs,n_sci_arrows,ierr)
    CALL xdrfint(ixdrs,n_sci_rock,ierr)
    CALL xdrfint(ixdrs,n_sci_pressure,ierr)
    CALL xdrfint(ixdrs,n_sci_fluxz,ierr)
    CALL xdrfint(ixdrs,n_sci_massfluxz,ierr)
    CALL xdrfint(ixdrs,n_sci_profile,ierr)
    CALL xdrfint(ixdrs,n_sci_profile_dump,ierr)
    CALL xdrfint(ixdrs,n_checkpoint,ierr)
    CALL xdrfint(ixdrs,g_accn_min, ierr)
    CALL xdrfint(ixdrs,g_accn_max, ierr)
    CALL xdrfint(ixdrs,g_accn_min_x, ierr)
    CALL xdrfint(ixdrs,g_accn_max_x, ierr)
    CALL xdrfint(ixdrs,g_accn_min_y, ierr)
    CALL xdrfint(ixdrs,g_accn_max_y, ierr)

    CALL xdrfclose(ixdrs,ierr)
  end if rank0

  CALL MPI_Bcast(fr,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(fg,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(fb,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(beta,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(amass_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(amass_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(amass_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_d,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(taubulk_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(taubulk_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(taubulk_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(taubulk_d,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_br,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_bs,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_ss,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_wr,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_wb,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_wr,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(tau_wb,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(shear_u,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(shear_omega,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_x,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_y,1,MPI_REAL8,0,comm_cart,ierr)

  CALL MPI_Bcast(s03_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s05_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s11_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s14_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s17_r,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s03_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s05_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s11_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s14_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s17_b,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s03_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s05_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s11_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s14_s,1,MPI_REAL8,0,comm_cart,ierr)
  CALL MPI_Bcast(s17_s,1,MPI_REAL8,0,comm_cart,ierr)

  CALL MPI_Bcast(n_sci_int,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_sur,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_od,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_wd,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_dir,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_vel,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_flo,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_arrows,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_rock,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_pressure,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_fluxz,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_massfluxz,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_profile,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_sci_profile_dump,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(n_checkpoint,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_min,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_max,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_min_x,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_max_x,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_min_y,1,MPI_INTEGER,0,comm_cart,ierr)
  CALL MPI_Bcast(g_accn_max_y,1,MPI_INTEGER,0,comm_cart,ierr)
end subroutine restore_params_xdr

!> Takes an XDR filehandle denoted by ixdrs, and attempts
!> to seek forwards by skipval bytes.
!> Calls \c Abend on error.
subroutine skip_xdr(ixdrs,skipval)
  integer*8              :: skipval,pos,newpos
  integer                :: ixdrs,ierr
  character(len=128)     :: msgstr
  CALL xdrfgetpos(ixdrs,pos)
  !if (ierr .eq. 0) then
  !        print*,'skip_xdr: xdrfgetpos returned error'
  !        call Abend
  !end if

  ! Seek to low-(xyz) corner of subdomain
  newpos = pos + skipval

  CALL xdrfsetpos(ixdrs,newpos,ierr)

  !if (ierr .eq. 0) then
  !        print*,'skip_xdr: xdrfsetpos returned error'
  !        call Abend
  !end if
  CALL xdrfgetpos(ixdrs,pos)

  if (pos /= newpos) then
    write(msgstr,"('FATAL ERROR: Position could not be set. pos. is: ',I0,' but should be: ', I0,'. Aborting...')") pos, newpos
    CALL log_msg(trim(msgstr),.true.)
    CALL Abend
  end if
end subroutine skip_xdr

subroutine skip_xdr_floats(ixdrs,skipval)
  integer   :: ixdrs
  integer*8 :: skipval
  CALL skip_xdr(ixdrs,4*skipval)
end subroutine skip_xdr_floats

subroutine skip_xdr_double(ixdrs,skipval)
  integer   :: ixdrs
  integer*8 :: skipval
  CALL skip_xdr(ixdrs,8*skipval)
end subroutine skip_xdr_double

!> Reads an XDR-float file into the rock matrix.
!> All CPUs open the file, but each CPU only loads its own
!> subsection of the rock.
subroutine read_rock_xdrf_par(filename,threshold,N)
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer          :: gx,gy,gz ! Global coords of low-(xyz) corner of our subdomain
  integer          :: ixdrs,ierr
  integer          :: x,y,z
  real*4           :: rockfloat !NB real*4 so equiv of C float!
  real(kind=rk)           :: threshold ! Values above threshold are considered rock.
  character(len=*) :: filename
  integer          :: rcount = 0
  integer*8        :: newskip
  integer*8        :: rocksize
  character(len=128)     :: msgstr

  ! NB: global coordinates are zero-based!
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz

  CALL xdrfopen(ixdrs,filename,"r",ierr)
  call xdrfdarrayset(nx*cdims(1),ny*cdims(2),nz*cdims(3),nx,ny,nz,1)

  rocksize = int(tnx,kind=8)*int(tny,kind=8)*int(tnz,kind=8)
  if (rocksize > HUGE(gx)/4-1) then
    CALL log_msg("ATTENTION: UNIX file operations used for seek calls. Check if long int is 64 bit.",.false.)
  end if

  if (ierr .eq. 0) then
    CALL log_msg("FATAL ERROR: read_rock_xdruchar_par: xdrfopen returned error. Aborting...",.true.)
    CALL Abend
  end if

  ! Seek to low-(xyz) corner of subdomain
  newskip =  gx + int(tnx,kind=8)*int(gy,kind=8) + &
                  int(tnx,kind=8)*int(tny,kind=8)* int(gz,kind=8)
  CALL skip_xdr_floats(ixdrs,newskip)

#ifdef XDRROCKWET
  CALL log_msg("  Using XDRROCKWET",.false.)
#else
  write(msgstr,"('  Not using XDRROCKWET, threshold = ',F16.10)") threshold
  CALL log_msg(trim(msgstr),.false.)
#endif

  do z = 1, nz
    do y = 1, ny
      ! Read in an X-span
      do x = 1, nx
        CALL xdrffloat(ixdrs,rockfloat,ierr)
        if (ierr .eq. 0) then
          CALL log_msg("FATAL ERROR: read_rock_xdruchar_par: xdrffloat returned error. Aborting...",.true.)
          !print*,'local coords: ',x,y,z
          CALL Abend
        end if
#ifdef XDRROCKWET
        ! Set rock state to 1 if value in xdr file 
        ! is larger than 0.
        ! Colour the rock according to value in xdr file.
        ! To be able to distinguish between rock/non-rock
        ! subtract an offset of 5.
        !
        if (rockfloat > 0.0) then
          rockfloat = rockfloat - 5.
          N(x,y,z)%rock_state = 1.d0
          rcount = rcount+1

          N(x,y,z)%n_r(:nnonrest)=0.d0
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:nnonrest)=0.d0
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:nnonrest)=0.d0
          N(x,y,z)%n_s(restvec)=0.d0
#endif
          
          if (rockfloat > 0.d0) then
#ifndef SINGLEFLUID
             N(x,y,z)%rock_colour_b = 0.d0
#endif
             N(x,y,z)%rock_colour_r = dble(rockfloat)
#ifdef WALLCONST
             N(x,y,z)%rock_colour_r = dble(rockfloat)/tau_r
#endif
          elseif (rockfloat < 0) then
#ifndef SINGLEFLUID
             N(x,y,z)%rock_colour_b = -dble(rockfloat)
#ifdef WALLCONST
             N(x,y,z)%rock_colour_b = -dble(rockfloat)/tau_b
#endif
#endif
             N(x,y,z)%rock_colour_r = 0.d0
          else
             N(x,y,z)%rock_colour_r = 0.d0
#ifndef SINGLEFLUID
             N(x,y,z)%rock_colour_b = 0.d0
#endif
          endif
       else
          N(x,y,z)%rock_state = 0.d0
        endif
#else
        ! Set rock state to 1 if value in xdr file
        ! is higher than threshold
        if (rockfloat > threshold) then
          N(x,y,z)%rock_state = 1.d0
          rcount = rcount+1
        else
          N(x,y,z)%rock_state = 0.d0
        endif
#endif
      end do ! x
      ! Skip tnx-nx sites to next span in XY slab 
      newskip = tnx-nx
      CALL skip_xdr_floats(ixdrs,newskip)
    end do ! y
    ! Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume 
    newskip = (tnx*tny)-(tnx*ny)
    CALL skip_xdr_floats(ixdrs, newskip)
  end do ! z

  CALL xdrfclose(ixdrs,ierr)

end subroutine read_rock_xdrf_par


!>Reads an XDR-integer file, parsed bitwise into the rock matrix.
!>All CPUs open the file, but each CPU only loads its own
!>subsection of the rock.
subroutine read_bit_rock_xdrf_par(filename,N)
  
  implicit none

  character(len=*) :: filename
  type(lbe_site),dimension(0:,0:,0:) :: N

  integer,parameter            :: local_ik=4
  integer                      :: ixdrs,ierr
  integer(kind=local_ik)             :: ival
  integer                      :: bitsize
 
  integer,dimension(3)         :: rockdims
  logical                      :: bigendian
  logical,dimension(0:local_ik*8-1)  :: bool
  integer                      :: i
  integer                      :: fields,bits,rockbits

  integer                      :: x,y,z
  integer                      :: gx,gy,gz

  integer                      :: skip

  bigendian=.false.
  bitsize=local_ik * 8

  ! open given file for reading
  call xdrfopen(ixdrs,filename,"r",ierr)

  ! read header (assumes 6 integers with coordinate information in every 2nd)
  ! Crude endianess check, assuming lattice length larger 16M are never used.
  ! If the dimension of one coordinate exceeds 65535, this will become gamble.
  call read_bit_rock_header(ixdrs,rockdims,bigendian)

  ! check matching of dimensions
  if (rockdims(1).ne.tnx.or.rockdims(2).ne.tny.or.rockdims(3).ne.tnz) then
     CALL log_msg("FATAL ERROR: read_bit_rock_xdr_par: bit rock dimension do not match system dimension. Aborting...",.true.)
     CALL Abend
  end if

  ! getting global coordinates of local point of origin
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz

  ! initialise local coordinates
  x=1;y=1;z=1

  rockbits=0
  bits=0
  fields=0

  ! read first field
  call xdrfint(ixdrs,ival,ierr)

  ! skip forward to local coordinate origin
  skip=gx+(tnx*gy)+tnx*tny*gz

  do while (ierr.eq.1)

     ! Read data to array of logicals
     if (bigendian.eqv..true.) then
        call read_bit_rock_int2boolArrBE(ival,bool)
     else
        call read_bit_rock_int2boolArr(ival,bool)
     end if

     do i=0,31
        if (skip.eq.0) then
           ! Write data
           if (bool(i).eqv..true.) then
              N(x,y,z)%rock_state=1.d0
              rockbits=rockbits+1
           else
              N(x,y,z)%rock_state=0.d0
           end if
           ! Iterate local coordinates
           call read_bit_rock_iterate(x,y,z,nx,ny,nz,tnx,tny,tnz,skip)
        else
           ! decrease skip count   
           skip=skip-1
        end if
        if (z.gt.nz) then
           ! end of data
           goto 10
        end if
        
        bits=bits + 1

     end do ! 32 bits
     ! update field (integer) count
     fields=fields + 1
     
     ! read next field
     call xdrfint(ixdrs,ival,ierr)

  end do ! while ierr.eq.1, reading integers 

  10 continue
 
  ! Close file handle 
  call xdrfclose(ixdrs,ierr)

end subroutine read_bit_rock_xdrf_par

! Since in read_bit_rock_
subroutine read_bit_rock_iterate(x,y,z,nx,ny,nz,tnx,tny,tnz,skip)

  implicit none

  integer               :: x,y,z
  integer               :: nx,ny,nz
  integer               :: tnx,tny,tnz
  integer               :: skip

  ! nested loops the other way round
  x=x+1
  if (x.gt.nx) then
     x=1
     y=y+1
     ! skip to next local block   
     skip=skip+tnx-nx     
  end if
  if (y.gt.ny) then
     y=1     
     z=z+1
     ! skip to next local block   
     skip = skip+(tnx*tny)-(tnx*ny)
  end if

  return

end subroutine read_bit_rock_iterate

subroutine read_bit_rock_header(ixdrs,rockdims,bigendian)
  
  implicit none
  
  integer                               :: ixdrs,ierr
  integer                               :: ival
  integer,dimension(3)                  :: rockdims
  logical,dimension(0:31)               :: bool
  integer                               :: i
  logical                               :: bigendian

  ! read header data
  do i=1,3
     call xdrfint(ixdrs,ival,ierr)
     call xdrfint(ixdrs,ival,ierr)
     rockdims(i)=ival
     ! check for bigendian, this is a kludge
     if (rockdims(i).gt.16777215) then
        bigendian=.true.
     end if
  end do
  
  if (bigendian.eqv..true.) then
      CALL log_msg("NOTICE: Found axis length larger 16M, assuming big endianess.",.false.)
     do i=1,3
        call read_bit_rock_int2boolArrBE(rockdims(i),bool)
        call read_bit_rock_boolArr2int(bool,rockdims(i))
     end do
  end if

  return

end subroutine read_bit_rock_header

subroutine read_bit_rock_int2boolArrBE(ival,bool)

  implicit none

  integer                               :: ival
  logical,dimension(0:sizeof(ival)*8-1) :: bool
  integer                               :: bytesize
  integer                               :: bitsize
  integer                               :: ierr,i

  bytesize=sizeof(ival)
  bitsize=bytesize * 8
  do i=0,bitsize-1
     bool(i)=btest(ival,8*(bytesize-1-int(i/8))+modulo(i,8))
  end do

  return

end subroutine read_bit_rock_int2boolArrBE

subroutine read_bit_rock_int2boolArr(ival,bool)

  implicit none

  integer                               :: ival 
  logical,dimension(0:sizeof(ival)*8-1) :: bool 
  integer                               :: bytesize 
  integer                               :: bitsize 
  integer                               :: i

  bytesize=sizeof(ival)
  bitsize=bytesize * 8
  do i=1,bitsize-1
     bool(i)=btest(ival,i)
  end do

  return

end subroutine read_bit_rock_int2boolArr

subroutine read_bit_rock_boolArr2int(bool,ival)

  implicit none

  integer                               :: ival
  logical,dimension(0:sizeof(ival)*8-1) :: bool
  integer                               :: bytesize
  integer                               :: bitsize
  integer                               :: i

  ival=0
  bytesize=sizeof(ival)
  bitsize=bytesize * 8
  do i=1,bitsize-1
     if (bool(i).eqv..true.) then
        ival=ival + 2**i
     end if
  end do

  return

end subroutine read_bit_rock_boolArr2int

!> Dumps a scalar in xdr format using the xdrf library.
subroutine dump_scalar_xdr(scalar,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:) :: scalar
  real*4 :: scalartmp
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: ixdrs,ierr,x,y,z
  integer :: nxi,nyi,nzi

  nxi = size(scalar,1)
  nyi = size(scalar,2)
  nzi = size(scalar,3)

  CALL lbe_make_filename_output(filename,trim(name),'.xdr',nt)
  CALL xdrfopen(ixdrs,filename,"w",ierr)

  if( ierr .eq. 0 ) then
    CALL log_msg("dump_scalar_xdr: xdrfopen returned error",.false.)
    CALL check_xdrfsloppy()
  endif

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (dump_double) then
          CALL xdrfdouble(ixdrs, scalar(x,y,z), ierr)
        else
          scalartmp = scalar(x,y,z)
          CALL xdrffloat(ixdrs, scalartmp, ierr)
        end if
      end do
    end do
  end do

  CALL xdrfclose(ixdrs,ierr)
  if (write_AVS_fld) CALL dump_avs_fld(trim(name),nxi,nyi,nzi,1)

end subroutine dump_scalar_xdr

!> Dumps an integer scalar field in xdr format using the xdrf library.
subroutine dump_iscalar_xdr(iscalar,name)
    integer,dimension(1:,1:,1:) :: iscalar
    character(LEN=*) :: name
    character(LEN=256) :: filename
    integer :: ixdrs,ierr,x,y,z
    integer :: nxi,nyi,nzi

    nxi = size(iscalar,1)
    nyi = size(iscalar,2)
    nzi = size(iscalar,3)

    call lbe_make_filename_output(filename,trim(name),'.xdr',nt)
    call xdrfopen(ixdrs,filename,"w",ierr)

    if (ierr .eq. 0) then
       call log_msg("dump_iscalar_xdr: xdrfopen returned error",.false.)
       call check_xdrfsloppy()
    end if

    do z=1,nzi
       do y=1,nyi
          do x=1,nxi
             call xdrfint(ixdrs,iscalar(x,y,z),ierr)
          end do
       end do
    end do
    call xdrfclose(ixdrs,ierr)

    if (write_AVS_fld) &
         &call error('dump_avs_fld() does not support integer fields yet')
end subroutine dump_iscalar_xdr

subroutine dump_scalar_xdr_parallel(scalar,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:) :: scalar
  real*4 :: scalartmp
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: ixdrs,ierr,x,y,z
  integer :: nxi,nyi,nzi
  logical :: temp_bool

  nxi = size(scalar,1)
  nyi = size(scalar,2)
  nzi = size(scalar,3)

  ! set post to true to coerse lbe_make_filename into producing a filename without pe
  temp_bool = post
  post = .true.
  CALL lbe_make_filename_output(filename,trim(name),'.xdr',nt)
  post = temp_bool

  CALL xdrfopen(ixdrs,filename,"w",ierr)
  CALL xdrfdarrayset( nx * cdims(1), ny * cdims(2), nz * cdims(3), nx, ny, nz, 1 )

  if( ierr .eq. 0 ) then
    CALL log_msg("dump_scalar_xdr_parallel: xdrfopen returned error",.false.)
    CALL check_xdrfsloppy()
  endif

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if( dump_double ) then
          CALL xdrfdarraydouble( ixdrs, scalar(x,y,z), ierr, x-1, y-1, z-1, ccoords(1), ccoords(2), ccoords(3), 0 )
        else
          scalartmp = scalar(x,y,z)
          CALL xdrfdarrayfloat( ixdrs, scalartmp, ierr, x-1, y-1, z-1, ccoords(1), ccoords(2), ccoords(3), 0 )
        end if
      end do
    end do
  end do
  CALL xdrfclose(ixdrs,ierr)
  ! if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,1)

end subroutine dump_scalar_xdr_parallel

!>Dumps a scalar in vtk format using the xdrf library.
subroutine dump_scalar_vtk(scalar,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:) :: scalar
  real*4, dimension(:,:,:), allocatable :: scalartmp
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: ixdrs,ierr,x,y,z
  integer :: nxi,nyi,nzi,itmp
  nxi = size(scalar,1)
  nyi = size(scalar,2)
  nzi = size(scalar,3)

  !! VTK output hard-wired to 'sample' for use in UNICORE
  !! framework, ARP, 30/8/2002.
  ! filename = " "
  ! filename = "sample"
  ! filename(7:7) = CHAR(0)
  !! End of ARP's changes
  CALL lbe_make_filename_output(filename,trim(name),'.vtk',nt)
  itmp = index(filename,' ')-1
  filename = filename(1:itmp)//CHAR(0)
  CALL vtkopen(ixdrs,filename,1,nxi,nyi,nzi,ierr)

  if(ierr .ne. 0)then
    CALL log_msg("dump_scalar_vtk: vtkopen failed",.false.)
    return
  endif

  allocate(scalartmp(1:nxi,1:nyi,1:nzi))
  scalartmp = scalar
  CALL vtkfloat(ixdrs, nxi*nyi*nzi, scalartmp, ierr)
  deallocate(scalartmp)

  CALL vtkclose(ixdrs,ierr)
end subroutine dump_scalar_vtk

!> Dump an array of vectors in xdr format.
subroutine dump_vector_xdr(vector,name)
  implicit none
  real(kind=rk),intent(inout),dimension(1:,1:,1:,1:) :: vector
  real*4 :: scalartmp
  character(LEN=*) :: name
  character(LEN=256) :: filename
  integer :: x,y,z,q
  integer :: ixdrs,ierr
  integer :: nxi,nyi,nzi
  integer :: l

  nxi = size(vector,2)
  nyi = size(vector,3)
  nzi = size(vector,4)

  CALL lbe_make_filename_output(filename,trim(name),'.xdr',nt)

  CALL xdrfopen(ixdrs,filename,"w",ierr)

  if( ierr .eq. 0 ) then
    CALL log_msg("dump_vector_xdr: xdrfopen returned error",.false.)
    CALL check_xdrfsloppy()
  endif

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        do q=1,size(vector,1)
          if (dump_double) then
            CALL xdrfdouble(ixdrs, vector(q,x,y,z), ierr)
          else
            scalartmp = vector(q,x,y,z)
            CALL xdrffloat(ixdrs, scalartmp, ierr)
          end if
        end do
      end do
    end do
  end do
  CALL xdrfclose(ixdrs,ierr)
  if (write_AVS_fld) CALL dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
end subroutine dump_vector_xdr

!> Dump an array of vectors in xdr format.
subroutine dump_vector_xdr_parallel(vector,name)
  implicit none
  real(kind=rk), dimension(1:,1:,1:,1:) :: vector
  real*4 :: scalartmp
  character(len=*) :: name
  character(len=256) :: filename
  integer :: x,y,z,q
  integer :: ixdrs,ierr 
  integer :: nxi,nyi,nzi
  logical :: temp_bool;

  nxi = size(vector,2)
  nyi = size(vector,3)
  nzi = size(vector,4)

  ! set post to true to coarse lbe_make_filename into producing a filename without pe
  temp_bool = post
  post = .true.
  call lbe_make_filename_output(filename,trim(name),'.xdr',nt)
  post = temp_bool

  CALL xdrfopen(ixdrs,filename,"w",ierr)
  CALL xdrfdarrayset( nx * cdims(1),ny * cdims(2),nz * cdims(3),nx,ny,nz,size(vector,1) )

  if( ierr .eq. 0 ) then
    CALL log_msg("dump_3scalar_xdr_parallel: xdrfopen returned error",.false.)
    CALL check_xdrfsloppy()
  endif

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        do q=1,size(vector,1)
          if (dump_double) then
            CALL xdrfdouble(ixdrs, vector(q,x,y,z), ierr)
            CALL xdrfdarraydouble( ixdrs, vector(q,x,y,z), ierr, x-1, y-1, z-1, ccoords(1), ccoords(2), ccoords(3), q-1 )
          else
            scalartmp = vector(q,x,y,z)
            CALL xdrfdarrayfloat( ixdrs, scalartmp, ierr, x-1, y-1, z-1, ccoords(1), ccoords(2), ccoords(3), q-1 )
          end if
        end do
      end do
    end do
  end do

  CALL xdrfclose(ixdrs,ierr)
  ! if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
end subroutine dump_vector_xdr_parallel

!>Builds a new lattice state by cloning copies of a smaller
!>saved system.
!>
!>NB: internally, this routine works with zero-based global
!>coordinates (spanning the whole lattice across all CPUs),
!>unless otherwise specified.
!>
!>Builds a new lattice state by copying a smaller
!>saved system somewhere into the new system.
!>
!>NB: internally, this routine works with zero-based global
!>coordinates (spanning the whole lattice across all CPUs),
!>unless otherwise specified.
subroutine restore_upscale(rnt,chkid,N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  logical :: intersects_p
  character(len=128) :: upscalefile = 'upscale.dat'
  integer :: rnt,chkid

  character(len=128) :: seedfile ! Name of upscale restart file
  integer :: seednx,seedny,seednz ! Size of upscale restart file
  logical :: pzero_p  ! Set if we need to skip the first 12 bytes
  namelist /upscale_input/ &
    seedfile, &
    pzero_p, &
    seednx,seedny,seednz
  integer :: gx,gy,gz ! Global coords of the start of our subdomain
  integer,dimension(3) :: seeddims,mysubdims
  integer :: junkint ! Dummy to read integers into

  integer :: repnx,repny,repnz ! No of repeats of seed in each dir
  integer :: repx,repy,repz    ! Repeat coordinates
  integer :: seedx,seedy,seedz ! Point within the seed file
  integer :: x,y,z ! Global coordinates of upscaled point

  integer :: npoints ! debug variable
  integer :: s ! Counter for lattice vector

  integer :: file_id,ierr ! for XDR file I/O
  type(lbe_site) :: seedsite ! Holds a single lattice site read from seed

  ! Read the upscale.dat file. This contains
  ! the total dimensions of the old lattice, and the 
  ! dimensions of the old CPU grid.

  call die_unless_exists(upscalefile)
  open(unit=17,file=upscalefile,status='UNKNOWN')
  read(unit=17,nml=upscale_input)
  close(unit=17)

  print*,'Read upscaler input:'
  print*,' seedfile=<',trim(seedfile),'> pzero_p=',pzero_p
  print*,' Seed size is (',seednx,',',seedny,',',seednz,')'

  seeddims(1) = seednx
  seeddims(2) = seedny
  seeddims(3) = seednz

  ! Check that the seed file exists.
  call die_unless_exists(seedfile)

  ! It does; open it.
  call xdrfopen(file_id,seedfile,"r",ierr)
  call check_xdrfopen(ierr,seedfile)

  ! If pzero_p is set, then the seed file was generated by saving
  ! a checkpoint from CPU 0; skip the three header integers in this case.

  if (pzero_p) then
    call xdrfint(file_id,junkint,ierr)
    call xdrfint(file_id,junkint,ierr)
    call xdrfint(file_id,junkint,ierr)
  end if

  !Find start of our subdomain in global coordinates.
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz
  mysubdims(1)=nx
  mysubdims(2)=ny
  mysubdims(3)=nz

  ! Recall that tnx,tny,tnz are the global lattice dimensions.
  ! Find out how many iterations of the seed span the global lattice.

  repnx = tnx/seednx
  repny = tny/seedny
  repnz = tnz/seednz

  print*,'Repeat dimensions ',repnx,',',repny,',',repnz

  ! Now, run through the seed file.
  ! For each point in the seed file, find the coordinates of
  ! each of the many points in global coordinates onto which
  ! it will map after upscaling.
  !
  ! For each of these points, calculate if it intersects our
  ! subvolume. If it does, set the lattice accordingly.
  !

  npoints = 0
  do seedx=0,seednx-1
    do seedy=0,seedny-1
      do seedz=0,seednz-1

    ! Read a new lattice site in from the seed.
    do s=1,size(seedsite%n_r)
        call xdrfdouble(file_id,seedsite%n_r(s),ierr)
        call check_xdrfdouble(ierr)
    end do
#ifndef SINGLEFLUID
    do s=1,size(seedsite%n_b)
        call xdrfdouble(file_id,seedsite%n_b(s),ierr)
        call check_xdrfdouble(ierr)
    end do
#endif
#ifndef NOSURFACTANT
    do s=1,size(seedsite%n_s)
        call xdrfdouble(file_id,seedsite%n_s(s),ierr)
        call check_xdrfdouble(ierr)
    end do
    do s=1,size(seedsite%da)
        call xdrfdouble(file_id,seedsite%da(s),ierr)
        call check_xdrfdouble(ierr)
    end do
    do s=1,size(seedsite%db)
        call xdrfdouble(file_id,seedsite%db(s),ierr)
        call check_xdrfdouble(ierr)
    end do
#endif
    call xdrfdouble(file_id,seedsite%rock_state,ierr)
    call check_xdrfdouble(ierr)

    do repx=0,repnx
      do repy=0,repny
        do repz=0,repnz

      !x,y,z store point location in global coords
      x=seedx + repx*seednx 
      y=seedy + repy*seedny 
      z=seedz + repz*seednz 

      if (point_in_cuboid_p(x,y,z,gx,gy,gz,mysubdims)) then
        ! The point is in our subdomain.
        N(x-gx+1,y-gy+1,z-gz+1) = seedsite

        npoints = npoints+1
      end if

        end do ! repz
      end do   ! repy
    end do     ! repx

      end do ! seedz
    end do   ! seedy
  end do     ! seedx

  ! We're done with the seed file; close it.
  call xdrfclose(file_id,ierr)

  if (npoints /= nx*ny*nz) then
    print*,'Internal inconsistency: read ',npoints,' points, ', &
    'expected ',nx*ny*nz
    call Abend
  end if

end subroutine restore_upscale

subroutine restore_cutout(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  logical :: intersects_p
  character(len=128) :: cutoutfile = 'cutout.dat'
  !integer :: rnt,chkid

  character(len=128) :: seedfile ! Name of upscale restart file
  integer :: seednx,seedny,seednz ! total size of seed file
  integer :: seedstartx,seedstarty,seedstartz ! Position of read
                                        ! starting point in seed file
  integer :: seedstopx,seedstopy,seedstopz ! End of read in seed file
  integer :: putx,puty,putz       ! Position of seed origin in
                                        ! global coordinates
  logical :: pzero_p  ! Set if we need to skip the first 12 bytes
  namelist /cutout_input/ &
    seedfile, &
    pzero_p, &
    seednx,seedny,seednz, &
    seedstartx,seedstarty,seedstartz, &
    seedstopx,seedstopy,seedstopz, &
                putx,puty,putz
  integer :: gx,gy,gz ! Global coords of the start of our subdomain
  integer,dimension(3) :: mysubdims
  integer :: junkint ! Dummy to read integers into

  integer :: seedx,seedy,seedz ! Point within the seed file
  integer :: x,y,z ! Global coordinates of upscaled point
        integer :: pnx,pny,pnz ! Size of seed chunk

  integer :: npoints ! debug variable
  integer :: s ! Counter for lattice vector

  integer :: file_id,ierr ! for XDR file I/O
  type(lbe_site) :: seedsite ! Holds a single lattice site read from seed

  ! Read the cutout.dat file. This contains
  ! the total dimensions of the old lattice, and the 
  ! the position where to start putting the old
        ! system within the new one.

  call die_unless_exists(cutoutfile)
  open(unit=17,file=cutoutfile,status='UNKNOWN')
  read(unit=17,nml=cutout_input)
  close(unit=17)
        
        if (myrankc == 0 ) then
  print*,'Read cutout input:'
  print*,' seedfile=<',trim(seedfile),'> pzero_p=',pzero_p
  print*,' Seed file size is (',seednx,',',seedny,',',seednz,')'
        print*,' Seed starts at (',putx,',',puty,',',putz,') in new system'
        endif

  pnx = seedstopx-seedstartx
  pny = seedstopy-seedstarty
  pnz = seedstopz-seedstartz

        if (myrankc == 0 ) then
        print*,' Size of seed chunk is (',pnx,',',pny,',',pnz,')'
        endif

  !Find start of our subdomain in global coordinates.
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz
  mysubdims(1)=nx
  mysubdims(2)=ny
  mysubdims(3)=nz

        ! Only open the seed file if it will be placed in our domain
        ! Cool looking if clause, eh?
        if ( ((gx.le.putx).AND.((gx+nx).ge.putx)) .OR. &
           ((gx.ge.putx).AND.((gx+nx).le.(putx+pnx))) .OR. &
           ((gx.le.(putx+pnx)).AND.((gx+nx).ge.(putx+pnx))) ) then
         if ( ((gy.le.puty).AND.((gy+ny).ge.puty)) .OR. &
           ((gy.ge.puty).AND.((gy+ny).le.(puty+pny))) .OR. &
           ((gy.le.(puty+pny)).AND.((gy+ny).ge.(puty+pny))) ) then
          if ( ((gz.le.putz).AND.((gz+nz).ge.putz)) .OR. &
           ((gz.ge.putz).AND.((gz+nz).le.(putz+pnz))) .OR. &
           ((gz.le.(putz+pnz)).AND.((gz+nz).ge.(putz+pnz))) ) then

        
        
    ! Check that the seed file exists.
    call die_unless_exists(seedfile)

    ! It does; open it.
    call xdrfopen(file_id,seedfile,"r",ierr)
    call check_xdrfopen(ierr,seedfile)

    ! If pzero_p is set, then the seed file was generated by saving
    ! a checkpoint from CPU 0; skip the three header integers in this case.

    if (pzero_p) then

    call xdrfint(file_id,junkint,ierr)
    call xdrfint(file_id,junkint,ierr)
    call xdrfint(file_id,junkint,ierr)

    end if

    ! Now, run through the seed file.
    ! For each point in the seed file, find the coordinates of
    ! each of the many points in global coordinates onto which
    ! it will map after upscaling.
    !
    ! For each of these points, calculate if it intersects our
    ! subvolume. If it does, set the lattice accordingly.
    !

    npoints = 0
    do seedx=0,seednx-1
      do seedy=0,seedny-1
        do seedz=0,seednz-1

      ! Read a new lattice site in from the seed.
      do s=1,size(seedsite%n_r)
        call xdrfdouble(file_id,seedsite%n_r(s),ierr)
        call check_xdrfdouble(ierr)
      end do
#ifndef SINGLEFLUID
      do s=1,size(seedsite%n_b)
        call xdrfdouble(file_id,seedsite%n_b(s),ierr)
        call check_xdrfdouble(ierr)
      end do
#endif
#ifndef NOSURFACTANT
      do s=1,size(seedsite%n_s)
        call xdrfdouble(file_id,seedsite%n_s(s),ierr)
        call check_xdrfdouble(ierr)
      end do
      do s=1,size(seedsite%da)
        call xdrfdouble(file_id,seedsite%da(s),ierr)
        call check_xdrfdouble(ierr)
      end do
      do s=1,size(seedsite%db)
        call xdrfdouble(file_id,seedsite%db(s),ierr)
        call check_xdrfdouble(ierr)
      end do
#endif
      call xdrfdouble(file_id,seedsite%rock_state,ierr)
        call check_xdrfdouble(ierr)

                  if ((seedx.ge.(seedstartx-1)).AND.(seedx.lt.seedstopx)) then
                  if ((seedy.ge.(seedstarty-1)).AND.(seedy.lt.seedstopy)) then
                  if ((seedz.ge.(seedstartz-1)).AND.(seedz.lt.seedstopz)) then

      !x,y,z store point location in global coords
      x=seedx - (seedstartx-1) + (putx-1) 
      y=seedy - (seedstarty-1) + (puty-1) 
      z=seedz - (seedstartz-1) + (putz-1) 

      if (point_in_cuboid_p(x,y,z,gx,gy,gz,mysubdims)) then
        ! The point is in our subdomain.
        N(x-gx+1,y-gy+1,z-gz+1) = seedsite

        npoints = npoints+1
      end if

                  endif
                  endif
                  endif

         end do ! seedz
       end do   ! seedy
     end do     ! seedx

     ! We're done with the seed file; close it.
     call xdrfclose(file_id,ierr)
        
           print*,'Read ',npoints,' points for CPU ',myrankc

          else ! Z direction matches
           print*,'Seed is not part of my domain', myrankc
          endif ! end Z direction matches
         endif ! end Y direction matches
        endif ! end X direction matches
end subroutine restore_cutout

!> Set retval to true if point (x,y,z) inside the cuboid starting at
!> (xc,yc,zc) of size cubedims
function point_in_cuboid_p(x,y,z,xc,yc,zc,cubedims)
  implicit none
  integer :: x,y,z,xc,yc,zc
  integer, dimension(3) :: cubedims
  logical :: point_in_cuboid_p

  point_in_cuboid_p = &
    ( x .ge. xc )     .and. &
    ( x .lt. (xc+cubedims(1)) ) .and. &
    ( y .ge. yc )     .and. &
    ( y .lt. (yc+cubedims(2)) ) .and. &
    ( z .ge. zc )     .and. &
    ( z .lt. (zc+cubedims(3)) )
end function point_in_cuboid_p

#ifdef DIST
subroutine read_dist_xdrf_par(filename,N)
  type(lbe_site),dimension(0:,0:,0:) :: N
        integer :: gx,gy,gz ! Global coords of low-(xyz) corner of our subdomain
        integer :: ixdrs,ierr2
        integer :: x,y,z
        real(kind=rk)  :: dist !
  character(LEN=80) :: path
  Integer :: slash,xdrpos
  character(LEN=*) :: filename
  integer :: rcount=0
  character(LEN=80) :: distfilename

slash = Index(filename,"/",.TRUE.)
xdrpos = Index(filename,".xdr",.TRUE.)
path = filename(:Index(filename,"/",.TRUE.))
distfilename=trim(path)//"rockdist"//filename(slash+1:xdrpos-1)//".xdr"



        ! NB: global coordinates are zero-based!
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz




        call xdrfopen(ixdrs,distfilename,"r",ierr2)

        if (ierr2 .eq. 0) then
                print*,'read_dist_xdruchar_par: xdrfopen returned error opening distance file:',distfilename
                call Abend
        end if


        ! Seek to low-(xyz) corner of subdomain

  !call skip_xdr_double(ixdrs, gx + tny*gy + tnx*tny*gz)
  ! Previous line was buggy. Kuni, 11.05.06
  call skip_xdr_double(ixdrs, gx + tnx*gy + tnx*tny*gz)

        do z=1,nz
                do y=1,ny
                        ! Read in an X-span
                        do x=1,nx

                                call xdrfdouble(ixdrs,dist,ierr2)


                                if (ierr2 .eq. 0) then
                    print*,'rank ',myrankc,'read_dist_xdruchar_par: xdrfdouble returned error at the local coordinates',x,y,z
                                        call Abend
                                end if

                                        N(x,y,z)%abst = dist
  !print*,myrankc,x,y,z,N(x,y,z)%abst
                        end do ! x

                        ! Skip tnx-nx sites to next span in XY slab 

                        call skip_xdr_double(ixdrs,tnx-nx)

                end do ! y
                ! Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume 

                call skip_xdr_double(ixdrs, (tnx*tny)-(tnx*ny) )
        end do ! z


        call xdrfclose(ixdrs,ierr2)

! print*,'Rank ',myrankc,'rock  count=',rcount

end subroutine read_dist_xdrf_par
! endif DIST
#endif

#ifdef RELTIME
subroutine read_rel_xdrf_par(filename,N)
  type(lbe_site),dimension(0:,0:,0:) :: N
        integer :: gx,gy,gz ! Global coords of low-(xyz) corner of our subdomain
        integer :: ixdrs,ierr
#ifndef SINGLEFLUID
  integer :: ixdrs2,ierr2 
#endif
#ifndef NOSURFACTANT
  integer ::ixdrs3,ierr3
#endif
        integer :: x,y,z
        real(kind=rk) :: rel
#ifndef SINGLEFLUID
  real(kind=rk) :: rel2
#endif
#ifndef NOSURFACTANT
  real(kind=rk) :: rel3
#endif
  character(LEN=*) :: filename
  character(LEN=80) :: path
  Integer :: slash,xdrpos
  integer :: rcount=0
  character(LEN=80) :: relfilename
#ifndef SINGLEFLUID
  character(LEN=80) :: relfilename2
#endif
#ifndef NOSURFACTANT
  character(LEN=80) :: relfilename3
#endif
  integer*8 :: skipval

slash = Index(filename,"/",.TRUE.)
xdrpos = Index(filename,".xdr",.TRUE.)
path = filename(:Index(filename,"/",.TRUE.))
relfilename=trim(path)//"taupos"//filename(slash+1:xdrpos-1)//"comp1.xdr"

#ifndef SINGLEFLUID
    relfilename2=trim(path)//"taupos"//filename(slash+1:xdrpos-1)//"comp2.xdr"
#endif
#ifndef NOSURFACTANT
    relfilename3=trim(path)//"taupos"//filename(slash+1:xdrpos-1)//"comp3.xdr"
#endif


        ! NB: global coordinates are zero-based!
  gx = ccoords(1)*nx
  gy = ccoords(2)*ny
  gz = ccoords(3)*nz




        call xdrfopen(ixdrs,relfilename,"r",ierr)

#ifndef SINGLEFLUID
  call xdrfopen(ixdrs2,relfilename2,"r",ierr2)  
#endif
#ifndef NOSURFACTANT
  call xdrfopen(ixdrs3,relfilename3,"r",ierr3)  
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef SINGLEFLUID
  
#endif
#ifndef NOSURFACTANT
  
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (ierr .eq. 0) then
                print*,'read_rel_xdruchar_par: xdrfopen returned error while opening local relaxationtimes', relfilename
                call Abend
        end if
#ifndef SINGLEFLUID
if (ierr2 .eq. 0) then
                print*,'read_rel_xdruchar_par: xdrfopen returned error while opening local relaxationtimes for component 2', relfilename2
                call Abend
        end if    
#endif
#ifndef NOSURFACTANT
  if (ierr3 .eq. 0) then
                print*,'read_rel_xdruchar_par: xdrfopen returned error while opening local relaxationtimes for component 3', relfilename3
                call Abend
        end if  
#endif

        ! Seek to low-(xyz) corner of subdomain

  !call skip_xdr_double(ixdrs, gx + tny*gy + tnx*tny*gz)
  ! Previous line was buggy. Kuni, 11.05.06
  skipval=gx + tnx*gy + tnx*tny*gz
  call skip_xdr_double(ixdrs, skipval)

        do z=1,nz
                do y=1,ny
                        ! Read in an X-span
                        do x=1,nx

                                call xdrfdouble(ixdrs,rel,ierr)

                                if (ierr .eq. 0) then
                    print*,'rank ',myrankc,'read_rel_xdruchar_par: xdrfdouble returned error at local coordinates',x,y,z
                                        call Abend
                                end if

                                        N(x,y,z)%taupos_r = rel


                        end do ! x

                        ! Skip tnx-nx sites to next span in XY slab 
                        skipval=tnx-nx
                        call skip_xdr_double(ixdrs,skipval)

                end do ! y
                ! Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume 
                skipval= (tnx*tny)-(tnx*ny) 
                call skip_xdr_double(ixdrs, skipval)
        end do ! z

        call xdrfclose(ixdrs,ierr)
! print*,'Rank ',myrankc,'rock  count=',rcount


#ifndef SINGLEFLUID
        skipval= gx + tnx*gy + tnx*tny*gz
    call skip_xdr_double(ixdrs2, skipval)

        do z=1,nz
                do y=1,ny
                        ! Read in an X-span
                        do x=1,nx
                                call xdrfdouble(ixdrs2,rel2,ierr2)

                                if (ierr2 .eq. 0) then
                    print*,'rank ',myrankc,'read_rel_xdruchar_par: xdrfdouble returned error at local coordinates',x,y,z,'for component 2'
                                        call Abend
                                end if

                                        N(x,y,z)%taupos_b = rel2


                        end do ! x

                        ! Skip tnx-nx sites to next span in XY slab 
                        skipval=tnx-nx
                        call skip_xdr_double(ixdrs2,skipval)

                end do ! y
                ! Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume 
                skipval= (tnx*tny)-(tnx*ny) 
                call skip_xdr_double(ixdrs2,skipval)
        end do ! z

        call xdrfclose(ixdrs2,ierr2)
#endif

#ifndef NOSURFACTANT
        skipval= gx + tnx*gy + tnx*tny*gz
    call skip_xdr_double(ixdrs3,skipval)

        do z=1,nz
                do y=1,ny
                        ! Read in an X-span
                        do x=1,nx
                                call xdrfdouble(ixdrs3,rel3,ierr3)

                                if (ierr3 .eq. 0) then
                    print*,'rank ',myrankc,'read_rel_xdruchar_par: xdrfdouble returned error at local coordinates',x,y,z,'for component 3'
                                        call Abend
                                end if

                                        N(x,y,z)%taupos_s = rel3


                        end do ! x

                        ! Skip tnx-nx sites to next span in XY slab 
                        skipval=tnx-nx
                        call skip_xdr_double(ixdrs3,skipval)

                end do ! y
                ! Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume 
                skipval=(tnx*tny)-(tnx*ny) 
                call skip_xdr_double(ixdrs3,skipval)
        end do ! z

        call xdrfclose(ixdrs3,ierr3)
#endif

end subroutine read_rel_xdrf_par
! endif RELTIME
#endif



#endif
end module lb3d_io_xdrf_module
