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

!> Contains subroutines for the sanity check
!> (Test of conserved properties and low mach limit as well as abend on NaN)
module lb3d_sanity_module

use lb3d_lattice_module
use lb3d_config_module, only: n_sanity_check, amass_r, amass_b, amass_s
use lb3d_log_module
use lb3d_mpi_module
use lb3d_analysis_module, only: norm, massflow
use lb3d_io_module, only: lb3d_io_write_data, lb3d_io_write_checkpoint
use lb3d_io_helper_module, only: check_dump_now
use lb3d_global_module, only:n_spec

#ifndef NOIEEEARITHMETIC
  use ieee_arithmetic
#endif

  implicit none

contains

  subroutine lb3d_sanity_check
    logical :: insanity

#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lb3d_sanity_check',.false.) 
#endif
    insanity=.false.

    if ( check_dump_now(.true.,n_sanity_check) ) then
       call lbe_sanity_check(N,insanity)
       if (insanity) then
          call lb3d_io_write_data
          call lb3d_io_write_checkpoint
          call abend
       end if
    end if
  end subroutine lb3d_sanity_check

!> performs checks on the physical sanity of the LB system and dumps a
!> short summary to standard output
!>
!> The checks currently encompass looking for NaNs, looking for
!> negative distribution functions or too large velocities (and the
!> frequency of both), maximum velocity (and its location in the
!> lattice), averages of velocities and densities, and total momentum
!> and densities.
!>
!> \param[in] N local lattice chunk with halo of depth 1
!>
!> \param[in,out] insanity set to \c .true. in case of issues, never
!> set to \c .false. here
!>
!> \note All checks should be put into the same routine in order to
!> have a single loop over the whole lattice only.
subroutine lbe_sanity_check(N,insanity)
    type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
    logical,intent(inout) :: insanity
    real(kind=rk),parameter :: max_vel=0.1_rk,small=1.0e-12_rk
    logical :: found_NaN,site_too_fast,tmps
    integer :: i,j,k,s,l,ierror
    real(kind=rk) :: tot
    real(kind=rk) :: nred,nblue,nsurf,rhored,rhoblue,rhosurf,sumred,sumblue,sumsurf&
         &,sumdens
    real(kind=rk) :: p(3),sum_p(3)
    real(kind=rk), dimension(3) :: velred,velredavg,velredavgsum
    real(kind=rk), dimension(3) :: velblue,velblueavg,velblueavgsum
    real(kind=rk), dimension(3) :: velsurf,velsurfavg,velsurfavgsum
    character(len=256)     :: msgstr
    real(kind=rk) :: maxvel(n_spec,3) ! maximum velocity for all species
    ! array dimensions: species, component of vel...max, component of
    ! lattice position where vel...max has its maximum value
    integer :: maxx(n_spec,3,3),buf_pos(3)
    ! counters might exceed 32bit integer range, MPI_INTEGER8 seems to
    ! be not available everywhere (?) so use double to be safe
    real(kind=rk) :: n_fs,n_nn,n_tf,sum_fs,sum_nn,sum_tf
    real(kind=rk) :: pair(2),rpair(2) ! for MPI_MAXLOC
    character(len=4) :: spec_name(3)=(/' red','blue','surf'/)

    CALL log_msg_ws("---------( Start LBE sanity check )----------",.false.)

    ! initialize
    nred=0.
    nblue=0.
    nsurf=0.
    rhored=0.
    rhoblue=0.
    rhosurf=0.
    sumred=0.
    sumblue=0.
    sumsurf=0.
    sumdens=0.
    velredavg=0.
    velredavgsum=0.
    velblueavg=0.
    velblueavgsum=0.
    velsurfavg=0.
    velsurfavgsum=0.
    p = 0.0_rk
    sum_p = 0.0_rk
    found_NaN = .false.
    n_fs = 0.0_rk
    n_nn = 0.0_rk
    n_tf = 0.0_rk
    maxvel=0.0_rk
    maxx = 0

    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             no_rock: if (N(i,j,k)%rock_state==0.0_rk) then
                n_fs = n_fs+1.0_rk

                ! look for NaNs
                directions: do s = 1, nvecs
                   if (check_NaN(N(i,j,k)%n_r(s))) found_NaN = .true.
#ifndef SINGLEFLUID
                   if (check_NaN(N(i,j,k)%n_b(s))) found_NaN = .true.
#ifndef NOSURFACTANT
                   if (check_NaN(N(i,j,k)%n_s(s))) found_NaN = .true.
#endif
#endif
                end do directions

                ! look for negative distribution functions
                if (any(N(i,j,k)%n_r<0.0_rk)&
#ifndef SINGLEFLUID
                     &.or.any(N(i,j,k)%n_b<0.0_rk)&
#ifndef NOSURFACTANT
                     &.or.any(N(i,j,k)%n_s<0.0_rk)&
#endif
#endif
                     &) n_nn = n_nn+1.0_rk

                ! accumulate total mass
                rhored  = sum(N(i,j,k)%n_r*g)
                nred  = nred  + rhored
#ifndef SINGLEFLUID
                rhoblue = sum(N(i,j,k)%n_b*g)
                nblue = nblue + rhoblue
#endif
#ifndef NOSURFACTANT
                rhosurf = sum(N(i,j,k)%n_s*g)
                nsurf = nsurf + rhosurf
#endif

      !          print*,cx,cy,cz,amass_r

                ! accumulate avg and max velocities
                site_too_fast = .false.
                if (abs(rhored)>small) then
                   velred(1) = sum(N(i,j,k)%n_r*g*cx)*amass_r/rhored
                   velred(2) = sum(N(i,j,k)%n_r*g*cy)*amass_r/rhored
                   velred(3) = sum(N(i,j,k)%n_r*g*cz)*amass_r/rhored
                   velredavg = velredavg+velred
                   do l=1,3
                      if (abs(velred(l))>=maxvel(1,l)) then
                         maxvel(1,l) = abs(velred(l))
                         maxx(1,l,:) = (/i,j,k/)+start-1
                      end if
                   end do
                   if (norm(velred)>max_vel) site_too_fast = .true.
                end if
#ifndef SINGLEFLUID
                if (abs(rhoblue)>small) then
                   velblue(1) = sum(N(i,j,k)%n_b*g*cx)*amass_b/rhoblue
                   velblue(2) = sum(N(i,j,k)%n_b*g*cy)*amass_b/rhoblue
                   velblue(3) = sum(N(i,j,k)%n_b*g*cz)*amass_b/rhoblue
                   velblueavg = velblueavg+velblue
                   do l=1,3
                      if (abs(velblue(l))>=maxvel(2,l)) then
                         maxvel(2,l) = abs(velblue(l))
                         maxx(2,l,:) = (/i,j,k/)+start-1
                      end if
                   end do
                   if (norm(velblue)>max_vel) site_too_fast = .true.
                end if
#endif
#ifndef NOSURFACTANT
                if (abs(rhosurf)>small) then
                   velsurf(1) = sum(N(i,j,k)%n_s*g*cx)*amass_s/rhosurf
                   velsurf(2) = sum(N(i,j,k)%n_s*g*cy)*amass_s/rhosurf
                   velsurf(3) = sum(N(i,j,k)%n_s*g*cz)*amass_s/rhosurf
                   velsurfavg = velsurfavg+velsurf
                   do l=1,3
                      if (abs(velsurf(l))>=maxvel(3,l)) then
                         maxvel(3,l) = abs(velsurf(l))
                         maxx(3,l,:) = (/i,j,k/)+start-1
                      end if
                   end do
                   if (norm(velsurf)>max_vel) site_too_fast = .true.
                end if
#endif
                if (site_too_fast) n_tf = n_tf+1.0_rk

                ! accumulate total momentum
                p = p+massflow(N(i,j,k))
             end if no_rock
          end do
       end do
    end do

    if (found_NaN) then
       CALL log_msg("Found NaN",.true.)
       insanity = .true.
    end if

    ! check mass conservation:
    ! Now, add up all the subdomains' values of nred,,nblue,nsurf and place
    ! the sum in rank 0's bigsum.
    CALL MPI_Reduce(nred,sumred,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
    msgstr = 'MPI_Reduce() of red failed'
    CALL checkmpi(ierror,msgstr)

#ifndef SINGLEFLUID
    CALL MPI_Reduce(nblue,sumblue,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
    msgstr = 'MPI_Reduce() of blue failed'
    CALL checkmpi(ierror,msgstr)
#endif

#ifndef NOSURFACTANT
    CALL MPI_Reduce(nsurf,sumsurf,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
    msgstr = 'MPI_Reduce() of surf failed'
    CALL checkmpi(ierror,msgstr)
#endif
    if (myrankc == 0) then
       tot=product(tsize)
       write(msgstr,"('Total densities (r,b,g)   : ',F16.4,' ',F16.4,' ',F16.4)") &
            max(0.d-30,sumred), max(0.d-15,sumblue), max(0.d-30,sumsurf)
       CALL log_msg(trim(msgstr),.false.)
       write(msgstr,"('Average densities (r,b,g) : ',F16.10,' ',F16.10,' ',F16.10)") &
            max(0.d-30,sumred/tot), max(0.d-30,sumblue/tot), max(0.d-30,sumsurf/tot)
       CALL log_msg(trim(msgstr),.false.)
       sumdens=sumred+sumblue+sumsurf
       write(msgstr,"('Total/average density     : ',F16.4,' ',F16.10)") &
            sumdens, max(0.d-30,sumdens/tot)
       CALL log_msg(trim(msgstr),.false.)
       CALL log_ws(.false.)
    endif

    CALL MPI_Reduce(velredavg,velredavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart&
         &,ierror)
    msgstr = 'MPI_Reduce() of avg red velocity failed'
    CALL checkmpi(ierror,msgstr)
#ifndef SINGLEFLUID
    CALL MPI_Reduce(velblueavg,velblueavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart&
         &,ierror)
    msgstr = 'MPI_Reduce() of avg blue velocity failed'
    CALL checkmpi(ierror,msgstr)
#endif
#ifndef NOSURFACTANT
    CALL MPI_Reduce(velsurfavg,velsurfavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart&
         &,ierror)
    msgstr = 'MPI_Reduce() of avg red velocity failed'
    CALL checkmpi(ierror,msgstr)
#endif

    ! determine position of maximum values of velocity components for
    ! all species present
    species: do i=1,n_spec
       vel_comp: do j=1,3
          ! If we want to use the predefined MPI_MAXLOC operation in
          ! Fortran we have to pack the rank into a floating point
          ! number here.
          pair = (/maxvel(i,j),real(myrankc,kind=rk)/)
          call MPI_Allreduce(pair,rpair,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC&
               &,comm_cart,ierror)
          ! rank holding the highest velocity sends it to rank 0
          if (myrankc==0) then
             maxvel(i,j) = rpair(1)
             if (nint(rpair(2))/=0) then ! else maxx(i,j,:) is on rank 0 already
                call MPI_Recv(buf_pos,3,LBE_REAL,nint(rpair(2)),0,comm_cart&
                     &,MPI_STATUS_IGNORE,ierror)
                maxx(i,j,:) = buf_pos
             end if
          else if (myrankc==nint(rpair(2))) then
             buf_pos = maxx(i,j,:)
             call MPI_Send(buf_pos,3,LBE_REAL,0,0,comm_cart,ierror)
          end if
       end do vel_comp
    end do species

    myrankc0: if (myrankc == 0) then
       call log_msg("Maximum absolute velocity : ",.false.)
       do i=1,n_spec
          write (msgstr,"('     ',A4,' (x,y,z)         : ',3(F16.10,:,' '))") &
               &spec_name(i),max(0.0e-30_rk,maxvel(i,:))
          call log_msg(trim(msgstr),.false.)
          write (msgstr&
               &,"('      at position',9X,': ',3('(',2(I4,','),I4,')',:,X))") &
               &maxx(i,1,:),maxx(i,2,:),maxx(i,3,:)
          call log_msg(trim(msgstr),.false.)
       end do

       CALL log_msg("Average velocity          : ",.false.)
       write(msgstr,"('      red (x,y,z)         : ',F16.10,' ',F16.10,' ',F16.10)") &
            velredavgsum(1)/tot, velredavgsum(2)/tot, velredavgsum(3)/tot
       CALL log_msg(trim(msgstr),.false.)
#ifndef SINGLEFLUID
       write(msgstr,"('     blue (x,y,z)         : ',F16.10,' ',F16.10,' ',F16.10)") &
            velblueavgsum(1)/tot, velblueavgsum(2)/tot, velblueavgsum(3)/tot
       CALL log_msg(trim(msgstr),.false.)
#endif
#ifndef NOSURFACTANT
       write(msgstr,"('     surf (x,y,z)         : ',F16.10,' ',F16.10,' ',F16.10)") &
            velsurfavgsum(1)/tot, velsurfavgsum(2)/tot, velsurfavgsum(3)/tot
       CALL log_msg(trim(msgstr),.false.)
#endif
       CALL log_ws(.false.)

       do i=1,n_spec
          if (any(maxvel(i,:)>0.1_rk)) then
             call log_msg("********  WARNING: MAXIMUM "//trim(spec_name(i))&
                  &//" VELOCITY IS TOO HIGH!  ********",.false.)
!!$             insanity = .true.
          end if

       end do
    end if myrankc0

    ! total number of fluid sites
    call MPI_Reduce(n_fs,sum_fs,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
    call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_fs failed')

    call MPI_Reduce(n_tf,sum_tf,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
    call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_tf failed')
    if (myrankc==0.and.sum_tf/=0.0_rk) then
       write (msgstr&
            &,"('********  WARNING: absolute species velocity higher than ',"&
            &//"F16.12,' on ',F16.2,' lattice sites ( ',F16.12,"&
            &//"' % of total fluid sites)  ********')") &
            &max_vel,sum_tf,100.0_rk*sum_tf/sum_fs
       CALL log_msg(trim(msgstr),.false.)
    end if

    call MPI_Reduce(n_nn,sum_nn,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
    call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_nn failed')
    if (myrankc==0.and.sum_nn/=0.0_rk) then
       write (msgstr&
            &,"('********  WARNING: negative distribution functions on '"&
            &//",F16.2,' lattice sites ( ',F16.12,"&
            &//"' % of total fluid sites)  ********')") &
            &sum_nn,100.0_rk*sum_nn/sum_fs
       CALL log_msg(trim(msgstr),.false.)
    end if

    ! total momentum
    call MPI_Reduce(p,sum_p,3,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
    call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of momentum failed')
    if (myrankc==0) then
       write (msgstr,"('Total momentum (x,y,z)    : ',3(F16.8,:,' '))") sum_p
       CALL log_msg(trim(msgstr),.false.)
    end if

    tmps = insanity
    CALL mpi_allreduce(tmps,insanity,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD&
         &,ierror)
    CALL checkmpi(ierror,'MPI_Allreduce() of insanity check failed')

    CALL log_msg_ws("----------( End LBE sanity check )-----------",.false.)
end subroutine lbe_sanity_check


logical function check_NaN(n)
  implicit none
  real(kind=rk) :: n
#ifdef NOIEEEARITHMETIC
#ifndef NOISNAN
  check_NaN = isnan(n)
#else
  ! UNABLE TO CHECK FOR NANs!
  check_NaN = .false.
#endif
#else
  check_NaN = ieee_is_nan(n)
#endif
  return
end function check_NaN

subroutine report_check_NaN_function(str)
  implicit none
  character(len=*), intent(out) :: str
#ifdef NOIEEEARITHMETIC
#ifndef NOISNAN
  str = "  NOIEEEARITHMETIC is set, NaNs are checked through function 'isnan'."
#else
  str = "  WARINING: NOIEEEARITHMETIC is set and NOISNAN is set, NaNs are not checked."
#endif
#else
  str = "  NOIEEEARITHMETIC is not set, NaNs are checked through function 'ieee_is_nan'."
#endif
end subroutine report_check_NaN_function


end module lb3d_sanity_module
