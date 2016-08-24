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

!> Calculation and file output of stress in the system.
module lb3d_io_stress_module

  use lb3d_config_module
  use lb3d_global_module
  use lb3d_mpi_module
  use lb3d_helper_module
  use lb3d_io_helper_module

  implicit none
  private
  public dump_stress, local_fluid_momentum_transfer, setup_stress 

contains

  subroutine dump_stress(N)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk) :: stress(2),stress_sum(2)
    character(len=256) :: msgstr
    character(len=256) cfg_file_name
    integer,parameter :: cfg_file_unit=12
    logical :: fexist
    integer ierror

    stress(:) = 0.0d0
    
    call local_fluid_momentum_transfer(N,stress)

    ! The momentum transfer through each of both planes is caused
    ! by forces from both sides. So the values for both planes are
    ! multiplied by one half here in order to let their sum be the
    ! force from both sides times one timestep.
    stress(:) = stress(:)*0.5_rk

    call MPI_Reduce(stress,stress_sum,2,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)

    rank0: if (myrankc==0) then

      call lbe_make_filename_append(cfg_file_name,'stress','.asc')

      inquire(file=cfg_file_name,exist=fexist)
      if (fexist) then
        open(unit=cfg_file_unit,file=cfg_file_name,status='OLD'&
             &,position='APPEND',recl=650)
      else
        open(unit=cfg_file_unit,file=cfg_file_name,status='NEW'&
             &,position='APPEND',recl=650)
      endif

      write (unit=cfg_file_unit, fmt='(I10.10,1X)',advance='no') nt
      write (unit=cfg_file_unit, fmt='(4(F20.12,X))',advance='no') (-stress_sum(1)+stress_sum(2))/(real(tnz)*real(tny)),-stress_sum(1)+stress_sum(2),stress_sum(1), stress_sum(2)
      write (unit=cfg_file_unit, fmt='()',advance='yes')
      close (cfg_file_unit)
    end if rank0

  end subroutine dump_stress

  !> Adds to \c lfmt(1:2) the momentum transfer through the layers
  !> \c x==dfz_minx and \c x==dfz_maxx (as far as part of the local
  !> chunk) caused by the fluid itself.
  !>
  !> \todo This does neigher work for
  !> \c amass_(r|b|s)/=1.0 nor without \c SINGLEFLUID!
  subroutine local_fluid_momentum_transfer(N,lfmt)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk),intent(out) :: lfmt(2)
    ! vectors that are advected in x direction AND carry momentum in z
    ! direction (p:positive x-direction, m:negative x-direction)
    integer,parameter,dimension(2) :: pvecs=(/9,10/),mvecs=(/13,14/)
    integer :: local_stress_minx,local_stress_maxx
    integer :: r,s,t,y,z
    character(len=256) :: msgstr

    lfmt = 0.d0

    local_stress_minx = stress_minx+1-start(1)
    if (1<=local_stress_minx.and.local_stress_minx<=nx) then
      ! write(msgstr,"('prev stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(trim(msgstr),.true.)
      ! loop through lattice nodes on plane x=dfz_minx
      ! call log_msg("Looping through plane 1",.true.)
      do y=1,ny
        do z=1,nz
          do t=1,size(pvecs)
            s = pvecs(t)
            r = bounce(s)
            ! both source and target node have to be fluid nodes
            if (N(local_stress_minx-1,y,z-cz(s))%rock_state==0.0_rk&
                 &.and.N(local_stress_minx,y,z)%rock_state==0.0_rk) then
              ! sum incoming and outgoing momentum
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
            end if
          end do
        end do
      end do
      ! write(msgstr,"('summed stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(trim(msgstr),.true.)
    end if

    local_stress_maxx = stress_maxx+1-start(1)
    if (1<=local_stress_maxx.and.local_stress_maxx<=nx) then
      ! write(msgstr,"('prev stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(trim(msgstr),.true.)
      ! loop through lattice nodes on plane x=dfz_maxx
      ! call log_msg("Looping through plane 2",.true.)
      do y=1,ny
        do z=1,nz
          do t=1,size(pvecs)
            s = mvecs(t)
            r = bounce(s)
            ! both source and target node have to be fluid nodes
            if (N(local_stress_maxx+1,y,z-cz(s))%rock_state==0.0_rk&
                 &.and.N(local_stress_maxx,y,z)%rock_state==0.0_rk) then
              ! sum incoming and outgoing momentum
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
            end if
          end do
        end do
      end do
      ! write(msgstr,"('summed stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(trim(msgstr),.true.)
    end if

  end subroutine local_fluid_momentum_transfer

  !> setup stress boundary coordinates
  subroutine setup_stress()
    implicit none
    if (stress_minx<0) stress_minx = 1
    if (stress_maxx<0) stress_maxx = tnx
  end subroutine setup_stress

end module lb3d_io_stress_module
