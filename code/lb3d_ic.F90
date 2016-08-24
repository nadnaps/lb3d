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

!> Contains subroutines implementing the initial conditions
module lb3d_ic_module
!      ==================================================================

  use lb3d_config_module,only:init_cond,restore,n_restore,nt,seed,fr,fb,fg,MRT,d_0
  use lb3d_global_module, only:rk,pi,tstart
  use lb3d_lattice_module
  use lb3d_bdist_module
  use lb3d_log_module
  use lb3d_init_rock_module, only: lbe_init_rock


  use lb3d_io_checkpoint_module, only:restore_checkpoint
#ifdef USEXDRF
  use lb3d_io_checkpoint_module,only: restore_upscale,restore_cutout
#endif

  implicit none

  type radial_model
     real(kind=rk) :: n_r, n_b, n_s
     integer :: dip
  end type radial_model
  
  !> \{
  !> \name indices into the radial model array for lbe_init_radial
  integer, parameter :: rad_inner=1,rad_middle=2,rad_outer=3
  !> \}

contains


  !> Calling subroutines to initialise the densities on the lattice \c N
  subroutine lb3d_ic_init (stage)


    integer :: stage

#ifdef LB3D_DEBUG_INFO    
    write(msgstr,"('In lb3d_ic_init stage',I0)") stage
    CALL log_msg(trim(msgstr),.false.) 
#endif

    select case (stage)

    case (4) ! post alloc

       ! Initialising the densities to zero
       call lb3d_lattice_clear(N)

       tstart = 0

       call lb3d_ic_apply(N,init_cond)

       nt = tstart

       ! Add perturbation
       if ( ( init_cond .ne. 7 ) .and. ( .not. restore ) ) then 
          ! Don't do this if restoring!
          ! Add a perturbation if required.
          if (perturbation .ne. 0.0) then
             write(msgstr,"('Perturbing by fraction ',F16.10)") perturbation
             CALL log_msg(trim(msgstr),.false.)
             CALL lbe_init_perturb(N)
          end if

          ! If shearing then remove any net x and z momentum
          if ( (boundary_cond .ne. 0) .and. &
               (inv_fluid .eq. 5 .or. inv_fluid .eq. 6) ) then
             CALL lbe_init_remove_momentum(N)
          end if
       end if
 
       call lbe_init_rock(N)
    case default

#ifdef LB3D_DEBUG_INFO
       call log_msg('nothing to do.',.false.)
#endif
    end select

  end subroutine lb3d_ic_init

  subroutine lb3d_ic_apply(N,init_cond)
    integer,intent(in) :: init_cond
    type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
    type(radial_model),dimension(3) :: model

    select case(init_cond)   
    case (-5)
        CALL log_msg("Initializing system using B_DIST with Z velocity.",.false.)
        CALL lbe_init_velz(N)

     case (-4)
        CALL log_msg("Initializing system using B_DIST.",.false.)
        CALL lbe_init_simple(N)
        
     case (-3)
        CALL log_msg("Initializing system using INIT_RATIO.",.false.)
        CALL lbe_init_ratio(N)
        
     case (-1)
        CALL log_msg("Initializing system using INIT_RAND.",.false.)
        CALL lbe_init_rand(N)
        
     case (0)
        CALL log_msg("Initializing system using INIT_FRAC.",.false.)
        CALL lbe_init_frac(N)
        
     case (1)
        CALL log_msg("Initializing system using INIT_LAMX.",.false.)
        CALL lbe_init_lamx(N)
        
     case (2)
        CALL log_msg("Initializing system using INIT_LAMY.",.false.)
        CALL lbe_init_lamy(N)
        
     case (3)
        CALL log_msg("Initializing system using INIT_LAMZ.",.false.)
        CALL lbe_init_lamz(N)

      case (7)
        CALL log_msg("Initializing system using INIT_RESTORE.",.false.)
        CALL lbe_init_restore(N)

      case (9)
        CALL log_msg("Initializing system using INIT_UPSCALE.",.false.)
        CALL lbe_init_upscale(N)

      case (10)
        CALL log_msg("Initializing system using INIT_RATIO and INIT_CUTOUT",.false.)
        CALL lbe_init_ratio(N)
#ifdef USEXDRF
        CALL restore_cutout(N)
#else
        CALL log_msg("XDRF is disabled, can't perform INIT_CUTOUT",.false.)
#endif

     case (11)
        CALL log_msg("Initializing system using INIT_BI_LAM_X.",.false.)
        CALL lbe_init_bi_lam_x(N)

     case (12)
        CALL log_msg("Initializing system using INIT_BI_LAM_Y.",.false.)
        CALL lbe_init_bi_lam_y(N)

     case (13)
        CALL log_msg("Initializing system using INIT_BI_LAM_Z.",.false.)
        CALL lbe_init_bi_lam_z(N)

    case (14) 
       CALL log_msg("Initializing system using INIT_RADIAL.",.false.)
       model(rad_inner)%n_r = fr
       model(rad_inner)%n_b = fb 
       model(rad_inner)%n_s = fg
       model(rad_inner)%dip = fd
       
       model(rad_middle)%n_r = qr 
       model(rad_middle)%n_b = qb 
       model(rad_middle)%n_s = qg
       model(rad_middle)%dip = qd

       model(rad_outer)%n_r = pr 
       model(rad_outer)%n_b = pb 
       model(rad_outer)%n_s = pg
       model(rad_outer)%dip = pd
       
       CALL lbe_init_radial(N,model)

    case (16)
       CALL log_msg("Initializing system using INIT_BI_SIN_LAM_X.",.false.)
       CALL lbe_init_bi_sin_lam_x(N)

    case default
      write(msgstr,"('FATAL ERROR: Unknown init_cond <',I0,'>. Aborting...')") init_cond
      CALL log_msg(trim(msgstr),.false.)
      CALL Abend
   end select     

  end subroutine lb3d_ic_apply


!>Restores the entire lattice state from one previously saved.
subroutine lbe_init_restore(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N

  CALL restore_checkpoint(n_restore, chk_uid, N)
  tstart = n_restore
 

end subroutine lbe_init_restore

  !>Builds a new lattice state by cloning copies of a smaller
  !>saved system. 
  subroutine lbe_init_upscale(N)
    implicit none
    type(lbe_site),dimension(0:,0:,0:) :: N
    
#ifdef USEXDRF
    CALL restore_upscale(n_restore, chk_uid, N)
#else
    CALL log_msg("XDRF is disabled, can't perform INIT_UPSCALE",.false.)
#endif
  end subroutine lbe_init_upscale

!> The mean x and z momentum are forced to be 0 - for Lees Edwards code to work.
subroutine lbe_init_remove_momentum_zx(N)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: s
  real(kind=rk) :: sumposx,sumnegx,sumposz,sumnegz
  
  sumposx = 0.0_rk
  sumnegx = 0.0_rk
  sumposz = 0.0_rk
  sumnegz = 0.0_rk
  
  do s = 1,nnp
     sumnegz = sumnegz + sum(N%n_r(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_r(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_r(posz(s)) = N%n_r(posz(s))*sumnegz/sumposz
  enddo
  
#ifndef NOSURFACTANT
  sumnegz = 0.0_rk
  sumposz = 0.0_rk
  
  do s = 1,nnp
     sumnegz = sumnegz + sum(N%n_s(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_s(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_s(posz(s)) = N%n_s(posz(s))*sumnegz/sumposz
  enddo
#endif
  
#ifndef SINGLEFLUID
  sumnegz = 0.0_rk
  sumposz = 0.0_rk
  
  do s=1,nnp
     sumnegz = sumnegz + sum(N%n_b(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_b(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_b(posz(s)) = N%n_b(posz(s))*sumnegz/sumposz
  enddo
#endif

  ! For Lees-Edwards x-momentum must sum to zero
  ! This is a mess - but it only runs once & I have other problems
  
  do s = 1,nnp
     sumnegx = sumnegx + sum(N%n_r(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_r(posx(s))*g(posx(s)))
  enddo
  ! Clip to avoid division by zero.
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_r(posx(s)) = N%n_r(posx(s))*sumnegx/sumposx
  enddo
  
  sumnegx = 0.d0
  sumposx = 0.d0
  
#ifndef NOSURFACTANT
  do s = 1,nnp
     sumnegx = sumnegx + sum(N%n_s(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_s(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_s(posx(s)) = N%n_s(posx(s))*sumnegx/sumposx
  enddo
#endif
  
#ifndef SINGLEFLUID
  sumnegx = 0.0_rk
  sumposx = 0.0_rk
  
  do s=1,nnp
     sumnegx = sumnegx + sum(N%n_b(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_b(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_b(posx(s)) = N%n_b(posx(s))*sumnegx/sumposx
  enddo
#endif
  ! It is nice if mean z-momentum is zero as then graph
  ! intercepts at +/- shear_u

end subroutine lbe_init_remove_momentum_zx

!> The mean x and z momentum are forced to be 0 - for Lees Edwards code to work.
subroutine lbe_init_remove_momentum(N)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: s
  real*8 :: sumposx,sumnegx,sumposz,sumnegz
  
  sumposx = 0.0_rk
  sumnegx = 0.0_rk
  sumposz = 0.0_rk
  sumnegz = 0.0_rk
  
  ! For Lees-Edwards x-momentum must sum to zero
  ! This is a mess - but it only runs once & I have other problems
  
  do s = 1,nnp
     sumnegx = sumnegx + sum(N%n_r(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_r(posx(s))*g(posx(s)))
  enddo
  ! Clip to avoid division by zero.
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_r(posx(s)) = N%n_r(posx(s))*sumnegx/sumposx
  enddo
  
  sumnegx = 0.d0
  sumposx = 0.d0
  
#ifndef NOSURFACTANT
  do s = 1,nnp
     sumnegx = sumnegx + sum(N%n_s(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_s(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_s(posx(s)) = N%n_s(posx(s))*sumnegx/sumposx
  enddo
#endif
  
#ifndef SINGLEFLUID
  sumnegx = 0.d0
  sumposx = 0.d0
  
  do s=1,nnp
     sumnegx = sumnegx + sum(N%n_b(negx(s))*g(negx(s)))
     sumposx = sumposx + sum(N%n_b(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
     N%n_b(posx(s)) = N%n_b(posx(s))*sumnegx/sumposx
  enddo
#endif
  ! It is nice if mean z-momentum is zero as then graph
  ! intercepts at +/- shear_u
  
  sumnegz = 0.d0
  sumposz = 0.d0
  
  do s = 1,nnp
     sumnegz = sumnegz + sum(N%n_r(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_r(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_r(posz(s)) = N%n_r(posz(s))*sumnegz/sumposz
  enddo
  
#ifndef NOSURFACTANT
  sumnegz = 0.d0
  sumposz = 0.d0
  
  do s = 1,nnp
     sumnegz = sumnegz + sum(N%n_s(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_s(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_s(posz(s)) = N%n_s(posz(s))*sumnegz/sumposz
  enddo
#endif
  
#ifndef SINGLEFLUID
  sumnegz = 0.d0
  sumposz = 0.d0
  
  do s=1,nnp
     sumnegz = sumnegz + sum(N%n_b(negz(s))*g(negz(s)))
     sumposz = sumposz + sum(N%n_b(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
     N%n_b(posz(s)) = N%n_b(posz(s))*sumnegz/sumposz
  enddo
#endif
end subroutine lbe_init_remove_momentum

!=========================================================================


subroutine lbe_init_velz(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r
  integer :: x, y, z
  
  do x = 1,nx
     do y = 1,ny
        do z = 1,nz
           if (COLLISIONTYPE_ID .eq. MRT) then
            !  CALL mrt_init_dist((/0.0_rk,0._rk,pr/),fr,N(x,y,z)%n_r)
           else
              CALL boltz_dist(0.0_rk,0.0_rk,pr,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,N(x,y,z)%n_r(:))
              N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*fr
           end if
#ifndef SINGLEFLUID
           if (COLLISIONTYPE_ID .eq. MRT) then
            !  CALL mrt_init_dist((/0.0_rk,0._rk,pb/),fb,N(x,y,z)%n_b)
           else
              CALL boltz_dist(0.0_rk,0.0_rk,pb,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,N(x,y,z)%n_b(:))
              N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*fb
           end if
#endif
#ifndef NOSURFACTANT
           if (COLLISIONTYPE_ID .eq. MRT) then
             ! CALL mrt_init_dist((/0.0_rk,0._rk,pg/),fg,N(x,y,z)%n_s)
           else
              CALL boltz_dist(0.0_rk,0.0_rk,pg,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,N(x,y,z)%n_s(:))
              N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)*fg
           end if

           if ( fd .eq. 0 ) then
              CALL random_number(theta)
              CALL random_number(phi)
              CALL random_number(r)
              theta = theta*pi ! 0 <= theta <= pi
              phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
              N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
              N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
              N(x,y,z)%da(3) = r*cos(theta)

              N(x,y,z)%db(:) =  N(x,y,z)%da(:)
           else
              N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              N(x,y,z)%db(:) =  N(x,y,z)%da(:)
           end if
#endif
        end do
     end do
  end do
end subroutine lbe_init_velz


subroutine lbe_init_simple(N)
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: x, y, z

! Now initialize my subdomain.

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     call boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F)
  end if

  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           if ( inv_fluid .eq. 11 ) then
              ! calculating a lineal interpolation of fr and pr for
              ! the initial condition for the density
              N(x,y,z)%n_r(:) = F(:)*( (pr - fr)*( start(3) + z - 2 )/(tnz - 1) + fr )
           else
              N(x,y,z)%n_r(:) = F(:)*fr
           end if
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
           N(x,y,z)%n_s(:) = F(:)*fg
           
           if ( fd .eq. 0 ) then
              CALL random_number(theta)
              CALL random_number(phi)
              CALL random_number(r)
              theta = theta*pi ! 0 <= theta <= pi
              phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
              N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
              N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
              N(x,y,z)%da(3) = r*cos(theta)

              N(x,y,z)%db(:) =  N(x,y,z)%da(:)
           else
              N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              N(x,y,z)%db(:) =  N(x,y,z)%da(:)
           end if
#endif
        end do
     end do
  end do
end subroutine lbe_init_simple


subroutine lbe_init_ratio(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: x, y, z
  
  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
           N(x,y,z)%n_s(:) = F(:)*fg

           CALL random_number(theta)
           CALL random_number(phi)
           CALL random_number(r)
           theta = theta*pi ! 0 <= theta <= pi
           phi = phi * 2.d0 * pi ! 0 <=  phi  <= 2pi
           N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
           N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
           N(x,y,z)%da(3) = r*cos(theta)
           
           N(x,y,z)%db(:) =  N(x,y,z)%da(:)
#endif
        end do
     end do
  end do
end subroutine lbe_init_ratio


subroutine lbe_init_rand(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r
  integer :: x, y, z
  real(kind=rk), dimension(nvecs) :: g_inv
  g_inv = 1./g
  
  CALL log_msg("Initializing subdomain.",.false.)
  do x = 1,nx
     do y = 1,ny
        do z = 1,nz
           CALL random_number(N(x,y,z)%n_r(:))
           N(x,y,z)%n_r(:) = fr*N(x,y,z)%n_r(:)*g_inv/nvecs
#ifndef SINGLEFLUID
           CALL random_number(N(x,y,z)%n_b(:))
           N(x,y,z)%n_b(:) = fb*N(x,y,z)%n_b(:)*g_inv/nvecs
#endif

#ifndef NOSURFACTANT
           CALL random_number(N(x,y,z)%n_s(:))
           N(x,y,z)%n_s(:) = fg*N(x,y,z)%n_s(:)*g_inv/nvecs

           CALL random_number(theta)
           CALL random_number(phi)
           CALL random_number(r)	! 0 <=   r   <= 1
           theta = theta*pi	! 0 <= theta <= pi
           phi = phi*2.d0*pi	! 0 <=  phi  <= 2pi

           N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
           N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
           N(x,y,z)%da(3) = r*cos(theta)
           
           N(x,y,z)%db(:) =  N(x,y,z)%da(:)
#endif
        end do
     end do
  end do
end subroutine lbe_init_rand


subroutine lbe_init_frac(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r, s, F(19)
  integer :: x, y, z
  real(kind=rk), dimension(nvecs) :: g_inv

  g_inv = 1.d0/g

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
         CALL random_number(N(x,y,z)%n_r(:))
         N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)/sum( N(x,y,z)%n_r(:)*g(:) )
         CALL random_number(s)
         N(x,y,z)%n_r(:) = fr*s*N(x,y,z)%n_r(:)
#ifndef SINGLEFLUID
         CALL random_number(N(x,y,z)%n_b(:))
         N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)/sum( N(x,y,z)%n_b(:)*g(:) )
         CALL random_number(s)
         N(x,y,z)%n_b(:) = fb*s*N(x,y,z)%n_b(:)
#ifndef NOSURFACTANT
         CALL random_number(N(x,y,z)%n_s(:))
         N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)/sum( N(x,y,z)%n_s(:)*g(:) )
         CALL random_number(s)
         N(x,y,z)%n_s(:) = fg*s*N(x,y,z)%n_s(:)

         CALL random_number(theta)
         CALL random_number(phi)
         CALL random_number(r)
         theta = theta*pi ! 0 <= theta <= pi
         phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

         N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
         N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
         N(x,y,z)%da(3) = r*cos(theta)
         
         N(x,y,z)%db(:) =  N(x,y,z)%da(:)

#endif
#endif
      end do
    end do
 end do
end subroutine lbe_init_frac


!> \{
!> \name Init subroutines for ternary systems
!> lbe_init_lamx
!> lbe_init_lamy
!> lbe_init_lamz
!> lbe_init_radial
!>
!> N.B. Be careful with the precompiler ifdefs if adding/removing
!> subroutines from here.

!> Sets up the system to have lamellae perpendicular to the X-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamx(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tx	! x-coordinate in global system.
 
  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.
  
  lamwidth = fr1 + fr2 + 2
  
  if (myrankc == 0) then
     if (mod(tnx,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif

  ! Note that tnx is zero-based.
  do x=1,nx
     tx = x-1 + ccoords(1)*nx
     i=mod(tx,lamwidth)
     
     if (i<fr1) then
        ! Make oil happen.
        do z = 1,nz
           do y = 1,ny
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg
              
              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
                 
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)

              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)
              end if

#endif
           end do
        end do
     else if ( (i == fr1) .or. (i == (fr1+fr2+1)) ) then
        ! Make surf happen.
        do z = 1,nz
          do y = 1,ny
             N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
             N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT            
             N(x,y,z)%n_s(:) = F(:)*qg
      
              if ( qd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
                 
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)

              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)
              end if

#endif
          end do
       end do
    elseif (i < (fr1 + 1 + fr2) ) then
       ! Make water happen.
       do z = 1,nz
          do y = 1,ny
             N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
             N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%n_s(:) = F(:)*pg

              if ( pd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
                 
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                 N(x,y,z)%db(:) =  N(x,y,z)%da(:)
              end if

#endif
          end do
       end do
    else
       print*,"Can't happen in lam-x!"
       call Abend
    end if
 end do
end subroutine lbe_init_lamx

!> Sets up the system to have lamellae perpendicular to the Y-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamy(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y ,z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: ty	! y-coordinate in global system.
  
  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.
  
  lamwidth = fr1 + fr2 + 2
  
  if (myrankc == 0) then
     if (mod(tny,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif
  
  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

 ! Note that tny is zero-based.
  do y = 1,ny
     ty = y - 1 + ccoords(2)*ny
     i = mod(ty,lamwidth)
     
     if ( i<fr1 ) then
        ! Make oil happen.
        do x = 1,nx
           do z = 1,nz
                N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
                N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%n_s(:) = F(:)*fg
                
                if ( fd .eq. 0 ) then
                   CALL random_number(theta)
                   CALL random_number(phi)
                   CALL random_number(r)
                   theta = theta*pi ! 0 <= theta <= pi
                   phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

                   N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                   N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                   N(x,y,z)%da(3) = r*cos(theta)
                   
                   N(x,y,z)%db(:) =  N(x,y,z)%da(:)
                else
                   N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                   N(x,y,z)%db(:) =  N(x,y,z)%da(:)
                end if

#endif
             end do
          end do
       elseif ( (i == fr1) .or. (i==(fr1+fr2+1)) ) then
        ! Make surf happen.
          do x = 1,nx
             do z = 1,nz
                N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
                N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%n_s(:) = F(:)*qg

                if ( qd .eq. 0 ) then
                   CALL random_number(theta)
                   CALL random_number(phi)
                   CALL random_number(r)
                   theta = theta*pi ! 0 <= theta <= pi
                   phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi

                   N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                   N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                   N(x,y,z)%da(3) = r*cos(theta)
                   
                   N(x,y,z)%db(:) =  N(x,y,z)%da(:)
                else
                   N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                   N(x,y,z)%db(:) =  N(x,y,z)%da(:)
                end if

#endif
             end do
          end do
       elseif ( i < (fr1 + 1 + fr2) ) then
          ! Make water happen.
          do x = 1,nx
             do z = 1,nz
                N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
                N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%n_s(:) = F(:)*pg

                if ( pd .eq. 0 ) then
                   CALL random_number(theta)
                   CALL random_number(phi)
                   CALL random_number(r)
                   theta = theta*pi ! 0 <= theta <= pi
                   phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                   N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                   N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                   N(x,y,z)%da(3) = r*cos(theta)
                else
                   N(x,y,z)%da(:) = (/ 0., 1., 0. /)
                end if
                N(x,y,z)%db = N(x,y,z)%da
#endif
             end do
          end do
       else
          print*,"Can't happen in lam-y!"
          call Abend
       end if
    end do
  end subroutine lbe_init_lamy

!> 
!> Sets up the system to have lamellae perpendicular to the Y-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamz(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tz	! z-coordinate in global system.
  
  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.
  
  lamwidth = fr1 + fr2 + 2
  
  if (myrankc == 0) then
     if (mod(tnz,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif
 
  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if
  
  ! Note that tnz is zero-based.
  do z=1,nz
     tz = z-1 + ccoords(3)*nz
     i=mod(tz,lamwidth)
     print*,'z=',z,'tz=',tz,'i=',i
     
     if ( i<fr1 ) then
        ! Make oil happen.
        print*,' oil'
        do x = 1,nx
           do y = 1,ny
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg

              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da

#endif
           end do
        end do
     elseif ((i == fr1) .or. (i==(fr1+fr2+1))) then
        ! Make surf happen.
        print*,' surf'
        do x = 1,nx
           do y = 1,ny
              N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*qg

              if ( qd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da
#endif
           end do
        end do
     elseif (i < (fr1 + 1 + fr2) ) then
        ! Make water happen.
        print*,' water'
        do x=1,nx
           do y=1,ny
              N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*pg

              if ( pd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da
#endif
           end do
        end do
     else
        print*,"Can't happen in lam-z!"
        call Abend
     end if
  end do
end subroutine lbe_init_lamz


!> This is a general routine for initialising droplike or vesicular
!> initial systems.  \c model contains an \c lbe_site type for each of
!> the three regions xE<lt>fr1, fr1E<lt>xE<gt>fr2, and xE<gt>fr2.
!>
!> Each site is initialised to the appropriate value.
!>
!>However, the dipole moment is ignored.
!>
!>The dipole moment of each site is set to have unit value, and
!>its orientation depends on the rock_state value for the region.
!>If \c rock_state is negative, the dipole points towards the centre of
!>the system (ie -r); if positive, it points away from the centre.
!>If zero, its orientation is chosen randomly.
!>
!>The ME3D initial conditions are all special cases of this model.
!>
!>Currently, the values are averaged over the eight points of the unit
!>cube. This sucks.
subroutine lbe_init_radial(N,model)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: x, y, z, s
  integer, dimension(3) :: centre     ! Global coords of centre of system.
  integer, dimension(3) :: base	      ! Global coords of corner of subdomain
  integer, dimension(3) :: offset
  integer, dimension(8,3) :: vertices ! Vertices of a unit cube
  real(kind=rk) :: if1, if2	              ! Real radii
 
  real(kind=rk), dimension(3) :: r	! Global coords of a given point
  real(kind=rk) :: rad			! Radius of point
  integer :: i		! Loop variable over vertices
  
  real(kind=rk) :: n_r
#ifndef SINGLEFLUID
  real(kind=rk) :: n_b
#endif
#ifndef NOSURFACTANT
  real(kind=rk) :: n_s
  real(kind=rk), dimension(3) :: d, dd
#endif
  real(kind=rk), parameter :: A = 1./8.	
  real(kind=rk) :: F(19)
  type(radial_model) :: model(3)  
  integer :: rad_index
#ifndef NOSURFACTANT
  real(kind=rk) :: theta, phi
#endif
  character(len=128) :: msgstr
  
  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! All coordinates in this routine are in the global
  ! coordinate system.

  ! Define the actual radii of the spheres.
  if1 = fr1*min(tnx, tny, tnz)/2.d0
  if2 = fr2*min(tnx, tny, tnz)/2.d0
  
  !Define the centre of the system.
  centre = (/ tnx, tny, tnz /)/2 + (/ drop_xshift, drop_yshift, drop_zshift /)

  ! Define the minimum xyz corner of my subdomain.
  base = ccoords*(/ nx, ny, nz /)
  ! Subtract the two.
  offset = base - centre

  ! Define the vectors of a unit cube.
  vertices(1,:) = (/ 0, 0, 0 /)
  vertices(2,:) = (/ 0, 0, 1 /)
  vertices(3,:) = (/ 0, 1, 0 /)
  vertices(4,:) = (/ 0, 1, 1 /)
  vertices(5,:) = (/ 1, 0, 0 /)
  vertices(6,:) = (/ 1, 0, 1 /)
  vertices(7,:) = (/ 1, 1, 0 /)
  vertices(8,:) = (/ 1, 1, 1 /)

  write(msgstr,"('if1 = ',F16.10,' if2 = ',F16.10)") if1, if2
  CALL log_msg(trim(msgstr),.false.)
  
  ! Note that r is zero-based.
  
  do x = 1,nx
     do y = 1,ny
        do z = 1,nz
           n_r = 0.d0
#ifndef SINGLEFLUID
           n_b = 0.d0
#endif
#ifndef NOSURFACTANT
           n_s = 0.d0
           d = 0.d0
#endif
           do i = 1,8 ! Average over each corner of the site cube.
              r = (/ x - 1, y - 1, z - 1 /) + offset + vertices(i,:) 
              rad = sqrt(sum(r*r))
              
              ! Find the appropriate index.
              if ( rad < if1 ) then
                 rad_index = rad_inner
              elseif ( (rad >= if1) .and. (rad <= if2) ) then
                 rad_index = rad_middle
              else
                 rad_index = rad_outer
              endif
              
              ! Set appropriate values.
              n_r = n_r + model(rad_index)%n_r
#ifndef SINGLEFLUID
              n_b = n_b + model(rad_index)%n_b
#endif
#ifndef NOSURFACTANT
              n_s = n_s + model(rad_index)%n_s
           
              ! Handle the dipole moment.
              ! (but only if there are actually any dipoles..)              
              if ( model(rad_index)%n_s > 0.d0 ) then
                 if ( model(rad_index)%dip < 0.d0 ) then                  
                    ! Dipole points towards centre.                
                    d = d - r/rad                    
                 elseif ( model(rad_index)%dip == 0) then
                    ! Make a random dipole of unit size.
                    call random_number(theta)
                    call random_number(phi)
                    theta = theta*pi	! 0 <= theta <= pi
                    phi = phi*2.0d0*pi	! 0 <= phi <= 2pi
                    dd(1) = sin(theta)*cos(phi)
                    dd(2) = sin(theta)*sin(phi)
                    dd(3) = cos(theta)
                    
                    ! Add it in to the averaging process.
                    d = d + dd
                 else
                    ! Dipole points away from centre
                    d = d + r/rad
                 endif
              endif
#endif
           end do

           ! n_r, n_s, n_b are the required average particle densities for one site.

           N(x,y,z)%n_r = A*n_r*F(:)
#ifndef SINGLEFLUID
           N(x,y,z)%n_b = A*n_b*F(:)
#endif
#ifndef NOSURFACTANT
           N(x,y,z)%n_s = A*n_s*F(:)
           N(x,y,z)%da = A*d

           N(x,y,z)%db = N(x,y,z)%da

#endif
           ! N(x,y,z)%rock_state=rad_index ! Useful for debugging.
        end do
     end do
  end do
end subroutine lbe_init_radial

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the y direction.
subroutine lbe_init_bi_lam_x(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tx	! y-coordinate in global system.


  lamwidth = fr1 + fr2
  
  if (myrankc == 0) then
     if (mod(tny,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! Note that ty is zero-based.
  do x=1,nx
     tx = x - 1 + ccoords(1)*nx
     i = mod(tx,lamwidth)
     
     if ( i<fr1 ) then
        ! Make oil happen.
        do y = 1,ny
           do z = 1,nz
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg

              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da

#endif
           end do
        end do
     else 
        ! Make water happen.
        do y = 1,ny
           do z = 1,nz
              N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*pg

              if ( qd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da

#endif
           end do
        end do
     end if
  end do
end subroutine lbe_init_bi_lam_x


!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the y direction.
subroutine lbe_init_bi_lam_y(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: ty	! y-coordinate in global system.

  lamwidth = fr1 + fr2
  
  if ( myrankc == 0 ) then
     if (mod(tny,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! Note that ty is zero-based.
  do y = 1,ny
     ty = y-1 + ccoords(2)*ny
     i = mod(ty,lamwidth)
     
     if (i<fr1) then
        ! Make oil happen.
        do x = 1,nx
           do z = 1,nz
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg

              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da

#endif
           end do
        end do
     else 
        ! Make water happen.
        do x = 1,nx
           do z = 1,nz
              N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*pg

              if ( pd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da             
#endif
           end do
        end do
     end if
  end do
end subroutine lbe_init_bi_lam_y


!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the z direction.
subroutine lbe_init_bi_lam_z(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk) :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tz	! z-coordinate in global system.
  
  lamwidth = fr1 + fr2

  if (myrankc == 0) then
     if (mod(tnz,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! Note that ty is zero-based.
  do z = 1,nz
     tz = z - 1 + ccoords(3)*nz
     i = mod(tz, lamwidth)

     if ( i < fr1 ) then
        ! Make oil happen.
        do x = 1,nx
           do y = 1,ny
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg

              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da
#endif
           end do
        end do
     else 
	! Make water happen.
        do x = 1,nx
           do y = 1,ny
              N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*pg

              if ( pd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da

#endif
           end do
        end do
     end if
  end do
end subroutine lbe_init_bi_lam_z

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the z direction.
subroutine lbe_init_bi_sin_lam_x(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real(kind=rk)  :: theta, phi, r, F(19) 
  integer :: lamwidth	! Width of one layer.
  integer :: tx, tz	! coordinates in global system.
  real(kind=rk)  :: sinAmplitude, sinPeriods ! parameters of sinusoidal disturbance

  sinAmplitude = 5
  sinPeriods = 1
    
  lamwidth = fr1 + fr2
  
  if (myrankc == 0) then
     if (mod(tnz,lamwidth) .ne. 0) then
        print*,'******* WARNING ******'
        print*,'* Truncated lamellae *'
        print*,'**********************'
     endif
  endif

  if (COLLISIONTYPE_ID .eq. MRT) then
     ! CALL mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  ! Note that ty is zero-based.
  do z = 1,nz 
     tz = z - 1 + ccoords(3)*nz
     i = mod(tz, lamwidth)
     do x = 1,nx
        tx = x - 1 + ccoords(1)*nx
        do y = 1,ny
           if ( i > ( fr1 + sinAmplitude*sin( 1.0d0/sinPeriods*(2.d0*pi/nx)*x ) ) ) then
              ! Make oil happen.
              N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*fg

              if ( fd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da
#endif
           else 
              N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
              N(x,y,z)%n_b(:) = F(:)*pb
#endif

#ifndef NOSURFACTANT
              N(x,y,z)%n_s(:) = F(:)*pg

              if ( pd .eq. 0 ) then
                 CALL random_number(theta)
                 CALL random_number(phi)
                 CALL random_number(r)
                 theta = theta*pi ! 0 <= theta <= pi
                 phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
                 N(x,y,z)%da(1) = r*sin(theta)*cos(phi)
                 N(x,y,z)%da(2) = r*sin(theta)*sin(phi)
                 N(x,y,z)%da(3) = r*cos(theta)
              else
                 N(x,y,z)%da(:) = (/ 0., 1., 0. /)
              end if
              N(x,y,z)%db = N(x,y,z)%da
#endif
           end if
        end do
     end do
  end do
end subroutine lbe_init_bi_sin_lam_x

!> This routine perturbs the system by up to a fraction \c perturbation.
!> It may be used to avoid metastable states.
!>
!> The occupation number on each site is multiplied by a random number
!> from a flat distrubution in the range
!> (1-perturbation)..(1+perturbation).
subroutine lbe_init_perturb(N)
	real(kind=rk), dimension(nvecs) :: wiggle
	integer :: x,y,z
	type(lbe_site),dimension(0:,0:,0:) :: N

	do z=1,nz
	 do y=1,ny
	  do x=1,nx
		call random_number(wiggle)
		wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
		N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*wiggle
#ifndef SINGLEFLUID
		call random_number(wiggle)
		wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
		N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*wiggle
#endif
		call random_number(wiggle)
#ifndef NOSURFACTANT
		wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
		N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)*wiggle
#endif
	  end do
	 end do
	end do
	
end subroutine lbe_init_perturb



end module lb3d_ic_module
