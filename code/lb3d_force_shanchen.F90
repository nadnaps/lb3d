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

!> Contains routines to calculate the Shan-Chen forces 
!> for single component multiphase, binary immiscible 
!> and ternary amphiphilic mixtures
!> \todo This should be more abstract, i.e. take nearest
!> neighbour densities as arguments rather than evaluating 
!> fixed ones.
module lb3d_force_shanchen_module

  use lb3d_global_module, only: rk
  use lb3d_lattice_module
  use lb3d_config_module

  implicit none

  interface lb3d_calc_sc_forces
     module procedure lb3d_calc_sc_force_single
#ifndef SINGLEFLUID
     module procedure lb3d_calc_sc_force_binary
#endif
#ifndef NOSURFACTANT
     module procedure lb3d_calc_sc_force_ternary
#endif

  end interface

contains


subroutine lb3d_calc_sc_force_single(N,x,y,z,f_r)

  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N

  integer,intent(in) :: x, y, z
  integer :: xa, ya, za, s

  real(kind=rk), dimension(1:,1:,1:,1:),intent(inout) :: f_r

  real(kind=rk),  dimension(3) :: f_rr, f_wr

  real(kind=rk) :: psi_local, psi_neighbour

  ! Green's Functions
  real(kind=rk) :: green_r

  ! Force weights
  integer :: g_sc(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)

#ifdef LB3D_DEBUG_INFO    
 ! call log_msg('In lb3d_calc_sc_force_single',.false.) 
#endif
  if (SCMP.eqv..true.) then
     psi_local = 0.d0
     psi_neighbour = 0.d0
     f_rr = 0.d0
     f_wr = 0.d0

     do s=1,nnonrest
        xa = x + cx(s)
        ya = y + cy(s)
        za = z + cz(s)

        if (N(xa,ya,za)%rock_state.ne.0.) then 
           ! Handle wall interaction

           if (ZEROFORCEOFFSET.eqv..true.) then
              ! Interpolate densities to make the wall force-free
              call zero_force_offset(N,x,y,z)
           end if

           psi_neighbour = get_psi( N(xa,ya,za)%rho_r+N(xa,ya,za)%rock_colour_r )

           green_r = psi_neighbour*g_sc(s)

           f_wr = f_wr + green_r*c(s,:)

        else

           psi_neighbour = get_psi( N(xa,ya,za)%rho_r )

           green_r = psi_neighbour*g_sc(s)

           f_rr = f_rr + green_r*c(s,:)

        end if

     end do ! s over nonrest

     ! Self interaction
     f_rr = - g_rr*f_rr

     ! Wall interaction
     f_wr = - g_wr*f_wr

     f_r(x,y,z,:) = (f_rr + (f_wr/tau_r*tau_wr) ) * (get_psi(N(x,y,z)%rho_r))
  end if

end subroutine lb3d_calc_sc_force_single


#ifndef SINGLEFLUID
subroutine lb3d_calc_sc_force_binary(N,x,y,z,f_b,f_r)

  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N

  integer,intent(in) :: x, y, z
  integer :: xa, ya, za, s

  real(kind=rk), dimension(1:,1:,1:,1:),intent(out) :: f_b
  real(kind=rk), dimension(1:,1:,1:,1:),intent(out) :: f_r

  real(kind=rk),  dimension(3) :: f_br, f_rb, f_wr, f_wb

  real(kind=rk) :: psi_b, psi_r

  ! Green's Functions
  real(kind=rk) :: green_b, green_r
  
  ! Force weights
  integer :: g_sc(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)

#ifdef LB3D_DEBUG_INFO    
 ! call log_msg('In lb3d_calc_sc_force_binary',.false.) 
#endif

  psi_b = 0.d0
  psi_r = 0.d0
  f_br = 0.d0
  f_rb = 0.d0
  f_wr = 0.d0
  f_wb = 0.d0
  
  do s=1,nnonrest
     xa = x + cx(s)
     ya = y + cy(s)
     za = z + cz(s)
     
     if (N(xa,ya,za)%rock_state.ne.0.) then 
        ! Handle wall interaction

        if (ZEROFORCEOFFSET) then
           ! Interpolate densities to make the wall force-free
           call zero_force_offset(N,x,y,z)
        end if

        psi_r = get_psi( N(xa,ya,za)%rho_r+N(xa,ya,za)%rock_colour_r )
        psi_b = get_psi( N(xa,ya,za)%rho_b+N(xa,ya,za)%rock_colour_b )

        green_r = psi_b*g_sc(s)
        green_b = psi_r*g_sc(s)
        
        f_wb = f_wb + green_r*c(s,:)
        f_wr = f_wr + green_b*c(s,:)

     else

        psi_r = get_psi( N(xa,ya,za)%rho_r )
        psi_b = get_psi( N(xa,ya,za)%rho_b )

        green_r = psi_b*g_sc(s)
        green_b = psi_r*g_sc(s)
        
        f_rb = f_rb + green_r*c(s,:)
        f_br = f_br + green_b*c(s,:)

     end if

  end do ! s over nonrest
  
  ! Two-component interactions
  f_br = - g_br*f_br
  f_rb = - g_br*f_rb

  ! Wall interactions
  f_wr = - g_wr*f_wr
  f_wb = - g_wb*f_wb

  f_b(x,y,z,:) = (f_br+(f_wb/tau_b*tau_wb)) * (get_psi(N(x,y,z)%rho_b))
  f_r(x,y,z,:) = (f_rb+(f_wr/tau_r*tau_wr)) * (get_psi(N(x,y,z)%rho_r))

end subroutine lb3d_calc_sc_force_binary
#endif


#ifndef NOSURFACTANT
subroutine lb3d_calc_sc_force_ternary(N,x,y,z,f_b,f_r,f_s)

  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N

  integer,intent(in) :: x, y, z
  integer :: xa, ya, za, s

  real(kind=rk), dimension(1:,1:,1:,1:),intent(out) :: f_b
  real(kind=rk), dimension(1:,1:,1:,1:),intent(out) :: f_r
  real(kind=rk), dimension(1:,1:,1:,1:),intent(out) :: f_s

  real(kind=rk), dimension(3) :: da
  real(kind=rk) :: d_dotc, di_dot_d, di_dotc

  real(kind=rk),  dimension(3) :: f_br, f_rb
  real(kind=rk),  dimension(3) :: f_cs, f_sc
  real(kind=rk),  dimension(3) :: f_ss
  real(kind=rk) :: psi_s

  real(kind=rk) :: psi_b, psi_r

  ! Green's Functions
  real(kind=rk) :: green_b, green_r
  real(kind=rk) :: greenc_s, green_ss
  
  ! Independent wall interactions
  real(kind=rk) :: green_wr, green_wb  

  integer :: g_sc(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)

#ifdef LB3D_DEBUG_INFO    
 ! call log_msg('In lb3d_calc_sc_force_ternary',.false.) 
#endif

  f_ss = 0.d0
  psi_s = 0.d0
  psi_b = 0.d0
  psi_r = 0.d0
  f_br = 0.d0
  f_rb = 0.d0
  f_cs = 0.d0
  f_sc = 0.d0
  
  do s=1,nnonrest
     xa = x + cx(s)
     ya = y + cy(s)
     za = z + cz(s)
     
!FIXME do not forget te restvec-case for obstacle-sites
     psi_r = get_psi( N(xa,ya,za)%rho_r )
     psi_b = get_psi( N(xa,ya,za)%rho_b )

     green_r = psi_b*g_sc(s)
     green_b = psi_r*g_sc(s)


  !   psi_s = get_psi( N(xa,ya,za)%rho_s )

     f_rb = f_rb + green_r*c(s,:)
     f_br = f_br + green_b*c(s,:)

     ! Now think about the force on colours due to surfactants.
     ! da = surf direction vector = d^eq(x + c_i)*psi_s(x + c_i)
     
     da = N(xa,ya,za)%db * get_psi( N(xa,ya,za)%rho_s )

!     da = N(xa,ya,za)%db * get_psi( N(xa,ya,za)%rho_s )

     ! Dot with the vector c_i
     ! Remember that d (c c) == (d.c) c
     !               ^  ^ ^      ^ ^  ^

     d_dotc = sum(da*c(s,:))
     f_cs = f_cs + (da - D/2.d0 *d_dotc*c(s,:) )*g_sc(s)

     ! That's the force on colours due to surfs done.
     ! Now think about the force on surfs due to colours.
     ! Dot c_i with the direction of the surfactants.
     
     di_dotc = sum( N(x,y,z)%db*c(s,:) )
     greenc_s = (psi_b - psi_r )*g_sc(s)
     f_sc = f_sc + greenc_s*(N(x,y,z)%db-D/2.*di_dotc*c(s,:))
!     f_sc = f_sc + greenc_s*(N(x,y,z)%db-D/2.*di_dotc*c(s,:))
     ! Now think about the surf-surf forces.

     di_dot_d = sum(da*N(x,y,z)%db)
     green_ss = 2.d0*D*g_sc(s)

     f_ss = f_ss + green_ss*( &
          + (di_dot_d - D/2.d0*di_dotc*d_dotc)*c(s,:)&
          + da*di_dotc + d_dotc*N(x,y,z)%db )

!     f_ss = f_ss + green_ss*( &
!          + (di_dot_d - D/2.d0*di_dotc*d_dotc)*c(s,:)&
!          + da*di_dotc + d_dotc*N(x,y,z)%db )

     

  end do ! s over nonrest
  
  f_br = - g_br*f_br
  f_rb = - g_br*f_rb

  f_sc = 2.d0*g_bs*f_sc
  f_cs = - 2.d0*g_bs*f_cs
  f_ss = - g_ss*f_ss

  f_b(x,y,z,:) = (f_br + f_cs) * (get_psi(N(x,y,z)%rho_b))
  f_r(x,y,z,:) = (f_rb - f_cs) * (get_psi(N(x,y,z)%rho_r))
  f_s(x,y,z,:) = (f_sc + f_ss) * (get_psi(N(x,y,z)%rho_s))

end subroutine lb3d_calc_sc_force_ternary
#endif


  pure function get_psi (rho)
    real(kind=rk),intent(in) :: rho
    real(kind=rk) :: get_psi
#ifdef COMPAT_PSIFUN
    select case (psifunc)
    case (0)
       get_psi = min(1.0_rk, real(rho,rk))
    case (1)
       get_psi = rho
    case (2)
#endif
    get_psi = 1.0_rk - exp(-rho)
#ifdef COMPAT_PSIFUN
    case default
!       call log_msg("ERROR: Unknown psi functional.",.false.)
!       call abend
    end select
#endif
    return
  end function get_psi

  !> Interpolating the densities surrounding a rock-site to make it
  !> Shan-Chen-force free.
  subroutine zero_force_offset(N,x,y,z)
    integer :: x, y, z, xa, ya, za, s
    real(kind=rk) :: f_r_zero_offset, f_b_zero_offset, f_s_zero_offset, nfc
    type(lbe_site),dimension(1-halo_extent:,1-halo_extent:,1-halo_extent:) :: N
    
    if (N(x,y,z)%rock_state.ne.0.d0) then  
       f_r_zero_offset = 0.d0
#ifndef SINGLEFLUID
       f_b_zero_offset = 0.d0
#ifndef NOSURFACTANT
       f_s_zero_offset = 0.d0
#endif
#endif
         nfc = 0.0d0
         if (ZFOSdiag.eqv..true.) then
            do s=1,nnonrest
               xa = x + cx(s)
               ya = y + cy(s)
               za = z + cz(s)
               if (N(xa,ya,za)%rock_state.eq.0.d0) then
                  f_r_zero_offset = f_r_zero_offset + sum(N(xa,ya,za)%n_r(:)*g(:))
#ifndef SINGLEFLUID
                  f_b_zero_offset = f_b_zero_offset + sum(N(xa,ya,za)%n_b(:)*g(:))
!!$#ifndef NOSURFACTANT
!!$                  f_s_zero_offset = f_s_zero_offset + sum(N(xa,ya,za)%n_s(:)*g(:))
!!$#endif
#endif
                  nfc = nfc + 1.0d0
               end if
            end do
         else
            do s=1,6
               xa = x + cx(s)
               ya = y + cy(s)
               za = z + cz(s)
               if (N(xa,ya,za)%rock_state.eq.0.d0) then
                  f_r_zero_offset = f_r_zero_offset + sum(N(xa,ya,za)%n_r(:)*g(:))
#ifndef SINGLEFLUID
                  f_b_zero_offset = f_b_zero_offset + sum(N(xa,ya,za)%n_b(:)*g(:))
!!$#ifndef NOSURFACTANT
!!$                  f_s_zero_offset = f_s_zero_offset + sum(N(xa,ya,za)%n_s(:)*g(:))
!!$#endif
#endif
                  nfc = nfc + 1.0d0
               end if
            end do
         end if
         
         if ( nfc .ne. 0 ) then
            N(x,y,z)%n_r(restvec) = f_r_zero_offset/nfc
#ifndef SINGLEFLUID
            N(x,y,z)%n_b(restvec) = f_b_zero_offset/nfc
!!$#ifndef NOSURFACTANT
!!$            N(x,y,z)%n_s(restvec) = f_s_zero_offset/nfc
!!$#endif
#endif
         end if
      end if
    end subroutine zero_force_offset

end module lb3d_force_shanchen_module
