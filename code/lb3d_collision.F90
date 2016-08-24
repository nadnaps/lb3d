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

!> Functions needed for a separated collision-portion of the 
!> lattice Boltzmann algorithm
module lb3d_collision_module
!========================================================================
  use lb3d_global_module, only: rk
  use lb3d_config_module
  use lb3d_bc_leesedwards_module
  use lb3d_bc_periodic_module
  use lb3d_timer_module
  use lb3d_force_shanchen_module
  use lb3d_lattice_module
  use lb3d_bdist_module
  use lb3d_analysis_module

  implicit none

contains

  subroutine lb3d_collision
#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lb3d_collision',.false.) 
#endif
    call lb3d_collision_step(N)
    
  end subroutine lb3d_collision
     

  !> Takes a subdomain array \c N and collides and redistributes 
  !> the particles at each site.
  !> If this is ever run in a situation where a Fortran 95 compiler
  !> is available, it could very possibly benefit from being rephrased
  !> in terms of the \c forall statement in F95.
  !>
  !> Contains the loops for the force calculation of Shan-Chen and Body-forces
  !> applied via an update of the equilibrium velocity as well.
  !> \todo The force calculation calls should be in \c lb3d_time_loop
  subroutine lb3d_collision_step(N)

    type(lbe_site),dimension(0:,0:,0:) :: N

    integer :: x, y, z, s
    integer :: tz, tx, ty
    integer :: xa, ya, za, ss
    logical :: minz, maxz
    integer :: xb, yb, zb
    integer :: xaa, yaa, zaa
    integer :: xbb, ybb, zbb

    ! Nett forces
    real*8, dimension(nx,ny,nz,3) :: f_r
#ifndef SINGLEFLUID
    real*8, dimension(nx,ny,nz,3) :: f_b
#endif
#ifndef NOSURFACTANT
    real*8, dimension(nx,ny,nz,3) :: f_s
#endif

    ! Velocities
    real*8, dimension(3) :: p_r, p_tilde, u_tilde
#ifndef SINGLEFLUID
    real*8, dimension(3) :: p_b
#endif
#ifndef NOSURFACTANT
    real*8, dimension(3) :: p_s
    real*8 :: sN_s, rho_s
    real*8, dimension(nvecs) :: boltz_s
    real*8, dimension(3) :: sN_eq(nvecs)
    real*8, dimension(3) :: du_s
    real(kind=rk), dimension(3) :: dd_s
#endif
    ! Number densities
    real*8 :: bN_s, rN_s
    
    ! Weighted masses
    real*8 :: rho_r, rho_tilde
#ifndef SINGLEFLUID
    real*8 :: rho_b
#endif

    
    ! Boltzmann distributions
    real*8, dimension(nvecs) :: boltz_b, boltz_r
    
    ! Equilibrium distributions at a site
    real*8, dimension(3) :: bN_eq(nvecs), rN_eq(nvecs)
    real*8, dimension(3) :: du_r, du_b
    real(kind=rk), dimension(3) :: dd_r, dd_b
    
    ! temporary variable for force calculation
    logical :: applyforce
    
    call start_timer(ti_intf)


#ifdef FORCE_SHANCHEN

#ifndef NOSURFACTANT
    ! Loop over all lattice sites (icky)

    do z=1,nz
       do y=1,ny
          do x=1,nx
             call lb3d_dipole_collision (N,x,y,z)
          end do
       end do
    end do
   
    if ( ( boundary_cond .ne. 0 ) .and. &
         ( inv_fluid .eq. 5 .or. inv_fluid .eq. 6 ) ) then

       call le_halo_exchange(N)
       call lb3d_update_density_loop (N)
       call lb3d_update_velocity_loop (N)
    else

       call lbe_halo_exchange(whole_N)

    endif
#endif ifndef NOSURFACTANT

#endif ifdef FORCE_SHANCHEN
    ! The equilibrium values and Hamiltonians have now been calculated.
    ! Now redistribute the particles.
    ! Loop over all lattice sites (yick).
    
    ! If we have only a single fluid, we do not need to calculate forces.
    ! #ifndef SINGLEFLUID

    do z=1,nz
       do y=1,ny
          do x=1,nx
             f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
             f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
             f_s(x,y,z,:) = 0.d0
#endif
#ifdef FORCE_SHANCHEN
             ! only calculate forces on fluid at non-rock sites
             if (N(x,y,z)%rock_state == 0.d0) then
#ifdef SINGLEFLUID
                call lb3d_calc_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
                call lb3d_calc_sc_forces(N,x,y,z,f_b,f_r)
#else
                call lb3d_calc_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
             endif
#endif ifdef FORCE_SHANCHEN
          end do
       end do
    end do

    call stop_timer(ti_intf)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Collision stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Loop over all lattice sites again (yawn..)
    call start_timer(ti_intc)
    do z=1,nz
       do y=1,ny
          do x=1,nx


#ifdef RELTIME
             ! If flag RELTIME is defined, read the tau for the local site
             ! from an array (from file read at init)
             call lb3d_get_local_tau (N,x,y,z)
#endif

 
#ifdef LBEFORCE
             ! Add a force contribution from the "lbe_force" a per lattice site 
             ! defined force used for coupling
             call lb3d_get_lbe_force (x,y,z,f_r,use_lbe_force)
#ifndef SINGLEFLUID
             call lb3d_get_lbe_force (x,y,z,f_b,use_lbe_force)
#endif
#ifndef NOSURFACTANT
             call lb3d_get_lbe_force (x,y,z,f_s,use_lbe_force)
#endif
#endif 

#ifdef FORCE_BODY
             ! Evaluate BC and input values to find out if body force is to
             ! be applied at the current site
             call lb3d_get_applyforce (x,y,z,applyforce)
             
             ! Add the force-contribution of acceleration g_accn
             ! if applyforce .eqv. .true
             call lb3d_get_body_force (x,y,z,amass_r*N(x,y,z)%rho_r,f_r,applyforce)
#ifndef SINGLEFLUID
             call lb3d_get_body_force (x,y,z,amass_b*N(x,y,z)%rho_b,f_b,applyforce)
#endif
#ifndef NOSURFACTANT
             call lb3d_get_body_force (x,y,z,amass_s*N(x,y,z)%rho_s,f_s,applyforce)
#endif
#ifdef FORCE_GUO             
             ! dd_s are neccessary for the force implementation
             ! according to Phys. Rev. E 65,(2002): Guo
             ! if applyforce .eqv. .true
             call lb3d_get_guo_force (tau_r,dd_r,applyforce)                
#ifndef SINGLEFLUID
             call lb3d_get_guo_force (tau_b,dd_b,applyforce)    
#endif
#ifndef NOSURFACTANT            
             call lb3d_get_guo_force (tau_s,dd_s,applyforce)                
#endif

#endif FORCE_GUO
#endif FORCE_BODY


             ! Calculate the combined equilibrium velocity of all components
             ! e.g. Shan Chen PRE 47
             call lb3d_get_sc_u_eq (N,x,y,z,u_tilde)

             ! Calculate the difference in equilibrium velocity for the individual
             ! components from the pseudo-potential forces, SC PRE 47
             ! FIXME This argument forrest is sub-optimal
             call lb3d_get_sc_du (x,y,z,tau_r,N(x,y,z)%rho_r,amass_r,f_r(x,y,z,:),du_r)
#ifndef SINGLEFLUID
             call lb3d_get_sc_du (x,y,z,tau_b,N(x,y,z)%rho_b,amass_b,f_b(x,y,z,:),du_b)
#endif
#ifndef NOSURFACTANT 
             call lb3d_get_sc_du (x,y,z,tau_s,N(x,y,z)%rho_s,amass_s,f_s(x,y,z,:),du_s)
#endif


                call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3), &
                     du_r(1), du_r(2), du_r(3), &
                     dd_r(1), dd_r(2), dd_r(3), boltz_r)
                rN_eq = N(x,y,z)%rho_r*boltz_r
                
                N(x,y,z)%n_r = N(x,y,z)%n_r - omega_r*( N(x,y,z)%n_r - rN_eq )
  

#ifndef SINGLEFLUID

                call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3), &
                     du_b(1), du_b(2), du_b(3), &
                     dd_b(1), dd_b(2), dd_b(3), boltz_b)
                bN_eq = N(x,y,z)%rho_b*boltz_b
                
                N(x,y,z)%n_b = N(x,y,z)%n_b - omega_b*( N(x,y,z)%n_b - bN_eq )
  

#endif                        
#ifndef NOSURFACTANT

                call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3), &
                     du_s(1), du_s(2), du_s(3), &
                     dd_s(1), dd_s(2), dd_s(3), boltz_s)
                sN_eq = N(x,y,z)%rho_s*boltz_s
                
                N(x,y,z)%n_s = N(x,y,z)%n_s - omega_s*( N(x,y,z)%n_s - sN_eq )
   
#endif
     end do
  end do
end do
call stop_timer(ti_intc)
end subroutine lb3d_collision_step

#ifndef NOSURFACTANT
!> Calculates the post-collision value of the dipole orientation 
!> for a given lattice site
subroutine lb3d_dipole_collision (N,x,y,z)
  implicit none 
  type(lbe_site),dimension(0:,0:,0:) :: N  

  integer :: x, y, z
  integer :: s, ss

  integer :: xa, ya, za  

  real(kind=rk), dimension(3) :: h_c
  real(kind=rk) :: hc_m ! Magnitude of h_c
  ! h_s is a vector involving the surfactant part of the Hamiltonian.
  real(kind=rk), dimension(3) :: h_s
  
  real(kind=rk) :: hs_dot
  
  ! h is the sum of h_s and h_c
  real(kind=rk), dimension(3) :: h
  real(kind=rk) :: h_m  ! Magnitude of h, times beta
  
  ! sd_eq is the equilibrium surfactant density.
  real(kind=rk), dimension(3) :: sd_eq
  
  ! These three are used in some more of Hudong's code
  ! which I don't quite follow.
  real(kind=rk), dimension(3) :: h_bar
  real(kind=rk) :: sd_mag
  real(kind=rk) :: exh

  ! sN is the surfactant density at a site.
  real(kind=rk) :: sN,sN_tmp
  real(kind=rk), dimension(3) :: sd_tmp

  real(kind=rk) :: psi_b, psi_r

  real(kind=rk), dimension(3) :: pas = (/ 0.d0, 0.d0, 0.d0 /)

  integer :: g_sc(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)

  h_c = 0.d0
  h_s = 0.d0
  sN = 0.d0

  h_s = N(x,y,z)%da

  sN_tmp = max(N(x,y,z)%rho_s,dble(10.e-5))

  ! sd_tmp is post-advect dipole
  sd_tmp = h_s/sN_tmp

  ! Now loop over only the nonrest vectors.

  do s=1,nnonrest
     psi_b = 0.d0
     psi_r = 0.d0
     pas = 0.d0

     xa = x + cx(s)
     ya = y + cy(s)
     za = z + cz(s)

     psi_r = psi_r + N(xa,ya,za)%rho_r
     psi_b = psi_b + N(xa,ya,za)%rho_b
     pas = N(xa,ya,za)%da

     ! psi_b, psi_r give the blue, red densities at (xa,ya,za).
     ! Calculate magnitude of h_c                 

     !FIXME why are psi here equal rho? forgotten when introducing psifunc 2?
     hc_m = - sign(dble(1.0),g_bs)*(psi_b - psi_r)*g_sc(s)

     ! Make it into a vector.
     h_c = h_c + hc_m * c(s,:)

     !hs_dot = D/2.*( cx(s)*pas(1) + cy(s)*pas(2) + cz(s)*pas(3) )
     hs_dot = D/2.d0*sum(c(s,:)*pas)

     ! Make h_s into a vector
     h_s = h_s + (pas - c(s,:)*hs_dot) * g_sc(s)
  end do ! s over non-rest vectors

  h = h_c + h_s

  ! h_m = beta * magnitude of h
  h_m = beta*sqrt(sum(h*h))

  sd_eq = 1.0_rk/3.0_rk*beta*h

  ! This is a case distinction. d_eq is only d0*(1/3*beta*h)*h_hat in the limit of
  ! Small beta*h (why is the lim 10^-3?), below evaluating the full integral:
  !  d_eq = d0*(coth(beta h) - 1/(beta*h))*h_hat
  if (h_m > 10.d-3) then
     exh = 1.d0/tanh(h_m)
     h_bar = beta*h/h_m
     sd_mag = exh - 1.d0/h_m
     sd_eq = h_bar*sd_mag
  endif

  ! This line was added in a hurry to allow d_0 to be varied
  sd_eq = sd_eq*d_0

  N(x,y,z)%db = omega_d*sd_eq + ( 1.d0 - omega_d )*sd_tmp
!  N(x,y,z)%db = N(x,y,z)%d
end subroutine lb3d_dipole_collision
#endif


#ifdef RELTIME
!> Reads the local relaxation time for a given lattice site 
!> from \c N. Requires compilation with -DRELTIME
subroutine lb3d_get_local_tau (N,x,y,z)
  type(lbe_site),dimension(0:,0:,0:) :: N  
  integer :: x, y, z
  

  tau_r = N(x,y,z)%taupos_r
  omega_r = 1.0d0/N(x,y,z)%taupos_r
#ifndef SINGLEFLUID
  tau_b = N(x,y,z)%taupos_b
  omega_b = 1.0d0/N(x,y,z)%taupos_b
#endif
#ifndef NOSURFACTANT
  tau_s = N(x,y,z)%taupos_s
  omega_s = 1.0d0/N(x,y,z)%taupos_s
#endif
  
end subroutine lb3d_get_local_tau
#endif

!> Evaluates whether a body-force is to be applied at a 
!> given lattice site, returns logical applyforce
subroutine lb3d_get_applyforce (x,y,z,applyforce)
  integer,intent(in) :: x, y, z
  logical,intent(out) :: applyforce

  integer :: tx, ty, tz

  ! Positions on the global lattice
  tz = z + ccoords(3)*nz
  tx = x + ccoords(1)*nx
  ty = y + ccoords(2)*ny

  applyforce = .false.

  if (inv_fluid.eq.9) then 
     if ( (tz.eq.1) .or. (tz.eq.tnz) ) then
        applyforce = .true.
     end if
  elseif ( (inv_fluid.eq.20) .or. (inv_fluid.eq.21) ) then
     if (tz.eq.1) then     
        applyforce = .true.      
     end if
  elseif ( &
       (tz.ge.g_accn_min) .and. (tz.le.g_accn_max) .and. &
       (tx.ge.g_accn_min_x) .and. (tx.le.g_accn_max_x) .and. &
       (ty.ge.g_accn_min_y) .and. (ty.le.g_accn_max_y) &
       ) then
     applyforce = .true.
  endif

end subroutine lb3d_get_applyforce


!> Calculates a force from the acceleration input parameter and 
!> adds it to a given force-array
subroutine lb3d_get_body_force (x,y,z,mass,force,applyforce)

!===============================================================================
  integer,intent(in) :: x, y, z
  real(kind=rk),intent(in)  :: mass
  real(kind=rk), dimension(1:,1:,1:,1:),intent(inout) :: force

  logical,intent(in) :: applyforce  

  
  if ( applyforce .eqv. .true. ) then

     force(x,y,z,3) = force(x,y,z,3) + mass * g_accn
     force(x,y,z,2) = force(x,y,z,2) + mass * g_accn_y
     force(x,y,z,1) = force(x,y,z,1) + mass * g_accn_x

  endif
  !===============================================================================

end subroutine lb3d_get_body_force

!> Read a force from an external (global) force-data array 
!> and adds it to a given force-array.
!> This is an entrypoint for coupled methods.
subroutine lb3d_get_lbe_force (x,y,z,force,use_lbe_force)

  integer,intent(in) :: x, y, z
  logical,intent(in) :: use_lbe_force

  real(kind=rk), dimension(1:,1:,1:,1:),intent(inout) :: force

  if (use_lbe_force .eqv. .true.)then  
     force(x,y,z,:) = force(x,y,z,:) + lbe_force(:,1,x,y,z)
     lbe_force(:,1,x,y,z) = 0.0_rk
  end if
end subroutine lb3d_get_lbe_force

!> Calculates terms of a Guo-style force implementation
!>  that are additional to Shan-Chen-style forces.
subroutine lb3d_get_guo_force (tau,dd,applyforce)

  logical,intent(in) :: applyforce
  real(kind=rk),intent(in) :: tau

  real(kind=rk),dimension(3),intent(out) :: dd

  if (applyforce.eqv..true.) then

     dd(3) = (tau - 0.5d0)*g_accn
     dd(2) = (tau - 0.5d0)*g_accn_y
     dd(1) = (tau - 0.5d0)*g_accn_x

  end if

end subroutine lb3d_get_guo_force

!> Calculates the velocity-shift to the equilibrium distribution
!> imposed through a Shan-Chen-style force
subroutine lb3d_get_sc_du (x,y,z,tau,n,amass,force,du)

  integer,intent(in) :: x, y, z
  real(kind=rk),intent(in)  :: amass, tau, n

  real(kind=rk), dimension(3),intent(in) :: force

  real(kind=rk),dimension(3),intent(out) :: du

  real(kind=rk) :: mass, pon

#ifdef COMPAT_PSIFUN                
  if (psifunc .eq. 0) then
     ! This expression is only for the original
     ! funny clipped version of the code
     ! Evaluates to 1 unless densities go over 1

     pon = min(dble(1.)/max(n,dble(10.e-9)), dble(1.))

  else
#endif

     ! sane code that does what the literature says
     pon = 1.0_rk / max(n,real(10.e-9,rk))

#ifdef COMPAT_PSIFUN             
  end if
#endif
 
  du = tau * pon * force / amass



end subroutine lb3d_get_sc_du

!> Calculates the momentum coupled, common equilibrium velocity 
!> for a Shan-Chen-style multi-component system
subroutine lb3d_get_sc_u_eq (N,x,y,z,u_tilde)
  type(lbe_site),intent(in),dimension(0:,0:,0:) :: N  
  integer,intent(in) :: x, y, z
  real(kind=rk),intent(out),dimension(3) :: u_tilde

  real(kind=rk),dimension(3) :: p_tilde
  real(kind=rk) :: rho_tilde

  p_tilde = amass_r*omega_r* N(x,y,z)%u_r(:) 
#ifndef SINGLEFLUID
  p_tilde = p_tilde + amass_b*omega_b* N(x,y,z)%u_b(:)
#endif
#ifndef NOSURFACTANT
  p_tilde = p_tilde + amass_s*omega_s* N(x,y,z)%u_s(:) 
#endif 
                
  ! rho_tilde = sum of masses over relaxation times.
  rho_tilde = amass_r*N(x,y,z)%rho_r*omega_r 
#ifndef SINGLEFLUID
  rho_tilde = rho_tilde + amass_b*N(x,y,z)%rho_b*omega_b 
#endif
#ifndef NOSURFACTANT
  rho_tilde = rho_tilde + amass_s*N(x,y,z)%rho_s*omega_s 
#endif
  
  ! Calculate averaged velocity
  u_tilde = p_tilde/max(rho_tilde,dble(10.e-9))
end subroutine lb3d_get_sc_u_eq

end module lb3d_collision_module

