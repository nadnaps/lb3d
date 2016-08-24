#include "lb3d.h"
!===================================================================
!======
!
! Copyright 1999-2012, Owners retain copyrights to their respective
!  works.
!
! This file is part of lb3d.
!
! lb3d is free software: you can redistribute it and/or modify it
!  under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
!  (at
! your option) any later version.
!
! lb3d is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!  or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
! License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with lb3d. If not, see <http://www.gnu.org/licenses/>.
!
!===================================================================
!======

!> Contains subroutines to calculate derived properties and interpolations
module lb3d_analysis_module
!========================================================================

  use lb3d_config_module
  use lb3d_lattice_module
  use lb3d_mpi_module, only: local_coordinates
  use lb3d_global_module, only: rk

  implicit none

  !> calculates local fluid densities
  interface density
     module procedure density_REAL_NVECS,density_LBE_SITE
  end interface

  !> interpolates scalar or vector fields
  interface interpolate
     module procedure interpolate_r,interpolate_r3
  end interface

  !> calculates local mass flow rates
  interface massflow
     module procedure massflow_REAL_NVECS,massflow_LBE_SITE
  end interface

  !> calculates local velocity vectors
  interface velocity
     module procedure velocity_REAL_NVECS,velocity_LBE_SITE
  end interface

contains

  !> Updates the local densities stored in \c N%rho_x
  !> from the distributions \c N%n_x for a given lattice site 
  subroutine lb3d_update_density (N,x,y,z)
  type(lbe_site),dimension(0:,0:,0:) :: N  
  integer :: x, y, z

             N(x,y,z)%rho_r = sum( N(x,y,z)%n_r(:) * g(:) )
#ifndef SINGLEFLUID
             N(x,y,z)%rho_b = sum( N(x,y,z)%n_b(:) * g(:) )
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%rho_s = sum( N(x,y,z)%n_s(:) * g(:) )
#endif

  end subroutine lb3d_update_density

  !> Updates the local velocities stored in \c N%u_x
  !> from the distributions \c N%n_x for a given lattice site 
  subroutine lb3d_update_velocity (N,x,y,z)
  type(lbe_site),dimension(0:,0:,0:) :: N  
  integer :: x, y, z
  integer :: s

             N(x,y,z)%u_r(:) = 0.0_rk
#ifndef SINGLEFLUID
             N(x,y,z)%u_b(:) = 0.0_rk
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%u_s(:) = 0.0_rk
#endif
             do s=1,nvecs
                N(x,y,z)%u_r(:) = N(x,y,z)%u_r(:) + ( N(x,y,z)%n_r(s) * g(s) * c(s,:) )
#ifndef SINGLEFLUID
                N(x,y,z)%u_b(:) = N(x,y,z)%u_b(:) + ( N(x,y,z)%n_b(s) * g(s) * c(s,:) )
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%u_s(:) = N(x,y,z)%u_s(:) + ( N(x,y,z)%n_s(s) * g(s) * c(s,:) )
#endif
             end do

  end subroutine lb3d_update_velocity
  
  !> Updates the local densities stored in \c N%rho_x
  !> from the distributions \c N%n_x for the whole lattice
  subroutine lb3d_update_density_loop (N)
  type(lbe_site),dimension(0:,0:,0:) :: N  
  integer :: x, y, z

    do z=0,nz+1
       do y=0,ny+1
          do x=0,nx+1
             N(x,y,z)%rho_r = sum( N(x,y,z)%n_r(:) * g(:) )
#ifndef SINGLEFLUID
             N(x,y,z)%rho_b = sum( N(x,y,z)%n_b(:) * g(:) )
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%rho_s = sum( N(x,y,z)%n_s(:) * g(:) )
#endif
          end do
       end do
    end do
  end subroutine lb3d_update_density_loop

  !> Updates the local velocity stored in \c N%u_x
  !> from the distributions \c N%n_x for the whole lattice
  subroutine lb3d_update_velocity_loop (N)
  type(lbe_site),dimension(0:,0:,0:) :: N  
  integer :: x, y, z
  integer :: s
    do z=0,nz+1
       do y=0,ny+1
          do x=0,nx+1
             N(x,y,z)%u_r(:) = 0.0_rk
#ifndef SINGLEFLUID
             N(x,y,z)%u_b(:) = 0.0_rk
#endif
#ifndef NOSURFACTANT
             N(x,y,z)%u_s(:) = 0.0_rk
#endif
             do s=1,nvecs
                N(x,y,z)%u_r(:) = N(x,y,z)%u_r(:) + ( N(x,y,z)%n_r(s) * g(s) * c(s,:) )
#ifndef SINGLEFLUID
                N(x,y,z)%u_b(:) = N(x,y,z)%u_b(:) + ( N(x,y,z)%n_b(s) * g(s) * c(s,:) )
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%u_s(:) = N(x,y,z)%u_s(:) + ( N(x,y,z)%n_s(s) * g(s) * c(s,:) )
#endif
             end do
          end do
       end do
    end do
  end subroutine lb3d_update_velocity_loop


!---Migrated from helper

!> Returns the cross (vector) product of two three-dimensional vectors
!>
!> \param[in] a vector \c a
!>
!> \param[in] b vector \c b
!>
!> \returns \f[\mathbf{a}\times\mathbf{b}\f]
pure function cross_product(a,b)
  real(kind=rk),dimension(3) :: cross_product
  real(kind=rk),intent(in) :: a(3),b(3)

  cross_product(1) = a(2)*b(3)-a(3)*b(2)
  cross_product(2) = a(3)*b(1)-a(1)*b(3)
  cross_product(3) = a(1)*b(2)-a(2)*b(1)
end function cross_product

!> trilinear interpolation of a discrete scalar field
!>
!> \param[in] Nr floating point scalar field (position indices are
!> local as in N(x,y,z); there is a halo of depth 1)
!>
!> \param[in] x global position. Must be at some place that is covered
!> by the local \c N, however, this subroutine is smart enough to deal
!> with pbc.
!>
!> \param[out] r interpolated value
subroutine interpolate_r(Nr,x,r)
  real(kind=rk),intent(in) :: Nr(0:,0:,0:)
  real(kind=rk),intent(in) :: x(3)
  real(kind=rk),intent(out) :: r
  integer l
  integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
  real(kind=rk) :: xx(3),weight

  CALL local_coordinates(x(:),xx(:))
  lx(:) = floor(xx(:))

  r = 0.0_rk
  lattice_points: do l=1,n_lp
    ! weight eac lattice point with the volume of the hypercube
    ! spanned by the particle position and the opposite lattice
    ! point
    lp(:) = lx(:) + lp_sur(:,l)
    opp_lp(:) = lx(:) + opp_lp_sur(:,l)
    weight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))
    r = r + weight*Nr(lp(1),lp(2),lp(3))
  end do lattice_points
end subroutine interpolate_r

!> trilinear interpolation of a discrete 3-vector field
!>
!> \param[in] Nr3 floating point 3-vector field (position indices are
!> local as in N(x,y,z); there is a halo of depth 1)
!>
!> \param[in] x global position. Must be at some place that is covered
!> by the local \c N, however, this subroutine is smart enough to deal
!> with pbc.
!>
!> \param[out] r3 interpolated value
subroutine interpolate_r3(Nr3,x,r3)
  real(kind=rk),intent(in) :: Nr3(0:,0:,0:,1:)
  real(kind=rk),intent(in) :: x(3)
  real(kind=rk),intent(out) :: r3(3)
  integer l
  integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
  real(kind=rk) :: xx(3),weight

  CALL local_coordinates(x(:),xx(:))
  lx(:) = floor(xx(:))

  r3(:) = 0.0_rk
  lattice_points: do l=1,n_lp
    ! weight each lattice point with the volume of the hypercube
    ! spanned by the particle position and the opposite lattice
    ! point
    lp(:) = lx(:) + lp_sur(:,l)
    opp_lp(:) = lx(:) + opp_lp_sur(:,l)
    weight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))
    r3 = r3 + weight*Nr3(lp(1),lp(2),lp(3),:)
  end do lattice_points
end subroutine interpolate_r3



!> Returns the length of a vector
!>
!> \param[in] v vector (3-dimensional)
!>
!> \returns \f[|\mathbf{v}|\f]
pure function norm(v)
  real(kind=rk) :: norm
  real(kind=rk),intent(in) :: v(3)
  norm = sqrt(dot_product(v,v))
end function norm

!> Normalizes a vector to unit length
!>
!> \param[in] v vector (3-dimensional)
!>
!> \returns \f[\mathbf{v}/|\mathbf{v}|\f]
pure function unit_vector(v)
  real(kind=rk),dimension(3) :: unit_vector
  real(kind=rk),intent(in) :: v(3)
  unit_vector = v/sqrt(dot_product(v,v))
end function unit_vector

!> Calculates the density for a given set of distributions.
!>
!> \param[in] n local population vector
!>
!> \returns density (not scaled with \c amass_{rbs})
function density_REAL_NVECS(n)
  real(kind=rk) :: density_REAL_NVECS
  real(kind=rk),intent(in) :: n(nvecs)

  density_REAL_NVECS = 2.0_rk*sum(n(1:6))+sum(n(7:19))
end function density_REAL_NVECS

!> Calculates the total density for a given lattice site.
!>
!> \param[in] s lattice site
!>
!> \returns density (not scaled with \c amass_{rbs})
function density_LBE_SITE(s)
  real(kind=rk) :: density_LBE_SITE
  type(lbe_site),intent(in) :: s

  density_LBE_SITE = &
#ifndef SINGLEFLUID
       &density(s%n_b)+&
#ifndef NOSURFACTANT
       &density(s%n_s)+&
#endif
#endif
       &density(s%n_r)
end function density_LBE_SITE

!> Calculates the mass flow (that is the momentum) for a given set of
!> distributions, a possible rescaling of masses using \c amass_[rbs]
!> is not taken into account.
!>
!> \param[in] n local population vector
!>
!> \returns mass flow vector
function massflow_REAL_NVECS(n)
  real(kind=rk) :: massflow_REAL_NVECS(3)
  real(kind=rk),intent(in) :: n(nvecs)

  massflow_REAL_NVECS(1) = &
       &2.0_rk*(n(1)-n(2))+n(7)+n(8)+n(9)+n(10)-n(11)-n(12)-n(13)-n(14)
  massflow_REAL_NVECS(2) = &
       &2.0_rk*(n(3)-n(4))+n(7)+n(11)+n(15)+n(16)-n(8)-n(12)-n(17)-n(18)
  massflow_REAL_NVECS(3) = &
       &2.0_rk*(n(5)-n(6))+n(9)+n(13)+n(15)+n(17)-n(10)-n(14)-n(16)-n(18)
end function massflow_REAL_NVECS

!> Calculates the total mass flow (that is the momentum) for a given
!> lattice site taking into account \c amass_[rbs].
!>
!> \param[in] s lattice site
!>
!> \returns mass flow vector
function massflow_LBE_SITE(s)
  real(kind=rk) :: massflow_LBE_SITE(3)
  type(lbe_site),intent(in) :: s

  massflow_LBE_SITE = &
#ifndef SINGLEFLUID
       &amass_b*massflow(s%n_b)+&
#ifndef NOSURFACTANT
       &amass_s*massflow(s%n_s)+&
#endif
#endif
       &amass_r*massflow(s%n_r)
end function massflow_LBE_SITE

!> Calculates the velocity for a given lattice site.
!>
!> \param[in] s lattice site
!>
!> \returns velocity vector
function velocity_LBE_SITE(s)
  real(kind=rk) :: velocity_LBE_SITE(3)
  type(lbe_site),intent(in) :: s

  velocity_LBE_SITE = massflow_LBE_SITE(s)/density_LBE_SITE(s)
end function velocity_LBE_SITE

!> Calculates the velocity for a given set of distributions.
!>
!> \param[in] n local population vector
!>
!> \returns velocity vector
function velocity_REAL_NVECS(n)
  real(kind=rk) :: velocity_REAL_NVECS(3)
  real(kind=rk),intent(in) :: n(nvecs)

  velocity_REAL_NVECS = massflow_REAL_NVECS(n)/density_REAL_NVECS(n)
end function velocity_REAL_NVECS

end module lb3d_analysis_module

