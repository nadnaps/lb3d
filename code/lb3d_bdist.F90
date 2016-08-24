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

!>This contains the function to calculate the equilibrium distribution.
!>
!>The bdist parameter in the input file chooses the expression that is
!>used.
!FIXME This needs clean-up
module lb3d_bdist_module

  use lb3d_mpi_module
  use lb3d_log_module
  use lb3d_config_module!, only: bdist
  implicit none

    real(kind=rk), dimension(19), save :: S_hat_r
#ifndef SINGLEFLUID  
    real(kind=rk), dimension(19), save :: S_hat_b
#ifndef NOSURFACTANT
    real(kind=rk), dimension(19), save :: S_hat_s
#endif
#endif

contains

#ifndef COMPAT_BDIST
    !> Equilibrium distribution function
    !>
    !> Takes a velocity specified by \c u, an increment \c du, and
    !> and array \c boltz of 19 elements; returns with \c boltz filled in with
    !> the equilibrium distribution of velocities.
    !>
    !> Uses either the Chen operator or the Luo operator, depending on the
    !> value of the input variable \c bdist.
    !>
    !> \todo vectorize this please.
  subroutine boltz_dist(ux, uy, uz, dux, duy, duz, ddx, ddy, ddz, boltz)
    real(kind=rk) :: ux, uy, uz
    real(kind=rk) :: dux, duy, duz, ddx, ddy, ddz
    real(kind=rk),intent(out) :: boltz(nvecs)

    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter,dimension(nvecs) :: weight = (/ &
!!$         1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,& 
!!$         1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,&
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/3.0_rk /)
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T

    integer :: s, ifactor, nni

    real(kind=rk) :: const,const1,const2,even,odd,tux,tuy,tuz,tu1x,tu1y,tu1z
    real(kind=rk) :: tuxt,tuyt,tuzt,uke,u1ke,dke
    real(kind=rk) :: cdotu,cdot1u,cdotd,cdotdu_s,udotdu_s
    real(kind=rk) :: ux_f,uy_f,uz_f,uke_eq,u_mag,delta
    real(kind=rk) :: expux,expuy,expuz,expek,a1eq
    real(kind=rk) :: n_iter
    real(kind=rk) :: rho1_tt,rho_tt,rho_eq
    real(kind=rk) :: px_tt,py_tt,pz_tt
    real(kind=rk) :: delt_ux,delt_uy,delt_uz
    real(kind=rk) :: ux_eq,uy_eq,uz_eq
    real(kind=rk) :: uholdit, choldit, cdholdit, ch2
    real(kind=rk) :: weight1by2,weight2by2,weight1by6,weight2by6

       tux = ux+dux
       tuy = uy+duy
       tuz = uz+duz

       tuxt = tux*holdit
       tuyt = tuy*holdit
       tuzt = tuz*holdit

       uke = 0.5_rk*(tux*tuxt + tuy*tuyt + tuz*tuzt)

       const = 1.0_rk-uke
       
       const1 = weight(1)*const
       const2 = weight(7)*const
  
       weight1by2 = weight(1)/2.0_rk
       weight2by2 = weight(7)/2.0_rk

       weight1by6 = weight(1)/6.0_rk
       weight2by6 = weight(7)/6.0_rk

       cdotu = tuxt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(1) = even+odd
       boltz(2) = even-odd
       cdotu = tuyt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(3) = even+odd
       boltz(4) = even-odd
       cdotu = tuzt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(5) = even+odd
       boltz(6) = even-odd
       cdotu = tuxt+tuyt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(7) = even+odd
       boltz(12) = even-odd
       cdotu = tuxt-tuyt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(8) = even+odd
       boltz(11) = even-odd
       cdotu = tuxt+tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(9) = even+odd
       boltz(14) = even-odd
       cdotu = tuxt-tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(10) = even+odd
       boltz(13) = even-odd
       cdotu = tuyt+tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(15) = even+odd
       boltz(18) = even-odd
       cdotu = tuyt-tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(16) = even+odd
       boltz(17) = even-odd

       boltz(restvec) = weight(restvec)*const

  end subroutine boltz_dist

#else
!ifdef COMPAT_BDIST

    !> Equilibrium distribution function
    !>
    !> Takes a velocity specified by \c u, an increment \c du, and
    !> and array \c boltz of 19 elements; returns with \c boltz filled in with
    !> the equilibrium distribution of velocities.
    !>
    !> Uses either the Chen operator or the Luo operator, depending on the
    !> value of the input variable \c bdist.
    !>
    !> \todo vectorize this please.
  subroutine boltz_dist(ux, uy, uz, dux, duy, duz, ddx, ddy, ddz, boltz)
    real(kind=rk) :: ux, uy, uz
    real(kind=rk) :: dux, duy, duz, ddx, ddy, ddz
    real(kind=rk),intent(out) :: boltz(nvecs)

    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter,dimension(nvecs) :: weight = (/ &
!!$         1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,& 
!!$         1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,&
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         1.0_rk/3.0_rk /)
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T

    integer :: s, ifactor, nni

    real(kind=rk) :: const,const1,const2,even,odd,tux,tuy,tuz,tu1x,tu1y,tu1z
    real(kind=rk) :: tuxt,tuyt,tuzt,uke,u1ke,dke
    real(kind=rk) :: cdotu,cdot1u,cdotd,cdotdu_s,udotdu_s
    real(kind=rk) :: ux_f,uy_f,uz_f,uke_eq,u_mag,delta
    real(kind=rk) :: expux,expuy,expuz,expek,a1eq
    real(kind=rk) :: n_iter
    real(kind=rk) :: rho1_tt,rho_tt,rho_eq
    real(kind=rk) :: px_tt,py_tt,pz_tt
    real(kind=rk) :: delt_ux,delt_uy,delt_uz
    real(kind=rk) :: ux_eq,uy_eq,uz_eq
    real(kind=rk) :: uholdit, choldit, cdholdit, ch2
    real(kind=rk) :: weight1by2,weight2by2,weight1by6,weight2by6


    select case (bdist)

       ! Luo operator
    case (0)
       uke = 0.5d0*(ux*ux + uy*uy + uz*uz)
       udotdu_s = ux*dux + uy*duy + uz*duz
       uholdit = uke*holdit

       do s = 1, nnonrest
          cdotu = cx(s)*ux + cy(s)*uy + cz(s)*uz
          cdotdu_s = cx(s)*dux + cy(s)*duy + cz(s)*duz
          choldit = cdotu*holdit
          cdholdit = cdotdu_s*holdit
          ch2 = choldit*choldit

          !Elena: Replace division by T
          boltz(s) = weight(s)*( 1.d0 + choldit + 0.5d0*ch2 - uholdit &
               + (choldit*ch2)/6.d0 - uholdit*choldit &
               + cdholdit + choldit*cdholdit - udotdu_s*holdit &
               + 0.5d0*ch2*cdholdit - choldit*udotdu_s*holdit &
               - uholdit*cdholdit )
       end do
       boltz(restvec) = weight(restvec)*(1.d0 - uke/T - udotdu_s/T)
       return

       ! Hudong operator.
    case(1)
       ! This is a third order expansion of f(u+delta_u) with no terms missing
       ! as described in Proc. R. Soc. Lond. A (2000) 456, 2043-2057
       ! PRL 81 1618 (1998) and PRE 58 6855 (1998) both claim that this is
       ! inconsistent with a systematic derivation of how force terms should be
       ! introduced. However, it is desirably more stable than the other
       ! offerings here.
       ! Phys. Rev. E 75 036712 use it as well...
       ux_f = ux + dux
       uy_f = uy + duy
       uz_f = uz + duz

       uke_eq = 0.5d0*(ux_f*ux_f + uy_f*uy_f + uz_f*uz_f)

       u_mag = sqrt(2.d0*uke_eq)
       ifactor = 200.d0*u_mag

       n_iter = max(2,ifactor)
       delta = sqrt(1.d0/n_iter)

       ux_eq = ux_f
       uy_eq = uy_f
       uz_eq = uz_f

       rho_eq = 1.d0

       do nni = 1, int(n_iter)
          expux = exp(ux_eq/T)
          expuy = exp(uy_eq/T)
          expuz = exp(uz_eq/T)
          expek = 1.d0/exp(uke_eq/T)

          do s = 1, nnonrest
             a1eq = (expux**cx(s))*(expuy**cy(s))*(expuz**cz(s))
             boltz(s) = rho_eq*weight(s)*a1eq*expek
          end do
          boltz(restvec) = rho_eq*weight(restvec)*expek

          rho1_tt = 0.d0
          px_tt = 0.d0
          py_tt = 0.d0
          pz_tt = 0.d0

          do s = 1, nnonrest
             rho1_tt = rho1_tt + boltz(s)
             px_tt = px_tt + boltz(s)*cx(s)
             py_tt = py_tt + boltz(s)*cy(s)
             pz_tt = pz_tt + boltz(s)*cz(s)
          end do

          rho_tt = rho1_tt + boltz(nvecs)

          delt_ux = ux_f - px_tt
          delt_uy = uy_f - py_tt
          delt_uz = uz_f - pz_tt

          rho_eq = rho_eq/(1.d0 + delta*(rho_tt - 1.d0)/rho_eq)
          ux_eq = ux_eq + delta*delt_ux
          uy_eq = uy_eq + delta*delt_uy
          uz_eq = uz_eq + delta*delt_uz
          uke_eq = 0.5d0*(ux_eq*ux_eq + uy_eq*uy_eq + uz_eq*uz_eq)
       end do
       return

    case (2)

       tux = ux+dux
       tuy = uy+duy
       tuz = uz+duz

       tuxt = tux*holdit
       tuyt = tuy*holdit
       tuzt = tuz*holdit

       uke = 0.5_rk*(tux*tuxt + tuy*tuyt + tuz*tuzt)

       const = 1.0_rk-uke
       
       const1 = weight(1)*const
       const2 = weight(7)*const
  
       weight1by2 = weight(1)/2.0_rk
       weight2by2 = weight(7)/2.0_rk

       weight1by6 = weight(1)/6.0_rk
       weight2by6 = weight(7)/6.0_rk

       cdotu = tuxt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(1) = even+odd
       boltz(2) = even-odd
       cdotu = tuyt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(3) = even+odd
       boltz(4) = even-odd
       cdotu = tuzt
       even = const1+weight1by2*cdotu**2
       odd = const1*cdotu+weight1by6*cdotu**3
       boltz(5) = even+odd
       boltz(6) = even-odd
       cdotu = tuxt+tuyt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(7) = even+odd
       boltz(12) = even-odd
       cdotu = tuxt-tuyt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(8) = even+odd
       boltz(11) = even-odd
       cdotu = tuxt+tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(9) = even+odd
       boltz(14) = even-odd
       cdotu = tuxt-tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(10) = even+odd
       boltz(13) = even-odd
       cdotu = tuyt+tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(15) = even+odd
       boltz(18) = even-odd
       cdotu = tuyt-tuzt
       even = const2+weight2by2*cdotu**2
       odd = const2*cdotu+weight2by6*cdotu**3
       boltz(16) = even+odd
       boltz(17) = even-odd

       boltz(restvec) = weight(restvec)*const

    case (3)
       ! This is a second order expansion of f(u+delta_u) with no terms missing
       ! PRL 81 1618 (1998) and PRE 58 6855 (1998) both claim that this is
       ! inconsistent with a systematic derivation of how force terms should be
       ! introduced.
       uke = 0.5d0*(ux*ux + uy*uy + uz*uz)
       udotdu_s = ux*dux + uy*duy + uz*duz
       
       do s = 1, nnonrest
          cdotu = cx(s)*ux + cy(s)*uy + cz(s)*uz
          cdotdu_s = cx(s)*dux + cy(s)*duy + cz(s)*duz
          
          boltz(s) = weight(s)*(1.d0 + cdotu/T &
               + 0.5d0*(cdotu/T)**2 - uke/T &
               + cdotdu_s/T + cdotu*cdotdu_s/T**2 - udotdu_s/T)
       end do
       boltz(restvec) = weight(restvec)*(1.d0 - uke/T - udotdu_s/T)
       return
       
       
    case (4)
       ! Second order expansion for the equilibrium distribution
       ! and force implemmentation according to
       ! Phys. Rev. E 65,(2002): Guo
       
       ! u = u_tilde
       ! u_tilde + f_{Shan Chen}_r + tau_r*g_accn
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz
       
       uke = 0.5d0*(tux*tux + tuy*tuy + tuz*tuz)
       
       ! dd#  = (tau - 0.5)*g_accn_#
       dke = 0.5d0*(ddx*ddx + ddy*ddy + ddz*ddz)
       
       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz
          cdotd = cx(s)*ddx + cy(s)*ddy + cz(s)*ddz 
          
          boltz(s) = weight(s)*( 1.d0 + cdotu/T &
               + 0.5d0*((cdotu/T)**2 - (cdotd/T)**2) &
               - (uke-dke)/T )
       end do
       boltz(restvec) = weight(restvec)*(1.d0 - (uke-dke)/T)
       return
       
    case (5)
       ! Third order expansion for the equilibrium distribution
       ! and force implemmentation according to
       ! Phys. Rev. E 65,(2002): Guo
       
       ! ux = u_tilde
       ! u_tile + f_{Shan Chen}_r + tau_r*g_accn
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz
       
       uke = 0.5d0*(tux*tux + tuy*tuy + tuz*tuz)
       
       ! dd#  = (tau - 0.5)*g_accn
       dke = 0.5d0*(ddx*ddx + ddy*ddy + ddz*ddz)
       
       ! u_tile + f_{Shan Chen}_r + 0.5*g_accn
       tu1x = ux + dux - ddx
       tu1y = uy + duy - ddy
       tu1z = uz + duz - ddz
       
       u1ke = 0.5d0*(tu1x*tu1x + tu1y*tu1y + tu1z*tu1z)
       
       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz
          cdot1u = cx(s)*tu1x + cy(s)*tu1y + cz(s)*tu1z
          cdotd = cx(s)*ddx + cy(s)*ddy + cz(s)*ddz 
          
          boltz(s) = weight(s)*( 1.d0 + cdotu/T &
               + 0.5d0*((cdotu/T)**2 - (cdotd/T)**2) &
               - (uke-dke)/T &
               - uke*cdot1u/T**2 + (cdot1u/T)**3/6.d0 &
               + 0.5d0*(u1ke/T)**2 - 0.5*(cdot1u**2*u1ke)/T**3 &
               + (cdot1u/T)**4/24.d0 )
       end do
       boltz(restvec) = weight(restvec)*(1.d0 - (uke-dke)/T + 0.5d0*(uke/T)**2)
       return
       
       ! Corrected Luo 'operator'
       ! Extended to fourth order (as a test)
       ! Particle number and momentum are conserved - however
       ! I believe that other 'hydronamic moments' are not.
       ! Do not use/trust this code!
       ! However, it is interesting that this code agrees
       ! very closely with the Hudong code (bdist=0)
    case (6)
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz
       
       uke = 0.5d0*(tux*tux + tuy*tuy + tuz*tuz)
       
       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz
          
          boltz(s) = weight(s)*(1.d0 + cdotu/T &
               + 0.5d0*(cdotu/T)**2 - uke/T &
               - uke*cdotu/T**2 + (cdotu/T)**3/6.d0 &
               + 0.5d0*(uke/T)**2 - 0.5*(cdotu**2*uke)/T**3 &
               + (cdotu/T)**4/24.d0)
       end do
       boltz(restvec) = weight(restvec)*(1.d0 - uke/T + 0.5d0*(uke/T)**2)
       return

       ! This is a second order expansion of f(u+delta_u) with no
       ! terms missing PRL 81 1618 (1998) and PRE 58 6855 (1998) both
       ! claim that this is inconsistent with a systematic derivation
       ! of how force terms should be introduced. So this is the same
       ! as bdist==4, but optimized by manually unrolling the loop
       ! and exploiting zero-valued components of the lattice
       ! vectors.
    case (7)
       ! tux = ux+dux
       ! tuy = uy+duy
       ! tuz = uz+duz

       ! tuxt = tux*holdit
       ! tuyt = tuy*holdit
       ! tuzt = tuz*holdit

       ! tuxta = tuxt*aN1
       ! tuyta = tuyt*aN1
       ! tuzta = tuzt*aN1

       ! uke = 0.5d0*(tux*tuxt + tuy*tuyt + tuz*tuzt)

       ! const = 1.0_rk-uke
       ! const1 = aN1*const

       ! cdotu = tuxt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta
       ! boltz(1) = even+odd
       ! boltz(2) = even-odd

       ! cdotu = tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta
       ! boltz(3) = even+odd
       ! boltz(4) = even-odd

       ! cdotu = tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuzta
       ! boltz(5) = even+odd
       ! boltz(6) = even-odd

       ! cdotu = tuxt+tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta+tuyta
       ! boltz(7) = even+odd
       ! boltz(12) = even-odd

       ! cdotu = tuxt-tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta-tuyta
       ! boltz(8) = even+odd
       ! boltz(11) = even-odd

       ! cdotu = tuxt+tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta+tuzta
       ! boltz(9) = even+odd
       ! boltz(14) = even-odd

       ! cdotu = tuxt-tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta-tuzta
       ! boltz(10) = even+odd
       ! boltz(13) = even-odd

       ! cdotu = tuyt+tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta-tuzta
       ! boltz(15) = even+odd
       ! boltz(18) = even-odd

       ! cdotu = tuyt-tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta-tuzta
       ! boltz(16) = even+odd
       ! boltz(17) = even-odd

       ! boltz(restvec) = aN0*const
       return
    end select

end subroutine boltz_dist

#endif

!> returns only component \c s of \c boltz_dist() as defined above
subroutine boltz_dist_s(ux,uy,uz,dux,duy,duz,ddx,ddy,ddz,boltz_s,s)
    real(kind=rk),intent(in) :: ux,uy,uz,dux,duy,duz, ddx, ddy, ddz
    real(kind=rk),intent(out) :: boltz_s
    integer,intent(in) :: s
    real(kind=rk),parameter,dimension(nvecs) :: weight = (/ &
!!$         & 1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,& 
!!$         & 1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,&
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/3.0_rk /)
    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T
    real(kind=rk) :: cdotu,const,tux,tuy,tuz

    tux = ux+dux
    tuy = uy+duy
    tuz = uz+duz
    const = 1.0_rk-0.5_rk*(tux*tux + tuy*tuy + tuz*tuz)*holdit

    select case (s)
       case (1)
           cdotu = tux*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (2)
           cdotu = tux*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (3)
           cdotu = tuy*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (4)
           cdotu = tuy*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (5)
           cdotu = tuz*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (6)
           cdotu = tuz*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (7)
           cdotu = (tux+tuy)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (8)
           cdotu = (tux-tuy)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (9)
           cdotu = (tux+tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (10)
           cdotu = (tux-tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (11)
           cdotu = (tux-tuy)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (12)
           cdotu = (tux+tuy)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (13)
           cdotu = (tux-tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (14)
           cdotu = (tux+tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (15)
           cdotu = (tuy+tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (16)
           cdotu = (tuy-tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk+cdotu)+(weight(s)/2.0_rk)*cdotu**2+(weight(s)/6.0_rk)*cdotu**3
       case (17)
           cdotu = (tuy-tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (18)
           cdotu = (tuy+tuz)*holdit
           boltz_s = weight(s)*const*(1.0_rk-cdotu)+(weight(s)/2.0_rk)*cdotu**2-(weight(s)/6.0_rk)*cdotu**3
       case (19)
           boltz_s = weight(s)*const
    end select

end subroutine boltz_dist_s

!> returns only 0th and 1st order terms of \c boltz_dist_s() with
!> \c bdist=2 (3rd order bdist).
!>
!> \param[in] tux x velocity component
!>
!> \param[in] tux y velocity component
!>
!> \param[in] tux z velocity component
!>
!> \param[out] boltz_s component \c s of the distribution function
!>
!> \param[in] s specifying which component of the distribution
!> function to return in \c boltz_s
!>
!> If necessary, \c (/dux,duy,duz/) should be added to \c (/ux,uy,uz/)
!> and then passed to this subroutine as sum.
subroutine boltz_dist1_s(tux,tuy,tuz,boltz_s,s)
    real(kind=rk),intent(in) :: tux,tuy,tuz
    real(kind=rk),intent(out) :: boltz_s
    integer,intent(in) :: s
    real(kind=rk),parameter,dimension(nvecs) :: weight = (/ &
         & 1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,& 
         & 1.0_rk/18.0_rk, 1.0_rk/18.0_rk, 1.0_rk/18.0_rk,&
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/36.0_rk, 1.0_rk/36.0_rk, 1.0_rk/36.0_rk,& 
         & 1.0_rk/3.0_rk /)
    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T

    select case (s)
       case (1)
           boltz_s = weight(s)*(1.0_rk+tux*holdit)
       case (2)
           boltz_s = weight(s)*(1.0_rk-tux*holdit)
       case (3)
           boltz_s = weight(s)*(1.0_rk+tuy*holdit)
       case (4)
           boltz_s = weight(s)*(1.0_rk-tuy*holdit)
       case (5)
           boltz_s = weight(s)*(1.0_rk+tuz*holdit)
       case (6)
           boltz_s = weight(s)*(1.0_rk-tuz*holdit)
       case (7)
           boltz_s = weight(s)*(1.0_rk+(tux+tuy)*holdit)
       case (8)
           boltz_s = weight(s)*(1.0_rk+(tux-tuy)*holdit)
       case (9)
           boltz_s = weight(s)*(1.0_rk+(tux+tuz)*holdit)
       case (10)
           boltz_s = weight(s)*(1.0_rk+(tux-tuz)*holdit)
       case (11)
           boltz_s = weight(s)*(1.0_rk-(tux-tuy)*holdit)
       case (12)
           boltz_s = weight(s)*(1.0_rk-(tux+tuy)*holdit)
       case (13)
           boltz_s = weight(s)*(1.0_rk-(tux-tuz)*holdit)
       case (14)
           boltz_s = weight(s)*(1.0_rk-(tux+tuz)*holdit)
       case (15)
           boltz_s = weight(s)*(1.0_rk+(tuy+tuz)*holdit)
       case (16)
           boltz_s = weight(s)*(1.0_rk+(tuy-tuz)*holdit)
       case (17)
           boltz_s = weight(s)*(1.0_rk-(tuy-tuz)*holdit)
       case (18)
           boltz_s = weight(s)*(1.0_rk-(tuy+tuz)*holdit)
       case (19)
           boltz_s = weight(s)
    end select
end subroutine boltz_dist1_s

end module lb3d_bdist_module
