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

!> Contains lattice specific variables and functions
!> \todo This is not yet exhaustive / fully abstract
module lb3d_lattice_module

    use lb3d_global_module!,only:msgstr,rk
    use lb3d_config_module, only: nx,ny,nz!,bcsel
    use lb3d_log_module

    implicit none

    
    type lbe_site
       sequence	! Contiguous in memory, please..

       real(kind=rk), dimension(19) :: n_r  ! Red fluid velocity distribution  
       real(kind=rk) :: rho_r               ! Density of red fluid             
       real(kind=rk), dimension(3) :: u_r   ! Mean velocity of red fluid       
       real(kind=rk)	:: rock_colour_r    ! Wetting parameter                

#ifdef RELTIME
     real(kind=rk)	:: taupos_r
#endif


#ifndef SINGLEFLUID
       real(kind=rk), dimension(19) :: n_b  ! Blue fluid velocity distribution 
       real(kind=rk) :: rho_b               ! Density of blue fluid            
       real(kind=rk), dimension(3) :: u_b   ! Mean velocity of blue fluid      
       real(kind=rk)	:: rock_colour_b    ! Wetting parameter                

#ifdef RELTIME
       real(kind=rk) :: taupos_b
#endif


#endif
#ifndef NOSURFACTANT
       real(kind=rk), dimension(19) :: n_s  ! Surfactant fluid velocity distribution 
       real(kind=rk) :: rho_s               ! Density of surfactant fluid            
       real(kind=rk), dimension(3):: u_s    ! Mean velocity of surfactant fluid      
       real(kind=rk)	:: rock_colour_s    ! Wetting parameter                

       !FIXME cleanup needed, this is just to quick-check prior to halo-exchange re-write
 !      real(kind=rk), dimension(3)  :: d    ! Dipole moment vector
       real(kind=rk), dimension(3)  :: da    ! Dipole moment vector  
       real(kind=rk), dimension(3)  :: db    ! Dipole moment vector

#ifdef RELTIME
       real(kind=rk) :: taupos_s
#endif
  
#endif       

       real(kind=rk)	:: rock_state  ! (0,1) for rockp.

#ifdef DIST
     real(kind=rk)	:: abst
#endif
       
    end type lbe_site



  !> \{
  !> \name These parameters are a result solely of the lattice chosen.

    integer, parameter :: D = 4       !< No dimensions of lattice

  !> Number of different lattice vectors.
  !>
  !> Contains the number of independent vectors (19), which
  !> corresponds to the size of the \c n_r, \c n_b, and
  !> \c n_s arrays of each site. When operating on all elements
  !> of such an array, it is strongly advised that the loop run
  !> from \c 1 to \c nvecs, so that the code will not need
  !> changing should another lattice be used.
  integer, parameter :: nvecs = 19    !< No of vectors
  integer, parameter :: ngvecs = 25   !< No of vectors (inc. degen.)
  integer, parameter :: nnonrest = 18 !< No of nonrest vectors
  real(kind=rk) , parameter :: T = 1.d0 / 3.d0 !< Lattice temp.
!  integer, parameter :: nd = 3    !< Number of dimensions of the system.
  !> The index of the rest vector.
  !>
  !> When referring explicitly to the rest particles, it is strongly
  !> advised that an index of \c restvec rather than 19 is used, to
  !> ensure code readability and minimise problems should the lattice
  !> be changed.
  integer, parameter :: restvec = 19
  !> \}
  !>
  !> \{
  !> \name Lattice vectors and associated degeneracies.
  !>
  !> \verbatim
  !> the vectors are the columns of the following matrix
  !>
  !> 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19             (index)
  !>
  !> x,-x, y,-y, z,-z, x+-y, x+-z,-x+-y,-x+-z, y+-z,-y+-z, 0    (short notation)
  !>
  !> 1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0           (vectors)
  !> 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0
  !> 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1, 0
  !>
  !> 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1(projection weights)
  !>
  !> the weights here are not the lattice weights, but the weights
  !> used for the downprojection from a D3Q25 to a D3Q19 lattice
  !>
  !> 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12  (lattice weights)
  !> \endverbatim

  !> the X components of the lattice vectors
  integer, parameter :: cx(nvecs) = &
    (/ 1,-1, 0, 0, 0, 0, &
       1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, &
       0 &
    /)

  !> the Y components of the lattice vectors
  integer, parameter :: cy(nvecs) = &
    (/ 0, 0, 1,-1, 0, 0, &
       1,-1, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, &
       0 &
    /)

  !> the Z components of the lattice vectors
  integer, parameter :: cz(nvecs) = &
    (/ 0, 0, 0, 0, 1,-1, &
       0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1, &
       0 &
    /)

  !> degeneracies of the lattice vectors
  !> lb3d_bdist make use of 'proper' 3d-weights, so g is unity
  !FIXME removal of all g by undefining this and follow compiler panic trails
  integer, parameter :: g(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)
!!$  integer, parameter :: g(nvecs) = &
!!$    (/ 1, 1, 1, 1, 1, 1, &
!!$       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
!!$       1 &
!!$    /)

  !> contains the index of the vector pointing in the opposite
  !> direction to c(i).
  integer, parameter :: bounce(nvecs) = &
    (/  2, 1, 4, 3, 6, 5, &
       12,11,14,13, 8, 7, 10, 9, 18,17,16,15, &
       19 &
    /)

  integer, parameter :: nnp = 5 !< size of the \c negx,\c negy, etc arrays

  !> \{
  !> \name Arrays used in invasive flow.
  !>
  !> \c negx contains a list of the indices of each vector which
  !> has a negative X componentl \c posy contains a list of the
  !> indices of each vector which has a positive Y component, etc.
  integer, parameter :: negx(nnp) = (/ 1 , 7 , 8 , 9 , 10 /)
  integer, parameter :: negy(nnp) = (/ 3 , 7 , 11 , 15 , 16 /)
  integer, parameter :: negz(nnp) = (/ 5 , 9 , 13 , 15 , 17 /)
  integer, parameter :: posx(nnp) = (/ 2 , 11 , 12 , 13 , 14 /)
  integer, parameter :: posy(nnp) = (/ 4 , 8 , 12 , 17 , 18 /)
  integer, parameter :: posz(nnp) = (/ 6 , 10 , 14 , 16 , 18 /)
  !> \}

  !> Array containing the lattice vectors.
  !>
  !> contains a vector corresponding to the ith lattice vector, but has
  !> to be initialised in code, since you can't initialise an array of
  !> rank greater than 1.  Fortran sucks.
  integer, save :: c(nvecs,3)
  !> \}

!!$#ifdef MD
!!$  !> depth of halo region in each direction
!!$  !>
!!$  !> meant to be enlarged on demand. Currently, only in the MD branch
!!$  !> the value is enlarged and currently only in the MD branch the
!!$  !> enlargement has some effect.
!!$  !>
!!$  !> \note However, also there currently only the rock state in the
!!$  !> additional halo is communicated and this is done only once during
!!$  !> initialization. This is for performance reason. The code to
!!$  !> completely replace lbe_halo_exchange() was implemented already
!!$  !> and can be found in lbe_md_fluid.F90.
!!$#ifndef RWALK
!!$  integer,save :: halo_extent = 1
!!$#else
!!$  integer,save :: halo_extent = 4
!!$#endif
!!$
!!$#else
!!$  ! NOMD
  integer,save :: halo_extent = 1 
  integer,save :: force_halo_extent = 1
!!$#endif

  !> # of surrounding lattice nodes for an arbitrary off-lattice position
  integer,parameter :: n_lp=2**3
  !> relative coordinates of these nodes
  integer,save :: lp_sur(3,n_lp)
  !> relative coordinates of the opposing nodes
  integer,save :: opp_lp_sur(3,n_lp)

    ! This might be a queer hack but allows to widen the halo region just
    ! for md without need to touch anything outside the head of this file.
    ! Now  whole_N  is the whole lattice...
    type(lbe_site),dimension(:,:,:),allocatable,target :: whole_N
    ! ... N  is just a pointer to a part of  whole_N .
    type(lbe_site),pointer :: N(:,:,:)

  contains

  !> Initialise lattice related properties and allocate the lattice
    subroutine lb3d_lattice_init (stage)
      !      ==================================================================

      integer :: stage
      integer :: i,j,k
      character(len=256) :: msgstr
      integer            :: ierror
      integer,parameter          :: rk = 8

      select case (stage)

      case (3) ! mem alloc

         do i = 1, nvecs
            c(i,:) = (/ cx(i), cy(i), cz(i) /)
         end do

         ! Calculate relative coordinates for lattice points
         ! surrounding an arbitrary point
         combinations: do i = 1, n_lp
            dimensions: do k = 1, 3
               lp_sur(k,i) = mod((i-1)/2**(3-k),2)
            end do dimensions
         end do combinations
         opp_lp_sur(:,:) = 1 - lp_sur(:,:)

         ! allocate the arena.
         write (unit=msgstr,fmt='("Allocating arena with halo_extent=",I0,"...")') &
              &halo_extent
         call log_msg(trim(msgstr),.false.)
         allocate(whole_N(&
              &1-halo_extent:nx+halo_extent,&
              &1-halo_extent:ny+halo_extent,&
              &1-halo_extent:nz+halo_extent),stat=ierror)

         !> \warning I found no way to set \c lbound(N) to 0. Instead they
         !> are 1. However, this imposes no problem as long as nowhere in the
         !> main program there is a direct access to an element of \c
         !> N. Subroutines that receive \c N as an argument specify the

         !> \todo Not sure if pre-allocating helps, post is the way to go according
         !> to IBM Fortran, but seems to break things, so its left like this for now
         !allocate(N(0:nx+1,0:ny+1,0:nz+1),stat=ierror)
         N => whole_N(0:nx+1,0:ny+1,0:nz+1)
  
         !> starting index anyway.

         if (ierror .ne. 0) then
            CALL log_msg("FATAL ERROR: Unable to allocate arena. Aborting...",.true.)
            CALL MPI_Abort
         end if

!!!
         ! This might be important to avoid random prebouncing at
         ! last_periodic_z+1 regarding uninitialized rock_state at
         ! last_periodic_z if last_periodic_z is in a lower z halo layer ???
         whole_N%rock_state = 0.0_rk
!!!
      case default

#ifdef LB3D_DEBUG_INFO
       call log_msg('nothing to do.',.false.) 
#endif

      end select
    end subroutine lb3d_lattice_init

 
  subroutine lb3d_lattice_clear(N)
    integer :: x,y,z
    type(lbe_site),dimension(0:,0:,0:) :: N
    do x = 0, nx+1
       do y = 0, ny+1
          do z = 0, nz+1
             N(x,y,z)%n_r = 0.d0
#ifndef SINGLEFLUID
             N(x,y,z)%n_b = 0.d0
#endif
             N(x,y,z)%rock_state = 0.d0
#ifndef NOSURFACTANT
             N(x,y,z)%n_s = 0.d0
             N(x,y,z)%da = 0.d0
             N(x,y,z)%db = 0.d0
#endif
          end do
       end do
    end do
    
  end subroutine lb3d_lattice_clear

  subroutine lb3d_lattice_clear_site(N,x,y,z)
    integer :: x,y,z
    type(lbe_site),dimension(0:,0:,0:) :: N

             N(x,y,z)%n_r = 0.d0
#ifndef SINGLEFLUID
             N(x,y,z)%n_b = 0.d0
#endif
             N(x,y,z)%rock_state = 0.d0
#ifndef NOSURFACTANT
             N(x,y,z)%n_s = 0.d0
             N(x,y,z)%da = 0.d0
             N(x,y,z)%db = 0.d0
#endif
  end subroutine lb3d_lattice_clear_site

 end module lb3d_lattice_module
