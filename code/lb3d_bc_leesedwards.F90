#include "lb3d.h"

!> Lees-Edwards boundary conditions, or other boundary conditions
!> which involve changing the velocity and/or position of particles
!>
!> In the case of LE to produce shear.
!>
!>When special boundary conditions are in force, code here replaces
!>the halo exchange code in lbe_parallel_module.
module lb3d_bc_leesedwards_module

    use lb3d_bdist_module
    use lb3d_mpi_parameters_module
    use lb3d_config_module, only: amass_b,amass_r,amass_s,omega_b,omega_r,omega_s&
         &,shear_omega,shear_u
    use lb3d_log_module
    implicit none

    !FIXME not sure where to put the dimensionality, took it out, because it
    !is creating difficulties with the dependencies and is not implemented 
    !consistently anyway
   ! integer, parameter :: nd = 3

    integer :: lelo            ! Lattice offset
    real(kind=rk)  :: lefrac          ! Fraction used in interpolation
    real(kind=rk)  :: shear_sum = 0.d0

    !Array of ranks of processors will need to write to in LE scheme
    !Two arrays are needed as a processor can be both on the top and
    !bottom.
    integer, dimension(:), allocatable :: topprocs
    integer, dimension(:), allocatable :: botprocs

    ! Ranks of Lees Edwards neighbours (p|m,chunk1|chunk2)
    integer :: leprocs(2,3)=-1

    ! Unique tags to identify which part of a halo transfer
    ! data belongs to.
    integer, parameter :: tag_le_px1 = 9
    integer, parameter :: tag_le_px2 = 10
    integer, parameter :: tag_le_px3 = 11
    integer, parameter :: tag_le_px4 = 12
    integer, parameter :: tag_le_mx1 = 13
    integer, parameter :: tag_le_mx2 = 14
    integer, parameter :: tag_le_mx3 = 15
    integer, parameter :: tag_le_mx4 = 16
contains

!>Perform initial tasks for Lees Edwards code:
!>Allocate and fill topprocs and botprocs - arrays holding ranks
!>of processors this rank will need to exchange information with.
    subroutine le_init
        integer :: ierror
        integer :: i              ! looping variable

        ! If this processor is on minimum x-plane, remember the
        ! choice of top processors.
        if (ccoords(1) .eq. 0) then
           allocate(topprocs(cdims(3)))
           do i = 1, cdims(3)
              call MPI_Cart_Rank(Comm_Cart, &
                   (/ cdims(1)-1, ccoords(2), i-1 /), &
                   topprocs(i), ierror )
           enddo
        endif

        ! and/or if processor is on max x-plane...
        if ((ccoords(1)+1) .eq. cdims(1)) then
           allocate(botprocs(cdims(3)))
           do i = 1, cdims(3)
              call MPI_Cart_Rank(Comm_Cart, &
                   (/ 0, ccoords(2), i-1 /), &
                   botprocs(i),ierror )
           enddo
        endif

        print *, myrankc, ccoords
        print *,topprocs
        print *,botprocs
    end subroutine le_init

    !>Set up Lees Edwards for current timestep. leprocs stores the
    !>ranks of this processor's neighbours for interfaces where
    !>boundary conditions are strange or changing. It replaces nnprocs
    !>when not equal to -1.
    subroutine le_neighbours()
        character(len=128) :: msgstr
        integer :: proco = 0

	! If doing rotated boundary conditions
	if (inv_fluid .eq. 6) then
           shear_u = 0.d0
           lelo = nz/2 + 1
           lefrac = 0.d0
           proco = cdims(3)/2
        else ! inv_fluid == 5
           shear_sum = shear_sum + shear_u
           lelo = mod(int(2.*shear_sum),nz)

           lefrac = 2.d0*shear_sum
           lefrac = abs(lefrac - aint(lefrac))
           proco = int(2.d0*shear_sum/nz)

           if (lelo .lt. 0) then
              write(msgstr,"('Lelo = ',I0,' , shear_sum = ',F16.10,"&
                   &//"' , shear_u = ',F16.10,' , shear_omega = ',F16.10)") & 
                   &lelo, shear_sum, shear_u, shear_omega
              call log_msg(trim(msgstr),.false.)         
              write(msgstr,"('lefrac = ',F16.10,' , proco = ',I0)") &
                   &lefrac, proco
              call log_msg(trim(msgstr),.false.)
         
           endif

        endif

        ! If this processor is on minimum x-plane, find the
        ! two processors to write to.
        if (ccoords(1) .eq. 0) then
           !leprocs changed order, one added, as suggested by Jens Oct 2003
           leprocs(1,1) = topprocs(modulo(ccoords(3)+proco-1,cdims(3))+1)
           leprocs(1,2) = topprocs(modulo(ccoords(3)+proco  ,cdims(3))+1)
           leprocs(1,3) = topprocs(modulo(ccoords(3)+proco+1,cdims(3))+1)
        endif

        ! and/or if processor is on max x-plane...
        if (ccoords(1) .eq. cdims(1) - 1) then
           !leproc(2,3) added as suggested by Jens, Oct 2003
           leprocs(2,1) = botprocs(modulo(ccoords(3)-proco-1,cdims(3))+1)
           leprocs(2,2) = botprocs(modulo(ccoords(3)-proco  ,cdims(3))+1)
           leprocs(2,3) = botprocs(modulo(ccoords(3)-proco+1,cdims(3))+1)
        endif

#ifdef DEBUG_LE
 ! Who I am and where I'm going:
	print *, 'Rank ',myrankc,' writing to', &
             leprocs(1,1),leprocs(1,2),leprocs(1,3),leprocs(2,1),leprocs(2,2),leprocs(2,3), &
             nnprocs(1,1),nnprocs(1,2)
#endif
    end subroutine le_neighbours


    !> Performs a halo exchange subject to Lees-Edwards Boundary
    !> conditions on the maximum and minimum x-planes
    !>
    !> This is an extension of the normal halo exchange function.
    subroutine le_halo_exchange(N)
        type(lbe_site),dimension(0:,0:,0:) :: N
        real(kind=rk),dimension(:,:,:),allocatable :: pinbuf
        real(kind=rk),dimension(:,:,:),allocatable :: poutbuf
        real(kind=rk),dimension(:,:,:),allocatable :: minbuf
        real(kind=rk),dimension(:,:,:),allocatable :: moutbuf

        ! These are some defined request numbers
        integer,parameter :: pin=1, min=2, pout=3, mout=4
        ! xreq allows extra requests - stores number of next
        ! undefined request.
        integer :: xreq=5

        ! FIXME: Should allocate this properly
        ! - as it is we'll not get more than 8 requests
        integer, dimension(10) :: requests

        integer, dimension(MPI_STATUS_SIZE,10) :: statuses
        integer :: ierror

        integer :: x,y,z
        integer :: ibuf      ! variable, depending on nosurfactant or surfactant
        integer :: i,j,k,s          ! looping variables
        character(len=80) :: bla

        ! These used only for velocity adjusting part ------------------------
        real(kind=rk),dimension(nvecs) :: boltz_du, boltz_u, boltz_u_du
        real(kind=rk),dimension(3) :: p_tilde, u_tilde
        real(kind=rk) :: rho_tilde, rho_tmp

        real(kind=rk) :: rho_r,rN_s,p_r(3),df_r(nvecs),fnew_r(nvecs)
#ifndef SINGLEFLUID
        real(kind=rk) :: rho_b,bN_s,p_b(3),df_b(nvecs),fnew_b(nvecs)
#ifndef NOSURFACTANT
        real(kind=rk) :: rho_s,sN_s,p_s(3),df_s(nvecs),fnew_s(nvecs)
#endif
#endif

        real(kind=rk) :: df_scale
        real(kind=rk),parameter :: zero=0.0_rk
        !----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Swap in the X direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!
        ! Begin asynchronous receives
!!!!!!

        xreq = 5

#ifndef NOSURFACTANT
        ibuf = 67
#else
        ibuf = 39
#endif
#ifdef SINGLEFLUID
        ibuf = 20
#endif

	! Receive from top into bottom halo
	if (leprocs(1,2) .eq. -1) then
	   allocate(pinbuf(ibuf,ny+2,nz+2))
	   call MPI_Irecv(	pinbuf,		& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,1),	& ! source
                tag_px,		& ! tag
                Comm_Cart,	& ! communicator
                requests(1),	& ! request
                ierror)
	   call MPI_Irecv(	xreq,		& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,1),	& ! source
                1,		& ! tag
                Comm_Cart,	& ! communicator
                requests(2),	& ! request
                ierror)
           bla = 'Bad +x async receive'
	   call checkmpi(ierror,bla)
	else
           allocate(pinbuf(ibuf,ny+2,nz+3))
           if(lelo.lt.0)then
              call MPI_Irecv(	pinbuf(1,1,3-lelo), &
                   ibuf*(ny+2)*(nz+lelo+1),          &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! source
                   tag_le_px3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(1),	& ! request
                   ierror)
              bla = 'Bad +x_le1 async receive'
              call checkmpi(ierror,bla)
              call MPI_Irecv(	pinbuf(1,1,1), &
                   ibuf*(ny+2)*(-lelo+2),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(1,1),	& ! source
                   tag_le_px4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(2),	& ! request
                   ierror)

           else !lelo >=0
              call MPI_Irecv(	pinbuf(1,1,nz-lelo+2), &
                   ibuf*(ny+2)*(lelo+2),          &
                   LBE_REAL,	& ! datatype
                   leprocs(1,3),	& ! source
                   tag_le_px1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(1),	& ! request
                   ierror)
              bla = 'Bad +x_le1 async receive'
              call checkmpi(ierror,bla)
              call MPI_Irecv(	pinbuf(1,1,1), &
                   ibuf*(ny+2)*(nz-lelo+1),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! source
                   tag_le_px2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(2),	& ! request
                   ierror)

              bla = 'Bad +x_le2 async receive'
              call checkmpi(ierror,bla)
              xreq = xreq + 1
           endif !lelo >= or < 0
	endif

 ! Receive from bottom into top halo
	if (leprocs(2,1) .eq. -1) then
	   allocate(minbuf(ibuf,ny+2,nz+2))
	   call MPI_Irecv(	minbuf,		& ! bufh
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,2),	& ! source
                tag_mx,		& ! tag
                Comm_Cart,	& ! communicator
                requests(3),	& ! request
                ierror)
	   call MPI_Irecv(	xreq,		& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,2),	& ! source
                2,		& ! tag
                Comm_Cart,	& ! communicator
                requests(4),	& ! request
                ierror)
           bla = 'Bad -x async receive'
	   call checkmpi(ierror,bla)
	else
	   allocate(minbuf(ibuf,ny+2,nz+3))
           if(lelo.lt.0)then

              call MPI_Irecv(	minbuf(1,1,nz+lelo+2),  &
                   ibuf*(ny+2)*(2-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,3),	& ! source
                   tag_le_mx3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(3),	& ! request
                   ierror)
              bla = 'Bad -x_le1 async receive'
              call checkmpi(ierror,bla)
              call MPI_Irecv(	minbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+1+lelo),    &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! source
                   tag_le_mx4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(4),	& ! request
                   ierror)

           else !lelo >= 0

              call MPI_Irecv(	minbuf(1,1,lelo+3),  &
                   ibuf*(ny+2)*(nz+1-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! source
                   tag_le_mx1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(3),	& ! request
                   ierror)
              bla = 'Bad -x_le1 async receive'
              call checkmpi(ierror,bla)
              call MPI_Irecv(	minbuf(1,1,1), &
                   ibuf*(ny+2)*(lelo+2),    &
                   LBE_REAL,	& ! datatype
                   leprocs(2,1),	& ! source
                   tag_le_mx2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(4),	& ! request
                   ierror)

              bla = 'Bad -x_le2 async receive'
              call checkmpi(ierror,bla)
              xreq = xreq + 1
           endif !lelo < or >= 0
	endif

!!!!!!
 ! Begin asynchronous sends
!!!!!!

 ! Send from top buffer destined for bottom halo
	if (leprocs(2,1) .eq. -1) then
	   allocate(poutbuf(ibuf,ny+2,nz+2))

           do k=1,nz+2
              do j=1,ny+2
                 poutbuf(1:19,j,k) = N(nx,j-1,k-1)%n_r(:)
#ifndef SINGLEFLUID
                 poutbuf(20:38,j,k) = N(nx,j-1,k-1)%n_b(:)
#endif
#ifndef NOSURFACTANT
                 poutbuf(39:57,j,k) = N(nx,j-1,k-1)%n_s(:)
                 poutbuf(58:60,j,k) = N(nx,j-1,k-1)%da(:)
                 poutbuf(61:63,j,k) = N(nx,j-1,k-1)%db(:)
#endif
                 poutbuf(ibuf,j,k) = N(nx,j-1,k-1)%rock_state
              enddo
           enddo
	   call MPI_ISend(	poutbuf,	& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,2),	& ! dest
                tag_px,		& ! tag
                Comm_Cart,	& ! communicator
                requests(5),	& ! request
                ierror)
	   call MPI_ISend(	2,	& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,2),	& ! dest
                1,		& ! tag
                Comm_Cart,	& ! communicator
                requests(6),	& ! request
                ierror)

           bla = 'Bad +x async send'
	   call checkmpi(ierror,bla)
	else
           allocate(poutbuf(ibuf,ny+2,nz+3))
           if(lelo.lt.0)then
              do k=1,nz+1+lelo
                 do j=1,ny+2
                    poutbuf(1:19,j,k)  = N(nx,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    poutbuf(20:38,j,k) = N(nx,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    poutbuf(39:57,j,k) = N(nx,j-1,k)%n_s(:)
                    poutbuf(58:60,j,k) = N(nx,j-1,k)%da(:)
                    poutbuf(61:63,j,k) = N(nx,j-1,k)%db(:)

#endif
                    poutbuf(ibuf,j,k) = N(nx,j-1,k)%rock_state
                 enddo
              enddo

              do k=nz+2+lelo,nz+3
                 do j=1,ny+2 
                    poutbuf(1:19,j,k)  = N(nx,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    poutbuf(20:38,j,k) = N(nx,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    poutbuf(39:57,j,k) = N(nx,j-1,k-3)%n_s(:)
                    poutbuf(58:60,j,k) = N(nx,j-1,k-3)%da(:)
                    poutbuf(61:63,j,k) = N(nx,j-1,k-3)%db(:)
#endif
                    poutbuf(ibuf,j,k) = N(nx,j-1,k-3)%rock_state
                 enddo
              enddo
              call MPI_ISend(	poutbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+lelo+1),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! dest
                   tag_le_px3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(5),	& ! request
                   ierror)
              bla = 'Bad +x async send'
              call checkmpi(ierror,bla)
              call MPI_ISend(	poutbuf(1,1,nz+2+lelo), &
                   ibuf*(ny+2)*(-lelo+2),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,3),	& ! dest
                   tag_le_px4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(6),	& ! request
                   ierror)           
           else  !lelo >= 0

              do k=1,lelo+2
                 do j=1,ny+2           
                    poutbuf(1:19,j,k)  = N(nx,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    poutbuf(20:38,j,k) = N(nx,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    poutbuf(39:57,j,k) = N(nx,j-1,k)%n_s(:)
                    poutbuf(58:60,j,k) = N(nx,j-1,k)%da(:)
                    poutbuf(61:63,j,k) = N(nx,j-1,k)%db(:)
#endif
                    poutbuf(ibuf,j,k) = N(nx,j-1,k)%rock_state
                 enddo
              enddo

              do k=lelo+3,nz+3
                 do j=1,ny+2 
                    poutbuf(1:19,j,k)  = N(nx,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    poutbuf(20:38,j,k) = N(nx,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    poutbuf(39:57,j,k) = N(nx,j-1,k-3)%n_s(:)
                    poutbuf(58:60,j,k) = N(nx,j-1,k-3)%da(:)
                    poutbuf(61:63,j,k) = N(nx,j-1,k-3)%db(:)

#endif
                    poutbuf(ibuf,j,k) = N(nx,j-1,k-3)%rock_state
                 enddo
              enddo
              call MPI_ISend(	poutbuf(1,1,1), &
                   ibuf*(ny+2)*(lelo+2),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(2,1),	& ! dest
                   tag_le_px1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(5),	& ! request
                   ierror)
              bla = 'Bad +x async send'
              call checkmpi(ierror,bla)
              call MPI_ISend(	poutbuf(1,1,lelo+3), &
                   ibuf*(ny+2)*(nz-lelo+1),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! dest
                   tag_le_px2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(6),	& ! request
                   ierror)
              bla = 'Bad +x async send'
              call checkmpi(ierror,bla)
              xreq = xreq + 1
           endif ! lelo>= or < 0
	endif


 ! Send from bottom buffer destined for top halo
	if (leprocs(1,2) .eq. -1) then
	   allocate(moutbuf(ibuf,ny+2,nz+2))

           do k=1,nz+2
              do j=1,ny+2
                 moutbuf(1:19,j,k) = N(1,j-1,k-1)%n_r(:)
#ifndef SINGLEFLUID
                 moutbuf(20:38,j,k) = N(1,j-1,k-1)%n_b(:)
#endif
#ifndef NOSURFACTANT
                 moutbuf(39:57,j,k) = N(1,j-1,k-1)%n_s(:)
                 moutbuf(58:60,j,k) = N(1,j-1,k-1)%da(:)
                 moutbuf(61:63,j,k) = N(1,j-1,k-1)%db(:)

#endif
                 moutbuf(ibuf,j,k) = N(1,j-1,k-1)%rock_state
              enddo
           enddo

	   call MPI_ISend(	moutbuf,	& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,1),	& ! dest
                tag_mx,		& ! tag
                Comm_Cart,	& ! communicator
                requests(7),	& ! request
                ierror)
	   call MPI_ISend(	2,	& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,1),	& ! dest
                2,		& ! tag
                Comm_Cart,	& ! communicator
                requests(8),	& ! request
                ierror)
           bla = 'Bad -x async send'
	   call checkmpi(ierror,bla)
	else
           allocate(moutbuf(ibuf,ny+2,nz+3))
           if(lelo.lt.0)then

              do k=1,-lelo+2
                 do j=1,ny+2
                    moutbuf(1:19,j,k) = N(1,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    moutbuf(20:38,j,k) = N(1,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    moutbuf(39:57,j,k) = N(1,j-1,k)%n_s(:)
                    moutbuf(58:60,j,k) = N(1,j-1,k)%da(:)
                    moutbuf(61:63,j,k) = N(1,j-1,k)%db(:)

#endif
                    moutbuf(ibuf,j,k) = N(1,j-1,k)%rock_state
                 enddo
              enddo

              do k=3-lelo,nz+3
                 do j=1,ny+2
                    moutbuf(1:19,j,k) = N(1,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    moutbuf(20:38,j,k) = N(1,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    moutbuf(39:57,j,k) = N(1,j-1,k-3)%n_s(:)
                    moutbuf(58:60,j,k) = N(1,j-1,k-3)%da(:)
                    moutbuf(61:63,j,k) = N(1,j-1,k-3)%db(:)

#endif
                    moutbuf(ibuf,j,k) = N(1,j-1,k-3)%rock_state
                 enddo
              enddo
              call MPI_ISend(	moutbuf(1,1,1), &
                   ibuf*(ny+2)*(2-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(1,1),	& ! dest
                   tag_le_mx3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(7),	& ! request
                   ierror)
              bla = 'Bad +x async send'
              call checkmpi(ierror,bla)
              call MPI_ISend(	moutbuf(1,1,3-lelo), &
                   ibuf*(ny+2)*(nz+1+lelo), &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! dest
                   tag_le_mx4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(8),	& ! request
                   ierror)

           else !lelo >= 0

              do k=1,nz+1-lelo
                 do j=1,ny+2
                    moutbuf(1:19,j,k) = N(1,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    moutbuf(20:38,j,k) = N(1,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    moutbuf(39:57,j,k) = N(1,j-1,k)%n_s(:)
                    moutbuf(58:60,j,k) = N(1,j-1,k)%da(:)
                    moutbuf(61:63,j,k) = N(1,j-1,k)%db(:)

#endif
                    moutbuf(ibuf,j,k) = N(1,j-1,k)%rock_state
                 enddo
              enddo

              do k=nz+2-lelo,nz+3
                 do j=1,ny+2
                    moutbuf(1:19,j,k) = N(1,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    moutbuf(20:38,j,k) = N(1,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    moutbuf(39:57,j,k) = N(1,j-1,k-3)%n_s(:)
                    moutbuf(58:60,j,k) = N(1,j-1,k-3)%da(:)
                    moutbuf(61:63,j,k) = N(1,j-1,k-3)%db(:)

#endif
                    moutbuf(ibuf,j,k) = N(1,j-1,k-3)%rock_state
                 enddo
              enddo
              call MPI_ISend(	moutbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+1-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! dest
                   tag_le_mx1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(7),	& ! request
                   ierror)
              bla = 'Bad +x async send'
              call checkmpi(ierror,bla)
              call MPI_ISend(	moutbuf(1,1,nz+2-lelo), &
                   ibuf*(ny+2)*(lelo+2), &
                   LBE_REAL,	& ! datatype
                   leprocs(1,3),	& ! dest
                   tag_le_mx2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(8),	& ! request
                   ierror)
              xreq = xreq + 1
           endif !lwlo < or >= 0
	endif

!!!!!!
 ! Now wait for all I/O to complete
!!!!!!

	call MPI_Waitall(8,requests(1:8),statuses,ierror)

        bla = 'MPI_Waitall() failed in X direction'
	call checkmpi(ierror,bla)

!!!!!!
 ! Interpolate points and write halos:
!!!!!!

 ! Write bottom halo

	if (leprocs(1,2) .ne. -1) then
    ! FIXME: So far I've failed to get a terse array version of this
    ! with an overloaded * operator working
	   do y=0,ny+1
	      do z=0,nz+1
		 N(0,y,z)%n_r(:) = (1-lefrac)*pinbuf(1:19,y+1,z+1) + &
                      lefrac*pinbuf(1:19,y+1,z+2)
#ifndef SINGLEFLUID
		 N(0,y,z)%n_b(:) = (1-lefrac)*pinbuf(20:38,y+1,z+1) + &
                      lefrac*pinbuf(20:38,y+1,z+2)
#endif
#ifndef NOSURFACTANT
		 N(0,y,z)%n_s(:) = (1-lefrac)*pinbuf(39:57,y+1,z+1) + & 
                      lefrac*pinbuf(39:57,y+1,z+2)
		 N(0,y,z)%da(:) = (1-lefrac)*pinbuf(58:60,y+1,z+1) +  &
                      lefrac*pinbuf(58:60,y+1,z+2)
		 N(0,y,z)%db(:) = (1-lefrac)*pinbuf(61:63,y+1,z+1) +  &
                      lefrac*pinbuf(61:63,y+1,z+2)
#endif
		 N(0,y,z)%rock_state = 0.0_rk

   !d ----- under construction - change velocity at boundary ------------

   ! bN_s = density of blue particles at this site, etc.
   ! p_b = velocity of blue particles at this site, etc.

		 p_r = 0.d0
		 rN_s = 0.d0
#ifndef SINGLEFLUID
		 p_b = 0.d0
		 bN_s = 0.d0
#ifndef NOSURFACTANT
		 p_s = 0.d0
		 sN_s = 0.d0
#endif
#endif

		 do s=1,nvecs
		    rN_s = rN_s + N(0,y,z)%n_r(s) * g(s)
		    p_r = p_r + N(0,y,z)%n_r(s)*g(s)*c(s,:)
#ifndef SINGLEFLUID
		    p_b = p_b + N(0,y,z)%n_b(s)*g(s)*c(s,:)
		    bN_s = bN_s + N(0,y,z)%n_b(s) * g(s)
#endif
#ifndef NOSURFACTANT
		    sN_s = sN_s + N(0,y,z)%n_s(s) * g(s)
		    p_s = p_s + N(0,y,z)%n_s(s)*g(s)*c(s,:)
#endif
		 end do

   ! Calculate weighted total momentum  p_tilde ,
   ! rho_i = mass of species i at this site, and
   ! rho_tilde = sum of masses over relaxation times.

                 p_tilde = amass_r*omega_r*p_r
		 rho_r = amass_r*rN_s
                 rho_tilde = rho_r*omega_r
#ifndef SINGLEFLUID
                 p_tilde = p_tilde + amass_b*omega_b*p_b
		 rho_b = amass_b*bN_s
                 rho_tilde = rho_tilde + rho_b*omega_b
#ifndef NOSURFACTANT
                 p_tilde = p_tilde + amass_s*omega_s*p_s
		 rho_s = amass_s*sN_s
                 rho_tilde = rho_tilde + rho_s*omega_s
#endif
#endif

		 ! rho_tmp = clipped version of this.

		 rho_tmp = max(rho_tilde,dble(10.e-9))

   ! Calculate averaged velocity

		 u_tilde = p_tilde / rho_tmp


   ! FIXME This should be optimized.

		 call boltz_dist(u_tilde(1),u_tilde(2),&
                      &u_tilde(3)-shear_u-shear_u,&
                      &zero,zero,zero,0.0_rk,0.0_rk,0.0_rk,boltz_u_du)

		 call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3),&
                      &zero,zero,zero,0.0_rk,0.0_rk,0.0_rk,boltz_u)

		 boltz_du = boltz_u_du - boltz_u

   !Should be zero (or very little) change in flux in x direction
   !print *,sum(boltz_du*cx*g)

		 df_r = rN_s * boltz_du
		 fnew_r = N(0,y,z)%n_r + df_r
#ifndef SINGLEFLUID
		 df_b = bN_s * boltz_du
		 fnew_b = N(0,y,z)%n_b + df_b
#ifndef NOSURFACTANT
		 df_s = sN_s * boltz_du
		 fnew_s = N(0,y,z)%n_s + df_s
#endif
#endif

   ! If the velocity adjustment would make the site density
   ! negative anywhere, then scale down the velocity adjustment
   ! so that it will not
		 if (minval(fnew_r)<0.0_rk&
#ifndef SINGLEFLUID
                    &.or.minval(fnew_b)<0.0_rk&
#ifndef NOSURFACTANT
                    &.or.minval(fnew_s)<0.0_rk&
#endif
#endif
                    &) then

      ! Avoid dividing by zero
                    if (minval((/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))>0.0_rk) then
                       df_scale = maxval(-(/df_r&
#ifndef SINGLEFLUID
                       &,df_b&
#ifndef NOSURFACTANT
                       &,df_s&
#endif
#endif
                       &/)/(/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))
                    else
                       df_scale = 0.0_rk
                    end if
                    if (abs(df_scale).gt.1.0_rk) then
                       fnew_r = N(0,y,z)%n_r + df_r/df_scale
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b + df_b/df_scale
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s + df_s/df_scale
#endif
#endif
                    else
                       fnew_r = N(0,y,z)%n_r
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s
#endif
#endif
                    endif
                 endif
                 N(0,y,z)%n_r = fnew_r
#ifndef SINGLEFLUID
                 N(0,y,z)%n_b = fnew_b
#ifndef NOSURFACTANT
                 N(0,y,z)%n_s = fnew_s
#endif
#endif
              enddo
           enddo
	else
           do k=1,nz+2
              do j=1,ny+2
                 N(0,j-1,k-1)%n_r(:) = pinbuf(1:19,j,k)
#ifndef SINGLEFLUID
                 N(0,j-1,k-1)%n_b(:) = pinbuf(20:38,j,k)
#endif
#ifndef NOSURFACTANT
                 N(0,j-1,k-1)%n_s(:) = pinbuf(39:57,j,k)
                 N(0,j-1,k-1)%da(:) = pinbuf(58:60,j,k)
                 N(0,j-1,k-1)%db(:) = pinbuf(61:63,j,k)

#endif
                 N(0,j-1,k-1)%rock_state = pinbuf(ibuf,j,k)
              enddo
           enddo
	endif

 ! Write top halo
	if (leprocs(2,1) .ne. -1) then
    ! FIXME: So far I've failed to get a terse array version of this
    ! with an overloaded * operator working
	   do y=0,ny+1
	      do z=0,nz+1
		 N(nx+1,y,z)%n_r(:) = lefrac*minbuf(1:19,y+1,z+1) &
                      &+ (1-lefrac)*minbuf(1:19,y+1,z+2)
#ifndef SINGLEFLUID
		 N(nx+1,y,z)%n_b(:) = lefrac*minbuf(20:38,y+1,z+1) &
                      &+ (1-lefrac)*minbuf(20:38,y+1,z+2)
#ifndef NOSURFACTANT
		 N(nx+1,y,z)%n_s(:) = lefrac*minbuf(39:57,y+1,z+1) &
                      &+ (1-lefrac)*minbuf(39:57,y+1,z+2)
		 N(nx+1,y,z)%da(:) = lefrac*minbuf(58:60,y+1,z+1) &
                      &+ (1-lefrac)*minbuf(58:60,y+1,z+2)
		 N(nx+1,y,z)%db(:) = lefrac*minbuf(61:63,y+1,z+1) &
                      &+ (1-lefrac)*minbuf(61:63,y+1,z+2)

#endif
#endif
		 N(nx+1,y,z)%rock_state = 0.d0

   ! bN_s = density of blue particles at this site, etc.
   ! p_b = velocity of blue particles at this site, etc.

                 p_r = 0.d0
                 rN_s = 0.d0
#ifndef SINGLEFLUID
                 p_b = 0.d0
                 bN_s = 0.d0
#ifndef NOSURFACTANT
                 p_s = 0.d0
                 sN_s = 0.d0
#endif
#endif
                 do s=1,nvecs
                    rN_s = rN_s + N(nx+1,y,z)%n_r(s) * g(s)
                    p_r = p_r + N(nx+1,y,z)%n_r(s)*g(s)*c(s,:)
#ifndef SINGLEFLUID
                    bN_s = bN_s + N(nx+1,y,z)%n_b(s) * g(s)
                    p_b = p_b + N(nx+1,y,z)%n_b(s)*g(s)*c(s,:)
#ifndef NOSURFACTANT
                    sN_s = sN_s + N(nx+1,y,z)%n_s(s) * g(s)
                    p_s = p_s + N(nx+1,y,z)%n_s(s)*g(s)*c(s,:)
#endif
#endif
                 end do

                 ! Calculate weighted total momentum
                 ! rho_i = mass of species i at this site.
                 ! rho_tilde = sum of masses over relaxation times.
                 p_tilde = amass_r*omega_r*p_r
                 rho_r = amass_r*rN_s
                 rho_tilde = rho_r*omega_r
#ifndef SINGLEFLUID
                 p_tilde = p_tilde + amass_b*omega_b*p_b
                 rho_b = amass_b*bN_s
                 rho_tilde = rho_tilde + rho_b*omega_b
#ifndef NOSURFACTANT
                 p_tilde = p_tilde + amass_s*omega_s*p_s
                 rho_s = amass_s*sN_s
                 rho_tilde = rho_tilde + rho_s*omega_s
#endif
#endif

                 ! rho_tmp = clipped version of this.
                 rho_tmp = max(rho_tilde,dble(10.e-9))

                 ! Calculate averaged velocity
                 u_tilde = p_tilde / rho_tmp

		 ! FIXME This should be optimized.

		 call boltz_dist(u_tilde(1),u_tilde(2)&
                      &,u_tilde(3)+shear_u+shear_u&
                      &,zero,zero,zero,0.0_rk,0.0_rk,0.0_rk,boltz_u_du)
		 call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3)&
                      &,zero,zero,zero,0.0_rk,0.0_rk,0.0_rk,boltz_u)

		 boltz_du = boltz_u_du - boltz_u

		 df_r = rN_s * boltz_du
		 fnew_r = N(nx+1,y,z)%n_r + df_r
#ifndef SINGLEFLUID
		 df_b = bN_s * boltz_du
		 fnew_b = N(nx+1,y,z)%n_b + df_b
#ifndef NOSURFACTANT
		 df_s = sN_s * boltz_du
		 fnew_s = N(nx+1,y,z)%n_s + df_s
#endif
#endif

   ! If the velocity adjustment would make the site density
   ! negative anywhere, then scale down the velocity adjustment
   ! so that it will not
		 if (minval(fnew_r)<0.0_rk&
#ifndef SINGLEFLUID
                    &.or.minval(fnew_b)<0.0_rk&
#ifndef NOSURFACTANT
                    &.or.minval(fnew_s)<0.0_rk&
#endif
#endif
                    &) then

      ! Avoid dividing by zero
                    if (minval((/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))>0.0_rk) then
                       df_scale = maxval(-(/df_r&
#ifndef SINGLEFLUID
                       &,df_b&
#ifndef NOSURFACTANT
                       &,df_s&
#endif
#endif
                       &/)/(/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))
                    else
                       df_scale = 0.0_rk
                    end if
                    if (abs(df_scale).gt.1.0_rk) then
                       fnew_r = N(0,y,z)%n_r + df_r/df_scale
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b + df_b/df_scale
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s + df_s/df_scale
#endif
#endif
                    else
                       fnew_r = N(0,y,z)%n_r
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s
#endif
#endif
                    endif
                 endif

                 N(nx+1,y,z)%n_r = fnew_r
#ifndef SINGLEFLUID
                 N(nx+1,y,z)%n_b = fnew_b
#ifndef NOSURFACTANT
                 N(nx+1,y,z)%n_s = fnew_s
#endif
#endif
              enddo
           enddo
        else
           do k=1,nz+2
              do j=1,ny+2
                 N(nx+1,j-1,k-1)%n_r = minbuf(1:19,j,k)
#ifndef SINGLEFLUID
                 N(nx+1,j-1,k-1)%n_b = minbuf(20:38,j,k)
#endif
#ifndef NOSURFACTANT
                 N(nx+1,j-1,k-1)%n_s = minbuf(39:57,j,k)
                 N(nx+1,j-1,k-1)%da = minbuf(58:60,j,k)
                 N(nx+1,j-1,k-1)%db = minbuf(61:63,j,k)

#endif
                 N(nx+1,j-1,k-1)%rock_state = minbuf(ibuf,j,k)
              enddo
           enddo
        endif

#ifdef DEBUG_LE
	do s=1,nvecs
	   if (minval( N(nx+1,:,:)%n_r(s) ) .lt. 0 ) then
	      print *,'Ack! Negative density alert stage'
	   endif
	   if (minval( N(0,:,:)%n_r(s) ) .lt. 0 ) then
	      print *,'Ack! Negative density alert stage'
	   endif
	enddo
#endif

	if(allocated(pinbuf))deallocate(pinbuf)
	if(allocated(poutbuf))deallocate(poutbuf)
	if(allocated(minbuf))deallocate(minbuf)
	if(allocated(moutbuf))deallocate(moutbuf)

!!!!
 ! X swap done.
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Swap in the Y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(ibuf,nx+2,nz+2))
	allocate(poutbuf(ibuf,nx+2,nz+2))
	allocate(moutbuf(ibuf,nx+2,nz+2))
	allocate(minbuf(ibuf,nx+2,nz+2))

!!!!!!
 ! Begin asynchronous receives
!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,1),	& ! source
             tag_py,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pin),	& ! request
             ierror)
        bla = 'Bad +y async receive'
	call checkmpi(ierror,bla)

	call MPI_Irecv(		minbuf,		& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,2),	& ! source
             tag_my,		& ! tag
             Comm_Cart,	& ! communicator
             requests(min),	& ! request
             ierror)
        bla = 'Bad -y async receive'
	call checkmpi(ierror,bla)

        do k=1,nz+2
           do i=1,nx+2
              poutbuf(1:19,i,k) = N(i-1,ny,k-1)%n_r
#ifndef SINGLEFLUID
              poutbuf(20:38,i,k) = N(i-1,ny,k-1)%n_b
#endif
#ifndef NOSURFACTANT
              poutbuf(39:57,i,k) = N(i-1,ny,k-1)%n_s
              poutbuf(58:60,i,k) = N(i-1,ny,k-1)%da
              poutbuf(61:63,i,k) = N(i-1,ny,k-1)%db

#endif
              poutbuf(ibuf,i,k) = N(i-1,ny,k-1)%rock_state

              moutbuf(1:19,i,k) = N(i-1,1,k-1)%n_r
#ifndef SINGLEFLUID
              moutbuf(20:38,i,k) = N(i-1,1,k-1)%n_b
#endif
#ifndef NOSURFACTANT
              moutbuf(39:57,i,k) = N(i-1,1,k-1)%n_s
              moutbuf(58:60,i,k) = N(i-1,1,k-1)%da
              moutbuf(61:63,i,k) = N(i-1,1,k-1)%db

#endif
              moutbuf(ibuf,i,k) = N(i-1,1,k-1)%rock_state
           enddo
        enddo

!!!!!!
	! Begin asynchronous sends
!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,2),	& ! dest
             tag_py,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pout),	& ! request
             ierror)
        bla = 'Bad +y async send'
	call checkmpi(ierror,bla)

	call MPI_ISend(		moutbuf,	& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,1),	& ! dest
             tag_my,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mout),	& ! request
             ierror)
        bla = 'Bad -y async send'
	call checkmpi(ierror,bla)

!!!!!!
 ! Now wait for all I/O to complete
!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
 !	print*,'Y statuses:',statuses
        bla = 'MPI_Waitall() failed in Y direction'
	call checkmpi(ierror,bla)

        do k=1,nz+2
           do i=1,nx+2
              N(i-1,ny+1,k-1)%n_r = minbuf(1:19,i,k)
#ifndef SINGLEFLUID
              N(i-1,ny+1,k-1)%n_b = minbuf(20:38,i,k)
#endif
#ifndef NOSURFACTANT
              N(i-1,ny+1,k-1)%n_s = minbuf(39:57,i,k)
              N(i-1,ny+1,k-1)%da= minbuf(58:60,i,k)
              N(i-1,ny+1,k-1)%db = minbuf(61:63,i,k)

#endif
              N(i-1,ny+1,k-1)%rock_state = minbuf(ibuf,i,k)

              N(i-1,0,k-1)%n_r = pinbuf(1:19,i,k)
#ifndef SINGLEFLUID
              N(i-1,0,k-1)%n_b = pinbuf(20:38,i,k)
#endif
#ifndef NOSURFACTANT
              N(i-1,0,k-1)%n_s = pinbuf(39:57,i,k)
              N(i-1,0,k-1)%da= pinbuf(58:60,i,k)
              N(i-1,0,k-1)%db = pinbuf(61:63,i,k)

#endif
              N(i-1,0,k-1)%rock_state = pinbuf(ibuf,i,k)
           enddo
        enddo

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

 !print*,'Rank ',myrankc,'completed Y exchange.'

!!!!
 ! Y swap done.
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Swap in the Z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(ibuf,nx+2,ny+2))
	allocate(poutbuf(ibuf,nx+2,ny+2))
	allocate(moutbuf(ibuf,nx+2,ny+2))
	allocate(minbuf(ibuf,nx+2,ny+2))

!!!!!!
 ! Begin asynchronous receives
!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,1),	& ! source
             tag_pz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pin),	& ! request
             ierror)
        bla = 'Bad +z async receive'
	call checkmpi(ierror,bla)

	call MPI_Irecv(		minbuf,		& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,2),	& ! source
             tag_mz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(min),	& ! request
             ierror)
        bla = 'Bad -z async receive'
	call checkmpi(ierror,bla)

        do j=1,ny+2
           do i=1,nx+2
              poutbuf(1:19,i,j) = N(i-1,j-1,nz)%n_r
#ifndef SINGLEFLUID
              poutbuf(20:38,i,j) = N(i-1,j-1,nz)%n_b
#endif
#ifndef NOSURFACTANT
              poutbuf(39:57,i,j) = N(i-1,j-1,nz)%n_s
              poutbuf(58:60,i,j) = N(i-1,j-1,nz)%da
              poutbuf(61:63,i,j) = N(i-1,j-1,nz)%db
#endif
              poutbuf(ibuf,i,j) = N(i-1,j-1,nz)%rock_state

              moutbuf(1:19,i,j) = N(i-1,j-1,1)%n_r
#ifndef SINGLEFLUID
              moutbuf(20:38,i,j) = N(i-1,j-1,1)%n_b
#endif
#ifndef NOSURFACTANT
              moutbuf(39:57,i,j) = N(i-1,j-1,1)%n_s
              moutbuf(58:60,i,j) = N(i-1,j-1,1)%da
              moutbuf(61:63,i,j) = N(i-1,j-1,1)%db

#endif
              moutbuf(ibuf,i,j) = N(i-1,j-1,1)%rock_state
           enddo
        enddo

!!!!!!
	! Begin asynchronous sends
!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,2),	& ! dest
             tag_pz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pout),	& ! request
             ierror)
        bla = 'Bad +z async send'
	call checkmpi(ierror,bla)

	call MPI_ISend(		moutbuf,	& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,1),	& ! dest
             tag_mz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mout),	& ! request
             ierror)
        bla = 'Bad -z async send'
	call checkmpi(ierror,bla)

!!!!!!
 ! Now wait for all I/O to complete
!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
 !	print*,'Z statuses:',statuses
        bla = 'MPI_Waitall() failed in Z direction'
	call checkmpi(ierror,bla)

        do j=1,ny+2
           do i=1,nx+2
              N(i-1,j-1,nz+1)%n_r = minbuf(1:19,i,j)
#ifndef SINGLEFLUID
              N(i-1,j-1,nz+1)%n_b = minbuf(20:38,i,j)
#endif
#ifndef NOSURFACTANT
              N(i-1,j-1,nz+1)%n_s = minbuf(39:57,i,j)
              N(i-1,j-1,nz+1)%da  = minbuf(58:60,i,j)
              N(i-1,j-1,nz+1)%db   = minbuf(61:63,i,j)

#endif
              N(i-1,j-1,nz+1)%rock_state = minbuf(ibuf,i,j)

              N(i-1,j-1,0)%n_r = pinbuf(1:19,i,j)
#ifndef SINGLEFLUID
              N(i-1,j-1,0)%n_b = pinbuf(20:38,i,j)
#endif
#ifndef NOSURFACTANT
              N(i-1,j-1,0)%n_s = pinbuf(39:57,i,j)
              N(i-1,j-1,0)%da  = pinbuf(58:60,i,j)
              N(i-1,j-1,0)%db   = pinbuf(61:63,i,j)

#endif
              N(i-1,j-1,0)%rock_state = pinbuf(ibuf,i,j)
           enddo
        enddo

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

 !print*,'Rank ',myrankc,'completed Z exchange.'

!!!!
 ! Z swap done.
!!!!
    end subroutine le_halo_exchange


!>Holds any code relating to Lees Edwards that needs running
!>at the end of the simulation.
subroutine le_cleanup

      implicit none

      if (ccoords(1) .eq. 0) then
         deallocate(topprocs)
      endif

      if ((ccoords(1)+1) .eq. cdims(1)) then
         deallocate(botprocs)
      endif

end subroutine le_cleanup

end module lb3d_bc_leesedwards_module
