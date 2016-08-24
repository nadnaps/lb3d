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

!      ==================================================================
module lb3d_force_module
!      ==================================================================

  use lb3d_global_module
  use lb3d_config_module, only: force,nx,ny,nz
  use lb3d_log_module

  use lb3d_lattice_module, only: lbe_site,halo_extent,force_halo_extent
  use lb3d_mpi_module!, only: check_allocate,comm_cart,error,nnprocs
  
  use lb3d_force_constant_module, only: force_apply_constant&
       &,force_init_constant,force_input_constant
  use lb3d_force_kolmogorov_module, only: force_apply_kolmogorov&
       &,force_init_kolmogorov,force_input_kolmogorov
 
  implicit none

  public lbe_add_force_halo,lbe_force_apply,lb3d_force_init,lbe_force_input

 !> \name buffers to recv halo forces into
    !> \{
    !> indices: component, species, x, y, z
    real*8,save,allocatable,dimension(:,:,:,:,:) :: x_rbuf,y_rbuf,z_rbuf
    !> \}

    !> mpi datatype representing the whole \c lbe_force data at a lattice
    !> point (in fortran something like \c real*8,dimension(3,n_spec) )
    integer,save :: forcedata_mpitype

    !> mpi datatypes representing \c [xyz]_rbuf .
    integer,save :: r_mpitype(3)

    !> \name more custom MPI data types
    !> \{
    !> MPI datatypes representing the lower (\c l)/upper (\c u) part
    !> of the force halo that is sent (\c s) during the exchange for
    !> each direction (indices 1/2/3 represent x/y/z). Datatype
    !> definitions are relative to the memory address of \c lbe_force
    !> as a whole.
    integer,save :: ls_f_mpitype(3) !< for lower halos
    integer,save :: us_f_mpitype(3) !< for upper halos
    !> \}


contains

!      ==================================================================
  subroutine lb3d_force_init (stage)
!      ==================================================================

    integer :: stage

#ifdef LB3D_DEBUG_INFO
    write(msgstr,"('In lb3d_force_init stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif

    select case (stage)
    case (0) ! pre mpi init
       
       call lbe_force_input()

    case (4) ! post mem alloc

       call lbe_force_init()

    case default
#ifdef LB3D_DEBUG_INFO
       call log_msg('nothing to do.',.false.) 
#endif
    end select    
  
  end subroutine lb3d_force_init

  subroutine lb3d_force_apply (stage)
!==================================================================

    integer :: stage

#ifdef LB3D_DEBUG_INFO
    write(msgstr,"('In lb3d_force_apply stage',I0)") stage
    call log_msg(trim(msgstr),.false.) 
#endif

    select case (stage)

    case (6) ! pre advection

    case (8) ! post advection / pre collision

    case (10) ! post collision / pre advection / default dump

    case default
       !FIXME Error

    end select

  end subroutine lb3d_force_apply


    !> Builds a custom mpi datatype that represents the whole force
    !> data on a lattice point (in fortran: \c
    !> real(kind=rk),dimension(3,n_spec) ) and stores it in \c fdmt.
    subroutine build_forcedata_mpitype(fdmt)
        integer,intent(out) :: fdmt ! type to be built
        real(kind=rk) :: sample(3,n_spec) ! sample of a force lattice point
        integer s_mt ! temporary mpi data type
        integer stride,ierror
        integer(kind=MPI_ADDRESS_KIND) adr1,adr2

        ! mpi datatype for all vector components of a species ( sample(:,s) )
        call mpi_get_address(sample(1,1),adr1,ierror)
        call mpi_get_address(sample(2,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_hvector(3,1,stride,MPI_REAL8,s_mt,ierror)

        ! mpi datatype for whole point data ( sample(:,:) )
#ifdef SINGLEFLUID
        ! for  SINGLEFLUID ,  n_spec==1 , so the type is already complete
        ! (in fact the #else branch would not compile on JUMP because the
        ! compiler complains "Subscript is out of bounds." in  sample(1,2) .
        fdmt = s_mt
#else
        ! for not  SINGLEFLUID , we need  n_spec  (2 or 3) times  s_mt :
        call mpi_get_address(sample(1,1),adr1,ierror)
        call mpi_get_address(sample(1,2),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_hvector(n_spec,1,stride,s_mt,fdmt,ierror)
#endif

        call mpi_type_commit(fdmt,ierror)
    end subroutine build_forcedata_mpitype

    !> Builds a custom mpi datatype that represents the part of \c f
    !> specified by the coordinate intervals \c xr(2),yr(2),zr(2) and
    !> stores it in \c fcmt . The starting indices of f are \c si(:).
    !>
    !> \param[in] si starting indices of \c f
    !>
    !> \param[in] f sample force lattice
    !>
    !> \param[in] xr chunk range in x-direction
    !>
    !> \param[in] yr chunk range in y-direction
    !>
    !> \param[in] zr chunk range in z-direction
    !>
    !> \param[out] fmct type to be built
    subroutine build_force_chunk_mpitype(si,f,xr,yr,zr,fcmt)
        integer,intent(in) :: si(3)
        real(kind=rk),intent(in) :: f(1:,1:,si(1):,si(2):,si(3):)
        integer,intent(in) :: xr(2),yr(2),zr(2)
        integer,intent(out) :: fcmt
        integer xrow_mt,xyplane_mt,xyzchunk_mt ! temporary mpi data types
        integer blocks,stride,ierror
        integer(kind=MPI_ADDRESS_KIND) adr1,adr2,base,offset
        integer lengths(1),displs(1) ! for  mpi_type_hindexed()

        ! mpi datatype for slices of  f  like  f(:,:,xr(1):xr(2),y,z)
        blocks = 1 + xr(2) - xr(1)
        call mpi_get_address(f(1,1,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,1,2,1,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_hvector(blocks,1,stride,forcedata_mpitype,xrow_mt,ierror)

        ! mpi datatype for slices of  f  like  f(:,:,xr(1):xr(2),yr(1):yr(2),z)
        blocks = 1 + yr(2) - yr(1)
        call mpi_get_address(f(1,1,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,1,1,2,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_hvector(blocks,1,stride,xrow_mt,xyplane_mt,ierror)

        ! mpi datatype for whole chunk
        ! f(:,:,xr(1):xr(2),yr(1):yr(2),zr(1):zr(2))
        blocks = 1 + zr(2) - zr(1)
        call mpi_get_address(f(1,1,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,1,1,1,2),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_hvector(blocks,1,stride,xyplane_mt,xyzchunk_mt,ierror)

        ! position of the beginning of the chunk relative to the beginning of  f
        call mpi_get_address(f(1,1,0,0,0),base,ierror)
        call mpi_get_address(f(1,1,xr(1),yr(1),zr(1)),offset,ierror)
        offset = offset - base

        ! fcmt becomes a datatype like  xyzchunk_mt  but relative to  base
        lengths = (/1/)
        displs = (/int(offset)/)
        call mpi_type_hindexed(1,lengths,displs,xyzchunk_mt,fcmt,ierror)
        call mpi_type_commit(fcmt,ierror)
    end subroutine build_force_chunk_mpitype

    !> initializes buffers and custom mpi types concerning \c lbe_force
    subroutine init_buffers_and_types
        integer stat
        integer :: h1,he

        h1=force_halo_extent-1
        he=force_halo_extent

        ! allocate and initialize arrays
        allocate (lbe_force(3,n_spec,-h1:nx+he,-h1:ny+he,-h1:nz+he)&
             &,x_rbuf(3,n_spec,     1:he,-h1:ny+he,-h1:nz+he)&
             &,y_rbuf(3,n_spec,-h1:nx+he,     1:he,-h1:nz+he)&
             &,z_rbuf(3,n_spec,-h1:nx+he,-h1:ny+he,     1:he),stat=stat)
        call check_allocate(stat&
             &,'init_buffers_and_types(): lbe_force,[xyz]_rbuf')
        lbe_force(:,:,:,:,:) = 0.0_rk

        ! create custom mpi types for communication of force halo

        call build_forcedata_mpitype(forcedata_mpitype)

        ! types used for sending
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,    0/),(/ -h1,ny+he/),(/ -h1,nz+he/),ls_f_mpitype(1))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/nx+1,nx+he/),(/ -h1,ny+he/),(/ -h1,nz+he/),us_f_mpitype(1))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,    0/),(/ -h1,nz+he/),ls_f_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ny+1,ny+he/),(/ -h1,nz+he/),us_f_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,ny+he/),(/ -h1,    0/),ls_f_mpitype(3))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,ny+he/),(/nz+1,nz+he/),us_f_mpitype(3))

        ! types for receiving
        call build_force_chunk_mpitype((/  1,-h1,-h1/),x_rbuf&
             &,(/  1,   he/),(/-h1,ny+he/),(/-h1,nz+he/),r_mpitype(1))
        call build_force_chunk_mpitype((/-h1,  1,-h1/),y_rbuf&
             &,(/-h1,nx+he/),(/  1,   he/),(/-h1,nz+he/),r_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,  1/),z_rbuf&
             &,(/-h1,nx+he/),(/-h1,ny+he/),(/  1,   he/),r_mpitype(3))
    end subroutine init_buffers_and_types

    !> Add halo forces from neighbor processes to own \c
    !> lbe_force. Send own force halo to neighbors. Reset force halo
    !> to zero.
    subroutine lbe_add_force_halo
        integer ierror,status(MPI_STATUS_SIZE)
        integer :: h1,he

        if (.not. use_lbe_force) return

        h1=force_halo_extent-1
        he=force_halo_extent

        ! send "downward" in x direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,ls_f_mpitype(1),nnprocs(1,1),0&
             &,x_rbuf(1,1,0,0,0),1,r_mpitype(1),nnprocs(1,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,nx-h1:nx,-h1:ny+he,-h1:nz+he)&
             & = lbe_force(:,:,nx-h1:nx,-h1:ny+he,-h1:nz+he)&
             & +    x_rbuf(:,:,    1:he,-h1:ny+he,-h1:nz+he)

        ! send "upward" in x direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,us_f_mpitype(1),nnprocs(1,2),0&
             &,x_rbuf(1,1,0,0,0),1,r_mpitype(1),nnprocs(1,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,    1:he,-h1:ny+he,-h1:nz+he)&
             & = lbe_force(:,:,    1:he,-h1:ny+he,-h1:nz+he)&
             & +    x_rbuf(:,:,    1:he,-h1:ny+he,-h1:nz+he)

        ! send "downward" in y direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,ls_f_mpitype(2),nnprocs(2,1),0&
             &,y_rbuf(1,1,0,0,0),1,r_mpitype(2),nnprocs(2,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,ny-h1:ny,-h1:nz+he)&
             & = lbe_force(:,:,-h1:nx+he,ny-h1:ny,-h1:nz+he)&
             & +    y_rbuf(:,:,-h1:nx+he,    1:he,-h1:nz+he)

        ! send "upward" in y direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,us_f_mpitype(2),nnprocs(2,2),0&
             &,y_rbuf(1,1,0,0,0),1,r_mpitype(2),nnprocs(2,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,    1:he,-h1:nz+he)&
             & = lbe_force(:,:,-h1:nx+he,    1:he,-h1:nz+he)&
             & +    y_rbuf(:,:,-h1:nx+he,    1:he,-h1:nz+he)

        ! send "downward" in z direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,ls_f_mpitype(3),nnprocs(3,1),0&
             &,z_rbuf(1,1,0,0,0),1,r_mpitype(3),nnprocs(3,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,-h1:ny+he,nz-h1:nz)&
             & = lbe_force(:,:,-h1:nx+he,-h1:ny+he,nz-h1:nz)&
             & +    z_rbuf(:,:,-h1:nx+he,-h1:ny+he,    1:he)

        ! send "upward" in z direction
        call mpi_sendrecv&
             &(lbe_force(1,1,0,0,0),1,us_f_mpitype(3),nnprocs(3,2),0&
             &,z_rbuf(1,1,0,0,0),1,r_mpitype(3),nnprocs(3,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,-h1:ny+he,    1:he)&
             & = lbe_force(:,:,-h1:nx+he,-h1:ny+he,    1:he)&
             & +    z_rbuf(:,:,-h1:nx+he,-h1:ny+he,    1:he)

        ! clear halo
        lbe_force(:,:,     -h1:0, -h1:ny+he, -h1:nz+he) = 0.0_rk
        lbe_force(:,:,nx+1:nx+he, -h1:ny+he, -h1:nz+he) = 0.0_rk
        lbe_force(:,:, -h1:nx+he,     -h1:0, -h1:nz+he) = 0.0_rk
        lbe_force(:,:, -h1:nx+he,ny+1:ny+he, -h1:nz+he) = 0.0_rk
        lbe_force(:,:, -h1:nx+he, -h1:ny+he,     -h1:0) = 0.0_rk
        lbe_force(:,:, -h1:nx+he, -h1:ny+he,nz+1:nz+he) = 0.0_rk
    end subroutine lbe_add_force_halo

    !> apply force depending on the chosen implementation
    !>
    !> \param[in] lbe_N local lattice chunk with halo 1
    !>
    !> \param[in] whole_N local lattice chunk with full halo
    subroutine lbe_force_apply(lbe_N,whole_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (.not.use_lbe_force) return

        select case (force)
        case ('constant')
           call force_apply_constant(lbe_N,whole_N)
        case ('kolmogorov')
           call force_apply_kolmogorov
        case ('none')
        end select
    end subroutine lbe_force_apply

    !> initialize common force data structures and branch off into
    !> actual force implementations
    subroutine lbe_force_init
        select case (force)
        case ('constant')
           call force_init_constant
        case ('kolmogorov')
           call force_init_kolmogorov
        case ('none')
        end select

        if (use_lbe_force) call init_buffers_and_types
    end subroutine lbe_force_init

    !> read in force namelists and initialize things that are required
    !> before other initialization routines are called
    subroutine lbe_force_input
        select case (force)
        case ('constant')
           call force_input_constant
        case ('kolmogorov')
           call force_input_kolmogorov
        case ('none')
        case default
           call error('unknown type of force: force="'//force//'"')
        end select
    end subroutine lbe_force_input



end module lb3d_force_module
