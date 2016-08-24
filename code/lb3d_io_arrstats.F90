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

!> Dumping of velocity field statistics
module lb3d_io_arrstats_module

    use lb3d_global_module
    use lb3d_lattice_module, only: lbe_site,halo_extent
    use lb3d_config_module, only: arrstats_intervals,arrstats_name&
         &,n_sci_arrstats_dump,nt,nx,ny,nz,sci_arrstats

    use lb3d_log_module
    use lb3d_mpi_module

    use lb3d_analysis_module, only: velocity

    use lb3d_io_helper_module, only: lbe_make_filename_restore,every_n_time_steps
    use lb3d_io_dump_data_module, only: dump_iscalar,dump_vector

#ifdef USEHDF
    use lb3d_io_hdf5_module, only: read_iscalar_phdf5,read_vector_phdf5
#endif

#ifdef USEXDRF
    use lb3d_io_xdrf_module, only: check_xdrfopen
#endif

    implicit none
    private
    public dump_arrstats,dump_arrstats_checkpoint,lb3d_io_setup_arrstats&
         &,restore_arrstats_checkpoint

    !> type representing data accumulated at one sampling interval
    type arrstats_data_type
       integer :: n_samples     !< number of sampled time steps
       !> number of samples per lattice site
       integer,allocatable,dimension(:,:,:) :: n
       !> accumulation buffer, contains 9 elements per lattice, the accumulated
       !> \c v(1:3), \c v(1:3)*v(1:3), and the quadratic cross terms \c v(1)*v(
       real(kind=rk),allocatable,dimension(:,:,:,:) :: m
    end type arrstats_data_type

    !> data accumulated at different sampling intervals
    type(arrstats_data_type),save :: arrstats_data(arrstats_intervals)

contains
    !> accumulates data at all sampling intervals, eventually dumps them
    !>
    !> \param[in] N local lattice chunk with halo of size 1
    subroutine dump_arrstats(N)
        type(lbe_site),intent(in) :: N(0:,0:,0:)
        integer i

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call sample_arrstats(N,i,arrstats_data(i))

           if (every_n_time_steps(n_sci_arrstats_dump(i))) then
              call log_msg(" -now dumping...",.false.)
              call dump_arrstats_data(i,arrstats_data(i))
              call reset_arrstats_data(arrstats_data(i))
           end if

        end do
    end subroutine dump_arrstats

    !> initialize \c arrstats output
    subroutine lb3d_io_setup_arrstats
        integer i

        if (.not.sci_arrstats) return

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call init_arrstats_data(i,arrstats_data(i))
        end do
    end subroutine lb3d_io_setup_arrstats

    !> initialize data objects related to a specific sampling interval
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine init_arrstats_data(i,asd)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        integer :: stat

        if (trim(arrstats_name(i))=='') &
             &write (arrstats_name(i),fmt='("d",I0)') n_sci_arrstats_dump(i)

        allocate (asd%n(1:nx,1:ny,1:nz),asd%m(1:9,1:nx,1:ny,1:nz),stat=stat)
        call check_allocate(stat&
             &,'init_arrstats(): asd%n,asd%m (arrstats_name <'&
             &//trim(arrstats_name(i))//'>)')

        call reset_arrstats_data(asd)
    end subroutine init_arrstats_data

    !> samples the current velocity field into one accumulation buffer
    !>
    !> \param[in] N local lattice chunk (halo extent 1)
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine sample_arrstats(N,i,asd)
        type(lbe_site),intent(in) :: N(0:,0:,0:)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        integer :: x,y,z
        real(kind=rk) :: v(3)

        call log_msg(' -"'//trim(arrstats_name(i))//'"',.false.)

        do x=1,nx
           do y=1,ny
              do z=1,nz
                 if (N(x,y,z)%rock_state==0.0_rk) then
                    asd%n(x,y,z) = asd%n(x,y,z)+1
                    v = velocity(N(x,y,z))
                    asd%m(1:3,x,y,z) = asd%m(1:3,x,y,z)+v ! velocity
                    asd%m(4:6,x,y,z) = asd%m(4:6,x,y,z)+v*v ! square terms
                    asd%m(7,x,y,z) = asd%m(7,x,y,z)+v(2)*v(1) ! cross terms...
                    asd%m(8,x,y,z) = asd%m(8,x,y,z)+v(3)*v(1)
                    asd%m(9,x,y,z) = asd%m(9,x,y,z)+v(3)*v(2)
                 end if
              end do
           end do
        end do
        asd%n_samples = asd%n_samples+1
    end subroutine sample_arrstats

    !> Writes the complete data related to one sampling interval to disk
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    !>
    !> \param[in] pprefix prefix to add to the filenames (optional)
    subroutine dump_arrstats_data(i,asd,pprefix)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        character(len=*),intent(in),optional :: pprefix
        character(len=32) :: prefix

        if (present(pprefix)) then
           prefix = pprefix
        else
           prefix = ''
        end if

        call dump_arrstats_desc(i,asd,prefix)
        call dump_iscalar(asd%n&
             &,trim(prefix)//'arrstats-n-'//trim(arrstats_name(i)))
        call dump_vector(asd%m&
             &,trim(prefix)//'arrstats-m-'//trim(arrstats_name(i)))
    end subroutine dump_arrstats_data

    !> Writes additional descriptive data regarding a state of
    !> accumulated data to disk
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in] asd corresponding accumulation buffer
    !>
    !> \param[in] pprefix prefix to add to the filenames (can be \c '')
    !>
    !> At the moment the file format is always XDR. The data written
    !> encompasses the number of sampled time steps and the sampling
    !> interval. From the first, the time interval covered by the data
    !> can be obtained if \c n_sci_arrstats is known.
    subroutine dump_arrstats_desc(i,asd,prefix)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(in) :: asd
        character(len=32),intent(in) :: prefix
        character(len=256) :: filename
        integer :: file_id,ierror

        if (myrankc/=0) return

        call lbe_make_filename_output(filename&
             &,trim(prefix)//'arrstats-desc-'//trim(arrstats_name(i)),'.xdr',nt)
#ifndef USEXDRF
        call error('dump_arrstats_desc() requires still XDRF')
#else
        call xdrfopen(file_id,filename,"w",ierror)
        call check_xdrfopen(ierror,filename)

        call xdrfint(file_id,asd%n_samples,ierror)
        call xdrfint(file_id,n_sci_arrstats_dump(i),ierror)

        call xdrfclose(file_id,ierror)
#endif
    end subroutine dump_arrstats_desc

    !> Restores the state accumulated data for one sampling interval
    !> from a checkpoint
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine restore_arrstats_desc(i,asd)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        character(len=256) :: filename,msgstr
        integer :: dummy_n_sci_arrstats_dump,file_id,ierror

        if (myrankc/=0) return

        call lbe_make_filename_restore(filename&
             &,'checkarrstats-desc-'//trim(arrstats_name(i)),'.xdr')
#ifndef USEXDRF
        call error('restore_arrstats_desc() requires still XDRF')
#else
        call xdrfopen(file_id,filename,'r',ierror)
        call check_xdrfopen(ierror,filename)

        call xdrfint(file_id,asd%n_samples,ierror)
        call xdrfint(file_id,dummy_n_sci_arrstats_dump,ierror)

        call xdrfclose(file_id,ierror)
#endif

        if (dummy_n_sci_arrstats_dump/=n_sci_arrstats_dump(i)) then
           write (unit=msgstr,fmt='("WARNING: read n_sci_arrstats_dump==",I0,'&
                &//'" from <",A,"> (which is ignored) but the input files '&
                &//'set n_sci_arrstats_dump==",I0," already.")') &
                &dummy_n_sci_arrstats_dump,filename,n_sci_arrstats_dump(i)
           call log_msg(trim(msgstr),.false.)
        end if
    end subroutine restore_arrstats_desc

    !> Resets an instance of \c arrstats_data_type for further accumulation
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine reset_arrstats_data(asd)
        type(arrstats_data_type),intent(inout) :: asd

        asd%n_samples = 0
        asd%n = 0
        asd%m = 0.0_rk
    end subroutine reset_arrstats_data

    !> Writes checkpoints for all active sampling intervals
    subroutine dump_arrstats_checkpoint
        integer i,n_wrote
        character(len=256) :: msgstr

        if (.not.sci_arrstats) return

        n_wrote = 0
        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle
           call dump_arrstats_data(i,arrstats_data(i),'check')
           n_wrote = n_wrote+1
        end do

        write (unit=msgstr,fmt='("  Wrote ",I0," arrstats checkpoints.")') &
             &n_wrote
        call log_msg(trim(msgstr),.false.)
    end subroutine dump_arrstats_checkpoint

    !> Restore all accumulated data from a checkpoint
    subroutine restore_arrstats_checkpoint
        character(len=256) :: filename,msgstr
        logical :: exists
        integer :: i

        if (.not.sci_arrstats) return

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call lbe_make_filename_restore(filename&
                &,'checkarrstats-desc-'//trim(arrstats_name(i)),'.xdr')
           inquire (file=filename,exist=exists)
           file_exists: if (exists) then
              call restore_arrstats_desc(i,arrstats_data(i))

              call lbe_make_filename_restore(filename&
                   &,'checkarrstats-n-'//trim(arrstats_name(i)),'.h5')
#ifdef USEHDF
              call read_iscalar_phdf5(arrstats_data(i)%n,filename)
#else
              call log_msg("Need HDF5 to restore arrstats",.false.)
              call Abend
#endif

              call lbe_make_filename_restore(filename&
                   &,'checkarrstats-m-'//trim(arrstats_name(i)),'.h5')
#ifdef USEHDF
              call read_vector_phdf5(arrstats_data(i)%m,filename&
                   &,size(arrstats_data(i)%m,1))
#else
              call log_msg("Need HDF5 to restore arrstats",.false.)
              call Abend
#endif

              write (unit=msgstr&
                   &,fmt='("Restored ",I0," samples for arrstats <",A,">.")') &
                   &arrstats_data(i)%n_samples,trim(arrstats_name(i))
              call log_msg(trim(msgstr),.false.)
           else file_exists
              call log_msg('arrstats checkpoint file <'//trim(filename)&
                   &//'> not found---resetting arrstats <'&
                   &//trim(arrstats_name(i))//'> to zero!',.false.)
              call reset_arrstats_data(arrstats_data(i))
           end if file_exists
        end do
    end subroutine restore_arrstats_checkpoint

end module lb3d_io_arrstats_module
