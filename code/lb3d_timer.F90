
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


!> easy to use timer module for code profiling

!> The first timer created has a special meaning: It should be started
!> and stopped in a way that the total runtime is detected, including
!> the time measured by all other timers. It is used to determine the
!> count for the timer "Other" that is automatically added at the end
!> of the list of timers and contains the time spent in all unassessed
!> parts of the code.

!> \par
!> Please note that the plain time counts are averaged over all
!> processes but the avg/min/max and histogram output might give a
!> rough idea of load balance.
module lb3d_timer_module
  use lb3d_global_module, only: rk, timesteps_count,start_time_t, end_time_t
  use lb3d_mpi_module, only:check_allocate,stats_rk
  use lb3d_mpi_parameters_module
  use lb3d_log_module
  
  implicit none
 ! include 'mpif.h'

  !> \name indices referencing timers
  !> \{
  integer,save :: ti_total,ti_adv,ti_intf,ti_intc,ti_halo,ti_inv,ti_dump
  !> \}

  !> specifies one timer
  type lbe_md_timer
     real(kind=rk) start      !< start time of the current counting period
     real(kind=rk) total      !< cumulative timer count
     character(len=16) name   !< short human-readable description
  end type lbe_md_timer

  !> contains all timers
  type(lbe_md_timer),save,allocatable,dimension(:) :: timers

contains
  !> register timers used in main LB code without MD
  subroutine lb3d_timer_init (stage)
    integer, intent(in) :: stage
#ifdef LB3D_DEBUG_INFO    
    call log_msg('In lb3d_timer_init.',.false.)
#endif
    call register_timer('Total',ti_total)
    call register_timer('Advection',ti_adv)
    call register_timer('Int:Force',ti_intf)
    call register_timer('Int:Coll',ti_intc)
    call register_timer('HaloExchange',ti_halo)
    call register_timer('Invasion',ti_inv)
    call register_timer('Dump',ti_dump)
  end subroutine lb3d_timer_init

  !> Create a new timer with name \c name. Its id is returned in \c ti.
  subroutine register_timer(name,ti)
    character(len=*),intent(in) :: name
    integer,intent(out) :: ti
    integer stat
    type(lbe_md_timer),allocatable,dimension(:) :: tmp


    if (.not.allocated(timers)) then
       allocate (timers(0),stat=stat)
       call check_allocate(stat,'register_timer(): timers')
    end if

    allocate (tmp(size(timers)),stat=stat)
    call check_allocate(stat,'register_timer(): tmp')
    tmp(:) = timers(:)
    deallocate (timers)
    allocate (timers(size(tmp)+1),stat=stat)
    call check_allocate(stat,'register_timer(): timers')
    timers(1:size(tmp)) = tmp(:)
    deallocate (tmp)

    ti = size(timers)
    timers(ti)%total = 0.0
    timers(ti)%name = name
  end subroutine register_timer

  !> let timer \c ti start counting
  subroutine start_timer(ti)
    integer,intent(in) :: ti

    timers(ti)%start = mpi_wtime()
  end subroutine start_timer

  !> let timer \c ti stop counting
  subroutine stop_timer(ti)
    integer,intent(in) :: ti

    timers(ti)%total = timers(ti)%total + mpi_wtime() - timers(ti)%start
  end subroutine stop_timer

  !> Writes a summary of all timer data to all units specified in  units .
  subroutine sum_timer(units)
    integer,intent(in) :: units(:)
    integer i,j,nt,ierror,stat
    integer ihisto(10),ihistotmp(10)
    real(kind=rk),allocatable :: avg(:),aavg(:)
    real(kind=rk) tavg,ave,xmax,xmin

    ! add another timer that will hold the differenc between the "Total"
    ! timer (always the first) and the sum of all other timers
    call register_timer('Other',nt)

    timers(nt)%total = timers(1)%total - sum(timers(2:nt-1)%total)

    allocate(avg(nt),aavg(nt),stat=stat)
    call check_allocate(stat,'sum_timer(): avg,aavg')

    aavg(:) = timers(:)%total/nprocs
    call mpi_reduce(aavg,avg,size(timers),MPI_REAL8,MPI_SUM,0,comm_cart,&
         &ierror)

    tavg = avg(1)           ! average of "Total" timer

    if (myrankc==0) then
       do i=1,nt
          do j=1,size(units)
             write (unit=units(j),fmt='(A," time/%",F15.6,F13.4)') &
                  &timers(i)%name,avg(i),100*avg(i)/tavg
          end do
       end do

       do j=1,size(units)
          write (unit=units(j),fmt='()')
       end do
    end if

    do i=1,nt
       call stats_rk(timers(i)%total,ave,xmax,xmin,ihisto,ihistotmp,10)
       if (myrankc==0) then
          do j=1,size(units)
             write (unit=units(j),fmt=&
                  &'(A," time: ",F13.4," ave",F13.4," max",F13.4," min")')&
                  &timers(i)%name,ave,xmax,xmin
             write (unit=units(j),fmt='("  Histogram: ",10(I6,:,X))') ihisto
          end do
       end if
    end do
  end subroutine sum_timer

!> Sets the \c inp_file variable to contain the name of the input file.
!> If the -f option is passed on the command line, then
!> this is used. Otherwise, .input-file is checked -- if it contains the
!> word "INTERACTIVE", then the input file name is read from stdin.
!> Otherwise, the file named in .input-file is read.
subroutine calculate_performance()
  integer*8              :: nsites, nupdates
  integer                :: delta_t
  real(kind=rk)                 :: updaterate, ratepercpu
  character(len=128)     :: msgstr


  delta_t = (end_time_t - start_time_t)
  nsites = int(tnx,kind=8)*int(tny,kind=8)*int(tnz,kind=8)
  nupdates = timesteps_count * nsites
  updaterate = real(nupdates) / real(delta_t)
  ratepercpu = updaterate / real(nprocs)

  write(msgstr,"('Updated ',i0,' lattice sites for ',i0,' time steps -> ',i0,' site updates.')") &
    nsites, timesteps_count, nupdates
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('This took ',i0,' seconds ->')") delta_t
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('    ',F16.2,' updates per second')") updaterate
  CALL log_msg(trim(msgstr),.false.)
  write(msgstr,"('    ',F16.2,' updates per second per CPU')") ratepercpu
  CALL log_msg(trim(msgstr),.false.)

end subroutine calculate_performance


end module lb3d_timer_module
