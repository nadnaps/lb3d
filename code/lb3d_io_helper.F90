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

!> helper routines for IO
module lb3d_io_helper_module

  use lb3d_mpi_module, only:abend
  use lb3d_config_module!, only: chk_uid,dump_double,dump_format,folder,cpfolder&
       !&,gr_out_file,nt,restore_string,srccpfolder
  implicit none

  !> Number of flags set - counted in lbe_detect_flags() and used for
  !> HDF metadata
  integer :: nflags

  !> \{
  !> \name not in namelist
  character(len=32), save :: restore_fmt_t = '("t",i8.8,"-",i10.10)'
  character(len=32), save :: restore_fmt_p = '("p",i8.8,"-",i10.10)'
  !> \}

contains

!> Wrapper for file deletion.
subroutine lbe_delete_file(filename)
  implicit none
  character(len=*)   :: filename
  character(len=512) :: msgstr, syscmd

! CALL system() not available on Blue Gene/P
#ifdef NOCALLSYSTEM
  write(msgstr,"('Calling truncate <',A,'>.')") trim(filename)
  DEBUG_CHECKPOINT_MSG_ALL(trim(msgstr))
  CALL truncate(trim(filename)//'\0',0)
  write(msgstr,"('Calling unlink <',A,'>.')") trim(filename)
  DEBUG_CHECKPOINT_MSG_ALL(trim(msgstr))
  CALL unlink(trim(filename)//'\0')
#else
  syscmd='rm -f '//trim(filename)
  write(msgstr,"('Executing command <',A,'>.')") trim(syscmd)
  DEBUG_CHECKPOINT_MSG_ALL(trim(msgstr))
  CALL system(trim(syscmd))
#endif

end subroutine lbe_delete_file

!> Makes a filename by appending the given suffix and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_output(buffer, prefix, suffix, time)
  implicit none

  character(len=*)   :: buffer, prefix, suffix
  integer            :: time
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) time, chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'_',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine lbe_make_filename_output

!> Makes a filename by appending the given suffix, plus CPU rank and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_cp_rank(buffer, prefix, suffix, time, rank)
  character(len=*)   :: buffer, prefix, suffix
  integer            :: time, rank
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) time, chk_uid

  write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
    trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(chkuidstr), rank, trim(suffix)
end subroutine lbe_make_filename_cp_rank

!> Makes a filename by appending the given suffix, plus CPU rank and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_restore_rank(buffer, prefix, suffix, rank)
  implicit none

  character(len=*)   :: buffer, prefix, suffix
  integer            :: rank

  if (trim(srccpfolder) .ne. "") then
    write(buffer,"(A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
      trim(srccpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), rank, trim(suffix)
  else
    write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
      trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), rank, trim(suffix)
  endif

end subroutine lbe_make_filename_restore_rank

!> Makes a filename by appending the given suffix and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_cp(buffer, prefix, suffix, time)
  character(len=*)   :: buffer, prefix, suffix
  integer            :: time
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) time, chk_uid

  write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,A)") & 
    trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)
end subroutine lbe_make_filename_cp


subroutine lbe_make_filename_restore(buffer, prefix, suffix)
  implicit none

  character(len=*)   :: buffer, prefix, suffix

  if (trim(srccpfolder) .ne. "") then
    write(buffer,"(A,'/',A,'_',A,'_',A,A)") & 
      trim(srccpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), trim(suffix)
  else
    write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,A)") & 
      trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), trim(suffix)
  endif
end subroutine lbe_make_filename_restore

subroutine lbe_make_filename_append(buffer, prefix, suffix)
  implicit none

  character(len=*)   :: buffer, prefix, suffix
  character(len=10)  :: chkuidstr

  write(chkuidstr,"(I0)") chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'-',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine lbe_make_filename_append

!>  This routine writes an AVS field file
!>  Variable veclen determines fiedltype:
!>  1=scalar,2=2scalar,3=3scalar,4=vector
!> 
!>  Added 21.06.02 by Jens
subroutine dump_avs_fld(prefix,nxi,nyi,nzi,vectmp)
  implicit none
  character(LEN=256) :: filename,fldname
  character(LEN=256) :: tmpstring
  character(LEN=*) :: prefix
  character(LEN=1)   :: vecstr,countstr,stridestr
  character(LEN=7)   :: skipstr
  integer :: nxi,nyi,nzi,veclen,i,i2,i3,vectmp

  ! Stupid compiler bug workaround on IRIX
  veclen = vectmp

  call lbe_make_filename_output(fldname,prefix,'.fld',nt)
  open(10,file=fldname)
  write(10,'(a)') '# AVS field file'
  write(10,'(a)') 'ndim=3'
  write(10,'(a,i3)') 'dim1=',nxi
  write(10,'(a,i3)') 'dim2=',nyi
  write(10,'(a,i3)') 'dim3=',nzi
  write(10,'(a)') 'nspace=3'
  write(10,'(a)') 'field=uniform'

  if (index(dump_format,'all').gt.0) then
    call lbe_make_filename_output(filename,prefix,'.all',nt)
    if (veclen.eq.4) veclen=3
    write(stridestr,'(i1.1)') veclen+3
    write(10,'(a)') 'data=float'
    write(vecstr,'(i1.1)') veclen
    write(10,'(a)') 'veclen='//vecstr
      do i=1,veclen
        write(countstr,'(i1.1)') i 
        write(skipstr,'(i3.3)') i+2
        tmpstring='variable '//trim(countstr)//' file='//      &
            trim(filename)//' filetype=ascii skip=0 offset='// &
            trim(skipstr)//' stride='//trim(stridestr)
        write(10,'(a,a,a)') trim(tmpstring)
      end do
  else
    if (index(dump_format,'xdr').gt.0) then
      if (veclen.eq.4) veclen=3
      call lbe_make_filename_output(filename,prefix,'.xdr',nt)
      write(vecstr,'(i1.1)') veclen
      write(10,'(a)') 'veclen='//vecstr
      if (dump_double) then
        write(10,'(a)') 'data=xdr_double'
      else
        write(10,'(a)') 'data=xdr_float'
      end if
      do i=1,veclen
        write(countstr,'(i1.1)') i
        write(skipstr,'(i7.7)') 8*(i-1)
        tmpstring='variable '//countstr//' file='//trim(filename)// &
          ' filetype=binary skip='//skipstr//' stride='//vecstr
        write(10,'(a,a,a)') trim(tmpstring)
      end do
    else
      call lbe_make_filename_output(filename,prefix,'.bin',nt)

      if (dump_double) then
        write(10,'(a)') 'data=double'
        i3 = 8
      else
        write(10,'(a)') 'data=float'
        i3 = 4
      end if

      i2 = veclen
      if (veclen.eq.4) veclen=3
      write(vecstr,'(i1.1)') veclen
      write(10,'(a)') 'veclen='//vecstr

      do i=1,veclen
        write(countstr,'(i1.1)') i
          if (i2.lt.4) then
          write(skipstr,'(i7.7)') i3*nxi*nyi*nzi*(i-1)
          stridestr='1'
          else
          write(skipstr,'(i7.7)') i3*(i-1)
          write(stridestr,'(i1.1)') veclen
          end if
        tmpstring='variable '//countstr//' file='//trim(filename)// &
        ' filetype=unformatted skip='//skipstr//&
        ' stride='//trim(stridestr)
        write(10,'(a,a,a)') trim(tmpstring)
      end do
    endif
  endif

  close(10)
end subroutine dump_avs_fld

!> initialize default names for fluxz regions
subroutine lbe_setup_fluxz()
  integer :: i

#ifdef MD
  provide_uid2i = .true.
#endif

  do i=1,fluxz_regions
    write (fluxz_name(i),fmt='(I0)') i
  end do
end subroutine lbe_setup_fluxz

!>Checks to see if the given file exists. If not,
!>prints its name and an error message, and aborts.
!>Otherwise, silently returns.
subroutine die_unless_exists(filename)
  character(len=128) :: filename
  logical :: exists_p
  inquire(file=filename ,exist = exists_p)
  if (.not. exists_p) then
    print*,'File <',trim(filename),'> not found. Aborting...'
    call Abend
  end if
  print*,'File <',trim(filename),'> found OK.'
end subroutine die_unless_exists

subroutine die_unless_nonzero(ierr,string)
  integer :: ierr
  character(len=*),intent(in) :: string

  if (ierr .eq. 0) then
    print*,string
    call abend
  end if
end subroutine die_unless_nonzero

!> This function makes sure that the boolean is checked before the mod.
!>
!> \param[in] sci disables or enables check of \c n_sci
!>
!> \param[in] n_sci if checked, returns \c .true. every \c n_sci time steps
!>
!> \returns \c .false. if \c .not.sci, \c every_n_time_steps(n_sci) otherwise
!>
!> 2010-05-06. Added by Stefan. Putting this in a function avoids
!> nasty crashes, and keeps the nested if-loops confined to one place.
logical function check_dump_now(sci, n_sci)
  logical,intent(in) :: sci
  integer,intent(in) :: n_sci

  if ( sci ) then
     check_dump_now = every_n_time_steps(n_sci)
  else
     check_dump_now = .false.
  endif
end function check_dump_now

!> checks whether the n-ths time step is reached while safely
!> returning \c .false. if \c n==0
!>
!> \param[in] n interval against which to check \c nt
!>
!> \returns \c .false. if \c n==0, \c mod(nt.n)==0 otherwise
!>
!> This can be used in querying md-style time interval parameters
!> where a value of 0 is valid and means that the respective feature
!> is disabled.
logical function every_n_time_steps(n)
    integer,intent(in) :: n

    if (n==0) then
       every_n_time_steps = .false.
    else
       every_n_time_steps = mod(nt,n)==0
    end if
end function every_n_time_steps


end module lb3d_io_helper_module
