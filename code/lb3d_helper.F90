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

!> General helper routines
module lb3d_helper_module

  implicit none
  
!  private

!  public check_dump_now,check_NaN,cross_product,density,every_n_time_steps&
!       &,interpolate,local_coordinates,makedatestr,massflow,norm&
!       &,report_check_NaN_function,unit_vector,unixtime,velocity

contains

subroutine makedatestr(datestr)
  integer,dimension(8) :: datevalues
  character(len=24)    :: datestr
  CALL date_and_time(values=datevalues)
  write(datestr,"(I0.4,'-',I0.2,'-',I0.2,X,I0.2,':',I0.2,':',I0.2)") &
    datevalues(1), datevalues(2), datevalues(3), datevalues(5), datevalues(6), datevalues(7)
end subroutine makedatestr

!> Some code to calculate UNIX time_ts
function is_leap_year(year)
  integer, intent(in) :: year
  integer :: is_leap_year

  ! http://www.timeanddate.com/date/leapyear.html says:
  !
  ! In the Gregorian calendar, which is the calendar used by
  ! most modern countries, the following rules decides which
  ! years are leap years:
  ! 
  ! 1. Every year divisible by 4 is a leap year.
  ! 
  ! 2. But every year divisible by 100 is NOT a leap year
  ! 
  ! 3. Unless the year is also divisible by 400, then it is
  ! still a leap year
  !
  ! The original normative standard is the papal bull
  ! "Inter Gravissimas", of 1582.
  !

  if (modulo(year,400) == 0) then
    is_leap_year = 1
  else
    if (modulo(year,100) == 0) then
      is_leap_year = 0
    else
      if (modulo(year,4) == 0) then
        is_leap_year = 1
      else
        is_leap_year = 0
      endif
    endif
  endif

end function is_leap_year

function days_in_year(year)
  integer, intent(in) :: year
  integer :: days_in_year
  if (1 == is_leap_year(year)) then
    days_in_year = 366
  else
    days_in_year = 365
  end if
end function days_in_year

!> Return the number of days between January 1 1970 and January 1 of
!> given year
function ydays_since_epoch(year)
  integer :: i
  integer :: total
  integer :: year
  integer :: ydays_since_epoch
  total = 0
  do i=1970,(year-1)
    total = total + days_in_year(i)
  end do
  ydays_since_epoch = total
end function ydays_since_epoch

!> Return the time_t, i.e. the number of seconds elapsed
!> since 00:00:00 UTC, January 1, 1970.
!>
!> Does not correct for leap seconds, consistent with POSIX.1
function unixtime()
  implicit none
  integer :: unixtime
  character(len=10) :: datestr,timestr,zonestr
  integer, dimension(8) :: values
  integer, dimension(12) :: monthlengths = &
          (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer :: days_since_epoch

  integer :: i
  integer :: delta_seconds ! Timezone offset

  call date_and_time(datestr,timestr,zonestr,values)
  ! values returns with:
  ! values(1) = year
  ! values(2) = month of year 
  ! values(3) = day of month
  ! values(4) = time different in minutes wrt utc
  ! values(5) = hour of day
  ! values(6) = minutes of hour
  ! values(7) = seconds of minute
  ! values(8) = milliseconds of second

  ! Figure out if this is a leap year, and
  ! adjust accordingly.

  if (1 == is_leap_year(values(1))) then
    monthlengths(2) = 29
  end if

  ! Days between Jan 1 this year and Jan 1 1970
  days_since_epoch = ydays_since_epoch(values(1))

  ! Now add up days between start of this month and Jan 1
  do i=1,(values(2)-1)
    days_since_epoch = days_since_epoch + monthlengths(i)
  end do
  days_since_epoch = days_since_epoch + values(3)-1

  ! days_since_epoch is now number of complete days since Jan 1 1970
  ! now correct for timezone

  delta_seconds = values(4)*60

  ! Note that the Intel compiler *will* adjust the time when you
  ! set the TZ environment variable, but it won't set values(4),
  ! so you can't tell whether you're in UTC or not.
  !
  !  FORTRAN! Because your blood pressure is too low, and there
  !  aren't enough head-shaped dents in the wall!
  !

  unixtime = values(7)                    &
            + 60*( values(6)               &
            + 60*( values(5)               &
            + 24* days_since_epoch ))      &
            - delta_seconds

end function unixtime




end module lb3d_helper_module
