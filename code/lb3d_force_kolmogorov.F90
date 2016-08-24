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

!> forcing that is sinusoidally modulated in x-direction, producing
!> Kolmogorov flow
module lb3d_force_kolmogorov_module

    use lb3d_global_module

    use lb3d_mpi_module
    use lb3d_log_module
    use lb3d_config_module, only: arg_input_dfile,arg_input_dfile_p,inp_file,nx

    implicit none

    private

    public force_init_kolmogorov,force_input_kolmogorov,force_apply_kolmogorov

    !> buffer to avoid force re-calculation at each time step and lattice site
    real(kind=rk),save,allocatable :: kolmogorov_flow_forces(:)

    real(kind=rk),save :: amplitude !< maximum absolute value of force

    namelist /force_kolmogorov/ amplitude

contains

    !> apply sinusoidally modulated forces for Kolmogorov flow
    subroutine force_apply_kolmogorov
        integer x

        do x=1,nx
           lbe_force(3,1,x,:,:) = kolmogorov_flow_forces(x)
        end do
    end subroutine force_apply_kolmogorov

    !> mainly pre-calculate position-dependent force
    subroutine force_init_kolmogorov
        integer stat,x

#ifndef SINGLEFLUID
        call log_msg("WARNING: force_apply_kolmogorov() not implemented for"&
             &//" multi-component case",.false.)

#endif
        use_lbe_force = .true.

        allocate (kolmogorov_flow_forces(1:nx),stat=stat)
        call check_allocate(stat&
             &,'force_init_kolmogorov():kolmogorov_flow_forces')

        do x=1,nx
           kolmogorov_flow_forces(x) &
                &= sin(2.0_rk*pi*(real(start(1)+x-1,kind=rk)-0.5_rk)/tsize(1))&
                &*amplitude
        end do
    end subroutine force_init_kolmogorov

    !> read namelist
    subroutine force_input_kolmogorov
        integer ierror
        character(len=128)     :: msgstr

        call log_ws(.false.)
        call log_msg("----( Reading FORCE KOLMOGOROV input )-----",.false.)
        call log_ws(.false.)

        if (myrankc==0) then
           open (unit=input_file_unit,file=trim(inp_file),err=100)
           read (unit=input_file_unit,nml=force_kolmogorov,err=100)
           close (unit=input_file_unit,err=100)
        end if

        if (arg_input_dfile_p > 0) then
           call log_msg("  Getting differential input...",.false.)
           open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
           read (unit=input_dfile_unit,nml=force_kolmogorov,iostat=ierror)
           if (ierror/=0) then
              call log_msg("    WARNING: Differential namelist not found "&
                   &//'or errors encountered.',.false.)
           end if
           close (unit=input_dfile_unit)
           call log_ws(.false.)
        end if

        write (msgstr,"('amplitude = ',F16.10)") amplitude
        call log_msg(trim(msgstr),.false.)
        call log_ws(.false.)

        call MPI_Bcast(amplitude,1,LBE_REAL,0,comm_cart,ierror)

        return
100     continue
        call error('Error reading input file "'//trim(inp_file))
    end subroutine force_input_kolmogorov

end module lb3d_force_kolmogorov_module
