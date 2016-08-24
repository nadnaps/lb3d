! Variable type

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

#define LBE_REAL MPI_REAL8

! Single fluid implies no surfactant
#ifdef SINGLEFLUID
#ifndef NOSURFACTANT
#define NOSURFACTANT
#endif
#endif

! Define defaults
#ifndef FORCE_BODY
#define FORCE_BODY
#endif
#ifndef FORCE_GUO
#define FORCE_GUO
#endif
#ifndef FORCE_SHANCHEN
#define FORCE_SHANCHEN
#endif
#ifndef COMPAT_BDIST
#define COMPAT_BDIST
#endif


! Electrolytes implies MD (for now)
!#ifdef ELEC
!#ifndef MD
!#define MD
!#endif
!#endif

! Debug macro
#ifdef DEBUG_CHECKPOINT
#define DEBUG_CHECKPOINT_MSG(msg) CALL log_msg(msg,.false.)
#define DEBUG_CHECKPOINT_MSG_WS(msg) CALL log_msg_ws(msg,.false.)
#define DEBUG_CHECKPOINT_MSG_ALL(msg) CALL log_msg(msg,.true.)
#else
#define DEBUG_CHECKPOINT_MSG(msg)
#define DEBUG_CHECKPOINT_MSG_WS(msg)
#define DEBUG_CHECKPOINT_MSG_ALL(msg)
#endif

! With the DEBUG_HDF5 flag set, show the debug output.
! Otherwise, replace the macros with nothing.
#ifdef DEBUG_HDF5
#define DEBUG_HDF5_MSG(msg) CALL log_msg(msg,.false.)
#define DEBUG_HDF5_MSG_WS(msg) CALL log_msg_ws(msg,.false.)
#define DEBUG_HDF5_MSG_ALL(msg) CALL log_msg(msg,.true.)
#else
#define DEBUG_HDF5_MSG(msg)
#define DEBUG_HDF5_MSG_WS(msg)
#define DEBUG_HDF5_MSG_ALL(msg)
#endif

#ifdef DEBUG_MPI
#define DEBUG_MPI_MSG(msg) CALL log_mpi(msg)
#else
#define DEBUG_MPI_MSG(msg)
#endif

