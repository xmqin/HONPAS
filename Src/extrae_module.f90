!*****************************************************************************
!*                        ANALYSIS PERFORMANCE TOOLS                         *
!*                                   Extrae                                  *
!*              Instrumentation package for parallel applications            *
!*****************************************************************************
!*     ___     This library is free software; you can redistribute it and/or *
!*    /  __         modify it under the terms of the GNU LGPL as published   *
!*   /  /  _____    by the Free Software Foundation; either version 2.1      *
!*  /  /  /     \   of the License, or (at your option) any later version.   *
!* (  (  ( B S C )                                                           *
!*  \  \  \_____/   This library is distributed in hope that it will be      *
!*   \  \__         useful but WITHOUT ANY WARRANTY; without even the        *
!*    \___          implied warranty of MERCHANTABILITY or FITNESS FOR A     *
!*                  PARTICULAR PURPOSE. See the GNU LGPL for more details.   *
!*                                                                           *
!* You should have received a copy of the GNU Lesser General Public License  *
!* along with this library; if not, write to the Free Software Foundation,   *
!* Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA          *
!* The GNU LEsser General Public License is contained in the file COPYING.   *
!*                                 ---------                                 *
!*   Barcelona Supercomputing Center - Centro Nacional de Supercomputacion   *
!*****************************************************************************

!* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
!| @file: $HeadURL: https://svn.bsc.es/repos/ptools/extrae/branches/2.4/include/extrae_module.f90 $
!| @last_commit: $Date: 2013-06-06 12:48:56 +0200 (Thu, 06 Jun 2013) $
!| @version:     $Revision: 1795 $
!* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

      module EXTRAE_MODULE
      integer, parameter :: INT4 = selected_int_kind(9)
      integer, parameter :: INT8 = selected_int_kind(18)

      integer(INT4) extrae_disable_all_options
      parameter (extrae_disable_all_options=0)

      integer(INT4) extrae_caller_option 
      parameter (extrae_caller_option=1)

      integer(INT4) extrae_hwc_option
      parameter (extrae_hwc_option=2)  

      integer(INT4) extrae_mpi_hwc_option
      parameter (extrae_mpi_hwc_option=4)

      integer(INT4) extrae_mpi_option
      parameter (extrae_mpi_option=8)

      integer(INT4) extrae_omp_option
      parameter (extrae_omp_option=16)

      integer(INT4) extrae_omp_hwc_option
      parameter (extrae_omp_hwc_option=32)

      integer(INT4) extrae_uf_hwc_option
      parameter (extrae_uf_hwc_option=64)

      integer(INT4) extrae_sampling_option
      parameter (extrae_sampling_option=128)

      integer(INT4) extrae_enable_all_options
      parameter (extrae_enable_all_options=255)

      interface

         subroutine extrae_init
         end subroutine extrae_init

         subroutine extrae_fini
         end subroutine extrae_fini

         subroutine extrae_version (major, minor, revision)
           integer, parameter :: INT4 = selected_int_kind(9)
         integer(INT4), intent (out) :: major
         integer(INT4), intent (out) :: minor
         integer(INT4), intent (out) :: revision
         end subroutine extrae_version

         subroutine extrae_shutdown
         end subroutine extrae_shutdown

         subroutine extrae_restart
         end subroutine extrae_restart

         subroutine extrae_event (extrae_type, extrae_value)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer, parameter :: INT8 = selected_int_kind(18)
         integer(INT4), intent(in) :: extrae_type
         integer(INT8), intent(in) :: extrae_value
         end subroutine extrae_event

         subroutine extrae_eventandcounters (extrae_type, extrae_value)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer, parameter :: INT8 = selected_int_kind(18)
         integer(INT4), intent(in) :: extrae_type
         integer(INT8), intent(in) :: extrae_value
         end subroutine extrae_eventandcounters


         subroutine extrae_nevent (num, extrae_types, extrae_values)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer, parameter :: INT8 = selected_int_kind(18)
         integer(INT4), intent(in) :: num
         integer(INT4), intent(in) :: extrae_types
         integer(INT8), intent(in) :: extrae_values
         end subroutine extrae_nevent

         subroutine extrae_neventandcounters (num, extrae_types, &
           extrae_value)
         integer, parameter :: INT8 = selected_int_kind(18)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer(INT4):: num
         integer(INT4):: extrae_types
         integer(INT8):: extrae_value
         end subroutine extrae_neventandcounters

         subroutine extrae_counters
         end subroutine

         subroutine extrae_previous_hwc_set
         end subroutine extrae_previous_hwc_set

         subroutine extrae_next_hwc_set
         end subroutine extrae_next_hwc_set

         subroutine extrae_set_options (options)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer(INT4), intent(in) :: options
         end subroutine extrae_set_options

         subroutine extrae_set_tracing_tasks (task_from, task_to)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer(INT4), intent(in) :: task_from
         integer(INT4), intent(in) :: task_to
         end subroutine extrae_set_tracing_tasks

         subroutine extrae_user_function (enter)
         integer, parameter :: INT4 = selected_int_kind(9)
         integer(INT4), intent(in) :: enter
         end subroutine extrae_user_function

      end interface

      end module EXTRAE_MODULE

