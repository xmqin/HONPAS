!> \brief YAML output of SIESTA variables for external post-processors.
!!
!! This module provides facilities to output SIESTA variables into a YAML
!! file at the end of a run.
!!
!! Variables are dumped in YAML format to provide structured data that
!! can be read partially or fully, in any order, at the option of the user.
!!
!! The resulting output can be processed in C/C++ using LibYAML (found at
!! https://github.com/yaml/libyaml), in Python through the ruamel.yaml
!! package (found at https://pypi.python.org/pypi/ruamel.yaml), and in Perl
!! thanks to the YAML::XS module (found at
!! http://search.cpan.org/~tinita/YAML-LibYAML-0.69/lib/YAML/XS.pod).
!!
!! \author Yann Pouillon
!! \date 2017-2018
!! \copyright GNU General Public License version 3
!!
!! \note The interest of using the YAML format is that output can be achieved
!!       without introducing new external dependencies.
!!
module m_io_yaml

  use precision, only: dp

  implicit none

  private

  ! Ensure compliance with the YAML 1.2 file format
  character(len=*), parameter :: CH10 = achar(10)   !< New line
  character(len=*), parameter :: CH34 = achar(34)   !< Double quote
  character(len=*), parameter :: YAML_HEADER = "%YAML 1.2"//CH10//"---"//CH10
  character(len=*), parameter :: YAML_FOOTER = CH10//"..."

  public :: siesta_write_yaml

contains

  !> \brief Creates a YAML file containing SIESTA build parameters and
  !!        final energy values.
  !!
  !! This routine calls io_assign() to get a free unit number and outputs
  !! valid YAML data structures into the associated OUTVARS.yml file. The
  !! data consists in two dictionaries:
  !!   - siesta: build parameters and status of optional features;
  !!   - energies: full decomposition of the total energy, using exactly
  !!     the same naming conventions as in the \ref m_energies module.
  !!
  !! The resulting file can easily be parsed using the ruamel.yaml Python
  !! module.
  !!
  !! \todo Output forces and stress tensor.
  subroutine siesta_write_yaml()

    use m_energies
    use version_info

    implicit none

    logical :: trigger
    integer :: ierr, yaml_fd

    ! Open YAML document
    call io_assign(yaml_fd)
    open(unit=yaml_fd, file='OUTVARS.yml', status='new', action='write', &
&     access='sequential', form='formatted', iostat=ierr)
    ! FIXME: Find out why the system sometines reports a failure while things
    !        have gone perfectly well.
    !if ( ierr .ne. 0 ) call die('could not open OUTVARS.yml')
    write(unit=yaml_fd, fmt='(A)') YAML_HEADER

    ! Dump SIESTA information
    write(unit=yaml_fd, fmt='(A)') "siesta:"
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "version", &
&     CH34, trim(adjustl(version_str)), CH34
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "arch", &
&     CH34, trim(adjustl(siesta_arch)), CH34
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "compiler", &
&     CH34, trim(adjustl(compiler_version)), CH34
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "fflags", &
&     CH34, trim(adjustl(fflags)), CH34
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "fppflags", &
&     CH34, trim(adjustl(fppflags)), CH34
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "libs", &
&     CH34, trim(adjustl(libs)), CH34
#ifdef MPI
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "build", &
&     CH34, "mpi", CH34
#else
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "build", &
&     CH34, "serial", CH34
#endif
#ifdef _OPENMP
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "openmp", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef USE_GEMM3M
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "gemm3m", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef CDF
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "netcdf", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef NCDF_4
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "netcdf4", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef NCDF_PARALLEL
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "netcdf4_mpi", &
&     CH34, trim(yesno(trigger)), CH34
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "metis", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef SIESTA__CHESS
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "chess", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef SIESTA__ELPA
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "elpa", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef SIESTA__FLOOK
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "flook", &
&     CH34, trim(yesno(trigger)), CH34
#ifdef SIESTA__PEXSI
    trigger = .true.
#else
    trigger = .false.
#endif
    write(unit=yaml_fd, fmt='(2X,A,":",1X,3(A))') "pexsi", &
&     CH34, trim(yesno(trigger)), CH34

    ! Dump energies
    write(unit=yaml_fd, fmt='(A,A)') CH10, "energies:"
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Ebs", Ebs
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Eions", Eions
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Ena", Ena
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Ekin", Ekin
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Enl", Enl
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Eso", Eso
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Edftu", Edftu
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "DEna", DEna
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "DUscf", DUscf
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "DUext", DUext
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Exc", Exc
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Ecorrec", Ecorrec
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Emadel", Emad
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Emeta", Emeta
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Emolmec", Emm
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Ekinion", Ekinion
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Eharris", Eharrs+Ekinion
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "Etot", Etot+Ekinion
    write(unit=yaml_fd, fmt='(2X,A,":",1X,E24.8)') "FreeEng", FreeE+Ekinion

    ! Close YAML document
    write(unit=yaml_fd, fmt='(A)') YAML_FOOTER
    call io_close(yaml_fd)

  end subroutine siesta_write_yaml

  !> \brief Internal function to translate booleans into "yes"/"no" strings.
  !!
  !! This function takes a boolean condition as input and returns a string
  !! corresponding to the boolean value.
  !!
  !! \param[in] cond: boolean condition
  !! \return string equal to "yes" for .true., and "no" for .false.
  function yesno(cond) result(word)

    logical, intent(in) :: cond

    character(len=3) :: word

    if ( cond ) then
      word = "yes"
    else
      word = "no "
    end if

  end function yesno

end module m_io_yaml
