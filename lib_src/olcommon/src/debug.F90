!******************************************************************************!
! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
!                                                                              !
! This file is part of OpenLoops.                                              !
!                                                                              !
! OpenLoops is free software: you can redistribute it and/or modify            !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! OpenLoops is distributed in the hope that it will be useful,                 !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
! GNU General Public License for more details.                                 !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.           !
!******************************************************************************!

module ol_debug
  ! precision independent message and error routines
  use, intrinsic :: iso_c_binding, only:  c_int
  use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, &
                                            stderr=>error_unit
  implicit none
  private
  public :: set_verbose, get_verbose, get_error, ol_msg, ol_error, ol_fatal, ol_fatal_update
  public :: ol_write_msg
  public :: error, verbose, do_not_stop
  public :: olodebug_unit

  integer, save :: verbose = 0
  integer, save :: error = 0
  logical, save :: do_not_stop = .false.
  integer, save :: olodebug_unit = stdout ! if >= 0, used in CutTools to print olo() calls which caused an error

  interface ol_msg
    module procedure ol_print_msg_level, ol_print_msg
  end interface ol_msg

  interface ol_error
    module procedure ol_error_msg, ol_error_level
  end interface ol_error

  contains

  subroutine set_verbose(level)
    implicit none
    integer, intent(in) :: level
    verbose = level
  end subroutine set_verbose

  subroutine get_verbose(level)
    implicit none
    integer, intent(out) :: level
    level = verbose
  end subroutine get_verbose

  subroutine ol_print_msg(msg)
    implicit none
    character(len=*), intent(in) :: msg
    call ol_msg(0,msg)
  end subroutine ol_print_msg

  subroutine ol_print_msg_level(level, msg)
    implicit none
    integer, intent(in) :: level
    character(len=*), intent(in) :: msg
    if (verbose >= level) call ol_write_msg("[OpenLoops] " // trim(msg))
  end subroutine ol_print_msg_level

  subroutine ol_error_level(err, msg)
    implicit none
    integer, intent(in) :: err
    character(len=*), intent(in), optional :: msg
    character(len=100) :: err_format
    error = err
    if (present(msg)) call ol_write_msg("[OpenLoops] Error: " // trim(msg), stderr)
  end subroutine ol_error_level

  subroutine ol_error_msg(msg)
    implicit none
    character(len=*), intent(in), optional :: msg
    call ol_error_level(1,msg)
  end subroutine ol_error_msg

  function get_error()
    implicit none
    integer :: get_error
    get_error = error
  end function get_error

  function get_error_c() bind(c,name="ol_get_error")
    implicit none
    integer(c_int) :: get_error_c
    get_error_c = error
  end function get_error_c

  subroutine ol_fatal(msg, fatal_err)
    implicit none
    character(len=*), optional, intent(in) :: msg
    integer, optional, intent(out) :: fatal_err
    error = 2
    if (present(msg)) call ol_write_msg("[OpenLoops] ERROR: " // trim(msg), stderr)
    if (present(fatal_err)) then
      fatal_err = 1
    else if (do_not_stop) then
      if (verbose > 0) call ol_write_msg("[OpenLoops] FATAL ERROR.", stderr)
    else
      if (verbose > 0) call ol_write_msg("[OpenLoops] STOP.", stderr)
      stop
    end if
  end subroutine ol_fatal

  subroutine ol_fatal_update(fatal_err)
    implicit none
    integer, optional, intent(out) :: fatal_err
    if (present(fatal_err)) then
      fatal_err = 1
    else
      call ol_fatal("Process library not up-to-date. Please run: ./openloops update --processes")
    end if
  end subroutine ol_fatal_update

  subroutine ol_write_msg(msg, unit)
    implicit none
    integer, intent(in), optional :: unit
    character(len=*), intent(in) :: msg
    character(len=100) :: format_msg
    write (format_msg, '("(a", I4, ")")' )  len_trim(msg)
    if (present(unit)) then
      write(unit,format_msg) trim(msg)
    else
      write(stdout,format_msg) trim(msg)
    end if
  end subroutine ol_write_msg

end module ol_debug
