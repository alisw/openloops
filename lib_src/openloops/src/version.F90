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


module ol_version
  implicit none
  character(16) :: version = VERSION
  integer :: process_api = PROCESSAPI
  character(7) :: revision = REVISION
  logical :: splash_todo = .true.
  integer, parameter :: welcome_length = 1300

  contains

  subroutine welcome(outstring)
    implicit none
    character(welcome_length), intent(out) :: outstring
    character(16) :: print_version
    character, parameter :: LF = char(10)
    if (len_trim(version) == 0) then
      print_version=""
    else
      print_version="version " // version
    end if
    outstring = &
      & " #####################################################################" // LF // &
      & " #        ___                                       " // adjustr(print_version) // " #" // LF // &
      & " #       /   \ ___  ____  _  _  |     __   __  ___   __    ___       #" // LF // &
      & " #       |   | |__| |__   |\ |  |    /  \ /  \ |__| /__    __/       #" // LF // &
      & " #       \___/ |    |___  | \|  |___ \__/ \__/ |    __/   /__        #" // LF // &
      & " #                                                                   #" // LF // &
      & " #####################################################################" // LF // &
      & " #       You are using OpenLoops 2 to evaluate loop amplitudes       #" // LF // &
      & " #                             Authors:                              #" // LF // &
      & " #       F. Buccioni, J.-N. Lang, J. Lindert, P. Maierhoefer,        #" // LF // &
      & " #                S. Pozzorini, M. Zoller, H. Zhang                  #" // LF // &
      & " #                                                                   #" // LF // &
      & " #         Please cite Eur.Phys.J. C79 (2019) no.10, 866             #" // LF // &
      & " #                     Phys. Rev. Lett. 108 (2012) 111601            #" // LF // &
      & " #                     Eur.Phys.J. C78 (2018) no.1, 70               #" // LF // &
      & " #                                                                   #" // LF // &
      & " #####################################################################" // LF
    splash_todo = .false.
  end subroutine welcome

  subroutine welcome_c(outstring) bind(c,name="ol_welcome")
    use, intrinsic :: iso_c_binding, only: c_char
    implicit none
    ! DO NOT CHANGE THE LENGTH OF THE WELCOME STRING,
    ! because it is hard-coded to 700 characters in Sherpa.
    character(kind=c_char), intent(out) :: outstring(700)
    character, parameter :: LF = char(10)
    integer :: i
    character(700) :: welcome_str = &
      & "####################################################################" // LF // &
      & "#        ___  ___  ____  _  _       __   __  ___   __    ___       #" // LF // &
      & "#       /   \ |__| |__   |\ |  |   /  \ /  \ |__| /__    __/       #" // LF // &
      & "#       \___/ |    |___  | \|  |__ \__/ \__/ |    __/   /__        #" // LF // &
      & "#       You are using OpenLoops 2 to evaluate loop amplitudes      #" // LF // &
      & "#        F. Buccioni, J.-N. Lang, J. Lindert, P. Maierhoefer,      #" // LF // &
      & "#               S. Pozzorini, M. Zoller, H. Zhang                  #" // LF // &
      & "#               Eur.Phys.J. C79 (2019) no.10, 866                  #" // LF // &
      & "#                        arXiv:1907.13071                          #" // LF // &
      & "####################################################################"
    if (len(trim(welcome_str)) >= 699) then
      ! Again, so that nobody does anything stupid without noticing.
      print *, "ol_welcome(): welcome string is too long"
      stop
    end if
    do i = 1, len(trim(welcome_str))
      outstring(i) = welcome_str(i:i)
    end do
    outstring(i) = char(0)
    splash_todo = .false.
  end subroutine welcome_c

  subroutine print_welcome()
    use ol_debug, only: ol_write_msg
    implicit none
    character(welcome_length) :: welcome_string
    call welcome(welcome_string)
    call ol_write_msg(trim(welcome_string))
  end subroutine print_welcome


  subroutine openloops_version_string(outstring)
    character(16) :: outstring
    outstring = trim(version)
  end subroutine openloops_version_string

  subroutine openloops_version_string_c(outstring) bind(c,name="ol_version_string")
    use, intrinsic :: iso_c_binding, only: c_char
    implicit none
    character(kind=c_char), intent(out) :: outstring(17)
    character, parameter :: LF = char(10)
    integer i
    do i = 1, len(trim(version))
      outstring(i) = version(i:i)
    end do
    outstring(i) = char(0)
  end subroutine openloops_version_string_c

end module ol_version


module openloops_version
  use ol_version, only: version, process_api, revision, splash_todo
end module openloops_version
