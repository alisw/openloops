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


module ol_helicity_init

contains

! **********************************************************************
subroutine heltable(n_in, n, tab)
! **********************************************************************
  use KIND_TYPES, only: intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_fatal
  implicit none
  integer,           intent(in)  :: n_in(:)
  integer(intkind2), intent(out) :: n(:), tab(:,:)
  integer :: i1, i2, i3, i4, i5, i6, i7

  n = n_in

  if (size(n) == 3) then
    i3 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        tab(1,i3) = i1
        tab(2,i3) = i2
        i3 = i3 + 1 ! <=> i3 = i2 + n(2) * (i1-1)
      end do
    end do
  else if (size(n) == 4) then
    i4 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          tab(1,i4) = i1
          tab(2,i4) = i2
          tab(3,i4) = i3
          i4 = i4 + 1 ! <=> i4 = i3 + n(3) * (i2-1 + n(2) * (i1-1))
        end do
      end do
    end do
  else if (size(n) == 5) then
    i5 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          do i4 = 1, n(4)
            tab(1,i5) = i1
            tab(2,i5) = i2
            tab(3,i5) = i3
            tab(4,i5) = i4
            i5 = i5 + 1
          end do
        end do
      end do
    end do
  else if (size(n) == 6) then
    i6 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          do i4 = 1, n(4)
            do i5 = 1, n(5)
              tab(1,i6) = i1
              tab(2,i6) = i2
              tab(3,i6) = i3
              tab(4,i6) = i4
              tab(5,i6) = i5
              i6 = i6 + 1
            end do
          end do
        end do
      end do
    end do
  else if (size(n) == 7) then
    i7 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          do i4 = 1, n(4)
            do i5 = 1, n(5)
              do i6 = 1, n(6)
                tab(1,i7) = i1
                tab(2,i7) = i2
                tab(3,i7) = i3
                tab(4,i7) = i4
                tab(5,i7) = i5
                tab(6,i7) = i6
                i7 = i7 + 1
              end do
            end do
          end do
        end do
      end do
    end do
  else
    call ol_fatal("heltable: " // trim(to_string(size(n))) // " point vertices are not supported.")
  end if
end subroutine heltable



! **********************************************************************
subroutine helbookkeeping_flip(hel, k, shift, eflip, exthel, firstpol)
! ----------------------------------------------------------------------
! THIS SUBROUTINE IS NOT USED
! ----------------------------------------------------------------------
! determines helicity state h = exthel(e+1,k) = 1 or 2 of ext. particle k
! corresponding to global helicity label e
! and global helicity label e -> eflip(e+1,k) resulting
! from flip of first two helicity states, hel(1:2), of k-th external-particle
! (eflip = e for particles with /= 2 helicity states)
! eflip is needed for gluon-helicity flips in dipole subtraction
! **********************************************************************
  use KIND_TYPES, only: intkind2
  implicit none
  integer,            intent(in)    :: hel(:), k, shift
  integer(intkind2),  intent(inout) :: eflip(:,:)
  integer,            intent(inout) :: exthel(:,:), firstpol(:)
  integer :: e, h, i, nhelmax, de_flip

  firstpol(k) = hel(1)    ! first helicity state of particle k

  nhelmax = size(eflip,1) ! max number of global helicity states
  e = 0                   ! global helicity label

  do while (e <= nhelmax - 1)
    do h = 1, size(hel) ! helicity states of particle k
      if (size(hel) == 2) then
        select case (h)
        case (1)
          de_flip =  shift ! 1 -> 2 flip
        case (2)
          de_flip = -shift ! 2 -> 1 flip
        end select
      else
        de_flip = 0        ! no flip
      end if
      do i = 1, shift ! helicity states of particles 1,...,k-1
        eflip(e+1,k)  = e + de_flip ! flips of helicty state e stored in (e+1)th row of eflip array (0 <= e <= nhelmax-1)
        exthel(e+1,k) = h
        e = e+1
      end do
    end do
  end do

end subroutine helbookkeeping_flip



! **********************************************************************
subroutine helsync_flip(nsync, nhel, hel, eflip, exthel)
! ----------------------------------------------------------------------
! nsync
! nhel                = number of non-vanishing helicity configurations
! hel(1:nhel)         = %e helicity labels of nonvanishing helicity configurations
! hel(nhel+1:nhelmax) = -1 (vanishing helicity config)
! eflip(e+1,k) in     = %e helicity label resulting from flip of particle-k helicity
! eflip(h,k)   out    = same with e-> non-vanishing config index h, hel(h)=e, i.e.
!                       hel(eflipout(h,k)) = eflip(hel(h)+1,k)
!                       when flip gives vanishing config => eflipout = nhel+1 (convention)
! exthel     in -> out: analogous restriction to non-vanishing helicity configurations
! nhelmax             = maximal number of helicity configurations
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer(intkind1), intent(in)    :: nsync
  integer(intkind2), intent(in)    :: nhel, hel(:)
  integer(intkind2), intent(inout) :: eflip(:,:)
  integer,           intent(inout) :: exthel(:,:)
  integer(intkind2) :: helinv(0:size(hel)-1)
  integer :: nhelmax, e, h, k

  if (nsync /= 1) then
    call ol_error(2, 'in subroutine helsync_flip:')
    call ol_error(2, 'nsync = ' // to_string(nsync) // ' not allowed')
    call ol_fatal()
  end if

  nhelmax = size(hel) ! maximal number of helicity configurations
  if (nhel < nhelmax) then
    do e = 0, nhelmax - 1
      helinv(e) = nhel + 1 ! spurious helicity state for vanishing configurations (initialisation)
    end do
  end if

  do h = 1, nhel
    helinv(hel(h)) = h ! non-vanishing helicity states
  end do

  do h = 1, nhel ! non-vanishing helicity states
    do k = 1, size(eflip,2) ! external particles
      eflip(h,k)  = helinv(eflip(hel(h)+1,k))
      exthel(h,k) = exthel(hel(h)+1,k)
    end do
  end do

  if (nhel < nhelmax) then
    do k = 1, size(eflip,2) ! external particle
      eflip(nhel+1:nhelmax,k)  = nhel + 1 ! spurious helicity state for vanishing configurations
      exthel(nhel+1:nhelmax,k) = 0        ! idem
    end do
  end if

end subroutine helsync_flip

end module ol_helicity_init
