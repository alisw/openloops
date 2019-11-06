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

module ol_helicity_bookkeeping_/**/REALKIND
  implicit none
  contains

! **********************************************************************
subroutine checkzero_scalar(S)
! ----------------------------------------------------------------------
! sets to zero irrelevant components of scalar wfun in order to
! find vanishing scalar subtrees via generic S(h)%j == 0 test
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(inout) :: S(:)
  integer :: h
  do h = 1, size(S)
    S(h)%j(2:4) = 0
  end do
end subroutine checkzero_scalar



! **********************************************************************
subroutine helbookkeeping_wf(hel, ex, shift)
! ----------------------------------------------------------------------
! attributes additive global labels (0, 1,..., n_k-1)*shift_k
! to helicity states (1,2,...,n_k) of kth external particle, where
! shift_1 = 1 and shift_k+1 = n_k*shift_k
! shift = shift_k + 1 is returned as output to fix next particles labels
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer,           intent(in)    :: hel(:)
  type(wfun),        intent(inout) :: ex(size(hel))
  integer,           intent(inout) :: shift
  integer                          :: h

  do h = 1, size(hel)
    if (all(ex(h)%j == 0)) then
      ex(h)%e = -1 ! marks vanishing helicity configurations
    else
      ex(h)%e = (h-1) * shift
    end if
  end do
  shift = shift * size(Hel)

end subroutine helbookkeeping_wf



! **********************************************************************
subroutine helbookkeeping_prop(ntry, WF1, WF2, n)
! ----------------------------------------------------------------------
! WF1(1:n) = input  wfun array
! WF2(1:n) = output wfun array
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  type(wfun),        intent(in)    :: WF1(n)
  type(wfun),        intent(out)   :: WF2(n)
  integer(intkind2) :: h1, i

  if (ntry > 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_prop:')
    call ol_fatal('ntry =' // to_string(ntry) // 'not allowed')
  end if

  ! sets n = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n
  do i = 1, n
    if (WF1(i)%e == -1_intkind2) then
      h1 = i - 1
      exit
    end if
    WF2(i)%e = WF1(i)%e ! sets WF2 helicity label
  end do

  do i = h1 + 1, n
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_prop:')
      call ol_error(2,'i, h1, n, WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " // &
                &   to_string(n) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
    WF2(i)%e = -1_intkind2  ! mark vanishing helicity configurations
  end do
  n = h1

end subroutine helbookkeeping_prop



! **********************************************************************
subroutine helbookkeeping_vert3(ntry, WF1, WF2, WF3, n, t)
! ----------------------------------------------------------------------
! WF1(1:n(1)), WF(1:n(2)) = input wfun arrays
! WF3(1:n(3))             = output wfun array
! vanishing components moved at the end of WF3 array and marked via %e=-1
! non-zero WF3 components ordered according to global helicity label
! WF1(h1), WF2(h2) <-> WF3(h3) connection stored in table hi= t(i,h3) with i=1,2
! array sizes n(1),n(2),n(3) restricted to non-vanishing components
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2))
  type(wfun),        intent(out)   :: WF3(n(3))
  integer(intkind2) :: h1, h2, h3, i, n2_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(2)
  type(wfun)        :: WFaux

  if (ntry /= 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_vert3:')
    call ol_error(2,' ntry =' // to_string(ntry) // ' not allowed')
    call ol_fatal()
  end if

  ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n(1)
  do i = 1, n(1)
    if (WF1(i)%e == -1_intkind2) then
      h1 = i-1
      exit
    end if
  end do
  do i = h1+1, n(1)
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert3:')
      call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                  &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
  end do
  n(1) = h1

  ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
  h2 = n(2)
  do i = 1, n(2)
    if (WF2(i)%e == -1_intkind2) then
      h2 = i - 1
      exit
    end if
  end do
  do i = h2+1, n(2)
    if (WF2(i)%e /= -1_intkind2) then
      call ol_error(2, 'in subroutine helbookkeeping_vert3:')
      call ol_error(2,'i, h1, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
      call ol_fatal()
    end if
  end do
  n2_in = n(2)
  n(2) = h2

  ! restrict I/O table to non-vanishing helicity configurations
  i  = 0 ! index of WF3 states
  h3 = 0 ! index of non-zero WF3 states
  do h1 = 1, n(1)
    do h2 = 1, n(2)
      i = (h1-1)*n2_in + h2         ! assumes input table in standard inititalisation form
      if (all(WF3(i)%j == 0)) cycle ! skips vanishing WF3 components
      h3 = h3 + 1
      t(1,h3) = h1
      t(2,h3) = h2
      WF3(h3)%e = WF1(h1)%e + WF2(h2)%e ! determines additive helicity label for outgoing states
      if (h3 == i) cycle
      WF3(h3)%j = WF3(i)%j              ! shifts non-zero components to first part of WF3 array
      WF3(h3)%h = WF3(i)%h              ! shifts non-zero components to first part of WF3 array
    end do
  end do

  ! put vanishing components at the end of WF3 array
  do i = h3 + 1, n(3)
    WF3(i)%j = 0
    WF3(i)%h = B"00"
    WF3(i)%e = -1_intkind2 ! marker for vanising helicity states
  end do

  ! sets n(3) = # of non-vanishing I/O helicity configurations
  n(3) = h3

  ! sort non-vanishing WF3(1:h3) components and adapt table accordingly
  do i = 1, n(3)-1
    emin = WF3(i)%e
    imin = i
    do iq = i+1, n(3)             ! look for component with helicity < emin
      if (WF3(iq)%e >= emin) cycle
      emin = WF3(iq)%e
      imin = iq
    end do
    if (imin > i) then
       WFaux       = WF3(i)          ! flip WF3(i) <-> WF3(imin) components
       WF3(i)      = WF3(imin)
       WF3(imin)   = WFaux
       taux        = t(1:2,i)      ! same for table
       t(1:2,i)    = t(1:2,imin)
       t(1:2,imin) = taux
    end if
  end do

end subroutine helbookkeeping_vert3



! **********************************************************************
subroutine helbookkeeping_vert4(ntry, WF1, WF2, WF3, WF4, n, t)
! ----------------------------------------------------------------------
! WF1(1:n(1)),WF2(1:n(2)),WF3(1:n(3)) = input wfun arrays
! WF4(1:n(4))                         = output wfun array
! vanishing components moved at the end of WF4 array and marked via %e = -1
! non-zero WF4 components ordered according to global helicity label
! WF1(h1), WF2(h2), WF3(h3) <-> WF4(h4) connection stored in table hi = t(i,h4) with i=1,2,3
! array sizes n(1), n(2), n(3), n(4) restricted to non-vanishing components
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2)), WF3(n(3))
  type(wfun),        intent(out)   :: WF4(n(4))
  integer(intkind2) :: h1, h2, h3, h4, i, n2_in, n3_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(3)
  type(wfun)        :: WFaux

  if(ntry /= 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_vert4:')
    call ol_error(2,'ntry =' // to_string(ntry) // ' not allowed')
    call ol_fatal()
  end if

  ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n(1)
  do i = 1, n(1)
    if (WF1(i)%e == -1_intkind2) then
      h1 = i - 1
      exit
    end if
  end do
  do i = h1 + 1, n(1)
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert4:')
      call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                  &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
  end do
  n(1) = h1

  ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
  h2 = n(2)
  do i = 1, n(2)
    if (WF2(i)%e == -1_intkind2) then
      h2 = i - 1
      exit
    end if
  end do
  do i = h2 + 1, n(2)
    if (WF2(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert4:')
      call ol_error(2,'i, h1, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
      call ol_fatal()
    end if
  end do
  n2_in = n(2)
  n(2) = h2

  ! sets n(3) = # of non-zero WF3 components and check that all zeros are at the end
  h3 = n(3)
  do i = 1, n(3)
    if (WF3(i)%e == -1_intkind2) then
      h3 = i - 1
      exit
    end if
  end do
  do i = h3 + 1, n(3)
    if (WF3(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert4:')
      call ol_error(2,'i, h3, n(3), WF3(i)%e =' // to_string(i) // " " // to_string(h3) // " " //  &
                  &  to_string(n(3)) // " " //  to_string(WF3(i)%e))
      call ol_fatal()
    end if
  end do
  n3_in = n(3)
  n(3) = h3


  ! restrict I/O table to non-vanishing helicity configurations
  i  = 0 ! index of WF4 states
  h4 = 0 ! index of non-zero WF4 states
  do h1 = 1, n(1)
    do h2 = 1, n(2)
      do h3 = 1, n(3)
        i  = h3 + n3_in * (h2-1 + n2_in * (h1-1)) ! assumes input table in standard inititalisation form
        if (all(WF4(i)%j == 0)) cycle ! skips vanishing WF4 components
        h4 = h4 + 1
        t(1,h4) = h1
        t(2,h4) = h2
        t(3,h4) = h3
        WF4(h4)%e = WF1(h1)%e + WF2(h2)%e + WF3(h3)%e ! additive helicity label for outgoing states
        if (h4 == i) cycle
        WF4(h4)%j = WF4(i)%j ! shifts non-zero components to first part of WF4 array
        WF4(h4)%h = WF4(i)%h ! shifts non-zero components to first part of WF4 array
      end do
    end do
  end do

  ! put vanishing components at the end of WF4 array
  do i = h4 + 1, n(4)
    WF4(i)%j = 0
    WF4(i)%h = B"00"
    WF4(i)%e = -1_intkind2 ! marker for vanising helicity states
  end do

  ! sets n(4) = # of non-vanishing I/O helicity configurations
  n(4) = h4

  ! sort non-vanishing WF4(1:h4) components and adapt table accordingly
  do i = 1, n(4) - 1
    emin = WF4(i)%e
    imin = i
    do iq = i + 1, n(4)              ! look for component with helicity < emin
      if (WF4(iq)%e >= emin) cycle
      emin = WF4(iq)%e
      imin = iq
    end do
    if (imin > i) then
       WFaux     = WF4(i) ! flip WF4(i) <-> WF4(imin) components
       WF4(i)    = WF4(imin)
       WF4(imin) = WFaux
       taux        = t(1:3,i) ! same for table
       t(1:3,i)    = t(1:3,imin)
       t(1:3,imin) = taux
    end if
  end do

end subroutine helbookkeeping_vert4


! **********************************************************************
subroutine helbookkeeping_vert5(ntry, WF1, WF2, WF3, WF4, WF5, n, t)
! ----------------------------------------------------------------------
! WF1(1:n(1)),WF2(1:n(2)),WF3(1:n(3)), WF3(1:n(4)) = input wfun arrays
! WF5(1:n(4))                                      = output wfun array
! vanishing components moved at the end of WF5 array and marked via %e = -1
! non-zero WF5 components ordered according to global helicity label
! WF1(h1), WF2(h2), WF3(h3), WF(h4) <-> WF5(h5) connection stored in
! table hi = t(i,h5) with i=1,2,3,4
! array sizes n(1), n(2), n(3), n(4), n(5) restricted to non-vanishing components
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2)), WF3(n(3)), WF4(n(4))
  type(wfun),        intent(out)   :: WF5(n(5))
  integer(intkind2) :: h1, h2, h3, h4, h5, i, n2_in, n3_in, n4_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(4)
  type(wfun)        :: WFaux

  if(ntry /= 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_vert5:')
    call ol_error(2,'ntry =' // to_string(ntry) // ' not allowed')
    call ol_fatal()
  end if

  ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n(1)
  do i = 1, n(1)
    if (WF1(i)%e == -1_intkind2) then
      h1 = i - 1
      exit
    end if
  end do
  do i = h1 + 1, n(1)
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                  &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
  end do
  n(1) = h1

  ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
  h2 = n(2)
  do i = 1, n(2)
    if (WF2(i)%e == -1_intkind2) then
      h2 = i - 1
      exit
    end if
  end do
  do i = h2 + 1, n(2)
    if (WF2(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h2, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
      call ol_fatal()
    end if
  end do
  n2_in = n(2)
  n(2) = h2

  ! sets n(3) = # of non-zero WF3 components and check that all zeros are at the end
  h3 = n(3)
  do i = 1, n(3)
    if (WF3(i)%e == -1_intkind2) then
      h3 = i - 1
      exit
    end if
  end do
  do i = h3 + 1, n(3)
    if (WF3(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h3, n(3), WF3(i)%e =' // to_string(i) // " " // to_string(h3) // " " //  &
                  &  to_string(n(3)) // " " //  to_string(WF3(i)%e))
      call ol_fatal()
    end if
  end do
  n3_in = n(3)
  n(3) = h3


  ! sets n(4) = # of non-zero WF4 components and check that all zeros are at the end
  h4 = n(4)
  do i = 1, n(4)
    if (WF4(i)%e == -1_intkind2) then
      h4 = i - 1
      exit
    end if
  end do
  do i = h4 + 1, n(4)
    if (WF4(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h4, n(4), WF4(i)%e =' // to_string(i) // " " // to_string(h4) // " " //  &
                  &  to_string(n(4)) // " " //  to_string(WF4(i)%e))
      call ol_fatal()
    end if
  end do
  n4_in = n(4)
  n(4) = h4

  ! restrict I/O table to non-vanishing helicity configurations
  i  = 0 ! index of WF5 states
  h5 = 0 ! index of non-zero WF5 states
  do h1 = 1, n(1)
    do h2 = 1, n(2)
      do h3 = 1, n(3)
        do h4 = 1, n(4)
          i  = h4 + n4_in * (h3-1 + n3_in * (h2-1 + n2_in * (h1-1))) ! assumes input table in standard inititalisation form
          if (all(WF5(i)%j == 0)) cycle ! skips vanishing WF5 components
          h5 = h5 + 1
          t(1,h5) = h1
          t(2,h5) = h2
          t(3,h5) = h3
          t(4,h5) = h4
          WF5(h5)%e = WF1(h1)%e + WF2(h2)%e + WF3(h3)%e + WF4(h4)%e ! additive helicity label for outgoing states
          if (h5 == i) cycle
          WF5(h5)%j = WF5(i)%j ! shifts non-zero components to first part of WF5 array
          WF5(h5)%h = WF5(i)%h ! shifts non-zero components to first part of WF5 array
        end do
      end do
    end do
  end do

  ! put vanishing components at the end of WF4 array
  do i = h5 + 1, n(5)
    WF5(i)%j = 0
    WF5(i)%h = B"00"
    WF5(i)%e = -1_intkind2 ! marker for vanising helicity states
  end do

  ! sets n(5) = # of non-vanishing I/O helicity configurations
  n(5) = h5

  ! sort non-vanishing WF5(1:h5) components and adapt table accordingly
  do i = 1, n(5) - 1
    emin = WF5(i)%e
    imin = i
    do iq = i + 1, n(5)              ! look for component with helicity < emin
      if (WF5(iq)%e >= emin) cycle
      emin = WF5(iq)%e
      imin = iq
    end do
    if (imin > i) then
       WFaux     = WF5(i) ! flip WF5(i) <-> WF5(imin) components
       WF5(i)    = WF5(imin)
       WF5(imin) = WFaux
       taux        = t(1:4,i) ! same for table
       t(1:4,i)    = t(1:4,imin)
       t(1:4,imin) = taux
    end if
  end do

end subroutine helbookkeeping_vert5


! **********************************************************************
subroutine helbookkeeping_vert6(ntry, WF1, WF2, WF3, WF4, WF5, WF6, n, t)
! ----------------------------------------------------------------------
! WF1(1:n(1)),WF2(1:n(2)),WF3(1:n(3)), WF3(1:n(4)) = input wfun arrays
! WF5(1:n(4))                                      = output wfun array
! vanishing components moved at the end of WF5 array and marked via %e = -1
! non-zero WF5 components ordered according to global helicity label
! WF1(h1), WF2(h2), WF3(h3), WF(h4) <-> WF5(h5) connection stored in
! table hi = t(i,h5) with i=1,2,3,4
! array sizes n(1), n(2), n(3), n(4), n(5) restricted to non-vanishing components
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(6), t(5,n(6))
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2)), WF3(n(3)), WF4(n(4)), WF5(n(5))
  type(wfun),        intent(out)   :: WF6(n(6))
  integer(intkind2) :: h1, h2, h3, h4, h5, h6, i, n2_in, n3_in, n4_in, n5_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(5)
  type(wfun)        :: WFaux

  if(ntry /= 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_vert6:')
    call ol_error(2,'ntry =' // to_string(ntry) // ' not allowed')
    call ol_fatal()
  end if

  ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n(1)
  do i = 1, n(1)
    if (WF1(i)%e == -1_intkind2) then
      h1 = i - 1
      exit
    end if
  end do
  do i = h1 + 1, n(1)
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                  &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
  end do
  n(1) = h1

  ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
  h2 = n(2)
  do i = 1, n(2)
    if (WF2(i)%e == -1_intkind2) then
      h2 = i - 1
      exit
    end if
  end do
  do i = h2 + 1, n(2)
    if (WF2(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h2, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
      call ol_fatal()
    end if
  end do
  n2_in = n(2)
  n(2) = h2

  ! sets n(3) = # of non-zero WF3 components and check that all zeros are at the end
  h3 = n(3)
  do i = 1, n(3)
    if (WF3(i)%e == -1_intkind2) then
      h3 = i - 1
      exit
    end if
  end do
  do i = h3 + 1, n(3)
    if (WF3(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h3, n(3), WF3(i)%e =' // to_string(i) // " " // to_string(h3) // " " //  &
                  &  to_string(n(3)) // " " //  to_string(WF3(i)%e))
      call ol_fatal()
    end if
  end do
  n3_in = n(3)
  n(3) = h3


  ! sets n(4) = # of non-zero WF4 components and check that all zeros are at the end
  h4 = n(4)
  do i = 1, n(4)
    if (WF4(i)%e == -1_intkind2) then
      h4 = i - 1
      exit
    end if
  end do
  do i = h4 + 1, n(4)
    if (WF4(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h4, n(4), WF4(i)%e =' // to_string(i) // " " // to_string(h4) // " " //  &
                  &  to_string(n(4)) // " " //  to_string(WF4(i)%e))
      call ol_fatal()
    end if
  end do
  n4_in = n(4)
  n(4) = h4


  ! sets n(5) = # of non-zero WF5 components and check that all zeros are at the end
  h5 = n(5)
  do i = 1, n(5)
    if (WF5(i)%e == -1_intkind2) then
      h5 = i - 1
      exit
    end if
  end do
  do i = h5 + 1, n(5)
    if (WF5(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert6:')
      call ol_error(2,'i, h5, n(5), WF5(i)%e =' // to_string(i) // " " // to_string(h5) // " " //  &
                  &  to_string(n(5)) // " " //  to_string(WF5(i)%e))
      call ol_fatal()
    end if
  end do
  n5_in = n(5)
  n(5) = h5


  ! restrict I/O table to non-vanishing helicity configurations
  i  = 0 ! index of WF6 states
  h6 = 0 ! index of non-zero WF5 states
  do h1 = 1, n(1)
    do h2 = 1, n(2)
      do h3 = 1, n(3)
        do h4 = 1, n(4)
          do h5 = 1, n(5)
            i  = h5 + n5_in * (h4-1 + n4_in * (h3-1 + n3_in * (h2-1 + n2_in * (h1-1)))) ! assumes input table in standard inititalisation form
            if (all(WF6(i)%j == 0)) cycle ! skips vanishing WF5 components
            h6 = h6 + 1
            t(1,h6) = h1
            t(2,h6) = h2
            t(3,h6) = h3
            t(4,h6) = h4
            t(5,h6) = h5
            WF6(h6)%e = WF1(h1)%e + WF2(h2)%e + WF3(h3)%e + WF4(h4)%e + WF5(h5)%e ! additive helicity label for outgoing states
            if (h6 == i) cycle
            WF6(h6)%j = WF6(i)%j ! shifts non-zero components to first part of WF6 array
            WF6(h6)%h = WF6(i)%h ! shifts non-zero components to first part of WF6 array
          end do
        end do
      end do
    end do
  end do

  ! put vanishing components at the end of WF4 array
  do i = h6 + 1, n(6)
    WF6(i)%j = 0
    WF6(i)%h = B"00"
    WF6(i)%e = -1_intkind2 ! marker for vanising helicity states
  end do

  ! sets n(6) = # of non-vanishing I/O helicity configurations
  n(6) = h6

  ! sort non-vanishing WF6(1:h6) components and adapt table accordingly
  do i = 1, n(6) - 1
    emin = WF6(i)%e
    imin = i
    do iq = i + 1, n(6)              ! look for component with helicity < emin
      if (WF6(iq)%e >= emin) cycle
      emin = WF6(iq)%e
      imin = iq
    end do
    if (imin > i) then
       WFaux     = WF6(i) ! flip WF6(i) <-> WF6(imin) components
       WF6(i)    = WF6(imin)
       WF6(imin) = WFaux
       taux        = t(1:5,i) ! same for table
       t(1:5,i)    = t(1:5,imin)
       t(1:5,imin) = taux
    end if
  end do

end subroutine helbookkeeping_vert6


! **********************************************************************
subroutine helbookkeeping_vert7(ntry, WF1, WF2, WF3, WF4, WF5, WF6, WF7, n, t)
! ----------------------------------------------------------------------
! WF1(1:n(1)),..., WF6(1:n(4))                     = input wfun arrays
! WF7(1:n(4))                                      = output wfun array
! vanishing components moved at the end of WF6 array and marked via %e = -1
! non-zero WF5 components ordered according to global helicity label
! WF1(h1), WF2(h2), WF3(h3), WF4(h4), WF(h5) <-> WF6(h6) connection stored in
! table hi = t(i,h6) with i=1,2,3,4,5
! array sizes n(1), n(2), n(3), n(4), n(5), n(6) restricted to non-vanishing components
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(7), t(6,n(7))
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2)), WF3(n(3)), WF4(n(4)), WF5(n(5)), WF6(n(6))
  type(wfun),        intent(out)   :: WF7(n(7))
  integer(intkind2) :: h1, h2, h3, h4, h5, h6, h7, i, n2_in, n3_in, n4_in, n5_in, n6_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(6)
  type(wfun)        :: WFaux

  if(ntry /= 1) then ! the following operations input table t in initialisation form
    call ol_error(2,'in subroutine helbookkeeping_vert7:')
    call ol_error(2,'ntry =' // to_string(ntry) // ' not allowed')
    call ol_fatal()
  end if
  ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
  h1 = n(1)
  do i = 1, n(1)
    if (WF1(i)%e == -1_intkind2) then
      h1 = i - 1
      exit
    end if
  end do
  do i = h1 + 1, n(1)
    if (WF1(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                  &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
      call ol_fatal()
    end if
  end do
  n(1) = h1

  ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
  h2 = n(2)
  do i = 1, n(2)
    if (WF2(i)%e == -1_intkind2) then
      h2 = i - 1
      exit
    end if
  end do
  do i = h2 + 1, n(2)
    if (WF2(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h2, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
      call ol_fatal()
    end if
  end do
  n2_in = n(2)
  n(2) = h2

  ! sets n(3) = # of non-zero WF3 components and check that all zeros are at the end
  h3 = n(3)
  do i = 1, n(3)
    if (WF3(i)%e == -1_intkind2) then
      h3 = i - 1
      exit
    end if
  end do
  do i = h3 + 1, n(3)
    if (WF3(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h3, n(3), WF3(i)%e =' // to_string(i) // " " // to_string(h3) // " " //  &
                  &  to_string(n(3)) // " " //  to_string(WF3(i)%e))
      call ol_fatal()
    end if
  end do
  n3_in = n(3)
  n(3) = h3


  ! sets n(4) = # of non-zero WF4 components and check that all zeros are at the end
  h4 = n(4)
  do i = 1, n(4)
    if (WF4(i)%e == -1_intkind2) then
      h4 = i - 1
      exit
    end if
  end do
  do i = h4 + 1, n(4)
    if (WF4(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert5:')
      call ol_error(2,'i, h4, n(4), WF4(i)%e =' // to_string(i) // " " // to_string(h4) // " " //  &
                  &  to_string(n(4)) // " " //  to_string(WF4(i)%e))
      call ol_fatal()
    end if
  end do
  n4_in = n(4)
  n(4) = h4


  ! sets n(5) = # of non-zero WF5 components and check that all zeros are at the end
  h5 = n(5)
  do i = 1, n(5)
    if (WF5(i)%e == -1_intkind2) then
      h5 = i - 1
      exit
    end if
  end do
  do i = h5 + 1, n(5)
    if (WF5(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert7:')
      call ol_error(2,'i, h5, n(5), WF5(i)%e =' // to_string(i) // " " // to_string(h5) // " " //  &
                  &  to_string(n(5)) // " " //  to_string(WF5(i)%e))
      call ol_fatal()
    end if
  end do
  n5_in = n(5)
  n(5) = h5


  ! sets n(6) = # of non-zero WF6 components and check that all zeros are at the end
  h6 = n(6)
  do i = 1, n(6)
    if (WF6(i)%e == -1_intkind2) then
      h6 = i - 1
      exit
    end if
  end do
  do i = h6 + 1, n(6)
    if (WF6(i)%e /= -1_intkind2) then
      call ol_error(2,'in subroutine helbookkeeping_vert7:')
      call ol_error(2,'i, h6, n(6), WF6(i)%e =' // to_string(i) // " " // to_string(h6) // " " //  &
                  &  to_string(n(6)) // " " //  to_string(WF6(i)%e))
      call ol_fatal()
    end if
  end do
  n6_in = n(6)
  n(6) = h6

  ! restrict I/O table to non-vanishing helicity configurations
  i  = 0 ! index of WF7 states
  h7 = 0 ! index of non-zero WF6 states
  do h1 = 1, n(1)
    do h2 = 1, n(2)
      do h3 = 1, n(3)
        do h4 = 1, n(4)
          do h5 = 1, n(5)
            do h6 = 1, n(6)
              i  = h6 + n6_in * (h5-1 + n5_in * (h4-1 + n4_in * (h3-1 + n3_in * (h2-1 + n2_in * (h1-1))))) ! assumes input table in standard inititalisation form
              if (all(WF7(i)%j == 0)) cycle ! skips vanishing WF5 components
              h7 = h7 + 1
              t(1,h7) = h1
              t(2,h7) = h2
              t(3,h7) = h3
              t(4,h7) = h4
              t(5,h7) = h5
              t(6,h7) = h6
              WF7(h7)%e = WF1(h1)%e + WF2(h2)%e + WF3(h3)%e + WF4(h4)%e + WF5(h5)%e + WF6(h6)%e ! additive helicity label for outgoing states
              if (h7 == i) cycle
              WF7(h7)%j = WF7(i)%j ! shifts non-zero components to first part of WF7 array
              WF7(h7)%h = WF7(i)%h ! shifts non-zero components to first part of WF7 array
            end do
          end do
        end do
      end do
    end do
  end do

  ! put vanishing components at the end of WF7 array
  do i = h7 + 1, n(7)
    WF7(i)%j = 0
    WF7(i)%h = B"00"
    WF7(i)%e = -1_intkind2 ! marker for vanising helicity states
  end do

  ! sets n(7) = # of non-vanishing I/O helicity configurations
  n(7) = h7

  ! sort non-vanishing WF5(1:h5) components and adapt table accordingly
  do i = 1, n(7) - 1
    emin = WF7(i)%e
    imin = i
    do iq = i + 1, n(7)              ! look for component with helicity < emin
      if (WF7(iq)%e >= emin) cycle
      emin = WF7(iq)%e
      imin = iq
    end do
    if (imin > i) then
       WFaux     = WF7(i) ! flip WF7(i) <-> WF7(imin) components
       WF7(i)    = WF7(imin)
       WF7(imin) = WFaux
       taux        = t(1:6,i) ! same for table
       t(1:6,i)    = t(1:6,imin)
       t(1:6,imin) = taux
    end if
  end do

end subroutine helbookkeeping_vert7




! **********************************************************************
subroutine helbookkeeping_cont(nsync, WF1, WF2, cont, n, t, nhel)
! ----------------------------------------------------------------------
! WF1(1:n(1)), WF(1:n(2))  = input wfun arrays
! cont(1:nhel)             = contraction array
! i3(i1,i2) helicity table stored in  ik= t(k,i3) with k=1,2
! nsync = 1 :  filter and order n ,t, and cont according to WF3 entries
!     (n3 = max # hels  -> n3 = # nonvanishing hels for actual diag)
! nsync = 2 :  syncronise n, t, and cont using global helsync table cont(i)%s
!     (n3 = # nonvanishing hels for actual diag => n3 = global # nonvanishing hels)
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: wfun, polcont
  implicit none
  integer(intkind1), intent(in)    :: nsync
  integer(intkind2), intent(inout) :: n(3), nhel, t(2,nhel)
  type(wfun),        intent(in)    :: WF1(n(1)), WF2(n(2))
  type(polcont),     intent(inout) :: cont(nhel)
  integer(intkind2) :: h1, h2, h3, i, n2_in
  integer(intkind2) :: iq, imin, emin
  integer(intkind2) :: taux(2)
  type(polcont)     :: contaux

  if (n(3) > nhel) then
    call ol_error(2,'in subroutine helbookkeeping_cont:')
    call ol_error(2,'n(3) = ' // to_string(n(3)) // ' > nhel = ' // to_string(nhel) // ' not allowed')
    call ol_fatal()
  end if

  select case (nsync)

  case (1_intkind1) ! filter and order n, t, and cont according to WF3 entries
    ! sets n(1) = # of non-zero WF1 components and check that all zeros are at the end
    h1 = n(1)
    do i = 1, n(1)
      if (WF1(i)%e == -1_intkind2) then
        h1 = i-1
        exit
      end if
    end do
    do i = h1+1, n(1)
      if (WF1(i)%e /= -1_intkind2) then
        call ol_error(2,'in subroutine helbookkeeping_cont:')
        call ol_error(2,'i, h1, n(1), WF1(i)%e =' // to_string(i) // " " // to_string(h1) // " " //  &
                    &  to_string(n(1)) // " " //  to_string(WF1(i)%e))
        call ol_fatal()
      end if
    end do
    n(1) = h1

    ! sets n(2) = # of non-zero WF2 components and check that all zeros are at the end
    h2 = n(2)
    do i = 1, n(2)
      if (WF2(i)%e == -1_intkind2) then
        h2 = i-1
        exit
      end if
    end do
    do i = h2+1, n(2)
      if (WF2(i)%e /= -1_intkind2) then
        call ol_error(2,'in subroutine helbookkeeping_cont:')
        call ol_error(2,'i, h2, n(2), WF2(i)%e =' // to_string(i) // " " // to_string(h2) // " " //  &
                  &  to_string(n(2)) // " " //  to_string(WF2(i)%e))
        call ol_fatal()
      end if
    end do
    n2_in = n(2)
    n(2) = h2

    ! restrict I/O table to non-vanishing helicity configurations
    i  = 0 ! index of cont states
    h3 = 0 ! index of non-zero cont states
    do h1 = 1, n(1)
      do h2 = 1, n(2)
        i = (h1-1)*n2_in + h2              ! assumes input table in standard inititalisation form
        if (cont(i)%j == 0) cycle          ! skips vanishing cont components
        h3 = h3 + 1
        t(1,h3) = h1
        t(2,h3) = h2
        cont(h3)%e = WF1(h1)%e + WF2(h2)%e ! determines additive helicity label for outgoing states
        if (h3 == i) cycle
        cont(h3)%j = cont(i)%j             ! shifts non-zero components to first part of cont array
      end do
    end do

    ! put vanishing components at the end of cont array
    do i = h3+1, n(3)
      cont(i)%j = 0
      cont(i)%e = -1_intkind2      ! marker for vanising helicity states
    end do

    ! sets n(3) = # of non-vanishing I/O helicity configurations
    n(3) = h3

    ! sort non-vanishing cont(1:h3) components and adapt table accordingly
    do i = 1, n(3)-1
      emin = cont(i)%e
      imin = i
      do iq = i+1, n(3)            ! look for component with helicity < emin
        if (cont(iq)%e >= emin) cycle
        emin = cont(iq)%e
        imin = iq
      end do
      if (imin > i) then
        contaux     = cont(i)      ! flip cont(i) <-> cont(imin) components
        cont(i)     = cont(imin)
        cont(imin)  = contaux
        taux        = t(1:2,i)     ! same for table
        t(1:2,i)    = t(1:2,imin)
        t(1:2,imin) = taux
      end if
    end do

  case (2_intkind1) ! syncronise n, t, and cont using global helsync table cont(i)%s
    n(3) = nhel
    do i = Nhel, 1 , -1
      if (cont(i)%s == 0) then
        cont(i)%j = 0
        t(1:2,i)  = 0_intkind2
      else
        cont(i)%j = cont(cont(i)%s)%j
        t(1:2,i)  = t(1:2,cont(i)%s)
      end if
    end do

  case default
    call ol_error(2, 'in subroutine helbookkeeping_cont:')
    call ol_error(2, 'nsync = ' // to_string(nsync) // ' not allowed')
    call ol_fatal()
  end select

end subroutine helbookkeeping_cont



! **********************************************************************
subroutine helsync(nsync, cont, Nhel, Hel)
! ----------------------------------------------------------------------
! cont(1:Nhelmax,1:Ndiag) = amplitudes of all individual polarised diagrams
! Nhelmax                 = maximal number of helicity configurations
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  use ol_data_types_/**/REALKIND, only: polcont
  implicit none
  integer(intkind1), intent(in)    :: nsync
  type(polcont),     intent(inout) :: cont(:,:)
  integer(intkind2), intent(out)   :: Hel(size(cont,1))
  integer(intkind2), intent(inout) :: NHel
  logical :: nonvanishing(size(cont,1))
  integer :: Nhelmax, h, h0, i, ishift, Ndiag, n

  if (nsync /= 1) then
    call ol_error(2,'in subroutine helsync:')
    call ol_error(2,'nsync = ' // to_string(nsync) // ' not allowed')
    call ol_fatal()
  end if

  Ndiag = size(cont,2)   ! number of diagrams
  Nhelmax = size(cont,1) ! maximal number of helicity configurations
  nonvanishing(1:Nhelmax) = .false.

  ! sets nonvanishing(h) = .true. if configuration with helicity Hel = h-1 non vanishing
  ! Hel = physical helicty    :   0 <= Hel <= Nhelmax - 1
  ! h = Hel+1 = helicty index :   1 <=  h  <= Nhelmax

  nexthel: do h = 1, Nhelmax
    do h0 = 1, Nhelmax
      do n = 1, Ndiag
        if (cont(h0,n)%e == h-1) then
          nonvanishing(h) = .true.
          cycle nexthel
        end if
      end do
    end do
  end do nexthel

  ! extract subset Hel(1:Nhel) of non-vanishing helicity configurations
  Nhel = 0
  do h = 1, Nhelmax
    if (nonvanishing(h)) then
      NHel      = Nhel + 1
      Hel(NHel) = h - 1
    end if
  end do
  Hel(Nhel+1:Nhelmax) = -1 ! put vanishing helicity configurations at the end

  ! helicity-index table cont(i,n)%s to determine position of helicity Hel(i) in diagram n
  !   cont(cont(i,n)%s,n)%e = Hel(i)   when diagram n non-vanisging for helicity state Hel(i)
  !        cont(i,n)%s      = 0        when diagram n vanishes for helicity state Hel(i)
  ! NB: construction below exploits pre-ordering of cont array, i.e. cont(h2,n)%e > cont(h1,n)%e for h2 > h1

  do n = 1, Ndiag
    ishift = 0        ! incremental counter of absent Hel(:) elements in cont(:,n)%e
    do i = 1, Nhel    ! index to scan cont(1:Nhel,n)%e (sub)row
      if (cont(i-ishift,n)%e == Hel(i)) then
        cont(i,n)%s = i - ishift    ! position of Hel(i) in n-th row
      else
        cont(i,n)%s = 0             ! mark absence of Hel(i) in n-th row
        ishift = ishift + 1         ! avoids incrementing i - ishift
      end if
    end do
    cont(Nhel+1:Nhelmax,n)%s = -1  ! irrelevant positions <-> helicity states that never contribute
  end do

end subroutine helsync



! **********************************************************************
subroutine flip_phase(P, pol, MOM, omega)
! ----------------------------------------------------------------------
! P     = emitter 4-momentum (contravariant standard representation, real)
! pol   = +1/-1 gluon-emitter helicity
! MOM   = auxiliary momentum (contravariant standard representation, real)
! omega = (<epsilon(pol),MOM>/|<epsilon(pol),MOM>|)^2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3)
  integer,           intent(in)  :: pol
  real(REALKIND),    intent(in)  :: MOM(0:3)
  complex(REALKIND), intent(out) :: omega(2)
  complex(REALKIND) :: eps(4), MOM_LC(4)
  call wf_V_Std(P, 0._/**/REALKIND, pol, eps) ! light-cone polarisation vector
  call Std2LC_Rep(MOM, MOM_LC)

  ! Do not use h_contractions::cont_PP(eps,MOM_LC) to avoid cyclic dependencies
  omega(1) = eps(1)*MOM_LC(2) + eps(2)*MOM_LC(1) - eps(3)*MOM_LC(4) - eps(4)*MOM_LC(3)
  omega(1) = omega(1)/abs(omega(1))
  omega(1) = omega(1)*omega(1)
  omega(2) = conjg(omega(1))

end subroutine flip_phase

end module ol_helicity_bookkeeping_/**/REALKIND
