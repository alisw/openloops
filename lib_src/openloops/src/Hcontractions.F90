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


! **********************************************************************
module ol_hel_contractions_/**/REALKIND
! Routines to contract two wave functions
! - cont_PP: two 4-vectors which are not stored in a wfun object
! - cont_SS: two scalars
! - cont_CD: two ghosts (this is the same as cont_SS)
! - Hcont_VV: two vector bosons
! - Hcont_QA: fermion and anti-fermion
! - cont_EE: two sigma particles, expressed as four gluons
! **********************************************************************
  implicit none
  contains

! **********************************************************************
function cont_PP(A, B)
! Contraction of (complex) Lorentz vectors in light-cone representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_PP
  complex(REALKIND), intent(in) :: A(4), B(4)
  cont_PP = A(1)*B(2) + A(2)*B(1) - A(3)*B(4) - A(4)*B(3)
  cont_PP = 0.5_/**/REALKIND * cont_PP
end function cont_PP

!
! !!! To be tested
!
! ! **********************************************************************
! subroutine cont_CD(nsync, A, B, cont, n, t, nhel, den)
! ! Contraction of two ghosts (this is the same as cont_SS for two scalars).
! ! ----------------------------------------------------------------------
! ! A(1:n(1)), B(1:n(2)) = wfun type (%h component ignored)
! ! cont(1:nhel)         = contractions
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND, intkind1, intkind2
!   use ol_data_types_/**/REALKIND
!   use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
!   implicit none
!   integer(intkind1), intent(in)    :: nsync
!   integer(intkind2), intent(inout) :: nhel, n(3), t(2,nhel)
!   type(wfun),        intent(in)    :: A(n(1)), B(n(2))
!   type(polcont),     intent(inout) :: cont(nhel)
!   complex(REALKIND),       intent(in)    :: den
!   integer :: h
!
!   do h = 1, n(3)
!     if (t(1,h) == 0_intkind2) then
!       cont(h)%j = 0
!     else
!       cont(h)%j = A(t(1,h))%j(1) * B(t(2,h))%j(1) * den
!     end if
!   end do
!
!   if(nsync <= 2) call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)
!
! end subroutine cont_CD


! **********************************************************************
subroutine Hcont_SS(nsync, A, B, cont, n, t, nhel, den)
! Contraction of two scalar particles.
! ----------------------------------------------------------------------
! A(1:n(1)), B(1:n(2)) = wfun type (%h component ignored)
! cont(1:nhel)         = contractions
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)    :: nsync
  integer(intkind2), intent(inout) :: nhel, n(3), t(2,nhel)
  type(wfun),        intent(in)    :: A(n(1)), B(n(2))
  type(Hpolcont),    intent(inout):: cont(nhel)
  complex(REALKIND),       intent(in)    :: den
  integer :: h

  do h = 1, n(3)
    if (t(1,h) == 0_intkind2) then
      cont(h)%j = 0
    else
      cont(h)%j = A(t(1,h))%j(1) * B(t(2,h))%j(1) * den
    end if
  end do

  if(nsync <= 2) then
    do h=1, n(3)
      if (A(t(1,h))%hf==-1_intkind2 .OR. B(t(2,h))%hf==-1_intkind2) then
        cont(h)%hf = -1_intkind2
      else
        cont(h)%hf = A(t(1,h))%hf + B(t(2,h))%hf
      end if
    end do
    call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)
  end if

end subroutine Hcont_SS

! **********************************************************************
subroutine Hcont_VV(nsync, A, B, cont, n, t, nhel, den)
! Contraction of (complex) Lorentz vectors in light-cone representation
! as part of a wave function object
! **********************************************************************
! A(1:n(1)),B(1:n(2)) = wfun type (%h component ignored)
! cont(1:nhel)        = contractions
!************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)   :: nsync
  integer(intkind2), intent(inout):: n(3), nhel, t(2,nhel)
  type(wfun),        intent(in)   :: A(n(1)), B(n(2))
  type(Hpolcont),    intent(inout):: cont(nhel)
  complex(REALKIND),       intent(in)   :: den
  complex(REALKIND) :: denhalf
  integer     :: h

  denhalf = 0.5_/**/REALKIND * den

  do h = 1, n(3)
    if (t(1,h) == 0_intkind2) then
      cont(h)%j = 0
    else
      cont(h)%j = A(t(1,h))%j(1)*B(t(2,h))%j(2) + A(t(1,h))%j(2)*B(t(2,h))%j(1) &
                - A(t(1,h))%j(3)*B(t(2,h))%j(4) - A(t(1,h))%j(4)*B(t(2,h))%j(3)
      cont(h)%j = cont(h)%j * denhalf
    end if
  end do

  if(nsync <= 2) then
    do h=1, n(3)
      if (A(t(1,h))%hf==-1_intkind2 .OR. B(t(2,h))%hf==-1_intkind2) then
        cont(h)%hf = -1_intkind2
      else
        cont(h)%hf = A(t(1,h))%hf + B(t(2,h))%hf
      end if
    end do
    call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)
  end if

end subroutine Hcont_VV

! **********************************************************************
subroutine Hcont_QA(nsync, Q, A, cont, n, t, nhel, den)
! Contraction of quark-antiquark currents
! (standard scalar product of complex 4-vectors)
! **********************************************************************
! Q(1:n(1))      = quark current & helicity
! A(1:n(2))      = anti-quark current & helicity
! cont(1:nhel)   = contractions
!************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)   :: nsync
  integer(intkind2), intent(inout):: n(3), nhel, t(2,nhel)
  type(wfun),        intent(in)   :: Q(n(1)), A(n(2))
  type(Hpolcont),     intent(inout):: cont(nhel)
  complex(REALKIND),       intent(in)   :: den
  integer        :: h
  integer, save  :: counter =1

  do h = 1, n(3)

    if (t(1,h) == 0_intkind2) then
      cont(h)%j = 0
    else
    select case (iand(Q(t(1,h))%h,A(t(2,h))%h))
    case(B01)
      cont(h)%j = Q(t(1,h))%j(3)*A(t(2,h))%j(3) + Q(t(1,h))%j(4)*A(t(2,h))%j(4)
      cont(h)%j = cont(h)%j * den
    case(B10)
      cont(h)%j = Q(t(1,h))%j(1)*A(t(2,h))%j(1) + Q(t(1,h))%j(2)*A(t(2,h))%j(2)
      cont(h)%j = cont(h)%j * den
    case(B00)
      cont(h)%j = 0
    case default
      cont(h)%j = Q(t(1,h))%j(1)*A(t(2,h))%j(1) + Q(t(1,h))%j(2)*A(t(2,h))%j(2) &
                + Q(t(1,h))%j(3)*A(t(2,h))%j(3) + Q(t(1,h))%j(4)*A(t(2,h))%j(4)
      cont(h)%j = cont(h)%j * den
    end select
    end if
  end do

  if(nsync <= 2) then
    do h = 1, n(3)
      if (Q(t(1,h))%hf==-1_intkind2 .OR. A(t(2,h))%hf==-1_intkind2) then
        cont(h)%hf = -1_intkind2
      else
        cont(h)%hf = Q(t(1,h))%hf + A(t(2,h))%hf
      end if
  end do
    call helbookkeeping_cont(nsync, Q, A, cont, n, t, nhel)
  end if

end subroutine Hcont_QA


end module ol_hel_contractions_/**/REALKIND
