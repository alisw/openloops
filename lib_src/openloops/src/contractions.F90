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


module ol_contractions_/**/REALKIND
  implicit none
  contains

! **********************************************************************
function cont_LL(A, B)
! Contraction of (real) Lorentz vectors in standard representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  real(REALKIND) :: cont_LL
  real(REALKIND), intent(in) :: A(0:3), B(0:3)
  cont_LL = A(0)*B(0) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
end function cont_LL

! **********************************************************************
function cont_L(A)
! Contraction of a (real) Lorentz vector in standard representation with itself
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  real(REALKIND) :: cont_L
  real(REALKIND), intent(in) :: A(0:3)
  cont_L = A(0)*A(0) - A(1)*A(1) - A(2)*A(2) - A(3)*A(3)
end function cont_L

! **********************************************************************
function cont_VV(A, B)
! Contraction of (complex) Lorentz vectors in light-cone representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_VV
  complex(REALKIND), intent(in) :: A(4), B(4)
  cont_VV = A(1)*B(2) + A(2)*B(1) - A(3)*B(4) - A(4)*B(3)
  cont_VV = 0.5_/**/REALKIND * cont_VV
end function cont_VV

! **********************************************************************
function cont_V(A)
! Contraction of a (complex) Lorentz vector in light-cone representation with itself
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_V
  complex(REALKIND), intent(in) :: A(4)
  cont_V = A(1)*A(2) + A(2)*A(1) - A(3)*A(4) - A(4)*A(3)
  cont_V = 0.5_/**/REALKIND * cont_V
end function cont_V

! **********************************************************************
function contsum_VVV(A1, A2, B)
! Contraction of (complex) Lorentz vectors in light-cone representation
! contsum_VVV = [A1+A2](mu) * B(mu)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: contsum_VVV
  complex(REALKIND), intent(in) :: A1(4), A2(4), B(4)
  contsum_VVV = (A1(1)+A2(1))*B(2) + (A1(2)+A2(2))*B(1) - (A1(3)+A2(3))*B(4) - (A1(4)+A2(4))*B(3)
  contsum_VVV = 0.5_/**/REALKIND * contsum_VVV
end function contsum_VVV

! **********************************************************************
function cont_EE(A, B, C, D)
! Contraction of two sigma particles in a four gluon vertex.
! Each sigma is represented by its two incoming gluon (complex) Lorentz
! vectors in light-cone representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_EE
  complex(REALKIND), intent(in) :: A(4), B(4), C(4), D(4)
  cont_EE = cont_VV(A,C)*cont_VV(B,D) - cont_VV(B,C)*cont_VV(A,D)
end function cont_EE

! **********************************************************************
function cont_SS(A, B)
! Contraction of two scalar particles.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_SS
  complex(REALKIND), intent(in) :: A(4), B(4)
  cont_SS = A(1) * B(1)
end function cont_SS

! **********************************************************************
function cont_CD(A, B)
! Contraction of two scalar particles.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_CD
  complex(REALKIND), intent(in) :: A(4), B(4)
  cont_CD = A(1) * B(1)
end function cont_CD

! **********************************************************************
function cont_QA(A, B)
! Contraction of quark-antiquark currents
! (standard scalar product of complex 4-vectors)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_QA
  complex(REALKIND), intent(in) :: A(4), B(4)
  cont_QA = A(1)*B(1) + A(2)*B(2) + A(3)*B(3) + A(4)*B(4)
end function cont_QA

! ! **********************************************************************
! subroutine adjoint_DIRAC(J_Q, J_A)
! ! Dirac conjugation in the Weyl representation
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   implicit none
!   complex(REALKIND) :: J_Q(4), J_A(4)
!   J_A(1) = -conjg(J_Q(3))
!   J_A(2) = -conjg(J_Q(4))
!   J_A(3) = -conjg(J_Q(1))
!   J_A(4) = -conjg(J_Q(2))
! end subroutine adjoint_DIRAC


! **********************************************************************
subroutine cont_EpVVV(B, C, D, Aout)
! Aout^i = I * gLC^i_a * ep_{abcd} * B^b*C^c*D^d
! where gLC is the metric in light-cone representation.
! Aout, B, C, D are in contravariant ligh-cone representation.
! Factors are so that X.gLC.Aout = cont_VV(X,Aout) = ep_{ijkl} * x^i*b^j*c^k*d^l
! where x,b,c,d are the vectors X,B,C,D in standard representation.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: B(4), C(4), D(4)
  complex(REALKIND), intent(out) :: Aout(4)
  complex(REALKIND)              :: C1D4, C1D3, C1D2, C2D3, C2D4, C3D4
  ! all vectors in standard representation
  ! aout(1) = b(2) * (c(3)*d(4) - c(4)*d(3)) + b(3) * (c(4)*d(2) - c(2)*d(4)) + b(4) * (c(2)*d(3) - c(3)*d(2))
  ! aout(2) = b(1) * (c(3)*d(4) - c(4)*d(3)) + b(3) * (c(4)*d(1) - c(1)*d(4)) + b(4) * (c(1)*d(3) - c(3)*d(1))
  ! aout(3) = b(1) * (c(4)*d(2) - c(2)*d(4)) + b(2) * (c(1)*d(4) - c(4)*d(1)) + b(4) * (c(2)*d(1) - c(1)*d(2))
  ! aout(4) = b(1) * (c(2)*d(3) - c(3)*d(2)) + b(2) * (c(3)*d(1) - c(1)*d(3)) + b(3) * (c(1)*d(2) - c(2)*d(1))
  ! light-cone representation
  ! Aout(1) = i/2 * (B(1) * (C(4)*D(3) - C(3)*D(4)) + B(3) * (C(1)*D(4) - C(4)*D(1)) + B(4) * (C(3)*D(1) - C(1)*D(3)))
  ! Aout(2) = i/2 * (B(2) * (C(3)*D(4) - C(4)*D(3)) + B(3) * (C(4)*D(2) - C(2)*D(4)) + B(4) * (C(2)*D(3) - C(3)*D(2)))
  ! Aout(3) = i/2 * (B(1) * (C(2)*D(3) - C(3)*D(2)) + B(2) * (C(3)*D(1) - C(1)*D(3)) + B(3) * (C(1)*D(2) - C(2)*D(1)))
  ! Aout(4) = i/2 * (B(1) * (C(4)*D(2) - C(2)*D(4)) + B(2) * (C(1)*D(4) - C(4)*D(1)) + B(4) * (C(2)*D(1) - C(1)*D(2)))
  C1D2    = C(1)*D(2) - C(2)*D(1)
  C1D3    = C(1)*D(3) - C(3)*D(1)
  C1D4    = C(1)*D(4) - C(4)*D(1)
  C2D3    = C(2)*D(3) - C(3)*D(2)
  C2D4    = C(2)*D(4) - C(4)*D(2)
  C3D4    = C(3)*D(4) - C(4)*D(3)
  Aout(1) = - B(1) * C3D4 - B(4) * C1D3 + B(3) * C1D4
  Aout(2) =   B(2) * C3D4 - B(3) * C2D4 + B(4) * C2D3
  Aout(3) =   B(3) * C1D2 + B(1) * C2D3 - B(2) * C1D3
  Aout(4) = - B(4) * C1D2 + B(2) * C1D4 - B(1) * C2D4
  Aout    = (0._/**/REALKIND, 0.5_/**/REALKIND) * Aout
end subroutine cont_EpVVV

end module ol_contractions_/**/REALKIND



! **********************************************************************
module ol_s_contractions_/**/REALKIND
! Routines to contract two wave functions
! - cont_PP: two 4-vectors which are not stored in a wfun object
! - cont_SS: two scalars
! - cont_CD: two ghosts (this is the same as cont_SS)
! - cont_VV: two vector bosons
! - cont_QA: fermion and anti-fermion
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


! **********************************************************************
function cont_SS(A, B)
! Contraction of two scalar particles.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND
  implicit none
  complex(REALKIND) :: cont_SS
  type(wfun), intent(in)  :: A, B
  cont_SS = A%j(1) * B%j(1)
end function cont_SS


! **********************************************************************
function cont_CD(A, B)
! Contraction of two ghosts (this is the same as cont_SS for two scalars).
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND
  implicit none
  complex(REALKIND) :: cont_CD
  type(wfun), intent(in)  :: A, B
  cont_CD = A%j(1) * B%j(1)
end function cont_CD


! **********************************************************************
function cont_VV(A, B)
! Contraction of (complex) Lorentz vectors in light-cone representation
! as part of a wave function oblect
! **********************************************************************
! A,B = wfun type (%h component ignored)
!************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND
  implicit none
  complex(REALKIND)             :: cont_VV
  type(wfun), intent(in)  :: A, B
  cont_VV = A%j(1)*B%j(2) + A%j(2)*B%j(1) - A%j(3)*B%j(4) - A%j(4)*B%j(3)
  cont_VV = 0.5_/**/REALKIND * cont_VV
end function cont_VV


! **********************************************************************
function cont_EE(A, B, C, D)
! Contraction of two sigma particles in a four gluon vertex.
! Each sigma is represented by its two incoming gluon (complex) Lorentz
! vectors in light-cone representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND
  implicit none
  complex(REALKIND) :: cont_EE
  type(wfun), intent(in)  :: A, B, C, D
  cont_EE = cont_VV(A,C)*cont_VV(B,D) - cont_VV(B,C)*cont_VV(A,D)
end function cont_EE


! **********************************************************************
function cont_QA(Q, A)
! Contraction of quark-antiquark currents
! (standard scalar product of complex 4-vectors)
! **********************************************************************
! Q      = quark current & helicity
! A      = anti-quark current & helicity
!-----------------------------------------------------------------------
! X%j    = (0, 0, 0, 0)  (0, 0, J3,J4)  (J1,J2, 0, 0)  (J1,J2,J3,J4)
! X%h    =    B"00"          B"01"         B"10"          B"11"
!************************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND
  implicit none
  complex(REALKIND) :: cont_QA
  type(wfun), intent(in)  :: Q, A

  select case (iand(Q%h,A%h))
  case(B01)
    cont_QA = Q%j(3)*A%j(3) + Q%j(4)*A%j(4)
  case(B10)
    cont_QA = Q%j(1)*A%j(1) + Q%j(2)*A%j(2)
  case(B00)
    cont_QA = 0
  case default
    cont_QA = Q%j(1)*A%j(1) + Q%j(2)*A%j(2) + Q%j(3)*A%j(3) + Q%j(4)*A%j(4)
  end select

end function cont_QA

end module ol_s_contractions_/**/REALKIND



! **********************************************************************
module ol_h_contractions_/**/REALKIND
! Routines to contract two wave functions
! - cont_PP: two 4-vectors which are not stored in a wfun object
! - cont_SS: two scalars
! - cont_CD: two ghosts (this is the same as cont_SS)
! - cont_VV: two vector bosons
! - cont_QA: fermion and anti-fermion
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


! **********************************************************************
subroutine cont_SS(nsync, A, B, cont, n, t, nhel, den)
! Contraction of two scalar particles.
! ----------------------------------------------------------------------
! A(1:n(1)), B(1:n(2)) = wfun type (%h component ignored)
! cont(1:nhel)         = contractions
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)    :: nsync
  integer(intkind2), intent(inout) :: nhel, n(3), t(2,nhel)
  type(wfun),        intent(in)    :: A(n(1)), B(n(2))
  type(polcont),     intent(inout) :: cont(nhel)
  complex(REALKIND),       intent(in)    :: den
  integer :: h

  do h = 1, n(3)
    if (t(1,h) == 0_intkind2) then
      cont(h)%j = 0
    else
      cont(h)%j = A(t(1,h))%j(1) * B(t(2,h))%j(1) * den
    end if
  end do

  if(nsync <= 2) call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)

end subroutine cont_SS


! **********************************************************************
subroutine cont_CD(nsync, A, B, cont, n, t, nhel, den)
! Contraction of two ghosts (this is the same as cont_SS for two scalars).
! ----------------------------------------------------------------------
! A(1:n(1)), B(1:n(2)) = wfun type (%h component ignored)
! cont(1:nhel)         = contractions
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)    :: nsync
  integer(intkind2), intent(inout) :: nhel, n(3), t(2,nhel)
  type(wfun),        intent(in)    :: A(n(1)), B(n(2))
  type(polcont),     intent(inout) :: cont(nhel)
  complex(REALKIND),       intent(in)    :: den
  integer :: h

  do h = 1, n(3)
    if (t(1,h) == 0_intkind2) then
      cont(h)%j = 0
    else
      cont(h)%j = A(t(1,h))%j(1) * B(t(2,h))%j(1) * den
    end if
  end do

  if(nsync <= 2) call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)

end subroutine cont_CD


! **********************************************************************
subroutine cont_VV(nsync, A, B, cont, n, t, nhel, den)
! Contraction of (complex) Lorentz vectors in light-cone representation
! as part of a wave function object
! **********************************************************************
! A(1:n(1)),B(1:n(2)) = wfun type (%h component ignored)
! cont(1:nhel)        = contractions
!************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)   :: nsync
  integer(intkind2), intent(inout):: n(3), nhel, t(2,nhel)
  type(wfun),        intent(in)   :: A(n(1)), B(n(2))
  type(polcont),     intent(inout):: cont(nhel)
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

  if(nsync <= 2) call helbookkeeping_cont(nsync, A, B, cont, n, t, nhel)

end subroutine cont_VV


! **********************************************************************
subroutine cont_QA(nsync, Q, A, cont, n, t, nhel, den)
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
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_cont
  implicit none
  integer(intkind1), intent(in)   :: nsync
  integer(intkind2), intent(inout):: n(3), nhel, t(2,nhel)
  type(wfun),        intent(in)   :: Q(n(1)), A(n(2))
  type(polcont),     intent(inout):: cont(nhel)
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

  if(nsync <= 2) call helbookkeeping_cont(nsync, Q, A, cont, n, t, nhel)

end subroutine cont_QA




! **********************************************************************
subroutine cont_EpPPP(B, C, D, Aout)
! Aout^i = I * gLC^i_a * ep_{abcd} * B^b*C^c*D^d
! where gLC is the metric in light-cone representation.
! Aout, B, C, D are in contravariant ligh-cone representation.
! Factors are so that X.gLC.Aout = cont_VV(X,Aout) = ep_{ijkl} * x^i*b^j*c^k*d^l
! where x,b,c,d are the vectors X,B,C,D in standard representation.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: B(4), C(4), D(4)
  complex(REALKIND), intent(out) :: Aout(4)
  complex(REALKIND)              :: C1D4, C1D3, C1D2, C2D3, C2D4, C3D4
  C1D2    = C(1)*D(2) - C(2)*D(1)
  C1D3    = C(1)*D(3) - C(3)*D(1)
  C1D4    = C(1)*D(4) - C(4)*D(1)
  C2D3    = C(2)*D(3) - C(3)*D(2)
  C2D4    = C(2)*D(4) - C(4)*D(2)
  C3D4    = C(3)*D(4) - C(4)*D(3)
  Aout(1) = - B(1) * C3D4 - B(4) * C1D3 + B(3) * C1D4
  Aout(2) =   B(2) * C3D4 - B(3) * C2D4 + B(4) * C2D3
  Aout(3) =   B(3) * C1D2 + B(1) * C2D3 - B(2) * C1D3
  Aout(4) = - B(4) * C1D2 + B(2) * C1D4 - B(1) * C2D4
  Aout    = (0._/**/REALKIND, 0.5_/**/REALKIND) * Aout
end subroutine cont_EpPPP



end module ol_h_contractions_/**/REALKIND
