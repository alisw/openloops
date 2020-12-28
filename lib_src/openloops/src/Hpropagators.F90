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


! *****************************************************************************
module ol_hel_propagators_/**/REALKIND
! Routines for dressing a wave function with a propagator
! - prop_Q_A: fermion
! - prop_A_Q: anti-fermion
! - prop_W_W: massive vector boson in unitary gauge
! *****************************************************************************
  use KIND_TYPES_BW
  implicit none
  contains

! *****************************************************************************
subroutine prop_Q_A(ntry, Q, mom, M, mQ, Q_out, n)
! dressing quark current with propagator
! -----------------------------------------------------------------------------
! Q(1:n)       = quark wfun
! K(4)         = incoming momentum (light-cone rep)
! M            = complex mass
! mQ           = 0/1 for massless/massive propagator
! Q_out(1:n)   = outgoing dressed quark current without I/(K^2-M^2)
! Q_out(h)%j(i)= [slash(K)+M](i,j) * Q%j(j)
! *****************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4, get_LC_mass2
  use ol_momenta_decl_/**/REALKIND, only: collconf, softconf
  use ol_parameters_decl_/**/DREALKIND, only: ir_hacks
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom
  integer(intkind2), intent(inout) :: n
  type(wfun),        intent(in)    :: Q(n)
  complex(REALKIND), intent(in)    :: M
  integer(intkind1), intent(in)    :: mQ
  type(wfun),        intent(inout) :: Q_out(n)
  complex(REALKIND) :: K(4)
  integer :: h,momh

  ! store wavefunction for later use
  if (ir_hacks .and. (softconf .ne. 0 .or. collconf .ne. 0)) then
    ! check if a soft particle is attached to an external particle
    if ((mom .eq. collconf) .or. &
        (iand(mom,softconf) .ne. 0 .and. iand(mom-softconf,mom-softconf-1) .eq. 0)) then
       !currently only supported for massless particles
      !if (get_LC_mass2(mom-softconf) .eq. 0) then
      do h = 1, n
        if (associated(Q_out(h)%j_prev)) deallocate(Q_out(h)%j_prev)
        allocate(Q_out(h)%j_prev(4))
        Q_out(h)%j_prev(:) = Q(h)%j(:)
      end do
      !end if
    end if
  end if

  K = get_LC_4(mom)
  do h = 1, n

    select case (Q(h)%h)

    case (B01)
      Q_out(h)%j(1) = - K(2)*Q(h)%j(3) + K(4)*Q(h)%j(4)
      Q_out(h)%j(2) = - K(1)*Q(h)%j(4) + K(3)*Q(h)%j(3)
      if (mQ == 0) then
        Q_out(h)%j(3:4) = 0
        Q_out(h)%h      = B10
      else
        Q_out(h)%j(3:4) = M*Q(h)%j(3:4)
        Q_out(h)%h      = B11
      end if

    case (B10)
      Q_out(h)%j(3) = - K(1)*Q(h)%j(1) - K(4)*Q(h)%j(2)
      Q_out(h)%j(4) = - K(2)*Q(h)%j(2) - K(3)*Q(h)%j(1)
      if (mQ == 0) then
        Q_out(h)%j(1:2) = 0
        Q_out(h)%h      = B01
      else
        Q_out(h)%j(1:2) = M*Q(h)%j(1:2)
        Q_out(h)%h      = B11
      end if

    case (B00)
      Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      Q_out(h)%h = B00

    case default
      if (mQ == 0) then
        Q_out(h)%j(1) = - K(2)*Q(h)%j(3) + K(4)*Q(h)%j(4)
        Q_out(h)%j(2) = - K(1)*Q(h)%j(4) + K(3)*Q(h)%j(3)
        Q_out(h)%j(3) = - K(1)*Q(h)%j(1) - K(4)*Q(h)%j(2)
        Q_out(h)%j(4) = - K(2)*Q(h)%j(2) - K(3)*Q(h)%j(1)
      else
        Q_out(h)%j(1) = - K(2)*Q(h)%j(3) + K(4)*Q(h)%j(4) + M*Q(h)%j(1)
        Q_out(h)%j(2) = - K(1)*Q(h)%j(4) + K(3)*Q(h)%j(3) + M*Q(h)%j(2)
        Q_out(h)%j(3) = - K(1)*Q(h)%j(1) - K(4)*Q(h)%j(2) + M*Q(h)%j(3)
        Q_out(h)%j(4) = - K(2)*Q(h)%j(2) - K(3)*Q(h)%j(1) + M*Q(h)%j(4)
      end if
      Q_out(h)%h      = B11

    end select

  end do

  if (ntry == 1) then
    Q_out(:)%n_part = Q(1)%n_part
    Q_out(:)%t = Q(1)%t
    Q_out(:)%hf = Q(:)%hf
    call helbookkeeping_prop(ntry, Q, Q_out, n)
  end if

end subroutine prop_Q_A


! *****************************************************************************
subroutine prop_A_Q(ntry, A, mom, M, mA, A_out, n)
! dressing anti-quark current with propagator
! -----------------------------------------------------------------------------
! A(1:n)    = anti-quark wfun
! K(4)      = incoming momentum (light-cone rep)
! M         = complex mass
! mA        = 0/1 for massless/massive propagator
! A_out(1:n = outgoing dressed anti-quark current without I/(K^2-M^2)
! A_out(h)%j(i)= A(h)%j(j) * [-slash(K)+M](j,i)
! *****************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4, get_LC_mass2
  use ol_momenta_decl_/**/REALKIND, only: collconf, softconf
  use ol_parameters_decl_/**/DREALKIND, only: ir_hacks
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom
  integer(intkind2), intent(inout) :: n
  type(wfun),        intent(in)    :: A(n)
  complex(REALKIND), intent(in)    :: M
  integer(intkind1), intent(in)    :: mA
  type(wfun),        intent(inout) :: A_out(n)
  complex(REALKIND) :: K(4)
  integer :: h,momh

  ! store wavefunction for later use
  if (ir_hacks .and. (softconf .ne. 0 .or. collconf .ne. 0)) then
    ! check if a soft particle is attached to an external particle
    if ((mom .eq. collconf) .or. &
        (iand(mom,softconf) .ne. 0 .and. iand(mom-softconf,mom-softconf-1) .eq. 0)) then
       !currently only supported for massless particles
      !if (get_LC_mass2(mom-softconf) .eq. 0) then
      do h = 1, n
        if (associated(A_out(h)%j_prev)) deallocate(A_out(h)%j_prev)
        allocate(A_out(h)%j_prev(4))
        A_out(h)%j_prev(:) = A(h)%j(:)
      end do
      !end if
    end if
  end if

  K = get_LC_4(mom)
  do h = 1, n
    select case (A(h)%h)

    case (B01)
      A_out(h)%j(1) = + K(1)*A(h)%j(3) + K(3)*A(h)%j(4)
      A_out(h)%j(2) = + K(2)*A(h)%j(4) + K(4)*A(h)%j(3)
      if (mA == 0) then
        A_out(h)%j(3:4) = 0
        A_out(h)%h      = B10
      else
        A_out(h)%j(3:4) = M*A(h)%j(3:4)
        A_out(h)%h      = B11
      end if

    case (B10)
      A_out(h)%j(3) = + K(2)*A(h)%j(1) - K(3)*A(h)%j(2)
      A_out(h)%j(4) = + K(1)*A(h)%j(2) - K(4)*A(h)%j(1)
      if (mA == 0) then
        A_out(h)%j(1:2) = 0
        A_out(h)%h      = B01
      else
        A_out(h)%j(1:2) = M*A(h)%j(1:2)
        A_out(h)%h      = B11
      end if

    case (B00)
      A_out(h)%j = 0  ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    case default
      if (mA == 0) then
        A_out(h)%j(1) = + K(1)*A(h)%j(3) + K(3)*A(h)%j(4)
        A_out(h)%j(2) = + K(2)*A(h)%j(4) + K(4)*A(h)%j(3)
        A_out(h)%j(3) = + K(2)*A(h)%j(1) - K(3)*A(h)%j(2)
        A_out(h)%j(4) = + K(1)*A(h)%j(2) - K(4)*A(h)%j(1)
      else
        A_out(h)%j(1) = + K(1)*A(h)%j(3) + K(3)*A(h)%j(4) + M*A(h)%j(1)
        A_out(h)%j(2) = + K(2)*A(h)%j(4) + K(4)*A(h)%j(3) + M*A(h)%j(2)
        A_out(h)%j(3) = + K(2)*A(h)%j(1) - K(3)*A(h)%j(2) + M*A(h)%j(3)
        A_out(h)%j(4) = + K(1)*A(h)%j(2) - K(4)*A(h)%j(1) + M*A(h)%j(4)
      end if
      A_out(h)%h = B11

    end select

  end do

  if (ntry == 1) then
    A_out(:)%n_part  = A(1)%n_part
    A_out(:)%t  = A(1)%t
    A_out(:)%hf = A(:)%hf
    call helbookkeeping_prop(ntry, A, A_out, n)
  end if

end subroutine prop_A_Q


! *****************************************************************************
subroutine prop_W_W(ntry, V, mom, M, mW, V_out, n)
! Unitary-gauge propagator for massive vector bosons
! -----------------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! V(:)          = incoming vector boson wavefunctions (light-cone
!                 representation),
!                 one for each helicity configuration
! K(4)          = incoming vector boson momentum (light-cone rep)
! M             = complex mass
! mW            = 0/1 for massless/massive propagator; not used because
!                 this subroutine is only defined for massive vector bosons.
! V_out(:)      = outgoing dressed vector boson wavefunctions (light-cone
!                 representation)
! V_out(h)%j(i) = V(h)%j(i) - (V(h)%j.K) * K(i)
! *****************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom
  integer(intkind2), intent(inout) :: n
  type(wfun),        intent(in)    :: V(1:n)
  complex(REALKIND), intent(in)    :: M
  integer(intkind1), intent(in)    :: mW
  type(wfun),        intent(out)   :: V_out(1:n)
  complex(REALKIND)                :: M2,K(4)
  integer                          :: h

  K = get_LC_4(mom)
  M2 = M * M
  do h = 1, n
    V_out(h)%j = V(h)%j - (cont_PP(V(h)%j, K) / M2) * K
  end do

  if (ntry == 1) then
    V_out(:)%n_part  = V(1)%n_part
    V_out(:)%t  = V(1)%t
    V_out(:)%hf = V(:)%hf
    call helbookkeeping_prop(ntry, V, V_out, n)
  end if

end subroutine prop_W_W

end module ol_hel_propagators_/**/REALKIND
