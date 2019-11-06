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


module ol_loop_momentum_/**/REALKIND
  implicit none
  contains
  subroutine loop_mom_tens(qloop, qarray)
    ! qloop  = real Lorentz 4-momentum, contravariant
    ! qarray = direct product of qloop in contravariant light cone representation
    use KIND_TYPES, only: REALKIND
    use ol_kinematics_/**/REALKIND, only: Std2LC_cmplx
    use ol_tensor_bookkeeping, only: qproductrules
    implicit none
    complex(REALKIND), intent(in)  :: qloop(0:3)
    complex(REALKIND), intent(out) :: qarray(:)
    complex(REALKIND) :: qloop_lc(4)
    integer :: l
    call Std2LC_cmplx(qloop, qloop_lc)
    qarray(1) = 1
    do l = 2, size(qarray)
      qarray(l) = qloop_lc(qproductrules(1,l)) * qarray(qproductrules(2,l))
    end do
  end subroutine loop_mom_tens
end module ol_loop_momentum_/**/REALKIND
