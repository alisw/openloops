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


module ol_Std2LC_converter_/**/REALKIND
  implicit none
  contains

  subroutine lorentz2lc_tensor(rank, lor, lc)
    ! convert tensor integrals
    ! lor = Lorentz contravariant -> lc = light-cone contravariant
    ! see init_tensorbookkeeping() for details
    use KIND_TYPES, only: REALKIND
    use ol_parameters_decl_/**/DREALKIND, only: CI
    use ol_tensor_bookkeeping, only: rank_to_size, l2lc
    implicit none
    integer, intent(in) :: rank
    complex(REALKIND), intent(in) :: lor(:)
    complex(REALKIND), intent(out) :: lc(:)
    integer :: pos, r, i, j
    complex(REALKIND) :: tmp
    pos = 1
    lc(1) = lor(1)
    do r = 1, rank
      do i = 1, rank_to_size(r)
        pos = pos + 1
        tmp = 0
        do j = 1, size(l2lc(r)%arr(i)%c,2)
          tmp = tmp + l2lc(r)%arr(i)%c(2,j) * lor(l2lc(r)%arr(i)%c(1,j))
        end do
        tmp = CI * tmp
        do j = 1, size(l2lc(r)%arr(i)%r,2)
          tmp = tmp + l2lc(r)%arr(i)%r(2,j) * lor(l2lc(r)%arr(i)%r(1,j))
        end do
        lc(pos) = tmp
      end do
    end do
  end subroutine lorentz2lc_tensor

end module ol_Std2LC_converter_/**/REALKIND
