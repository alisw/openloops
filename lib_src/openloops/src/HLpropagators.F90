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


module ol_h_prop_interface_/**/REALKIND
  use ol_h_vert_interface_/**/REALKIND, only: valid_hol
  use ol_parameters_decl_/**/REALKIND, only: hp_switch
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: req_qp_cmp
#endif
  implicit none
  contains

!***********************************************************************
subroutine Hloop_A_Q(ntry, G_Q, momid, M, Gout_Q, n)
!-----------------------------------------------------------------------
! OpenLoops step: anti-fermion -> fermion propagator
! ----------------------------------------------------------------------
! ntry      = 1 (2) for 1st (subsequent) PS points
! G_Q(:)    = input open-loop
! K         = momentum flowing in the loop propagator
! M         = mass of the propagating particle
! Gout_Q(:) = output open loop with n helicity states
!                 corresp. to "unattached" external legs
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2, QREALKIND
  use ol_data_types_/**/REALKIND, only: hol
  use ol_prop_interface_/**/REALKIND, only: loop_A_Q
  use hel_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_5,get_mass
#ifdef PRECISION_dp
  use ol_prop_interface_/**/QREALKIND, only: loop_A_Q_qp => loop_A_Q
  use ol_kinematics_/**/QREALKIND, only: get_LC_5_qp=>get_LC_5, &
                                         get_mass_qp=>get_mass
#endif
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: momid
  integer(intkind2), intent(inout) :: n
  type(hol),         intent(in)    :: G_Q
  type(hol),         intent(inout) :: Gout_Q
  integer,           intent(in)    :: M
  integer :: h

  if (ntry == 1) call helbookkeeping_prop(ntry, G_Q, Gout_Q, n)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j=0
  do h = 1, n  ! recursion step
    call loop_A_Q(G_Q%j(:,:,:,h),get_LC_5(momid),get_mass(M),Gout_Q%j(:,:,:,h))
  end do

#ifdef PRECISION_dp
  if (req_qp_cmp(G_Q)) then
    Gout_Q%j_qp = 0
    do h = 1, n  ! recursion step
      call loop_A_Q_qp(G_Q%j_qp(:,:,:,h),get_LC_5_qp(momid),get_mass_qp(M), &
                       Gout_Q%j_qp(:,:,:,h))
    end do
  end if
#endif


end subroutine Hloop_A_Q

!***********************************************************************
subroutine Hloop_Q_A(ntry, G_Q, momid, M, Gout_Q, n)
!-----------------------------------------------------------------------
! OpenLoops step: fermion -> anti-fermion propagator
! ----------------------------------------------------------------------
! ntry      = 1 (2) for 1st (subsequent) PS points
! G_Q(:)    = input open-loop
! K         = momentum flowing in the loop propagator
! M         = mass of the propagating particle
! Gout_Q(:) = output open loop with n helicity states
!             corresp. to "unattached" external legs
!**********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2, QREALKIND
  use ol_data_types_/**/REALKIND, only: hol
  use ol_prop_interface_/**/REALKIND, only: loop_Q_A
  use hel_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_5,get_mass
#ifdef PRECISION_dp
  use ol_prop_interface_/**/QREALKIND, only: loop_Q_A_qp => loop_Q_A
  use ol_kinematics_/**/QREALKIND, only: get_LC_5_qp=>get_LC_5, &
                                         get_mass_qp=>get_mass
#endif
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: momid
  integer(intkind2), intent(inout) :: n
  type(hol),         intent(in)    :: G_Q
  type(hol),         intent(inout) :: Gout_Q
  integer,           intent(in)    :: M
  integer :: h

  if (ntry == 1) call helbookkeeping_prop(ntry, G_Q, Gout_Q, n)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j=0
  do h = 1, n  ! recursion step
    call loop_Q_A(G_Q%j(:,:,:,h),get_LC_5(momid),get_mass(M),Gout_Q%j(:,:,:,h))
  end do

#ifdef PRECISION_dp
  if (req_qp_cmp(G_Q)) then
    Gout_Q%j_qp = 0
    do h = 1, n  ! recursion step
      call loop_Q_A_qp(G_Q%j_qp(:,:,:,h),get_LC_5_qp(momid),get_mass_qp(M), &
                       Gout_Q%j_qp(:,:,:,h))
    end do
  end if
#endif

end subroutine Hloop_Q_A

end module ol_h_prop_interface_/**/REALKIND
