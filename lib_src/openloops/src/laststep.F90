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


module ol_last_step_/**/REALKIND
  implicit none
  contains

! **********************************************************************
subroutine check_last_AQ_V(switch, G_A, J_Q, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare AQ -> V gluon-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_AQ_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  complex(REALKIND) :: Gout_V(4,size(G_A,2),4)

  if (switch == 0) then
    call loop_AQ_V(G_A, J_Q, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_AQ_V(G_A, J_Q, Gtensor)

  else if (switch == 2) then
    call loop_AQ_V(G_A, J_Q, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_AQ_V


! **********************************************************************
subroutine check_last_QA_V(switch, G_Q, J_A, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare QA -> V gluon-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_QA_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  complex(REALKIND) :: Gout_V(4,size(G_Q,2),4)

  if (switch == 0) then
    call loop_QA_V(G_Q, J_A, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_QA_V(G_Q, J_A, Gtensor)

  else if (switch == 2) then
    call loop_QA_V(G_Q, J_A, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_QA_V


! **********************************************************************
subroutine check_last_AQ_Z(switch, G_A, J_Q, Gtensor, g_RL)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare AQ -> Z Z-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_AQ_Z
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4), g_RL(2)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  complex(REALKIND) :: Gout_Z(4,size(G_A,2),4)

  if (switch == 0) then
    call loop_AQ_Z(G_A, J_Q, Gout_Z, g_RL)
    call loop_cont_VV(Gout_Z, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_AQ_Z(g_RL, G_A, J_Q, Gtensor)

  else if (switch == 2) then
    call loop_AQ_Z(G_A, J_Q, Gout_Z, g_RL)
    call loop_trace(Gout_Z, Gtensor)

  end if
end subroutine check_last_AQ_Z


! **********************************************************************
subroutine check_last_QA_Z(switch, G_Q, J_A, Gtensor, g_RL)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare QA -> Z Z-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_QA_Z
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4), g_RL(2)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  complex(REALKIND) :: Gout_Z(4,size(G_Q,2),4)

  if (switch == 0) then
    call loop_QA_Z(G_Q, J_A, Gout_Z, g_RL)
    call loop_cont_VV(Gout_Z, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_QA_Z(g_RL, G_Q, J_A, Gtensor)

  else if (switch == 2) then
    call loop_QA_Z(G_Q, J_A, Gout_Z, g_RL)
    call loop_trace(Gout_Z, Gtensor)

  end if
end subroutine check_last_QA_Z


! **********************************************************************
subroutine check_last_AQ_W(switch, G_A, J_Q, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare AQ -> W W-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_AQ_W
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  complex(REALKIND) :: Gout_W(4,size(G_A,2),4)

  if (switch == 0) then
    call loop_AQ_W(G_A, J_Q, Gout_W)
    call loop_cont_VV(Gout_W, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_AQ_W(G_A, J_Q, Gtensor)

  else if (switch == 2) then
    call loop_AQ_W(G_A, J_Q, Gout_W)
    call loop_trace(Gout_W, Gtensor)

  end if
end subroutine check_last_AQ_W


! **********************************************************************
subroutine check_last_QA_W(switch, G_Q, J_A, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare QA -> W gluon-like interaction
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_QA_W
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  complex(REALKIND) :: Gout_W(4,size(G_Q,2),4)

  if (switch == 0) then
    call loop_QA_W(G_Q, J_A, Gout_W)
    call loop_cont_VV(Gout_W, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_QA_W(G_Q, J_A, Gtensor)

  else if (switch == 2) then
    call loop_QA_W(G_Q, J_A, Gout_W)
    call loop_trace(Gout_W, Gtensor)

  end if
end subroutine check_last_QA_W


! **********************************************************************
subroutine check_last_A_Q(switch, G_A, K, M, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! dressing anti-quark current with propagator
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_prop_interface_/**/REALKIND, only: loop_A_Q
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_QA
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), K(5), M
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_A(:,:,:)

  if (switch == 0) then
    allocate(Gout_A(4,HR(4,size(G_A,2)),4))
    call loop_A_Q(G_A, K, M, Gout_A)
    call loop_cont_QA(Gout_A, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_A)

  else if (switch == 1) then
    call last_A_Q(G_A, K, M, Gtensor)

  else if (switch == 2) then
    allocate(Gout_A(4,HR(4,size(G_A,2)),4))
    call loop_A_Q(G_A, K, M, Gout_A)
    call loop_trace(Gout_A, Gtensor)
    deallocate(Gout_A)

  end if
end subroutine check_last_A_Q


! **********************************************************************
subroutine check_last_Q_A(switch, G_Q, K, M, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! dressing quark current with propagator
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_prop_interface_/**/REALKIND, only: loop_Q_A
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_QA
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), K(5), M
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_Q(:,:,:)

  if (switch == 0) then
    allocate(Gout_Q(4,HR(4,size(G_Q,2)),4))
    call loop_Q_A(G_Q, K, M, Gout_Q)
    call loop_cont_QA(Gout_Q, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_Q)

  else if (switch == 1) then
    call last_Q_A(G_Q, K, M, Gtensor)

  else if (switch == 2) then
    allocate(Gout_Q(4,HR(4,size(G_Q,2)),4))
    call loop_Q_A(G_Q, K, M, Gout_Q)
    call loop_trace(Gout_Q, Gtensor)
    deallocate(Gout_Q)

  end if
end subroutine check_last_Q_A


!************************************************************
subroutine check_last_UV_W(switch, Gin_V, Ploop, J_V, Ptree, Gtensor)
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare VV -> V vertex
!************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_UV_W
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), Ploop(4), J_V(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)

  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(Gin_V,2)),4))
    call loop_UV_W(Gin_V, Ploop, J_V, Ptree, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)

  else if (switch == 1) then
    call last_VV_V(Gin_V, Ploop, J_V, Ptree, Gtensor)

  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(Gin_V,2)),4))
    call loop_UV_W(Gin_V, Ploop, J_V, Ptree, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)

  end if
end subroutine check_last_UV_W


! *******************************************************************
subroutine check_last_UW_V(switch, Gin_V, Ploop, J_V, Ptree, Gtensor)
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), Ploop(4), J_V(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  call check_last_UV_W(switch, Gin_V, Ploop, J_V, Ptree, Gtensor)
  Gtensor = -Gtensor
end subroutine check_last_UW_V


!*******************************************************
subroutine check_last_VE_V(switch, Gin, J1, J2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare 4-gluon vertex when the sigma particle enters the
! loop as a tree wavefunction
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_VE_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout(4,size(Gin,2),4)

  if (switch == 0) then
    call loop_VE_V(Gin, J1, J2, Gout)
    call loop_cont_VV(Gout, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_VE_V(Gin, J1, J2, Gtensor)

  else if (switch == 2) then
    call loop_VE_V(Gin, J1, J2, Gout)
    call loop_trace(Gout, Gtensor)

  end if
end subroutine check_last_VE_V


!*******************************************************
subroutine check_last_GGG_G_23(switch, Gin, J1, J2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare 4-gluon vertex when the sigma particle enters the
! loop as a tree wavefunction
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GGG_G_23
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout(4,size(Gin,2),4)

  if (switch == 0) then
    call loop_GGG_G_23(Gin, J1, J2, Gout)
    call loop_cont_VV(Gout, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_GGG_G_23(Gin, J1, J2, Gtensor)

  else if (switch == 2) then
    call loop_GGG_G_23(Gin, J1, J2, Gout)
    call loop_trace(Gout, Gtensor)

  end if
end subroutine check_last_GGG_G_23


!*******************************************************
subroutine check_last_EV_V(switch, Gin, J1, J2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare 4-gluon vertex when the sigma particle is in the loop
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_EV_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout(4,size(Gin,2),4)

  if (switch == 0) then
    call loop_EV_V(Gin, J1, J2, Gout)
    call loop_cont_VV(Gout, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_EV_V(Gin, J1, J2, Gtensor)

  else if (switch == 2) then
    call loop_EV_V(Gin, J1, J2, Gout)
    call loop_trace(Gout, Gtensor)

  end if
end subroutine check_last_EV_V


!*******************************************************
subroutine check_last_GGG_G_12(switch, Gin, J1, J2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare 4-gluon vertex when the sigma particle is in the loop
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GGG_G_12
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout(4,size(Gin,2),4)

  if (switch == 0) then
    call loop_GGG_G_12(Gin, J1, J2, Gout)
    call loop_cont_VV(Gout, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_GGG_G_12(Gin, J1, J2, Gtensor)

  else if (switch == 2) then
    call loop_GGG_G_12(Gin, J1, J2, Gout)
    call loop_trace(Gout, Gtensor)

  end if
end subroutine check_last_GGG_G_12


!*******************************************************
subroutine check_last_CV_D(switch, Gin, Ploop, J_V, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare ghost-gluon -> ghost vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_V(4), Ploop(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_CV_D(Gin, Ploop, J_V, Ptree, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_CV_D(Gin, Ploop, J_V, Ptree, Gtensor)

  else if (switch == 2) then
    call last_CV_D(Gin, Ploop, J_V, Ptree, Gtensor)

  end if
end subroutine check_last_CV_D


!*******************************************************
subroutine check_last_DV_C(switch, Gin, Ploop, J_V, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare antighost-gluon -> antighost vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_V(4), Ploop(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_DV_C(Gin, Ploop, J_V, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_DV_C(Gin, Ploop, J_V, Gtensor)

  else if (switch == 2) then
    call last_DV_C(Gin, Ploop, J_V, Gtensor)

  end if
end subroutine check_last_DV_C


!*******************************************************
subroutine check_last_QA_S(switch, G_Q, J_A, Gtensor, g_RL)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare quark-antiquark -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4), g_RL(2)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_QA_S(g_RL, G_Q, J_A, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_QA_S(g_RL, G_Q, J_A, Gtensor)

  else if (switch == 2) then
    call last_QA_S(g_RL, G_Q, J_A, Gtensor)

  end if
end subroutine check_last_QA_S


!*******************************************************
subroutine check_last_AQ_S(switch, G_A, J_Q, Gtensor, g_RL)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare antiquark-quark -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4), g_RL(2)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_AQ_S(g_RL, G_A, J_Q, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_AQ_S(g_RL, G_A, J_Q, Gtensor)

  else if (switch == 2) then
    call last_AQ_S(g_RL, G_A, J_Q, Gtensor)

  end if
end subroutine check_last_AQ_S


!*******************************************************
subroutine check_last_VV_S(switch, G_V, J_V, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-vector -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_VV_S(G_V, J_V, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_VV_S(G_V, J_V, Gtensor)

  else if (switch == 2) then
    call last_VV_S(G_V, J_V, Gtensor)

  end if
end subroutine check_last_VV_S


!*******************************************************
subroutine check_last_VS_V(switch, G_V, J_S, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-scalar -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_VS_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND)              :: Gout_V(4,size(G_V,2),4)

  if (switch == 0) then
    call loop_VS_V(G_V, J_S, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_VS_V(G_V, J_S, Gtensor)

  else if (switch == 2) then
    call loop_VS_V(G_V, J_S, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_VS_V


!*******************************************************
subroutine check_last_SV_V(switch, G_S, J_V, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-vector -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_SV_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_V(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND)              :: Gout_V(4,size(G_S,2),4)

  if (switch == 0) then
    call loop_SV_V(G_S, J_V, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_SV_V(G_S, J_V, Gtensor)

  else if (switch == 2) then
    call loop_SV_V(G_S, J_V, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_SV_V


!*******************************************************
subroutine check_last_SV_T(switch, G_S, Ploop, J_V, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-vector -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_V(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_SV_T(G_S, Ploop, J_V, Ptree, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_SV_T(G_S, Ploop, J_V, Ptree, Gtensor)

  else if (switch == 2) then
    call last_SV_T(G_S, Ploop, J_V, Ptree, Gtensor)

  end if
end subroutine check_last_SV_T


!*******************************************************
subroutine check_last_TV_S(switch, G_S, Ploop, J_V, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-vector -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_V(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_TV_S(G_S, Ploop, J_V, Ptree, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_TV_S(G_S, Ploop, J_V, Ptree, Gtensor)

  else if (switch == 2) then
    call last_TV_S(G_S, Ploop, J_V, Ptree, Gtensor)

  end if
end subroutine check_last_TV_S


!*******************************************************
subroutine check_last_VS_T(switch, G_V, Ploop, J_S, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-scalar -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), Ploop(4), J_S(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_VS_T(G_V, Ploop, J_S, Ptree, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_VS_T(G_V, Ploop, J_S, Ptree, Gtensor)

  else if (switch == 2) then
    call last_VS_T(G_V, Ploop, J_S, Ptree, Gtensor)

  end if
end subroutine check_last_VS_T


!*******************************************************
subroutine check_last_VT_S(switch, G_V, Ploop, J_S, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-scalar -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), Ploop(4), J_S(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_VT_S(G_V, Ploop, J_S, Ptree, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_VT_S(G_V, Ploop, J_S, Ptree, Gtensor)

  else if (switch == 2) then
    call last_VT_S(G_V, Ploop, J_S, Ptree, Gtensor)

  end if
end subroutine check_last_VT_S


!*******************************************************
subroutine check_last_ST_V(switch, G_S, Ploop, J_S, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-scalar -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_ST_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_S(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)

  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(G_S,2)),4))
    call loop_ST_V(G_S, Ploop, J_S, Ptree, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)

  else if (switch == 1) then
    call last_ST_V(G_S, Ploop, J_S, Ptree, Gtensor)

  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(G_S,2)),4))
    call loop_ST_V(G_S, Ploop, J_S, Ptree, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)

  end if
end subroutine check_last_ST_V


!*******************************************************
subroutine check_last_TS_V(switch, G_S, Ploop, J_S, Ptree, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-scalar -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_TS_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_S(4), Ptree(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)

  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(G_S,2)),4))
    call loop_TS_V(G_S, Ploop, J_S, Ptree, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)

  else if (switch == 1) then
    call last_TS_V(G_S, Ploop, J_S, Ptree, Gtensor)

  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(G_S,2)),4))
    call loop_TS_V(G_S, Ploop, J_S, Ptree, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)

  end if
end subroutine check_last_TS_V


!*******************************************************
subroutine check_last_SS_S(switch, G_S, J_S, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-scalar -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_SS_S(G_S, J_S, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_SS_S(G_S, J_S, Gtensor)

  else if (switch == 2) then
    call last_SS_S(G_S, J_S, Gtensor)

  end if
end subroutine check_last_SS_S


!*******************************************************
subroutine check_last_SSS_S(switch, G_S, J_S1, J_S2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-scalar-scalar -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S1(4), J_S2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_SSS_S(G_S, J_S1, J_S2, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_SSS_S(G_S, J_S1, J_S2, Gtensor)

  else if (switch == 2) then
    call last_SSS_S(G_S, J_S1, J_S2, Gtensor)

  end if
end subroutine check_last_SSS_S


!*******************************************************
subroutine check_last_VVS_S(switch, G_V, J_V, J_S, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-vector-scalar -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V(4), J_S(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_VVS_S(G_V, J_V, J_S, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_VVS_S(G_V, J_V, J_S, Gtensor)

  else if (switch == 2) then
    call last_VVS_S(G_V, J_V, J_S, Gtensor)

  end if
end subroutine check_last_VVS_S


!*******************************************************
subroutine check_last_SSV_V(switch, G_S, J_S, J_V, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-scalar-vector -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_SSV_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S(4), J_V(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND)              :: Gout_V(4,size(G_S,2),4)

  if (switch == 0) then
    call loop_SSV_V(G_S, J_S, J_V, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_SSV_V(G_S, J_S, J_V, Gtensor)

  else if (switch == 2) then
    call loop_SSV_V(G_S, J_S, J_V, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_SSV_V


!*******************************************************
subroutine check_last_VSS_V(switch, G_V, J_S1, J_S2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare vector-scalar-scalar -> vector vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_VSS_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND)              :: Gout_V(4,size(G_V,2),4)

  if (switch == 0) then
    call loop_VSS_V(G_V, J_S1, J_S2, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_VSS_V(G_V, J_S1, J_S2, Gtensor)

  else if (switch == 2) then
    call loop_VSS_V(G_V, J_S1, J_S2, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_VSS_V


!*******************************************************
subroutine check_last_SVV_S(switch, G_S, J_V1, J_V2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare scalar-vector-vector -> scalar vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_V1(4), J_V2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)

  if (switch == 0) then
    call last_SVV_S(G_S, J_V1, J_V2, Gtensor)
    Gtensor = exloop(1,2)*exloop(1,1)*Gtensor

  else if (switch == 1) then
    call last_SVV_S(G_S, J_V1, J_V2, Gtensor)

  else if (switch == 2) then
    call last_SVV_S(G_S, J_V1, J_V2, Gtensor)

  end if
end subroutine check_last_SVV_S


!*******************************************************
subroutine check_last_WWV_V(switch, G_V, JV1, JV2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare V V V -> V vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_WWV_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), JV1(4), JV2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout_V(4,size(G_V,2),4)

  if (switch == 0) then
    call loop_WWV_V(G_V, JV1, JV2, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_WWV_V(G_V, JV1, JV2, Gtensor)

  else if (switch == 2) then
    call loop_WWV_V(G_V, JV1, JV2, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_WWV_V


!*******************************************************
subroutine check_last_VWW_V(switch, G_V, JV1, JV2, Gtensor)
!-------------------------------------------------------
! 0 = tree; 1 = loop in one step; 2 = loop in two steps
! bare V V V -> V vertex
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_VWW_V
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), JV1(4), JV2(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: Gout_V(4,size(G_V,2),4)

  if (switch == 0) then
    call loop_VWW_V(G_V, JV1, JV2, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)

  else if (switch == 1) then
    call last_VWW_V(G_V, JV1, JV2, Gtensor)

  else if (switch == 2) then
    call loop_VWW_V(G_V, JV1, JV2, Gout_V)
    call loop_trace(Gout_V, Gtensor)

  end if
end subroutine check_last_VWW_V


! ======================================================================
! Last vertices  for HEFT.
! ======================================================================

!*******************************************************
subroutine check_last_GH_G(switch, G_V, plin, J_S, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GH_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GH_G(G_V, plin, J_S, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GH_G(G_V, plin, J_S, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GH_G(G_V, plin, J_S, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GH_G


!*******************************************************
subroutine check_last_GHG_G(switch, G_V, plin, J_S, J_V3, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHG_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S(4), J_V3(4), p3(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHG_G(G_V, plin, J_S, J_V3, p3, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHG_G(G_V, plin, J_S, J_V3, p3, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHG_G(G_V, plin, J_S, J_V3, p3, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHG_G


!*******************************************************
subroutine check_last_GHGG_G_12(switch, G_V, J_S, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHGG_G_12
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHGG_G_12(G_V, J_S, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHGG_G_12(G_V, J_S, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHGG_G_12(G_V, J_S, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHGG_G_12


!*******************************************************
subroutine check_last_GHGG_G_23(switch, G_V, J_S, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHGG_G_23
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHGG_G_23(G_V, J_S, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHGG_G_23(G_V, J_S, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHGG_G_23(G_V, J_S, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHGG_G_23

! ======================================================================
! Last vertices  for HHEFT.
! ======================================================================

!*******************************************************
subroutine check_last_GHH_G(switch, G_V, plin, J_S1, J_S2, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHH_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S1(4), J_S2(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GHH_G(G_V, plin, J_S1, J_S2, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHH_G(G_V, plin, J_S1, J_S2, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GHH_G(G_V, plin, J_S1, J_S2, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHH_G


!*******************************************************
subroutine check_last_GHHG_G(switch, G_V, plin, J_S1, J_S2, J_V3, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHG_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S1(4), J_S2(4), J_V3(4), p3(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHHG_G(G_V, plin, J_S1, J_S2, J_V3, p3, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHG_G(G_V, plin, J_S1, J_S2, J_V3, p3, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHHG_G(G_V, plin, J_S1, J_S2, J_V3, p3, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHG_G


!*******************************************************
subroutine check_last_GHHGG_G_12(switch, G_V, J_S1, J_S2, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHGG_G_12
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHGG_G_12(G_V, J_S1, J_S2, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHGG_G_12(G_V, J_S1, J_S2, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHGG_G_12(G_V, J_S1, J_S2, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHGG_G_12


!*******************************************************
subroutine check_last_GHHGG_G_23(switch, G_V, J_S1, J_S2, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHGG_G_23
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHGG_G_23(G_V, J_S1, J_S2, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHGG_G_23(G_V, J_S1, J_S2, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHGG_G_23(G_V, J_S1, J_S2, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHGG_G_23


! ======================================================================
! Last vertices  for HHHEFT.
! ======================================================================

!*******************************************************
subroutine check_last_GHHH_G(switch, G_V, plin, J_S1, J_S2, J_S3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHH_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S1(4), J_S2(4), J_S3(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GHHH_G(G_V, plin, J_S1, J_S2, J_S3, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHH_G(G_V, plin, J_S1, J_S2, J_S3, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,HR(4,size(G_V,2))),4))
    call loop_GHHH_G(G_V, plin, J_S1, J_S2, J_S3, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHH_G


!*******************************************************
subroutine check_last_GHHHG_G(switch, G_V, plin, J_S1, J_S2, J_S3, J_V3, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHHG_G
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), plin(4), J_S1(4), J_S2(4), J_S3(4), J_V3(4), p3(4), plout(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHHHG_G(G_V, plin, J_S1, J_S2, J_S3, J_V3, p3, Gout_V, plout)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHHG_G(G_V, plin, J_S1, J_S2, J_S3, J_V3, p3, Gtensor, plout)
  else if (switch == 2) then
    allocate(Gout_V(4,HR(4,size(G_V,2)),4))
    call loop_GHHHG_G(G_V, plin, J_S1, J_S2, J_S3, J_V3, p3, Gout_V, plout)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHHG_G


!*******************************************************
subroutine check_last_GHHHGG_G_12(switch, G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHHGG_G_12
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4), J_S3(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHHGG_G_12(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHHGG_G_12(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHHGG_G_12(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHHGG_G_12


!*******************************************************
subroutine check_last_GHHHGG_G_23(switch, G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: loop_GHHHGG_G_23
  use ol_loop_routines_/**/REALKIND, only: loop_trace, loop_cont_VV
  use ol_pseudotree_/**/REALKIND, only: exloop
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4), J_S3(4), J_V3(4), J_V4(4)
  integer,           intent(in)  :: switch
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND), allocatable :: Gout_V(:,:,:)
  if (switch == 0) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHHGG_G_23(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gout_V)
    call loop_cont_VV(Gout_V, exloop(:,2), exloop(:,1), Gtensor)
    deallocate(Gout_V)
  else if (switch == 1) then
    call last_GHHHGG_G_23(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
  else if (switch == 2) then
    allocate(Gout_V(4,size(G_V,2),4))
    call loop_GHHHGG_G_23(G_V, J_S1, J_S2, J_S3, J_V3, J_V4, Gout_V)
    call loop_trace(Gout_V, Gtensor)
    deallocate(Gout_V)
  end if
end subroutine check_last_GHHHGG_G_23


! **********************************************************************
!                          LAST STEP ROUTINES
! **********************************************************************

! **********************************************************************
subroutine last_AQ_V(G_A, J_Q, Gtensor)
! bare AQ -> V gluon-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4)
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  integer :: l
  do l = 1, size(G_A,2)
    Gtensor(l) = (( - G_A(1,l,1)*J_Q(3) - G_A(4,l,1)*J_Q(2))  & ! alpha = beta = 1
               +  ( - G_A(2,l,2)*J_Q(4) - G_A(3,l,2)*J_Q(1))  & ! alpha = beta = 2
               +  ( - G_A(1,l,3)*J_Q(4) + G_A(3,l,3)*J_Q(2))  & ! alpha = beta = 3
               +  ( - G_A(2,l,4)*J_Q(3) + G_A(4,l,4)*J_Q(1)))   ! alpha = beta = 4
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_AQ_V


! **********************************************************************
subroutine last_QA_V(G_Q, J_A, Gtensor)
! bare QA -> V gluon-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4)
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  integer :: l
  do l = 1, size(G_Q,2)
    Gtensor(l) = (( - G_Q(3,l,1)*J_A(1) - G_Q(2,l,1)*J_A(4))  & ! alpha = beta = 1
               +  ( - G_Q(4,l,2)*J_A(2) - G_Q(1,l,2)*J_A(3))  & ! alpha = beta = 2
               +  ( - G_Q(4,l,3)*J_A(1) + G_Q(2,l,3)*J_A(3))  & ! alpha = beta = 3
               +  ( - G_Q(3,l,4)*J_A(2) + G_Q(1,l,4)*J_A(4)))   ! alpha = beta = 4
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_QA_V


! **********************************************************************
subroutine last_AQ_Z(g_RL, G_A, J_Q, Gtensor)
! bare AQ -> Z Z-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4), g_RL(2)
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  integer :: l
  do l = 1, size(G_A,2)
    Gtensor(l) = (- g_RL(1)*(G_A(4,l,1)*J_Q(2)  &
                           + G_A(3,l,2)*J_Q(1)  &
                           - G_A(3,l,3)*J_Q(2)  &
                           - G_A(4,l,4)*J_Q(1)) &
                  - g_RL(2)*(G_A(1,l,1)*J_Q(3)  &
                           + G_A(2,l,2)*J_Q(4)  &
                           + G_A(1,l,3)*J_Q(4)  &
                           + G_A(2,l,4)*J_Q(3)))
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_AQ_Z


! **********************************************************************
subroutine last_QA_Z(g_RL, G_Q, J_A, Gtensor)
! bare QA -> Z Z-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4), g_RL(2)
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  integer :: l
  do l = 1, size(G_Q,2)
    Gtensor(l) = (- g_RL(1)*(G_Q(2,l,1)*J_A(4)  &
                           + G_Q(1,l,2)*J_A(3)  &
                           - G_Q(2,l,3)*J_A(3)  &
                           - G_Q(1,l,4)*J_A(4)) &
                  - g_RL(2)*(G_Q(3,l,1)*J_A(1)  &
                           + G_Q(4,l,2)*J_A(2)  &
                           + G_Q(4,l,3)*J_A(1)  &
                           + G_Q(3,l,4)*J_A(2)))
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_QA_Z


! **********************************************************************
subroutine last_AQ_W(G_A, J_Q, Gtensor)
! bare AQ -> W W-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4)
  complex(REALKIND), intent(out) :: Gtensor(size(G_A,2))
  integer :: l
  do l = 1, size(G_A,2)
    Gtensor(l) = (- G_A(1,l,1)*J_Q(3)  &
                  - G_A(2,l,2)*J_Q(4)  &
                  - G_A(1,l,3)*J_Q(4)  &
                  - G_A(2,l,4)*J_Q(3))
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_AQ_W


! **********************************************************************
subroutine last_QA_W(G_Q, J_A, Gtensor)
! bare QA -> W W-like interaction
! ----------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4)
  complex(REALKIND), intent(out) :: Gtensor(size(G_Q,2))
  integer :: l
  do l = 1, size(G_Q,2)
    Gtensor(l) = (- G_Q(3,l,1)*J_A(1)  &
                  - G_Q(4,l,2)*J_A(2)  &
                  - G_Q(4,l,3)*J_A(1)  &
                  - G_Q(3,l,4)*J_A(2))
  end do
  Gtensor = Gtensor + Gtensor
end subroutine last_QA_W


! **********************************************************************
subroutine last_A_Q(G_A,K,M,Gtensor)
! dressing anti-quark current with propagator
! ----------------------------------------------------------------------
! dress fermion, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), K(5), M
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: length_in, l

  length_in = size(G_A,2)
  Gtensor = 0

  do l = 1, length_in
    ! dressing proportional to loop momentum, raise the rank
    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + G_A(3,l,1) + G_A(2,l,4)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) + G_A(4,l,2) + G_A(1,l,3)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) + G_A(4,l,1) - G_A(2,l,3)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) + G_A(3,l,2) - G_A(1,l,4)
    ! dressing proportional to tree level momentum, same rank
    Gtensor(l) = Gtensor(l) + K(1)*G_A(3,l,1) + K(3)*G_A(4,l,1) + M*G_A(1,l,1) & ! alpha = beta = 1
                            + K(2)*G_A(4,l,2) + K(4)*G_A(3,l,2) + M*G_A(2,l,2) & ! alpha = beta = 2
                            + K(2)*G_A(1,l,3) - K(3)*G_A(2,l,3) + M*G_A(3,l,3) & ! alpha = beta = 3
                            + K(1)*G_A(2,l,4) - K(4)*G_A(1,l,4) + M*G_A(4,l,4)   ! alpha = beta = 4
  end do
end subroutine last_A_Q


! **********************************************************************
subroutine last_Q_A(G_Q,K,M,Gtensor)
! dressing quark current with propagator
! ----------------------------------------------------------------------
! dress fermion, calculate only diagonal elements and take trace
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), K(5), M
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: length_in, l

  length_in = size(G_Q,2)
  Gtensor = 0

  do l = 1, length_in
    ! dressing proportional to loop momentum, raise the rank
    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) - G_Q(4,l,2) - G_Q(1,l,3)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) - G_Q(3,l,1) - G_Q(2,l,4)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) + G_Q(3,l,2) - G_Q(1,l,4)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) + G_Q(4,l,1) - G_Q(2,l,3)
    ! dressing proportional to tree level momentum, same rank
    Gtensor(l) = Gtensor(l) - K(2)*G_Q(3,l,1) + K(4)*G_Q(4,l,1) + M*G_Q(1,l,1) & ! alpha = beta = 1
                            - K(1)*G_Q(4,l,2) + K(3)*G_Q(3,l,2) + M*G_Q(2,l,2) & ! alpha = beta = 2
                            - K(1)*G_Q(1,l,3) - K(4)*G_Q(2,l,3) + M*G_Q(3,l,3) & ! alpha = beta = 3
                            - K(2)*G_Q(2,l,4) - K(3)*G_Q(1,l,4) + M*G_Q(4,l,4)   ! alpha = beta = 4
  end do
end subroutine last_Q_A


!************************************************************
subroutine last_VV_V(Gin_V, Ploop, J_V, Ptree, Gtensor)
! bare VV -> V vertex
!------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
!************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), Ploop(4), J_V(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: A(4), B(4), C, Jhalf(4), Jtwo(4), Ptmp(4,3)
  integer           :: l, length_in

  length_in = size(Gin_V,2)
  Gtensor = 0
  Ptmp(:,1) = Ploop + 2*Ptree
  Ptmp(:,2) = Ptree + 2*Ploop
  Ptmp(:,3) = Ploop - Ptree
  C = cont_VV(Ptmp(:,2),J_V)
  Jhalf = 0.5_/**/REALKIND * J_V
! covariant components of 2*J_V_mu = 2 * g_(mu,nu)*J_V^nu, factor 2 simplifies in metric tensor
  Jtwo(1) =  J_V(2)
  Jtwo(2) =  J_V(1)
  Jtwo(3) = -J_V(4)
  Jtwo(4) = -J_V(3)

  do l = 1, length_in
    A(1) = cont_VV(Gin_V(:,l,1), J_V)
    A(2) = cont_VV(Gin_V(:,l,2), J_V)
    A(3) = cont_VV(Gin_V(:,l,3), J_V)
    A(4) = cont_VV(Gin_V(:,l,4), J_V)
    B(1) = cont_VV(Gin_V(:,l,1), Ptmp(:,1))
    B(2) = cont_VV(Gin_V(:,l,2), Ptmp(:,1))
    B(3) = cont_VV(Gin_V(:,l,3), Ptmp(:,1))
    B(4) = cont_VV(Gin_V(:,l,4), Ptmp(:,1))

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + A(1) + Gin_V(2,l,1)*Jhalf(1) - Gin_V(1,l,1)*Jtwo(1) & ! alpha = beta = 1
                                               + Gin_V(2,l,2)*Jhalf(2) - Gin_V(2,l,2)*Jtwo(1) & ! alpha = beta = 2
                                               + Gin_V(2,l,3)*Jhalf(3) - Gin_V(3,l,3)*Jtwo(1) & ! alpha = beta = 3
                                               + Gin_V(2,l,4)*Jhalf(4) - Gin_V(4,l,4)*Jtwo(1)   ! alpha = beta = 4

    Gtensor(HR(2,l)) = Gtensor(HR(2,l))        + Gin_V(1,l,1)*Jhalf(1) - Gin_V(1,l,1)*Jtwo(2) & ! alpha = beta = 1
                                        + A(2) + Gin_V(1,l,2)*Jhalf(2) - Gin_V(2,l,2)*Jtwo(2) & ! alpha = beta = 2
                                               + Gin_V(1,l,3)*Jhalf(3) - Gin_V(3,l,3)*Jtwo(2) & ! alpha = beta = 3
                                               + Gin_V(1,l,4)*Jhalf(4) - Gin_V(4,l,4)*Jtwo(2)   ! alpha = beta = 4

    Gtensor(HR(3,l)) = Gtensor(HR(3,l))        - Gin_V(4,l,1)*Jhalf(1) - Gin_V(1,l,1)*Jtwo(3) & ! alpha = beta = 1
                                               - Gin_V(4,l,2)*Jhalf(2) - Gin_V(2,l,2)*Jtwo(3) & ! alpha = beta = 2
                                        + A(3) - Gin_V(4,l,3)*Jhalf(3) - Gin_V(3,l,3)*Jtwo(3) & ! alpha = beta = 3
                                               - Gin_V(4,l,4)*Jhalf(4) - Gin_V(4,l,4)*Jtwo(3)   ! alpha = beta = 4

    Gtensor(HR(4,l)) = Gtensor(HR(4,l))        - Gin_V(3,l,1)*Jhalf(1) - Gin_V(1,l,1)*Jtwo(4) & ! alpha = beta = 1
                                               - Gin_V(3,l,2)*Jhalf(2) - Gin_V(2,l,2)*Jtwo(4) & ! alpha = beta = 2
                                               - Gin_V(3,l,3)*Jhalf(3) - Gin_V(3,l,3)*Jtwo(4) & ! alpha = beta = 3
                                        + A(4) - Gin_V(3,l,4)*Jhalf(4) - Gin_V(4,l,4)*Jtwo(4)   ! alpha = beta = 4

    Gtensor(l) = Gtensor(l) + A(1)*Ptmp(1,3) + B(1)*J_V(1) - C*Gin_V(1,l,1) & ! alpha = beta = 1
                            + A(2)*Ptmp(2,3) + B(2)*J_V(2) - C*Gin_V(2,l,2) & ! alpha = beta = 2
                            + A(3)*Ptmp(3,3) + B(3)*J_V(3) - C*Gin_V(3,l,3) & ! alpha = beta = 3
                            + A(4)*Ptmp(4,3) + B(4)*J_V(4) - C*Gin_V(4,l,4)   ! alpha = beta = 4
  end do
end subroutine last_VV_V


! *******************************************************************
subroutine last_VE_V(Gin, J1, J2, Gtensor)
!--------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle enters the loop
! as a tree wavefunction
!--------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! *******************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J2(1)+Gin(:,l,2)*J2(2)+Gin(:,l,3)*J2(3)+Gin(:,l,4)*J2(4), J1) &
               - cont_VV(Gin(:,l,1)*J1(1)+Gin(:,l,2)*J1(2)+Gin(:,l,3)*J1(3)+Gin(:,l,4)*J1(4), J2)
  end do
end subroutine last_VE_V


! *******************************************************************
subroutine last_GGG_G_23(Gin, J1, J2, Gtensor)
! *******************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J2(1)+Gin(:,l,2)*J2(2)+Gin(:,l,3)*J2(3)+Gin(:,l,4)*J2(4), J1) &
               - cont_VV(Gin(:,l,1)*J1(1)+Gin(:,l,2)*J1(2)+Gin(:,l,3)*J1(3)+Gin(:,l,4)*J1(4), J2)
  end do
end subroutine last_GGG_G_23


! *******************************************************************
subroutine last_EV_V(Gin, J1, J2, Gtensor)
!--------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle is in the loop
!--------------------------------------------------------------------
! attach last vertex, calculate only diagonal elements and take trace
! *******************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: J1J2
  integer :: l

  J1J2 = cont_VV(J1,J2)
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J1(1)+Gin(:,l,2)*J1(2)+Gin(:,l,3)*J1(3)+Gin(:,l,4)*J1(4), J2) &
               - J1J2 * (Gin(1,l,1) + Gin(2,l,2) + Gin(3,l,3) + Gin(4,l,4))
  end do
end subroutine last_EV_V


! *******************************************************************
subroutine last_GGG_G_12(Gin, J1, J2, Gtensor)
! *******************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J1(4), J2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: J1J2
  integer :: l

  J1J2 = cont_VV(J1,J2)
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J1(1)+Gin(:,l,2)*J1(2)+Gin(:,l,3)*J1(3)+Gin(:,l,4)*J1(4), J2) &
               - J1J2 * (Gin(1,l,1) + Gin(2,l,2) + Gin(3,l,3) + Gin(4,l,4))
  end do
end subroutine last_GGG_G_12


!****************************************************
subroutine last_CV_D(Gin, Ploop, J_V, Ptree, Gtensor)
!----------------------------------------------------
! bare ghost-gluon -> ghost vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_CV_D
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_V(4), Ploop(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(Gin,2)
  rank_out = size(Gtensor)

  call vert_loop_CV_D(rank_in, rank_out, Gin(1,:,1), Ploop, J_V, Ptree, Gtensor)

end subroutine last_CV_D


!****************************************************
subroutine last_DV_C(Gin, Ploop, J_V, Gtensor)
!----------------------------------------------------
! bare antighost-gluon -> antighost vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_DV_C
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_V(4), Ploop(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(Gin,2)
  rank_out = size(Gtensor)

  call vert_loop_DV_C(rank_in, rank_out, Gin(1,:,1), Ploop, J_V, Gtensor)

end subroutine last_DV_C


!****************************************************
subroutine last_QA_S(g_RL, G_Q, J_A, Gtensor)
!----------------------------------------------------
! bare quark-antiquark -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_QA_S
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), J_A(4), g_RL(2)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_Q,2)
  rank_out = size(Gtensor)

  call vert_loop_QA_S(rank_in, rank_out, g_RL, G_Q(:,:,1), J_A, Gtensor)

end subroutine last_QA_S


!****************************************************
subroutine last_AQ_S(g_RL, G_A, J_Q, Gtensor)
!----------------------------------------------------
! bare antiquark-quark -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_AQ_S
  implicit none
  complex(REALKIND), intent(in)  :: G_A(:,:,:), J_Q(4), g_RL(2)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_A,2)
  rank_out = size(Gtensor)

  call vert_loop_AQ_S(rank_in, rank_out, g_RL, G_A(:,:,1), J_Q, Gtensor)

end subroutine last_AQ_S


!****************************************************
subroutine last_VV_S(G_V, J_V, Gtensor)
!----------------------------------------------------
! bare vector-vector -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_VV_S
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_V,2)
  rank_out = size(Gtensor)

  call vert_loop_VV_S(rank_in, rank_out, G_V(:,:,1), J_V, Gtensor)

end subroutine last_VV_S


!****************************************************
subroutine last_VS_V(G_V, J_S, Gtensor)
!----------------------------------------------------
! bare vector-scalar -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l
  do l = 1, size(G_V,2)
    Gtensor(l) = (G_V(1,l,1) + G_V(2,l,2) + G_V(3,l,3) + G_V(4,l,4))*J_S(1)
  end do
end subroutine last_VS_V


!****************************************************
subroutine last_SV_V(G_S, J_V, Gtensor)
!----------------------------------------------------
! bare scalar-vector -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_V(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l
  do l = 1, size(G_S,2)
    Gtensor(l) = G_S(1,l,1)*J_V(1) + G_S(1,l,2)*J_V(2) + G_S(1,l,3)*J_V(3) + G_S(1,l,4)*J_V(4)
  end do
end subroutine last_SV_V


!****************************************************
subroutine last_SV_T(G_S, Ploop, J_V, Ptree, Gtensor)
!----------------------------------------------------
! bare scalar-vector -> scalar vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_SV_T
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_V(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_S,2)
  rank_out = size(Gtensor)

  call vert_loop_SV_T(rank_in, rank_out, G_S(1,:,1), Ploop, J_V, Ptree, Gtensor)
end subroutine last_SV_T


!****************************************************
subroutine last_TV_S(G_S, Ploop, J_V, Ptree, Gtensor)
!----------------------------------------------------
! bare scalar-vector -> scalar vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_TV_S
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_V(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_S,2)
  rank_out = size(Gtensor)

  call vert_loop_TV_S(rank_in, rank_out, G_S(1,:,1), Ploop, J_V, Ptree, Gtensor)
end subroutine last_TV_S


!****************************************************
subroutine last_VS_T(G_V, Ploop, J_S, Ptree, Gtensor)
!----------------------------------------------------
! bare vector-scalar -> scalar vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_VS_T
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), Ploop(4), J_S(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_V,2)
  rank_out = size(Gtensor)

  call vert_loop_VS_T(rank_in, rank_out, G_V(:,:,1), Ploop, J_S(1), Ptree, Gtensor)
end subroutine last_VS_T


!****************************************************
subroutine last_VT_S(G_V, Ploop, J_S, Ptree, Gtensor)
!----------------------------------------------------
! bare vector-scalar -> scalar vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_VT_S
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), Ploop(4), J_S(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_V,2)
  rank_out = size(Gtensor)

  call vert_loop_VT_S(rank_in, rank_out, G_V(:,:,1), Ploop, J_S(1), Ptree, Gtensor)
end subroutine last_VT_S


!****************************************************
subroutine last_ST_V(G_S, Ploop, J_S, Ptree, Gtensor)
!----------------------------------------------------
! bare scalar-scalar -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_S(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  Gtensor = 0

  do l = 1, size(G_S,2)

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + G_S(1,l,1)*J_S(1)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) + G_S(1,l,2)*J_S(1)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) + G_S(1,l,3)*J_S(1)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) + G_S(1,l,4)*J_S(1)

    Gtensor(l) = Gtensor(l) + G_S(1,l,1)*J_S(1)*(Ploop(1)-Ptree(1)) + G_S(1,l,2)*J_S(1)*(Ploop(2)-Ptree(2)) &
                            + G_S(1,l,3)*J_S(1)*(Ploop(3)-Ptree(3)) + G_S(1,l,4)*J_S(1)*(Ploop(4)-Ptree(4))
  end do
end subroutine last_ST_V


!****************************************************
subroutine last_TS_V(G_S, Ploop, J_S, Ptree, Gtensor)
!----------------------------------------------------
! bare scalar-scalar -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), Ploop(4), J_S(4), Ptree(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  Gtensor = 0

  do l = 1, size(G_S,2)

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) - G_S(1,l,1)*J_S(1)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) - G_S(1,l,2)*J_S(1)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) - G_S(1,l,3)*J_S(1)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) - G_S(1,l,4)*J_S(1)

    Gtensor(l) = Gtensor(l) + G_S(1,l,1)*J_S(1)*(-Ploop(1)+Ptree(1)) + G_S(1,l,2)*J_S(1)*(-Ploop(2)+Ptree(2)) &
                            + G_S(1,l,3)*J_S(1)*(-Ploop(3)+Ptree(3)) + G_S(1,l,4)*J_S(1)*(-Ploop(4)+Ptree(4))
  end do
end subroutine last_TS_V


!****************************************************
subroutine last_SS_S(G_S, J_S, Gtensor)
!----------------------------------------------------
! bare scalar-scalar -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_SS_S
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_S,2)
  rank_out = size(Gtensor)

  call vert_loop_SS_S(rank_in, rank_out, G_S(1,:,1), J_S(1), Gtensor)

end subroutine last_SS_S


!****************************************************
subroutine last_SSS_S(G_S, J_S1, J_S2, Gtensor)
!----------------------------------------------------
! bare scalar-scalar-scalar -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_SSS_S
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S1(4), J_S2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_S,2)
  rank_out = size(Gtensor)

  call vert_loop_SSS_S(rank_in, rank_out, G_S(1,:,1), J_S1(1), J_S2(1), Gtensor)

end subroutine last_SSS_S


!****************************************************
subroutine last_VVS_S(G_V, J_V, J_S, Gtensor)
!----------------------------------------------------
! bare vector-vector-scalar -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_VVS_S
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V(4), J_S(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_V,2)
  rank_out = size(Gtensor)

  call vert_loop_VVS_S(rank_in, rank_out, G_V(:,:,1), J_V, J_S(1), Gtensor)

end subroutine last_VVS_S


!****************************************************
subroutine last_SSV_V(G_S, J_S, J_V, Gtensor)
!----------------------------------------------------
! bare scalar-scalar-vector -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_S(4), J_V(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  do l = 1, size(G_S,2)
    Gtensor(l) = J_S(1) * (G_S(1,l,1) * J_V(1) + G_S(1,l,2) * J_V(2) + G_S(1,l,3) * J_V(3) + G_S(1,l,4) * J_V(4))
  end do
end subroutine last_SSV_V


!****************************************************
subroutine last_VSS_V(G_V, J_S1, J_S2, Gtensor)
!----------------------------------------------------
! bare vector-scalar-scalar -> vector vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_S1(4), J_S2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: l

  do l = 1, size(G_V,2)
    Gtensor(l) = J_S1(1) * J_S2(1) * (G_V(1,l,1) + G_V(2,l,2) + G_V(3,l,3) + G_V(4,l,4))
  end do
end subroutine last_VSS_V


!****************************************************
subroutine last_SVV_S(G_S, J_V1, J_V2, Gtensor)
!----------------------------------------------------
! bare scalar-vector-vector -> scalar vertex
!------------------------------------------------------------
! attach last vertex, same as intermediate vertex (scalar coeff)
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_vert_interface_/**/REALKIND, only: vert_loop_SVV_S
  implicit none
  complex(REALKIND), intent(in)  :: G_S(:,:,:), J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer :: rank_in, rank_out
  rank_in  = size(G_S,2)
  rank_out = size(Gtensor)

  call vert_loop_SVV_S(rank_in, rank_out, G_S(1,:,1), J_V1, J_V2, Gtensor)

end subroutine last_SVV_S


!****************************************************
subroutine last_WWV_V(G_V, J_V1, J_V2, Gtensor)
!----------------------------------------------------
! bare V V V -> V vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer           :: l
  complex(REALKIND) :: J1J2

  J1J2 = cont_VV(J_V1, J_V2)
  do l = 1, size(G_V,2)
    Gtensor(l) = 2*cont_VV(G_V(:,l,1) * J_V2(1) + G_V(:,l,2) * J_V2(2) + G_V(:,l,3) * J_V2(3) + G_V(:,l,4) * J_V2(4),J_V1) &
                 - J1J2 * (G_V(1,l,1) + G_V(2,l,2) + G_V(3,l,3) + G_V(4,l,4))                                              &
                 - cont_VV(G_V(:,l,1) * J_V1(1) + G_V(:,l,2) * J_V1(2) + G_V(:,l,3) * J_V1(3) + G_V(:,l,4) * J_V1(4),J_V2)
  end do
end subroutine last_WWV_V


!****************************************************
subroutine last_VWW_V(G_V, J_V1, J_V2, Gtensor)
!----------------------------------------------------
! bare V V V -> V vertex
!****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: G_V(:,:,:), J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  integer           :: l
  complex(REALKIND) :: J1J2

  J1J2 = cont_VV(J_V1, J_V2)
  J1J2 = J1J2 + J1J2
  do l = 1, size(G_V,2)
    Gtensor(l) = J1J2 * (G_V(1,l,1) + G_V(2,l,2) + G_V(3,l,3) + G_V(4,l,4)) &
               - cont_VV(G_V(:,l,1) * J_V1(1) + G_V(:,l,2) * J_V1(2) + G_V(:,l,3) * J_V1(3) + G_V(:,l,4) * J_V1(4),J_V2) &
               - cont_VV(G_V(:,l,1) * J_V2(1) + G_V(:,l,2) * J_V2(2) + G_V(:,l,3) * J_V2(3) + G_V(:,l,4) * J_V2(4),J_V1)
  end do
end subroutine last_VWW_V


!*******************************************************
subroutine last_GH_G(Gin_V, pi, J_S, Gtensor, po)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), pi(4), J_S(4), po(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: s2, pi2(4), po2(4), poS(4), pio2(4), pipoS, Gpo(4), Gin2(10), GinS
  integer :: l

  Gtensor = 0
  s2 = (0.5_/**/REALKIND * J_S(1))
  pi2 = s2 * pi
  po2 = s2 * po
  poS = J_S(1) * po
  pio2(1) = - (pi2(2)+po2(2)) ! -S/2*(pi+pi)_rho, covariant
  pio2(2) = - (pi2(1)+po2(1))
  pio2(3) =   (pi2(4)+po2(4))
  pio2(4) =   (pi2(3)+po2(3))
  pipoS = cont_VV(pi,poS)

  do l = 1, size(Gin_V,2)
    Gpo(1) = cont_VV(Gin_V(:,l,1),poS)
    Gpo(2) = cont_VV(Gin_V(:,l,2),poS)
    Gpo(3) = cont_VV(Gin_V(:,l,3),poS)
    Gpo(4) = cont_VV(Gin_V(:,l,4),poS)
    ! Gout(rank+0) = Gin.po*pi^nu - pi.po*Gin^nu
    Gtensor(l) = Gtensor(l) + Gpo(1)*pi(1) + Gpo(2)*pi(2) + Gpo(3)*pi(3) + Gpo(4)*pi(4) &
                            - pipoS * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    ! Gout(rank+1) = [ Gin_rho*pi^nu + Gin.po*g^nu_rho - (pi+po)_rho*Gin^nu ] * l^rho
    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + Gin_V(2,l,1) * pi2(1) + Gpo(1) + Gin_V(1,l,1) * pio2(1) &
                                        + Gin_V(2,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(1) &
                                        + Gin_V(2,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(1) &
                                        + Gin_V(2,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(1)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) + Gin_V(1,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(2) &
                                        + Gin_V(1,l,2) * pi2(2) + Gpo(2) + Gin_V(2,l,2) * pio2(2) &
                                        + Gin_V(1,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(2) &
                                        + Gin_V(1,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(2)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) - Gin_V(4,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(3) &
                                        - Gin_V(4,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(3) &
                                        - Gin_V(4,l,3) * pi2(3) + Gpo(3) + Gin_V(3,l,3) * pio2(3) &
                                        - Gin_V(4,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(3)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) - Gin_V(3,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(4) &
                                        - Gin_V(3,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(4) &
                                        - Gin_V(3,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(4) &
                                        - Gin_V(3,l,4) * pi2(4) + Gpo(4) + Gin_V(4,l,4) * pio2(4)
    ! Gout(rank+2) = [ Gin_rho*g^nu_sigma - Gin^nu*g_rho_sigma ] l^rho l^sigma
    Gin2(1) =   s2 *  Gin_V(2,l,1)
    Gin2(2) =   s2 * (Gin_V(2,l,2) + Gin_V(1,l,1))
    Gin2(3) =   s2 * (Gin_V(2,l,3) - Gin_V(4,l,1))
    Gin2(4) =   s2 * (Gin_V(2,l,4) - Gin_V(3,l,1))
    Gin2(5) =   s2 *  Gin_V(1,l,2)
    Gin2(6) =   s2 * (Gin_V(1,l,3) - Gin_V(4,l,2))
    Gin2(7) =   s2 * (Gin_V(1,l,4) - Gin_V(3,l,2))
    Gin2(8) = - s2 *  Gin_V(4,l,3)
    Gin2(9) = - s2 * (Gin_V(4,l,4) + Gin_V(3,l,3))
    Gin2(10)= - s2 *  Gin_V(3,l,4)
    GinS = J_S(1) * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    Gtensor(HR(1,HR(1,l))) = Gtensor(HR(1,HR(1,l))) + Gin2(1)
    Gtensor(HR(2,HR(1,l))) = Gtensor(HR(2,HR(1,l))) + Gin2(2) - GinS
    Gtensor(HR(3,HR(1,l))) = Gtensor(HR(3,HR(1,l))) + Gin2(3)
    Gtensor(HR(4,HR(1,l))) = Gtensor(HR(4,HR(1,l))) + Gin2(4)
    Gtensor(HR(2,HR(2,l))) = Gtensor(HR(2,HR(2,l))) + Gin2(5)
    Gtensor(HR(3,HR(2,l))) = Gtensor(HR(3,HR(2,l))) + Gin2(6)
    Gtensor(HR(4,HR(2,l))) = Gtensor(HR(4,HR(2,l))) + Gin2(7)
    Gtensor(HR(3,HR(3,l))) = Gtensor(HR(3,HR(3,l))) + Gin2(8)
    Gtensor(HR(4,HR(3,l))) = Gtensor(HR(4,HR(3,l))) + Gin2(9) + GinS
    Gtensor(HR(4,HR(4,l))) = Gtensor(HR(4,HR(4,l))) + Gin2(10)
  end do
end subroutine last_GH_G


!*******************************************************
subroutine last_GHG_G(Gin_V, plin, J_S, J_V, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), plin(4), J_S(4), J_V(4), p3(4), plout(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV(4), A(4), B(4), C, Jhalf(4), Jtwo(4), Ptmp(4,3)
  integer           :: l

  Gtensor = 0
  SV = J_S(1) * J_V
  Ptmp(:,1) = p3 + plout
  Ptmp(:,2) = plin + plout
  Ptmp(:,3) = plin - p3
  C = -cont_VV(Ptmp(:,2),SV)
  Jhalf = 0.5_/**/REALKIND * SV
  ! covariant components of -2*S*J_V_mu = -2*S*g_(mu,nu)*J_V^nu, factor 2 simplifies in metric tensor
  Jtwo(1) = -SV(2)
  Jtwo(2) = -SV(1)
  Jtwo(3) =  SV(4)
  Jtwo(4) =  SV(3)

  do l = 1, size(Gin_V,2)
    A(1) = cont_VV(Gin_V(:,l,1), SV)
    A(2) = cont_VV(Gin_V(:,l,2), SV)
    A(3) = cont_VV(Gin_V(:,l,3), SV)
    A(4) = cont_VV(Gin_V(:,l,4), SV)
    B(1) = cont_VV(Gin_V(:,l,1), Ptmp(:,1))
    B(2) = cont_VV(Gin_V(:,l,2), Ptmp(:,1))
    B(3) = cont_VV(Gin_V(:,l,3), Ptmp(:,1))
    B(4) = cont_VV(Gin_V(:,l,4), Ptmp(:,1))

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + A(1) + Gin_V(2,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(1) & ! alpha = beta = 1
                                               + Gin_V(2,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(1) & ! alpha = beta = 2
                                               + Gin_V(2,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(1) & ! alpha = beta = 3
                                               + Gin_V(2,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(1)   ! alpha = beta = 4

    Gtensor(HR(2,l)) = Gtensor(HR(2,l))        + Gin_V(1,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(2) & ! alpha = beta = 1
                                        + A(2) + Gin_V(1,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(2) & ! alpha = beta = 2
                                               + Gin_V(1,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(2) & ! alpha = beta = 3
                                               + Gin_V(1,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(2)   ! alpha = beta = 4

    Gtensor(HR(3,l)) = Gtensor(HR(3,l))        - Gin_V(4,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(3) & ! alpha = beta = 1
                                               - Gin_V(4,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(3) & ! alpha = beta = 2
                                        + A(3) - Gin_V(4,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(3) & ! alpha = beta = 3
                                               - Gin_V(4,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(3)   ! alpha = beta = 4

    Gtensor(HR(4,l)) = Gtensor(HR(4,l))        - Gin_V(3,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(4) & ! alpha = beta = 1
                                               - Gin_V(3,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(4) & ! alpha = beta = 2
                                               - Gin_V(3,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(4) & ! alpha = beta = 3
                                        + A(4) - Gin_V(3,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(4)   ! alpha = beta = 4

    Gtensor(l) = Gtensor(l) + A(1)*Ptmp(1,3) + A(2)*Ptmp(2,3) + A(3)*Ptmp(3,3) + A(4)*Ptmp(4,3) &
                            + B(1)*SV(1) + B(2)*SV(2) + B(3)*SV(3) + B(4)*SV(4) &
                            + C * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
  end do
end subroutine last_GHG_G


!*******************************************************
subroutine last_GHGG_G_12(Gin, J_S, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV4(4), SV3V4
  integer :: l
  SV4 = J_S(1) * J_V4
  SV3V4 = cont_VV(J_V3,SV4)
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V3(1)+Gin(:,l,2)*J_V3(2)+Gin(:,l,3)*J_V3(3)+Gin(:,l,4)*J_V3(4), SV4) &
               - SV3V4 * (Gin(1,l,1) + Gin(2,l,2) + Gin(3,l,3) + Gin(4,l,4))
  end do
end subroutine last_GHGG_G_12


!*******************************************************
subroutine last_GHGG_G_23(Gin, J_S, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV3(4)
  integer :: l
  SV3 = J_S(1) * J_V3
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V4(1)+Gin(:,l,2)*J_V4(2)+Gin(:,l,3)*J_V4(3)+Gin(:,l,4)*J_V4(4),  SV3) &
               - cont_VV(Gin(:,l,1)*SV3(1) +Gin(:,l,2)* SV3(2)+Gin(:,l,3)* SV3(3)+Gin(:,l,4)* SV3(4), J_V4)
  end do
end subroutine last_GHGG_G_23


!*******************************************************
subroutine last_GHH_G(Gin_V, pi, J_S1, J_S2, Gtensor, po)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), pi(4), J_S1(4), J_S2(4), po(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: s2, pi2(4), po2(4), poS(4), pio2(4), pipoS, Gpo(4), Gin2(10), GinS
  integer :: l

  Gtensor = 0
  s2 = (0.5_/**/REALKIND * J_S1(1) * J_S2(1))
  pi2 = s2 * pi
  po2 = s2 * po
  poS = J_S1(1) * J_S2(1) * po
  pio2(1) = - (pi2(2)+po2(2)) ! -S/2*(pi+pi)_rho, covariant
  pio2(2) = - (pi2(1)+po2(1))
  pio2(3) =   (pi2(4)+po2(4))
  pio2(4) =   (pi2(3)+po2(3))
  pipoS = cont_VV(pi,poS)

  do l = 1, size(Gin_V,2)
    Gpo(1) = cont_VV(Gin_V(:,l,1),poS)
    Gpo(2) = cont_VV(Gin_V(:,l,2),poS)
    Gpo(3) = cont_VV(Gin_V(:,l,3),poS)
    Gpo(4) = cont_VV(Gin_V(:,l,4),poS)
    ! Gout(rank+0) = Gin.po*pi^nu - pi.po*Gin^nu
    Gtensor(l) = Gtensor(l) + Gpo(1)*pi(1) + Gpo(2)*pi(2) + Gpo(3)*pi(3) + Gpo(4)*pi(4) &
                            - pipoS * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    ! Gout(rank+1) = [ Gin_rho*pi^nu + Gin.po*g^nu_rho - (pi+po)_rho*Gin^nu ] * l^rho
    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + Gin_V(2,l,1) * pi2(1) + Gpo(1) + Gin_V(1,l,1) * pio2(1) &
                                        + Gin_V(2,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(1) &
                                        + Gin_V(2,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(1) &
                                        + Gin_V(2,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(1)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) + Gin_V(1,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(2) &
                                        + Gin_V(1,l,2) * pi2(2) + Gpo(2) + Gin_V(2,l,2) * pio2(2) &
                                        + Gin_V(1,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(2) &
                                        + Gin_V(1,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(2)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) - Gin_V(4,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(3) &
                                        - Gin_V(4,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(3) &
                                        - Gin_V(4,l,3) * pi2(3) + Gpo(3) + Gin_V(3,l,3) * pio2(3) &
                                        - Gin_V(4,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(3)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) - Gin_V(3,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(4) &
                                        - Gin_V(3,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(4) &
                                        - Gin_V(3,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(4) &
                                        - Gin_V(3,l,4) * pi2(4) + Gpo(4) + Gin_V(4,l,4) * pio2(4)
    ! Gout(rank+2) = [ Gin_rho*g^nu_sigma - Gin^nu*g_rho_sigma ] l^rho l^sigma
    Gin2(1) =   s2 *  Gin_V(2,l,1)
    Gin2(2) =   s2 * (Gin_V(2,l,2) + Gin_V(1,l,1))
    Gin2(3) =   s2 * (Gin_V(2,l,3) - Gin_V(4,l,1))
    Gin2(4) =   s2 * (Gin_V(2,l,4) - Gin_V(3,l,1))
    Gin2(5) =   s2 *  Gin_V(1,l,2)
    Gin2(6) =   s2 * (Gin_V(1,l,3) - Gin_V(4,l,2))
    Gin2(7) =   s2 * (Gin_V(1,l,4) - Gin_V(3,l,2))
    Gin2(8) = - s2 *  Gin_V(4,l,3)
    Gin2(9) = - s2 * (Gin_V(4,l,4) + Gin_V(3,l,3))
    Gin2(10)= - s2 *  Gin_V(3,l,4)
    GinS = J_S1(1) * J_S2(1) * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    Gtensor(HR(1,HR(1,l))) = Gtensor(HR(1,HR(1,l))) + Gin2(1)
    Gtensor(HR(2,HR(1,l))) = Gtensor(HR(2,HR(1,l))) + Gin2(2) - GinS
    Gtensor(HR(3,HR(1,l))) = Gtensor(HR(3,HR(1,l))) + Gin2(3)
    Gtensor(HR(4,HR(1,l))) = Gtensor(HR(4,HR(1,l))) + Gin2(4)
    Gtensor(HR(2,HR(2,l))) = Gtensor(HR(2,HR(2,l))) + Gin2(5)
    Gtensor(HR(3,HR(2,l))) = Gtensor(HR(3,HR(2,l))) + Gin2(6)
    Gtensor(HR(4,HR(2,l))) = Gtensor(HR(4,HR(2,l))) + Gin2(7)
    Gtensor(HR(3,HR(3,l))) = Gtensor(HR(3,HR(3,l))) + Gin2(8)
    Gtensor(HR(4,HR(3,l))) = Gtensor(HR(4,HR(3,l))) + Gin2(9) + GinS
    Gtensor(HR(4,HR(4,l))) = Gtensor(HR(4,HR(4,l))) + Gin2(10)
  end do
end subroutine last_GHH_G


!*******************************************************
subroutine last_GHHG_G(Gin_V, plin, J_S1, J_S2, J_V, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), plin(4), J_S1(4), J_S2(4), J_V(4), p3(4), plout(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV(4), A(4), B(4), C, Jhalf(4), Jtwo(4), Ptmp(4,3)
  integer           :: l

  Gtensor = 0
  SV = J_S1(1) * J_S2(1) * J_V
  Ptmp(:,1) = p3 + plout
  Ptmp(:,2) = plin + plout
  Ptmp(:,3) = plin - p3
  C = -cont_VV(Ptmp(:,2),SV)
  Jhalf = 0.5_/**/REALKIND * SV
  ! covariant components of -2*S*J_V_mu = -2*S*g_(mu,nu)*J_V^nu, factor 2 simplifies in metric tensor
  Jtwo(1) = -SV(2)
  Jtwo(2) = -SV(1)
  Jtwo(3) =  SV(4)
  Jtwo(4) =  SV(3)

  do l = 1, size(Gin_V,2)
    A(1) = cont_VV(Gin_V(:,l,1), SV)
    A(2) = cont_VV(Gin_V(:,l,2), SV)
    A(3) = cont_VV(Gin_V(:,l,3), SV)
    A(4) = cont_VV(Gin_V(:,l,4), SV)
    B(1) = cont_VV(Gin_V(:,l,1), Ptmp(:,1))
    B(2) = cont_VV(Gin_V(:,l,2), Ptmp(:,1))
    B(3) = cont_VV(Gin_V(:,l,3), Ptmp(:,1))
    B(4) = cont_VV(Gin_V(:,l,4), Ptmp(:,1))

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + A(1) + Gin_V(2,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(1) & ! alpha = beta = 1
                                               + Gin_V(2,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(1) & ! alpha = beta = 2
                                               + Gin_V(2,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(1) & ! alpha = beta = 3
                                               + Gin_V(2,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(1)   ! alpha = beta = 4

    Gtensor(HR(2,l)) = Gtensor(HR(2,l))        + Gin_V(1,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(2) & ! alpha = beta = 1
                                        + A(2) + Gin_V(1,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(2) & ! alpha = beta = 2
                                               + Gin_V(1,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(2) & ! alpha = beta = 3
                                               + Gin_V(1,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(2)   ! alpha = beta = 4

    Gtensor(HR(3,l)) = Gtensor(HR(3,l))        - Gin_V(4,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(3) & ! alpha = beta = 1
                                               - Gin_V(4,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(3) & ! alpha = beta = 2
                                        + A(3) - Gin_V(4,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(3) & ! alpha = beta = 3
                                               - Gin_V(4,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(3)   ! alpha = beta = 4

    Gtensor(HR(4,l)) = Gtensor(HR(4,l))        - Gin_V(3,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(4) & ! alpha = beta = 1
                                               - Gin_V(3,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(4) & ! alpha = beta = 2
                                               - Gin_V(3,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(4) & ! alpha = beta = 3
                                        + A(4) - Gin_V(3,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(4)   ! alpha = beta = 4

    Gtensor(l) = Gtensor(l) + A(1)*Ptmp(1,3) + A(2)*Ptmp(2,3) + A(3)*Ptmp(3,3) + A(4)*Ptmp(4,3) &
                            + B(1)*SV(1) + B(2)*SV(2) + B(3)*SV(3) + B(4)*SV(4) &
                            + C * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
  end do
end subroutine last_GHHG_G


!*******************************************************
subroutine last_GHHGG_G_12(Gin, J_S1, J_S2, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S1(4), J_S2(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV4(4), SV3V4
  integer :: l
  SV4 = J_S1(1) * J_S2(1) * J_V4
  SV3V4 = cont_VV(J_V3,SV4)
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V3(1)+Gin(:,l,2)*J_V3(2)+Gin(:,l,3)*J_V3(3)+Gin(:,l,4)*J_V3(4), SV4) &
               - SV3V4 * (Gin(1,l,1) + Gin(2,l,2) + Gin(3,l,3) + Gin(4,l,4))
  end do
end subroutine last_GHHGG_G_12


!*******************************************************
subroutine last_GHHGG_G_23(Gin, J_S1, J_S2, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S1(4), J_S2(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV3(4)
  integer :: l
  SV3 = J_S1(1) * J_S2(1) * J_V3
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V4(1)+Gin(:,l,2)*J_V4(2)+Gin(:,l,3)*J_V4(3)+Gin(:,l,4)*J_V4(4),  SV3) &
               - cont_VV(Gin(:,l,1)*SV3(1) +Gin(:,l,2)* SV3(2)+Gin(:,l,3)* SV3(3)+Gin(:,l,4)* SV3(4), J_V4)
  end do
end subroutine last_GHHGG_G_23


!*******************************************************
subroutine last_GHHH_G(Gin_V, pi, J_S1, J_S2, J_S3, Gtensor, po)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), pi(4), J_S1(4), J_S2(4), J_S3(4), po(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: s2, pi2(4), po2(4), poS(4), pio2(4), pipoS, Gpo(4), Gin2(10), GinS
  integer :: l

  Gtensor = 0
  s2 = (0.5_/**/REALKIND * J_S1(1) * J_S2(1) * J_S3(1))
  pi2 = s2 * pi
  po2 = s2 * po
  poS = J_S1(1) * J_S2(1) * J_S3(1) * po
  pio2(1) = - (pi2(2)+po2(2)) ! -S/2*(pi+pi)_rho, covariant
  pio2(2) = - (pi2(1)+po2(1))
  pio2(3) =   (pi2(4)+po2(4))
  pio2(4) =   (pi2(3)+po2(3))
  pipoS = cont_VV(pi,poS)

  do l = 1, size(Gin_V,2)
    Gpo(1) = cont_VV(Gin_V(:,l,1),poS)
    Gpo(2) = cont_VV(Gin_V(:,l,2),poS)
    Gpo(3) = cont_VV(Gin_V(:,l,3),poS)
    Gpo(4) = cont_VV(Gin_V(:,l,4),poS)
    ! Gout(rank+0) = Gin.po*pi^nu - pi.po*Gin^nu
    Gtensor(l) = Gtensor(l) + Gpo(1)*pi(1) + Gpo(2)*pi(2) + Gpo(3)*pi(3) + Gpo(4)*pi(4) &
                            - pipoS * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    ! Gout(rank+1) = [ Gin_rho*pi^nu + Gin.po*g^nu_rho - (pi+po)_rho*Gin^nu ] * l^rho
    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + Gin_V(2,l,1) * pi2(1) + Gpo(1) + Gin_V(1,l,1) * pio2(1) &
                                        + Gin_V(2,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(1) &
                                        + Gin_V(2,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(1) &
                                        + Gin_V(2,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(1)
    Gtensor(HR(2,l)) = Gtensor(HR(2,l)) + Gin_V(1,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(2) &
                                        + Gin_V(1,l,2) * pi2(2) + Gpo(2) + Gin_V(2,l,2) * pio2(2) &
                                        + Gin_V(1,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(2) &
                                        + Gin_V(1,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(2)
    Gtensor(HR(3,l)) = Gtensor(HR(3,l)) - Gin_V(4,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(3) &
                                        - Gin_V(4,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(3) &
                                        - Gin_V(4,l,3) * pi2(3) + Gpo(3) + Gin_V(3,l,3) * pio2(3) &
                                        - Gin_V(4,l,4) * pi2(4)          + Gin_V(4,l,4) * pio2(3)
    Gtensor(HR(4,l)) = Gtensor(HR(4,l)) - Gin_V(3,l,1) * pi2(1)          + Gin_V(1,l,1) * pio2(4) &
                                        - Gin_V(3,l,2) * pi2(2)          + Gin_V(2,l,2) * pio2(4) &
                                        - Gin_V(3,l,3) * pi2(3)          + Gin_V(3,l,3) * pio2(4) &
                                        - Gin_V(3,l,4) * pi2(4) + Gpo(4) + Gin_V(4,l,4) * pio2(4)
    ! Gout(rank+2) = [ Gin_rho*g^nu_sigma - Gin^nu*g_rho_sigma ] l^rho l^sigma
    Gin2(1) =   s2 *  Gin_V(2,l,1)
    Gin2(2) =   s2 * (Gin_V(2,l,2) + Gin_V(1,l,1))
    Gin2(3) =   s2 * (Gin_V(2,l,3) - Gin_V(4,l,1))
    Gin2(4) =   s2 * (Gin_V(2,l,4) - Gin_V(3,l,1))
    Gin2(5) =   s2 *  Gin_V(1,l,2)
    Gin2(6) =   s2 * (Gin_V(1,l,3) - Gin_V(4,l,2))
    Gin2(7) =   s2 * (Gin_V(1,l,4) - Gin_V(3,l,2))
    Gin2(8) = - s2 *  Gin_V(4,l,3)
    Gin2(9) = - s2 * (Gin_V(4,l,4) + Gin_V(3,l,3))
    Gin2(10)= - s2 *  Gin_V(3,l,4)
    GinS = J_S1(1) * J_S2(1) * J_S3(1) * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
    Gtensor(HR(1,HR(1,l))) = Gtensor(HR(1,HR(1,l))) + Gin2(1)
    Gtensor(HR(2,HR(1,l))) = Gtensor(HR(2,HR(1,l))) + Gin2(2) - GinS
    Gtensor(HR(3,HR(1,l))) = Gtensor(HR(3,HR(1,l))) + Gin2(3)
    Gtensor(HR(4,HR(1,l))) = Gtensor(HR(4,HR(1,l))) + Gin2(4)
    Gtensor(HR(2,HR(2,l))) = Gtensor(HR(2,HR(2,l))) + Gin2(5)
    Gtensor(HR(3,HR(2,l))) = Gtensor(HR(3,HR(2,l))) + Gin2(6)
    Gtensor(HR(4,HR(2,l))) = Gtensor(HR(4,HR(2,l))) + Gin2(7)
    Gtensor(HR(3,HR(3,l))) = Gtensor(HR(3,HR(3,l))) + Gin2(8)
    Gtensor(HR(4,HR(3,l))) = Gtensor(HR(4,HR(3,l))) + Gin2(9) + GinS
    Gtensor(HR(4,HR(4,l))) = Gtensor(HR(4,HR(4,l))) + Gin2(10)
  end do
end subroutine last_GHHH_G


!*******************************************************
subroutine last_GHHHG_G(Gin_V, plin, J_S1, J_S2, J_S3, J_V, p3, Gtensor, plout)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin_V(:,:,:), plin(4), J_S1(4), J_S2(4), J_S3(4), J_V(4), p3(4), plout(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV(4), A(4), B(4), C, Jhalf(4), Jtwo(4), Ptmp(4,3)
  integer           :: l

  Gtensor = 0
  SV = J_S1(1) * J_S2(1) * J_S3(1) * J_V
  Ptmp(:,1) = p3 + plout
  Ptmp(:,2) = plin + plout
  Ptmp(:,3) = plin - p3
  C = -cont_VV(Ptmp(:,2),SV)
  Jhalf = 0.5_/**/REALKIND * SV
  ! covariant components of -2*S*J_V_mu = -2*S*g_(mu,nu)*J_V^nu, factor 2 simplifies in metric tensor
  Jtwo(1) = -SV(2)
  Jtwo(2) = -SV(1)
  Jtwo(3) =  SV(4)
  Jtwo(4) =  SV(3)

  do l = 1, size(Gin_V,2)
    A(1) = cont_VV(Gin_V(:,l,1), SV)
    A(2) = cont_VV(Gin_V(:,l,2), SV)
    A(3) = cont_VV(Gin_V(:,l,3), SV)
    A(4) = cont_VV(Gin_V(:,l,4), SV)
    B(1) = cont_VV(Gin_V(:,l,1), Ptmp(:,1))
    B(2) = cont_VV(Gin_V(:,l,2), Ptmp(:,1))
    B(3) = cont_VV(Gin_V(:,l,3), Ptmp(:,1))
    B(4) = cont_VV(Gin_V(:,l,4), Ptmp(:,1))

    Gtensor(HR(1,l)) = Gtensor(HR(1,l)) + A(1) + Gin_V(2,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(1) & ! alpha = beta = 1
                                               + Gin_V(2,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(1) & ! alpha = beta = 2
                                               + Gin_V(2,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(1) & ! alpha = beta = 3
                                               + Gin_V(2,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(1)   ! alpha = beta = 4

    Gtensor(HR(2,l)) = Gtensor(HR(2,l))        + Gin_V(1,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(2) & ! alpha = beta = 1
                                        + A(2) + Gin_V(1,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(2) & ! alpha = beta = 2
                                               + Gin_V(1,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(2) & ! alpha = beta = 3
                                               + Gin_V(1,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(2)   ! alpha = beta = 4

    Gtensor(HR(3,l)) = Gtensor(HR(3,l))        - Gin_V(4,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(3) & ! alpha = beta = 1
                                               - Gin_V(4,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(3) & ! alpha = beta = 2
                                        + A(3) - Gin_V(4,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(3) & ! alpha = beta = 3
                                               - Gin_V(4,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(3)   ! alpha = beta = 4

    Gtensor(HR(4,l)) = Gtensor(HR(4,l))        - Gin_V(3,l,1)*Jhalf(1) + Gin_V(1,l,1)*Jtwo(4) & ! alpha = beta = 1
                                               - Gin_V(3,l,2)*Jhalf(2) + Gin_V(2,l,2)*Jtwo(4) & ! alpha = beta = 2
                                               - Gin_V(3,l,3)*Jhalf(3) + Gin_V(3,l,3)*Jtwo(4) & ! alpha = beta = 3
                                        + A(4) - Gin_V(3,l,4)*Jhalf(4) + Gin_V(4,l,4)*Jtwo(4)   ! alpha = beta = 4

    Gtensor(l) = Gtensor(l) + A(1)*Ptmp(1,3) + A(2)*Ptmp(2,3) + A(3)*Ptmp(3,3) + A(4)*Ptmp(4,3) &
                            + B(1)*SV(1) + B(2)*SV(2) + B(3)*SV(3) + B(4)*SV(4) &
                            + C * (Gin_V(1,l,1) + Gin_V(2,l,2) + Gin_V(3,l,3) + Gin_V(4,l,4))
  end do
end subroutine last_GHHHG_G


!*******************************************************
subroutine last_GHHHGG_G_12(Gin, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S1(4), J_S2(4), J_S3(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV4(4), SV3V4
  integer :: l
  SV4 = J_S1(1) * J_S2(1) * J_S3(1) * J_V4
  SV3V4 = cont_VV(J_V3,SV4)
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V3(1)+Gin(:,l,2)*J_V3(2)+Gin(:,l,3)*J_V3(3)+Gin(:,l,4)*J_V3(4), SV4) &
               - SV3V4 * (Gin(1,l,1) + Gin(2,l,2) + Gin(3,l,3) + Gin(4,l,4))
  end do
end subroutine last_GHHHGG_G_12


!*******************************************************
subroutine last_GHHHGG_G_23(Gin, J_S1, J_S2, J_S3, J_V3, J_V4, Gtensor)
!*******************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), J_S1(4), J_S2(4), J_S3(4), J_V3(4), J_V4(4)
  complex(REALKIND), intent(out) :: Gtensor(:)
  complex(REALKIND) :: SV3(4)
  integer :: l
  SV3 = J_S1(1) * J_S2(1) * J_S3(1) * J_V3
  do l = 1, size(Gin,2)
    Gtensor(l) = cont_VV(Gin(:,l,1)*J_V4(1)+Gin(:,l,2)*J_V4(2)+Gin(:,l,3)*J_V4(3)+Gin(:,l,4)*J_V4(4),  SV3) &
               - cont_VV(Gin(:,l,1)*SV3(1) +Gin(:,l,2)* SV3(2)+Gin(:,l,3)* SV3(3)+Gin(:,l,4)* SV3(4), J_V4)
  end do
end subroutine last_GHHHGG_G_23




end module ol_last_step_/**/REALKIND
