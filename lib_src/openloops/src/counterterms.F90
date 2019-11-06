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


module ol_counterterms_/**/REALKIND

  ! TODO:  <10-11-18, J.-N. Lang> !
  ! make copies of other 2point ct structures

  interface counter_Q_A
    module procedure counter_Q_A_orig, counter_Q_A_pid
  end interface

  interface counter_A_Q
    module procedure counter_A_Q_orig, counter_A_Q_pid
  end interface

  interface counter_V_V
    module procedure counter_V_V_orig, counter_V_V_pid
  end interface

contains

! **********************************************************************
subroutine counter_Q_A_orig(ctQA, J_Q, K, Jout_Q)
! Q -> Q counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(2), J_Q(4), K(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = ctQA(1) * ( - K(2)*J_Q(3) + K(4)*J_Q(4)) - ctQA(2) * J_Q(1)
  Jout_Q(2) = ctQA(1) * ( - K(1)*J_Q(4) + K(3)*J_Q(3)) - ctQA(2) * J_Q(2)
  Jout_Q(3) = ctQA(1) * ( - K(1)*J_Q(1) - K(4)*J_Q(2)) - ctQA(2) * J_Q(3)
  Jout_Q(4) = ctQA(1) * ( - K(2)*J_Q(2) - K(3)*J_Q(1)) - ctQA(2) * J_Q(4)
end subroutine counter_Q_A_orig
! **********************************************************************
subroutine counter_Q_A_pid(ctQA, pid, J_Q, K, Jout_Q)
! Q -> Q counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(2), J_Q(4), K(4)
  integer,           intent(in)  :: pid
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = ctQA(1) * ( - K(2)*J_Q(3) + K(4)*J_Q(4)) - ctQA(2) * J_Q(1)
  Jout_Q(2) = ctQA(1) * ( - K(1)*J_Q(4) + K(3)*J_Q(3)) - ctQA(2) * J_Q(2)
  Jout_Q(3) = ctQA(1) * ( - K(1)*J_Q(1) - K(4)*J_Q(2)) - ctQA(2) * J_Q(3)
  Jout_Q(4) = ctQA(1) * ( - K(2)*J_Q(2) - K(3)*J_Q(1)) - ctQA(2) * J_Q(4)
end subroutine counter_Q_A_pid


! **********************************************************************
subroutine counter_A_Q_orig(ctQA, J_A, K, Jout_A)
! A -> A counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(2), J_A(4), K(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = ctQA(1) * ( + K(1)*J_A(3) + K(3)*J_A(4)) - ctQA(2) * J_A(1)
  Jout_A(2) = ctQA(1) * ( + K(2)*J_A(4) + K(4)*J_A(3)) - ctQA(2) * J_A(2)
  Jout_A(3) = ctQA(1) * ( + K(2)*J_A(1) - K(3)*J_A(2)) - ctQA(2) * J_A(3)
  Jout_A(4) = ctQA(1) * ( + K(1)*J_A(2) - K(4)*J_A(1)) - ctQA(2) * J_A(4)
end subroutine counter_A_Q_orig

! **********************************************************************
subroutine counter_A_Q_pid(ctQA, pid, J_A, K, Jout_A)
! A -> A counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(2), J_A(4), K(4)
  integer,           intent(in)  :: pid
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = ctQA(1) * ( + K(1)*J_A(3) + K(3)*J_A(4)) - ctQA(2) * J_A(1)
  Jout_A(2) = ctQA(1) * ( + K(2)*J_A(4) + K(4)*J_A(3)) - ctQA(2) * J_A(2)
  Jout_A(3) = ctQA(1) * ( + K(2)*J_A(1) - K(3)*J_A(2)) - ctQA(2) * J_A(3)
  Jout_A(4) = ctQA(1) * ( + K(1)*J_A(2) - K(4)*J_A(1)) - ctQA(2) * J_A(4)
end subroutine counter_A_Q_pid


! **********************************************************************
subroutine counter_V_V_orig(ctVV, J_V, K, Jout_V)
! V -> V counter term
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: ctVV(3), J_V(4), K(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = (ctVV(1)*cont_VV(K,K) + ctVV(2)) * J_V + ctVV(3) * (cont_VV(J_V,K)) * K
end subroutine counter_V_V_orig

subroutine counter_V_V_pid(ctVV, pid, J_V, K, Jout_V)
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: ctVV(3), J_V(4), K(4)
  integer,           intent(in) :: pid
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = (ctVV(1)*cont_VV(K,K) + ctVV(2)) * J_V + ctVV(3) * (cont_VV(J_V,K)) * K
end subroutine counter_V_V_pid


! **********************************************************************
subroutine counter_S_S(ctSS, J_S, K, Jout_S)
! S -> S counter term
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_V
  implicit none
  complex(REALKIND), intent(in)  :: ctSS(2), J_S(4), K(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = (ctSS(1)*cont_V(K) - ctSS(2)) * J_S(1)
end subroutine counter_S_S

! **********************************************************************
subroutine counter_S_V(ctSV, J_S, K, Jout_V)
! S -> V counter term
! i k
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_V
  implicit none
  complex(REALKIND), intent(in)  :: ctSV, J_S(4), K(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V =  - ctSV * K * J_S(1)
end subroutine counter_S_V


! **********************************************************************
subroutine counter_V_S(ctSV, J_V, K, Jout_S)
! S -> V counter term
! -i k
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_V
  implicit none
  complex(REALKIND), intent(in)  :: ctSV, J_V(4), K(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) =  ctSV * cont_V(K)
end subroutine counter_V_S


! **********************************************************************
subroutine counter_Q_A_LR(ctQA, J_Q, K, Jout_Q)
! Q -> Q counter term with left/right splitting
! ctQA(1)*slash(k)*P_R+ctQA(2)*slash(k)*P_L - ctQA(3) P_R - ctQA(4) P_L
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(4), J_Q(4), K(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = ctQA(2) * ( - K(2)*J_Q(3) + K(4)*J_Q(4)) - ctQA(3) * J_Q(1)
  Jout_Q(2) = ctQA(2) * ( - K(1)*J_Q(4) + K(3)*J_Q(3)) - ctQA(3) * J_Q(2)
  Jout_Q(3) = ctQA(1) * ( - K(1)*J_Q(1) - K(4)*J_Q(2)) - ctQA(4) * J_Q(3)
  Jout_Q(4) = ctQA(1) * ( - K(2)*J_Q(2) - K(3)*J_Q(1)) - ctQA(4) * J_Q(4)
end subroutine counter_Q_A_LR



! **********************************************************************
subroutine counter_A_Q_LR(ctQA, J_A, K, Jout_A)
! A -> A counter term with left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: ctQA(4), J_A(4), K(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = ctQA(1) * ( + K(1)*J_A(3) + K(3)*J_A(4)) - ctQA(3) * J_A(1)
  Jout_A(2) = ctQA(1) * ( + K(2)*J_A(4) + K(4)*J_A(3)) - ctQA(3) * J_A(2)
  Jout_A(3) = ctQA(2) * ( + K(2)*J_A(1) - K(3)*J_A(2)) - ctQA(4) * J_A(3)
  Jout_A(4) = ctQA(2) * ( + K(1)*J_A(2) - K(4)*J_A(1)) - ctQA(4) * J_A(4)
end subroutine counter_A_Q_LR



! ======================================================================
! Vertex counter terms.
! For QCD corrections, these are just copies of the vertex routines.
! ======================================================================


! **********************************************************************
subroutine counter_ZQ_A(g_RL, J_Z, J_Q, Jout_Q)
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! J_Q(4)     = incoming quark current
! J_Z(4)     = incoming Z current ("light-cone" rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_Q(4)  = outgoing quark current
! Jout_Q(i)  = J_Z(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_Q(4), J_Z(4), Jout_Q(4)
  complex(REALKIND) :: g_RL(2)
  Jout_Q(1) = g_RL(2) * ( - J_Z(2)*J_Q(3) + J_Z(4)*J_Q(4))
  Jout_Q(2) = g_RL(2) * ( - J_Z(1)*J_Q(4) + J_Z(3)*J_Q(3))
  Jout_Q(3) = g_RL(1) * ( - J_Z(1)*J_Q(1) - J_Z(4)*J_Q(2))
  Jout_Q(4) = g_RL(1) * ( - J_Z(2)*J_Q(2) - J_Z(3)*J_Q(1))
end subroutine counter_ZQ_A


! **********************************************************************
subroutine counter_AZ_Q(g_RL, J_A, J_Z, Jout_A)
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! J_A(4)     = incoming anti-quark current
! J_Z(4)     = incoming Z current (light-cone rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_A(4)  = outgoing anti-quark current
! Jout_A(i)  = J_A(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * J_Z(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_A(4), J_Z(4), Jout_A(4)
  complex(REALKIND) :: g_RL(2)
  Jout_A(1) = g_RL(1) * ( - J_Z(1)*J_A(3) - J_Z(3)*J_A(4))
  Jout_A(2) = g_RL(1) * ( - J_Z(2)*J_A(4) - J_Z(4)*J_A(3))
  Jout_A(3) = g_RL(2) * ( - J_Z(2)*J_A(1) + J_Z(3)*J_A(2))
  Jout_A(4) = g_RL(2) * ( - J_Z(1)*J_A(2) + J_Z(4)*J_A(1))
end subroutine counter_AZ_Q


! **********************************************************************
subroutine counter_QA_Z(g_RL, J_Q, J_A, Jout_Z)
! bare QA -> Z Z-like interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! g_RL(1)   = right-handed coupling gR
! g_RL(2)   = left-handed coupling gL
! Jout_Z(4) = outgoing Z current (light-cone rep.)
! Jout_Z(A) = J_A(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_Q(4), J_A(4), Jout_Z(4)
  complex(REALKIND) :: g_RL(2)
  Jout_Z(1) = - g_RL(2)*J_A(1)*J_Q(3) - g_RL(1)*J_A(4)*J_Q(2)
  Jout_Z(2) = - g_RL(2)*J_A(2)*J_Q(4) - g_RL(1)*J_A(3)*J_Q(1)
  Jout_Z(3) = - g_RL(2)*J_A(1)*J_Q(4) + g_RL(1)*J_A(3)*J_Q(2)
  Jout_Z(4) = - g_RL(2)*J_A(2)*J_Q(3) + g_RL(1)*J_A(4)*J_Q(1)
  Jout_Z = Jout_Z + Jout_Z
end subroutine counter_QA_Z


! **********************************************************************
subroutine counter_WQ_A(J_W, J_Q, Jout_Q)
! bare WQ -> Q W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = incoming quark current
! J_W(4)    = incoming W current ("light-cone" rep.)
! Jout_Q(4) = outgoing quark current
! Jout_Q(i) = J_W(A) * [gamma_A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_Q(4), J_W(4), Jout_Q(4)
  Jout_Q(1)   = - J_W(2)*J_Q(3) + J_W(4)*J_Q(4)
  Jout_Q(2)   = - J_W(1)*J_Q(4) + J_W(3)*J_Q(3)
  Jout_Q(3:4) = 0
end subroutine counter_WQ_A


! **********************************************************************
subroutine counter_AW_Q(J_A, J_W, Jout_A)
! bare AW -> A W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_W(4)    = incoming W current (light-cone rep.)
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * [gamma_A*w_L](j,i) * J_W(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_A(4), J_W(4), Jout_A(4)
  Jout_A(1:2) = 0
  Jout_A(3)   = - J_W(2)*J_A(1) + J_W(3)*J_A(2)
  Jout_A(4)   = - J_W(1)*J_A(2) + J_W(4)*J_A(1)
end subroutine counter_AW_Q


! **********************************************************************
subroutine counter_QA_W(J_Q, J_A, Jout_W)
! bare QA -> W W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_W(4) = outgoing W current (light-cone rep.)
! Jout_W(A) = J_A(i) * [gamma^A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_Q(4), J_A(4), Jout_W(4)
  Jout_W(1) = - J_A(1)*J_Q(3)
  Jout_W(2) = - J_A(2)*J_Q(4)
  Jout_W(3) = - J_A(1)*J_Q(4)
  Jout_W(4) = - J_A(2)*J_Q(3)
  Jout_W = Jout_W + Jout_W
end subroutine counter_QA_W


! **********************************************************************
subroutine counter_VQ_A(J_V, J_Q, Jout_Q)
! VQ -> Q counter term; without left/right splitting
! Factorised wrt. vert_VQ_A: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_V(2)*J_Q(3)+J_V(4)*J_Q(4)
  Jout_Q(2) = - J_V(1)*J_Q(4)+J_V(3)*J_Q(3)
  Jout_Q(3) = - J_V(1)*J_Q(1)-J_V(4)*J_Q(2)
  Jout_Q(4) = - J_V(2)*J_Q(2)-J_V(3)*J_Q(1)
end subroutine counter_VQ_A


! **********************************************************************
subroutine counter_AV_Q(J_A, J_V, Jout_A)
! AV -> A counter term; without left/right splitting
! Factorised wrt. vert_AV_Q: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = - J_V(1)*J_A(3) - J_V(3)*J_A(4)
  Jout_A(2) = - J_V(2)*J_A(4) - J_V(4)*J_A(3)
  Jout_A(3) = - J_V(2)*J_A(1) + J_V(3)*J_A(2)
  Jout_A(4) = - J_V(1)*J_A(2) + J_V(4)*J_A(1)
end subroutine counter_AV_Q


! **********************************************************************
subroutine counter_QA_V(J_Q, J_A, Jout_V)
! QA -> V counter term; without left/right splitting
! Factorised wrt. vert_QA_V: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = Jout_V + Jout_V
end subroutine counter_QA_V


! **********************************************************************
subroutine counter_VQ_A_LR(ctVFF,J_V, J_Q, Jout_Q)
! VQ -> Q counter term; with left/right splitting -> ZQ_A
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_V(4), ctVFF(2)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  call counter_ZQ_A(ctVFF,J_V, J_Q, Jout_Q)
end subroutine counter_VQ_A_LR

! **********************************************************************
subroutine counter_AV_Q_LR(ctVFF,J_A, J_V, Jout_A)
! AV -> A counter term; with left/right splitting -> AZ_Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4),ctVFF(2)
  complex(REALKIND), intent(out) :: Jout_A(4)
  call counter_AZ_Q(ctVFF,J_A, J_V, Jout_A)
end subroutine counter_AV_Q_LR

! **********************************************************************
subroutine counter_QA_V_LR(ctVFF,J_Q, J_A, Jout_V)
! QA -> V counter term; with left/right splitting -> QA_Z
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_A(4),ctVFF(2)
  complex(REALKIND), intent(out) :: Jout_V(4)
  call counter_QA_Z(ctVFF,J_Q, J_A, Jout_V)
end subroutine counter_QA_V_LR



! ! **********************************************************************
! subroutine counter_UV_W(J_V1, P1, J_V2, P2, Jout_V)
! ! VV -> V counter term
! ! Factorised wrt. vert_UV_W: ctVVV = dlnG*gQCD + 3/2 * dZg
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   use ol_contractions_/**/REALKIND, only: cont_VV
!   implicit none
!   complex(REALKIND), intent(in)  :: J_V1(4), P1(4), J_V2(4), P2(4)
!   complex(REALKIND), intent(out) :: Jout_V(4)
!   Jout_V = cont_VV(J_V1,J_V2) * (P1 - P2) + cont_VV(P1+P2+P2,J_V1) * J_V2 - cont_VV(P1+P1+P2,J_V2) * J_V1
! end subroutine counter_UV_W
subroutine counter_UV_W(J_V1, P1, J_V2, P2, Jout_V)
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND) :: P1(4), P2(4)
  complex(REALKIND) :: J_V1(4), J_V2(4), Jout_V(4)
  complex(REALKIND) :: J1J2, P1J2, P2J1
  J1J2 = cont_VV(J_V1,J_V2)
  P1J2 = cont_VV(P1+P1+P2,J_V2)
  P2J1 = cont_VV(P1+P2+P2,J_V1)
  Jout_V = J1J2 * (P1 - P2) + P2J1 * J_V2 - P1J2 * J_V1
end subroutine counter_UV_W


! **********************************************************************
subroutine counter_EV_V(J_V1, J_V2, J_V3, Jout_V)
! sigma vertex counter term, where the sigma wave function is replaced by two gluon wave functions J_V1 and J_V2
! Jout_V(d) = (g(a,c)*g(b,d) + g(1,4)*g(2,3)) * J_V1(a) * J_V1(b) * J_V1(c)
!           = J_V1.J_V3 * J_V2(d) + J_V2.J_V3 * J_V1(d)
! Factorised wrt. vert_EV_V: ctVVVV = 1/2*dlnG*gQCD + dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = cont_VV(J_V1,J_V3) * J_V2 - cont_VV(J_V2,J_V3) * J_V1
end subroutine counter_EV_V


! **********************************************************************
subroutine counter_AQ_S(g_RL, J_A, J_Q, Jout_S)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming anti-fermion current: J_A(4)
! Incoming fermion current:      J_Q(4)
! Outgoing scalar current:       Jout_S = gR*J_A.P_R.J_Q + gL*J_A.P_L.J_Q
!   with the right- and left-handed projectors P_R = (1+y5)/2 and P_L = (1-y5)/2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_A(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = g_RL(1) * (J_A(1)*J_Q(1) + J_A(2)*J_Q(2)) + g_RL(2) * (J_A(3)*J_Q(3) + J_A(4)*J_Q(4))
end subroutine counter_AQ_S


! **********************************************************************
subroutine counter_QS_A(g_RL, J_Q, J_S, Jout_A)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming fermion current:      J_Q(4)
! Incoming scalar current:       J_S
! Outgoing anti-fermion current: Jout_A(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_Q(4), J_S(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = g_RL(1) * J_Q(1) * J_S(1)
  Jout_A(2) = g_RL(1) * J_Q(2) * J_S(1)
  Jout_A(3) = g_RL(2) * J_Q(3) * J_S(1)
  Jout_A(4) = g_RL(2) * J_Q(4) * J_S(1)
end subroutine counter_QS_A


! **********************************************************************
subroutine counter_SA_Q(g_RL, J_S, J_A, Jout_Q)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming scalar current:       J_S
! Incoming anti-fermion current: J_A(4)
! Outgoing fermion current:      Jout_Q(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_S(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = g_RL(1) * J_A(1) * J_S(1)
  Jout_Q(2) = g_RL(1) * J_A(2) * J_S(1)
  Jout_Q(3) = g_RL(2) * J_A(3) * J_S(1)
  Jout_Q(4) = g_RL(2) * J_A(4) * J_S(1)
end subroutine counter_SA_Q


! **********************************************************************
subroutine counter_VG_G(J_V, Jin_G, p2, Jout_G, p3)
! Z-gluon-gluon vertex for R2
! Jout_G(c) = ep(a,b,c,d) * J_V(a) * Jin_G(b) * (p2-p3)(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_EpVVV
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), Jin_G(4), p2(4), p3(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  call cont_EpVVV(J_V, Jin_G, p2-p3, Jout_G)
end subroutine counter_VG_G


! **********************************************************************
subroutine counter_GG_V(J_G1, p1, J_G2, p2, Jout_V)
! Z-gluon-gluon vertex for R2
! Jout_V(a) = ep(a,b,c,d) * J_G1(b) * J_G2(c) * (p1-p2)(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_EpVVV
  implicit none
  complex(REALKIND), intent(in)  :: J_G1(4), p1(4), J_G2(4), p2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  call cont_EpVVV(J_G1, J_G2, p1-p2, Jout_V)
end subroutine counter_GG_V


! **********************************************************************
subroutine counter_SG_G(J_S, Jin_G, Jout_G)
! Higgs-gluon-gluon vertex for R2
! Jout_G = J_S * Jin_G
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), Jin_G(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  Jout_G = J_S(1) * Jin_G
end subroutine counter_SG_G


! **********************************************************************
subroutine counter_GG_S(J_G1, J_G2, Jout_S)
! Higgs-gluon-gluon vertex for R2
! Jout_S = J_G1.J_G2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_G1(4), J_G2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = cont_VV(J_G1,J_G2)
end subroutine counter_GG_S


! **********************************************************************
subroutine counter_VVG_G(J_V1, J_V2, Jin_G, Jout_G)
! Vector-vector-gluon-gluon vertex for R2
! Jout_G = J_V1.J_V2 * Jin_G + J_V1.Jin_G * J_V2 + J_V2.Jin_G * J_V1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), Jin_G(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  Jout_G = cont_VV(J_V1,J_V2) * Jin_G + cont_VV(J_V1,Jin_G) * J_V2 + cont_VV(J_V2,Jin_G) * J_V1
  ! Jout_G = Jout_G + 4*(1-g5s)*sum_ij(a_V1QiQj*aV2QiQj) * cont_VV(J_V1,J_V2) * Jin_G
end subroutine counter_VVG_G


! **********************************************************************
subroutine counter_SSG_G(J_S1, J_S2, Jin_G, Jout_G)
! Scalar-scalar-gluon-gluon vertex for R2
! Jout_G = J_S1 * J_S2 * Jin_G
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4), Jin_G(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  Jout_G = J_S1(1) * J_S2(1) * Jin_G
end subroutine counter_SSG_G


! **********************************************************************
subroutine counter_GGS_S(J_G1, J_G2, Jin_S, Jout_S)
! Scalar-scalar-gluon-gluon vertex for R2
! Jout_S = Jin_S * J_G1.J_G2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_G1(4), J_G2(4), Jin_S(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = Jin_S(1) * cont_VV(J_G1,J_G2)
end subroutine counter_GGS_S


! **********************************************************************
subroutine counter_GGG_V(gVA, J_G1, J_G2, J_G3, Jout_V)
! Vector-gluon-gluon-gluon vertex for R2
! Jout_V(a) = [v * (g(a,b)*g(c,d) + g(a,c)*g(b,d) + g(a,d)*g(b,c)) + a * ep(a,b,c,d)] * J_G1(b) * J_G2(c) * J_G3(d)
! v = 2/3*vector_coupling; a = -6*i*axial_coupling
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_EpVVV
  implicit none
  complex(REALKIND), intent(in)  :: gVA(2), J_G1(4), J_G2(4), J_G3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  if (gVA(2) /= 0) then
    call cont_EpVVV(J_G1, J_G2, J_G3, Jout_V)
    Jout_V = gVA(2) * Jout_V
  else
    Jout_V = 0
  end if
  Jout_V = Jout_V + gVA(1) * (cont_VV(J_G1,J_G2)*J_G3 + cont_VV(J_G2,J_G3)*J_G1 + cont_VV(J_G3,J_G1)*J_G2)
end subroutine counter_GGG_V


! **********************************************************************
subroutine counter_VGG_G(gVA, J_V, J_G1, J_G2, Jout_G)
! Vector-gluon-gluon-gluon vertex for R2
! Jout_G(d) = [v * (g(a,b)*g(c,d) + g(a,c)*g(b,d) + g(a,d)*g(b,c)) + a * ep(a,b,c,d)] * J_V(a) * J_G1(b) * J_G2(c)
! v = 2/3*vector_coupling; a = -6*i*axial_coupling
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_EpVVV
  implicit none
  complex(REALKIND), intent(in)  :: gVA(2), J_V(4), J_G1(4), J_G2(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  if (gVA(2) /= 0) then
    call cont_EpVVV(J_V, J_G1, J_G2, Jout_G)
    Jout_G = - gVA(2) * Jout_G
  else
    Jout_G = 0
  end if
  Jout_G = Jout_G + gVA(1) * (cont_VV(J_V,J_G1)*J_G2 + cont_VV(J_G1,J_G2)*J_V + cont_VV(J_G2,J_V)*J_G1)
end subroutine counter_VGG_G


! **********************************************************************
subroutine counter_GGG_G(J_G1, J_G2, J_G3, Jout_G)
! gluon-gluon-gluon-gluon vertex for R2, factorised Lorentz monomials g(a,b)*g(c,d)
! Jout_G(d) = g(a,b)*g(c,d) * J_G1(a) * J_G2(b) * J_G3(c) = J_G1.J_G2 * J_G3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_G1(4), J_G2(4), J_G3(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  Jout_G = cont_VV(J_G1,J_G2) * J_G3
end subroutine counter_GGG_G





! ======================================================================
! Additional vertex counter terms.
! for EW corrections, these are just copies of the vertex routines.
! ======================================================================

! **********************************************************************
subroutine counter_SS_S(J_S1, J_S2, Jout_S)
! Three scalar vertex
! Incoming scalar currents: J_S1, J_S2
! Outgoing scalar current:  Jout_S = J_S1 * J_S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = J_S1(1) * J_S2(1)
end subroutine counter_SS_S

! **********************************************************************
! subroutine vert_VS_S(J_V, J_S, P1, Jout_S)
subroutine counter_VS_T(J_V, P1, J_S, P2, Jout_S)
! Vector boson + two scalars vertex
! Incoming vector current: J_V(4), incoming momentum P1(4) (light-cone rep.)
! Incoming scalar current: J_S,    incoming momentum P2(4) (light-cone rep.)
! Outgoing scalar current: Jout_S = J_V.(2*P2+P1) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), P1(4), J_S(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = cont_VV(P1+P2+P2, J_V) * J_S(1)
end subroutine counter_VS_T


! **********************************************************************
! subroutine vert_SV_S(J_S, P2, J_V, Jout_S)
subroutine counter_TV_S(J_S, P1, J_V, P2, Jout_S)
! Vector boson + two scalars vertex
! Incoming scalar current: J_S, incoming momentum P2(4) (light-cone rep.)
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing scalar current: Jout_S = J_V.(-2*P1-P2) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), P1(4), J_V(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = - cont_VV(P1+P1+P2, J_V) * J_S(1)
end subroutine counter_TV_S


! **********************************************************************
! subroutine vert_SS_V(J_S1, P1, J_S2, P2, Jout_V)
subroutine counter_ST_V(J_S1, P1, J_S2, P2, Jout_V)
! Vector boson + two scalars vertex
! Incoming scalar currents: J_S1, J_S2, incoming momenta P1(4), P2(4) (light-cone rep.)
! Outgoing vector current: Jout_V(4) = J_S1 * J_S2 * (P1 - P2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), P1(4), J_S2(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = (J_S1(1) * J_S2(1)) * (P1 - P2)
! end subroutine vert_SS_V
end subroutine counter_ST_V


! **********************************************************************
subroutine counter_VV_S(J_V1, J_V2, Jout_S)
! Two vector boson + scalar vertex
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = cont_VV(J_V1,J_V2)
end subroutine counter_VV_S


! **********************************************************************
subroutine counter_VS_V(J_V, J_S, Jout_V)
! Two vector boson + scalar vertex
! Incoming vector current: J_V(4) (light-cone rep.)
! Incoming scalar current: J_S
! Outgoing vector current: Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), J_S(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = J_V * J_S(1)
end subroutine counter_VS_V


! **********************************************************************
subroutine counter_SV_V(J_S, J_V, Jout_V)
! Two vector boson + scalar vertex
! Incoming scalar current: J_S
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing vector current: Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = J_S(1) * J_V
end subroutine counter_SV_V




! **********************************************************************
subroutine counter_SSS_S(J_S1, J_S2, J_S3, Jout_S)
! Four scalar vertex
! Incoming scalar currents: J_S1, J_S2, J_S3
! Outgoing scalar current:  Jout_S = J_S1 * J_S2 * J_S3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4), J_S3(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = J_S1(1) * J_S2(1) * J_S3(1)
end subroutine counter_SSS_S



! **********************************************************************
subroutine counter_VVV_V(J_V1, J_V2, J_V3, Jout_V)
! neutral Vector-vector-vector-vector pure R2
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [g(a1,a2)*g(a3,a4) + g(a2,a3)*g(a1,a4)
!               + g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  complex(REALKIND) :: J1J2, J1J3, J2J3

  J1J2 = cont_VV(J_V1, J_V2)
  J1J3 = cont_VV(J_V1, J_V3)
  J2J3 = cont_VV(J_V2, J_V3)
  Jout_V = J1J2 * J_V3 + J2J3 * J_V1 + J1J3 * J_V2

end subroutine counter_VVV_V


! **********************************************************************
subroutine counter_WWV_V(ctWWVV, J_V1, J_V2, J_V3, Jout_V)
! W-W-vector-vector UV & pure R2
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [ctVVVV(1) * g(a1,a2)*g(a3,a4) + ctVVVV(2) * g(a2,a3)*g(a1,a4)
!               + ctVVVV(2) * g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4), ctWWVV(2)
  complex(REALKIND), intent(out) :: Jout_V(4)
  complex(REALKIND) :: J1J2, J1J3, J2J3

  J1J2 = cont_VV(J_V1, J_V2)
  J1J3 = cont_VV(J_V1, J_V3)
  J2J3 = cont_VV(J_V2, J_V3)
  Jout_V = J1J2 * J_V3 * ctWWVV(1) + J2J3 * J_V1 * ctWWVV(2) + J1J3 * J_V2 * ctWWVV(2)

end subroutine counter_WWV_V

! **********************************************************************
subroutine counter_VWW_V(ctWWVV,J_V1,J_V2,J_V3,Jout_V)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4), ctWWVV(2)
  complex(REALKIND), intent(out) :: Jout_V(4)
  call counter_WWV_V(ctWWVV,J_V2,J_V3,J_V1,Jout_V)
end subroutine counter_VWW_V


! ! **********************************************************************
! subroutine counter_WVV_W(ctWWVV, J_V1, J_V2, J_V3, Jout_V)
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   use ol_contractions_/**/REALKIND, only: cont_VV
!   implicit none
!
!   complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4), ctWWVV(3)
!   complex(REALKIND), intent(out) :: Jout_V(4)
!   complex(REALKIND) :: J1J2, J1J3, J2J3
!
!   J1J2 = cont_VV(J_V1, J_V2)
!   J1J3 = cont_VV(J_V1, J_V3)
!   J2J3 = cont_VV(J_V2, J_V3)
!   Jout_V = J1J2 * J_V3 * ctWWVV(2) + J2J3 * J_V1 * ctWWVV(1) + J1J3 * J_V2 * ctWWVV(2)
!
! end subroutine counter_WVV_W
!

! ! **********************************************************************
! subroutine counter_WVW_V(ctWWVV, J_V1, J_V2, J_V3, Jout_V)
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   use ol_contractions_/**/REALKIND, only: cont_VV
!   implicit none
!
!   complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4), ctWWVV(2)
!   complex(REALKIND), intent(out) :: Jout_V(4)
!   complex(REALKIND) :: J1J2, J1J3, J2J3
!
!   J1J2 = cont_VV(J_V1, J_V2)
!   J1J3 = cont_VV(J_V1, J_V3)
!   J2J3 = cont_VV(J_V2, J_V3)
!   Jout_V = J1J2 * J_V3 * ctWWVV(2) + J2J3 * J_V1 * ctWWVV(1) + J1J3 * J_V2 * ctWWVV(1)
!
! end subroutine counter_WVW_V
!

! **********************************************************************
subroutine counter_VVS_S(J_V1, J_V2, J_S, Jout_S)
! Two vector boson + two scalars vertex
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Incoming scalar current:  J_S
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_S(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = cont_VV(J_V1,J_V2) * J_S(1)
end subroutine counter_VVS_S


! **********************************************************************
subroutine counter_SSV_V(J_S1, J_S2, J_V, Jout_V)
! Two vector boson + two scalars vertex
! Incoming scalar currents: J_S1, J_S2
! Incoming vector current:  J_V(4) (light-cone rep.)
! Outgoing vector current:  Jout_V (light-cone rep.)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = (J_S1(1) * J_S2(1)) * J_V
end subroutine counter_SSV_V


! **********************************************************************
subroutine counter_VSS_V(J_V, J_S1, J_S2, Jout_V)
! Two vector boson + two scalars vertex
! Incoming vector current:  J_V(4) (light-cone rep.)
! Incoming scalar currents: J_S1, J_S2
! Outgoing vector current:  Jout_V (light-cone rep.)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), J_S1(4), J_S2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = (J_S1(1) * J_S2(1)) * J_V
end subroutine counter_VSS_V


! **********************************************************************
subroutine counter_SVV_S(J_S, J_V1, J_V2, Jout_S)
! Two vector boson + two scalars vertex
! Incoming scalar current:  J_S
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = J_S(1) * cont_VV(J_V1,J_V2)
end subroutine counter_SVV_S



! ======================================================================
! Vertex counter terms for HEFT.
! ======================================================================

! **********************************************************************
subroutine counter_HG_G(f, S, V2, P2, V_out, P3)
! Higgs-gluon-gluon vertex: general tensor structure
! Jout_G = J_S*(P3*(-f(1)*P2.V2+f(3)*P3.V2)+P2*(-f(2)*P3.V2+f(3)*P2.V2)
!           +V2*(f(4)*P2.P2+f(4)*P3.P3-f(5)*P2.P3))
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), P3(4)
  complex(REALKIND), intent(out) :: V_out(4)
  complex(REALKIND) :: P3V2, P2V2  ,  P2P3, P2P2, P3P3
  P3V2 = cont_VV(P3,V2)
  P2V2 = cont_VV(P2,V2)
!   V_out = (P3*(-f(1)*P2V2+f(3)*P3V2) + P2*(-f(2)*P3V2+f(3)*P2V2) &
!       & +  V2*(f(4)*cont_V(P2)+f(4)*cont_V(P3)-f(5)*cont_VV(P2,P3))) * S(1)

  P2P3 = cont_VV(P2,P3)
  P2P2 = cont_V(P2)
  P3P3 = cont_V(P3)
  V_out = P3*(-f(1)*P2V2+f(3)*P3V2)+P2*(-f(2)*P3V2+f(3)*P2V2)+V2*(f(4)*P2P2+f(4)*P3P3-f(5)*P2P3)
  V_out = S(1)*V_out
end subroutine counter_HG_G


! **********************************************************************
subroutine counter_GG_H(f, V1, P1, V2, P2, S_out)
! Higgs-gluon-gluon vertex with general form factor contributions
! S_out = (f(1)*P1(1)*P2(2)+f(2)*P1(2)*P2(1)+f(3)*(P1(1)*P1(2)+P2(1)*P2(2))
!          +g(1,2)*(f(4)*(P1.P1+P2.P2)+f(5)*P1.P2))*V1(1)*V2(2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P1V1, P2V2, P1V2, P2V1  ,  V1V2, P1P2, P1P1, P2P2
  P1V1 = cont_VV(P1,V1)
  P2V2 = cont_VV(P2,V2)
  P1V2 = cont_VV(P1,V2)
  P2V1 = cont_VV(P2,V1)
!   S_out(1) = f(1)*P1V1*P2V2 + f(2)*P1V2*P2V1 + f(3)*(P1V1*P1V2 + P2V1*P2V2) &
!          & + (f(4)*(cont_V(P1)+cont_V(P2)) + f(5)*cont_VV(P1,P2)) * cont_VV(V1,V2)

  V1V2 = cont_VV(V1,V2)
  P1P2 = cont_VV(P1,P2)
  P1P1 = cont_V(P1)
  P2P2 = cont_V(P2)
  S_out(1) = f(1)*P1V1*P2V2+f(2)*P1V2*P2V1+f(3)*P1V1*P1V2+f(3)*P2V1*P2V2
  S_out(1) = S_out(1) + V1V2*(f(4)*P1P1 + f(4)*P2P2 + f(5)*P1P2)

end subroutine counter_GG_H


! **********************************************************************
subroutine counter_GGG_H(V1, P1, V2, P2, V3, P3, S_out)
! gluon gluon gluon -> Higgs vertex for HEFT
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(P1-P2,V3) + cont_VV(V2,V3)*cont_VV(P2-P3,V1) + cont_VV(V3,V1)*cont_VV(P3-P1,V2)
end subroutine counter_GGG_H


! **********************************************************************
subroutine counter_HGG_G(S, V2, P2, V3, P3, V_out, P4)
! Higgs gluon gluon -> gluon vertex for HEFT
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S*V2(2)*V3(3)
!          = S*V2.V3*(P2-P3)(4) + S*(P3+P4).V2*V3(4) - S*(P4+P2).V3*V2(4)
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S(1) * cont_VV(V2,V3)) * (P2-P3) + (S(1) * cont_VV(P3+P4,V2)) * V3 - (S(1) * cont_VV(P4+P2,V3)) * V2
end subroutine counter_HGG_G


! **********************************************************************
subroutine counter_HGGG_G(S, V2, V3, V4, V_out)
! Effective Higgs gluon gluon gluon -> gluon vertex
! for factorised Lorentz monomials (g(2,3)*g(4,5))*S*V(2)*V(3)*V(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S(1)*cont_VV(V2,V3)) * V4
end subroutine counter_HGGG_G


! **********************************************************************
subroutine counter_GGGG_H(V1, V2, V3, V4, S_out)
! gluon gluon gluon gluon -> Higgs vertex for HEFT
! for factorised Lorentz monomials (g(1,2)*g(3,4))*V(1)*V(2)*V(3)*V(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(V3,V4)
end subroutine counter_GGGG_H


! **********************************************************************
subroutine counter_AQ_H(JA, PA, JQ, PQ, S_out)
! Fermion-scalar-vertex for R2 in HEFT
! Incoming anti-fermion current: JA(4)
! Incoming fermion current:      JQ(4)
! Outgoing scalar current:       S_out = A.slash(PQ-PA).Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: JA(4), PA(4), JQ(4), PQ(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P(4)
  P = PQ-PA
  S_out(1) = (-JA(3)*P(1)-JA(4)*P(3))*JQ(1) + (-JA(4)*P(2)-JA(3)*P(4))*JQ(2) &
         & + (-JA(1)*P(2)+JA(2)*P(3))*JQ(3) + (-JA(2)*P(1)+JA(1)*P(4))*JQ(4)
end subroutine counter_AQ_H


! **********************************************************************
subroutine counter_QH_A(JQ, PQ, S, JA_out, PA)
! Fermion-scalar-vertex for R2 in HEFT
! Incoming fermion current:      JQ(4)
! Incoming scalar current:       S
! Outgoing anti-fermion current: JA_out_i = slash(PQ+PA)_ij * JQ_j * S
! note that PA is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: JQ(4), PQ(4), S(4), PA(4)
  complex(REALKIND), intent(out) :: JA_out(4)
  complex(REALKIND) :: P(4)
  P = PQ+PA
  JA_out(1) = - P(2)*JQ(3) + P(4)*JQ(4)
  JA_out(2) = - P(1)*JQ(4) + P(3)*JQ(3)
  JA_out(3) = - P(1)*JQ(1) - P(4)*JQ(2)
  JA_out(4) = - P(2)*JQ(2) - P(3)*JQ(1)
  JA_out = JA_out*S(1)
end subroutine counter_QH_A


! **********************************************************************
subroutine counter_HA_Q(S, JA, PA, JQ_out, PQ)
! Fermion-scalar-vertex for R2 in HEFT
! Incoming scalar current:       S
! Incoming anti-fermion current: JA(4)
! Outgoing fermion current:      JQ_out_i = S * JA_j * slash(-PQ-PA)_ji
! note that PQ is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S(4), JA(4), PA(4), PQ(4)
  complex(REALKIND), intent(out) :: JQ_out(4)
  complex(REALKIND) :: P(4)
  P = -PQ-PA
  JQ_out(1) = -JA(3)*P(1)-JA(4)*P(3)
  JQ_out(2) = -JA(4)*P(2)-JA(3)*P(4)
  JQ_out(3) = -JA(1)*P(2)+JA(2)*P(3)
  JQ_out(4) = -JA(2)*P(1)+JA(1)*P(4)
  JQ_out = JQ_out*S(1)
end subroutine counter_HA_Q


! **********************************************************************
subroutine counter_HQA_V(S, J_Q, J_A, Jout_V)
! Fermion-scalar-gluon-vertex for R2 in HEFT
! extension of tree-vertex vert_QA_V
! ----------------------------------------------------------------------
! S(4)      = Incoming scalar current
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_V(4) = outgoing gluon current (light-cone rep.)
! Outgoing gluon current:        V_out(4)
! V_out(A) = JA(i) * gamma^A(i,j) * JQ(j) * S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S(4), J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = (Jout_V + Jout_V)*S(1)
end subroutine counter_HQA_V


! **********************************************************************
subroutine counter_VHQ_A(J_V, S, J_Q, Jout_Q)
! Fermion-scalar-gluon-vertex for R2 in HEFT
! extension of tree-vertex vert_VQ_A
! ----------------------------------------------------------------------
! J_V(4)    = incoming gluon current (light-cone rep.)
! S(4)      = Incoming scalar current
! J_Q(4)    = incoming quark current
! Jout_Q(4) = outgoing quark current
! Jout_Q(4) = J_V(1) * gamma_1(4,3) * J_Q(3) * S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), S(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_V(2)*J_Q(3)+J_V(4)*J_Q(4)
  Jout_Q(2) = - J_V(1)*J_Q(4)+J_V(3)*J_Q(3)
  Jout_Q(3) = - J_V(1)*J_Q(1)-J_V(4)*J_Q(2)
  Jout_Q(4) = - J_V(2)*J_Q(2)-J_V(3)*J_Q(1)
  Jout_Q = Jout_Q*S(1)
end subroutine counter_VHQ_A


! **********************************************************************
subroutine counter_AVH_Q(J_A, J_V, S, Jout_A)
! Fermion-scalar-gluon-vertex for R2 in HEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S(4)      = Incoming scalar current
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * gamma_A(j,i) * J_V(A) * S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4), S(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = - J_V(1)*J_A(3) - J_V(3)*J_A(4)
  Jout_A(2) = - J_V(2)*J_A(4) - J_V(4)*J_A(3)
  Jout_A(3) = - J_V(2)*J_A(1) + J_V(3)*J_A(2)
  Jout_A(4) = - J_V(1)*J_A(2) + J_V(4)*J_A(1)
  Jout_A = Jout_A*S(1)
end subroutine counter_AVH_Q


! **********************************************************************
subroutine counter_QAV_H(J_Q, J_A, J_V, S_out)
! Fermion-scalar-gluon-vertex for R2 in HEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S(4)      = Incoming scalar current
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * gamma_A(j,i) * J_V(A) * S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_A(4), J_V(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = (-J_A(3)*J_V(1)-J_A(4)*J_V(3))*J_Q(1) + (-J_A(4)*J_V(2)-J_A(3)*J_V(4))*J_Q(2) &
         & + (-J_A(1)*J_V(2)+J_A(2)*J_V(3))*J_Q(3) + (-J_A(2)*J_V(1)+J_A(1)*J_V(4))*J_Q(4)
end subroutine counter_QAV_H


! ======================================================================
! Vertex counter terms for HHEFT.
! ======================================================================

! **********************************************************************
subroutine counter_HHG_G(f, S1, S2, V2, P2, V_out, P3)
! Higgs-Higgs-gluon-gluon vertex: general tensor structure
! Jout_G = J_S1*J_S2*(P3*(-f(1)*P2.V2+f(3)*P3.V2)+P2*(-f(2)*P3.V2+f(3)*P2.V2)
!           +V2*(f(4)*P2.P2+f(4)*P3.P3-f(5)*P2.P3))
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V2(4), P2(4), P3(4)
  complex(REALKIND), intent(out) :: V_out(4)
  complex(REALKIND) :: P3V2, P2V2
  P3V2 = cont_VV(P3,V2)
  P2V2 = cont_VV(P2,V2)
  V_out = P3*(-f(1)*P2V2+f(3)*P3V2) + P2*(-f(2)*P3V2+f(3)*P2V2) &
      & + V2*(f(4)*cont_V(P2)+f(4)*cont_V(P3)-f(5)*cont_VV(P2,P3))
  V_out = S1(1)*S2(1)*V_out
end subroutine counter_HHG_G


! **********************************************************************
subroutine counter_HGG_H(f, S_in, V1, P1, V2, P2, S_out)
! Higgs-Higgs-gluon-gluon vertex with general form factor contributions
! S_out = (f(1)*P1(1)*P2(2)+f(2)*P1(2)*P2(1)+f(3)*(P1(1)*P1(2)+P2(1)*P2(2))
!          +g(1,2)*(f(4)*(P1.P1+P2.P2)+f(5)*P1.P2))*V1(1)*V2(2)*S_in
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: S_in(4), V1(4), P1(4), V2(4), P2(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P1V1, P2V2, P1V2, P2V1
  P1V1 = cont_VV(P1,V1)
  P2V2 = cont_VV(P2,V2)
  P1V2 = cont_VV(P1,V2)
  P2V1 = cont_VV(P2,V1)
  S_out(1) = f(1)*P1V1*P2V2 + f(2)*P1V2*P2V1 + f(3)*(P1V1*P1V2 + P2V1*P2V2) &
         & + (f(4)*(cont_V(P1)+cont_V(P2)) + f(5)*cont_VV(P1,P2)) * cont_VV(V1,V2)
  S_out(1) = S_out(1)*S_in(1)
end subroutine counter_HGG_H


! **********************************************************************
subroutine counter_HGGG_H(S_in, V1, P1, V2, P2, V3, P3, S_out)
! Higgs gluon gluon gluon -> Higgs vertex for HHEFT
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)*S_in
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S_in(4), V1(4), P1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(P1-P2,V3) + cont_VV(V2,V3)*cont_VV(P2-P3,V1) + cont_VV(V3,V1)*cont_VV(P3-P1,V2)
  S_out(1) = S_out(1)*S_in(1)
end subroutine counter_HGGG_H


! **********************************************************************
subroutine counter_HHGG_G(S1, S2, V2, P2, V3, P3, V_out, P4)
! Higgs gluon gluon -> gluon vertex for HEFT
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S1*S2*V2(2)*V3(3)
!          = S1*S2 * (V2.V3*(P2-P3)(4) + (P3+P4).V2*V3(4) - (P4+P2).V3*V2(4))
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V2(4), P2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = cont_VV(V2,V3)*(P2-P3) + cont_VV(P3+P4,V2)*V3 - cont_VV(P4+P2,V3)*V2
  V_out = S1(1)*S2(1)*V_out
end subroutine counter_HHGG_G


! **********************************************************************
subroutine counter_HHGGG_G(S1, S2, V2, V3, V4, V_out)
! Effective Higgs Higgs gluon gluon gluon -> gluon vertex
! for factorised Lorentz monomials (g(2,3)*g(4,5))*S1*S2*V(2)*V(3)*V(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S1(1)*S2(1)*cont_VV(V2,V3)) * V4
end subroutine counter_HHGGG_G


! **********************************************************************
subroutine counter_HGGGG_H(S1, V1, V2, V3, V4, S_out)
! Higgs gluon gluon gluon gluon -> Higgs vertex for HEFT
! for factorised Lorentz monomials (g(1,2)*g(3,4))*V(1)*V(2)*V(3)*V(4)*S1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), V1(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(V3,V4)*S1(1)
end subroutine counter_HGGGG_H


! **********************************************************************
subroutine counter_HAQ_H(S_in, JA, PA, JQ, PQ, S_out)
! Fermion-scalar-vertex for R2 in HHEFT
! Incoming scalar current:       S_in(4)
! Incoming anti-fermion current: JA(4)
! Incoming fermion current:      JQ(4)
! Outgoing scalar current:       S_out = A.slash(PQ-PA).Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S_in(4), JA(4), PA(4), JQ(4), PQ(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P(4)
  P = PQ-PA
  S_out(1) = (-JA(3)*P(1)-JA(4)*P(3))*JQ(1) + (-JA(4)*P(2)-JA(3)*P(4))*JQ(2) &
         & + (-JA(1)*P(2)+JA(2)*P(3))*JQ(3) + (-JA(2)*P(1)+JA(1)*P(4))*JQ(4)
  S_out(1) = S_out(1)*S_in(1)
end subroutine counter_HAQ_H


! **********************************************************************
subroutine counter_QHH_A(JQ, PQ, S1, S2, JA_out, PA)
! Fermion-scalar-vertex for R2 in HHEFT
! Incoming fermion current:      JQ(4)
! Incoming scalar current:       S1(4), S2(4)
! Outgoing anti-fermion current: JA_out_i = slash(PQ+PA)_ij * JQ_j * S1 * S2
! note that PA is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: JQ(4), PQ(4), S1(4), S2(4), PA(4)
  complex(REALKIND), intent(out) :: JA_out(4)
  complex(REALKIND) :: P(4)
  P = PQ+PA
  JA_out(1) = - P(2)*JQ(3) + P(4)*JQ(4)
  JA_out(2) = - P(1)*JQ(4) + P(3)*JQ(3)
  JA_out(3) = - P(1)*JQ(1) - P(4)*JQ(2)
  JA_out(4) = - P(2)*JQ(2) - P(3)*JQ(1)
  JA_out = JA_out*S1(1)*S2(1)
end subroutine counter_QHH_A


! **********************************************************************
subroutine counter_HHA_Q(S1, S2, JA, PA, JQ_out, PQ)
! Fermion-scalar-scalar-vertex for R2 in HHEFT
! Incoming scalar currents:      S1(4), S2(4)
! Incoming anti-fermion current: JA(4)
! Outgoing fermion current:      JQ_out_i = S * JA_j * slash(-PQ-PA)_ji
! note that PQ is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), JA(4), PA(4), PQ(4)
  complex(REALKIND), intent(out) :: JQ_out(4)
  complex(REALKIND) :: P(4)
  P = -PQ-PA
  JQ_out(1) = -JA(3)*P(1)-JA(4)*P(3)
  JQ_out(2) = -JA(4)*P(2)-JA(3)*P(4)
  JQ_out(3) = -JA(1)*P(2)+JA(2)*P(3)
  JQ_out(4) = -JA(2)*P(1)+JA(1)*P(4)
  JQ_out = JQ_out*S1(1)*S2(1)
end subroutine counter_HHA_Q


! **********************************************************************
subroutine counter_HHQA_V(S1, S2, J_Q, J_A, Jout_V)
! Fermion--scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_QA_V
! ----------------------------------------------------------------------
! S(4)      = Incoming scalar current
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_V(4) = outgoing gluon current (light-cone rep.)
! Outgoing gluon current:        V_out(4)
! V_out(A) = JA(i) * gamma^A(i,j) * JQ(j) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = (Jout_V + Jout_V)*S1(1)*S2(1)
end subroutine counter_HHQA_V


! **********************************************************************
subroutine counter_VHHQ_A(J_V, S1, S2, J_Q, Jout_Q)
! Fermion-scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_VQ_A
! ----------------------------------------------------------------------
! J_V(4)    = incoming gluon current (light-cone rep.)
! S1(4)      = Incoming scalar current
! S2(4)      = Incoming scalar current
! J_Q(4)    = incoming quark current
! Jout_Q(4) = outgoing quark current
! Jout_Q(4) = J_V(1) * gamma_1(4,3) * J_Q(3) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), S1(4), S2(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_V(2)*J_Q(3)+J_V(4)*J_Q(4)
  Jout_Q(2) = - J_V(1)*J_Q(4)+J_V(3)*J_Q(3)
  Jout_Q(3) = - J_V(1)*J_Q(1)-J_V(4)*J_Q(2)
  Jout_Q(4) = - J_V(2)*J_Q(2)-J_V(3)*J_Q(1)
  Jout_Q = Jout_Q*S1(1)*S2(1)
end subroutine counter_VHHQ_A


! **********************************************************************
subroutine counter_AVHH_Q(J_A, J_V, S1, S2, Jout_A)
! Fermion-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S1(4)      = Incoming scalar current
! S2(4)      = Incoming scalar current
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * gamma_A(j,i) * J_V(A) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4), S1(4), S2(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = - J_V(1)*J_A(3) - J_V(3)*J_A(4)
  Jout_A(2) = - J_V(2)*J_A(4) - J_V(4)*J_A(3)
  Jout_A(3) = - J_V(2)*J_A(1) + J_V(3)*J_A(2)
  Jout_A(4) = - J_V(1)*J_A(2) + J_V(4)*J_A(1)
  Jout_A = Jout_A*S1(1)*S2(1)
end subroutine counter_AVHH_Q


! **********************************************************************
subroutine counter_HQAV_H(S_in, J_Q, J_A, J_V, S_out)
! Fermion-scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! S_in(4)   = Incoming scalar current
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S(4)      = Outgoing scalar current
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S_in(4), J_Q(4), J_A(4), J_V(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = (-J_A(3)*J_V(1)-J_A(4)*J_V(3))*J_Q(1) + (-J_A(4)*J_V(2)-J_A(3)*J_V(4))*J_Q(2) &
         & + (-J_A(1)*J_V(2)+J_A(2)*J_V(3))*J_Q(3) + (-J_A(2)*J_V(1)+J_A(1)*J_V(4))*J_Q(4)
  S_out(1) = S_out(1)*S_in(1)
end subroutine counter_HQAV_H


! ======================================================================
! Vertex counter terms for HHHEFT.
! ======================================================================

! **********************************************************************
subroutine counter_HHHG_G(f, S1, S2, S3, V2, P2, V_out, P3)
! Higgs-Higgs-Higgs-gluon-gluon vertex: general tensor structure
! Jout_G = J_S1*J_S2*J_S3*(P3*(-f(1)*P2.V2+f(3)*P3.V2)+P2*(-f(2)*P3.V2+f(3)*P2.V2)
!           +V2*(f(4)*P2.P2+f(4)*P3.P3-f(5)*P2.P3))
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V2(4), P2(4), P3(4)
  complex(REALKIND), intent(out) :: V_out(4)
  complex(REALKIND) :: P3V2, P2V2
  P3V2 = cont_VV(P3,V2)
  P2V2 = cont_VV(P2,V2)
  V_out = P3*(-f(1)*P2V2+f(3)*P3V2) + P2*(-f(2)*P3V2+f(3)*P2V2) &
      & + V2*(f(4)*cont_V(P2)+f(4)*cont_V(P3)-f(5)*cont_VV(P2,P3))
  V_out = S1(1)*S2(1)*S3(1)*V_out
end subroutine counter_HHHG_G


! **********************************************************************
subroutine counter_HHGG_H(f, S1, S2, V1, P1, V2, P2, S_out)
! Higgs-Higgs-gluon-gluon vertex with general form factor contributions
! S_out = (f(1)*P1(1)*P2(2)+f(2)*P1(2)*P2(1)+f(3)*(P1(1)*P1(2)+P2(1)*P2(2))
!          +g(1,2)*(f(4)*(P1.P1+P2.P2)+f(5)*P1.P2))*V1(1)*V2(2)*S_in
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_V
  implicit none
  real(REALKIND), intent(in)     :: f(5)
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V1(4), P1(4), V2(4), P2(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P1V1, P2V2, P1V2, P2V1
  P1V1 = cont_VV(P1,V1)
  P2V2 = cont_VV(P2,V2)
  P1V2 = cont_VV(P1,V2)
  P2V1 = cont_VV(P2,V1)
  S_out(1) = f(1)*P1V1*P2V2 + f(2)*P1V2*P2V1 + f(3)*(P1V1*P1V2 + P2V1*P2V2) &
         & + (f(4)*(cont_V(P1)+cont_V(P2)) + f(5)*cont_VV(P1,P2)) * cont_VV(V1,V2)
  S_out(1) = S_out(1)*S1(1)*S2(1)
end subroutine counter_HHGG_H


! **********************************************************************
subroutine counter_HHGGG_H(S1, S2, V1, P1, V2, P2, V3, P3, S_out)
! Higgs Higgs gluon gluon gluon -> Higgs vertex for HHEFT
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)*S_in
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V1(4), P1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(P1-P2,V3) + cont_VV(V2,V3)*cont_VV(P2-P3,V1) + cont_VV(V3,V1)*cont_VV(P3-P1,V2)
  S_out(1) = S_out(1)*S1(1)*S2(1)
end subroutine counter_HHGGG_H


! **********************************************************************
subroutine counter_HHHGG_G(S1, S2, S3, V2, P2, V3, P3, V_out, P4)
! Higgs Higgs gluon gluon -> gluon vertex for HEFT
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S1*S2*V2(2)*V3(3)
!          = S1*S2 * (V2.V3*(P2-P3)(4) + (P3+P4).V2*V3(4) - (P4+P2).V3*V2(4))
! note: copy from tree vertices
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V2(4), P2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = cont_VV(V2,V3)*(P2-P3) + cont_VV(P3+P4,V2)*V3 - cont_VV(P4+P2,V3)*V2
  V_out = S1(1)*S2(1)*S3(1)*V_out
end subroutine counter_HHHGG_G


! **********************************************************************
subroutine counter_HHHGGG_G(S1, S2, S3, V2, V3, V4, V_out)
! Effective Higgs Higgs gluon gluon gluon -> gluon vertex
! for factorised Lorentz monomials (g(2,3)*g(4,5))*S1*S2*V(2)*V(3)*V(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S1(1)*S2(1)*S3(1)*cont_VV(V2,V3)) * V4
end subroutine counter_HHHGGG_G


! **********************************************************************
subroutine counter_HHGGGG_H(S1, S2, V1, V2, V3, V4, S_out)
! Higgs gluon gluon gluon gluon -> Higgs vertex for HEFT
! for factorised Lorentz monomials (g(1,2)*g(3,4))*V(1)*V(2)*V(3)*V(4)*S1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V1(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(V3,V4)*S1(1)*S2(1)
end subroutine counter_HHGGGG_H


! **********************************************************************
subroutine counter_HHAQ_H(S1, S2, JA, PA, JQ, PQ, S_out)
! Fermion-scalar-vertex for R2 in HHEFT
! Incoming scalar current:       S_in(4)
! Incoming anti-fermion current: JA(4)
! Incoming fermion current:      JQ(4)
! Outgoing scalar current:       S_out = A.slash(PQ-PA).Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), JA(4), PA(4), JQ(4), PQ(4)
  complex(REALKIND), intent(out) :: S_out(4)
  complex(REALKIND) :: P(4)
  P = PQ-PA
  S_out(1) = (-JA(3)*P(1)-JA(4)*P(3))*JQ(1) + (-JA(4)*P(2)-JA(3)*P(4))*JQ(2) &
         & + (-JA(1)*P(2)+JA(2)*P(3))*JQ(3) + (-JA(2)*P(1)+JA(1)*P(4))*JQ(4)
  S_out(1) = S_out(1)*S1(1)*S2(4)
end subroutine counter_HHAQ_H


! **********************************************************************
subroutine counter_QHHH_A(JQ, PQ, S1, S2, S3, JA_out, PA)
! Fermion-scalar-vertex for R2 in HHEFT
! Incoming fermion current:      JQ(4)
! Incoming scalar current:       S1(4), S2(4)
! Outgoing anti-fermion current: JA_out_i = slash(PQ+PA)_ij * JQ_j * S1 * S2
! note that PA is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: JQ(4), PQ(4), S1(4), S2(4), S3(4), PA(4)
  complex(REALKIND), intent(out) :: JA_out(4)
  complex(REALKIND) :: P(4)
  P = PQ+PA
  JA_out(1) = - P(2)*JQ(3) + P(4)*JQ(4)
  JA_out(2) = - P(1)*JQ(4) + P(3)*JQ(3)
  JA_out(3) = - P(1)*JQ(1) - P(4)*JQ(2)
  JA_out(4) = - P(2)*JQ(2) - P(3)*JQ(1)
  JA_out = JA_out*S1(1)*S2(1)*S3(1)
end subroutine counter_QHHH_A


! **********************************************************************
subroutine counter_HHHA_Q(S1, S2, S3, JA, PA, JQ_out, PQ)
! Fermion-scalar-scalar-vertex for R2 in HHEFT
! Incoming scalar currents:      S1(4), S2(4), S3(4)
! Incoming anti-fermion current: JA(4)
! Outgoing fermion current:      JQ_out_i = S * JA_j * slash(-PQ-PA)_ji
! note that PQ is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), JA(4), PA(4), PQ(4)
  complex(REALKIND), intent(out) :: JQ_out(4)
  complex(REALKIND) :: P(4)
  P = -PQ-PA
  JQ_out(1) = -JA(3)*P(1)-JA(4)*P(3)
  JQ_out(2) = -JA(4)*P(2)-JA(3)*P(4)
  JQ_out(3) = -JA(1)*P(2)+JA(2)*P(3)
  JQ_out(4) = -JA(2)*P(1)+JA(1)*P(4)
  JQ_out = JQ_out*S1(1)*S2(1)*S3(1)
end subroutine counter_HHHA_Q


! **********************************************************************
subroutine counter_HHHQA_V(S1, S2, S3, J_Q, J_A, Jout_V)
! Fermion--scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_QA_V
! ----------------------------------------------------------------------
! S(4)      = Incoming scalar current
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_V(4) = outgoing gluon current (light-cone rep.)
! Outgoing gluon current:        V_out(4)
! V_out(A) = JA(i) * gamma^A(i,j) * JQ(j) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = (Jout_V + Jout_V)*S1(1)*S2(1)*S3(1)
end subroutine counter_HHHQA_V


! **********************************************************************
subroutine counter_VHHHQ_A(J_V, S1, S2, S3, J_Q, Jout_Q)
! Fermion-scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_VQ_A
! ----------------------------------------------------------------------
! J_V(4)    = incoming gluon current (light-cone rep.)
! S1(4)      = Incoming scalar current
! S2(4)      = Incoming scalar current
! J_Q(4)    = incoming quark current
! Jout_Q(4) = outgoing quark current
! Jout_Q(4) = J_V(1) * gamma_1(4,3) * J_Q(3) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), S1(4), S2(4), S3(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_V(2)*J_Q(3)+J_V(4)*J_Q(4)
  Jout_Q(2) = - J_V(1)*J_Q(4)+J_V(3)*J_Q(3)
  Jout_Q(3) = - J_V(1)*J_Q(1)-J_V(4)*J_Q(2)
  Jout_Q(4) = - J_V(2)*J_Q(2)-J_V(3)*J_Q(1)
  Jout_Q = Jout_Q*S1(1)*S2(1)*S3(1)
end subroutine counter_VHHHQ_A


! **********************************************************************
subroutine counter_AVHHH_Q(J_A, J_V, S1, S2, S3, Jout_A)
! Fermion-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S1(4)      = Incoming scalar current
! S2(4)      = Incoming scalar current
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * gamma_A(j,i) * J_V(A) * S1 * S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4), S1(4), S2(4), S3(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = - J_V(1)*J_A(3) - J_V(3)*J_A(4)
  Jout_A(2) = - J_V(2)*J_A(4) - J_V(4)*J_A(3)
  Jout_A(3) = - J_V(2)*J_A(1) + J_V(3)*J_A(2)
  Jout_A(4) = - J_V(1)*J_A(2) + J_V(4)*J_A(1)
  Jout_A = Jout_A*S1(1)*S2(1)*S3(1)
end subroutine counter_AVHHH_Q


! **********************************************************************
subroutine counter_HHQAV_H(S1, S2, J_Q, J_A, J_V, S_out)
! Fermion-scalar-scalar-gluon-vertex for R2 in HHEFT
! extension of tree-vertex vert_AV_Q
! ----------------------------------------------------------------------
! S_in(4)   = Incoming scalar current
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! S(4)      = Outgoing scalar current
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), J_Q(4), J_A(4), J_V(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = (-J_A(3)*J_V(1)-J_A(4)*J_V(3))*J_Q(1) + (-J_A(4)*J_V(2)-J_A(3)*J_V(4))*J_Q(2) &
         & + (-J_A(1)*J_V(2)+J_A(2)*J_V(3))*J_Q(3) + (-J_A(2)*J_V(1)+J_A(1)*J_V(4))*J_Q(4)
  S_out(1) = S_out(1)*S1(1)*S2(1)
end subroutine counter_HHQAV_H



! ======================================================================
! Vertex counter terms for HiggsPO.
! ======================================================================

! **********************************************************************
subroutine counter_AZS_Q(g_RL, J_A, J_Z, S, Jout_A)
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! J_A(4)     = incoming anti-quark current
! J_Z(4)     = incoming Z current (light-cone rep.)
! S(1)       = incoming scalar
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_A(4)  = outgoing anti-quark current
! Jout_A(i)  = J_A(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * J_Z(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_A(4), J_Z(4), S(4), Jout_A(4)
  complex(REALKIND) :: g_RL(2)
  Jout_A(1) = g_RL(1) * S(1) * ( - J_Z(1)*J_A(3) - J_Z(3)*J_A(4))
  Jout_A(2) = g_RL(1) * S(1) * ( - J_Z(2)*J_A(4) - J_Z(4)*J_A(3))
  Jout_A(3) = g_RL(2) * S(1) * ( - J_Z(2)*J_A(1) + J_Z(3)*J_A(2))
  Jout_A(4) = g_RL(2) * S(1) * ( - J_Z(1)*J_A(2) + J_Z(4)*J_A(1))
end subroutine counter_AZS_Q


! **********************************************************************
subroutine counter_SQA_Z(g_RL, S, J_Q, J_A, Jout_Z)
! bare SQA -> Z Z-like interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! g_RL(1)   = right-handed coupling gR
! g_RL(2)   = left-handed coupling gL
! Jout_Z(4) = outgoing Z current (light-cone rep.)
! Jout_Z(A) = J_A(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: S(4), J_Q(4), J_A(4), Jout_Z(4)
  complex(REALKIND) :: g_RL(2)
  Jout_Z(1) = - g_RL(2)*J_A(1)*J_Q(3) - g_RL(1)*J_A(4)*J_Q(2)
  Jout_Z(2) = - g_RL(2)*J_A(2)*J_Q(4) - g_RL(1)*J_A(3)*J_Q(1)
  Jout_Z(3) = - g_RL(2)*J_A(1)*J_Q(4) + g_RL(1)*J_A(3)*J_Q(2)
  Jout_Z(4) = - g_RL(2)*J_A(2)*J_Q(3) + g_RL(1)*J_A(4)*J_Q(1)
  Jout_Z = (Jout_Z + Jout_Z)*S(1)
end subroutine counter_SQA_Z


! **********************************************************************
subroutine counter_ZSQ_A(g_RL, J_Z, S, J_Q, Jout_Q)
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! J_Q(4)     = incoming quark current
! S(4)       = incoming scalar
! J_Z(4)     = incoming Z current ("light-cone" rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_Q(4)  = outgoing quark current
! Jout_Q(i)  = J_Z(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: J_Q(4), J_Z(4), S(4), Jout_Q(4)
  complex(REALKIND) :: g_RL(2)
  Jout_Q(1) = g_RL(2) * S(1) * ( - J_Z(2)*J_Q(3) + J_Z(4)*J_Q(4))
  Jout_Q(2) = g_RL(2) * S(1) * ( - J_Z(1)*J_Q(4) + J_Z(3)*J_Q(3))
  Jout_Q(3) = g_RL(1) * S(1) * ( - J_Z(1)*J_Q(1) - J_Z(4)*J_Q(2))
  Jout_Q(4) = g_RL(1) * S(1) * ( - J_Z(2)*J_Q(2) - J_Z(3)*J_Q(1))
end subroutine counter_ZSQ_A


! **********************************************************************
subroutine counter_SQA_V(S, J_Q, J_A, Jout_V)
! QA -> V counter term; without left/right splitting
! Factorised wrt. vert_QA_V: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: S(4), J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = (Jout_V + Jout_V)*S(1)
end subroutine counter_SQA_V


end module ol_counterterms_/**/REALKIND
