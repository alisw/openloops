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


module ol_vertices_/**/REALKIND
  implicit none

  ! interface to support explicit momenta and momentum integers
  interface vert_UV_W
    module procedure vert_UV_W_mexpl, vert_UV_W_mids
  end interface


  contains

! **********************************************************************
subroutine vert_ZQ_A(g_RL, J_Z, J_Q, Jout_Q)
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
  complex(REALKIND), intent(in)  :: g_RL(2), J_Z(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = g_RL(2) * ( - J_Z(2)*J_Q(3) + J_Z(4)*J_Q(4))
  Jout_Q(2) = g_RL(2) * ( - J_Z(1)*J_Q(4) + J_Z(3)*J_Q(3))
  Jout_Q(3) = g_RL(1) * ( - J_Z(1)*J_Q(1) - J_Z(4)*J_Q(2))
  Jout_Q(4) = g_RL(1) * ( - J_Z(2)*J_Q(2) - J_Z(3)*J_Q(1))
end subroutine vert_ZQ_A


! **********************************************************************
subroutine vert_AZ_Q(g_RL, J_A, J_Z, Jout_A)
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
  complex(REALKIND), intent(in)  :: g_RL(2), J_A(4), J_Z(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = g_RL(1) * ( - J_Z(1)*J_A(3) - J_Z(3)*J_A(4))
  Jout_A(2) = g_RL(1) * ( - J_Z(2)*J_A(4) - J_Z(4)*J_A(3))
  Jout_A(3) = g_RL(2) * ( - J_Z(2)*J_A(1) + J_Z(3)*J_A(2))
  Jout_A(4) = g_RL(2) * ( - J_Z(1)*J_A(2) + J_Z(4)*J_A(1))
end subroutine vert_AZ_Q


! **********************************************************************
subroutine vert_QA_Z(g_RL, J_Q, J_A, Jout_Z)
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
  complex(REALKIND), intent(in)  :: g_RL(2), J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_Z(4)
  Jout_Z(1) = - g_RL(2)*J_A(1)*J_Q(3) - g_RL(1)*J_A(4)*J_Q(2)
  Jout_Z(2) = - g_RL(2)*J_A(2)*J_Q(4) - g_RL(1)*J_A(3)*J_Q(1)
  Jout_Z(3) = - g_RL(2)*J_A(1)*J_Q(4) + g_RL(1)*J_A(3)*J_Q(2)
  Jout_Z(4) = - g_RL(2)*J_A(2)*J_Q(3) + g_RL(1)*J_A(4)*J_Q(1)
  Jout_Z = Jout_Z + Jout_Z
end subroutine vert_QA_Z


! **********************************************************************
subroutine vert_WQ_A(J_W, J_Q, Jout_Q)
! bare WQ -> Q W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = incoming quark current
! J_W(4)    = incoming W current ("light-cone" rep.)
! Jout_Q(4) = outgoing quark current
! Jout_Q(i) = J_W(A) * [gamma_A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_W(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_W(2)*J_Q(3) + J_W(4)*J_Q(4)
  Jout_Q(2) = - J_W(1)*J_Q(4) + J_W(3)*J_Q(3)
  Jout_Q(3) = 0
  Jout_Q(4) = 0
end subroutine vert_WQ_A


! **********************************************************************
subroutine vert_AW_Q(J_A, J_W, Jout_A)
! bare AW -> A W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_W(4)    = incoming W current (light-cone rep.)
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * [gamma_A*w_L](j,i) * J_W(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_W(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = 0
  Jout_A(2) = 0
  Jout_A(3) = - J_W(2)*J_A(1) + J_W(3)*J_A(2)
  Jout_A(4) = - J_W(1)*J_A(2) + J_W(4)*J_A(1)
end subroutine vert_AW_Q


! **********************************************************************
subroutine vert_QA_W(J_Q, J_A, Jout_W)
! bare QA -> W W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_W(4) = outgoing W current (light-cone rep.)
! Jout_W(A) = J_A(i) * [gamma^A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_A(4)
  complex(REALKIND), intent(out) :: Jout_W(4)
  Jout_W(1) = - J_A(1)*J_Q(3)
  Jout_W(2) = - J_A(2)*J_Q(4)
  Jout_W(3) = - J_A(1)*J_Q(4)
  Jout_W(4) = - J_A(2)*J_Q(3)
  Jout_W = Jout_W + Jout_W
end subroutine vert_QA_W


! **********************************************************************
subroutine vert_VQ_A(J_V, J_Q, Jout_Q)
! bare VQ -> Q gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = incoming quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! Jout_Q(4) = outgoing quark current
! Jout_Q(i) = J_V(A) * gamma_A(i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), J_Q(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = - J_V(2)*J_Q(3)+J_V(4)*J_Q(4)
  Jout_Q(2) = - J_V(1)*J_Q(4)+J_V(3)*J_Q(3)
  Jout_Q(3) = - J_V(1)*J_Q(1)-J_V(4)*J_Q(2)
  Jout_Q(4) = - J_V(2)*J_Q(2)-J_V(3)*J_Q(1)
end subroutine vert_VQ_A


! **********************************************************************
subroutine vert_AV_Q(J_A, J_V, Jout_A)
! bare AV -> A gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! J_A(4)    = incoming anti-quark current
! J_V(4)    = incoming gluon current (light-cone rep.)
! Jout_A(4) = outgoing anti-quark current
! Jout_A(i) = J_A(j) * gamma_A(j,i) * J_V(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_A(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = - J_V(1)*J_A(3) - J_V(3)*J_A(4)
  Jout_A(2) = - J_V(2)*J_A(4) - J_V(4)*J_A(3)
  Jout_A(3) = - J_V(2)*J_A(1) + J_V(3)*J_A(2)
  Jout_A(4) = - J_V(1)*J_A(2) + J_V(4)*J_A(1)
end subroutine vert_AV_Q


! **********************************************************************
subroutine vert_QA_V(J_Q, J_A, Jout_V)
! bare QA -> V gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_V(4) = outgoing gluon current (light-cone rep.)
! Jout_V(A) = J_A(i) * gamma^A(i,j) * J_Q(j)
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
end subroutine vert_QA_V


! **********************************************************************
! subroutine vert_VV_V(J_V1, P1, J_V2, P2, Jout_V)
! TODO:  <09-11-18, J.-N. Lang> !

subroutine vert_UV_W_mexpl(J_V1, P1, J_V2, P2, Jout_V)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming gluon currents (light-cone rep.)
! Pi(4)      = incoming J_Vi momentum  (light-cone rep.)
! Jout_V(4)  = outgoing gluon current  (light-cone rep.)
! Jout_V(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * J_V1(a1) * J_V2(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  complex(REALKIND), intent(in)  :: J_V1(4), P1(4), J_V2(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  complex(REALKIND) :: J1J2, P1J2, P2J1

  J1J2 = cont_VV(J_V1,J_V2)
  P1J2 = cont_VV(P1+P1+P2,J_V2)
  P2J1 = cont_VV(P1+P2+P2,J_V1)
  Jout_V = J1J2 * (P1 - P2) + P2J1 * J_V2 - P1J2 * J_V1

end subroutine vert_UV_W_mexpl

! **********************************************************************
subroutine vert_UV_W_mids(J_V1, mom1, J_V2, mom2, Jout_V)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming gluon currents (light-cone rep.)
! Pi(4)      = incoming J_Vi momentum  (light-cone rep.)
! Jout_V(4)  = outgoing gluon current  (light-cone rep.)
! Jout_V(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * J_V1(a1) * J_V2(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none

  integer,           intent(in)  :: mom1, mom2
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  complex(REALKIND) :: J1J2, P1J2, P2J1, P1(4), P2(4)

  P1 = get_LC_4(mom1)
  P2 = get_LC_4(mom2)
  J1J2 = cont_VV(J_V1,J_V2)
  P1J2 = cont_VV(P1+P1+P2,J_V2)
  P2J1 = cont_VV(P1+P2+P2,J_V1)
  Jout_V = J1J2 * (P1 - P2) + P2J1 * J_V2 - P1J2 * J_V1

end subroutine vert_UV_W_mids


! **********************************************************************
subroutine vert_EV_V(J_V1, J_V2, J_V3, Jout_V)
! sigma vertex, where the sigma wave function is replaced by two gluon wave functions J_V1 and J_V2
! Jout_V(d) = (g(a,c)*g(b,d) - g(1,4)*g(2,3)) * J_V1(a) * J_V1(b) * J_V1(c)
!           = J_V1.J_V3 * J_V2(d) - J_V2.J_V3 * J_V1(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = cont_VV(J_V1,J_V3) * J_V2 - cont_VV(J_V2,J_V3) * J_V1
end subroutine vert_EV_V


! **********************************************************************
subroutine vert_GGG_G(J_G1, J_G2, J_G3, Jout_G)
! Four-gluon vertex: factorised colour monomials f(a,b,x)*f(c,d,x) (same as EV_V)
! Jout_G(d) = (g(a,c)*g(b,d) - g(1,4)*g(2,3)) * J_G1(a) * J_G1(b) * J_G1(c)
!           = J_G1.J_G3 * J_G2(d) - J_G2.J_G3 * J_G1(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_G1(4), J_G2(4), J_G3(4)
  complex(REALKIND), intent(out) :: Jout_G(4)
  Jout_G = cont_VV(J_G1,J_G3)*J_G2 - cont_VV(J_G2,J_G3)*J_G1
end subroutine vert_GGG_G


! ! **********************************************************************
! subroutine vert_WWG_G(J1, J2, J3, Jout)
! ! Four-gluon vertex: factorised Lorentz monomials g(a,b)*g(c,d)
! ! Jout(d) = g(a,b)*g(c,d) * J1(a) * J2(b) * J3(c) = J1.J2 * J3(d)
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   use ol_contractions_/**/REALKIND, only: cont_VV
!   implicit none
!   complex(REALKIND), intent(in)  :: J1(4), J2(4), J3(4)
!   complex(REALKIND), intent(out) :: Jout(4)
!   Jout = cont_VV(J1,J2) * J3
! end subroutine vert_WWG_G


! **********************************************************************
subroutine vert_WWV_V(J_V1, J_V2, J_V3, Jout_V)
! bare W+ W- A/Z A/Z vertex
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [2*g(a1,a2)*g(a3,a4) - g(a2,a3)*g(a1,a4)
!               - g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  complex(REALKIND) :: J1J2, J1J3, J2J3

  J1J2 = cont_VV(J_V1, J_V2)
  J1J2 = J1J2 + J1J2
  J1J3 = cont_VV(J_V1, J_V3)
  J2J3 = cont_VV(J_V2, J_V3)
  Jout_V = J1J2 * J_V3 - J2J3 * J_V1 - J1J3 * J_V2

end subroutine vert_WWV_V

subroutine vert_VWW_V(J_V1,J_V2,J_V3,Jout_V)
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_V1(4), J_V2(4), J_V3(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  call vert_WWV_V(J_V2,J_V3,J_V1,Jout_V)
end subroutine vert_VWW_V


! **********************************************************************
subroutine vert_SS_S(J_S1, J_S2, Jout_S)
! Three scalar vertex
! Incoming scalar currents: J_S1, J_S2
! Outgoing scalar current:  Jout_S = J_S1 * J_S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = J_S1(1) * J_S2(1)
end subroutine vert_SS_S


! **********************************************************************
subroutine vert_SSS_S(J_S1, J_S2, J_S3, Jout_S)
! Four scalar vertex
! Incoming scalar currents: J_S1, J_S2, J_S3
! Outgoing scalar current:  Jout_S = J_S1 * J_S2 * J_S3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S1(4), J_S2(4), J_S3(4)
  complex(REALKIND), intent(out) :: Jout_S(4)
  Jout_S(1) = J_S1(1) * J_S2(1) * J_S3(1)
end subroutine vert_SSS_S


! **********************************************************************
! subroutine vert_VS_S(J_V, J_S, P1, Jout_S)
subroutine vert_VS_T(J_V, P1, J_S, P2, Jout_S)
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
! end subroutine vert_VS_S
end subroutine vert_VS_T


! **********************************************************************
! subroutine vert_SV_S(J_S, P2, J_V, Jout_S)
subroutine vert_TV_S(J_S, P1, J_V, P2, Jout_S)
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
! end subroutine vert_SV_S
end subroutine vert_TV_S


! **********************************************************************
! subroutine vert_SS_V(J_S1, P1, J_S2, P2, Jout_V)
subroutine vert_ST_V(J_S1, P1, J_S2, P2, Jout_V)
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
end subroutine vert_ST_V


! **********************************************************************
subroutine vert_VV_S(J_V1, J_V2, Jout_S)
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
end subroutine vert_VV_S


! **********************************************************************
subroutine vert_VS_V(J_V, J_S, Jout_V)
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
end subroutine vert_VS_V


! **********************************************************************
subroutine vert_SV_V(J_S, J_V, Jout_V)
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
end subroutine vert_SV_V


! **********************************************************************
subroutine vert_VVS_S(J_V1, J_V2, J_S, Jout_S)
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
end subroutine vert_VVS_S


! **********************************************************************
subroutine vert_SSV_V(J_S1, J_S2, J_V, Jout_V)
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
end subroutine vert_SSV_V


! **********************************************************************
subroutine vert_VSS_V(J_V, J_S1, J_S2, Jout_V)
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
end subroutine vert_VSS_V


! **********************************************************************
subroutine vert_SVV_S(J_S, J_V1, J_V2, Jout_S)
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
end subroutine vert_SVV_S


! **********************************************************************
subroutine vert_AQ_S(g_RL, J_A, J_Q, Jout_S)
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
end subroutine vert_AQ_S


! **********************************************************************
subroutine vert_QS_A(g_RL, J_Q, J_S, Jout_A)
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
end subroutine vert_QS_A


! **********************************************************************
subroutine vert_SA_Q(g_RL, J_S, J_A, Jout_Q)
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
end subroutine vert_SA_Q


! **********************************************************************
subroutine vert_CD_V(J_C, J_D, P2, Jout_V)
! Ghost/anti-ghost/gluon vertex
! Incoming ghost current:      J_C
! Incoming anti-ghost current: J_D, momentum P2
! Outgoing gluon current:      Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_C(4), J_D(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = - J_C(1) * J_D(1) * P2
end subroutine vert_CD_V


! **********************************************************************
subroutine vert_DV_C(J_D, P1, J_V, Jout_D)
! Anti-ghost/gluon/anti-ghost vertex
! Incoming anti-ghost current: J_D, momentum P1
! Incoming gluon current:      J_V(4)
! Outgoing anti-ghost current: Jout_D
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_D(4), P1(4), J_V(4)
  complex(REALKIND), intent(out) :: Jout_D(4)
  Jout_D(1) = - J_D(1) * cont_VV(P1, J_V)
end subroutine vert_DV_C


! **********************************************************************
subroutine vert_VC_D(J_V, P1, J_C, P2, Jout_C)
! Ghost/gluon/ghost vertex
! Incoming gluon current: J_V(4), momentum P1
! Incoming ghost current: J_C,    momentum P2
! Outgoing ghost current: Jout_C
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_V(4), P1(4), J_C(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_C(4)
  Jout_C(1) = J_C(1) * cont_VV(P1+P2, J_V)
end subroutine vert_VC_D


! ======================================================================
! HEFT vertices
! ======================================================================

! **********************************************************************
subroutine vert_GG_H(V1, P1, V2, P2, S_out)
! Effective gluon gluon -> Higgs vertex
! S_out = (P1.P2*g(1,2) - P1(2)*P2(1))*V1(1)*V2(2) = V1.V2*P1.P2 - P1.V2*P2.V1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(P1,P2) - cont_VV(P1,V2)*cont_VV(P2,V1)
end subroutine vert_GG_H


! **********************************************************************
subroutine counter_GG_H_vert(V1, P1, V2, P2, S_out)
! Effective gluon gluon -> Higgs vertex for double operator insertion
! S_out = (P1.P2*g(1,2) - P1(2)*P2(1))*V1(1)*V2(2) = V1.V2*P1.P2 - P1.V2*P2.V1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4)
  complex(REALKIND), intent(out) :: S_out(4)
  call vert_GG_H(V1, P1, V2, P2, S_out)
end subroutine counter_GG_H_vert


! **********************************************************************
subroutine vert_HG_G(S, V2, P2, V_out, P3)
! Effective Higgs gluon -> gluon vertex
! V_out(3) = (P2(3)*P3(2) - P2.P3*g(2,3))*S*V2(2) = S*P3.V2*P2(3) - S*P2.P3*V2(3)
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), P3(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S(1) * cont_VV(P3,V2)) * P2 - (S(1) * cont_VV(P2,P3)) * V2
end subroutine vert_HG_G


! **********************************************************************
subroutine counter_HG_G_vert(S, V2, P2, V_out, P3)
! Effective Higgs gluon -> gluon vertex for double operator insertion
! V_out(3) = (P2(3)*P3(2) - P2.P3*g(2,3))*S*V2(2) = S*P3.V2*P2(3) - S*P2.P3*V2(3)
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), P3(4)
  complex(REALKIND), intent(out) :: V_out(4)
  call vert_HG_G(S, V2, P2, V_out, P3)
end subroutine counter_HG_G_vert


! **********************************************************************
subroutine vert_GGG_H(V1, P1, V2, P2, V3, P3, S_out)
! Effective gluon gluon gluon -> Higgs vertex
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V2)*cont_VV(P1-P2,V3) + cont_VV(V2,V3)*cont_VV(P2-P3,V1) + cont_VV(V3,V1)*cont_VV(P3-P1,V2)
end subroutine vert_GGG_H


! **********************************************************************
subroutine counter_GGG_H_vert(V1, P1, V2, P2, V3, P3, S_out)
! Effective gluon gluon gluon -> Higgs vertex for double operator insertion
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), P1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  call vert_GGG_H(V1, P1, V2, P2, V3, P3, S_out)
end subroutine counter_GGG_H_vert


! **********************************************************************
subroutine vert_HGG_G(S, V2, P2, V3, P3, V_out, P4)
! Effective Higgs gluon gluon -> gluon vertex
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S*V2(2)*V3(3)
!          = S*V2.V3*(P2-P3)(4) + S*(P3+P4).V2*V3(4) - S*(P4+P2).V3*V2(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S(1) * cont_VV(V2,V3)) * (P2-P3) + (S(1) * cont_VV(P3+P4,V2)) * V3 - (S(1) * cont_VV(P4+P2,V3)) * V2
end subroutine vert_HGG_G


! **********************************************************************
subroutine counter_HGG_G_vert(S, V2, P2, V3, P3, V_out, P4)
! Effective Higgs gluon gluon -> gluon vertex for double operator insertion
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S*V2(2)*V3(3)
!          = S*V2.V3*(P2-P3)(4) + S*(P3+P4).V2*V3(4) - S*(P4+P2).V3*V2(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), P2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  call vert_HGG_G(S, V2, P2, V3, P3, V_out, P4)
end subroutine counter_HGG_G_vert


! **********************************************************************
subroutine vert_GGGG_H(V1, V2, V3, V4, S_out)
! Effective gluon gluon gluon gluon -> Higgs vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(1,3)*g(2,4) - g(1,4)*g(2,3)) * V1(1) * V2(2) * V3(3) * V4(4)
!       = V1.V3*V2.V4 - V1.V4*V2.V3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = cont_VV(V1,V3)*cont_VV(V2,V4) - cont_VV(V1,V4)*cont_VV(V2,V3)
end subroutine vert_GGGG_H


! **********************************************************************
subroutine counter_GGGG_H_vert(V1, V2, V3, V4, S_out)
! Effective gluon gluon gluon gluon -> Higgs vertex for double operator insertion
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(1,3)*g(2,4) - g(1,4)*g(2,3)) * V1(1) * V2(2) * V3(3) * V4(4)
!       = V1.V3*V2.V4 - V1.V4*V2.V3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: V1(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  call vert_GGGG_H(V1, V2, V3, V4, S_out)
end subroutine counter_GGGG_H_vert


! **********************************************************************
subroutine vert_HGGG_G(S, V2, V3, V4, V_out)
! Effective Higgs gluon gluon gluon -> gluon vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(5) = (g(2,4)*g(3,5) - g(2,5)*g(3,4)) * S * V2(2) * V3(3) * V4(4)
!          = S*V2.V4*V3(5) - S*V3.V4*V2(5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = (S(1)*cont_VV(V2,V4)) * V3 - (S(1)*cont_VV(V3,V4)) * V2
end subroutine vert_HGGG_G


! **********************************************************************
subroutine counter_HGGG_G_vert(S, V2, V3, V4, V_out)
! Effective Higgs gluon gluon gluon -> gluon vertex for double operator insertion
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(5) = (g(2,4)*g(3,5) - g(2,5)*g(3,4)) * S * V2(2) * V3(3) * V4(4)
!          = S*V2.V4*V3(5) - S*V3.V4*V2(5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S(4), V2(4), V3(4), V4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  call vert_HGGG_G(S, V2, V3, V4, V_out)
end subroutine counter_HGGG_G_vert


! ======================================================================
! HHEFT vertices
! ======================================================================


! **********************************************************************
subroutine vert_HGG_H(S1, V2, P2, V3, P3,  S_out)
! Effective Higgs gluon gluon -> Higgs vertex
! S_out = (P2.P3*g(2,3) - P2(3)*P3(2))*V2(2)*V3(3)*S1 = (V2.V3*P2.P3 - P2.V3*P3.V2)*S1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), V2(4), P2(4), V3(4), P3(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = S1(1)*(cont_VV(V2,V3)*cont_VV(P2,P3) - cont_VV(P2,V3)*cont_VV(P3,V2))
end subroutine vert_HGG_H


! **********************************************************************
subroutine vert_HHG_G(S1, S2, V3, P3, V_out, P4)
! Effective Higgs Higgs gluon -> gluon vertex
! V_out(4) = (P3(4)*P4(3) - P3.P4*g(3,4))*S1*S2*V3(3) = S1*S2*(P3.V2*P3(4) - P3.P4*V3(4))
! note that P4 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*(cont_VV(P4,V3) * P3 - cont_VV(P3,P4) * V3)
end subroutine vert_HHG_G


! **********************************************************************
subroutine counter_HHG_G_vert(S1, S2, V3, P3, V_out, P4)
! Effective Higgs Higgs gluon -> gluon vertex
! V_out(4) = (P3(4)*P4(3) - P3.P4*g(3,4))*S1*S2*V3(3) = S1*S2*(P3.V2*P3(4) - P3.P4*V3(4))
! note that P4 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), P3(4), P4(4)
  complex(REALKIND), intent(out) :: V_out(4)
  call vert_HHG_G(S1, S2, V3, P3, V_out, P4)
end subroutine counter_HHG_G_vert


! **********************************************************************
subroutine vert_HHHG_G(S1, S2, S3, V4, P4, V_out, P5)
! Effective Higgs Higgs Higgs gluon -> gluon vertex
! V_out(5) = (P4(5)*P5(4) - P4.P5*g(4,5))*S1*S2*S3*V4(4) = S1*S2*S3*(P4.V3*P4(5) - P4.P5*V4(5))
! note that P4 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V4(4), P4(4), P5(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*S3(1)*(cont_VV(P5,V4) * P4 - cont_VV(P4,P5) * V4)
end subroutine vert_HHHG_G


! **********************************************************************
subroutine vert_HGGG_H(S1, V2, P2, V3, P3, V4, P4, S_out)
! Effective Higgs gluon gluon gluon -> Higgs vertex
! S_out = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3-P4)(2) + g(4,2)*(P4-P2)(3)) * V2(2)*V3(3)*V4(4)*S1
!       = (V2.V3*(P2-P3).V4 + V3.V4*(P3-P4).V2 + V4.V2*(P4-P2).V3)*S1
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), V2(4), P2(4), V3(4), P3(4), V4(4), P4(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = S1(1)*(cont_VV(V2,V3)*cont_VV(P2-P3,V4) + cont_VV(V3,V4)*cont_VV(P3-P4,V2) + cont_VV(V4,V2)*cont_VV(P4-P2,V3))
end subroutine vert_HGGG_H


! **********************************************************************
subroutine vert_HHGGG_H(S1, S2, V3, P3, V4, P4, V5, P5, S_out)
! Effective Higgs Higgs gluon gluon gluon -> Higgs vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), P3(4), V4(4), P4(4), V5(4), P5(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = S1(1)*S2(1)*(cont_VV(V3,V4)*cont_VV(P3-P4,V5) + cont_VV(V4,V5)*cont_VV(P4-P5,V3) + cont_VV(V5,V3)*cont_VV(P5-P3,V4))
end subroutine vert_HHGGG_H


! **********************************************************************
subroutine vert_HHGG_G(S1, S2, V3, P3, V4, P4, V_out, P5)
! Effective Higgs Higgs gluon gluon -> gluon vertex
! V_out(5) = (g(3,4)*(P3-P4)(5) + g(4,5)*(P4+P5)(3) + g(5,3)*(-P5-P3)(4))*S1*S2*V4(3)*V4(4)
!          = S1*S2*(V3.V4*(P3-P4)(5) + (P4+P5).V3*V4(5) - (P5+P3).V4*V3(5))
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), P3(4), V4(4), P4(4), P5(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*(cont_VV(V3,V4) * (P3-P4) + cont_VV(P4+P5,V3) * V4 - cont_VV(P5+P3,V4) * V3)
end subroutine vert_HHGG_G


! **********************************************************************
subroutine counter_HHGG_G_vert(S1, S2, V3, P3, V4, P4, V_out, P5)
! Effective Higgs Higgs gluon gluon -> gluon vertex
! V_out(5) = (g(3,4)*(P3-P4)(5) + g(4,5)*(P4+P5)(3) + g(5,3)*(-P5-P3)(4))*S1*S2*V4(3)*V4(4)
!          = S1*S2*(V3.V4*(P3-P4)(5) + (P4+P5).V3*V4(5) - (P5+P3).V4*V3(5))
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), P3(4), V4(4), P4(4), P5(4)
  complex(REALKIND), intent(out) :: V_out(4)
  call vert_HHGG_G(S1, S2, V3, P3, V4, P4, V_out, P5)
end subroutine counter_HHGG_G_vert



! **********************************************************************
subroutine vert_HHHGG_G(S1, S2, S3, V4, P4, V5, P5, V_out, P6)
! Effective Higgs Higgs Higgs gluon gluon -> gluon vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V4(4), P4(4), V5(4), P5(4), P6(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*S3(1)*(cont_VV(V4,V5) * (P4-P5) + cont_VV(P5+P6,V4) * V5 - cont_VV(P6+P4,V5) * V4)
end subroutine vert_HHHGG_G


! **********************************************************************
subroutine vert_HGGGG_H(S1, V2, V3, V4, V5, S_out)
! Effective Higgs gluon gluon gluon gluon -> Higgs vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(2,4)*g(3,5) - g(2,5)*g(3,4)) * S1 * V2(2) * V3(3) * V4(4) * V5(5)
!       = S1*(V2.V4*V3.V5 - V2.V5*V3.V4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), V2(4), V3(4), V4(4), V5(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = S1(1)*(cont_VV(V2,V4)*cont_VV(V3,V5) - cont_VV(V2,V5)*cont_VV(V3,V4))
end subroutine vert_HGGGG_H

! **********************************************************************
subroutine vert_HHGGGG_H(S1, S2, V3, V4, V5, V6, S_out)
! Effective Higgs Higgs gluon gluon gluon gluon -> Higgs vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), V4(4), V5(4), V6(4)
  complex(REALKIND), intent(out) :: S_out(4)
  S_out(1) = S1(1)*S2(1)*(cont_VV(V3,V5)*cont_VV(V4,V6) - cont_VV(V3,V6)*cont_VV(V4,V5))
end subroutine vert_HHGGGG_H



! **********************************************************************
subroutine vert_HHGGG_G(S1, S2, V3, V4, V5, V_out)
! Effective Higgs Higgs gluon gluon gluon -> gluon vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(6) = (g(3,5)*g(4,6) - g(3,6)*g(4,5)) * S1 * S2 * V3(3) * V4(4) * V5(5)
!          = S1*S2*(V3.V5*V4(6) - V4.V5*V3(6))
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), V3(4), V4(4), V5(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*(cont_VV(V3,V5) * V4 - cont_VV(V4,V5) * V3)
end subroutine vert_HHGGG_G


! **********************************************************************
subroutine vert_HHHGGG_G(S1, S2, S3, V4, V5, V6, V_out)
! Effective Higgs Higgs Higgs gluon gluon gluon -> gluon vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: S1(4), S2(4), S3(4), V4(4), V5(4), V6(4)
  complex(REALKIND), intent(out) :: V_out(4)
  V_out = S1(1)*S2(1)*S3(1)*(cont_VV(V4,V6) * V5 - cont_VV(V5,V6) * V4)
end subroutine vert_HHHGGG_G



! **********************************************************************
! HiggsPO Effective couplings
! **********************************************************************

! **********************************************************************
subroutine vert_SV_V_kap(J_S, J_V, P, Jout_V, P_out)
! Two vector boson + scalar vertex
! Incoming scalar current: J_S
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing vector current: Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), J_V(4), P(4), P_out(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = J_S(1) * J_V
end subroutine vert_SV_V_kap


! **********************************************************************
subroutine vert_SV_V_eps(J_S, J_V, P, Jout_V, P_out)
! Two vector boson + scalar vertex (
! Incoming scalar current: J_S
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing vector current: Jout_V = -(P.P_out)*V1*S - (P_out.V1)*P*S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), J_V(4), P(4), P_out(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V = -cont_VV(P,P_out)*J_V + cont_VV(P_out,J_V)*P
  Jout_V = J_S(1) * Jout_V
end subroutine vert_SV_V_eps


! **********************************************************************
subroutine vert_SV_V_aeps(J_S, J_V, P1, Jout_V, P2)
! Two vector boson + scalar vertex
! Incoming scalar current: J_S
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing vector current: Jout_V(a) = ep(a,b,c,d)*J_V(b)*P(c)*P2(d)*S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_EpVVV
  implicit none
  complex(REALKIND), intent(in)  :: J_S(4), J_V(4), P1(4), P2(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  call cont_EpVVV(J_V, P1, P2, Jout_V)
  Jout_V = J_S(1) * Jout_V
end subroutine vert_SV_V_aeps


! **********************************************************************
subroutine vert_ZSQ_A(g_RL, J_Z, S, J_Q, Jout_Q)
! bare ZSQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! S(4)       = incoming scalar current
! J_Q(4)     = incoming quark current
! J_Z(4)     = incoming Z current ("light-cone" rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_Q(4)  = outgoing quark current
! Jout_Q(i)  = J_Z(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_Z(4), J_Q(4), S(4)
  complex(REALKIND), intent(out) :: Jout_Q(4)
  Jout_Q(1) = g_RL(2) * ( - J_Z(2)*J_Q(3) + J_Z(4)*J_Q(4))
  Jout_Q(2) = g_RL(2) * ( - J_Z(1)*J_Q(4) + J_Z(3)*J_Q(3))
  Jout_Q(3) = g_RL(1) * ( - J_Z(1)*J_Q(1) - J_Z(4)*J_Q(2))
  Jout_Q(4) = g_RL(1) * ( - J_Z(2)*J_Q(2) - J_Z(3)*J_Q(1))
  Jout_Q = Jout_Q*S(1)
end subroutine vert_ZSQ_A


! **********************************************************************
subroutine vert_SAZ_Q(g_RL, S, J_A, J_Z, Jout_A)
! bare SAZ -> A Z-like interaction
! ----------------------------------------------------------------------
! S(4)       = incoming scalar current
! J_A(4)     = incoming anti-quark current
! J_Z(4)     = incoming Z current (light-cone rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_A(4)  = outgoing anti-quark current
! Jout_A(i)  = J_A(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * J_Z(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_A(4), J_Z(4), S(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = g_RL(1) * ( - J_Z(1)*J_A(3) - J_Z(3)*J_A(4))
  Jout_A(2) = g_RL(1) * ( - J_Z(2)*J_A(4) - J_Z(4)*J_A(3))
  Jout_A(3) = g_RL(2) * ( - J_Z(2)*J_A(1) + J_Z(3)*J_A(2))
  Jout_A(4) = g_RL(2) * ( - J_Z(1)*J_A(2) + J_Z(4)*J_A(1))
  Jout_A = Jout_A*S(1)
end subroutine vert_SAZ_Q


! **********************************************************************
subroutine vert_AZS_Q(g_RL, J_A, J_Z, S, Jout_A)
! bare AZS -> A Z-like interaction
! ----------------------------------------------------------------------
! J_A(4)     = incoming anti-quark current
! J_Z(4)     = incoming Z current (light-cone rep.)
! S(4)       = incoming scalar current
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_A(4)  = outgoing anti-quark current
! Jout_A(i)  = J_A(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * J_Z(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_A(4), J_Z(4), S(4)
  complex(REALKIND), intent(out) :: Jout_A(4)
  Jout_A(1) = g_RL(1) * ( - J_Z(1)*J_A(3) - J_Z(3)*J_A(4))
  Jout_A(2) = g_RL(1) * ( - J_Z(2)*J_A(4) - J_Z(4)*J_A(3))
  Jout_A(3) = g_RL(2) * ( - J_Z(2)*J_A(1) + J_Z(3)*J_A(2))
  Jout_A(4) = g_RL(2) * ( - J_Z(1)*J_A(2) + J_Z(4)*J_A(1))
  Jout_A = Jout_A*S(1)
end subroutine vert_AZS_Q


! **********************************************************************
subroutine vert_SQA_Z(g_RL, S, J_Q, J_A, Jout_Z)
! bare SQA -> Z Z-like interaction
! ----------------------------------------------------------------------
! S(4)       = incoming scalar current
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! g_RL(1)   = right-handed coupling gR
! g_RL(2)   = left-handed coupling gL
! Jout_Z(4) = outgoing Z current (light-cone rep.)
! Jout_Z(A) = J_A(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: g_RL(2), J_Q(4), J_A(4), S(4)
  complex(REALKIND), intent(out) :: Jout_Z(4)
  Jout_Z(1) = - g_RL(2)*J_A(1)*J_Q(3) - g_RL(1)*J_A(4)*J_Q(2)
  Jout_Z(2) = - g_RL(2)*J_A(2)*J_Q(4) - g_RL(1)*J_A(3)*J_Q(1)
  Jout_Z(3) = - g_RL(2)*J_A(1)*J_Q(4) + g_RL(1)*J_A(3)*J_Q(2)
  Jout_Z(4) = - g_RL(2)*J_A(2)*J_Q(3) + g_RL(1)*J_A(4)*J_Q(1)
  Jout_Z = (Jout_Z + Jout_Z)*S(1)
end subroutine vert_SQA_Z


! **********************************************************************
subroutine vert_SQA_V(S, J_Q, J_A, Jout_V)
! bare QA -> V gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! Jout_V(4) = outgoing gluon current (light-cone rep.)
! Jout_V(A) = J_A(i) * gamma^A(i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: J_Q(4), J_A(4), S(4)
  complex(REALKIND), intent(out) :: Jout_V(4)
  Jout_V(1) = - J_A(1)*J_Q(3) - J_A(4)*J_Q(2)
  Jout_V(2) = - J_A(2)*J_Q(4) - J_A(3)*J_Q(1)
  Jout_V(3) = - J_A(1)*J_Q(4) + J_A(3)*J_Q(2)
  Jout_V(4) = - J_A(2)*J_Q(3) + J_A(4)*J_Q(1)
  Jout_V = (Jout_V + Jout_V)*S(1)
end subroutine vert_SQA_V


end module ol_vertices_/**/REALKIND



module ol_s_vertices_/**/REALKIND
  implicit none

  ! interface to support explicit momenta and momentum integers
  interface vert_UV_W
    module procedure vert_UV_W_mexpl, vert_UV_W_mids
  end interface

  contains

!************************************************************************
subroutine vert_ZQ_A(g_RL, Z, Q, Q_out)
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! Q          = incoming quark
! Z          = incoming Z  ("light-cone" rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Q_out      = outgoing quark
! Q_out%j(i) = Z%j(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*Q%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: Q, Z
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: Q_out

  select case (Q%h)

  case (B01)
    Q_out%j(1) = g_RL(2) * ( - Z%j(2)*Q%j(3) + Z%j(4)*Q%j(4))
    Q_out%j(2) = g_RL(2) * ( - Z%j(1)*Q%j(4) + Z%j(3)*Q%j(3))
    Q_out%j(3:4) = 0
    Q_out%h    = B10

  case (B10)
    Q_out%j(1:2) = 0
    Q_out%j(3) = g_RL(1) * ( - Z%j(1)*Q%j(1) - Z%j(4)*Q%j(2))
    Q_out%j(4) = g_RL(1) * ( - Z%j(2)*Q%j(2) - Z%j(3)*Q%j(1))
    Q_out%h    = B01

  case (B00)
    Q_out%j = 0
    Q_out%h = B00

  case default
    Q_out%j(1) = g_RL(2) * ( - Z%j(2)*Q%j(3) + Z%j(4)*Q%j(4))
    Q_out%j(2) = g_RL(2) * ( - Z%j(1)*Q%j(4) + Z%j(3)*Q%j(3))
    Q_out%j(3) = g_RL(1) * ( - Z%j(1)*Q%j(1) - Z%j(4)*Q%j(2))
    Q_out%j(4) = g_RL(1) * ( - Z%j(2)*Q%j(2) - Z%j(3)*Q%j(1))
    Q_out%h    = B11

  end select

end subroutine vert_ZQ_A


! **********************************************************************
subroutine vert_AZ_Q(g_RL, A, Z, A_out)
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! A          = incoming anti-quark
! Z          = incoming Z  (light-cone rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! A_out      = outgoing anti-quark
! A_out%j(i) = A%j(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * Z%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: A, Z
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: A_out

  select case (A%h)

  case (B01)
    A_out%j(1) = g_RL(1) * ( - Z%j(1)*A%j(3) - Z%j(3)*A%j(4))
    A_out%j(2) = g_RL(1) * ( - Z%j(2)*A%j(4) - Z%j(4)*A%j(3))
    A_out%j(3:4) = 0
    A_out%h    = B10

  case (B10)
    A_out%j(1:2) = 0
    A_out%j(3) = g_RL(2) * ( - Z%j(2)*A%j(1) + Z%j(3)*A%j(2))
    A_out%j(4) = g_RL(2) * ( - Z%j(1)*A%j(2) + Z%j(4)*A%j(1))
    A_out%h    = B01

  case (B00)
    A_out%j = 0
    A_out%h = B00

  case default
    A_out%j(1) = g_RL(1) * ( - Z%j(1)*A%j(3) - Z%j(3)*A%j(4))
    A_out%j(2) = g_RL(1) * ( - Z%j(2)*A%j(4) - Z%j(4)*A%j(3))
    A_out%j(3) = g_RL(2) * ( - Z%j(2)*A%j(1) + Z%j(3)*A%j(2))
    A_out%j(4) = g_RL(2) * ( - Z%j(1)*A%j(2) + Z%j(4)*A%j(1))
    A_out%h    = B11

  end select

end subroutine vert_AZ_Q


! **********************************************************************
subroutine vert_QA_Z(g_RL, Q, A, Z_out)
! bare QA -> Z Z-like interaction
! ----------------------------------------------------------------------
! Q         = incoming quark
! A         = incoming anti-quark
! g_RL(1)   = right-handed coupling gR
! g_RL(2)   = left-handed coupling gL
! Z_out%j   = outgoing Z (light-cone rep.)
! Z_out%j(A)= A%j(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * Q%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: Q, A
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: Z_out
  complex(REALKIND) :: A_aux(4)

  select case (ishft(Q%h,2) + A%h)

  case (B1111)
    A_aux(1:2)   = g_RL(2)*A%j(1:2)
    A_aux(3:4)   = g_RL(1)*A%j(3:4)
    Z_out%j(1) = - A_aux(1)*Q%j(3) - A_aux(4)*Q%j(2)
    Z_out%j(2) = - A_aux(2)*Q%j(4) - A_aux(3)*Q%j(1)
    Z_out%j(3) = - A_aux(1)*Q%j(4) + A_aux(3)*Q%j(2)
    Z_out%j(4) = - A_aux(2)*Q%j(3) + A_aux(4)*Q%j(1)
    Z_out%j    =  Z_out%j + Z_out%j

  case (B0110, B0111, B1110)
    A_aux(1:2)   = g_RL(2)*A%j(1:2)
    Z_out%j(1) = - A_aux(1)*Q%j(3)
    Z_out%j(2) = - A_aux(2)*Q%j(4)
    Z_out%j(3) = - A_aux(1)*Q%j(4)
    Z_out%j(4) = - A_aux(2)*Q%j(3)
    Z_out%j    =  Z_out%j + Z_out%j

  case (B1001, B1101, B1011)
    A_aux(3:4)   = g_RL(1)*A%j(3:4)
    Z_out%j(1) = - A_aux(4)*Q%j(2)
    Z_out%j(2) = - A_aux(3)*Q%j(1)
    Z_out%j(3) = + A_aux(3)*Q%j(2)
    Z_out%j(4) = + A_aux(4)*Q%j(1)
    Z_out%j    =  Z_out%j + Z_out%j

  case default
    Z_out%j = 0

  end select

end subroutine vert_QA_Z


! **********************************************************************
subroutine vert_WQ_A(W, Q, Q_out)
! bare WQ -> Q W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! Q         = incoming quark
! W         = incoming W  ("light-cone" rep.)
! Q_out     = outgoing quark
! Q_out%j(i)= W%j(A) * [gamma_A*w_L](i,j) * Q%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: Q, W
  type(wfun), intent(out) :: Q_out

  select case (Q%h)

  case (B01, B11)
  Q_out%j(1) = - W%j(2)*Q%j(3) + W%j(4)*Q%j(4)
  Q_out%j(2) = - W%j(1)*Q%j(4) + W%j(3)*Q%j(3)
  Q_out%j(3:4) = 0
  Q_out%h    = B10

  case default
  Q_out%j = 0
  Q_out%h = B00

  end select

end subroutine vert_WQ_A


! **********************************************************************
subroutine vert_AW_Q(A, W, A_out)
! bare AW -> A W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! A         = incoming anti-quark
! W         = incoming W  (light-cone rep.)
! A_out     = outgoing anti-quark
! A_out%j(i)= A%j(j) * [gamma_A*w_L](j,i) * W%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: A, W
  type(wfun), intent(out) :: A_out

  select case (A%h)

  case (B10, B11)
    A_out%j(1:2) = 0
    A_out%j(3) = - W%j(2)*A%j(1) + W%j(3)*A%j(2)
    A_out%j(4) = - W%j(1)*A%j(2) + W%j(4)*A%j(1)
    A_out%h    = B01

  case default
    A_out%j = 0
    A_out%h = B00

  end select

end subroutine vert_AW_Q


! **********************************************************************
subroutine vert_QA_W(Q, A, W_out)
! bare QA -> W W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! Q         = incoming quark
! A         = incoming anti-quark
! W_out     = outgoing W  (light-cone rep.)
! W_out%j(A)= A%j(i) * [gamma^A*w_L](i,j) * Q%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: Q, A
  type(wfun), intent(out) :: W_out

  select case (ishft(Q%h,2) + A%h)

  case (B1111, B1110, B0111, B0110)
    W_out%j(1) = - A%j(1)*Q%j(3)
    W_out%j(2) = - A%j(2)*Q%j(4)
    W_out%j(3) = - A%j(1)*Q%j(4)
    W_out%j(4) = - A%j(2)*Q%j(3)
    W_out%j    = W_out%j + W_out%j

  case default
    W_out%j = 0

  end select

end subroutine vert_QA_W


! **********************************************************************
subroutine vert_VQ_A(V, Q, Q_out)
! bare VQ -> Q gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! Q          = incoming quark
! V          = incoming gluon (light-cone rep.)
! Q_out      = outgoing quark
! Q_out%j(i) = V%j(A) * gamma_A(i,j) * Q%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: V, Q
  type(wfun), intent(out) :: Q_out

  select case (Q%h)

  case (B01)
    Q_out%j(1) = - V%j(2)*Q%j(3)+V%j(4)*Q%j(4)
    Q_out%j(2) = - V%j(1)*Q%j(4)+V%j(3)*Q%j(3)
    Q_out%j(3:4) = 0
    Q_out%h    = B10

  case (B10)
    Q_out%j(1:2) = 0
    Q_out%j(3) = - V%j(1)*Q%j(1)-V%j(4)*Q%j(2)
    Q_out%j(4) = - V%j(2)*Q%j(2)-V%j(3)*Q%j(1)
    Q_out%h    = B01

  case (B00)
    Q_out%j = 0
    Q_out%h = B00

  case default
    Q_out%j(1) = - V%j(2)*Q%j(3)+V%j(4)*Q%j(4)
    Q_out%j(2) = - V%j(1)*Q%j(4)+V%j(3)*Q%j(3)
    Q_out%j(3) = - V%j(1)*Q%j(1)-V%j(4)*Q%j(2)
    Q_out%j(4) = - V%j(2)*Q%j(2)-V%j(3)*Q%j(1)
    Q_out%h    = B11

  end select

end subroutine vert_VQ_A


! **********************************************************************
subroutine vert_AV_Q(A, V, A_out)
! bare AV -> A gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! A          = incoming anti-quark
! V          = incoming gluon (light-cone rep.)
! A_out      = outgoing anti-quark wfun
! A_out%j(i) = A%j(j) * gamma_A(j,i) * V%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: A, V
  type(wfun), intent(out) :: A_out

  select case (A%h)

  case (B01)
    A_out%j(1) = - V%j(1)*A%j(3) - V%j(3)*A%j(4)
    A_out%j(2) = - V%j(2)*A%j(4) - V%j(4)*A%j(3)
    A_out%j(3:4) = 0
    A_out%h    = B10

  case (B10)
    A_out%j(1:2) = 0
    A_out%j(3) = - V%j(2)*A%j(1) + V%j(3)*A%j(2)
    A_out%j(4) = - V%j(1)*A%j(2) + V%j(4)*A%j(1)
    A_out%h    = B01

  case (B00)
    A_out%j = 0
    A_out%h = B00

  case default
    A_out%j(1) = - V%j(1)*A%j(3) - V%j(3)*A%j(4)
    A_out%j(2) = - V%j(2)*A%j(4) - V%j(4)*A%j(3)
    A_out%j(3) = - V%j(2)*A%j(1) + V%j(3)*A%j(2)
    A_out%j(4) = - V%j(1)*A%j(2) + V%j(4)*A%j(1)
    A_out%h    = B11

  end select

end subroutine vert_AV_Q


! **********************************************************************
subroutine vert_QA_V(Q, A, V_out)
! bare QA -> V gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! Q       = incoming quark
! A       = incoming anti-quark
! V_out   = outgoing gluon
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: Q, A
  type(wfun), intent(out) :: V_out

  select case (ishft(Q%h,2) + A%h)

  case (B1111)
    V_out%j(1) = - A%j(1)*Q%j(3) - A%j(4)*Q%j(2)
    V_out%j(2) = - A%j(2)*Q%j(4) - A%j(3)*Q%j(1)
    V_out%j(3) = - A%j(1)*Q%j(4) + A%j(3)*Q%j(2)
    V_out%j(4) = - A%j(2)*Q%j(3) + A%j(4)*Q%j(1)
    V_out%j = V_out%j + V_out%j

  case (B0110, B0111, B1110)
    V_out%j(1) = - A%j(1)*Q%j(3)
    V_out%j(2) = - A%j(2)*Q%j(4)
    V_out%j(3) = - A%j(1)*Q%j(4)
    V_out%j(4) = - A%j(2)*Q%j(3)
    V_out%j = V_out%j + V_out%j

  case (B1001, B1101, B1011)
    V_out%j(1) = - A%j(4)*Q%j(2)
    V_out%j(2) = - A%j(3)*Q%j(1)
    V_out%j(3) =   A%j(3)*Q%j(2)
    V_out%j(4) =   A%j(4)*Q%j(1)
    V_out%j = V_out%j + V_out%j

  case default
    V_out%j = 0

  end select

end subroutine vert_QA_V


! **********************************************************************
subroutine vert_UV_W_mexpl(V1, P1, V2, P2, V_out)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! Vi          = incoming gluon (light-cone rep.)
! Pi(4)       = incoming Vi momentum  (light-cone rep.)
! V_out(4)    = outgoing gluon (light-cone rep.)
! V_out%j(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * V1%j(a1) * V2%j(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun),        intent(in)  :: V1, V2
  complex(REALKIND), intent(in)  :: P1(4), P2(4)
  type(wfun),        intent(out) :: V_out
  complex(REALKIND) :: J1J2, P1J2, P2J1

  J1J2 = cont_PP(V1%j,V2%j)
  P1J2 = cont_PP(P1+P1+P2,V2%j)
  P2J1 = cont_PP(P1+P2+P2,V1%j)
  V_out%j = J1J2 * (P1 - P2) + P2J1 * V2%j - P1J2 * V1%j
end subroutine vert_UV_W_mexpl


! **********************************************************************
subroutine vert_UV_W_mids(V1, mom1, V2, mom2, V_out)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! Vi          = incoming gluon (light-cone rep.)
! Pi(4)       = incoming Vi momentum  (light-cone rep.)
! V_out(4)    = outgoing gluon (light-cone rep.)
! V_out%j(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * V1%j(a1) * V2%j(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  type(wfun),        intent(in)  :: V1, V2
  integer,           intent(in)  :: mom1, mom2
  type(wfun),        intent(out) :: V_out
  complex(REALKIND) :: J1J2, P1J2, P2J1, P1(4), P2(4)

  P1 = get_LC_4(mom1)
  P2 = get_LC_4(mom2)
  J1J2 = cont_PP(V1%j,V2%j)
  P1J2 = cont_PP(P1+P1+P2,V2%j)
  P2J1 = cont_PP(P1+P2+P2,V1%j)
  V_out%j = J1J2 * (P1 - P2) + P2J1 * V2%j - P1J2 * V1%j
end subroutine vert_UV_W_mids


! **********************************************************************
subroutine vert_EV_V(V1, V2, V3, V_out)
! sigma vertex, where the sigma wave function is replaced
! by two gluon wave functions V1%j and V2%j
! V_out%j(d) = (g(a,c)*g(b,d) + g(1,4)*g(2,3)) * V1%j(a) * V1%j(b) * V1%j(c)
!            = V1%j.V3%j * V2%j(d) + V2%j.V3%j * V1%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: V1, V2, V3
  type(wfun), intent(out) :: V_out
  V_out%j = cont_PP(V1%j,V3%j) * V2%j - cont_PP(V2%j,V3%j) * V1%j
end subroutine vert_EV_V


! **********************************************************************
subroutine vert_GGG_G(G1, G2, G3, G_out)
! Four-gluon vertex: factorised colour monomials f(a,b,x)*f(c,d,x) (same as EV_V)
! G_out%j(d) = (g(a,c)*g(b,d) + g(1,4)*g(2,3)) * G1%j(a) * G1%j(b) * G1%j(c)
!           = G1%j.G3%j * G2%j(d) + G2%j.G3%j * G1%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)   :: G1, G2, G3
  type(wfun), intent(out)  :: G_out
  G_out%j = cont_PP(G1%j,G3%j)*G2%j - cont_PP(G2%j,G3%j)*G1%j
end subroutine vert_GGG_G


! **********************************************************************
subroutine vert_WWG_G(V1, V2, V3, V_out)
! Four-gluon vertex: factorised Lorentz monomials g(a,b)*g(c,d)
! V_out%j(d) = g(a,b)*g(c,d) * V1%j(a) * V2%j(b) * V3%j(c)
!            = V1%j.V2%j * V3%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: V1, V2, V3
  type(wfun), intent(out) :: V_out
  V_out%j = cont_PP(V1%j,V2%j) * V3%j
end subroutine vert_WWG_G


! **********************************************************************
subroutine vert_WWV_V(V1, V2, V3, V_out)
! bare W+ W- A/Z A/Z vertex
! ----------------------------------------------------------------------
! Vi          = incoming vector bosons (light-cone rep.)
! V_out       = outgoing vector boson (light-cone rep.)
! V_out%j(a4) = [2*g(a1,a2)*g(a3,a4) - g(a2,a3)*g(a1,a4)
!               - g(a1,a3)*g(a2,a4)] * V1%j(a1) * V2%j(a2) * V3%j(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: V1, V2, V3
  type(wfun), intent(out) :: V_out
  complex(REALKIND) :: J1J2, J1J3, J2J3

  J1J2 = cont_PP(V1%j, V2%j)
  J1J2 = J1J2 + J1J2
  J1J3 = cont_PP(V1%j, V3%j)
  J2J3 = cont_PP(V2%j, V3%j)

  V_out%j = J1J2 * V3%j - J2J3 * V1%j - J1J3 * V2%j

end subroutine vert_WWV_V


subroutine vert_VWW_V(V1,V2,V3,V_out)
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: V1, V2, V3
  type(wfun), intent(out) :: V_out
  call vert_WWV_V(V2,V3,V1,V_out)
end subroutine vert_VWW_V


! **********************************************************************
subroutine vert_SS_S(J_S1, J_S2, Jout_S)
! Three scalar vertex
! Incoming scalars : J_S1, J_S2
! Outgoing scalar  : Jout_S = J_S1 * J_S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: J_S1, J_S2
  type(wfun), intent(out) :: Jout_S
  Jout_S%j(1) = J_S1%j(1) * J_S2%j(1)
end subroutine vert_SS_S


! **********************************************************************
subroutine vert_SSS_S(J_S1, J_S2, J_S3, Jout_S)
! Four scalar vertex
! Incoming scalars : J_S1, J_S2, J_S3
! Outgoing scalar  : Jout_S = J_S1 * J_S2 * J_S3
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: J_S1, J_S2, J_S3
  type(wfun), intent(out) :: Jout_S
  Jout_S%j(1) = J_S1%j(1) * J_S2%j(1) * J_S3%j(1)
end subroutine vert_SSS_S


! **********************************************************************
! subroutine vert_VS_S(V, J_S, P1, Jout_S)
subroutine vert_VS_T(V, P1, J_S, P2, Jout_S)
! Vector boson + two scalars vertex
! Incoming vector : V,      incoming momentum P1(4) (light-cone rep.)
! Incoming scalar : J_S,    incoming momentum P2(4) (light-cone rep.)
! Outgoing scalar : Jout_S = V%j.(2*P2+P1) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun),        intent(in)  :: V, J_S
  complex(REALKIND), intent(in)  :: P1(4), P2(4)
  type(wfun),        intent(out) :: Jout_S
  Jout_S%j(1) = cont_PP(P1+P2+P2, V%j) * J_S%j(1)
! end subroutine vert_VS_S
end subroutine vert_VS_T


! **********************************************************************
! subroutine vert_SV_S(J_S, P2, V, Jout_S)
subroutine vert_TV_S(J_S, P1, V, P2, Jout_S)
! Vector boson + two scalars vertex
! Incoming scalar : J_S, incoming momentum P2(4) (light-cone rep.)
! Incoming vector : V    (light-cone rep.)
! Outgoing scalar : Jout_S = V%j.(-2*P1-P2) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun),        intent(in)  :: V, J_S
  complex(REALKIND), intent(in)  :: P1(4), P2(4)
  type(wfun),        intent(out) :: Jout_S
  Jout_S%j(1) = - cont_PP(P1+P1+P2, V%j) * J_S%j(1)
! end subroutine vert_SV_S
end subroutine vert_TV_S


! **********************************************************************
! subroutine vert_SS_V(J_S1, P1, J_S2, P2, V_out)
subroutine vert_ST_V(J_S1, P1, J_S2, P2, V_out)
! Vector boson + two scalars vertex
! Incoming scalars: J_S1, J_S2, incoming momenta P1(4), P2(4) (light-cone rep.)
! Outgoing vector : V_out%j = J_S1 * J_S2 * (P1 - P2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: J_S1, J_S2
  complex(REALKIND), intent(in)  :: P1(4), P2(4)
  type(wfun),        intent(out) :: V_out
  V_out%j = (J_S1%j(1) * J_S2%j(1)) * (P1 - P2)
! end subroutine vert_SS_V
end subroutine vert_ST_V


! **********************************************************************
subroutine vert_VV_S(V1, V2, Jout_S)
! Two vector boson + scalar vertex
! Incoming vectors: V1, V2  (light-cone rep.)
! Outgoing scalar : Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: V1, V2
  type(wfun), intent(out) :: Jout_S
  Jout_S%j(1) = cont_PP(V1%j,V2%j)
end subroutine vert_VV_S


! **********************************************************************
subroutine vert_VS_V(V, J_S, V_out)
! Two vector boson + scalar vertex
! Incoming vector : V (light-cone rep.)
! Incoming scalar : J_S
! Outgoing vector : V_out
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: V, J_S
  type(wfun), intent(out) :: V_out
  V_out%j = V%j * J_S%j(1)
end subroutine vert_VS_V


! **********************************************************************
subroutine vert_SV_V(J_S, V, V_out)
! Two vector boson + scalar vertex
! Incoming scalar : J_S
! Incoming vector : V (light-cone rep.)
! Outgoing vector : V_out
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: J_S, V
  type(wfun), intent(out) :: V_out
  V_out%j = J_S%j(1) * V%j
end subroutine vert_SV_V


! **********************************************************************
subroutine vert_VVS_S(V1, V2, J_S, Jout_S)
! Two vector boson + two scalars vertex
! Incoming vectors:  V1, V2 (light-cone rep.)
! Incoming scalar :  J_S
! Outgoing scalar :  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: V1, V2, J_S
  type(wfun), intent(out) :: Jout_S
  Jout_S%j(1) = cont_PP(V1%j,V2%j) * J_S%j(1)
end subroutine vert_VVS_S


! **********************************************************************
subroutine vert_SSV_V(J_S1, J_S2, V, V_out)
! Two vector boson + two scalars vertex
! Incoming scalars: J_S1, J_S2
! Incoming vector : V (light-cone rep.)
! Outgoing vector : V_out (light-cone rep.)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: J_S1, J_S2, V
  type(wfun), intent(out) :: V_out
  V_out%j = (J_S1%j(1) * J_S2%j(1)) * V%j
end subroutine vert_SSV_V


! **********************************************************************
subroutine vert_VSS_V(V, J_S1, J_S2, V_out)
! Two vector boson + two scalars vertex
! Incoming vector :  V (light-cone rep.)
! Incoming scalar s: J_S1, J_S2
! Outgoing vector :  V_out (light-cone rep.)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(in)  :: V, J_S1, J_S2
  type(wfun), intent(out) :: V_out
  V_out%j = (J_S1%j(1) * J_S2%j(1)) * V%j
end subroutine vert_VSS_V


! **********************************************************************
subroutine vert_SVV_S(J_S, V1, V2, Jout_S)
! Two vector boson + two scalars vertex
! Incoming scalar :  J_S
! Incoming vectors:  V1, V2 (light-cone rep.)
! Outgoing scalar :  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun), intent(in)  :: J_S, V1, V2
  type(wfun), intent(out) :: Jout_S
  Jout_S%j(1) = J_S%j(1) * cont_PP(V1%j,V2%j)
end subroutine vert_SVV_S


! **********************************************************************
subroutine vert_AQ_S(g_RL, A, Q, Jout_S)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming anti-fermion : A
! Incoming fermion      : Q
! Outgoing scalar       : Jout_S = gR*A%j.P_R.Q%j + gL*A%j.P_L.Q%j
!   with the right- and left-handed projectors P_R = (1+y5)/2 and P_L = (1-y5)/2
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: A, Q
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: Jout_S

  select case (ishft(Q%h,2) + A%h)
  case (B1111)
    Jout_S%j(1) = g_RL(1) * (A%j(1)*Q%j(1) + A%j(2)*Q%j(2)) + g_RL(2) * (A%j(3)*Q%j(3) + A%j(4)*Q%j(4))
  case (B1010, B1110, B1011)
    Jout_S%j(1) = g_RL(1) * (A%j(1)*Q%j(1) + A%j(2)*Q%j(2))
  case (B0101, B1101, B0111)
    Jout_S%j(1) = g_RL(2) * (A%j(3)*Q%j(3) + A%j(4)*Q%j(4))
  case default
    Jout_S%j(1) = 0
  end select

end subroutine vert_AQ_S


! **********************************************************************
subroutine vert_QS_A(g_RL, Q, J_S, A_out)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming fermion      :  Q
! Incoming scalar       :  J_S
! Outgoing anti-fermion :  A_out
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: Q, J_S
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: A_out
  complex(REALKIND) :: g_aux(2)

  select case (Q%h)

  case (B01)
    g_aux(2)   = g_RL(2) * J_S%j(1)
    A_out%j(1:2) = 0
    A_out%j(3) = g_aux(2) * Q%j(3)
    A_out%j(4) = g_aux(2) * Q%j(4)
    A_out%h    = B01

  case (B10)
    g_aux(1)   = g_RL(1) * J_S%j(1)
    A_out%j(1) = g_aux(1) * Q%j(1)
    A_out%j(2) = g_aux(1) * Q%j(2)
    A_out%j(3:4) = 0
    A_out%h    = B10

  case (B00)
    A_out%j = 0
    A_out%h = B00

  case default
    g_aux        = g_RL * J_S%j(1)
    A_out%j(1:2) = g_aux(1) * Q%j(1:2)
    A_out%j(3:4) = g_aux(2) * Q%j(3:4)
    A_out%h      = B11

  end select

end subroutine vert_QS_A


! **********************************************************************
subroutine vert_SA_Q(g_RL, J_S, A, Q_out)
! Fermion-scalar-vertex
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming scalar       : J_S
! Incoming anti-fermion : A
! Outgoing fermion      : Q_out
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: J_S, A
  complex(REALKIND), intent(in)  :: g_RL(2)
  type(wfun),        intent(out) :: Q_out
  complex(REALKIND) :: g_aux(2)

  select case (A%h)

  case (B01)
    g_aux(2)   = g_RL(2) * J_S%j(1)
    Q_out%j(1:2) = 0
    Q_out%j(3) = g_aux(2) * A%j(3)
    Q_out%j(4) = g_aux(2) * A%j(4)
    Q_out%h    = B01

  case (B10)
    g_aux(1)   = g_RL(1) * J_S%j(1)
    Q_out%j(1) = g_aux(1) * A%j(1)
    Q_out%j(2) = g_aux(1) * A%j(2)
    Q_out%j(3:4) = 0
    Q_out%h    = B10

  case (B00)
    Q_out%j = 0
    Q_out%h = B00

  case default
    g_aux = g_RL * J_S%j(1)
    Q_out%j(1:2) = g_aux(1) * A%j(1:2)
    Q_out%j(3:4) = g_aux(2) * A%j(3:4)
    Q_out%h      = B11

  end select

end subroutine vert_SA_Q


! **********************************************************************
subroutine vert_CD_V(J_C, J_D, P2, V_out)
! Ghost/anti-ghost/gluon vertex
! Incoming ghost      : J_C
! Incoming anti-ghost : J_D, momentum P2
! Outgoing gluon      : V_out
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun),        intent(in)  :: J_C, J_D
  complex(REALKIND), intent(in)  :: P2(4)
  type(wfun),        intent(out) :: V_out
  V_out%j = - J_C%j(1) * J_D%j(1) * P2
end subroutine vert_CD_V


! **********************************************************************
subroutine vert_DV_C(J_D, P1, V, Jout_D)
! Anti-ghost/gluon/anti-ghost vertex
! Incoming anti-ghost : J_D, momentum P1
! Incoming gluon      : V
! Outgoing anti-ghost : Jout_D
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun),        intent(in)  :: J_D, V
  complex(REALKIND), intent(in)  :: P1(4)
  type(wfun),        intent(out) :: Jout_D
  Jout_D%j(1) = - J_D%j(1) * cont_PP(P1, V%j)
end subroutine vert_DV_C


! **********************************************************************
subroutine vert_VC_D(V, P1, J_C, P2, Jout_C)
! Ghost/gluon/ghost vertex
! Incoming gluon : V,      momentum P1
! Incoming ghost : J_C,    momentum P2
! Outgoing ghost : Jout_C
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_s_contractions_/**/REALKIND, only: cont_PP
  implicit none
  type(wfun),        intent(in)  :: V, J_C
  complex(REALKIND), intent(in)  :: P1(4), P2(4)
  type(wfun),        intent(out) :: Jout_C
  Jout_C%j(1) = J_C%j(1) * cont_PP(P1+P2, V%j)
end subroutine vert_VC_D

end module ol_s_vertices_/**/REALKIND



module ol_h_vertices_/**/REALKIND
  implicit none

  ! interface to support explicit momenta and momentum integers
  interface vert_UV_W
    module procedure vert_UV_W_mexpl, vert_UV_W_mids
  end interface

  interface vert_VS_T
    module procedure vert_VS_T_mexpl, vert_VS_T_mids
  end interface

  interface vert_TV_S
    module procedure vert_TV_S_mexpl, vert_TV_S_mids
  end interface

  interface vert_ST_V
    module procedure vert_ST_V_mexpl, vert_ST_V_mids
  end interface

  interface vert_DV_C
    module procedure vert_DV_C_mexpl, vert_DV_C_mids
  end interface

  interface vert_VC_D
    module procedure vert_VC_D_mexpl, vert_VC_D_mids
  end interface

  contains

! ! **********************************************************************
! subroutine standard_vert_UV_W(ntry, V1, P1, V2, P2, V_out, n, t)
! ! bare VV -> V vertex
! ! ----------------------------------------------------------------------
! ! ntry           = 1 (2) for 1st (subsequent) PS points
! ! Vi(1:n(i))     = incoming gluon (light-cone representation)
! ! Pi(4)          = incoming Vi momentum  (light-cone representation)
! ! V_out(1:n(3))  = outgoing gluon (light-cone representation)
! ! V_out(h)%j(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
! !                + g(a3,a1)*[-Pout-P1](a2)} * V1(t(1,h))%j(a1) * V2(t(2,h))%j(a2)
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND
!   use ol_data_types_/**/REALKIND, only: wfun
!   use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
!   use ol_h_contractions_/**/REALKIND, only: cont_PP
!   implicit none
!   integer(intkind1), intent(in)    :: ntry
!   integer(intkind2), intent(inout) :: n(3), t(2,n(3))
!   type(wfun),        intent(in)    :: V1(n(1)), V2(n(2))
!   complex(REALKIND),       intent(in)    :: P1(4), P2(4)
!   type(wfun),        intent(out)   :: V_out(n(3))
!   complex(REALKIND) :: J1J2, P1J2, P2J1
!   complex(REALKIND) :: P12(4), P112(4), P122(4)
!   integer :: h
!
!   P12  = P1 - P2
!   P112 = P1 + P1 + P2
!   P122 = P1 + P2 + P2
!   do h = 1, n(3)
!     J1J2 = cont_PP(V1(t(1,h))%j,V2(t(2,h))%j)
!     P1J2 = cont_PP(P112,V2(t(2,h))%j)
!     P2J1 = cont_PP(P122,V1(t(1,h))%j)
!     V_out(h)%j = J1J2 * P12 + P2J1 * V2(t(2,h))%j - P1J2 * V1(t(1,h))%j
!   end do
!
!   if (ntry == 1) call helbookkeeping_vert3(ntry, V1, V2, V_out, n, t)
!
! end subroutine standard_vert_UV_W

!************************************************************************
subroutine vert_ZQ_A(g_RL, ntry, Z, Q, Q_out, n, t)
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! Z(1:n(2)      = incoming Z (light-cone representation)
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling gL
! Q_out(1:n(3)) = outgoing quark
! Q_out(h)%j(i) = Z(t(1,h))%j(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*Q(t(2,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: Z(n(1)), Q(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: Q_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (Q(t(2,h))%h)

    case (B01)
      Q_out(h)%j(1) = g_RL(2) * ( - Z(t(1,h))%j(2)*Q(t(2,h))%j(3) + Z(t(1,h))%j(4)*Q(t(2,h))%j(4))
      Q_out(h)%j(2) = g_RL(2) * ( - Z(t(1,h))%j(1)*Q(t(2,h))%j(4) + Z(t(1,h))%j(3)*Q(t(2,h))%j(3))
      Q_out(h)%j(3:4) = 0
      Q_out(h)%h    = B10

    case (B10)
      Q_out(h)%j(1:2) = 0
      Q_out(h)%j(3) = g_RL(1) * ( - Z(t(1,h))%j(1)*Q(t(2,h))%j(1) - Z(t(1,h))%j(4)*Q(t(2,h))%j(2))
      Q_out(h)%j(4) = g_RL(1) * ( - Z(t(1,h))%j(2)*Q(t(2,h))%j(2) - Z(t(1,h))%j(3)*Q(t(2,h))%j(1))
      Q_out(h)%h    = B01

    case (B00)
      Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      Q_out(h)%h = B00

    case default
      Q_out(h)%j(1) = g_RL(2) * ( - Z(t(1,h))%j(2)*Q(t(2,h))%j(3) + Z(t(1,h))%j(4)*Q(t(2,h))%j(4))
      Q_out(h)%j(2) = g_RL(2) * ( - Z(t(1,h))%j(1)*Q(t(2,h))%j(4) + Z(t(1,h))%j(3)*Q(t(2,h))%j(3))
      Q_out(h)%j(3) = g_RL(1) * ( - Z(t(1,h))%j(1)*Q(t(2,h))%j(1) - Z(t(1,h))%j(4)*Q(t(2,h))%j(2))
      Q_out(h)%j(4) = g_RL(1) * ( - Z(t(1,h))%j(2)*Q(t(2,h))%j(2) - Z(t(1,h))%j(3)*Q(t(2,h))%j(1))
      Q_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, Z, Q, Q_out, n, t)

end subroutine vert_ZQ_A


! **********************************************************************
subroutine vert_AZ_Q(g_RL, ntry, A, Z, A_out, n, t)
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! A(1:n(1))     = incoming anti-quark
! Z(1:n(2))     = incoming Z (light-cone representation)
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling gL
! A_out(1:n(3)) = outgoing anti-quark
! A_out(h)%j(i) = A(t(1,h))%j(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * Z(t(2,h))%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: A(n(1)), Z(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: A_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (A(t(1,h))%h)

    case (B01)
      A_out(h)%j(1) = g_RL(1) * ( - Z(t(2,h))%j(1)*A(t(1,h))%j(3) - Z(t(2,h))%j(3)*A(t(1,h))%j(4))
      A_out(h)%j(2) = g_RL(1) * ( - Z(t(2,h))%j(2)*A(t(1,h))%j(4) - Z(t(2,h))%j(4)*A(t(1,h))%j(3))
      A_out(h)%j(3:4) = 0
      A_out(h)%h    = B10

    case (B10)
      A_out(h)%j(1:2) = 0
      A_out(h)%j(3) = g_RL(2) * ( - Z(t(2,h))%j(2)*A(t(1,h))%j(1) + Z(t(2,h))%j(3)*A(t(1,h))%j(2))
      A_out(h)%j(4) = g_RL(2) * ( - Z(t(2,h))%j(1)*A(t(1,h))%j(2) + Z(t(2,h))%j(4)*A(t(1,h))%j(1))
      A_out(h)%h    = B01

    case (B00)
      A_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    case default
      A_out(h)%j(1) = g_RL(1) * ( - Z(t(2,h))%j(1)*A(t(1,h))%j(3) - Z(t(2,h))%j(3)*A(t(1,h))%j(4))
      A_out(h)%j(2) = g_RL(1) * ( - Z(t(2,h))%j(2)*A(t(1,h))%j(4) - Z(t(2,h))%j(4)*A(t(1,h))%j(3))
      A_out(h)%j(3) = g_RL(2) * ( - Z(t(2,h))%j(2)*A(t(1,h))%j(1) + Z(t(2,h))%j(3)*A(t(1,h))%j(2))
      A_out(h)%j(4) = g_RL(2) * ( - Z(t(2,h))%j(1)*A(t(1,h))%j(2) + Z(t(2,h))%j(4)*A(t(1,h))%j(1))
      A_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, A, Z, A_out, n, t)

end subroutine vert_AZ_Q


! **********************************************************************
subroutine vert_QA_Z(g_RL, ntry, Q, A, Z_out, n, t)
! bare QA -> Z Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! A(1:n(2))     = incoming anti-quark
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling gL
! Z_out(1:n(3)) = outgoing Z (light-cone representation)
! Z_out(h)%j(A) = A(t(2,h))%j(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * Q(t(1,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: Q(n(1)), A(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: Z_out(n(3))
  complex(REALKIND) :: A_aux(4)
  integer :: h

  do h = 1, n(3)
    select case (ishft(Q(t(1,h))%h,2) + A(t(2,h))%h)

    case (B1111)
      A_aux(1:2)    = g_RL(2)*A(t(2,h))%j(1:2)
      A_aux(3:4)    = g_RL(1)*A(t(2,h))%j(3:4)
      Z_out(h)%j(1) = - A_aux(1)*Q(t(1,h))%j(3) - A_aux(4)*Q(t(1,h))%j(2)
      Z_out(h)%j(2) = - A_aux(2)*Q(t(1,h))%j(4) - A_aux(3)*Q(t(1,h))%j(1)
      Z_out(h)%j(3) = - A_aux(1)*Q(t(1,h))%j(4) + A_aux(3)*Q(t(1,h))%j(2)
      Z_out(h)%j(4) = - A_aux(2)*Q(t(1,h))%j(3) + A_aux(4)*Q(t(1,h))%j(1)
      Z_out(h)%j    =  Z_out(h)%j + Z_out(h)%j

    case (B0110, B0111, B1110)
      A_aux(1:2)    = g_RL(2)*A(t(2,h))%j(1:2)
      Z_out(h)%j(1) = - A_aux(1)*Q(t(1,h))%j(3)
      Z_out(h)%j(2) = - A_aux(2)*Q(t(1,h))%j(4)
      Z_out(h)%j(3) = - A_aux(1)*Q(t(1,h))%j(4)
      Z_out(h)%j(4) = - A_aux(2)*Q(t(1,h))%j(3)
      Z_out(h)%j    =  Z_out(h)%j + Z_out(h)%j

    case (B1001, B1101, B1011)
      A_aux(3:4)    = g_RL(1)*A(t(2,h))%j(3:4)
      Z_out(h)%j(1) = - A_aux(4)*Q(t(1,h))%j(2)
      Z_out(h)%j(2) = - A_aux(3)*Q(t(1,h))%j(1)
      Z_out(h)%j(3) = + A_aux(3)*Q(t(1,h))%j(2)
      Z_out(h)%j(4) = + A_aux(4)*Q(t(1,h))%j(1)
      Z_out(h)%j    =  Z_out(h)%j + Z_out(h)%j

    case default
      Z_out(h)%j = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, Q, A, Z_out, n, t)

end subroutine vert_QA_Z


! **********************************************************************
subroutine vert_WQ_A(ntry, W, Q, Q_out, n, t)
! bare WQ -> Q W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! W(1:n(2))     = incoming W  (light-cone representation)
! Q_out(1:n(3)) = outgoing quark
! Q_out(h)%j(i) = W(t(2,h))%j(A) * [gamma_A*w_L](i,j) * Q(t(1,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: W(n(1)), Q(n(2))
  type(wfun),        intent(out)   :: Q_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (Q(t(2,h))%h)

    case (B01, B11)
    Q_out(h)%j(1) = - W(t(1,h))%j(2)*Q(t(2,h))%j(3) + W(t(1,h))%j(4)*Q(t(2,h))%j(4)
    Q_out(h)%j(2) = - W(t(1,h))%j(1)*Q(t(2,h))%j(4) + W(t(1,h))%j(3)*Q(t(2,h))%j(3)
     Q_out(h)%j(3:4) = 0
    Q_out(h)%h    = B10

    case default
    Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
    Q_out(h)%h = B00

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, W, Q, Q_out, n, t)

end subroutine vert_WQ_A


! **********************************************************************
subroutine vert_AW_Q(ntry, A, W, A_out, n, t)
! bare AW -> A W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! A(1:n(1))     = incoming anti-quark
! W(1:n(2))     = incoming W  (light-cone representation)
! A_out(1:n(3)) = outgoing anti-quark
! A_out(h)%j(i) = A(t(1,h))%j(j) * [gamma_A*w_L](j,i) * W(t(2,h))%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: A(n(1)), W(n(2))
  type(wfun),        intent(out)   :: A_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (A(t(1,h))%h)

    case (B10, B11)
      A_out(h)%j(1:2) = 0
      A_out(h)%j(3) = - W(t(2,h))%j(2)*A(t(1,h))%j(1) + W(t(2,h))%j(3)*A(t(1,h))%j(2)
      A_out(h)%j(4) = - W(t(2,h))%j(1)*A(t(1,h))%j(2) + W(t(2,h))%j(4)*A(t(1,h))%j(1)
      A_out(h)%h    = B01

    case default
      A_out(h)%j = 0  ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, A, W, A_out, n, t)

end subroutine vert_AW_Q


! **********************************************************************
subroutine vert_QA_W(ntry, Q, A, W_out, n, t)
! bare QA -> W W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! A(1:n(2))     = incoming anti-quark
! W_out(1:n(3)) = outgoing W  (light-cone representation)
! W_out(h)%j(A) = A(t(2,h))%j(i) * [gamma^A*w_L](i,j) * Q(t(1,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: Q(n(1)), A(n(2))
  type(wfun),        intent(out)   :: W_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (ishft(Q(t(1,h))%h,2) + A(t(2,h))%h)

    case (B1111, B1110, B0111, B0110)
      W_out(h)%j(1) = - A(t(2,h))%j(1)*Q(t(1,h))%j(3)
      W_out(h)%j(2) = - A(t(2,h))%j(2)*Q(t(1,h))%j(4)
      W_out(h)%j(3) = - A(t(2,h))%j(1)*Q(t(1,h))%j(4)
      W_out(h)%j(4) = - A(t(2,h))%j(2)*Q(t(1,h))%j(3)
      W_out(h)%j    = W_out(h)%j + W_out(h)%j

    case default
      W_out(h)%j = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, Q, A, W_out, n, t)

end subroutine vert_QA_W


! **********************************************************************
subroutine vert_VQ_A(ntry, V, Q, Q_out, n, t)
! bare VQ -> Q gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! V(1:n(1))     = incoming gluon (light-cone representation)
! Q(1:n(2))     = incoming quark
! Q_out(1:n(3)) = outgoing quark
! Q_out(h)%j(i) = V(t(1,h))%j(A) * gamma_A(i,j) * Q(t(2,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), Q(n(2))
  type(wfun),        intent(out)   :: Q_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (Q(t(2,h))%h)

    case (B01)
      Q_out(h)%j(1) = - V(t(1,h))%j(2)*Q(t(2,h))%j(3)+V(t(1,h))%j(4)*Q(t(2,h))%j(4)
      Q_out(h)%j(2) = - V(t(1,h))%j(1)*Q(t(2,h))%j(4)+V(t(1,h))%j(3)*Q(t(2,h))%j(3)
      Q_out(h)%j(3:4) = 0
      Q_out(h)%h    = B10

    case (B10)
      Q_out(h)%j(1:2) = 0
      Q_out(h)%j(3) = - V(t(1,h))%j(1)*Q(t(2,h))%j(1)-V(t(1,h))%j(4)*Q(t(2,h))%j(2)
      Q_out(h)%j(4) = - V(t(1,h))%j(2)*Q(t(2,h))%j(2)-V(t(1,h))%j(3)*Q(t(2,h))%j(1)
      Q_out(h)%h    = B01

    case (B00)
      Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      Q_out(h)%h = B00

    case default
      Q_out(h)%j(1) = - V(t(1,h))%j(2)*Q(t(2,h))%j(3)+V(t(1,h))%j(4)*Q(t(2,h))%j(4)
      Q_out(h)%j(2) = - V(t(1,h))%j(1)*Q(t(2,h))%j(4)+V(t(1,h))%j(3)*Q(t(2,h))%j(3)
      Q_out(h)%j(3) = - V(t(1,h))%j(1)*Q(t(2,h))%j(1)-V(t(1,h))%j(4)*Q(t(2,h))%j(2)
      Q_out(h)%j(4) = - V(t(1,h))%j(2)*Q(t(2,h))%j(2)-V(t(1,h))%j(3)*Q(t(2,h))%j(1)
      Q_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, V, Q, Q_out, n, t)

end subroutine vert_VQ_A


! **********************************************************************
subroutine vert_AV_Q(ntry, A, V, A_out, n, t)
! bare AV -> A gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! A(1:n(1))     = incoming anti-quark
! V(1:n(2))     = incoming gluon (light-cone representation)
! A_out(1:n(3)) = outgoing anti-quark
! A_out(h)%j(i) = A(t(1,h))%j(j) * gamma_A(j,i) * V(t(2,h))%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: A(n(1)), V(n(2))
  type(wfun),        intent(out)   :: A_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (A(t(1,h))%h)

    case (B01)
      A_out(h)%j(1) = - V(t(2,h))%j(1)*A(t(1,h))%j(3) - V(t(2,h))%j(3)*A(t(1,h))%j(4)
      A_out(h)%j(2) = - V(t(2,h))%j(2)*A(t(1,h))%j(4) - V(t(2,h))%j(4)*A(t(1,h))%j(3)
      A_out(h)%j(3:4) = 0
      A_out(h)%h    = B10

    case (B10)
      A_out(h)%j(1:2) = 0
      A_out(h)%j(3) = - V(t(2,h))%j(2)*A(t(1,h))%j(1) + V(t(2,h))%j(3)*A(t(1,h))%j(2)
      A_out(h)%j(4) = - V(t(2,h))%j(1)*A(t(1,h))%j(2) + V(t(2,h))%j(4)*A(t(1,h))%j(1)
      A_out(h)%h    = B01

    case (B00)
      A_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    case default
      A_out(h)%j(1) = - V(t(2,h))%j(1)*A(t(1,h))%j(3) - V(t(2,h))%j(3)*A(t(1,h))%j(4)
      A_out(h)%j(2) = - V(t(2,h))%j(2)*A(t(1,h))%j(4) - V(t(2,h))%j(4)*A(t(1,h))%j(3)
      A_out(h)%j(3) = - V(t(2,h))%j(2)*A(t(1,h))%j(1) + V(t(2,h))%j(3)*A(t(1,h))%j(2)
      A_out(h)%j(4) = - V(t(2,h))%j(1)*A(t(1,h))%j(2) + V(t(2,h))%j(4)*A(t(1,h))%j(1)
      A_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, A, V, A_out, n, t)

end subroutine vert_AV_Q


! **********************************************************************
subroutine vert_QA_V(ntry, Q, A, V_out, n, t)
! bare QA -> V gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! A(1:n(2))     = incoming anti-quark
! V_out(1:n(3)) = outgoing gluon
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: Q(n(1)), A(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (ishft(Q(t(1,h))%h,2) + A(t(2,h))%h)

    case (B1111)
      V_out(h)%j(1) = - A(t(2,h))%j(1)*Q(t(1,h))%j(3) - A(t(2,h))%j(4)*Q(t(1,h))%j(2)
      V_out(h)%j(2) = - A(t(2,h))%j(2)*Q(t(1,h))%j(4) - A(t(2,h))%j(3)*Q(t(1,h))%j(1)
      V_out(h)%j(3) = - A(t(2,h))%j(1)*Q(t(1,h))%j(4) + A(t(2,h))%j(3)*Q(t(1,h))%j(2)
      V_out(h)%j(4) = - A(t(2,h))%j(2)*Q(t(1,h))%j(3) + A(t(2,h))%j(4)*Q(t(1,h))%j(1)
      V_out(h)%j = V_out(h)%j + V_out(h)%j

    case (B0110, B0111, B1110)
      V_out(h)%j(1) = - A(t(2,h))%j(1)*Q(t(1,h))%j(3)
      V_out(h)%j(2) = - A(t(2,h))%j(2)*Q(t(1,h))%j(4)
      V_out(h)%j(3) = - A(t(2,h))%j(1)*Q(t(1,h))%j(4)
      V_out(h)%j(4) = - A(t(2,h))%j(2)*Q(t(1,h))%j(3)
      V_out(h)%j = V_out(h)%j + V_out(h)%j

    case (B1001, B1101, B1011)
      V_out(h)%j(1) = - A(t(2,h))%j(4)*Q(t(1,h))%j(2)
      V_out(h)%j(2) = - A(t(2,h))%j(3)*Q(t(1,h))%j(1)
      V_out(h)%j(3) =   A(t(2,h))%j(3)*Q(t(1,h))%j(2)
      V_out(h)%j(4) =   A(t(2,h))%j(4)*Q(t(1,h))%j(1)
      V_out(h)%j = V_out(h)%j + V_out(h)%j

    case default
      V_out(h)%j = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, Q, A, V_out, n, t)

end subroutine vert_QA_V


! **********************************************************************
subroutine vert_UV_W_mexpl(ntry, V1, P1, V2, P2, V_out, n, t)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! ntry           = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))     = incoming gluon (light-cone representation)
! Pi(4)          = incoming Vi momentum  (light-cone representation)
! V_out(1:n(3))  = outgoing gluon (light-cone representation)
! V_out(h)%j(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * V1(t(1,h))%j(a1) * V2(t(2,h))%j(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
!   use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)   :: ntry
  integer(intkind2), intent(inout):: n(3), t(2,n(3))
  type(wfun),        intent(in)   :: V1(n(1)), V2(n(2))
  complex(REALKIND), intent(in)   :: P1(4), P2(4)
  type(wfun),        intent(out)  :: V_out(n(3))
  complex(REALKIND) :: J1J2, P1J2(n(2)), P2J1(n(1))
  complex(REALKIND) :: P1half(4), P2half(4), P12(4), P112(4), P122(4)
  integer :: h

  P1half = 0.5_/**/REALKIND * P1
  P2half = 0.5_/**/REALKIND * P2
  P12  = P1half - P2half
  P112 = P1half + P1half + P2half
  P122 = P1half + P2half + P2half
  do h = 1, n(1)
    P2J1(h) = P122(1) * V1(h)%j(2) + P122(2) * V1(h)%j(1) &
            - P122(3) * V1(h)%j(4) - P122(4) * V1(h)%j(3)
  end do
  do h = 1, n(2)
    P1J2(h) = P112(1) * V2(h)%j(2) + P112(2) * V2(h)%j(1) &
             -P112(3) * V2(h)%j(4) - P112(4) * V2(h)%j(3)
  end do
  do h = 1, n(3)
    J1J2 = V1(t(1,h))%j(1) * V2(t(2,h))%j(2) + V1(t(1,h))%j(2) * V2(t(2,h))%j(1) &
          -V1(t(1,h))%j(3) * V2(t(2,h))%j(4) - V1(t(1,h))%j(4) * V2(t(2,h))%j(3)

    V_out(h)%j = J1J2 * P12 + P2J1(t(1,h)) * V2(t(2,h))%j - P1J2(t(2,h)) * V1(t(1,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, V1, V2, V_out, n, t)

end subroutine vert_UV_W_mexpl

! **********************************************************************
subroutine vert_UV_W_mids(ntry, V1, mom1, V2, mom2, V_out, n, t)
! bare VV -> V vertex
! ----------------------------------------------------------------------
! ntry           = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))     = incoming gluon (light-cone representation)
! Pi(4)          = incoming Vi momentum  (light-cone representation)
! V_out(1:n(3))  = outgoing gluon (light-cone representation)
! V_out(h)%j(a3) = {  g(a1,a2)*[P1-P2](a3) + g(a2,a3)*[P2+Pout](a1)
!               + g(a3,a1)*[-Pout-P1](a2)} * V1(t(1,h))%j(a1) * V2(t(2,h))%j(a2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
!   use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)   :: ntry
  integer,           intent(in)   :: mom1, mom2
  integer(intkind2), intent(inout):: n(3), t(2,n(3))
  type(wfun),        intent(in)   :: V1(n(1)), V2(n(2))
  type(wfun),        intent(out)  :: V_out(n(3))
  complex(REALKIND) :: J1J2, P1J2(n(2)), P2J1(n(1))
  complex(REALKIND) :: P1half(4), P2half(4), P12(4), P112(4), P122(4)
  integer :: h

  P1half = 0.5_/**/REALKIND * get_LC_4(mom1)
  P2half = 0.5_/**/REALKIND * get_LC_4(mom2)
  P12  = P1half - P2half
  P112 = P1half + P1half + P2half
  P122 = P1half + P2half + P2half
  do h = 1, n(1)
    P2J1(h) = P122(1) * V1(h)%j(2) + P122(2) * V1(h)%j(1) &
            - P122(3) * V1(h)%j(4) - P122(4) * V1(h)%j(3)
  end do
  do h = 1, n(2)
    P1J2(h) = P112(1) * V2(h)%j(2) + P112(2) * V2(h)%j(1) &
             -P112(3) * V2(h)%j(4) - P112(4) * V2(h)%j(3)
  end do
  do h = 1, n(3)
    J1J2 = V1(t(1,h))%j(1) * V2(t(2,h))%j(2) + V1(t(1,h))%j(2) * V2(t(2,h))%j(1) &
          -V1(t(1,h))%j(3) * V2(t(2,h))%j(4) - V1(t(1,h))%j(4) * V2(t(2,h))%j(3)

    V_out(h)%j = J1J2 * P12 + P2J1(t(1,h)) * V2(t(2,h))%j - P1J2(t(2,h)) * V1(t(1,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, V1, V2, V_out, n, t)

end subroutine vert_UV_W_mids


! **********************************************************************
subroutine vert_EV_V(ntry, V1, V2, V3, V_out, n, t)
! sigma vertex, where the sigma wave function is replaced
! by two gluon wave functions V1%j and V2%j
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))    = incoming sigmas
! V_out(1:n(3)) = outgoing sigma
! V_out(h)%j(d) = (g(a,c)*g(b,d) + g(1,4)*g(2,3)) * V1(t(1,h))%j(a)
!                               * V2(t(2,h))%j(b) * V3(t(3,h))%j(c)
!               = V1(t(1,h))%j.V3(t(3,h))%j * V2(t(2,h))%j(d)
!               + V2(t(2,h))%j.V3(t(3,h))%j * V1(t(1,h))%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h

  do h = 1, n(4)
    V_out(h)%j = cont_PP(V1(t(1,h))%j,V3(t(3,h))%j) * V2(t(2,h))%j &
               - cont_PP(V2(t(2,h))%j,V3(t(3,h))%j) * V1(t(1,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, V1, V2, V3, V_out, n, t)

end subroutine vert_EV_V


! **********************************************************************
subroutine vert_GGG_G(ntry, G1, G2, G3, G_out, n, t)
! Four-gluon vertex: factorised colour monomials f(a,b,x)*f(c,d,x) (same as EV_V)
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Gi(1:n(i))    = incoming gluons
! G_out(1:n(4)) = outgoing gluon
! G_out(h)%j(d) = (g(a,c)*g(b,d) + g(1,4)*g(2,3)) * G1(t(1,h))%j(a)
!                               * G2(t(2,h))%j(b) * G3(t(3,h))%j(c)
!               = G1(t(1,h))%j.G3(t(3,h))%j * G2(t(2,h))%j(d)
!               + G2(t(2,h))%j.G3(t(3,h))%j * G1(t(1,h))%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1),  intent(in)    :: ntry
  integer(intkind2),  intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)   :: G1(n(1)), G2(n(2)), G3(n(3))
  type(wfun), intent(out)  :: G_out(n(4))
  integer :: h

  do h = 1, n(4)
    G_out(h)%j = cont_PP(G1(t(1,h))%j,G3(t(3,h))%j)*G2(t(2,h))%j &
               - cont_PP(G2(t(2,h))%j,G3(t(3,h))%j)*G1(t(1,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, G1, G2, G3, G_out, n, t)

end subroutine vert_GGG_G


! **********************************************************************
subroutine vert_WWG_G(ntry, V1, V2, V3, V_out, n, t)
! Four-gluon vertex: factorised Lorentz monomials g(a,b)*g(c,d)
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))    = incoming gluons
! V_out(1:n(4)) = outgoing gluon
! V_out(h)%j(d) = g(a,b)*g(c,d) * V1%j(a) * V2(t(2,h))%j(b) * V3(t(3,h))%j(c)
!               = V1(t(1,h))%j.V2(t(2,h))%j * V3(t(3,h))%j(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h

  do h = 1, n(4)
    V_out(h)%j = cont_PP(V1(t(1,h))%j,V2(t(2,h))%j) * V3(t(3,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, V1, V2, V3, V_out, n, t)

end subroutine vert_WWG_G


! **********************************************************************
subroutine vert_WWV_V(ntry, V1, V2, V3, V_out, n, t)
! bare W+ W- A/Z A/Z vertex
! ----------------------------------------------------------------------
! ntry           = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))     = incoming vector bosons (light-cone representation)
! V_out(1:n(4))  = outgoing vector boson (light-cone representation)
! V_out(h)%j(a4) = [2*g(a1,a2)*g(a3,a4) - g(a2,a3)*g(a1,a4) - g(a1,a3)*g(a2,a4)]
!                  * V1(t(1,h))%j(a1) * V2(t(2,h))%j(a2) * V3(t(3,h))%j(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  complex(REALKIND) :: J1J2, J1J3, J2J3
  integer :: h

  do h = 1, n(4)
    J1J2 = cont_PP(V1(t(1,h))%j, V2(t(2,h))%j)
    J1J2 = J1J2 + J1J2
    J1J3 = cont_PP(V1(t(1,h))%j, V3(t(3,h))%j)
    J2J3 = cont_PP(V2(t(2,h))%j, V3(t(3,h))%j)
    V_out(h)%j = J1J2 * V3(t(3,h))%j - J2J3 * V1(t(1,h))%j - J1J3 * V2(t(2,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, V1, V2, V3, V_out, n, t)

end subroutine vert_WWV_V


! **********************************************************************
subroutine vert_VWW_V(ntry, V1, V2, V3, V_out, n, t)
! ntry          = 1 (2) for 1st (subsequent) PS points
! Vi(1:n(i))    = incoming vector bosons
! V_out(1:n(4)) = outgoing vector boson
! **********************************************************************

  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  complex(REALKIND) :: J2J3, J2J1, J3J1
  integer :: h

  do h = 1, n(4)
    J2J3 = cont_PP(V2(t(2,h))%j, V3(t(3,h))%j)
    J2J3 = J2J3 + J2J3
    J2J1 = cont_PP(V2(t(2,h))%j, V1(t(1,h))%j)
    J3J1 = cont_PP(V3(t(3,h))%j, V1(t(1,h))%j)
    V_out(h)%j = J2J3 * V1(t(1,h))%j - J3J1 * V2(t(2,h))%j - J2J1 * V3(t(3,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, V1, V2, V3, V_out, n, t)

end subroutine vert_VWW_V


! **********************************************************************
subroutine vert_SS_S(ntry, S1, S2, S_out, n, t)
! Three scalar vertex
! ----------------------------------------------------------------------
! ntry             : 1 (2) for 1st (subsequent) PS points
! Incoming scalars : Si(1:n(i))
! Outgoing scalar  : S_out(1:n(3))
!                    S_out(h)%j(1) = S1(t(1,h))%j(1) * S2(t(2,h))%j(1)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2))
  type(wfun),        intent(out)   :: S_out(n(3))
  integer :: h

  do h = 1, n(3)
    S_out(h)%j(1) = S1(t(1,h))%j(1) * S2(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, S1, S2, S_out, n, t)
  end if

end subroutine vert_SS_S


! **********************************************************************
subroutine vert_SSS_S(ntry, S1, S2, S3, S_out, n, t)
! Four scalar vertex
! ----------------------------------------------------------------------
! ntry             : 1 (2) for 1st (subsequent) PS points
! Incoming scalars : Si(1:n(i))
! Outgoing scalar  : S_out(1:n(3))
!                    S_out(h)%j(1) = S1(t(1,h))%j(1) * S2(t(2,h))%j(1)
!                                  * S3(t(3,h))%j(1)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), S3(n(3))
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h

  do h = 1, n(4)
    S_out(h)%j(1) = S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1)
  end do

  if (ntry == 1) then
    call  checkzero_scalar(S_out)
    call helbookkeeping_vert4(ntry, S1, S2, S3, S_out, n, t)
  end if

end subroutine vert_SSS_S


! **********************************************************************
subroutine vert_VS_T_mexpl(ntry, V, P1, S, P2, S_out, n, t)
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vector : V(1:n(1)),    incoming momentum P1(4) (light-cone representation)
! Incoming scalar : S(1:n(2)),    incoming momentum P2(4) (light-cone representation)
! Outgoing scalar : S_out(1:n(3))
!                   S_out(h)%j(1) = V(t(1,h))%j.(2*P2+P1) * S(t(2,h))%j(1)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), S(n(2))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  type(wfun),        intent(out)   :: S_out(n(3))
  complex(REALKIND) :: P122(4)
  integer :: h

  P122 = P1 + P2 + P2
  do h = 1, n(3)
    S_out(h)%j(1) = cont_PP(P122, V(t(1,h))%j) * S(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    call  checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, V, S, S_out, n, t)
  end if

end subroutine vert_VS_T_mexpl

! **********************************************************************
subroutine vert_VS_T_mids(ntry, V, mom1, S, mom2, S_out, n, t)
! Vector boson + two scalars vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1, mom2
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), S(n(2))
  type(wfun),        intent(out)   :: S_out(n(3))
  complex(REALKIND) :: P122(4)
  integer :: h

  P122 = get_LC_4(mom1 + mom2) + get_LC_4(mom2)
  do h = 1, n(3)
    S_out(h)%j(1) = cont_PP(P122, V(t(1,h))%j) * S(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    call  checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, V, S, S_out, n, t)
  end if

end subroutine vert_VS_T_mids


! **********************************************************************
subroutine vert_TV_S_mexpl(ntry, S, P1, V, P2, S_out, n, t)
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1)), incoming momentum P2(4) (light-cone representation)
! Incoming vector : V(1:n(2))    (light-cone representation)
! Outgoing scalar : S_out(1:n(3))
!                   S_out(h)%j(1) = V(t(2,h))%j.(-2*P1-P2) * S(t(1,h))%j(1)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  type(wfun),        intent(out)   :: S_out(n(3))
  complex(REALKIND) :: P112(4)
  integer :: h

  P112 = P1 + P1 + P2
  do h = 1, n(3)
    S_out(h)%j(1) = - cont_PP(P112, V(t(2,h))%j) * S(t(1,h))%j(1)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, S, V, S_out, n, t)
  end if

end subroutine vert_TV_S_mexpl

! **********************************************************************
subroutine vert_TV_S_mids(ntry, S, mom1, V, mom2, S_out, n, t)
! Vector boson + two scalars vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1, mom2
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  type(wfun),        intent(out)   :: S_out(n(3))
  complex(REALKIND) :: P112(4)
  integer :: h

  P112 = get_LC_4(mom1) + get_LC_4(mom1 + mom2)
  do h = 1, n(3)
    S_out(h)%j(1) = - cont_PP(P112, V(t(2,h))%j) * S(t(1,h))%j(1)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, S, V, S_out, n, t)
  end if

end subroutine vert_TV_S_mids

! **********************************************************************
subroutine vert_ST_V_mexpl(ntry, S1, P1, S2, P2, V_out, n, t)
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalars: Si(1:n(i)),  incoming momenta Pi(4) (light-cone representation)
! Outgoing vector : V_out(1:n(3))
!                   V_out(h)%j = S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * (P1 - P2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND) :: P12(4)
  integer :: h

  P12 = P1 - P2
  do h = 1, n(3)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1)) * P12
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S1, S2, V_out, n, t)

end subroutine vert_ST_V_mexpl

! **********************************************************************
subroutine vert_ST_V_mids(ntry, S1, mom1, S2, mom2, V_out, n, t)
! Vector boson + two scalars vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1, mom2
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND) :: P12(4)
  integer :: h

  P12 = get_LC_4(mom1 - mom2)
  do h = 1, n(3)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1)) * P12
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S1, S2, V_out, n, t)

end subroutine vert_ST_V_mids


! **********************************************************************
subroutine vert_VV_S(ntry, V1, V2, S_out, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vectors: Vi(1:n(i))  (light-cone representation)
! Outgoing scalar : S_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2))
  type(wfun),        intent(out)   :: S_out(n(3))
  integer :: h

  do h = 1, n(3)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j,V2(t(2,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, V1, V2, S_out, n, t)
  end if

end subroutine vert_VV_S


! **********************************************************************
subroutine vert_VS_V(ntry, V, S, V_out, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vector : V(1:n(1)) (light-cone representation)
! Incoming scalar : S(1:n(2))
! Outgoing vector : V_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), S(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  integer :: h

  do h = 1, n(3)
    V_out(h)%j = V(t(1,h))%j * S(t(2,h))%j(1)
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, V, S, V_out, n, t)

end subroutine vert_VS_V


! **********************************************************************
subroutine vert_SV_V(ntry, S, V, V_out, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1))
! Incoming vector : V(1:n(2)) (light-cone representation)
! Outgoing vector : V_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  integer :: h

  do h = 1, n(3)
    V_out(h)%j = S(t(1,h))%j(1) * V(t(2,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)

end subroutine vert_SV_V


! **********************************************************************
subroutine vert_VVS_S(ntry, V1, V2, S, S_out, n, t)
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vectors: Vi(1:n(i)) (light-cone representation)
! Incoming scalar : S(1:n(3))
! Outgoing scalar : S_out(1:n(4))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), S(n(3))
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h

  do h = 1, n(4)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j,V2(t(2,h))%j) * S(t(3,h))%j(1)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert4(ntry, V1, V2, S, S_out, n, t)
  end if

end subroutine vert_VVS_S


! **********************************************************************
subroutine vert_SSV_V(ntry, S1, S2, V, V_out, n, t)
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalars: Si(1:n(i))
! Incoming vector : V(1:n(3)) (light-cone representation)
! Outgoing vector : V_out(1:n(4)) (light-cone representation)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h

  do h = 1, n(4)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1)) * V(t(3,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, S1, S2, V, V_out, n, t)

end subroutine vert_SSV_V


! **********************************************************************
subroutine vert_VSS_V(ntry, V, S1, S2, V_out, n, t)
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vector : V(1:n(1)) (light-cone representation)
! Incoming scalars: Si(1:n(i))
! Outgoing vector : V_out(1:n(4)) (light-cone representation)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V(n(1)), S1(n(2)), S2(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h

  do h = 1, n(4)
    V_out(h)%j = (S1(t(2,h))%j(1) * S2(t(3,h))%j(1)) * V(t(1,h))%j
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, V, S1, S2, V_out, n, t)

end subroutine vert_VSS_V


! **********************************************************************
subroutine vert_SVV_S(ntry, S, V1, V2, S_out, n, t)
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1))
! Incoming vectors: Vi(1:n(i)) (light-cone representation)
! Outgoing scalar : S_out(1:n(4))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S(n(1)), V1(n(2)), V2(n(3))
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h

  do h = 1, n(4)
    S_out(h)%j(1) = S(t(1,h))%j(1) * cont_PP(V1(t(2,h))%j,V2(t(3,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert4(ntry, S, V1, V2, S_out, n, t)
  end if

end subroutine vert_SVV_S


! **********************************************************************
subroutine vert_AQ_S(g_RL, ntry, A, Q, S_out, n, t)
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! ntry                  : 1 (2) for 1st (subsequent) PS points
! right-handed coupling : g_RL(1)
! left-handed coupling  : g_RL(2)
! Incoming anti-fermion : A(1:n(1))
! Incoming fermion      : Q(1:n(2))
! Outgoing scalar       : S_out(1:n(3))
!                         S_out(h)%j(1) = gR*A(t(1,h))%j.P_R.Q(t(2,h))%j
!                                       + gL*A(t(1,h))%j.P_L.Q(t(2,h))%j
! with the chiral projectors P_R = (1+y5)/2 and P_L = (1-y5)/2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: A(n(1)), Q(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: S_out(n(3))
  integer :: h

  do h = 1, n(3)
    select case (ishft(Q(t(2,h))%h,2) + A(t(1,h))%h)

    case (B1111)
      S_out(h)%j(1) = g_RL(1) * (A(t(1,h))%j(1)*Q(t(2,h))%j(1) + A(t(1,h))%j(2)*Q(t(2,h))%j(2)) &
                    + g_RL(2) * (A(t(1,h))%j(3)*Q(t(2,h))%j(3) + A(t(1,h))%j(4)*Q(t(2,h))%j(4))

    case (B1010, B1110, B1011)
      S_out(h)%j(1) = g_RL(1) * (A(t(1,h))%j(1)*Q(t(2,h))%j(1) + A(t(1,h))%j(2)*Q(t(2,h))%j(2))

    case (B0101, B1101, B0111)
      S_out(h)%j(1) = g_RL(2) * (A(t(1,h))%j(3)*Q(t(2,h))%j(3) + A(t(1,h))%j(4)*Q(t(2,h))%j(4))

    case default
      S_out(h)%j(1) = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, A, Q, S_out, n, t)
  end if

end subroutine vert_AQ_S


! **********************************************************************
subroutine vert_QS_A(g_RL, ntry, Q, S, Q_out, n, t)
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! ntry                  : 1 (2) for 1st (subsequent) PS points
! right-handed coupling : g_RL(1)
! left-handed coupling  : g_RL(2)
! Incoming fermion      : Q(1:n(1))
! Incoming scalar       : S(1:n(2))
! Outgoing fermion      : Q_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: Q(n(1)), S(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: Q_out(n(3))
  complex(REALKIND) :: g_aux(2)
  integer :: h

  do h = 1, n(3)
    select case (Q(t(1,h))%h)

    case (B01)
      g_aux(2)   = g_RL(2) * S(t(2,h))%j(1)
      Q_out(h)%j(1:2) = 0
      Q_out(h)%j(3:4) = g_aux(2) * Q(t(1,h))%j(3:4)
      Q_out(h)%h    = B01

    case (B10)
      g_aux(1)   = g_RL(1) * S(t(2,h))%j(1)
      Q_out(h)%j(1:2) = g_aux(1) * Q(t(1,h))%j(1:2)
      Q_out(h)%j(3:4) = 0
      Q_out(h)%h    = B10

    case (B00)
      Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      Q_out(h)%h = B00

    case default
      g_aux        = g_RL * S(t(2,h))%j(1)
      Q_out(h)%j(1:2) = g_aux(1) * Q(t(1,h))%j(1:2)
      Q_out(h)%j(3:4) = g_aux(2) * Q(t(1,h))%j(3:4)
      Q_out(h)%h      = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, Q, S, Q_out, n, t)

end subroutine vert_QS_A


! **********************************************************************
subroutine vert_SA_Q(g_RL, ntry, S, A, A_out, n, t)
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! ntry                  : 1 (2) for 1st (subsequent) PS points
! right-handed coupling : g_RL(1)
! left-handed coupling  : g_RL(2)
! Incoming scalar       : S(1:n(1))
! Incoming anti-fermion : A(1:n(2))
! Outgoing anti-fermion : A_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), A(n(2))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: A_out(n(3))
  complex(REALKIND) :: g_aux(2)
  integer :: h

  do h = 1, n(3)
    select case (A(t(2,h))%h)

    case (B01)
      g_aux(2)   = g_RL(2) * S(t(1,h))%j(1)
      A_out(h)%j(1:2) = 0
      A_out(h)%j(3:4) = g_aux(2) * A(t(2,h))%j(3:4)
      A_out(h)%h    = B01

    case (B10)
      g_aux(1)   = g_RL(1) * S(t(1,h))%j(1)
      A_out(h)%j(1:2) = g_aux(1) * A(t(2,h))%j(1:2)
      A_out(h)%j(3:4) = 0
      A_out(h)%h    = B10

    case (B00)
      A_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    case default
      g_aux = g_RL * S(t(1,h))%j(1)
      A_out(h)%j(1:2) = g_aux(1) * A(t(2,h))%j(1:2)
      A_out(h)%j(3:4) = g_aux(2) * A(t(2,h))%j(3:4)
      A_out(h)%h      = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S, A, A_out, n, t)

end subroutine vert_SA_Q


! **********************************************************************
subroutine vert_CD_V(ntry, C, D, P2, V_out, n, t)
! Ghost/anti-ghost/gluon vertex
! ----------------------------------------------------------------------
! ntry                : 1 (2) for 1st (subsequent) PS points
! Incoming ghost      : C(1:n(1))
! Incoming anti-ghost : D(1:n(2)), momentum P2
! Outgoing gluon      : V_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: C(n(1)), D(n(2))
  complex(REALKIND), intent(in)    :: P2(4)
  type(wfun),        intent(out)   :: V_out(n(3))
  integer :: h

  do h = 1, n(3)
    V_out(h)%j = - C(t(1,h))%j(1) * D(t(2,h))%j(1) * P2
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, C, D, V_out, n, t)

end subroutine vert_CD_V


! **********************************************************************
subroutine vert_DV_C_mexpl(ntry, D, P1, V, D_out, n, t)
! Anti-ghost/gluon/anti-ghost vertex
! ----------------------------------------------------------------------
! ntry                : 1 (2) for 1st (subsequent) PS points
! Incoming anti-ghost : D(1:n(1)), momentum P1
! Incoming gluon      : V(1:n(2))
! Outgoing anti-ghost : D_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: D(n(1)), V(n(2))
  complex(REALKIND), intent(in)    :: P1(4)
  type(wfun),        intent(out)   :: D_out(n(3))
  integer :: h

  do h = 1, n(3)
    D_out(h)%j(1) = - D(t(1,h))%j(1) * cont_PP(P1, V(t(2,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(D_out)
    call helbookkeeping_vert3(ntry, D, V, D_out, n, t)
  end if

end subroutine vert_DV_C_mexpl

! **********************************************************************
subroutine vert_DV_C_mids(ntry, D, mom1, V, D_out, n, t)
! Anti-ghost/gluon/anti-ghost vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: D(n(1)), V(n(2))
  type(wfun),        intent(out)   :: D_out(n(3))
  integer :: h

  do h = 1, n(3)
    D_out(h)%j(1) = - D(t(1,h))%j(1) * cont_PP(get_LC_4(mom1), V(t(2,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(D_out)
    call helbookkeeping_vert3(ntry, D, V, D_out, n, t)
  end if

end subroutine vert_DV_C_mids

! **********************************************************************
subroutine vert_VC_D_mexpl(ntry, V, P1, C, P2, C_out, n, t)
! Ghost/gluon/ghost vertex
! ----------------------------------------------------------------------
! ntry           : 1 (2) for 1st (subsequent) PS points
! Incoming gluon : V(1:n(1)),    momentum P1
! Incoming ghost : C(1:n(2)),    momentum P2
! Outgoing ghost : C_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), C(n(2))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  type(wfun),        intent(out)   :: C_out(n(3))
  complex(REALKIND) :: P12(4)
  integer :: h

  P12 = P1 + P2
  do h = 1, n(3)
    C_out(h)%j(1) = C(t(2,h))%j(1) * cont_PP(P12, V(t(1,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(C_out)
    call helbookkeeping_vert3(ntry, V, C, C_out, n, t)
  end if

end subroutine vert_VC_D_mexpl

! **********************************************************************
subroutine vert_VC_D_mids(ntry, V, mom1, C, mom2, C_out, n, t)
! Ghost/gluon/ghost vertex
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1, mom2
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), C(n(2))
  type(wfun),        intent(out)   :: C_out(n(3))
  complex(REALKIND) :: P12(4)
  integer :: h

  P12 = get_LC_4(mom1 + mom2)
  do h = 1, n(3)
    C_out(h)%j(1) = C(t(2,h))%j(1) * cont_PP(P12, V(t(1,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(C_out)
    call helbookkeeping_vert3(ntry, V, C, C_out, n, t)
  end if

end subroutine vert_VC_D_mids


! ======================================================================
! HEFT vertices
! ======================================================================


! **********************************************************************
subroutine vert_GG_H(ntry, V1, P1, V2, P2, S_out, n, t)
! Effective gluon gluon -> Higgs vertex
! S_out = (P1.P2*g(1,2) - P1(2)*P2(1))*V1(1)*V2(2) = V1.V2*P1.P2 - P1.V2*P2.V1
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  type(wfun),        intent(out)   :: S_out(n(3))
  integer :: h
  do h = 1, n(3)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j, V2(t(2,h))%j) * cont_PP(P1, P2) &
                & - cont_PP(P1, V2(t(2,h))%j) * cont_PP(P2, V1(t(1,h))%j)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, V1, V2, S_out, n, t)
  end if
end subroutine vert_GG_H


! **********************************************************************
subroutine vert_HG_G(ntry, S, V2, P2, V_out, P3, n, t)
! Effective Higgs gluon -> gluon vertex
! V_out(3) = (P2(3)*P3(2) - P2.P3*g(2,3))*S*V2(2) = S*P3.V2*P2(3) - S*P2.P3*V2(3)
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V2(n(2))
  complex(REALKIND), intent(in)    :: P2(4), P3(4)
  type(wfun),        intent(out)   :: V_out(n(3))
  integer :: h
  do h = 1, n(3)
    V_out(h)%j = (S(t(1,h))%j(1) * cont_PP(P3, V2(t(2,h))%j)) * P2 &
             & - (S(t(1,h))%j(1) * cont_PP(P2, P3)) * V2(t(2,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert3(ntry, S, V2, V_out, n, t)
end subroutine vert_HG_G


! **********************************************************************
subroutine vert_GGG_H(ntry, V1, P1, V2, P2, V3, P3, S_out, n, t)
! Effective gluon gluon gluon -> Higgs vertex
! S_out = (g(1,2)*(P1-P2)(3) + g(2,3)*(P2-P3)(1) + g(3,1)*(P3-P1)(2)) * V1(1)*V2(2)*V3(3)
!       = V1.V2*(P1-P2).V3 + V2.V3*(P2-P3).V1 + V3.V1*(P3-P1).V2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4), P3(4)
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h
  do h = 1, n(4)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j, V2(t(2,h))%j) * cont_PP(P1-P2, V3(t(3,h))%j) &
                & + cont_PP(V2(t(2,h))%j, V3(t(3,h))%j) * cont_PP(P2-P3, V1(t(1,h))%j) &
                & + cont_PP(V3(t(3,h))%j, V1(t(1,h))%j) * cont_PP(P3-P1, V2(t(2,h))%j)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert4(ntry, V1, V2, V3, S_out, n, t)
  end if
end subroutine vert_GGG_H


! **********************************************************************
subroutine vert_HGG_G(ntry, S, V2, P2, V3, P3, V_out, P4, n, t)
! Effective Higgs gluon gluon -> gluon vertex
! ! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3-P4)(2) + g(4,2)*(P4-P2)(3)) * S*V2(2)*V3(3)
! !          = S*V2.V3*(P2-P3)(4) + S*(P3-P4).V2*V3(4) + S*(P4-P2).V3*V2(4)
! V_out(4) = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3+P4)(2) + g(4,2)*(-P4-P2)(3)) * S*V2(2)*V3(3)
!          = S*V2.V3*(P2-P3)(4) + S*(P3+P4).V2*V3(4) - S*(P4+P2).V3*V2(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S(n(1)), V2(n(2)), V3(n(3))
  complex(REALKIND), intent(in)    :: P2(4), P3(4), P4(4)
  type(wfun),        intent(out)   :: V_out(n(4))
  complex(REALKIND) :: P34_V2(n(2)), P42_V3(n(3))
  integer :: h
  do h = 1, n(2)
    P34_V2(h) = cont_PP(P3+P4, V2(h)%j)
  end do
  do h = 1, n(3)
    P42_V3(h) = cont_PP(P4+P2, V3(h)%j)
  end do
  do h = 1, n(4)
    V_out(h)%j = (S(t(1,h))%j(1) * cont_PP(V2(t(2,h))%j, V3(t(3,h))%j)) * (P2-P3) &
             & + (S(t(1,h))%j(1) * P34_V2(t(2,h))) * V3(t(3,h))%j &
             & - (S(t(1,h))%j(1) * P42_V3(t(3,h))) * V2(t(2,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert4(ntry, S, V2, V3, V_out, n, t)
end subroutine vert_HGG_G


! **********************************************************************
subroutine vert_GGGG_H(ntry, V1, V2, V3, V4, S_out, n, t)
! Effective gluon gluon gluon gluon -> Higgs vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(1,3)*g(2,4) - g(1,4)*g(2,3)) * V1(1) * V2(2) * V3(3) * V4(4)
!       = V1.V3*V2.V4 - V1.V4*V2.V3
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2)), V3(n(3)), V4(n(4))
  type(wfun),        intent(out)   :: S_out(n(5))
  integer :: h
  do h = 1, n(5)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j, V3(t(3,h))%j) * cont_PP(V2(t(2,h))%j, V4(t(4,h))%j) &
                & - cont_PP(V1(t(1,h))%j, V4(t(4,h))%j) * cont_PP(V2(t(2,h))%j, V3(t(3,h))%j)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert5(ntry, V1, V2, V3, V4, S_out, n, t)
  end if
end subroutine vert_GGGG_H


! **********************************************************************
subroutine vert_HGGG_G(ntry, S, V2, V3, V4, V_out, n, t)
! Effective Higgs gluon gluon gluon -> gluon vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(5) = (g(2,4)*g(3,5) - g(2,5)*g(3,4)) * S * V2(2) * V3(3) * V4(4)
!          = S*V2.V4*V3(5) - S*V3.V4*V2(5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: S(n(1)), V2(n(2)), V3(n(3)), V4(n(4))
  type(wfun),        intent(out)   :: V_out(n(5))
  integer :: h
  do h = 1, n(5)
    V_out(h)%j = (S(t(1,h))%j(1) * cont_PP(V2(t(2,h))%j, V4(t(4,h))%j)) * V3(t(3,h))%j &
             & - (S(t(1,h))%j(1) * cont_PP(V3(t(3,h))%j, V4(t(4,h))%j)) * V2(t(2,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert5(ntry, S, V2, V3, V4, V_out, n, t)
end subroutine vert_HGGG_G


! ======================================================================
! HHEFT vertices
! =====================================================================


! **********************************************************************
subroutine vert_HGG_H(ntry, S1, V2, P2, V3, P3, S_out, n, t)
! Effective Higgs gluon gluon -> Higgs vertex
! S_out = (P2.P3*g(2,3) - P2(3)*P3(2))*V2(2)*V3(3)*S1 = S1*(V2.V3*P2.P3 - P2.V3*P3.V2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S1(n(1))
  type(wfun),        intent(in)    :: V2(n(2)), V3(n(3))
  complex(REALKIND), intent(in)    :: P2(4), P3(4)
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h
  do h = 1, n(4)
    S_out(h)%j(1) = S1(t(1,h))%j(1)*(cont_PP(V2(t(2,h))%j, V3(t(3,h))%j) * cont_PP(P2, P3) &
                & - cont_PP(P2, V3(t(3,h))%j) * cont_PP(P3, V2(t(2,h))%j))
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert4(ntry, S1, V2, V3, S_out, n, t)
  end if
end subroutine vert_HGG_H


! **********************************************************************
subroutine vert_HHG_G(ntry, S1, S2, V3, P3, V_out, P4, n, t)
! Effective Higgs Higgs gluon -> gluon vertex
! V_out(4) = (P3(4)*P4(3) - P3.P4*g(3,4))*S1*S2*V3(3) = S*P4.V3*P3(4) - S1*S2*P3.P4*V3(4)
! note that P3 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V3(n(3))
  complex(REALKIND), intent(in)    :: P3(4), P4(4)
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h
  do h = 1, n(4)
    V_out(h)%j = S1(t(1,h))%j(1) * S2(t(2,h))%j(1)* (cont_PP(P4, V3(t(3,h))%j) * P3 &
             & - cont_PP(P3, P4) * V3(t(3,h))%j)
  end do
  if (ntry == 1) call helbookkeeping_vert4(ntry, S1, S2, V3, V_out, n, t)
end subroutine vert_HHG_G


! **********************************************************************
subroutine vert_HGGG_H(ntry, S1, V2, P2, V3, P3, V4, P4, S_out, n, t)
! Effective Higgs gluon gluon gluon -> Higgs vertex
! S_out = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3-P4)(2) + g(4,2)*(P4-P2)(3)) * S1*V2(2)*V3(3)*V4(4)
!       = S1*(V2.V3*(P2-P3).V4 + V3.V4*(P3-P4).V2 + V4.V2*(P4-P2).V3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: S1(n(1)), V2(n(2)), V3(n(3)), V4(n(4))
  complex(REALKIND), intent(in)    :: P2(4), P3(4), P4(4)
  type(wfun),        intent(out)   :: S_out(n(5))
  integer :: h
  do h = 1, n(5)
    S_out(h)%j(1) = cont_PP(V2(t(2,h))%j, V3(t(3,h))%j) * cont_PP(P2-P3, V4(t(4,h))%j) * S1(t(1,h))%j(1) &
                & + cont_PP(V3(t(3,h))%j, V4(t(4,h))%j) * cont_PP(P3-P4, V2(t(2,h))%j) * S1(t(1,h))%j(1) &
                & + cont_PP(V4(t(4,h))%j, V2(t(2,h))%j) * cont_PP(P4-P2, V3(t(3,h))%j) * S1(t(1,h))%j(1)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert5(ntry, S1, V2, V3, V4, S_out, n, t)
  end if
end subroutine vert_HGGG_H


! **********************************************************************
subroutine vert_HHGG_G(ntry, S1, S2, V3, P3, V4, P4, V_out, P5, n, t)
! Effective Higgs Higgs gluon gluon -> gluon vertex
! ! V_out(5) = (g(3,4)*(P3-P4)(5) + g(4,5)*(P4-P5)(3) + g(5,3)*(P5-P3)(4)) * S1*S2*V3(3)*V4(4)
! !          = S1*S2*V3.V4*(P3-P4)(5) + S1*S2*(P4-P5).V3*V4(5) + S1*S2*(P5-P3).V4*V3(5)
! V_out(5) = (g(3,4)*(P3-P4)(5) + g(4,5)*(P4+P5)(3) + g(5,3)*(-P5-P3)(4)) * S1*S2*V3(3)*V4(4)
!          = S1*2*V3.V4*(P3-P4)(5) + S1*S2*(P4+P5).V3*V4(5) - S1*S2*(P5+P3).V4*V3(5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V3(n(3)), V4(n(4))
  complex(REALKIND), intent(in)    :: P3(4), P4(4), P5(4)
  type(wfun),        intent(out)   :: V_out(n(5))
  complex(REALKIND) :: P45_V3(n(3)), P53_V4(n(4))
  integer :: h
  do h = 1, n(3)
    P45_V3(h) = cont_PP(P4+P5, V3(h)%j)
  end do
  do h = 1, n(4)
    P53_V4(h) = cont_PP(P5+P3, V4(h)%j)
  end do
  do h = 1, n(5)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * cont_PP(V3(t(3,h))%j, V4(t(4,h))%j)) * (P3-P4) &
             & + (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * P45_V3(t(3,h))) * V4(t(4,h))%j &
             & - (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * P53_V4(t(4,h))) * V3(t(3,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert5(ntry, S1, S2, V3, V4, V_out, n, t)
end subroutine vert_HHGG_G


! **********************************************************************
subroutine vert_HGGGG_H(ntry, S1, V2, V3, V4, V5, S_out, n, t)
! Effective Higgs gluon gluon gluon gluon -> Higgs vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(2,4)*g(3,5) - g(2,5)*g(3,4)) * S_in * V2(2) * V3(3) * V4(4) * V5(5)
!       = S_in * (V2.V4*V3.V5 - V2.V5*V3.V4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert6, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(6), t(5,n(6))
  type(wfun),        intent(in)    :: S1(n(1)), V2(n(2)), V3(n(3)), V4(n(4)), V5(n(5))
  type(wfun),        intent(out)   :: S_out(n(6))
  integer :: h
  do h = 1, n(6)
    S_out(h)%j(1) = S1(t(1,h))%j(1) * cont_PP(V2(t(2,h))%j, V4(t(4,h))%j) * cont_PP(V3(t(3,h))%j, V5(t(5,h))%j) &
                & - S1(t(1,h))%j(1) * cont_PP(V2(t(2,h))%j, V5(t(5,h))%j) * cont_PP(V3(t(3,h))%j, V4(t(4,h))%j)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert6(ntry, S1, V2, V3, V4, V5, S_out, n, t)
  end if
end subroutine vert_HGGGG_H


! **********************************************************************
subroutine vert_HHGGG_G(ntry, S1, S2, V3, V4, V5, V_out, n, t)
! Effective Higgs Higgs gluon gluon gluon -> gluon vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(6) = (g(3,5)*g(4,6) - g(3,6)*g(3,5)) * S1 * S2 * V3(3) * V4(4) * V5(5)
!          = S1*S2*V2.V4*V3(5) - S1*S2*V3.V4*V2(5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert6
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(6), t(5,n(6))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V3(n(3)), V4(n(4)), V5(n(5))
  type(wfun),        intent(out)   :: V_out(n(6))
  integer :: h
  do h = 1, n(6)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * cont_PP(V3(t(3,h))%j, V5(t(5,h))%j)) * V4(t(4,h))%j &
             & - (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * cont_PP(V4(t(4,h))%j, V5(t(5,h))%j)) * V3(t(3,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert6(ntry, S1, S2, V3, V4, V5, V_out, n, t)
end subroutine vert_HHGGG_G


! ======================================================================
! HHHEFT vertices
! =====================================================================


! **********************************************************************
subroutine vert_HHGG_H(ntry, S1, S2, V3, P3, V4, P4, S_out, n, t)
! Effective Higgs Higgs gluon gluon -> Higgs vertex
! S_out = (P3.P4*g(3,4) - P3(4)*P4(3))*V3(3)*V4(4)*S1*S2 = S1*S2*(V3.V4*P3.P4 - P3.V4*P4.V3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(5))
  type(wfun),        intent(in)    :: S1(n(1))
  type(wfun),        intent(in)    :: S2(n(2))
  type(wfun),        intent(in)    :: V3(n(3)), V4(n(4))
  complex(REALKIND), intent(in)    :: P3(4), P4(4)
  type(wfun),        intent(out)   :: S_out(n(4))
  integer :: h
  do h = 1, n(5)
    S_out(h)%j(1) = S1(t(1,h))%j(1)*S2(t(2,h))%j(1)*(cont_PP(V3(t(3,h))%j, V4(t(4,h))%j) * cont_PP(P3, P4) &
                & - cont_PP(P3, V4(t(4,h))%j) * cont_PP(P4, V3(t(3,h))%j))
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert5(ntry, S1, S2, V3, V4, S_out, n, t)
  end if
end subroutine vert_HHGG_H


! **********************************************************************
subroutine vert_HHHG_G(ntry, S1, S2, S3, V4, P4, V_out, P5, n, t)
! Effective Higgs Higgs gluon -> gluon vertex
! V_out(4) = (P4(5)*P5(4) - P4.P3*g(4,5))*S1*S2*S3*V4(4) = S1*S2*S3*P5.V4*P4(3) - S1*S2*S3*P4.P5*V4(3)
! note that P4 is outgoing
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert5
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(5), t(4,n(3))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), S3(n(3)), V4(n(4))
  complex(REALKIND), intent(in)    :: P4(4), P5(4)
  type(wfun),        intent(out)   :: V_out(n(5))
  integer :: h
  do h = 1, n(5)
    V_out(h)%j = S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * (cont_PP(P5, V4(t(4,h))%j) * P4 &
             & - cont_PP(P4, P5) * V4(t(4,h))%j)
  end do
  if (ntry == 1) call helbookkeeping_vert5(ntry, S1, S2, S3, V4, V_out, n, t)
end subroutine vert_HHHG_G


! **********************************************************************
subroutine vert_HHGGG_H(ntry, S1, S2, V3, P3, V4, P4, V5, P5, S_out, n, t)
! Effective Higgs Higgs gluon gluon gluon -> Higgs vertex
! S_out = (g(2,3)*(P2-P3)(4) + g(3,4)*(P3-P4)(2) + g(4,2)*(P4-P2)(3)) * S1*V2(2)*V3(3)*V4(4)
!       = S1*(V2.V3*(P2-P3).V4 + V3.V4*(P3-P4).V2 + V4.V2*(P4-P2).V3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert6, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(6), t(5,n(4))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V3(n(3)), V4(n(4)), V5(n(5))
  complex(REALKIND), intent(in)    :: P3(4), P4(4), P5(4)
  type(wfun),        intent(out)   :: S_out(n(6))
  integer :: h
  do h = 1, n(6)
    S_out(h)%j(1) = cont_PP(V3(t(3,h))%j, V4(t(4,h))%j) * cont_PP(P3-P4, V5(t(5,h))%j) * S1(t(1,h))%j(1) * S2(t(2,h))%j(1) &
                & + cont_PP(V4(t(4,h))%j, V5(t(5,h))%j) * cont_PP(P4-P3, V3(t(3,h))%j) * S1(t(1,h))%j(1) * S2(t(2,h))%j(1) &
                & + cont_PP(V5(t(5,h))%j, V3(t(3,h))%j) * cont_PP(P5-P3, V4(t(4,h))%j) * S1(t(1,h))%j(1) * S2(t(2,h))%j(1)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert6(ntry, S1, S2, V3, V4, V5, S_out, n, t)
  end if
end subroutine vert_HHGGG_H


! **********************************************************************
subroutine vert_HHHGG_G(ntry, S1, S2, S3, V4, P4, V5, P5, V_out, P6, n, t)
! Effective Higgs Higgs Higgs gluon gluon -> gluon vertex
! V_out(6) = (g(4,5)*(P4-P5)(6) + g(5,6)*(P5+P6)(4) + g(6,4)*(-P6-P4)(5)) * S1*S2*S3*V4(4)*V5(5)
!          = S1*S2*S3*V4.V5*(P4-P5)(6) + S1*S2*S3*(P5+P6).V4*V5(6) - S1*S2*S3*(P6+P4).V5*V4(6)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert6
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(6), t(5,n(6))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), S3(n(3)), V4(n(4)), V5(n(5))
  complex(REALKIND), intent(in)    :: P4(4), P5(4), P6(4)
  type(wfun),        intent(out)   :: V_out(n(6))
  complex(REALKIND) :: P56_V4(n(4)), P64_V5(n(5))
  integer :: h
  do h = 1, n(4)
    P56_V4(h) = cont_PP(P5+P6, V4(h)%j)
  end do
  do h = 1, n(5)
    P64_V5(h) = cont_PP(P6+P4, V5(h)%j)
  end do
  do h = 1, n(6)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * cont_PP(V4(t(4,h))%j, V5(t(5,h))%j)) * (P4-P5) &
             & + (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * P56_V4(t(4,h))) * V5(t(5,h))%j &
             & - (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * P64_V5(t(5,h))) * V4(t(4,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert6(ntry, S1, S2, S3, V4, V5, V_out, n, t)
end subroutine vert_HHHGG_G


! **********************************************************************
subroutine vert_HHGGGG_H(ntry, S1, S2, V3, V4, V5, V6, S_out, n, t)
! Effective Higgs Higgs gluon gluon gluon gluon -> Higgs vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! S_out = (g(3,5)*g(4,6) - g(3,6)*g(4,5)) * S1 * S2 * S3 * V3(3) * V4(4) * V5(5) * V6(6)
!       = S1 * S2 * S3 * (V3.V5*V4.V6 - V3.V6*V4.V5)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert7, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(7), t(6,n(7))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), V3(n(3)), V4(n(4)), V5(n(5)), V6(n(6))
  type(wfun),        intent(out)   :: S_out(n(7))
  integer :: h
  do h = 1, n(7)
    S_out(h)%j(1) = S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * cont_PP(V3(t(3,h))%j, V5(t(5,h))%j) * cont_PP(V4(t(4,h))%j, V6(t(6,h))%j) &
                & - S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * cont_PP(V3(t(3,h))%j, V6(t(6,h))%j) * cont_PP(V4(t(4,h))%j, V5(t(5,h))%j)
  end do
  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert7(ntry, S1, S2, V3, V4, V5, V6, S_out, n, t)
  end if
end subroutine vert_HHGGGG_H


! **********************************************************************
subroutine vert_HHHGGG_G(ntry, S1, S2, S3, V4, V5, V6, V_out, n, t)
! Effective Higgs Higgs Higgs gluon gluon gluon -> gluon vertex
! for factorised colour monomials f(a,b,x)*f(c,d,x)
! V_out(7) = (g(4,6)*g(5,7) - g(4,7)*g(4,6)) * S1 * S2 * S3 * V4(4) * V5(5) * V6(6)
!          = S1*S2*S3*V3.V5*V4(6) - S1*S2*S3*V4.V5*V3(6)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert7
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(7), t(6,n(7))
  type(wfun),        intent(in)    :: S1(n(1)), S2(n(2)), S3(n(3)), V4(n(4)), V5(n(5)), V6(n(6))
  type(wfun),        intent(out)   :: V_out(n(7))
  integer :: h
  do h = 1, n(7)
    V_out(h)%j = (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * cont_PP(V4(t(4,h))%j, V6(t(6,h))%j)) * V5(t(5,h))%j &
             & - (S1(t(1,h))%j(1) * S2(t(2,h))%j(1) * S3(t(3,h))%j(1) * cont_PP(V5(t(5,h))%j, V6(t(6,h))%j)) * V4(t(4,h))%j
  end do
  if (ntry == 1) call helbookkeeping_vert7(ntry, S1, S2, S3, V4, V5, V6, V_out, n, t)
end subroutine vert_HHHGGG_G



! **********************************************************************
! HiggsPO Effective couplings
! **********************************************************************

! **********************************************************************
subroutine vert_VV_S_kap(ntry, V1, P1, V2, P2, S_out, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vectors: Vi(1:n(i))  (light-cone representation)
! Outgoing scalar : S_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V1(n(1)), V2(n(2))
  type(wfun),        intent(out)   :: S_out(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  integer :: h

  do h = 1, n(3)
    S_out(h)%j(1) = cont_PP(V1(t(1,h))%j,V2(t(2,h))%j)
  end do

  if (ntry == 1) then
    call checkzero_scalar(S_out)
    call helbookkeeping_vert3(ntry, V1, V2, S_out, n, t)
  end if

end subroutine vert_VV_S_kap


! **********************************************************************
subroutine vert_VS_V_kap(ntry, V, P1, S, V_out, P2, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming vector : V(1:n(1)) (light-cone representation)
! Incoming scalar : S(1:n(2))
! Outgoing vector : V_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: V(n(1)), S(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  integer :: h

  do h = 1, n(3)
    V_out(h)%j = V(t(1,h))%j * S(t(2,h))%j(1)
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, V, S, V_out, n, t)

end subroutine vert_VS_V_kap


! **********************************************************************
subroutine vert_SV_V_kap(ntry, S, V, P1, V_out, P2, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1))
! Incoming vector : V(1:n(2)) (light-cone representation)
! Outgoing vector : V_out(1:n(3))
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  integer :: h
!
  do h = 1, n(3)
    V_out(h)%j = S(t(1,h))%j(1) * V(t(2,h))%j
  end do
!
  if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)
!
end subroutine vert_SV_V_kap
!

! ! **********************************************************************
! subroutine vert_SV_W_kap(ntry, S, V, P1, V_out, P2, n, t)
! ! Two vector boson + scalar vertex
! ! ----------------------------------------------------------------------
! ! ntry            : 1 (2) for 1st (subsequent) PS points
! ! Incoming scalar : S(1:n(1))
! ! Incoming vector : V(1:n(2)) (light-cone representation)
! ! Outgoing vector : V_out(1:n(3))
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND, intkind1, intkind2
!   use ol_data_types_/**/REALKIND, only: wfun
!   use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
!   implicit none
!   integer(intkind1), intent(in)    :: ntry
!   integer(intkind2), intent(inout) :: n(3), t(2,n(3))
!   type(wfun),        intent(in)    :: S(n(1)), V(n(2))
!   type(wfun),        intent(out)   :: V_out(n(3))
!   complex(REALKIND), intent(in)    :: P1(4), P2(4)
!   integer :: h
!
!   do h = 1, n(3)
!     V_out(h)%j = S(t(1,h))%j(1) * V(t(2,h))%j
!   end do
!
!   if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)
!
! end subroutine vert_SV_W_kap
!


! **********************************************************************
subroutine vert_SV_V_eps(ntry, S, V, P1, V_out, P2, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1))
! Incoming vector : V(1:n(2)) (light-cone representation)
! Outgoing vector : V_out(1:n(3)) = -(P1.P2)*V1*S - (P2.V1)*P1*S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_contractions_/**/REALKIND, only: cont_PP
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  complex(REALKIND) ::  P2V1(n(2))
  integer :: h

  do h = 1, n(2)
    P2V1(h) = cont_PP(P2, V(h)%j)
  end do

  do h = 1, n(3)
    V_out(h)%j = S(t(1,h))%j(1) * (-cont_PP(P1,P2)*V(t(2,h))%j+P2V1(h)*P1)
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)

end subroutine vert_SV_V_eps


! **********************************************************************
subroutine vert_SV_V_aeps(ntry, S, V, P1, V_out, P2, n, t)
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! ntry            : 1 (2) for 1st (subsequent) PS points
! Incoming scalar : S(1:n(1))
! Incoming vector : V(1:n(2)) (light-cone representation)
! Outgoing vector : V_out(1:n(3)) = ep(a,b,c,d)*J_V(b)*P(c)*P2(d)*S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_contractions_/**/REALKIND, only: cont_EpPPP
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: S(n(1)), V(n(2))
  type(wfun),        intent(out)   :: V_out(n(3))
  complex(REALKIND), intent(in)    :: P1(4), P2(4)
  integer :: h

  do h = 1, n(3)
    call cont_EpPPP(V(t(2,h))%j, P1, P2, V_out(h)%j)
    V_out(h)%j = S(t(1,h))%j(1) * V_out(h)%j
  end do

  if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)

end subroutine vert_SV_V_aeps

! ! **********************************************************************
! subroutine vert_VS_V_aeps(ntry, V, P1, S, V_out, P2, n, t)
! ! Two vector boson + scalar vertex
! ! ----------------------------------------------------------------------
! ! ntry            : 1 (2) for 1st (subsequent) PS points
! ! Incoming scalar : S(1:n(1))
! ! Incoming vector : V(1:n(2)) (light-cone representation)
! ! Outgoing vector : V_out(1:n(3)) = ep(a,b,c,d)*J_V(b)*P(c)*P2(d)
! ! **********************************************************************
!   use KIND_TYPES, only: REALKIND, intkind1, intkind2
!   use ol_data_types_/**/REALKIND, only: wfun
!   use ol_h_contractions_/**/REALKIND, only: cont_EpPPP
!   use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
!   implicit none
!   integer(intkind1), intent(in)    :: ntry
!   integer(intkind2), intent(inout) :: n(3), t(2,n(3))
!   type(wfun),        intent(in)    :: S(n(2)), V(n(1))
!   type(wfun),        intent(out)   :: V_out(n(3))
!   complex(REALKIND), intent(in)    :: P1(4), P2(4)
!   integer :: h
!
!   do h = 1, n(3)
!     call cont_EpPPP(V(t(1,h))%j, P1, P2, V_out(h)%j)
!     V_out(h)%j = S(t(2,h))%j(1) * V_out(h)%j
!   end do
!
!   if (ntry == 1) call helbookkeeping_vert3(ntry, S, V, V_out, n, t)
!
! end subroutine vert_VS_V_aeps
!

! **********************************************************************
subroutine vert_SQA_V(g, ntry, S, Q, A, V_out, n, t)
! bare SQA -> V gluon-like (i.e. vector-like) interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! A(1:n(2))     = incoming anti-quark
! V_out(1:n(3)) = outgoing gluon
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  real(REALKIND),    intent(in)    :: g
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S(n(1)), Q(n(2)), A(n(3))
  type(wfun),        intent(out)   :: V_out(n(4))
  integer :: h

  do h = 1, n(4)
    select case (ishft(Q(t(2,h))%h,2) + A(t(3,h))%h)

    case (B1111)
      V_out(h)%j(1) = - A(t(3,h))%j(1)*Q(t(2,h))%j(3) - A(t(3,h))%j(4)*Q(t(2,h))%j(2)
      V_out(h)%j(2) = - A(t(3,h))%j(2)*Q(t(2,h))%j(4) - A(t(3,h))%j(3)*Q(t(2,h))%j(1)
      V_out(h)%j(3) = - A(t(3,h))%j(1)*Q(t(2,h))%j(4) + A(t(3,h))%j(3)*Q(t(2,h))%j(2)
      V_out(h)%j(4) = - A(t(3,h))%j(2)*Q(t(2,h))%j(3) + A(t(3,h))%j(4)*Q(t(2,h))%j(1)
      V_out(h)%j = g*(V_out(h)%j + V_out(h)%j)*S(t(1,h))%j(1)

    case (B0110, B0111, B1110)
      V_out(h)%j(1) = - A(t(3,h))%j(1)*Q(t(2,h))%j(3)
      V_out(h)%j(2) = - A(t(3,h))%j(2)*Q(t(2,h))%j(4)
      V_out(h)%j(3) = - A(t(3,h))%j(1)*Q(t(2,h))%j(4)
      V_out(h)%j(4) = - A(t(3,h))%j(2)*Q(t(2,h))%j(3)
      V_out(h)%j = g*(V_out(h)%j + V_out(h)%j)*S(t(1,h))%j(1)

    case (B1001, B1101, B1011)
      V_out(h)%j(1) = - A(t(3,h))%j(4)*Q(t(2,h))%j(2)
      V_out(h)%j(2) = - A(t(3,h))%j(3)*Q(t(2,h))%j(1)
      V_out(h)%j(3) =   A(t(3,h))%j(3)*Q(t(2,h))%j(2)
      V_out(h)%j(4) =   A(t(3,h))%j(4)*Q(t(2,h))%j(1)
      V_out(h)%j = g*(V_out(h)%j + V_out(h)%j)*S(t(1,h))%j(1)

    case default
      V_out(h)%j = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, S, Q, A, V_out, n, t)

end subroutine vert_SQA_V


!************************************************************************
subroutine vert_ZSQ_A(g_RL, ntry, Z, S, Q, Q_out, n, t)
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Z(1:n(1)      = incoming Z (light-cone representation)
! S(1:n(2)      = incoming scalar
! Q(1:n(3))     = incoming quark
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling gL
! Q_out(1:n(3)) = outgoing quark
! Q_out(h)%j(i) = Z(t(1,h))%j(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*Q(t(2,h))%j(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: Z(n(1)), S(n(2)), Q(n(3))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: Q_out(n(4))
  integer :: h

  do h = 1, n(4)
    select case (Q(t(3,h))%h)

    case (B01)
      Q_out(h)%j(1) = g_RL(2) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(2)*Q(t(3,h))%j(3) + Z(t(1,h))%j(4)*Q(t(3,h))%j(4))
      Q_out(h)%j(2) = g_RL(2) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(1)*Q(t(3,h))%j(4) + Z(t(1,h))%j(3)*Q(t(3,h))%j(3))
      Q_out(h)%j(3:4) = 0
      Q_out(h)%h    = B10

    case (B10)
      Q_out(h)%j(1:2) = 0
      Q_out(h)%j(3) = g_RL(1) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(1)*Q(t(3,h))%j(1) - Z(t(1,h))%j(4)*Q(t(3,h))%j(2))
      Q_out(h)%j(4) = g_RL(1) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(2)*Q(t(3,h))%j(2) - Z(t(1,h))%j(3)*Q(t(3,h))%j(1))
      Q_out(h)%h    = B01

    case (B00)
      Q_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      Q_out(h)%h = B00

    case default
      Q_out(h)%j(1) = g_RL(2) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(2)*Q(t(3,h))%j(3) + Z(t(1,h))%j(4)*Q(t(3,h))%j(4))
      Q_out(h)%j(2) = g_RL(2) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(1)*Q(t(3,h))%j(4) + Z(t(1,h))%j(3)*Q(t(3,h))%j(3))
      Q_out(h)%j(3) = g_RL(1) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(1)*Q(t(3,h))%j(1) - Z(t(1,h))%j(4)*Q(t(3,h))%j(2))
      Q_out(h)%j(4) = g_RL(1) * S(t(2,h))%j(1) * ( - Z(t(1,h))%j(2)*Q(t(3,h))%j(2) - Z(t(1,h))%j(3)*Q(t(3,h))%j(1))
      Q_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, Z, S, Q, Q_out, n, t)

end subroutine vert_ZSQ_A


! **********************************************************************
subroutine vert_SQA_Z(g_RL, ntry, S, Q, A, Z_out, n, t)
! bare SQA -> Z Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! Q(1:n(1))     = incoming quark
! A(1:n(2))     = incoming anti-quark
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling g
! Z_out(1:n(3)) = outgoing Z
! Z_out(h)%j(A) = A(t(3,h))%j(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * Q(t(2,h))%j(j) * S(t(1,h))%j(1)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  complex(REALKIND),    intent(in)    :: g_RL(2)
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: S(n(1)), Q(n(2)), A(n(3))
  type(wfun),        intent(out)   :: Z_out(n(4))
  complex(REALKIND) :: A_aux(4)
  integer :: h

  do h = 1, n(4)
    select case (ishft(Q(t(2,h))%h,2) + A(t(3,h))%h)

    case (B1111)
      A_aux(1:2)    = g_RL(2)*A(t(3,h))%j(1:2)
      A_aux(3:4)    = g_RL(1)*A(t(3,h))%j(3:4)
      Z_out(h)%j(1) = - A_aux(1)*Q(t(2,h))%j(3) - A_aux(4)*Q(t(2,h))%j(2)
      Z_out(h)%j(2) = - A_aux(2)*Q(t(2,h))%j(4) - A_aux(3)*Q(t(2,h))%j(1)
      Z_out(h)%j(3) = - A_aux(1)*Q(t(2,h))%j(4) + A_aux(3)*Q(t(2,h))%j(2)
      Z_out(h)%j(4) = - A_aux(2)*Q(t(2,h))%j(3) + A_aux(4)*Q(t(2,h))%j(1)
      Z_out(h)%j    =  (Z_out(h)%j + Z_out(h)%j)*S(t(1,h))%j(1)

    case (B0110, B0111, B1110)
      A_aux(1:2)    = g_RL(2)*A(t(3,h))%j(1:2)
      Z_out(h)%j(1) = - A_aux(1)*Q(t(2,h))%j(3)
      Z_out(h)%j(2) = - A_aux(2)*Q(t(2,h))%j(4)
      Z_out(h)%j(3) = - A_aux(1)*Q(t(2,h))%j(4)
      Z_out(h)%j(4) = - A_aux(2)*Q(t(2,h))%j(3)
      Z_out(h)%j    =  (Z_out(h)%j + Z_out(h)%j)*S(t(1,h))%j(1)

    case (B1001, B1101, B1011)
      A_aux(3:4)    = g_RL(1)*A(t(3,h))%j(3:4)
      Z_out(h)%j(1) = - A_aux(4)*Q(t(2,h))%j(2)
      Z_out(h)%j(2) = - A_aux(3)*Q(t(2,h))%j(1)
      Z_out(h)%j(3) = + A_aux(3)*Q(t(2,h))%j(2)
      Z_out(h)%j(4) = + A_aux(4)*Q(t(2,h))%j(1)
      Z_out(h)%j    =  (Z_out(h)%j + Z_out(h)%j)*S(t(1,h))%j(1)

    case default
      Z_out(h)%j = 0 ! needed to detect vanishing helicity configurations

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, S, Q, A, Z_out, n, t)

end subroutine vert_SQA_Z


! **********************************************************************
subroutine vert_AZS_Q(g_RL, ntry, A, Z, S, A_out, n, t)
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! A(1:n(1))     = incoming anti-quark
! Z(1:n(2))     = incoming Z (light-cone representation)
! S(1:n(3)      = incoming scalar
! g_RL(1)       = right-handed coupling gR
! g_RL(2)       = left-handed coupling gL
! A_out(1:n(3)) = outgoing anti-quark
! A_out(h)%j(i) = A(t(1,h))%j(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * Z(t(2,h))%j(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun),        intent(in)    :: A(n(1)), Z(n(2)), S(n(3))
  complex(REALKIND), intent(in)    :: g_RL(2)
  type(wfun),        intent(out)   :: A_out(n(4))
  integer :: h

  do h = 1, n(4)
    select case (A(t(1,h))%h)

    case (B01)
      A_out(h)%j(1) = g_RL(1)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(1)*A(t(1,h))%j(3) - Z(t(2,h))%j(3)*A(t(1,h))%j(4))
      A_out(h)%j(2) = g_RL(1)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(2)*A(t(1,h))%j(4) - Z(t(2,h))%j(4)*A(t(1,h))%j(3))
      A_out(h)%j(3:4) = 0
      A_out(h)%h    = B10

    case (B10)
      A_out(h)%j(1:2) = 0
      A_out(h)%j(3) = g_RL(2)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(2)*A(t(1,h))%j(1) + Z(t(2,h))%j(3)*A(t(1,h))%j(2))
      A_out(h)%j(4) = g_RL(2)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(1)*A(t(1,h))%j(2) + Z(t(2,h))%j(4)*A(t(1,h))%j(1))
      A_out(h)%h    = B01

    case (B00)
      A_out(h)%j = 0 ! needed to detect vanishing helicity configurations
      A_out(h)%h = B00

    case default
      A_out(h)%j(1) = g_RL(1)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(1)*A(t(1,h))%j(3) - Z(t(2,h))%j(3)*A(t(1,h))%j(4))
      A_out(h)%j(2) = g_RL(1)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(2)*A(t(1,h))%j(4) - Z(t(2,h))%j(4)*A(t(1,h))%j(3))
      A_out(h)%j(3) = g_RL(2)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(2)*A(t(1,h))%j(1) + Z(t(2,h))%j(3)*A(t(1,h))%j(2))
      A_out(h)%j(4) = g_RL(2)*S(t(3,h))%j(1)*( - Z(t(2,h))%j(1)*A(t(1,h))%j(2) + Z(t(2,h))%j(4)*A(t(1,h))%j(1))
      A_out(h)%h    = B11

    end select
  end do

  if (ntry == 1) call helbookkeeping_vert4(ntry, A, Z, S, A_out, n, t)

end subroutine vert_AZS_Q


end module ol_h_vertices_/**/REALKIND
