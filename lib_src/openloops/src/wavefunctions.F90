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


module ol_wavefunctions_/**/REALKIND
  use KIND_TYPES, only: REALKIND, DREALKIND
  implicit none
  private
  public :: wf_S, wf_V, wf_V_Std, wf_Q, wf_A, wfIN_Q
  public :: pol_wf_S, pol_wf_V, pol_wf_Q, pol_wf_A
  real(REALKIND) :: small_real = 1.e-44_dp
  contains

! **********************************************************************
subroutine wf_S(P, M, POL, J_S)
! Wave function for a scalar particle. Just returns 1.
! (version without POLSEL kept for compatibility with old process code)
! **********************************************************************
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_S(4)
  J_S(1) = 1
  J_S(2:4) = 0
end subroutine wf_S

! **********************************************************************
subroutine pol_wf_S(P, M, POL, J_S, POLSEL)
! Wave function for a scalar particle. Just returns 1.
! (version without POLSEL kept for compatibility with old process code)
! **********************************************************************
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  integer, optional, intent(in)  :: POLSEL
  complex(REALKIND), intent(out) :: J_S(4)
  J_S(1) = 1
  J_S(2:4) = 0
end subroutine pol_wf_S

! **********************************************************************
subroutine wf_V(P, M, POL, J_V)
! vector boson wave function (incoming and outgoing)
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL:    -1|0|+1 polarisation
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0
! J_V(1:4) = EPS(P,POL)
!          = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0
! J_V(1:4) = EPS^*(-P,POL)
!          = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_V(4)
  integer :: i

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then
    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P, J_V)
    else
      ! normal wavefunction
      call wf_V_Std(P, M, POL, J_V)
    end if

  else
    call wf_V_Std(P, M, POL, J_V)
  end if

end subroutine wf_V


! **********************************************************************
subroutine pol_wf_V(P, M, POL, J_V, POLSEL)
! vector boson wave function (incoming and outgoing)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL:    -1|0|+1 polarisation
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0
! J_V(1:4) = EPS(P,POL)
!          = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0
! J_V(1:4) = EPS^*(-P,POL)
!          = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  integer, optional, intent(in)  :: POLSEL
  complex(REALKIND), intent(out) :: J_V(4)
  integer :: i

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then
    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P, J_V)
    else
      ! normal wavefunction
      call wf_V_Std(P, M, POL, J_V)
    end if

  else
    if (present(POLSEL)) then
      call wf_V_Std(P, M, POL, J_V, POLSEL)
    else
      call wf_V_Std(P, M, POL, J_V)
    end if
  end if

end subroutine pol_wf_V


! ! **********************************************************************
! subroutine wf_V(P, M, POL, J_V)
! ! vector boson wave function (incoming and outgoing)
! ! ----------------------------------------------------------------------
! ! P(0:3): incoming momentum (standard representation)
! ! POL:    -1|0|+1 polarisation
! ! M >= 0: real mass
! ! ----------------------------------------------------------------------
! ! if P(0) > 0
! ! J_V(1:4) = EPS(P,POL)
! !          = incoming vector boson wave function (light-cone representation)
! ! ----------------------------------------------------------------------
! ! if P(0) < 0
! ! J_V(1:4) = EPS^*(-P,POL)
! !          = outgoing vector boson wave function (light-cone representation)
! ! **********************************************************************
!   implicit none
!
!   real(REALKIND),    intent(in)  :: P(0:3), M
!   integer,     intent(in)  :: POL
!   complex(REALKIND), intent(out) :: J_V(4)
!   complex(REALKIND)  :: J_AUX(4)
!
!   if (P(0) >= 0) then ! incoming gluon -> EPS(P)
!     call wf_interface_V(P,M,POL,J_V)
! !     call wfIN_V(P,M,POL,J_V)
! !     call wfIN_V_MG(P,M,POL,J_V)
!   else if (P(0) < 0) then ! outgoing gluon -> [EPS(-P)]^*
!     call wf_interface_V(-P,M,POL,J_AUX)
! !     call wfIN_V(-P,M,POL,J_AUX)
! !     call wfIN_V_MG(P,M,POL,J_AUX)
!
!     J_V(1) = conjg(J_AUX(1))
!     J_V(2) = conjg(J_AUX(2))
!     J_V(3) = conjg(J_AUX(4)) ! light-cone conj: 3 <--> 4
!     J_V(4) = conjg(J_AUX(3))
!   end if
!
! end subroutine wf_V


! **********************************************************************
subroutine wf_V_Std(P, M, POL, J_V, POLSEL)
! wave function for IN/OUT vector boson
! ----------------------------------------------------------------------
! P(0:3) : incoming momentum P^mu (standard representation)
! POL    : -1|0|+1 polarisation
! M >= 0 : real mass
! ----------------------------------------------------------------------
! if P(0) > 0
! J_V(1:4) = EPS(P,POL)
!          = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0
! J_V(1:4) = EPS^*(-P,POL)
!          = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  integer, optional, intent(in)  :: POLSEL
  complex(REALKIND), intent(out) :: J_V(4)
  complex(REALKIND) :: J_AUX(4)

  if (P(0) >= 0) then ! incoming gluon -> EPS(P)
    call wfIN_V(P,M,POL,J_V, POLSEL)
!    call wf_interface_V(P,M,POL,J_V) ! gauge-fixing of Stefanos algebraic code
!    call wfIN_V_MG(P,M,POL,J_V) ! MadGraph convention
  else if (P(0) < 0) then ! outgoing gluon -> EPS^*(-P)
    call wfIN_V(-P,M,POL,J_AUX, POLSEL)
!    call wf_interface_V(-P,M,POL,J_AUX) ! gauge-fixing Stefanos algebraic code
!    call wfIN_V_MG(P,M,POL,J_AUX) ! MadGraph convention

    J_V(1) = conjg(J_AUX(1))
    J_V(2) = conjg(J_AUX(2))
    J_V(3) = conjg(J_AUX(4)) ! light-cone conj: 3 <--> 4
    J_V(4) = conjg(J_AUX(3))
  end if

end subroutine wf_V_Std



subroutine wf_interface_V(P,M,POL,J_V)
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_V(4)
  if (M == 0) then
    call wf_gf_V(P,POL,J_V)
  else
    call wfIN_V(P,M,POL,J_V)
  end if
end subroutine wf_interface_V



subroutine wf_gf_V(P,POL,J_V)
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: gf_array, inverse_crossing
  use ol_parameters_decl_/**/REALKIND, only: CI
  use ol_momenta_decl_/**/REALKIND, only: Q
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3)
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_V(4) ! Light-cone repr
  real(REALKIND)    :: Pmod, Pgfmod, cos_theta, sin_theta, Pgf(0:3), J_Lor(0:3), P_check(0:3,size(P)) ! Lorentz repr
  integer           :: i, part_i

  do i = 1, size(P_check)
  call LC2Std_Rep(Q(1:4,2**(i-1)), P_check(:,inverse_crossing(i)))
  end do

  part_i = 0

  do i = 1, size(P_check)
    if ( all(abs(P) - abs(P_check(:,i)) < 1d-10 ) ) then
      part_i = i ! identify the particle number to associate the gf_array(i)
      exit
    end if
  end do

  Pgf(1:3) = P_check(1:3,gf_array(part_i))
  Pgf(0) = sqrt(Pgf(1)*Pgf(1)+Pgf(2)*Pgf(2)+Pgf(3)*Pgf(3))! build it as light-like

  Pmod = sqrt(P(1)*P(1)+P(2)*P(2)+P(3)*P(3))
  Pgfmod = Pgf(0)
  cos_theta = (P(1)*Pgf(1)+P(2)*Pgf(2)+P(3)*Pgf(3))/(Pmod*Pgfmod)
  sin_theta = sqrt(1-cos_theta*cos_theta)

  if (sin_theta == 0) then
    if (POL == 1) then
      J_Lor(0)   = 0
      J_Lor(1)   = 1
      J_Lor(2:3) = 0
    else if (POL == -1) then
      J_Lor(0:1) = 0
      J_Lor(2)   = 1
      J_Lor(3)   = 0
    end if
  else
    if (POL == 1) then
      J_Lor(0) = 0
      J_Lor(1) = (P(2)*Pgf(3)-P(3)*Pgf(2))/(Pmod*Pgfmod)
      J_Lor(2) = (P(3)*Pgf(1)-P(1)*Pgf(3))/(Pmod*Pgfmod)
      J_Lor(3) = (P(1)*Pgf(2)-P(2)*Pgf(1))/(Pmod*Pgfmod)
    else if (POL == -1) then
      J_Lor(0) = 1 + cos_theta
      J_Lor(1:3) = P(1:3)/Pmod + Pgf(1:3)/Pgfmod
    end if
    J_Lor = J_Lor/sin_theta

  end if

  J_V(1) =   J_Lor(0) -      J_Lor(3)
  J_V(2) =   J_Lor(0) +      J_Lor(3)
  J_V(3) = - J_Lor(1) - CI * J_Lor(2) ! substitute for cmplx()
  J_V(4) = - J_Lor(1) + CI * J_Lor(2) ! substitute for cmplx()

end subroutine wf_gf_V


! **********************************************************************
subroutine wf_Q(P, M, POL, J_Q)
! wave function for an incoming quark or outgoing anti-quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M >= 0   : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! J_Q(1:4) = U(P,M,POL)
!          = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! J_Q(1:4) = V(-P,M,POL) = U(-P,-M,-POL)
!          = outgoing anti-quark wave function
! **********************************************************************
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_Q(4)

  if (P(0) >= 0) then ! in-quark -> U(P,M,POL)
    call wfIN_Q(P,M,POL,J_Q)
  else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL)=U(-P,-M,-POL)
    ! call wfIN_Q(-P,-M,-POL,J_Q)
    call wfIN_Q(-P,-M,POL,J_Q)
  end if
end subroutine wf_Q


! **********************************************************************
subroutine pol_wf_Q(P, M, POL, J_Q, POLSEL)
! wave function for an incoming quark or outgoing anti-quark
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M >= 0   : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! J_Q(1:4) = U(P,M,POL)
!          = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! J_Q(1:4) = V(-P,M,POL) = U(-P,-M,-POL)
!          = outgoing anti-quark wave function
! **********************************************************************
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_Q(4)
  integer, optional, intent(in)  :: POLSEL

  if (present(POLSEL)) then
    if (P(0) >= 0) then ! in-quark -> U(P,M,POL)
      call wfIN_Q(P,M,POL,J_Q, POLSEL)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL)=U(-P,-M,-POL)
      ! call wfIN_Q(-P,-M,-POL,J_Q)
      call wfIN_Q(-P,-M,POL,J_Q, POLSEL)
    end if
  else
    if (P(0) >= 0) then ! in-quark -> U(P,M,POL)
      call wfIN_Q(P,M,POL,J_Q)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL)=U(-P,-M,-POL)
      ! call wfIN_Q(-P,-M,-POL,J_Q)
      call wfIN_Q(-P,-M,POL,J_Q)
    end if
  end if

end subroutine pol_wf_Q


! **********************************************************************
subroutine wf_A(P, M, POL, J_A)
! wave function for incoming anti-quark or outgoing quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M>=0     : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0)>0
! J_A(1:4) = Vbar(P,M,POL) = Ubar(P,-M,-POL)
!          = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0)<0
! J_A(1:4) = Ubar(-P,M,POL)
!          = OUTGOING-antiquark wave function
! **********************************************************************
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: J_A(4)
  complex(REALKIND) :: J_AUX(4)

  if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL)=U(P,-M,-POL)
    call wfIN_Q(P,-M,-POL,J_AUX)
  else if(P(0) < 0) then ! out-quark -> U(-P,M,POL)
!     call wfIN_Q(-P,M,POL,J_AUX)
    call wfIN_Q(-P,M,-POL,J_AUX)
  end if

  ! Dirac conjugation of spinor
  J_A(1) = -conjg(J_AUX(3))
  J_A(2) = -conjg(J_AUX(4))
  J_A(3) = -conjg(J_AUX(1))
  J_A(4) = -conjg(J_AUX(2))

end subroutine wf_A


! **********************************************************************
subroutine pol_wf_A(P, M, POL, J_A, POLSEL)
! wave function for incoming anti-quark or outgoing quark
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M>=0     : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0)>0
! J_A(1:4) = Vbar(P,M,POL) = Ubar(P,-M,-POL)
!          = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0)<0
! J_A(1:4) = Ubar(-P,M,POL)
!          = OUTGOING-antiquark wave function
! **********************************************************************
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  integer,  optional,  intent(in)  :: POLSEL
  complex(REALKIND), intent(out) :: J_A(4)
  complex(REALKIND) :: J_AUX(4)

  if (present(POLSEL)) then
    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL)=U(P,-M,-POL)
      call wfIN_Q(P,-M,-POL,J_AUX,POLSEL)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL)
  !     call wfIN_Q(-P,M,POL,J_AUX)
      call wfIN_Q(-P,M,-POL,J_AUX,POLSEL)
    end if
  else
    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL)=U(P,-M,-POL)
      call wfIN_Q(P,-M,-POL,J_AUX)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL)
  !     call wfIN_Q(-P,M,POL,J_AUX)
      call wfIN_Q(-P,M,-POL,J_AUX)
    end if
  end if

  ! Dirac conjugation of spinor
  J_A(1) = -conjg(J_AUX(3))
  J_A(2) = -conjg(J_AUX(4))
  J_A(3) = -conjg(J_AUX(1))
  J_A(4) = -conjg(J_AUX(2))

end subroutine pol_wf_A


! **********************************************************************
subroutine wfIN_V(P, M, POL, EPS, POLSEL)
! wave function EPS(P,POL) for vector boson; P0 > 0
! as in (3.15) of hep-ph/9805445 (Dittmaier)
! I/O see subroutine wf_V
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND, only: CI, sqrt05
  use ol_debug, only: ol_fatal, ol_error
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: EPS(4)
  integer, optional, intent(in)  :: POLSEL
  real(REALKIND)    :: P2_T, P_T, P_MOD
  real(REALKIND)    :: SIN_THETA, ONEPCOS_THETA , ONEMCOS_THETA
  real(REALKIND)    :: COS_THETA, COS_PHI, SIN_PHI
  complex(REALKIND) :: EPHI_PLUS, EPHI_MINUS

  if (P(0) < 0) then
    call ol_fatal('in subroutine wfIN_V: P0 < 0 forbidden')
  end if

  if (present(POLSEL)) then
    if (POLSEL == 0) then
      continue
    else if (POL == 0 .and. POLSEL == 2) then
      continue
    else if (POL /= POLSEL) then
      EPS = small_real
      return
    end if
  end if

  P2_T  = P(1)**2 + P(2)**2
  P_T   = sqrt(P2_T)
  P_MOD = sqrt(P2_T + P(3)**2)

  if (P2_T == 0) then ! momentum along beam direction

    SIN_THETA = 0
    COS_PHI = 1
    SIN_PHI = 0
    if (P(3) > 0) then
      COS_THETA     = 1
      ONEPCOS_THETA = 2
      ONEMCOS_THETA = 0
    else if (P(3) <= 0) then
      COS_THETA     = -1
      ONEPCOS_THETA =  0
      ONEMCOS_THETA =  2
    end if
    EPHI_PLUS  = COS_PHI   ! COS_PHI+CI*SIN_PHI
    EPHI_MINUS = COS_PHI   ! COS_PHI-CI*SIN_PHI

  else if (P2_T > 0) then ! momentum not along beam direction

    SIN_THETA = P_T  / P_MOD
    COS_THETA = P(3) / P_MOD
    COS_PHI   = P(1) / P_T
    SIN_PHI   = P(2) / P_T
    if (P(3) > 0) then
      ONEPCOS_THETA = (P_MOD + P(3)) / P_MOD
      ONEMCOS_THETA = P2_T / (P_MOD * (P_MOD + P(3)))
    else if (P(3) <= 0) then
      ONEPCOS_THETA = P2_T / (P_MOD * (P_MOD - P(3)))
      ONEMCOS_THETA = (P_MOD - P(3)) / P_MOD
    end if
    EPHI_PLUS  = COS_PHI+CI*SIN_PHI ! (P(1) + CI*P(2)) / P_T
    EPHI_MINUS = COS_PHI-CI*SIN_PHI ! (P(1) - CI*P(2)) / P_T

  else
    call ol_error('in subroutine wfIN_V: P^2_T < 0 forbidden')
  end if

  if (POL == 1) then ! plus polarisation
    EPS(1) = - EPHI_MINUS * SIN_THETA * sqrt05
    EPS(2) = - EPS(1)
    EPS(3) = -ONEMCOS_THETA * sqrt05
    EPS(4) = (EPHI_MINUS**2) * ONEPCOS_THETA * sqrt05
  else if (POL == -1) then ! minus polarisation
    EPS(1) = - EPHI_PLUS * SIN_THETA * sqrt05
    EPS(2) = - EPS(1)
    EPS(3) = (EPHI_PLUS**2) * ONEPCOS_THETA * sqrt05
    EPS(4) = -ONEMCOS_THETA * sqrt05
  else if (POL == 0) then ! longitudinal polarisation
    EPS(1) = P(0)/M * (P_MOD/P(0) - COS_THETA)
    EPS(2) = P(0)/M * (P_MOD/P(0) + COS_THETA)
    EPS(3) = -P(0)/M * COS_PHI * SIN_THETA - CI * P(0)/M * SIN_PHI * SIN_THETA ! substitute for cmplx()
    EPS(4) = -P(0)/M * COS_PHI * SIN_THETA + CI * P(0)/M * SIN_PHI * SIN_THETA ! substitute for cmplx()
  end if

  ! workaround
  EPS = EPS + small_real

end subroutine wfIN_V


! **********************************************************************
subroutine wfIN_V_MG(P, M, POL, EPS)
! wave function EPS(P,POL) for vector boson; P0 > 0
! as in Appendix A.2 of KEK-91-11 (HELAS)
! I/O see subroutine wf_V
! Change wfIN_V into wfIN_V_MG in the subroutine wf_V to use this convention
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND, only: CI, sqrt05
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  complex(REALKIND), intent(out) :: EPS(4)
  real(REALKIND)    :: P2_T, P_T, P_MOD
  complex(REALKIND) :: ea(4), eb(4), epss(4)

  P2_T  = P(1)*P(1) + P(2)*P(2)
  P_T   = sqrt(P2_T)
  P_MOD = sqrt(P2_T + P(3)*P(3))

  if (POL == -1 .or. POL == 1) then

    if (P2_T == 0) then

      ea(1)   = 0
      ea(2)   = 1
      ea(3:4) = 0

      eb(1:2) = 0
      eb(3)   = P(3)/P_MOD
      eb(4)   = 0

    else

      ea(1)   =   0
      ea(2:3) =   P(1:2)*P(3)/(P_MOD*P_T)
      ea(4)   = - P_T/P_MOD

      eb(1) =   0
      eb(2) = - P(2)/P_T
      eb(3) =   P(1)/P_T
      eb(4) =   0

    end if

    if (POL == -1) then
      epss = - (ea + CI * eb) * sqrt05
    else if (POL == 1) then
      epss =   (ea - CI * eb) * sqrt05
    end if

  else if (POL == 0) then
    epss(1) =     P_MOD / M
    epss(2) = P(1)*P(0) / (M*P_MOD)
    epss(3) = P(2)*P(0) / (M*P_MOD)
    epss(4) = P(3)*P(0) / (M*P_MOD)
  end if

  EPS(1) =   epss(1) -      epss(4)
  EPS(2) =   epss(1) +      epss(4)
  EPS(3) = - epss(2) - CI * epss(3)
  EPS(4) = - epss(2) + CI * epss(3)

end subroutine wfIN_V_MG


! **********************************************************************
subroutine wfIN_Q(P, M, POL, WF, POLSEL)
! wave function U(P,M,POL) for incoming Quark; P(0) > 0; M = 0 or M > 0 or M < 0
! adapted from hep-ph/9805445 (Dittamier): m -> -m owing to Y.Z. conventions for Chiral rep.
! related to hep-ph/0002082 (HELAC) via Q -> POL*exp(POL*I*Phi)*Q
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND, only: CI
  use ol_debug, only: ol_fatal
  implicit none

  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL
  integer, optional,  intent(in)  :: POLSEL
  complex(REALKIND), intent(out) :: WF(4)
  real(REALKIND)    :: P2_T, P_T, P_MOD
  real(REALKIND)    :: COST, SINT, COSTHALF, SINTHALF, CHI, COSP, SINP
  real(REALKIND)    :: LAPLUS, LAMINUS
  complex(REALKIND) :: ZETA

  if (P(0) < 0) then
    call ol_fatal('in subroutine wfIN_Q: P0 < 0 forbidden')
  end if

  if (present(POLSEL)) then
    if (POLSEL == 0) then
      continue
    else if (POL /= POLSEL) then
      WF = 0
      return
    end if
  end if

  P2_T  = P(1)*P(1) + P(2)*P(2)
  P_T   = sqrt(P2_T)
  P_MOD = sqrt(P2_T + P(3)*P(3))

  if (P_MOD == 0) then
    COST = 1
    SINT = 0
  else
    COST = P(3) / P_MOD
    SINT = P_T  / P_MOD
  end if

  if (P2_T == 0) then
    COSP = 1
    SINP = 0
  else
    COSP = P(1) / P_T
    SINP = P(2) / P_T
  end if

  if (COST > 0) then
    COSTHALF = sqrt((1+COST)/2)
    SINTHALF = SINT / (2 * COSTHALF)
  else
    SINTHALF = sqrt((1-COST)/2)
    COSTHALF = SINT / (2 * SINTHALF)
  end if

  LAPLUS  = sqrt(P(0)+P_MOD)
  LAMINUS = M / LAPLUS

  if (POL == 1) then
    ZETA = COSTHALF*COSP - CI*COSTHALF*SINP ! substitute for cmplx()
    WF(1) =   LAPLUS  * ZETA     + small_real ! workaround
    WF(2) =   LAPLUS  * SINTHALF + small_real ! workaround
    WF(3) = - LAMINUS * ZETA
    WF(4) = - LAMINUS * SINTHALF
  else if (POL == -1) then
    ZETA  = COSTHALF*COSP + CI*COSTHALF*SINP ! substitute for cmplx()
    WF(1) = - LAMINUS * SINTHALF
    WF(2) =   LAMINUS * ZETA
    WF(3) =   LAPLUS  * SINTHALF + small_real ! workaround
    WF(4) = - LAPLUS  * ZETA     + small_real ! workaround
  end if

end subroutine wfIN_Q


end module ol_wavefunctions_/**/REALKIND



! **********************************************************************
module ol_s_wavefunctions_/**/REALKIND
! Routines to calculate external wave functions
! - wf_S: scalar
! - wf_Q: fermion
! - wf_A: anti-fermion
! - wf_V: vector boson
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  private
  public :: wf_S, wf_V, wf_Q, wf_A
  contains

! **********************************************************************
subroutine wf_S(P, M, POL, J_S)
! Wave function for a scalar particle.
! Just returns 1 as the first component of a wfun type.
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL
  type(wfun),     intent(out) :: J_S
  J_S%j(1) = 1
  J_S%j(2:4) = 0
end subroutine wf_S


! **********************************************************************
subroutine wf_Q(P, M, POL, Q)
! wave function for an incoming quark or outgoing anti-quark
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M >= 0   : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! Q%j(1:4) = U(P,M,POL)
!          = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! Q%j(1:4) = V(-P,M,POL) = U(-P,-M,-POL)
!          = outgoing anti-quark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL
  type(wfun),     intent(out) :: Q

  if (P(0) >= 0) then ! in-quark -> U(P,M,POL)
    call wfIN_Q(P,M,POL,Q%j)
  else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL)=U(-P,-M,-POL)
    call wfIN_Q(-P,-M,POL,Q%j)
  end if

  if(M /= 0)then
    Q%h = B"11"
  else
    if(POL == 1)then
      Q%h = B"10"
    else
      Q%h = B"01"
    end if
  end if

end subroutine wf_Q


! **********************************************************************
subroutine wf_A(P, M, POL, A)
! wave function for incoming anti-quark or outgoing quark
! ----------------------------------------------------------------------
! P(0:3)   : incoming momentum P^mu (standard representation)
! M>=0     : mass
! POL      : +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!            but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0)>0
! A%j(1:4) = Vbar(P,M,POL) = Ubar(P,-M,-POL)
!          = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0)<0
! A%j(1:4) = Ubar(-P,M,POL)
!          = OUTGOING-antiquark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL
  type(wfun),     intent(out) :: A
  complex(REALKIND) :: J_AUX(4)

  if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL)=U(P,-M,-POL)
    call wfIN_Q(P,-M,-POL,J_AUX)
  else if(P(0) < 0) then ! out-quark -> U(-P,M,POL)
    call wfIN_Q(-P,M,-POL,J_AUX)
  end if

  ! Dirac conjugation of spinor
  A%j(1) = -conjg(J_AUX(3))
  A%j(2) = -conjg(J_AUX(4))
  A%j(3) = -conjg(J_AUX(1))
  A%j(4) = -conjg(J_AUX(2))

  if(M /= 0)then
    A%h = B"11"
  else
    if(POL == 1)then
      A%h = B"10"
    else
      A%h = B"01"
    end if
  end if

end subroutine wf_A


! **********************************************************************
subroutine wf_V(P, M, POL, V)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL
  type(wfun),     intent(out) :: V
  integer :: i

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then

    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P,V%j)
    else
      ! normal wavefunction
      call wf_V_Std(P, M, POL, V%j)
    end if

  else

    call wf_V_Std(P, M, POL, V%j)

  end if

end subroutine wf_V

end module ol_s_wavefunctions_/**/REALKIND



! **********************************************************************
module ol_h_wavefunctions_/**/REALKIND
! Routines to calculate external wave functions:
! return an array of wave functions, one for each polarisation.
! - wf_S: scalar
! - wf_Q: fermion
! - wf_A: anti-fermion
! - wf_V: vector boson
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  private
  public :: wf_S, wf_V, wf_Q, wf_A
  public :: pol_wf_S, pol_wf_V, pol_wf_Q, pol_wf_A
  contains

! **********************************************************************
subroutine wf_S(P, M, POL, S)
! Wave function for a scalar particle.
! Just returns 1 in the component 1 of a 4 component wave function.
! (version without POLSEL kept for compatibility with old process code)
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(1) ! only 1 helicity state
  type(wfun),     intent(out) :: S(1)
  S(1)%j(1) = 1 ! S%j(2:4) components are not used
  S(1)%j(2:4) = 0
end subroutine wf_S


! **********************************************************************
subroutine pol_wf_S(P, M, POL, S, POLSEL)
! Wave function for a scalar particle.
! Just returns 1 in the component 1 of a 4 component wave function.
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)    :: P(0:3), M
  integer,        intent(in)    :: POL(1) ! only 1 helicity state
  integer, optional, intent(in) :: POLSEL
  type(wfun),     intent(out)   :: S(1)
  S(1)%j(1) = 1 ! S%j(2:4) components are not used
  S(1)%j(2:4) = 0
end subroutine pol_wf_S


! **********************************************************************
subroutine wf_Q(P, M, POL, Q)
! wave function for an incoming quark or outgoing anti-quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! Q(k)%j(1:4) = U(P,M,POL(k))
!             = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! Q(k)%j(1:4) = V(-P,M,POL(k)) = U(-P,-M,-POL(k))
!             = outgoing anti-quark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: Q(:)
  integer :: k

  do k = 1, size(POL)

    if (POL(k) == 99) then
      ! signal to ignore all remaining polarisation states
      Q(k:size(POL))%j(1) = 0
      Q(k:size(POL))%j(2) = 0
      Q(k:size(POL))%j(3) = 0
      Q(k:size(POL))%j(4) = 0
      Q(k:size(POL))%h = B"00"
      exit
    end if

    if (P(0) >= 0) then ! in-quark -> U(P,M,POL(k))
      call wfIN_Q(P, M, POL(k), Q(k)%j)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL(k)) = U(-P,-M,-POL(k))
      call wfIN_Q(-P, -M, POL(k), Q(k)%j)
    end if

    if(M /= 0)then
      Q(k)%h = B"11"
    else
      if(POL(k) == 1)then
        Q(k)%h = B"10"
      else
        Q(k)%h = B"01"
      end if
    end if

  end do

end subroutine wf_Q


! **********************************************************************
subroutine pol_wf_Q(P, M, POL, Q, POLSEL)
! wave function for an incoming quark or outgoing anti-quark
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! Q(k)%j(1:4) = U(P,M,POL(k))
!             = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! Q(k)%j(1:4) = V(-P,M,POL(k)) = U(-P,-M,-POL(k))
!             = outgoing anti-quark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  integer, optional,  intent(in)  :: POLSEL
  type(wfun),     intent(out) :: Q(:)
  integer :: k

  do k = 1, size(POL)

    if (POL(k) == 99) then
      ! signal to ignore all remaining polarisation states
      Q(k:size(POL))%j(1) = 0
      Q(k:size(POL))%j(2) = 0
      Q(k:size(POL))%j(3) = 0
      Q(k:size(POL))%j(4) = 0
      Q(k:size(POL))%h = B"00"
      exit
    end if

    if (present(POLSEL)) then
      if (POLSEL == 0) then
        continue
      else if (POL(k) /= POLSEL) then
        Q(k)%j(1) = 0
        Q(k)%j(2) = 0
        Q(k)%j(3) = 0
        Q(k)%j(4) = 0
        Q(k)%h = B"00"
        cycle
      end if
    end if

    if (P(0) >= 0) then ! in-quark -> U(P,M,POL(k))
      call wfIN_Q(P, M, POL(k), Q(k)%j)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL(k)) = U(-P,-M,-POL(k))
      call wfIN_Q(-P, -M, POL(k), Q(k)%j)
    end if

    if(M /= 0)then
      Q(k)%h = B"11"
    else
      if(POL(k) == 1)then
        Q(k)%h = B"10"
      else
        Q(k)%h = B"01"
      end if
    end if

  end do

end subroutine pol_wf_Q


! **********************************************************************
subroutine wf_A(P, M, POL, A)
! wave function for incoming anti-quark or outgoing quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! A(k)%j(1:4) = Vbar(P,M,POL(k)) = Ubar(P,-M,-POL(k))
!             = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! A(k)%j(1:4) = Ubar(-P,M,POL(k))
!             = OUTGOING-antiquark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none

  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: A(size(POL))
  complex(REALKIND) :: J_AUX(4)
  integer :: k

  do k = 1, size(POL)

    if(POL(k) == 99) then ! signal to ignore all remaining polarisation states
      A(k:size(POL))%j(1) = 0
      A(k:size(POL))%j(2) = 0
      A(k:size(POL))%j(3) = 0
      A(k:size(POL))%j(4) = 0
      A(k:size(POL))%h = B"00"
      exit
    end if

    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL(k))=U(P,-M,-POL(k))
      call wfIN_Q(P, -M, -POL(k), J_AUX)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL(k))
      call wfIN_Q(-P, M, -POL(k), J_AUX)
    end if
    ! Dirac conjugation of spinor
    A(k)%j(1) = -conjg(J_AUX(3))
    A(k)%j(2) = -conjg(J_AUX(4))
    A(k)%j(3) = -conjg(J_AUX(1))
    A(k)%j(4) = -conjg(J_AUX(2))

    if (M /= 0) then
      A(k)%h = B"11"
    else
      if (POL(k) == 1) then
        A(k)%h = B"10"
      else
        A(k)%h = B"01"
      end if
    end if

  end do

end subroutine wf_A


! **********************************************************************
subroutine pol_wf_A(P, M, POL, A, POLSEL)
! wave function for incoming anti-quark or outgoing quark
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! A(k)%j(1:4) = Vbar(P,M,POL(k)) = Ubar(P,-M,-POL(k))
!             = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! A(k)%j(1:4) = Ubar(-P,M,POL(k))
!             = OUTGOING-antiquark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none

  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  integer, optional, intent(in)  :: POLSEL
  type(wfun),     intent(out) :: A(size(POL))
  complex(REALKIND) :: J_AUX(4)
  integer :: k

  do k = 1, size(POL)

    if(POL(k) == 99) then ! signal to ignore all remaining polarisation states
      A(k:size(POL))%j(1) = 0
      A(k:size(POL))%j(2) = 0
      A(k:size(POL))%j(3) = 0
      A(k:size(POL))%j(4) = 0
      A(k:size(POL))%h = B"00"
      exit
    end if

    if (present(POLSEL)) then
      if (POLSEL == 0) then
        continue
      else if (-POL(k) /= POLSEL) then
        A(k)%j(1) = 0
        A(k)%j(2) = 0
        A(k)%j(3) = 0
        A(k)%j(4) = 0
        A(k)%h = B"00"
        cycle
      end if
    end if

    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL(k))=U(P,-M,-POL(k))
      call wfIN_Q(P, -M, -POL(k), J_AUX)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL(k))
      call wfIN_Q(-P, M, -POL(k), J_AUX)
    end if
    ! Dirac conjugation of spinor
    A(k)%j(1) = -conjg(J_AUX(3))
    A(k)%j(2) = -conjg(J_AUX(4))
    A(k)%j(3) = -conjg(J_AUX(1))
    A(k)%j(4) = -conjg(J_AUX(2))

    if (M /= 0) then
      A(k)%h = B"11"
    else
      if (POL(k) == 1) then
        A(k)%h = B"10"
      else
        A(k)%h = B"01"
      end if
    end if

  end do

end subroutine pol_wf_A


! **********************************************************************
subroutine wf_V(P, M, POL, V)
! vector boson wave function (incoming and outgoing)
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL(:): set of helicity states
! POL(k): +1|0|-1 as defined in subroutine wfIN_V
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0:
! V(k)%j(1:4) = EPS(P,POL(k))
!             = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0:
! V(k)%j(1:4) = EPS^*(-P,POL(k))
!             = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: V(:)
  integer :: i, k

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then

    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P,V(1)%j)
      do k = 2, size(POL)
        V(k)%j = 0
      end do
    else
      do k = 1, size(POL)
        if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          exit
        end if
        ! normal wavefunction
        call wf_V_Std(P, M, POL(k), V(k)%j)
      end do
    end if

  else

    do k = 1, size(POL)
      if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
        V(k:size(POL))%j(1) = 0
        V(k:size(POL))%j(2) = 0
        V(k:size(POL))%j(3) = 0
        V(k:size(POL))%j(4) = 0
        exit
      end if
      call wf_V_Std(P, M, POL(k), V(k)%j)
    end do

  end if

end subroutine wf_V


! **********************************************************************
subroutine pol_wf_V(P, M, POL, V, POLSEL)
! vector boson wave function (incoming and outgoing)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL(:): set of helicity states
! POL(k): +1|0|-1 as defined in subroutine wfIN_V
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0:
! V(k)%j(1:4) = EPS(P,POL(k))
!             = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0:
! V(k)%j(1:4) = EPS^*(-P,POL(k))
!             = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL(:)
  integer, optional, intent(in)  :: POLSEL
  type(wfun),        intent(out) :: V(:)
  integer :: i, k

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then

    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P,V(1)%j)
      do k = 2, size(POL)
        V(k)%j = 0
      end do
    else
      do k = 1, size(POL)
        if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          exit
        end if
        ! normal wavefunction
        call wf_V_Std(P, M, POL(k), V(k)%j)
      end do
    end if

  else

    do k = 1, size(POL)
      if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
        V(k:size(POL))%j(1) = 0
        V(k:size(POL))%j(2) = 0
        V(k:size(POL))%j(3) = 0
        V(k:size(POL))%j(4) = 0
        exit
      end if
      if (present(POLSEL)) then
        if (POLSEL == 0) then
          continue
        else if (POL(k) == 0 .and. POLSEL == 2) then
          continue
        else if (POL(k) /= POLSEL) then
          V(k)%j(1) = 0
          V(k)%j(2) = 0
          V(k)%j(3) = 0
          V(k)%j(4) = 0
          cycle
        end if
      end if
      call wf_V_Std(P, M, POL(k), V(k)%j)
    end do

  end if

end subroutine pol_wf_V


end module ol_h_wavefunctions_/**/REALKIND



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! The Following Module is used in the On-the-Fly Helicity Sum mode
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! **********************************************************************
module ol_hel_wavefunctions_/**/REALKIND
! Routines to calculate external wave functions:
! return an array of wave functions, one for each polarisation.
! - wf_S: scalar
! - wf_Q: fermion
! - wf_A: anti-fermion
! - wf_V: vector boson
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  private
  public :: wf_S, wf_V, wf_Q, wf_A
  public :: pol_wf_S, pol_wf_V, pol_wf_Q, pol_wf_A
  public :: init_hybrid_wf, init_hybrid_exwf
  contains

! **********************************************************************
subroutine wf_S(P, M, POL, S)
! Wave function for a scalar particle.
! Just returns 1 in the component 1 of a 4 component wave function.
! (version without POLSEL kept for compatibility with old process code)
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(1) ! only 1 helicity state
  type(wfun),     intent(out) :: S(1)
  S(1)%j(1) = 1 ! S%j(2:4) components are not used
  S(1)%j(2:4) = 0
end subroutine wf_S


! **********************************************************************
subroutine pol_wf_S(P, M, POL, S, POLSEL,w_ind)
! Wave function for a scalar particle.
! Just returns 1 in the component 1 of a 4 component wave function.
! w_ind = is the label for the external particle
! **********************************************************************
  use KIND_TYPES, only: intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(1) ! only 1 helicity state
  integer,        intent(in)  ::  w_ind
  integer, optional, intent(in) :: POLSEL
  type(wfun),     intent(out)   :: S(1)

  ! t = 2^(w_ind-1)
  S(1)%t = IBSET(0, w_ind-1)

  S(1)%j(1) = 1 ! S%j(2:4) components are not used
  S(1)%j(2:4) = 0

  S(1)%hf = (POL(1)+2)*IBSET(0,2*w_ind - 2)

end subroutine pol_wf_S


! **********************************************************************
subroutine wf_Q(P, M, POL, Q)
! wave function for an incoming quark or outgoing anti-quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! Q(k)%j(1:4) = U(P,M,POL(k))
!             = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! Q(k)%j(1:4) = V(-P,M,POL(k)) = U(-P,-M,-POL(k))
!             = outgoing anti-quark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: Q(:)
  integer :: k

  do k = 1, size(POL)

    if (POL(k) == 99) then
      ! signal to ignore all remaining polarisation states
      Q(k:size(POL))%j(1) = 0
      Q(k:size(POL))%j(2) = 0
      Q(k:size(POL))%j(3) = 0
      Q(k:size(POL))%j(4) = 0
      Q(k:size(POL))%h = B"00"
      exit
    end if

    if (P(0) >= 0) then ! in-quark -> U(P,M,POL(k))
      call wfIN_Q(P, M, POL(k), Q(k)%j)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL(k)) = U(-P,-M,-POL(k))
      call wfIN_Q(-P, -M, POL(k), Q(k)%j)
    end if

    if(M /= 0)then
      Q(k)%h = B"11"
    else
      if(POL(k) == 1)then
        Q(k)%h = B"10"
      else
        Q(k)%h = B"01"
      end if
    end if

  end do

end subroutine wf_Q


! **********************************************************************
subroutine pol_wf_Q(P, M, POL, Q, POLSEL, w_ind)
! wave function for an incoming quark or outgoing anti-quark
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing anti-quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! Q(k)%j(1:4) = U(P,M,POL(k))
!             = incoming quark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! Q(k)%j(1:4) = V(-P,M,POL(k)) = U(-P,-M,-POL(k))
!             = outgoing anti-quark wave function
! w_ind = is the label for the external particle
! **********************************************************************
  use KIND_TYPES, only: intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:), w_ind
  integer, optional,  intent(in)  :: POLSEL
  type(wfun),     intent(out) :: Q(:)
  integer :: k

  ! n_part = 1 is the number of external particles in the subtree
  Q(:)%n_part = 1_intkind2

  ! t = 2^(w_ind-1)
  Q(:)%t = IBSET(0, w_ind-1)

  do k = 1, size(POL)

    if (POL(k) == 99) then
      ! signal to ignore all remaining polarisation states
      Q(k:size(POL))%j(1) = 0
      Q(k:size(POL))%j(2) = 0
      Q(k:size(POL))%j(3) = 0
      Q(k:size(POL))%j(4) = 0
      Q(k:size(POL))%h = B"00"
      exit
    end if

    if (present(POLSEL)) then
      if (POLSEL == 0) then
        continue
      else if (POL(k) /= POLSEL) then
        Q(k)%j(1) = 0
        Q(k)%j(2) = 0
        Q(k)%j(3) = 0
        Q(k)%j(4) = 0
        Q(k)%h = B"00"
        cycle
      end if
    end if

    if (P(0) >= 0) then ! in-quark -> U(P,M,POL(k))
      call wfIN_Q(P, M, POL(k), Q(k)%j)
    else if (P(0) < 0) then ! out-antiquark -> V(-P,M,POL(k)) = U(-P,-M,-POL(k))
      call wfIN_Q(-P, -M, POL(k), Q(k)%j)
    end if

    if(M /= 0)then
      Q(k)%h = B"11"
    else
      if(POL(k) == 1)then
        Q(k)%h  = B"10"
      else
        Q(k)%h  = B"01"
      end if
    end if

  end do

  !initialisation of the global helicity label for the fermion wavefunction
  !hf = hel_i 4^(w_ind -1) with hel_i = {0,1,2,3}
  do k = 1, size(POL)
    if(present(POLSEL)) then
      if(POLSEL == 0 .or. POL(k) == POLSEL) then
        Q(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2)
      else if(POLSEL /= 0 .and. POL(k) /= POLSEL) then
        Q(k)%hf = -1_intkind2
      end if
    else
      Q(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2)
    end if
  end do

  if(present(POLSEL) .and. POLSEL /=0) call sort_hf_wf(Q)

end subroutine pol_wf_Q


! **********************************************************************
subroutine wf_A(P, M, POL, A)
! wave function for incoming anti-quark or outgoing quark
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! A(k)%j(1:4) = Vbar(P,M,POL(k)) = Ubar(P,-M,-POL(k))
!             = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! A(k)%j(1:4) = Ubar(-P,M,POL(k))
!             = OUTGOING-antiquark wave function
! **********************************************************************
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none

  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: A(size(POL))
  complex(REALKIND) :: J_AUX(4)
  integer :: k

  do k = 1, size(POL)

    if(POL(k) == 99) then ! signal to ignore all remaining polarisation states
      A(k:size(POL))%j(1) = 0
      A(k:size(POL))%j(2) = 0
      A(k:size(POL))%j(3) = 0
      A(k:size(POL))%j(4) = 0
      A(k:size(POL))%h = B"00"
      exit
    end if

    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL(k))=U(P,-M,-POL(k))
      call wfIN_Q(P, -M, -POL(k), J_AUX)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL(k))
      call wfIN_Q(-P, M, -POL(k), J_AUX)
    end if
    ! Dirac conjugation of spinor
    A(k)%j(1) = -conjg(J_AUX(3))
    A(k)%j(2) = -conjg(J_AUX(4))
    A(k)%j(3) = -conjg(J_AUX(1))
    A(k)%j(4) = -conjg(J_AUX(2))

    if (M /= 0) then
      A(k)%h = B"11"
    else
      if (POL(k) == 1) then
        A(k)%h = B"10"
      else
        A(k)%h = B"01"
      end if
    end if

  end do

end subroutine wf_A


! **********************************************************************
subroutine pol_wf_A(P, M, POL, A, POLSEL,w_ind)
! wave function for incoming anti-quark or outgoing quark
! ----------------------------------------------------------------------
! P(0:3): incoming momentum P^mu (standard representation)
! M >= 0: mass
! POL(:): set of helicity states to be summed
! POL(k): +1|-1 quark polarisation as (14,15) of hep-ph/0002082 (HELAC)
!         but with flipped polarisation for outgoing quarks
! ----------------------------------------------------------------------
! if P(0) > 0:
! A(k)%j(1:4) = Vbar(P,M,POL(k)) = Ubar(P,-M,-POL(k))
!             = INCOMING-antiquark wave function
! ----------------------------------------------------------------------
! if P(0) < 0:
! A(k)%j(1:4) = Ubar(-P,M,POL(k))
!             = OUTGOING-antiquark wave function
! **********************************************************************
  use KIND_TYPES, only: intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_wavefunctions_/**/REALKIND, only: wfIN_Q
  implicit none

  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:), w_ind
  integer, optional, intent(in)  :: POLSEL
  type(wfun),     intent(out) :: A(size(POL))
  complex(REALKIND) :: J_AUX(4)
  integer :: k

  ! n_part = 1 is the number of external particles in the subtree
  A(:)%n_part = 1_intkind2

  ! t = 2^(w_ind-1)
  A(:)%t = IBSET(0, w_ind-1)

  do k = 1, size(POL)

    if(POL(k) == 99) then ! signal to ignore all remaining polarisation states
      A(k:size(POL))%j(1) = 0
      A(k:size(POL))%j(2) = 0
      A(k:size(POL))%j(3) = 0
      A(k:size(POL))%j(4) = 0
      A(k:size(POL))%h = B"00"
      exit
    end if

    if (present(POLSEL)) then
      if (POLSEL == 0) then
        continue
      else if (-POL(k) /= POLSEL) then
        A(k)%j(1) = 0
        A(k)%j(2) = 0
        A(k)%j(3) = 0
        A(k)%j(4) = 0
        A(k)%h = B"00"
        cycle
      end if
    end if

    if (P(0) >= 0) then ! in-antiquark -> V(P,M,POL(k))=U(P,-M,-POL(k))
      call wfIN_Q(P, -M, -POL(k), J_AUX)
    else if(P(0) < 0) then ! out-quark -> U(-P,M,POL(k))
      call wfIN_Q(-P, M, -POL(k), J_AUX)
    end if
    ! Dirac conjugation of spinor
    A(k)%j(1) = -conjg(J_AUX(3))
    A(k)%j(2) = -conjg(J_AUX(4))
    A(k)%j(3) = -conjg(J_AUX(1))
    A(k)%j(4) = -conjg(J_AUX(2))

    if (M /= 0) then
      A(k)%h = B"11"
    else
      if (POL(k) == 1) then
        A(k)%h = B"10"
      else
        A(k)%h = B"01"
      end if
    end if

  end do

  do k = 1, size(POL)
    if(present(POLSEL)) then
      if(POLSEL == 0 .or. -POL(k) == POLSEL) then
        A(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2)
      else if(POLSEL /= 0 .and. -POL(k) /= POLSEL) then
        A(k)%hf = -1_intkind2
      end if
    else
      A(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2)
    end if
  end do

  if(present(POLSEL) .and. POLSEL /=0) call sort_hf_wf(A)

end subroutine pol_wf_A


! **********************************************************************
subroutine wf_V(P, M, POL, V)
! vector boson wave function (incoming and outgoing)
! (version without POLSEL kept for compatibility with old process code)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL(:): set of helicity states
! POL(k): +1|0|-1 as defined in subroutine wfIN_V
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0:
! V(k)%j(1:4) = EPS(P,POL(k))
!             = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0:
! V(k)%j(1:4) = EPS^*(-P,POL(k))
!             = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND), intent(in)  :: P(0:3), M
  integer,        intent(in)  :: POL(:)
  type(wfun),     intent(out) :: V(:)
  integer :: i, k

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then

    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P,V(1)%j)
      do k = 2, size(POL)
        V(k)%j = 0
      end do
    else
      do k = 1, size(POL)
        if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          exit
        end if
        ! normal wavefunction
        call wf_V_Std(P, M, POL(k), V(k)%j)
      end do
    end if

  else

    do k = 1, size(POL)
      if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
        V(k:size(POL))%j(1) = 0
        V(k:size(POL))%j(2) = 0
        V(k:size(POL))%j(3) = 0
        V(k:size(POL))%j(4) = 0
        exit
      end if
      call wf_V_Std(P, M, POL(k), V(k)%j)
    end do

  end if

end subroutine wf_V


! **********************************************************************
subroutine pol_wf_V(P, M, POL, V, POLSEL,w_ind)
! vector boson wave function (incoming and outgoing)
! ----------------------------------------------------------------------
! P(0:3): incoming momentum (standard representation)
! POL(:): set of helicity states
! POL(k): +1|0|-1 as defined in subroutine wfIN_V
! M >= 0: real mass
! ----------------------------------------------------------------------
! if P(0) > 0:
! V(k)%j(1:4) = EPS(P,POL(k))
!             = incoming vector boson wave function (light-cone representation)
! ----------------------------------------------------------------------
! if P(0) < 0:
! V(k)%j(1:4) = EPS^*(-P,POL(k))
!             = outgoing vector boson wave function (light-cone representation)
! **********************************************************************
  use KIND_TYPES, only: intkind2
  use ol_external_decl_/**/REALKIND, only: P_ex, Ward_array
  use ol_parameters_decl_/**/DREALKIND, only: Ward_tree, Ward_loop
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3), M
  integer,           intent(in)  :: POL(:), w_ind
  integer, optional, intent(in)  :: POLSEL
  type(wfun),        intent(out) :: V(:)
  integer :: i, k

  ! n_part = 1 is the number of external particles in the subtree
  V(:)%n_part = 1_intkind2

  ! t = 2^(w_ind-1)
  V(:)%t = IBSET(0,w_ind-1)

  if (Ward_tree /= 0 .or. Ward_loop /= 0) then

    do i = 1, size(P_ex,2)
      ! identify the particle number to associate the Ward_array(i)
      if ((P(0) >= 0 .and. all(P == P_ex(:,i))) .or. (P(0) < 0 .and. all(-P == P_ex(:,i)))) exit
    end do

    if (Ward_array(i) == 1) then
      call Std2LC_Rep(P,V(1)%j)
      do k = 2, size(POL)
        V(k)%j = 0
      end do
    else
      do k = 1, size(POL)
        if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          V(k:size(POL))%j(1) = 0
          exit
        end if
        ! normal wavefunction
        call wf_V_Std(P, M, POL(k), V(k)%j)
        V(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2) ! POL(k) = -1,0,1
      end do
    end if

  else

    do k = 1, size(POL)
      if (POL(k) == 99) then ! signal to ignore all remaining polarisation states
        V(k:size(POL))%j(1) = 0
        V(k:size(POL))%j(2) = 0
        V(k:size(POL))%j(3) = 0
        V(k:size(POL))%j(4) = 0
        exit
      end if
      if (present(POLSEL)) then
        if (POLSEL == 0) then
          continue
        else if (POL(k) == 0 .and. POLSEL == 2) then
          continue
        else if (POL(k) /= POLSEL) then
          V(k)%j(1) = 0
          V(k)%j(2) = 0
          V(k)%j(3) = 0
          V(k)%j(4) = 0
          cycle
        end if
      end if
      call wf_V_Std(P, M, POL(k), V(k)%j)
      V(k)%hf = (POL(k)+2)*IBSET(0,2*w_ind - 2) ! POL(k) = -1,0,1
    end do

  end if

end subroutine pol_wf_V


subroutine sort_hf_wf(wf)
  use KIND_TYPES, only: intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  type(wfun), intent(out) :: wf(:)
  type(wfun) :: wf_aux(size(wf))
  integer :: i, ntot

  ntot = 0
  wf_aux(:)%hf = -1_intkind2
  wf_aux(:)%h = B"00"
  wf_aux(:)%t = wf(:)%t
  wf_aux(:)%n_part = wf(:)%n_part
  wf_aux(:)%e = wf(:)%e

  do i = 1, size(wf)
    wf_aux(i)%j(:) = 0

    if(wf(i)%hf /= -1) then
      ntot = ntot + 1
      wf_aux(ntot)  = wf(i)
    end if
  end do

  wf = wf_aux

end subroutine sort_hf_wf

! **********************************************************************
! TODO:  <25-03-19, J.-N. Lang> !
! deprecated
subroutine init_hybrid_exwf(ex)
! **********************************************************************
  use KIND_TYPES, only: QREALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_loop_handling_/**/REALKIND, only: hp_switch
  implicit none
  type(wfun), intent(inout) :: ex(:)

#ifdef PRECISION_dp
  !if (hp_switch .eq. 1) then
    !ex(:)%j_qp(1) = cmplx(ex(:)%j(1),kind=qp)
    !ex(:)%j_qp(2) = cmplx(ex(:)%j(2),kind=qp)
    !ex(:)%j_qp(3) = cmplx(ex(:)%j(3),kind=qp)
    !ex(:)%j_qp(4) = cmplx(ex(:)%j(4),kind=qp)
  !end if
#endif
end subroutine init_hybrid_exwf

! **********************************************************************
! TODO:  <25-03-19, J.-N. Lang> !
! deprecated
subroutine init_hybrid_wf(wf)
! **********************************************************************
  use KIND_TYPES, only: QREALKIND
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_loop_handling_/**/REALKIND, only: hp_switch
  implicit none
  type(wfun), intent(inout) :: wf(:,:)
  integer :: swf1,swf2

#ifdef PRECISION_dp
  !if (hp_switch .eq. 1) then
    !swf1 = size(wf,1)
    !swf2 = size(wf,2)
    !wf(:swf1,:swf2)%j_qp(1) = cmplx(wf(:swf1,:swf2)%j(1),kind=qp)
    !wf(:swf1,:swf2)%j_qp(2) = cmplx(wf(:swf1,:swf2)%j(2),kind=qp)
    !wf(:swf1,:swf2)%j_qp(3) = cmplx(wf(:swf1,:swf2)%j(3),kind=qp)
    !wf(:swf1,:swf2)%j_qp(4) = cmplx(wf(:swf1,:swf2)%j(4),kind=qp)
  !end if
#endif
end subroutine init_hybrid_wf

end module ol_hel_wavefunctions_/**/REALKIND
