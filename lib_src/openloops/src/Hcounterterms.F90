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


module ol_h_counterterms_/**/REALKIND

  implicit none

  interface counter_Q_A
    module procedure counter_Q_A_orig,counter_Q_A_pid
  end interface

  interface counter_A_Q
    module procedure counter_A_Q_orig,counter_A_Q_pid
  end interface

  interface counter_V_V
    module procedure counter_V_V_orig,counter_V_V_pid
  end interface

contains

! **********************************************************************
subroutine counter_Q_A_orig(ctQA, ntry, J_Q, mom, Jout_Q, n)
! Q -> Q counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(2)
  type(wfun),        intent(inout) :: J_Q(n)
  type(wfun),        intent(out)   :: Jout_Q(n)

  call counter_Q_A(ctQA, 0, ntry, J_Q, mom, Jout_Q, n)
end subroutine counter_Q_A_orig

! **********************************************************************
subroutine counter_Q_A_pid(ctQA, pid, ntry, J_Q, mom, Jout_Q, n)
! Q -> Q counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: DREALKIND,REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_qcd_offshell_selfenergies/**/REALKIND, only: quark_ofsse
  use ol_momenta_decl_/**/REALKIND, only: softconf
  use ol_parameters_decl_/**/DREALKIND, only: bubble_vertex, ir_hacks
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_kinematics_/**/REALKIND, only: get_LC_5, get_LC_mass2
  use ol_debug, only: ol_error
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: pid,mom
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(2)
  type(wfun),        intent(inout) :: J_Q(n)
  type(wfun),        intent(out)   :: Jout_Q(n)
  complex(REALKIND) :: GG(2), K(5)
  integer :: h,momh
  logical :: use_soft_rel

  K = get_LC_5(mom)
  GG = 0
  if (bubble_vertex .eq. 1) GG = quark_ofsse(get_LC_mass2(mom),pid)

  ! TODO: clean  <27-11-18, J.-N. Lang> !
  ! hack to cure severe soft/collinear instabilities in 2 point fermionic propagator
  use_soft_rel = .false.
  if (ir_hacks .and. softconf .ne. 0) then
    momh = mom-softconf
    ! check if a soft particle is attached to an external particle
    if (iand(mom,softconf) .ne. 0 .and. iand(mom-softconf,mom-softconf-1) .eq. 0) then
      ! currently only supported for massless particles
      if (get_LC_mass2(momh) .eq. 0 .and. associated(J_Q(1)%j_prev)) then
      !if (get_LC_mass2(momh) .eq. 0) then
        use_soft_rel = .true.
      end if
    end if
  end if

  if (use_soft_rel) then
    do h = 1, n
      Jout_Q(h)%j(1) = (ctQA(1)+GG(1)) * K(5) *J_Q(h)%j_prev(1) - (ctQA(2)+GG(2)) * J_Q(h)%j(1)
      Jout_Q(h)%j(2) = (ctQA(1)+GG(1)) * K(5) *J_Q(h)%j_prev(2) - (ctQA(2)+GG(2)) * J_Q(h)%j(2)
      Jout_Q(h)%j(3) = (ctQA(1)+GG(1)) * K(5) *J_Q(h)%j_prev(3) - (ctQA(2)+GG(2)) * J_Q(h)%j(3)
      Jout_Q(h)%j(4) = (ctQA(1)+GG(1)) * K(5) *J_Q(h)%j_prev(4) - (ctQA(2)+GG(2)) * J_Q(h)%j(4)

      Jout_Q(h)%h = B11 ! default value for the label h

      deallocate(J_Q(h)%j_prev)
    end do
  else
    do h = 1, n

      Jout_Q(h)%j(1) = (ctQA(1)+GG(1)) * ( - K(2)*J_Q(h)%j(3) + K(4)*J_Q(h)%j(4)) - (ctQA(2)+GG(2)) * J_Q(h)%j(1)
      Jout_Q(h)%j(2) = (ctQA(1)+GG(1)) * ( - K(1)*J_Q(h)%j(4) + K(3)*J_Q(h)%j(3)) - (ctQA(2)+GG(2)) * J_Q(h)%j(2)
      Jout_Q(h)%j(3) = (ctQA(1)+GG(1)) * ( - K(1)*J_Q(h)%j(1) - K(4)*J_Q(h)%j(2)) - (ctQA(2)+GG(2)) * J_Q(h)%j(3)
      Jout_Q(h)%j(4) = (ctQA(1)+GG(1)) * ( - K(2)*J_Q(h)%j(2) - K(3)*J_Q(h)%j(1)) - (ctQA(2)+GG(2)) * J_Q(h)%j(4)

      Jout_Q(h)%h = B11 ! default value for the label h

    end do
  end if

  if (ntry == 1) then
    Jout_Q(:)%t = J_Q(1)%t
    Jout_Q(:)%n_part = J_Q(1)%n_part
    Jout_Q(:)%hf = J_Q(:)%hf
    call helbookkeeping_prop(ntry, J_Q, Jout_Q, n)
  end if

end subroutine counter_Q_A_pid

! **********************************************************************
subroutine counter_A_Q_orig(ctQA, ntry, J_A, mom, Jout_A, n)
! A -> A counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: DREALKIND, REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer          , intent(in)    :: mom
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(2)
  type(wfun),        intent(inout) :: J_A(n)
  type(wfun),        intent(out)   :: Jout_A(n)

  call counter_A_Q(ctQA, 0, ntry, J_A, mom, Jout_A, n)

end subroutine counter_A_Q_orig

! **********************************************************************
subroutine counter_A_Q_pid(ctQA, pid, ntry, J_A, mom, Jout_A, n)
! A -> A counter term without left/right splitting
! **********************************************************************
  use KIND_TYPES, only: DREALKIND, REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_qcd_offshell_selfenergies/**/REALKIND, only: quark_ofsse
  use ol_momenta_decl_/**/REALKIND, only: softconf
  use ol_parameters_decl_/**/DREALKIND, only: bubble_vertex, ir_hacks
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_debug, only: ol_fatal
  use ol_kinematics_/**/REALKIND, only: get_LC_5, get_LC_mass2
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer          , intent(in)    :: pid,mom
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(2)
  type(wfun),        intent(inout) :: J_A(n)
  type(wfun),        intent(out)   :: Jout_A(n)
  complex(REALKIND) :: GG(2),K(5)
  integer :: h,momh
  logical :: use_soft_rel

  K = get_LC_5(mom)
  GG = 0
  if (bubble_vertex .eq. 1) GG = quark_ofsse(get_LC_mass2(mom),pid)

  use_soft_rel = .false.
  if (ir_hacks .and. softconf .ne. 0) then
    momh = mom-softconf
    ! check if a soft particle is attached to an external particle
    if (iand(mom,softconf) .ne. 0 .and. iand(mom-softconf,mom-softconf-1) .eq. 0) then
      ! currently only supported for massless particles
      if (get_LC_mass2(momh) .eq. 0 .and. associated(J_A(1)%j_prev)) then
      !if (get_LC_mass2(momh) .eq. 0) then
        use_soft_rel = .true.
      end if
    end if
  end if

  if (use_soft_rel) then
    do h = 1, n
      Jout_A(h)%j(1) = (ctQA(1) + GG(1))*K(5) * J_A(h)%j_prev(1) - (ctQA(2)+GG(2)) * J_A(h)%j(1)
      Jout_A(h)%j(2) = (ctQA(1) + GG(1))*K(5) * J_A(h)%j_prev(2) - (ctQA(2)+GG(2)) * J_A(h)%j(2)
      Jout_A(h)%j(3) = (ctQA(1) + GG(1))*K(5) * J_A(h)%j_prev(3) - (ctQA(2)+GG(2)) * J_A(h)%j(3)
      Jout_A(h)%j(4) = (ctQA(1) + GG(1))*K(5) * J_A(h)%j_prev(4) - (ctQA(2)+GG(2)) * J_A(h)%j(4)

      Jout_A(h)%h = B11 ! default value for the label h

      deallocate(J_A(h)%j_prev)
    end do

  else
    do h = 1, n
      Jout_A(h)%j(1) = (ctQA(1) + GG(1)) * ( + K(1)*J_A(h)%j(3) + K(3)*J_A(h)%j(4)) - (ctQA(2)+GG(2)) * J_A(h)%j(1)
      Jout_A(h)%j(2) = (ctQA(1) + GG(1)) * ( + K(2)*J_A(h)%j(4) + K(4)*J_A(h)%j(3)) - (ctQA(2)+GG(2)) * J_A(h)%j(2)
      Jout_A(h)%j(3) = (ctQA(1) + GG(1)) * ( + K(2)*J_A(h)%j(1) - K(3)*J_A(h)%j(2)) - (ctQA(2)+GG(2)) * J_A(h)%j(3)
      Jout_A(h)%j(4) = (ctQA(1) + GG(1)) * ( + K(1)*J_A(h)%j(2) - K(4)*J_A(h)%j(1)) - (ctQA(2)+GG(2)) * J_A(h)%j(4)

      Jout_A(h)%h = B11 ! default value for the label h

    end do
  end if

  if (ntry == 1) then
    Jout_A(:)%t = J_A(1)%t
    Jout_A(:)%n_part = J_A(1)%n_part
    Jout_A(:)%hf = J_A(:)%hf
    call helbookkeeping_prop(ntry, J_A, Jout_A, n)
  end if

end subroutine counter_A_Q_pid

! **********************************************************************
subroutine counter_Q_A_LR(ctQA, ntry, J_Q, mom, Jout_Q, n)
! Q -> Q counter term with left/right splitting
! ctQA(1)*slash(k)*P_R+ctQA(2)*slash(k)*P_L - ctQA(3) P_R - ctQA(4) P_L
!************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(4)
  integer,           intent(in)    :: mom
  type(wfun),        intent(in)    :: J_Q(n)
  type(wfun),        intent(out)   :: Jout_Q(n)
  complex(REALKIND) :: K(4)
  integer :: h

  K = get_LC_4(mom)
  do h = 1, n
    Jout_Q(h)%j(1) = ctQA(2) * ( - K(2)*J_Q(h)%j(3) + K(4)*J_Q(h)%j(4)) - ctQA(3) * J_Q(h)%j(1)
    Jout_Q(h)%j(2) = ctQA(2) * ( - K(1)*J_Q(h)%j(4) + K(3)*J_Q(h)%j(3)) - ctQA(3) * J_Q(h)%j(2)
    Jout_Q(h)%j(3) = ctQA(1) * ( - K(1)*J_Q(h)%j(1) - K(4)*J_Q(h)%j(2)) - ctQA(4) * J_Q(h)%j(3)
    Jout_Q(h)%j(4) = ctQA(1) * ( - K(2)*J_Q(h)%j(2) - K(3)*J_Q(h)%j(1)) - ctQA(4) * J_Q(h)%j(4)

    Jout_Q(h)%h = B11 ! default value for the label h
  end do

  if (ntry == 1) then
    Jout_Q(:)%t = J_Q(1)%t
    Jout_Q(:)%n_part = J_Q(1)%n_part
    Jout_Q(:)%hf = J_Q(:)%hf
    call helbookkeeping_prop(ntry, J_Q, Jout_Q, n)
  end if

end subroutine counter_Q_A_LR

! **********************************************************************
subroutine counter_A_Q_LR(ctQA, ntry, J_A, mom, Jout_A, n)
! Q -> Q counter term with left/right splitting
! ctQA(1)*slash(k)*P_R+ctQA(2)*slash(k)*P_L - ctQA(3) P_R - ctQA(4) P_L
!************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctQA(4)
  integer,           intent(in)    :: mom
  type(wfun),        intent(in)    :: J_A(n)
  type(wfun),        intent(out)   :: Jout_A(n)
  complex(REALKIND) :: K(4)
  integer :: h

  K = get_LC_4(mom)
  do h = 1, n
    Jout_A(h)%j(1) = ctQA(1) * ( + K(1)*J_A(h)%j(3) + K(3)*J_A(h)%j(4)) - ctQA(3) * J_A(h)%j(1)
    Jout_A(h)%j(2) = ctQA(1) * ( + K(2)*J_A(h)%j(4) + K(4)*J_A(h)%j(3)) - ctQA(3) * J_A(h)%j(2)
    Jout_A(h)%j(3) = ctQA(2) * ( + K(2)*J_A(h)%j(1) - K(3)*J_A(h)%j(2)) - ctQA(4) * J_A(h)%j(3)
    Jout_A(h)%j(4) = ctQA(2) * ( + K(1)*J_A(h)%j(2) - K(4)*J_A(h)%j(1)) - ctQA(4) * J_A(h)%j(4)

    Jout_A(h)%h = B11 ! default value for the label h
  end do

  if (ntry == 1) then
    Jout_A(:)%t = J_A(1)%t
    Jout_A(:)%n_part = J_A(1)%n_part
    Jout_A(:)%hf = J_A(:)%hf
    call helbookkeeping_prop(ntry, J_A, Jout_A, n)
  end if

end subroutine counter_A_Q_LR

! **********************************************************************
subroutine counter_V_V_orig(ctVV, ntry, J_V, momid, Jout_V, n)
  use KIND_TYPES, only: DREALKIND, REALKIND, intkind1, intkind2
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: momid
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctVV(3)
  type(wfun),        intent(in)    :: J_V(n)
  type(wfun),        intent(out)   :: Jout_V(n)
  complex(REALKIND) :: GG(3)
  integer :: h

  call counter_V_V_pid(ctVV, 0, ntry, J_V, momid, Jout_V, n)
end subroutine counter_V_V_orig

! **********************************************************************
subroutine counter_V_V_pid(ctVV, pid, ntry, J_V, momid, Jout_V, n)
! V -> V counter term
! **********************************************************************
  use KIND_TYPES, only: DREALKIND, REALKIND, intkind1, intkind2
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_qcd_offshell_selfenergies/**/REALKIND, only: gluon_ofsse
  use ol_parameters_decl_/**/DREALKIND, only: bubble_vertex
  use ol_kinematics_/**/REALKIND, only: get_LC_mass2, get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: pid,momid
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctVV(3)
  type(wfun),        intent(in)    :: J_V(n)
  type(wfun),        intent(out)   :: Jout_V(n)
  complex(REALKIND) :: GG(3)
  integer :: h

  GG = 0
  if (bubble_vertex .eq. 1) GG = gluon_ofsse(get_LC_mass2(momid),pid)
  do h = 1, n
    Jout_V(h)%j = ((ctVV(1) + GG(1))*get_LC_mass2(momid) + ctVV(2) + GG(2)) * J_V(h)%j + &
                   (ctVV(3) + GG(3)) * (cont_VV(J_V(h)%j,get_LC_4(momid))) * get_LC_4(momid)
  end do

  if (ntry == 1) then
    Jout_V(:)%t       = J_V(1)%t
    Jout_V(:)%n_part  = J_V(1)%n_part
    Jout_V(:)%hf      = J_V(:)%hf
    call helbookkeeping_prop(ntry, J_V, Jout_V, n)
  end if

end subroutine counter_V_V_pid

! **********************************************************************
subroutine counter_S_S(ctSS, ntry, J_S, mom, Jout_S, n)
! ----------------------------------------------------------------------
! S -> S counter term
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop, checkzero_scalar
  use ol_contractions_/**/REALKIND, only: cont_V
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctSS(2)
  integer,           intent(in)    :: mom
  type(wfun),        intent(in)    :: J_S(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  complex(REALKIND) :: K(4)

  K = get_LC_4(mom)
  Jout_S(1:n)%j(1) = (ctSS(1)*cont_V(K) - ctSS(2)) * J_S(1:n)%j(1)

  if (ntry == 1) then
    Jout_S(:)%n_part = J_S(1)%n_part
    Jout_S(:)%t = J_S(1)%t
    Jout_S(:)%hf = J_S(:)%hf
    call checkzero_scalar(Jout_S)
    call helbookkeeping_prop(ntry, J_S, Jout_S, n)
  end if

end subroutine counter_S_S

! **********************************************************************
subroutine counter_S_V(ctSV,ntry,J_S,mom,Jout_V,n)
! ----------------------------------------------------------------------
! S -> V counter term
! i k
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctSV
  integer,           intent(in)    :: mom
  type(wfun),        intent(in)    :: J_S(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  complex(REALKIND) :: K(4)

  K = get_LC_4(mom)
  Jout_V(1:n)%j(1) = -ctSV * K(1) * J_S(1:n)%j(1)
  Jout_V(1:n)%j(2) = -ctSV * K(2) * J_S(1:n)%j(1)
  Jout_V(1:n)%j(3) = -ctSV * K(3) * J_S(1:n)%j(1)
  Jout_V(1:n)%j(4) = -ctSV * K(4) * J_S(1:n)%j(1)

  if (ntry == 1) then
    Jout_V(:)%n_part = J_S(1)%n_part
    Jout_V(:)%t = J_S(1)%t
    Jout_V(:)%hf = J_S(:)%hf
    call helbookkeeping_prop(ntry, J_S, Jout_V, n)
  end if

end subroutine counter_S_V

! **********************************************************************
subroutine counter_V_S(ctSV,ntry,J_V,mom,Jout_S,n)
! ----------------------------------------------------------------------
! S -> V counter term
! -i k
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_V
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_prop
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  complex(REALKIND), intent(in)    :: ctSV
  integer,           intent(in)    :: mom
  type(wfun),        intent(in)    :: J_V(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  complex(REALKIND) :: K(4)

  K = get_LC_4(mom)
  Jout_S(1:n)%j(1) =  ctSV * cont_V(K)

  if (ntry == 1) then
    Jout_S(:)%n_part = J_V(1)%n_part
    Jout_S(:)%t = J_V(1)%t
    Jout_S(:)%hf = J_V(:)%hf
    call helbookkeeping_prop(ntry, J_V, Jout_S, n)
  end if

end subroutine counter_V_S

! ======================================================================
! Vertex counter terms.
! ======================================================================

! **********************************************************************
subroutine counter_GGG_G(ntry, J_G1, J_G2, J_G3, Jout_G, n, t)
! gluon-gluon-gluon-gluon vertex for R2, factorised Lorentz monomials g(a,b)*g(c,d)
! Jout_G(d) = g(a,b)*g(c,d) * J_G1(a) * J_G2(b) * J_G3(c) = J_G1.J_G2 * J_G3
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none

  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_G1(n(1)), J_G2(n(2)), J_G3(n(3))
  type(wfun), intent(out) :: Jout_G(:)

  integer :: h

  do h = 1, n(4)
    Jout_G(h)%j = cont_VV(J_G1(t(1,h))%j,J_G2(t(2,h))%j) * J_G3(t(3,h))%j
  end do

  if (ntry == 1) then
    Jout_G(:)%t = J_G1(1)%t + J_G2(1)%t + J_G3(1)%t
    Jout_G(:)%n_part = J_G1(1)%n_part + J_G2(1)%n_part + J_G3(1)%n_part
    Jout_G(:)%hf = J_G1(t(1,:))%hf + J_G2(t(2,:))%hf + J_G3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_G1, J_G2, J_G3, Jout_G, n, t)
  end if

end subroutine counter_GGG_G

! **********************************************************************
subroutine counter_VQ_A(ntry, J_V, J_Q, Jout_Q, n, t)
! VQ -> Q counter term; without left/right splitting
! Factorised wrt. vert_VQ_A: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: J_V(n(1)), J_Q(n(2))
  type(wfun),        intent(out)   :: Jout_Q(n(3))
  integer :: h

  do h = 1, n(3)

    Jout_Q(h)%hf = J_V(t(1,h))%hf + J_Q(t(2,h))%hf

    Jout_Q(h)%j(1) = - J_V(t(1,h))%j(2)*J_Q(t(2,h))%j(3) + J_V(t(1,h))%j(4)*J_Q(t(2,h))%j(4)
    Jout_Q(h)%j(2) = - J_V(t(1,h))%j(1)*J_Q(t(2,h))%j(4) + J_V(t(1,h))%j(3)*J_Q(t(2,h))%j(3)
    Jout_Q(h)%j(3) = - J_V(t(1,h))%j(1)*J_Q(t(2,h))%j(1) - J_V(t(1,h))%j(4)*J_Q(t(2,h))%j(2)
    Jout_Q(h)%j(4) = - J_V(t(1,h))%j(2)*J_Q(t(2,h))%j(2) - J_V(t(1,h))%j(3)*J_Q(t(2,h))%j(1)

    Jout_Q(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_Q(:)%t = J_V(1)%t + J_Q(1)%t
    Jout_Q(:)%n_part = J_V(1)%n_part + J_Q(1)%n_part
    Jout_Q(:)%hf = J_V(t(1,:))%hf + J_Q(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_V, J_Q, Jout_Q, n, t)
  end if

end subroutine counter_VQ_A

! **********************************************************************
subroutine counter_AV_Q(ntry, J_A, J_V, Jout_A, n, t)
! AV -> A counter term; without left/right splitting
! Factorised wrt. vert_AV_Q: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:), J_V(:)
  type(wfun),        intent(out)   :: Jout_A(:)
  integer :: h

  do h = 1, n(3)

    Jout_A(h)%j(1) = - J_V(t(2,h))%j(1)*J_A(t(1,h))%j(3) - J_V(t(2,h))%j(3)*J_A(t(1,h))%j(4)
    Jout_A(h)%j(2) = - J_V(t(2,h))%j(2)*J_A(t(1,h))%j(4) - J_V(t(2,h))%j(4)*J_A(t(1,h))%j(3)
    Jout_A(h)%j(3) = - J_V(t(2,h))%j(2)*J_A(t(1,h))%j(1) + J_V(t(2,h))%j(3)*J_A(t(1,h))%j(2)
    Jout_A(h)%j(4) = - J_V(t(2,h))%j(1)*J_A(t(1,h))%j(2) + J_V(t(2,h))%j(4)*J_A(t(1,h))%j(1)

    Jout_A(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_A(:)%n_part = J_A(1)%n_part + J_V(1)%n_part
    Jout_A(:)%t = J_A(1)%t + J_V(1)%t
    Jout_A(:)%hf = J_A(t(1,:))%hf + J_V(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_A, J_V, Jout_A, n, t)
  end if

end subroutine counter_AV_Q

! **********************************************************************
subroutine counter_QA_V(ntry, J_Q, J_A, Jout_V, n, t)
! ----------------------------------------------------------------------
! QA -> V counter term; without left/right splitting
! ----------------------------------------------------------------------
! Factorised wrt. vert_QA_V: ctQAV = gQCD*dlnG + dZf + 1/2 * dZg
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun),        intent(in)    :: J_Q(n(1)), J_A(n(2))
  type(wfun),        intent(out)   :: Jout_V(n(3))
  integer :: h

  do h = 1, n(3)

    Jout_V(h)%j(1) = - J_A(t(2,h))%j(1)*J_Q(t(1,h))%j(3) - J_A(t(2,h))%j(4)*J_Q(t(1,h))%j(2)
    Jout_V(h)%j(2) = - J_A(t(2,h))%j(2)*J_Q(t(1,h))%j(4) - J_A(t(2,h))%j(3)*J_Q(t(1,h))%j(1)
    Jout_V(h)%j(3) = - J_A(t(2,h))%j(1)*J_Q(t(1,h))%j(4) + J_A(t(2,h))%j(3)*J_Q(t(1,h))%j(2)
    Jout_V(h)%j(4) = - J_A(t(2,h))%j(2)*J_Q(t(1,h))%j(3) + J_A(t(2,h))%j(4)*J_Q(t(1,h))%j(1)
    Jout_V(h)%j(:) = Jout_V(h)%j(:) + Jout_V(h)%j(:)

    Jout_V(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_V(:)%t = J_Q(1)%t + J_A(1)%t
    Jout_V(:)%n_part = J_Q(1)%t + J_A(1)%n_part
    Jout_V(:)%hf = J_Q(t(1,:))%hf + J_A(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_Q, J_A, Jout_V, n, t)
  end if

end subroutine counter_QA_V

! **********************************************************************
subroutine counter_QA_Z(g_RL, ntry, J_Q, J_A, Jout_Z, n, t)
! ----------------------------------------------------------------------
! bare QA -> Z Z-like interaction
! ----------------------------------------------------------------------
! J_Q(4)    = quark current
! J_A(4)    = anti-quark current
! g_RL(1)   = right-handed coupling gR
! g_RL(2)   = left-handed coupling gL
! Jout_Z(4) = outgoing Z current (light-cone rep.)
! Jout_Z(A) = J_A(i) * [gamma^A*(gR*w_R+gL*w_L)](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none

  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun), intent(in) :: J_Q(n(1)), J_A(n(2))
  type(wfun), intent(out) :: Jout_Z(n(3))
  complex(REALKIND), intent(in) :: g_RL(2)
  integer :: h

  do h = 1, n(3)
    Jout_Z(h)%j(1) = - g_RL(2)*J_A(t(2,h))%j(1)*J_Q(t(1,h))%j(3) - &
    g_RL(1)*J_A(t(2,h))%j(4)*J_Q(t(1,h))%j(2)

    Jout_Z(h)%j(2) = - g_RL(2)*J_A(t(2,h))%j(2)*J_Q(t(1,h))%j(4) - &
    g_RL(1)*J_A(t(2,h))%j(3)*J_Q(t(1,h))%j(1)

    Jout_Z(h)%j(3) = - g_RL(2)*J_A(t(2,h))%j(1)*J_Q(t(1,h))%j(4) + &
    g_RL(1)*J_A(t(2,h))%j(3)*J_Q(t(1,h))%j(2)

    Jout_Z(h)%j(4) = - g_RL(2)*J_A(t(2,h))%j(2)*J_Q(t(1,h))%j(3) + &
    g_RL(1)*J_A(t(2,h))%j(4)*J_Q(t(1,h))%j(1)

    Jout_Z(h)%j(:) = Jout_Z(h)%j(:) + Jout_Z(h)%j(:)

    Jout_Z(h)%h = B11 ! default value for the label h
  end do

  if (ntry == 1) then
    Jout_Z(:)%t = J_Q(1)%t + J_A(1)%t
    Jout_Z(:)%n_part = J_Q(1)%t + J_A(1)%n_part
    Jout_Z(:)%hf = J_Q(t(1,:))%hf + J_A(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_Q, J_A, Jout_Z, n, t)
  end if

end subroutine counter_QA_Z

! **********************************************************************
subroutine counter_AZ_Q(g_RL, ntry, J_A, J_Z, Jout_A, n, t)
! ----------------------------------------------------------------------
! bare AZ -> A Z-like interaction
! ----------------------------------------------------------------------
! J_A(4)     = incoming anti-quark current
! J_Z(4)     = incoming Z current (light-cone rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_A(4)  = outgoing anti-quark current
! Jout_A(i)  = J_A(j) * [gamma_A*(gR*w_R+gL*w_L)](j,i) * J_Z(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:), J_Z(:)
  type(wfun),        intent(out)   :: Jout_A(:)
  complex(REALKIND) :: g_RL(2)
  integer :: h
  integer(intkind2) :: h1, h2

  do h = 1, n(3)
    h1 = t(1,h)
    h2 = t(2,h)
    Jout_A(h)%j(1) = g_RL(1) * ( - J_Z(h2)%j(1)*J_A(h1)%j(3) - J_Z(h2)%j(3)*J_A(h1)%j(4))
    Jout_A(h)%j(2) = g_RL(1) * ( - J_Z(h2)%j(2)*J_A(h1)%j(4) - J_Z(h2)%j(4)*J_A(h1)%j(3))
    Jout_A(h)%j(3) = g_RL(2) * ( - J_Z(h2)%j(2)*J_A(h1)%j(1) + J_Z(h2)%j(3)*J_A(h1)%j(2))
    Jout_A(h)%j(4) = g_RL(2) * ( - J_Z(h2)%j(1)*J_A(h1)%j(2) + J_Z(h2)%j(4)*J_A(h1)%j(1))

    Jout_A(h)%h = B11 ! default value for the label h
  end do

  if (ntry == 1) then
    Jout_A(:)%n_part = J_A(1)%n_part + J_Z(1)%n_part
    Jout_A(:)%t = J_A(1)%t + J_Z(1)%t
    Jout_A(:)%hf = J_A(t(1,:))%hf + J_Z(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_A, J_Z, Jout_A, n, t)
  end if

end subroutine counter_AZ_Q

! **********************************************************************
subroutine counter_ZQ_A(g_RL, ntry, J_Z, J_Q, Jout_Q, n, t)
! ----------------------------------------------------------------------
! bare ZQ -> Q Z-like interaction
! ----------------------------------------------------------------------
! J_Q(:)     = incoming quark current
! J_Z(:)     = incoming Z current ("light-cone" rep.)
! g_RL(1)    = right-handed coupling gR
! g_RL(2)    = left-handed coupling gL
! Jout_Q(:)  = outgoing quark current
! Jout_Q(i)  = J_Z(A)*[gamma_A*(gR*w_R+gL*w_L)](i,j)*J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Z(:), J_Q(:)
  type(wfun),        intent(out)   :: Jout_Q(:)
  complex(REALKIND) :: g_RL(2)
  integer :: h
  integer(intkind2) :: h1, h2

  do h = 1, n(3)
    h1 = t(1,h)
    h2 = t(2,h)
    Jout_Q(h)%j(1) = g_RL(2) * ( - J_Z(h1)%j(2)*J_Q(h2)%j(3) + J_Z(h1)%j(4)*J_Q(h2)%j(4))
    Jout_Q(h)%j(2) = g_RL(2) * ( - J_Z(h1)%j(1)*J_Q(h2)%j(4) + J_Z(h1)%j(3)*J_Q(h2)%j(3))
    Jout_Q(h)%j(3) = g_RL(1) * ( - J_Z(h1)%j(1)*J_Q(h2)%j(1) - J_Z(h1)%j(4)*J_Q(h2)%j(2))
    Jout_Q(h)%j(4) = g_RL(1) * ( - J_Z(h1)%j(2)*J_Q(h2)%j(2) - J_Z(h1)%j(3)*J_Q(h2)%j(1))

    Jout_Q(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_Q(:)%n_part = J_Z(1)%n_part + J_Q(1)%n_part
    Jout_Q(:)%t = J_Z(1)%t + J_Q(1)%t
    Jout_Q(:)%hf = J_Z(t(1,:))%hf + J_Q(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_Z, J_Q, Jout_Q, n, t)
  end if

end subroutine counter_ZQ_A

! **********************************************************************
subroutine counter_UV_W(ntry, J_V1, mom1id, J_V2, mom2id, Jout_V, n, t)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer,           intent(in)    :: mom1id,mom2id
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(2,n(3))
  type(wfun) :: J_V1(n(1)), J_V2(n(2))
  type(wfun) :: Jout_V(n(3))
  complex(REALKIND) :: J1J2(n(3)), P1J2(n(2)), P2J1(n(1))
  complex(REALKIND) :: P1(4), P2(4)
  integer :: h

  P1 = get_LC_4(mom1id)
  P2 = get_LC_4(mom2id)

  Jout_V(:)%t = J_V1(1)%t + J_V2(1)%t
  Jout_V(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part

  do h = 1, n(3)

    Jout_V(h)%hf = J_V1(t(1,h))%hf + J_V2(t(2,h))%hf

    J1J2(h) = cont_VV(J_V1(t(1,h))%j, J_V2(t(2,h))%j )
    P1J2(t(2,h)) = cont_VV(P1+P1+P2,J_V2(t(2,h))%j)
    P2J1(t(1,h)) = cont_VV(P1+P2+P2,J_V1(t(1,h))%j)
    Jout_V(h)%j = J1J2(h) * (P1 - P2) + P2J1(t(1,h)) * J_V2(t(2,h))%j - P1J2(t(2,h)) * J_V1(t(1,h))%j

  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part
    Jout_V(:)%t = J_V1(1)%t + J_V2(1)%t
    Jout_V(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_V1, J_V2, Jout_V, n, t)
  end if

end subroutine counter_UV_W

! **********************************************************************
subroutine counter_EV_V(ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun) :: J_V1(n(1)), J_V2(n(2)), J_V3(n(3))
  type(wfun) :: Jout_V(n(4))
  integer :: h

  do h = 1, n(3)
    Jout_V(h)%j(:) = cont_VV(J_V1(t(1,h))%j(1:4),J_V3(t(3,h))%j(1:4)) * J_V2(t(2,h))%j(:) -&
                     cont_VV(J_V2(t(2,h))%j(1:4),J_V3(t(3,h))%j(1:4)) * J_V1(t(1,h))%j(:)
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part + J_V3(1)%n_part
    Jout_V(:)%t = J_V1(1)%t + J_V2(1)%t + J_V3(1)%t
    Jout_V(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf + J_V3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
  end if

end subroutine counter_EV_V


! **********************************************************************
subroutine counter_AW_Q(ntry, J_A, J_W, Jout_A, n, t)
! ----------------------------------------------------------------------
! bare AW -> A W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_A(:)    = incoming anti-quark current
! J_W(:)    = incoming W current (light-cone rep.)
! Jout_A(:) = outgoing anti-quark current
! Jout_A(:) = J_A(j) * [gamma_A*w_L](j,i) * J_W(A)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:), J_W(:)
  type(wfun),        intent(out)   :: Jout_A(:)
  integer :: h

  do h = 1, n(3)
    Jout_A(h)%j(1:2) = 0
    Jout_A(h)%j(3)= - J_W(t(2,h))%j(2)*J_A(t(1,h))%j(1) + J_W(t(2,h))%j(3)*J_A(t(1,h))%j(2)
    Jout_A(h)%j(4)= - J_W(t(2,h))%j(1)*J_A(t(1,h))%j(2) + J_W(t(2,h))%j(4)*J_A(t(1,h))%j(1)

    Jout_A(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_A(:)%n_part = J_A(1)%n_part + J_W(1)%n_part
    Jout_A(:)%t = J_A(1)%t + J_W(1)%t
    Jout_A(:)%hf = J_A(t(1,:))%hf + J_W(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_A, J_W, Jout_A, n, t)
  end if

end subroutine counter_AW_Q

! **********************************************************************
subroutine counter_WQ_A(ntry, J_W, J_Q, Jout_Q, n, t)
! ----------------------------------------------------------------------
! bare WQ -> Q W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(:)    = incoming quark current
! J_W(:)    = incoming W current ("light-cone" rep.)
! Jout_Q(:) = outgoing quark current
! Jout_Q(:) = J_W(A) * [gamma_A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_W(:), J_Q(:)
  type(wfun),        intent(out)   :: Jout_Q(:)
  integer :: h

  do h = 1, n(3)
    Jout_Q(h)%j(1)= - J_W(t(1,h))%j(2)*J_Q(t(2,h))%j(3) + J_W(t(1,h))%j(4)*J_Q(t(2,h))%j(4)
    Jout_Q(h)%j(2)= - J_W(t(1,h))%j(1)*J_Q(t(2,h))%j(4) + J_W(t(1,h))%j(3)*J_Q(t(2,h))%j(3)
    Jout_Q(h)%j(3:4) = 0

    Jout_Q(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_Q(:)%n_part = J_Q(1)%n_part + J_W(1)%n_part
    Jout_Q(:)%t = J_Q(1)%t + J_W(1)%t
    Jout_Q(:)%hf = J_W(t(1,:))%hf + J_Q(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_W, J_Q, Jout_Q, n, t)
  end if

end subroutine counter_WQ_A

! **********************************************************************
subroutine counter_QA_W(ntry, J_Q, J_A, Jout_W, n, t)
! ----------------------------------------------------------------------
! bare QA -> W W-like (i.e. left-handed) interaction
! ----------------------------------------------------------------------
! J_Q(:)    = incoming quark current
! J_W(:)    = incoming W current ("light-cone" rep.)
! Jout_Q(:) = outgoing quark current
! Jout_Q(:) = J_W(A) * [gamma_A*w_L](i,j) * J_Q(j)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:), J_A(:)
  type(wfun),        intent(out)   :: Jout_W(:)
  integer :: h, h1, h2

  do h = 1, n(3)
    h1 = t(1,h)
    h2 = t(2,h)
    Jout_W(h)%j(1) = - J_A(h2)%j(1)*J_Q(h1)%j(3)
    Jout_W(h)%j(2) = - J_A(h2)%j(2)*J_Q(h1)%j(4)
    Jout_W(h)%j(3) = - J_A(h2)%j(1)*J_Q(h1)%j(4)
    Jout_W(h)%j(4) = - J_A(h2)%j(2)*J_Q(h1)%j(3)
    Jout_W(h)%j = Jout_W(h)%j + Jout_W(h)%j

    Jout_W(h)%h = B11 ! default value for the label h

  end do

  if (ntry == 1) then
    Jout_W(:)%n_part = J_Q(1)%n_part + J_A(1)%n_part
    Jout_W(:)%t = J_Q(1)%t + J_A(1)%t
    Jout_W(:)%hf = J_Q(t(1,:))%hf + J_A(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_Q, J_A, Jout_W, n, t)
  end if

end subroutine counter_QA_W


! **********************************************************************
subroutine counter_SG_G(ntry, J_S, Jin_G, Jout_G, n, t)
! ----------------------------------------------------------------------
! Higgs-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_G = J_S * Jin_G
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:), Jin_G(:)
  type(wfun),        intent(out)   :: Jout_G(:)
  integer :: h

  do h = 1, n(3)
    Jout_G(h)%j(:) = J_S(t(1,h))%j(1) * Jin_G(t(2,h))%j(:)
  end do

  if (ntry == 1) then
    Jout_G(:)%n_part = J_S(1)%n_part + Jin_G(1)%n_part
    Jout_G(:)%t = J_S(1)%t + Jin_G(1)%t
    Jout_G(:)%hf = J_S(t(1,:))%hf + Jin_G(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S, Jin_G, Jout_G, n, t)
  end if

end subroutine counter_SG_G


! **********************************************************************
subroutine counter_GG_S(ntry, J_G1, J_G2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Higgs-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_S = J_G1.J_G2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_G1(:), J_G2(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  integer :: h

  do h = 1, n(3)
    Jout_S(h)%j(1) = cont_VV(J_G1(t(1,h))%j,J_G2(t(2,h))%j)
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_G1(1)%n_part + J_G2(1)%n_part
    Jout_S(:)%t = J_G1(1)%t + J_G2(1)%t
    Jout_S(:)%hf = J_G1(t(1,:))%hf + J_G2(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_G1, J_G2, Jout_S, n, t)
  end if

end subroutine counter_GG_S

! **********************************************************************
subroutine counter_GG_V(ntry, J_G1, mom1, J_G2, mom2, Jout_V, n, t)
! ----------------------------------------------------------------------
! Z-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_V(a) = ep(a,b,c,d) * J_G1(b) * J_G2(c) * (p1-p2)(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_EpVVV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom1,mom2
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_G1(:), J_G2(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  complex(REALKIND) :: p1(4), p2(4)
  integer :: h

  p1 = get_LC_4(mom1)
  p2 = get_LC_4(mom2)
  do h=1, n(3)
    call cont_EpVVV(J_G1(t(1,h))%j, J_G2(t(2,h))%j, p1-p2, Jout_V(h)%j)
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_G1(1)%n_part + J_G2(1)%n_part
    Jout_V(:)%t = J_G1(1)%t + J_G2(1)%t
    Jout_V(:)%hf = J_G1(t(1,:))%hf + J_G2(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_G1, J_G2, Jout_V, n, t)
  end if

end subroutine counter_GG_V

! **********************************************************************
subroutine counter_VG_G(ntry, J_V, Jin_G, mom2, Jout_G, mom3, n, t)
! ----------------------------------------------------------------------
! Z-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_G(c) = ep(a,b,c,d) * J_V(a) * Jin_G(b) * (p2-p3)(d)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_EpVVV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: mom2,mom3
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)  :: J_V(:), Jin_G(:)
  type(wfun),        intent(out) :: Jout_G(:)
  complex(REALKIND) :: p23(4)
  integer           :: h

  p23 = -get_LC_4(mom3-mom2)
  do h = 1, n(3)
    call cont_EpVVV(J_V(t(1,h))%j, Jin_G(t(2,h))%j, p23, Jout_G(h)%j)
  end do

   if (ntry == 1) then
    Jout_G(:)%n_part = J_V(1)%n_part + Jin_G(1)%n_part
    Jout_G(:)%t = J_V(1)%t + Jin_G(1)%t
    Jout_G(:)%hf = J_V(t(1,:))%hf + Jin_G(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_V, Jin_G, Jout_G, n, t)
  end if

end subroutine counter_VG_G


! **********************************************************************
subroutine counter_VVG_G(ntry, J_V1, J_V2, Jin_G, Jout_G, n, t)
! ----------------------------------------------------------------------
! Vector-vector-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_G = J_V1.J_V2 * Jin_G + J_V1.Jin_G * J_V2 + J_V2.Jin_G * J_V1
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),  intent(in) :: J_V1(:), J_V2(:), Jin_G(:)
  type(wfun),  intent(out) :: Jout_G(:)
  integer :: h
  integer(intkind2) :: h1, h2, h3

  do h = 1, n(4)
    h1 = t(1,h)
    h2 = t(2,h)
    h3 = t(3,h)
    Jout_G(h)%j(1:4) = cont_VV(J_V1(h1)%j,J_V2(h2)%j) * Jin_G(h3)%j(1:4) +&
    & cont_VV(J_V1(h1)%j,Jin_G(h3)%j) * J_V2(h2)%j(1:4) +&
    & cont_VV(J_V2(h2)%j,Jin_G(h3)%j) * J_V1(h1)%j(1:4)
  ! Jout_G = Jout_G + 4*(1-g5s)*sum_ij(a_V1QiQj*aV2QiQj) * cont_VV(J_V1,J_V2) * Jin_G
  end do

  if (ntry == 1) then
    Jout_G(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part + Jin_G(1)%n_part
    Jout_G(:)%t = J_V1(1)%t + J_V2(1)%t + Jin_G(1)%t
    Jout_G(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf + Jin_G(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_V1, J_V2, Jin_G, Jout_G, n, t)
  end if

end subroutine counter_VVG_G

! **********************************************************************
subroutine counter_SSG_G(ntry, J_S1, J_S2, Jin_G, Jout_G, n, t)
! ----------------------------------------------------------------------
! Scalar-scalar-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_G = J_S1 * J_S2 * Jin_G
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),  intent(in) :: J_S1(:), J_S2(:), Jin_G(:)
  type(wfun),  intent(out) :: Jout_G(:)
  integer :: h

  do h = 1, n(4)
    Jout_G(h)%j(:) = J_S1(t(1,h))%j(1) * J_S2(t(2,h))%j(1) * Jin_G(t(3,h))%j(:)
  end do

  if (ntry == 1) then
    Jout_G(:)%n_part = J_S1(1)%n_part + J_S2(1)%n_part + Jin_G(1)%n_part
    Jout_G(:)%t = J_S1(1)%t + J_S2(1)%t + Jin_G(1)%t
    Jout_G(:)%hf = J_S1(t(1,:))%hf + J_S2(t(2,:))%hf + Jin_G(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_S1, J_S2, Jin_G, Jout_G, n, t)
  end if

end subroutine counter_SSG_G

! **********************************************************************
subroutine counter_GGS_S(ntry, J_G1, J_G2, Jin_S, Jout_S, n, t)
! ----------------------------------------------------------------------
! Scalar-scalar-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_S = Jin_S * J_G1.J_G2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4, checkzero_scalar
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),  intent(in) :: J_G1(:), J_G2(:), Jin_S(:)
  type(wfun),  intent(out) :: Jout_S(:)
  integer :: h

  do h = 1, n(4)
    Jout_S(h)%j(1) = Jin_S(t(3,h))%j(1) * cont_VV(J_G1(t(1,h))%j(1:4),J_G2(t(2,h))%j(1:4))
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_G1(1)%n_part + J_G2(1)%n_part + Jin_S(1)%n_part
    Jout_S(:)%t = J_G1(1)%t + J_G2(1)%t + Jin_S(1)%t
    Jout_S(:)%hf = J_G1(t(1,:))%hf + J_G2(t(2,:))%hf + Jin_S(t(3,:))%hf
    call checkzero_scalar(Jout_S)
    call helbookkeeping_vert4(ntry, J_G1, J_G2, Jin_S, Jout_S, n, t)
  end if

end subroutine counter_GGS_S

! **********************************************************************
subroutine counter_GGG_V(gVA, ntry, J_G1, J_G2, J_G3, Jout_V, n, t)
! ----------------------------------------------------------------------
! Vector-gluon-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_V(a) = [v * (g(a,b)*g(c,d) + g(a,c)*g(b,d) + g(a,d)*g(b,c)) + a * ep(a,b,c,d)] * J_G1(b) * J_G2(c) * J_G3(d)
! v = 2/3*vector_coupling; a = -6*i*axial_coupling
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_EpVVV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  complex(REALKIND), intent(in)  :: gVA(2)
  type(wfun),  intent(in) :: J_G1(:), J_G2(:), J_G3(:)
  type(wfun),  intent(out) :: Jout_V(:)
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  integer(intkind2) :: h, h1, h2, h3

  do h = 1, n(4)
    h1 = t(1,h)
    h2 = t(2,h)
    h3 = t(3,h)
    if (gVA(2) /= 0) then
      call cont_EpVVV(J_G1(h1)%j, J_G2(h2)%j, J_G3(h3)%j, Jout_V(h)%j)
      Jout_V(h)%j(1:4) = gVA(2) * Jout_V(h)%j(1:4)
    else
      Jout_V(h)%j = 0
    end if

    Jout_V(h)%j(1:4) = Jout_V(h)%j(1:4) + &
      gVA(1) * (cont_VV(J_G1(h1)%j,J_G2(h2)%j)*J_G3(h3)%j(1:4) + &
      cont_VV(J_G2(h2)%j,J_G3(h3)%j)*J_G1(h1)%j(1:4) + &
      cont_VV(J_G3(h3)%j,J_G1(h1)%j)*J_G2(h2)%j(1:4))
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_G1(1)%n_part + J_G2(1)%n_part + J_G3(1)%n_part
    Jout_V(:)%t = J_G1(1)%t + J_G2(1)%t + J_G3(1)%t
    Jout_V(:)%hf = J_G1(t(1,:))%hf + J_G2(t(2,:))%hf + J_G3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_G1, J_G2, J_G3, Jout_V, n, t)
  end if

end subroutine counter_GGG_V

! **********************************************************************
subroutine counter_VGG_G(gVA, ntry, J_V, J_G1, J_G2, Jout_G, n, t)
! ----------------------------------------------------------------------
! Vector-gluon-gluon-gluon vertex for R2
! ----------------------------------------------------------------------
! Jout_G(d) = [v * (g(a,b)*g(c,d) + g(a,c)*g(b,d) + g(a,d)*g(b,c)) + a * ep(a,b,c,d)] * J_V(a) * J_G1(b) * J_G2(c)
! v = 2/3*vector_coupling; a = -6*i*axial_coupling
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_VV, cont_EpVVV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in) :: ntry
  complex(REALKIND), intent(in) :: gVA(2)
  type(wfun), intent(in)  :: J_V(:), J_G1(:), J_G2(:)
  type(wfun), intent(out) :: Jout_G(:)
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  integer(intkind2) :: h
  integer(intkind2) :: h1, h2, h3


  do h = 1, n(4)
    h1 = t(1,h)
    h2 = t(2,h)
    h3 = t(3,h)
    if (gVA(2) /= 0) then
      call cont_EpVVV(J_V(h1)%j, J_G1(h2)%j, J_G2(h3)%j, Jout_G(h)%j)
      Jout_G(h)%j(1:4) = - gVA(2) * Jout_G(h)%j(1:4)
    else
      Jout_G(h)%j = 0
    end if
    Jout_G(h)%j(1:4) = Jout_G(h)%j(1:4) + &
      gVA(1) * (cont_VV(J_V(h1)%j,J_G1(h2)%j)*J_G2(h3)%j(1:4) + &
                cont_VV(J_G1(h2)%j,J_G2(h3)%j)*J_V(h1)%j(1:4) + &
                cont_VV(J_G2(h3)%j,J_V(h1)%j) *J_G1(h2)%j(1:4))
  end do

  if (ntry == 1) then
   Jout_G(:)%n_part = J_V(1)%n_part + J_G1(1)%n_part + J_G2(1)%n_part
   Jout_G(:)%t = J_V(1)%t + J_G1(1)%t + J_G2(1)%t
   Jout_G(:)%hf = J_V(t(1,:))%hf + J_G1(t(2,:))%hf + J_G2(t(3,:))%hf
   call helbookkeeping_vert4(ntry, J_V, J_G1, J_G2, Jout_G, n, t)
 end if

end subroutine counter_VGG_G


! **********************************************************************
subroutine counter_QS_A(g_RL, ntry, J_Q, J_S, Jout_A, n, t)
! ----------------------------------------------------------------------
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming fermion current:      J_Q
! Incoming scalar current:       J_S
! Outgoing anti-fermion current: Jout_A(4)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:), J_S(:)
  type(wfun),        intent(out)   :: Jout_A(:)
  complex(REALKIND) :: g_RL(2)
  integer :: h

  do h = 1, n(3)
    Jout_A(h)%j(1:2) = g_RL(1) * J_S(t(2,h))%j(1) * J_Q(t(1,h))%j(1:2)
    Jout_A(h)%j(3:4) = g_RL(2) * J_S(t(2,h))%j(1) * J_Q(t(1,h))%j(3:4)

    Jout_A(h)%h = B11
  end do

  if (ntry == 1) then
    Jout_A(:)%n_part = J_Q(1)%n_part + J_S(1)%n_part
    Jout_A(:)%t = J_Q(1)%t + J_S(1)%t
    Jout_A(:)%hf = J_Q(t(1,:))%hf + J_S(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_Q, J_S, Jout_A, n, t)
  end if

end subroutine counter_QS_A

! **********************************************************************
subroutine counter_AQ_S(g_RL, ntry, J_A, J_Q, Jout_S, n, t)
! ----------------------------------------------------------------------
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming anti-fermion current: J_A
! Incoming fermion current:      J_Q
! Outgoing scalar current:       Jout_S = gR*J_A.P_R.J_Q + gL*J_A.P_L.J_Q
!   with the right- and left-handed projectors P_R = (1+y5)/2 and P_L = (1-y5)/2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:), J_Q(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  complex(REALKIND), intent(in) :: g_RL(2)
  integer :: h

  do h = 1, n(3)
    Jout_S(h)%j(1) = g_RL(1) * (J_A(t(1,h))%j(1)*J_Q(t(2,h))%j(1) + J_A(t(1,h))%j(2)*J_Q(t(2,h))%j(2)) +&
    g_RL(2) * (J_A(t(1,h))%j(3)*J_Q(t(2,h))%j(3) + J_A(t(1,h))%j(4)*J_Q(t(2,h))%j(4))
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_A(1)%n_part + J_Q(1)%n_part
    Jout_S(:)%t = J_A(1)%t + J_Q(1)%t
    Jout_S(:)%hf = J_A(t(1,:))%hf + J_Q(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_A, J_Q, Jout_S, n, t)
  end if

end subroutine counter_AQ_S

! **********************************************************************
subroutine counter_SA_Q(g_RL, ntry, J_S, J_A, Jout_Q, n, t)
! ----------------------------------------------------------------------
! Fermion-scalar-vertex
! ----------------------------------------------------------------------
! g_RL(1) = right-handed coupling
! g_RL(2) = left-handed coupling
! Incoming scalar current:       J_S
! Incoming anti-fermion current: J_A
! Outgoing fermion current:      Jout_Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use KIND_TYPES_BW
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:), J_A(:)
  type(wfun),        intent(out)   :: Jout_Q(:)
  complex(REALKIND), intent(in) :: g_RL(2)
  integer:: h

  do h=1, n(3)
    Jout_Q(h)%j(1:2) = g_RL(1) * J_S(t(1,h))%j(1) * J_A(t(2,h))%j(1:2)
    Jout_Q(h)%j(3:4) = g_RL(2) * J_S(t(1,h))%j(1) * J_A(t(2,h))%j(3:4)

    Jout_Q(h)%h = B11
  end do

  if (ntry == 1) then
    Jout_Q(:)%n_part = J_S(1)%n_part + J_A(1)%n_part
    Jout_Q(:)%t = J_S(1)%t + J_A(1)%t
    Jout_Q(:)%hf = J_S(t(1,:))%hf + J_A(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S, J_A, Jout_Q, n, t)
  end if

end subroutine counter_SA_Q

!******************************************************************************!
!                               EW Counterterms                                !
!******************************************************************************!


! **********************************************************************
subroutine counter_VQ_A_LR(ctVFF,ntry,J_V,J_Q,Jout_Q,n,t)
! ----------------------------------------------------------------------
! VQ -> Q counter term; with left/right splitting -> ZQ_A
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  complex(REALKIND), intent(in)    :: ctVFF(2)
  type(wfun),        intent(in)    :: J_V(:), J_Q(:)
  type(wfun),        intent(out)   :: Jout_Q(:)
  call counter_ZQ_A(ctVFF,ntry,J_V,J_Q,Jout_Q,n,t)
end subroutine counter_VQ_A_LR


! **********************************************************************
subroutine counter_AV_Q_LR(ctVFF,ntry,J_A,J_V,Jout_A,n,t)
! ----------------------------------------------------------------------
! AV -> A counter term; with left/right splitting -> AZ_Q
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  complex(REALKIND), intent(in)    :: ctVFF(2)
  type(wfun),        intent(in)    :: J_A(:), J_V(:)
  type(wfun),        intent(out)   :: Jout_A(:)
  call counter_AZ_Q(ctVFF,ntry,J_A,J_V,Jout_A,n,t)
end subroutine counter_AV_Q_LR


! **********************************************************************
subroutine counter_QA_V_LR(ctVFF,ntry,J_Q,J_A,Jout_V,n,t)
! ----------------------------------------------------------------------
! QA -> V counter term; with left/right splitting -> QA_Z
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  complex(REALKIND), intent(in)    :: ctVFF(2)
  type(wfun),        intent(in)    :: J_Q(:), J_A(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  call counter_QA_Z(ctVFF,ntry,J_Q,J_A,Jout_V,n,t)
end subroutine counter_QA_V_LR


! **********************************************************************
subroutine counter_SS_S(ntry, J_S1, J_S2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Three scalar vertex
! ----------------------------------------------------------------------
! Incoming scalar currents: J_S1, J_S2
! Outgoing scalar current:  Jout_S = J_S1 * J_S2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  integer:: h

  do h=1, n(3)
    Jout_S(h)%j(1) = J_S1(t(1,h))%j(1) * J_S2(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_S1(1)%n_part + J_S2(1)%n_part
    Jout_S(:)%t = J_S1(1)%t + J_S2(1)%t
    Jout_S(:)%hf = J_S1(t(1,:))%hf + J_S2(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S1, J_S2, Jout_S, n, t)
  end if

end subroutine counter_SS_S

! **********************************************************************
subroutine counter_VS_T(ntry, J_V, mom1, J_S, mom2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Incoming vector current: J_V(4), incoming momentum P1(4) (light-cone rep.)
! Incoming scalar current: J_S,    incoming momentum P2(4) (light-cone rep.)
! Outgoing scalar current: Jout_S = J_V.(2*P2+P1) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: mom1, mom2
  type(wfun),        intent(in)    :: J_V(:), J_S(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  complex(REALKIND) :: P1(4), P2(4)
  integer:: h

  P1 = get_LC_4(mom1)
  P2 = get_LC_4(mom2)

  do h=1, n(3)
    Jout_S(h)%j(1) = cont_VV(P1+P2+P2, J_V(t(1,h))%j(1:4)) * J_S(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_V(1)%n_part + J_S(1)%n_part
    Jout_S(:)%t = J_V(1)%t + J_S(1)%t
    Jout_S(:)%hf = J_V(t(1,:))%hf + J_S(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_V, J_S, Jout_S, n, t)
  end if

end subroutine counter_VS_T

! **********************************************************************
subroutine counter_TV_S(ntry, J_S, mom1, J_V, mom2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Incoming scalar current: J_S, incoming momentum P2(4) (light-cone rep.)
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing scalar current: Jout_S = J_V.(-2*P1-P2) * J_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_contractions_/**/REALKIND, only: cont_VV
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: mom1, mom2
  type(wfun),        intent(in)    :: J_S(:), J_V(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  complex(REALKIND) :: P1(4), P2(4)
  integer:: h

  P1 = get_LC_4(mom1)
  P2 = get_LC_4(mom2)
  do h=1, n(3)
    Jout_S(h)%j(1) = - cont_VV(P1+P1+P2, J_V(t(2,h))%j(1:4)) * J_S(t(1,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_S(1)%n_part + J_V(1)%n_part
    Jout_S(:)%t = J_S(1)%t + J_V(1)%t
    Jout_S(:)%hf = J_S(t(1,:))%hf + J_V(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S, J_V, Jout_S, n, t)
  end if

end subroutine counter_TV_S


! **********************************************************************
subroutine counter_ST_V(ntry, J_S1, mom1, J_S2, mom2, Jout_V, n, t)
! ----------------------------------------------------------------------
! Vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Vector boson + two scalars vertex
! Incoming scalar currents: J_S1, J_S2, incoming momenta P1(4), P2(4) (light-cone rep.)
! Outgoing vector current: Jout_V(4) = J_S1 * J_S2 * (P1 - P2)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: mom1, mom2
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  complex(REALKIND) :: P1(4), P2(4)
  integer:: h

  P1 = get_LC_4(mom1)
  P2 = get_LC_4(mom2)
  do h=1, n(3)
    Jout_V(h)%j = (J_S1(t(1,h))%j(1) * J_S2(t(2,h))%j(1)) * (P1 - P2)
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_S1(1)%n_part + J_S2(1)%n_part
    Jout_V(:)%t = J_S1(1)%t + J_S2(1)%t
    Jout_V(:)%hf = J_S1(t(1,:))%hf + J_S2(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S1, J_S2, Jout_V, n, t)
  end if

end subroutine counter_ST_V


! **********************************************************************
subroutine counter_VV_S(ntry, J_V1, J_V2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3, checkzero_scalar
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(wfun),        intent(out)   :: Jout_S(:)
  integer:: h

  do h=1, n(3)
    Jout_S(h)%j(1) = cont_VV(J_V1(t(1,h))%j(1:4),J_V2(t(2,h))%j(1:4))
  end do

  if (ntry == 1) then
    Jout_S(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part
    Jout_S(:)%t = J_V1(1)%t + J_V2(1)%t
    Jout_S(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf
    call checkzero_scalar(Jout_S)
    call helbookkeeping_vert3(ntry, J_V1, J_V2, Jout_S, n, t)
  end if

end subroutine counter_VV_S


! **********************************************************************
subroutine counter_VS_V(ntry, J_V, J_S, Jout_V, n, t)
! ----------------------------------------------------------------------
! Two vector boson + scalar vertex
! ----------------------------------------------------------------------
! Incoming vector current: J_V(4) (light-cone rep.)
! Incoming scalar current: J_S
! Outgoing vector current: Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:), J_S(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  integer:: h

  do h=1, n(3)
    Jout_V(h)%j(:) = J_V(t(1,h))%j(:) * J_S(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_V(1)%n_part + J_S(1)%n_part
    Jout_V(:)%t = J_V(1)%t + J_S(1)%t
    Jout_V(:)%hf = J_V(t(1,:))%hf + J_S(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_V, J_S, Jout_V, n, t)
  end if

end subroutine counter_VS_V


! **********************************************************************
subroutine counter_SV_V(ntry, J_S, J_V, Jout_V, n, t)
! ----------------------------------------------------------------------
! Two vector bosons + scalar vertex
! ----------------------------------------------------------------------
! Incoming scalar current: J_S
! Incoming vector current: J_V(4) (light-cone rep.)
! Outgoing vector current: Jout_V
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert3
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:), J_V(:)
  type(wfun),        intent(out)   :: Jout_V(:)
  integer:: h

  do h=1, n(3)
    Jout_V(h)%j(1:4) = J_S(t(1,h))%j(1) * J_V(t(2,h))%j(1:4)
  end do

  if (ntry == 1) then
    Jout_V(:)%n_part = J_S(1)%n_part + J_V(1)%n_part
    Jout_V(:)%t = J_S(1)%t + J_V(1)%t
    Jout_V(:)%hf = J_S(t(1,:))%hf + J_V(t(2,:))%hf
    call helbookkeeping_vert3(ntry, J_S, J_V, Jout_V, n, t)
  end if

end subroutine counter_SV_V

! **********************************************************************
subroutine counter_SSS_S(ntry, J_S1, J_S2, J_S3, Jout_S, n, t)
! ----------------------------------------------------------------------
! Four scalar vertex
! ----------------------------------------------------------------------
! Incoming scalar currents: J_S1, J_S2, J_S3
! Outgoing scalar current:  Jout_S = J_S1 * J_S2 * J_S3
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_S1(n(1)), J_S2(n(2)), J_S3(n(3))
  type(wfun), intent(out) :: Jout_S(:)
  integer :: h

  do h = 1, n(4)
    Jout_S(h)%j(1) = J_S1(t(1,h))%j(1) * J_S2(t(2,h))%j(1) * J_S3(t(2,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_S(:)%t = J_S1(1)%t + J_S2(1)%t + J_S3(1)%t
    Jout_S(:)%n_part = J_S1(1)%n_part + J_S2(1)%n_part + J_S3(1)%n_part
    Jout_S(:)%hf = J_S1(t(1,:))%hf + J_S2(t(2,:))%hf + J_S3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_S1, J_S2, J_S3, Jout_S, n, t)
  end if

end subroutine counter_SSS_S

! **********************************************************************
subroutine counter_VVV_V(ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
! ----------------------------------------------------------------------
! neutral Vector-vector-vector-vector pure R2
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [g(a1,a2)*g(a3,a4) + g(a2,a3)*g(a1,a4)
!               + g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_V1(n(1)), J_V2(n(2)), J_V3(n(3))
  type(wfun), intent(out) :: Jout_V(:)
  integer :: h
  complex(REALKIND) :: J1J2, J1J3, J2J3

  do h = 1, n(4)
    J1J2 = cont_VV(J_V1(t(1,h))%j(1:4), J_V2(t(2,h))%j(1:4))
    J1J3 = cont_VV(J_V1(t(1,h))%j(1:4), J_V3(t(3,h))%j(1:4))
    J2J3 = cont_VV(J_V2(t(2,h))%j(1:4), J_V3(t(3,h))%j(1:4))
    Jout_V(h)%j(:) = J2J3 * J_V1(t(1,h))%j(:) + J1J3 * J_V2(t(2,h))%j(:) +&
                     J1J2 * J_V3(t(3,h))%j(:)
  end do

  if (ntry == 1) then
    Jout_V(:)%t = J_V1(1)%t + J_V2(1)%t + J_V3(1)%t
    Jout_V(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part + J_V3(1)%n_part
    Jout_V(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf + J_V3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
  end if

end subroutine counter_VVV_V

! **********************************************************************
subroutine counter_WWV_V(ctWWVV, ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
! ----------------------------------------------------------------------
! W-W-vector-vector UV & pure R2
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [ctVVVV(1) * g(a1,a2)*g(a3,a4) + ctVVVV(2) * g(a2,a3)*g(a1,a4)
!               + ctVVVV(2) * g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_V1(n(1)), J_V2(n(2)), J_V3(n(3))
  complex(REALKIND), intent(in)  :: ctWWVV(2)
  type(wfun), intent(out) :: Jout_V(:)
  integer :: h
  complex(REALKIND) :: J1J2, J1J3, J2J3

  do h = 1, n(4)
    J1J2 = cont_VV(J_V1(t(1,h))%j(1:4), J_V2(t(2,h))%j(1:4))
    J1J3 = cont_VV(J_V1(t(1,h))%j(1:4), J_V3(t(3,h))%j(1:4))
    J2J3 = cont_VV(J_V2(t(2,h))%j(1:4), J_V3(t(3,h))%j(1:4))
    Jout_V(h)%j(:) = J2J3 * J_V1(t(1,h))%j(:) * ctWWVV(2) +&
                     J1J3 * J_V2(t(2,h))%j(:) * ctWWVV(2) +&
                     J1J2 * J_V3(t(3,h))%j(:) * ctWWVV(1)
  end do

  if (ntry == 1) then
    Jout_V(:)%t = J_V1(1)%t + J_V2(1)%t + J_V3(1)%t
    Jout_V(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part + J_V3(1)%n_part
    Jout_V(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf + J_V3(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
  end if

end subroutine counter_WWV_V

! **********************************************************************
subroutine counter_VWW_V(ctWWVV, ntry, J_V1, J_V2, J_V3, Jout_V, n, t)
! ----------------------------------------------------------------------
! W-W-vector-vector UV & pure R2
! ----------------------------------------------------------------------
! J_Vi(4)    = incoming vector boson currents (light-cone rep.)
! Jout_V(4)  = outgoing vector boson current  (light-cone rep.)
! Jout_V(a4) = [ctVVVV(1) * g(a1,a2)*g(a3,a4) + ctVVVV(2) * g(a2,a3)*g(a1,a4)
!               + ctVVVV(2) * g(a1,a3)*g(a2,a4)] * J_V1(a1) * J_V2(a2) * J_V3(a3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_V1(n(1)), J_V2(n(2)), J_V3(n(3))
  complex(REALKIND), intent(in)  :: ctWWVV(2)
  type(wfun), intent(out) :: Jout_V(:)
  integer :: h

  call counter_WWV_V(ctWWVV,ntry,J_V2,J_V3,J_V1,Jout_V,n,t)

end subroutine counter_VWW_V

! **********************************************************************
subroutine counter_VVS_S(ntry, J_V1, J_V2, J_S, Jout_S, n, t)
! ----------------------------------------------------------------------
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Incoming scalar current:  J_S
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_V1(n(1)), J_V2(n(2)), J_S(n(3))
  type(wfun), intent(out) :: Jout_S(:)
  integer :: h

  do h = 1, n(4)
    Jout_S(h)%j(1) = cont_VV(J_V1(t(1,h))%j(1:4),J_V2(t(2,h))%j(1:4)) * J_S(t(3,h))%j(1)
  end do

  if (ntry == 1) then
    Jout_S(:)%t = J_V1(1)%t + J_V2(1)%t + J_S(1)%t
    Jout_S(:)%n_part = J_V1(1)%n_part + J_V2(1)%n_part + J_S(1)%n_part
    Jout_S(:)%hf = J_V1(t(1,:))%hf + J_V2(t(2,:))%hf + J_S(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_V1, J_V2, J_S, Jout_S, n, t)
  end if

end subroutine counter_VVS_S

! **********************************************************************
subroutine counter_SSV_V(ntry, J_S1, J_S2, J_V, Jout_V, n, t)
! ----------------------------------------------------------------------
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Incoming scalar currents: J_S1, J_S2
! Incoming vector current:  J_V(4) (light-cone rep.)
! Outgoing vector current:  Jout_V (light-cone rep.)
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_S1(n(1)), J_S2(n(2)), J_V(n(3))
  type(wfun), intent(out) :: Jout_V(:)
  integer :: h

  do h = 1, n(4)
    Jout_V(h)%j(:) = (J_S1(t(1,h))%j(1) * J_S2(t(2,h))%j(1)) * J_V(t(3,h))%j(:)
  end do

  if (ntry == 1) then
    Jout_V(:)%t = J_S1(1)%t + J_S2(1)%t + J_V(1)%t
    Jout_V(:)%n_part = J_S1(1)%n_part + J_S2(1)%n_part + J_V(1)%n_part
    Jout_V(:)%hf = J_S1(t(1,:))%hf + J_S2(t(2,:))%hf + J_V(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_S1, J_S2, J_V, Jout_V, n, t)
  end if

end subroutine counter_SSV_V

! **********************************************************************
subroutine counter_SVV_S(ntry, J_S, J_V1, J_V2, Jout_S, n, t)
! ----------------------------------------------------------------------
! Two vector boson + two scalars vertex
! ----------------------------------------------------------------------
! Incoming scalar current:  J_S
! Incoming vector currents: J_V1(4), J_V2(4) (light-cone rep.)
! Outgoing scalar current:  Jout_S
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: helbookkeeping_vert4
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  integer(intkind1), intent(in) :: ntry
  integer(intkind2), intent(inout) :: n(4), t(3,n(4))
  type(wfun), intent(in)  :: J_S(n(1)), J_V1(n(2)), J_V2(n(3))
  type(wfun), intent(out) :: Jout_S(:)
  integer :: h

  do h = 1, n(4)
    Jout_S(h)%j(1) = J_S(t(1,h))%j(1) * cont_VV(J_V1(t(2,h))%j(1:4),J_V2(t(3,h))%j(1:4))
  end do

  if (ntry == 1) then
    Jout_S(:)%t = J_S(1)%t + J_V1(1)%t + J_V2(1)%t
    Jout_S(:)%n_part = J_S(1)%n_part + J_V1(1)%n_part + J_V2(1)%n_part
    Jout_S(:)%hf = J_S(t(1,:))%hf + J_V1(t(2,:))%hf + J_V2(t(3,:))%hf
    call helbookkeeping_vert4(ntry, J_S, J_V1, J_V2, Jout_S, n, t)
  end if

end subroutine counter_SVV_S

end module ol_h_counterterms_/**/REALKIND
