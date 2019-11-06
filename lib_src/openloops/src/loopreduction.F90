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


module ol_loop_reduction_/**/REALKIND
  use ol_debug, only: ol_fatal, ol_msg, ol_error
  use KIND_TYPES, only: REALKIND
  implicit none

  ! Thresholds used for the expansions
  real(REALKIND) :: m_thres = 1.E-9
  real(REALKIND) :: delta_thres = 1.E-2

  ! Thresholds for qp trigger in triangles
  real(REALKIND) :: coll_thres = 1d-4
  real(REALKIND) :: soft_thres = 1d-4

#ifdef PRECISION_dp
  ! In order to switch off the expansions, set DeltaExp = .false.
  logical :: DeltaExp = .true.
#else
  ! Expansions are switched on in quad-precision by default.
  ! Set DeltaExp = .false. in the following to siwtch them off in quad-precision
  logical :: DeltaExp = .true.
#endif

  interface TI_reduction
    module procedure TI_reduction_1
  end interface


contains


! ******************************************************************************
subroutine TI_bubble_red(Gin_A,momid,msq,Gout_A,M2R1,A0_0,A0_1)
! Reduction of tensor bubbles
! ------------------------------------------------------------------------------
! Gin_A   = input closed-loop two-point integral
! msq     = array containing the squared masses m0^2, m1^2
! Gout_A  = scalar output bubble
! M2R1    = R1 rational part of the amplitude
! A0_i    = tadpole with mass m_i
! ******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ofred_reduction_/**/REALKIND, only: twopoint_reduction
  use ofred_reduction_/**/QREALKIND, only: twopoint_reduction_qp=>twopoint_reduction
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hybrid_qp_mode,hp_alloc_mode
  use ol_data_types_/**/REALKIND, only: hcl, met
  use ol_kinematics_/**/REALKIND, only: get_mass2,get_LC_5
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: req_qp_cmp,hcl_alloc_hybrid
  use ol_kinematics_/**/QREALKIND, only: get_mass2_qp=>get_mass2
  use ol_kinematics_/**/DREALKIND, only: get_LC_5_qp
#endif
  type(hcl),           intent(in)    :: Gin_A
  integer,             intent(in)    :: momid
  integer,             intent(in)    :: msq(0:1)
  type(hcl),           intent(inout) :: Gout_A
  type(met),           intent(inout) :: M2R1
  type(hcl), optional, intent(inout) :: A0_0, A0_1
  complex(REALKIND)  :: p(1:5)
  complex(REALKIND) :: redcoeff(1:4)
#ifdef PRECISION_dp
  complex(QREALKIND) :: p_qp(1:5)
  complex(QREALKIND) :: redcoeff_qp(1:4)
#endif


  Gout_A%mode = Gin_A%mode
  M2R1%mode = Gin_A%mode
  if (present(A0_0) .AND. present(A0_1)) then
    A0_0%mode = Gin_A%mode
    A0_1%mode = Gin_A%mode
  else if(present(A0_0)) then
    A0_0%mode = Gin_A%mode
  end if

  Gout_A%error = Gin_A%error
  M2R1%error = Gin_A%error
  if (present(A0_0) .AND. present(A0_1)) then
    A0_0%error = Gin_A%error
    A0_1%error = Gin_A%error
    A0_0%ndrs = 0
    A0_0%nred = 0
    A0_1%ndrs = 0
    A0_1%nred = 0
#ifdef PRECISION_dp
    A0_0%ndrs_qp = 0
    A0_0%nred_qp = 0
    A0_1%ndrs_qp = 0
    A0_1%nred_qp = 0
#endif
  else if(present(A0_0)) then
    A0_0%error = Gin_A%error
    A0_0%ndrs = 0
    A0_0%nred = 0
#ifdef PRECISION_dp
    A0_0%ndrs_qp = 0
    A0_0%nred_qp = 0
#endif
  end if

#ifdef PRECISION_dp
  Gout_A%ndrs = Gin_A%ndrs
  Gout_A%nred = Gin_A%nred
  Gout_A%ndrs_qp = Gin_A%ndrs_qp
  Gout_A%nred_qp = Gin_A%nred_qp
  if (req_qp_cmp(Gin_A)) then
    Gout_A%nred_qp = Gout_A%nred_qp + 1
  else
    Gout_A%nred = Gout_A%nred + 1
  end if
#else
  Gout_A%ndrs = Gin_A%ndrs
  Gout_A%nred = Gin_A%nred + 1
#endif

  p = get_LC_5(momid)

#ifdef PRECISION_dp
  if (Gin_A%mode .ne. hybrid_qp_mode) then
#endif
    redcoeff = 0._/**/REALKIND

    call twopoint_reduction(Gin_A%cmp,p,get_mass2(msq),redcoeff)

    Gout_A%cmp(1) = redcoeff(1)
    M2R1%cmp = M2R1%cmp + redcoeff(4)

    if (present(A0_0) .AND. present(A0_1)) then
      A0_0%cmp(1) = get_mass2(msq(0))*redcoeff(2)
      A0_1%cmp(1) = get_mass2(msq(1))*redcoeff(3)
    else if(present(A0_0)) then
      if (msq(0) == msq(1)) then
        A0_0%cmp(1) = get_mass2(msq(0))*(redcoeff(2) + redcoeff(3))
      else if (msq(0) == 0) then
        A0_0%cmp(1) = get_mass2(msq(1))*redcoeff(3)
      else
        A0_0%cmp(1) = get_mass2(msq(0))*redcoeff(2)
      end if
    end if

#ifdef PRECISION_dp
  else
    Gout_A%cmp(1) = 0
    if (present(A0_0) .AND. present(A0_1)) then
      A0_0%cmp(1) = 0
      A0_1%cmp(1) = 0
    else if(present(A0_0)) then
      A0_0%cmp(1) = 0
    end if
  end if
#endif

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin_A)) then

    p_qp = get_LC_5_qp(momid)
    redcoeff_qp = 0._/**/QREALKIND

    call twopoint_reduction_qp(Gin_A%cmp_qp,p_qp,get_mass2_qp(msq),redcoeff_qp)
    if (hp_alloc_mode .gt. 1) call hcl_alloc_hybrid(Gout_A)

    Gout_A%cmp_qp(1) = redcoeff_qp(1)
    M2R1%cmp_qp = M2R1%cmp_qp + redcoeff_qp(4)

    if (present(A0_0) .AND. present(A0_1)) then
      if (hp_alloc_mode .gt. 1) then
        call hcl_alloc_hybrid(A0_0)
        call hcl_alloc_hybrid(A0_1)
      end if
      A0_0%cmp_qp(1) = get_mass2_qp(msq(0))*redcoeff_qp(2)
      A0_1%cmp_qp(1) = get_mass2_qp(msq(1))*redcoeff_qp(3)
    else if(present(A0_0)) then
      if (hp_alloc_mode .gt. 1) call hcl_alloc_hybrid(A0_0)
      if (msq(0) == msq(1)) then
        A0_0%cmp_qp(1) = get_mass2_qp(msq(0))*(redcoeff_qp(2) + redcoeff_qp(3))
      else if (msq(0) == 0) then
        A0_0%cmp_qp(1) = get_mass2_qp(msq(1))*redcoeff_qp(3)
      else
        A0_0%cmp_qp(1) = get_mass2_qp(msq(0))*redcoeff_qp(2)
      end if
    end if
  else if (hp_switch .eq. 1 .and. hp_alloc_mode .eq. 0) then
    Gout_A%cmp_qp(1) = 0
    if (present(A0_0) .AND. present(A0_1)) then
      A0_0%cmp_qp(1) = 0
      A0_1%cmp_qp(1) = 0
    else if(present(A0_0)) then
      A0_0%cmp_qp(1) = 0
    end if
  end if
#endif

end subroutine TI_bubble_red


! ******************************************************************************
subroutine TI_triangle_red(Gin_A,RedBasis,msq,Gout_A,Gout_A0,Gout_A1, &
                           Gout_A2,M2R1,A0msq,A0_0,A0_1,A0_2)
! Reduction of tensor triangles
! ------------------------------------------------------------------------------
! Gin_A    = input closed-loop three-point integral
! RedBasis = Reduction basis
! msq      = array of the squared masses m_0^2, m_1^2, m_2^2
! Gout_A   = output scalar triangle
! Gout_Ai  = output scalar bubbles originating from Di pinch
! M2R1     = R1 rational part of the amplitude
! A0msq    = squared masses in the tadpoles
! A0_i     = tadpole of mass m_i
! ******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_data_types_/**/REALKIND, only: basis
  use ofred_reduction_/**/REALKIND, only: otf_3pt_reduction_last
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,      &
                                              hybrid_qp_mode, &
                                              hp_irtri_trig,  &
                                              hp_metri_trig,  &
                                              hp_alloc_mode
  use ol_data_types_/**/REALKIND, only: hcl,met
  use ol_kinematics_/**/REALKIND, only: get_mass2
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: req_qp_cmp,upgrade_qp,hcl_alloc_hybrid
  use ol_loop_reduction_/**/QREALKIND, only: &
      tch_triangle_check_qp=>tch_triangle_check, &
      t_channel_triangle_reduction_qp=>t_channel_triangle_reduction
  use ofred_reduction_/**/QREALKIND, only: &
      otf_3pt_reduction_last_qp=>otf_3pt_reduction_last
  use ol_data_types_/**/QREALKIND, only: basis_qp=>basis
  use ofred_basis_construction_/**/REALKIND, only: normalize_gamma
  use ofred_basis_construction_/**/QREALKIND, only: &
    construct_RedBasis_qp => construct_RedBasis
  use ol_kinematics_/**/QREALKIND, only: get_mass2_qp=>get_mass2
#endif
  implicit none
  type(hcl),           intent(inout) :: Gin_A
  type(basis),         intent(in)    :: RedBasis
  integer,             intent(in)    :: msq(0:2)
  type(hcl),           intent(inout) :: Gout_A,Gout_A0,Gout_A1,Gout_A2
  type(met),           intent(inout) :: M2R1
  integer, optional,   intent(in)    :: A0msq(:)
  type(hcl), optional, intent(inout) :: A0_0, A0_1, A0_2
  logical :: unstable
#ifdef PRECISION_dp
  type(basis_qp)  :: Redbasis_qp
#endif

  complex(REALKIND) :: Gout_R1, momenta(5,3),p(0:3)
  complex(QREALKIND) :: Gout_R1_qp, momenta_qp(5,3)
  real(REALKIND) :: sdlt,gn
  real(QREALKIND) :: sdlt_qp
  integer :: perm(3), m_ind
  logical :: tch_top

  !! Check for a t-channel triangle topology with a massless external leg
  call tch_triangle_check(RedBasis%mom1, RedBasis%mom2, tch_top, perm, sdlt, momenta)

#ifdef PRECISION_dp
  unstable = .false.
  if (hp_switch .eq. 1 .and. .not. req_qp_cmp(Gin_A)) then
    if (hp_metri_trig) call tri_check_me(unstable)
    if (hp_irtri_trig .and. .not. unstable) call triangle_check_ir(RedBasis%mom1, RedBasis%mom2, unstable)
  end if
#endif

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. unstable) then
    call upgrade_qp(Gin_A)
  end if
#endif

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin_A)) then
    call tch_triangle_check_qp(RedBasis%mom1, RedBasis%mom2, tch_top, perm, sdlt_qp, momenta_qp)
  end if
#endif

  Gout_A%mode = Gin_A%mode
  Gout_A0%mode = Gin_A%mode
  Gout_A1%mode = Gin_A%mode
  Gout_A2%mode = Gin_A%mode
  M2R1%mode = Gin_A%mode

  Gout_A%error = Gin_A%error
  Gout_A0%error = Gin_A%error
  Gout_A1%error = Gin_A%error
  Gout_A2%error = Gin_A%error
  M2R1%error = Gin_A%error
  if (present(A0_2)) then
    A0_2%mode = Gin_A%mode
    A0_2%error = Gin_A%error
    A0_2%ndrs = 0
    A0_2%nred = 0
#ifdef PRECISION_dp
    A0_2%ndrs_qp = 0
    A0_2%nred_qp = 0
#endif
  end if
  if (present(A0_1)) then
    A0_1%mode = Gin_A%mode
    A0_1%error = Gin_A%error
    A0_1%ndrs = 0
    A0_1%nred = 0
#ifdef PRECISION_dp
    A0_1%ndrs_qp = 0
    A0_1%nred_qp = 0
#endif
  end if
  if (present(A0_0)) then
    A0_0%mode = Gin_A%mode
    A0_0%error = Gin_A%error
    A0_0%ndrs = 0
    A0_0%nred = 0
#ifdef PRECISION_dp
    A0_0%ndrs_qp = 0
    A0_0%nred_qp = 0
#endif
  end if


#ifdef PRECISION_dp
  Gout_A%ndrs = Gin_A%ndrs
  Gout_A0%ndrs = 0
  Gout_A1%ndrs = 0
  Gout_A2%ndrs = 0
  Gout_A%ndrs_qp = Gin_A%ndrs_qp
  Gout_A0%ndrs_qp = 0
  Gout_A1%ndrs_qp = 0
  Gout_A2%ndrs_qp = 0
  Gout_A%nred = Gin_A%nred
  Gout_A0%nred = 0
  Gout_A1%nred = 0
  Gout_A2%nred = 0
  Gout_A%nred_qp = Gin_A%nred_qp
  Gout_A0%nred_qp = 0
  Gout_A1%nred_qp = 0
  Gout_A2%nred_qp = 0
  if (req_qp_cmp(Gin_A)) then
    Gout_A%nred_qp = Gout_A%nred_qp + 1
  else
    Gout_A%nred = Gout_A%nred + 1
  end if
#else
  Gout_A%ndrs = Gin_A%ndrs
  Gout_A0%ndrs = 0
  Gout_A1%ndrs = 0
  Gout_A2%ndrs = 0
  Gout_A%nred = Gin_A%nred + 1
  Gout_A0%nred = 0
  Gout_A1%nred = 0
  Gout_A2%nred = 0
#endif

#ifdef PRECISION_dp
  if (Gin_A%mode .ne. hybrid_qp_mode) then
#endif
    Gout_R1 = 0._/**/REALKIND

    if(tch_top) then
    !! Exact reduction formulas or Gram-\Delta expansions
      if (present(A0_2)) then
        call t_channel_triangle_reduction(perm,sdlt,momenta,get_mass2(msq),&
                                          Gin_A%cmp(:),Gout_A%cmp(:),      &
                                          Gout_A0%cmp(:),Gout_A1%cmp(:),   &
                                          Gout_A2%cmp(:),Gout_R1,          &
                                          get_mass2(A0msq),A0_0%cmp(:),    &
                                          A0_1%cmp(:),A0_2%cmp(:))
      else if(present(A0_1)) then
        call t_channel_triangle_reduction(perm,sdlt,momenta,get_mass2(msq),&
                                          Gin_A%cmp(:),Gout_A%cmp(:),      &
                                          Gout_A0%cmp(:),Gout_A1%cmp(:),   &
                                          Gout_A2%cmp(:),Gout_R1,          &
                                          get_mass2(A0msq),A0_0%cmp(:),    &
                                          A0_1%cmp(:))
      else if(present(A0_0)) then
        call t_channel_triangle_reduction(perm,sdlt,momenta,get_mass2(msq),&
                                          Gin_A%cmp(:),Gout_A%cmp(:),      &
                                          Gout_A0%cmp(:),Gout_A1%cmp(:),   &
                                          Gout_A2%cmp(:),Gout_R1,          &
                                          get_mass2(A0msq),A0_0%cmp(:))
      else
        call t_channel_triangle_reduction(perm,sdlt,momenta,get_mass2(msq),&
                                          Gin_A%cmp(:),Gout_A%cmp(:),      &
                                          Gout_A0%cmp(:),Gout_A1%cmp(:),   &
                                          Gout_A2%cmp(:),Gout_R1)
      end if

    else

      !! On-the-fly reduction-like
      if(present(A0_2)) then
        call otf_3pt_reduction_last(Gin_A%cmp(:),RedBasis,get_mass2(msq),  &
                                    Gout_A%cmp(:),Gout_A0%cmp(:),          &
                                    Gout_A1%cmp(:),Gout_A2%cmp(:),Gout_R1, &
                                    get_mass2(A0msq),A0_0%cmp(:),          &
                                    A0_1%cmp(:),A0_2%cmp(:))
      else if(present(A0_1)) then
        call otf_3pt_reduction_last(Gin_A%cmp,RedBasis,get_mass2(msq),     &
                                    Gout_A%cmp(:),Gout_A0%cmp(:),          &
                                    Gout_A1%cmp(:),Gout_A2%cmp(:),Gout_R1, &
                                    get_mass2(A0msq),A0_0%cmp(:),A0_1%cmp(:))
      else if(present(A0_0)) then
        call otf_3pt_reduction_last(Gin_A%cmp(:),RedBasis,get_mass2(msq),  &
                                    Gout_A%cmp(:),Gout_A0%cmp(:),          &
                                    Gout_A1%cmp(:),Gout_A2%cmp(:),Gout_R1, &
                                    get_mass2(A0msq),A0_0%cmp(:))
      else
        call otf_3pt_reduction_last(Gin_A%cmp(:),RedBasis,get_mass2(msq), &
                                    Gout_A%cmp(:),Gout_A0%cmp(:),         &
                                    Gout_A1%cmp(:),Gout_A2%cmp(:),Gout_R1)
      end if

    end if

    M2R1%cmp = M2R1%cmp + Gout_R1
#ifdef PRECISION_dp
  else
    Gout_A%cmp = 0
    Gout_A0%cmp = 0
    Gout_A1%cmp = 0
    Gout_A2%cmp = 0
    if (present(A0_2)) A0_2%cmp = 0
    if (present(A0_1)) A0_1%cmp = 0
    if (present(A0_0)) A0_0%cmp = 0
  end if
#endif

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin_A)) then

  if (hp_alloc_mode .gt. 1) then
    call hcl_alloc_hybrid(Gout_A)
    call hcl_alloc_hybrid(Gout_A0)
    call hcl_alloc_hybrid(Gout_A1)
    call hcl_alloc_hybrid(Gout_A2)
  end if

  Gout_R1_qp = 0._/**/QREALKIND

  if(tch_top) then
  !! Exact reduction formulas or Gram-\Delta expansions
    if (present(A0_2)) then
      call hcl_alloc_hybrid(A0_0)
      call hcl_alloc_hybrid(A0_1)
      call hcl_alloc_hybrid(A0_2)
      call t_channel_triangle_reduction_qp(perm,sdlt_qp,momenta_qp,            &
                                           get_mass2_qp(msq),Gin_A%cmp_qp(:),  &
                                           Gout_A%cmp_qp(:),Gout_A0%cmp_qp(:), &
                                           Gout_A1%cmp_qp(:),Gout_A2%cmp_qp(:),&
                                           Gout_R1_qp,get_mass2_qp(A0msq),     &
                                           A0_0%cmp_qp(:),A0_1%cmp_qp(:),      &
                                           A0_2%cmp_qp(:))
    else if(present(A0_1)) then
      call hcl_alloc_hybrid(A0_0)
      call hcl_alloc_hybrid(A0_1)
      call t_channel_triangle_reduction_qp(perm,sdlt_qp,momenta_qp,            &
                                           get_mass2_qp(msq),Gin_A%cmp_qp(:),  &
                                           Gout_A%cmp_qp(:),Gout_A0%cmp_qp(:), &
                                           Gout_A1%cmp_qp(:),Gout_A2%cmp_qp(:),&
                                           Gout_R1_qp,get_mass2_qp(A0msq),     &
                                           A0_0%cmp_qp(:),A0_1%cmp_qp(:))
    else if(present(A0_0)) then
      call hcl_alloc_hybrid(A0_0)
      call t_channel_triangle_reduction_qp(perm,sdlt_qp,momenta_qp,            &
                                           get_mass2_qp(msq),Gin_A%cmp_qp(:),  &
                                           Gout_A%cmp_qp(:),Gout_A0%cmp_qp(:), &
                                           Gout_A1%cmp_qp(:),Gout_A2%cmp_qp(:),&
                                           Gout_R1_qp,get_mass2_qp(A0msq),     &
                                           A0_0%cmp_qp(:))
    else
      call t_channel_triangle_reduction_qp(perm,sdlt_qp,momenta_qp,            &
                                           get_mass2_qp(msq),Gin_A%cmp_qp(:),  &
                                           Gout_A%cmp_qp(:),Gout_A0%cmp_qp(:), &
                                           Gout_A1%cmp_qp(:),Gout_A2%cmp_qp(:),&
                                           Gout_R1_qp)
    end if

  else
    call construct_RedBasis_qp(RedBasis%mom1, &
                               RedBasis%mom2, &
                               Redbasis_qp)

    !! On-the-fly reduction-like
    if(present(A0_2)) then
      call hcl_alloc_hybrid(A0_0)
      call hcl_alloc_hybrid(A0_1)
      call hcl_alloc_hybrid(A0_2)
      call otf_3pt_reduction_last_qp(Gin_A%cmp_qp(:),Redbasis_qp,         &
                                     get_mass2_qp(msq),Gout_A%cmp_qp(:),  &
                                     Gout_A0%cmp_qp(:),Gout_A1%cmp_qp(:), &
                                     Gout_A2%cmp_qp(:),Gout_R1_qp,        &
                                     get_mass2_qp(A0msq),A0_0%cmp_qp(:),  &
                                     A0_1%cmp_qp(:),A0_2%cmp_qp(:))
    else if(present(A0_1)) then
      call hcl_alloc_hybrid(A0_0)
      call hcl_alloc_hybrid(A0_1)
      call otf_3pt_reduction_last_qp(Gin_A%cmp_qp,Redbasis_qp,            &
                                     get_mass2_qp(msq),Gout_A%cmp_qp(:),  &
                                     Gout_A0%cmp_qp(:),Gout_A1%cmp_qp(:), &
                                     Gout_A2%cmp_qp(:),Gout_R1_qp,        &
                                     get_mass2_qp(A0msq),A0_0%cmp_qp(:),  &
                                     A0_1%cmp_qp(:))
    else if(present(A0_0)) then
      call hcl_alloc_hybrid(A0_0)
      call otf_3pt_reduction_last_qp(Gin_A%cmp_qp(:),Redbasis_qp,         &
                                     get_mass2_qp(msq),Gout_A%cmp_qp(:),  &
                                     Gout_A0%cmp_qp(:),Gout_A1%cmp_qp(:), &
                                     Gout_A2%cmp_qp(:),Gout_R1_qp,        &
                                     get_mass2_qp(A0msq),A0_0%cmp_qp(:))
    else
      call otf_3pt_reduction_last_qp(Gin_A%cmp_qp(:),Redbasis_qp,         &
                                     get_mass2_qp(msq),Gout_A%cmp_qp(:),  &
                                     Gout_A0%cmp_qp(:),Gout_A1%cmp_qp(:), &
                                     Gout_A2%cmp_qp(:),Gout_R1_qp)
    end if

  end if


  M2R1%cmp_qp = M2R1%cmp_qp + Gout_R1_qp

  else if (hp_switch .eq. 1 .and. hp_alloc_mode .eq. 0) then
    Gout_A%cmp_qp = 0
    Gout_A0%cmp_qp = 0
    Gout_A1%cmp_qp = 0
    Gout_A2%cmp_qp = 0
    if (present(A0_2)) A0_2%cmp_qp = 0
    if (present(A0_1)) A0_1%cmp_qp = 0
    if (present(A0_0)) A0_0%cmp_qp = 0
  end if
#endif

  contains

  ! trigger for missing expansions
  subroutine tri_check_me(unstable)
    logical, intent(out) :: unstable
    integer              :: masses(0:2)

    unstable = .false.
    if (DeltaExp .and. tch_top) then
      if (sdlt < delta_thres) unstable = .true.
    end if

    if (unstable) then
      if(perm(1) == 2) then

        if(perm(2) == 0) then
          masses = (/msq(2),msq(0),msq(1)/)
        else
          masses = (/msq(2),msq(1),msq(0)/)
        end if

      else if(perm(1) == 0) then

        if(perm(2) == 2) then
          masses = (/msq(0),msq(2),msq(1)/)

        else
          masses = (/msq(0),msq(1),msq(2)/)

        end if

      else if(perm(1) == 1) then

        if(perm(2) == 2) then
          masses = (/msq(1),msq(2),msq(0)/)
        else
          masses = (/msq(1),msq(0),msq(2)/)
        end if
      end if

      ! the following cases are handled by trred
      if (get_mass2(masses(1))==get_mass2(masses(2)) .AND. &
          get_mass2(masses(1))==get_mass2(masses(0)) .AND. &
          get_mass2(masses(0))==0) then
        ! (0,0,0) masses configuration
        unstable = .false.
      else if (get_mass2(masses(1))==get_mass2(masses(2)) .AND. &
        get_mass2(masses(1))==0) then
        ! (m0,0,0) masses configuration
        unstable = .false.
      else if (get_mass2(masses(1))==get_mass2(masses(2)) .AND. &
        get_mass2(masses(1))==get_mass2(masses(0))) then
        ! (m,m,m) masses configuration
        unstable = .false.
      else if (get_mass2(masses(1))==get_mass2(masses(2)) .AND. &
        get_mass2(masses(0))==0) then
        ! (0,m,m) masses configuration
        unstable = .false.
      else if (get_mass2(masses(1))==get_mass2(masses(2))) then
        ! (m0,m1,m1) masses configuration
        unstable = .false.
      end if
    end if

  end subroutine tri_check_me

end subroutine TI_triangle_red


! ***************************************************************************************
!                            t-channel triangle topology
! ***************************************************************************************
!               -->                      <--
!               p1 _____________________ -p2
!                     \             /
!                      \           /
!   p1^2= -p^2          \         /
!                        \       /
!   p2^2= -p^2(1+dlt)     \     /
!                          \   /
!                           \ /
!                            |
!                            |
!                         (p2-p1)^2 = 0
!
! ***************************************************************************************


! ******************************************************************************
subroutine tch_triangle_check(mom1, mom2, tr_t_top, perm, sdlt, mom)
! ------------------------------------------------------------------------------
! Check for a t-channel triangle topology with one external massless leg
! TODO: this check might be cached
! ******************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_loop_parameters_decl_/**/REALKIND, only: mureg, zero
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx
  implicit none
  integer, intent(in) :: mom1, mom2
  logical, intent(out) :: tr_t_top
  integer, intent(out) :: perm(3)
  complex(REALKIND), intent(out) :: mom(5,3)
  real(REALKIND), intent(out) :: sdlt
  logical :: zero_mass(3)
  complex(REALKIND) :: k1(5), k2(5), k12(5), p(0:3)
  real(REALKIND) :: rone

  rone = 1._/**/REALKIND

  !! k1 and k2 are the external momenta flowing into the triangle
  k1(1:4)  = L(1:4,mom1)
  k1(5)    = L(5,mom1) + L(6,mom1)
  k2(1:4)  = L(1:4,mom2-mom1)
  k2(5)    = L(5,mom2-mom1) + L(6,mom2-mom1)
  k12(1:4) = L(1:4,mom2)
  k12(5)   = L(5,mom2) + L(6,mom2)

  tr_t_top  = .FALSE.
  zero_mass = .FALSE.

  !! Check for an external massless leg.
  zero_mass(1) = k1(5) .eq. zero
  zero_mass(2) = k2(5) .eq. zero
  zero_mass(3) = k12(5) .eq. zero

  if((zero_mass(1) .and. zero_mass(2)) .or. &
     (zero_mass(1) .and. zero_mass(3)) .or. &
     (zero_mass(2) .and. zero_mass(3)) ) return

  if(zero_mass(1) .AND. REAL(k12(5)*k2(5)) > 0) then
    tr_t_top = .TRUE.
    mom(1:5,1) =   k1(1:5)
    mom(1:4,2) = - k2(1:4)
    mom(5,2)   =   k2(5)
    mom(1:4,3) = - k12(1:4)
    mom(5,3)   =   k12(5)
    if(abs(k2(5)) > abs(k12(5))) then
      perm = [2,0,1]
      sdlt = abs(k2(5))/abs(k12(5)) - rone
    else
      perm = [2,1,0]
      sdlt = abs(k12(5))/abs(k2(5)) - rone
    end if

  else if(zero_mass(2) .AND. REAL(k1(5)*k12(5)) > 0) then
    tr_t_top = .TRUE.
    mom(1:5,1) = k1(1:5)
    mom(1:5,2) = k2(1:5)
    mom(1:5,3) = k12(1:5)
    if(abs(k1(5)) > abs(k12(5))) then
      perm = [0,2,1]
      sdlt = abs(k1(5))/abs(k12(5)) - rone
    else
      perm = [0,1,2]
      sdlt = abs(k12(5))/abs(k1(5)) - rone
    end if

  else if(zero_mass(3) .AND. REAL(k1(5)*k2(5)) > 0) then
    tr_t_top = .TRUE.
    mom(1:4,1) = - k1(1:4)
    mom(5,1)   =   k1(5)
    mom(1:5,2) =   k2(1:5)
    mom(1:5,3) =   k12(1:5)
    if(abs(k1(5)) > abs(k2(5))) then
      perm = [1,2,0]
      sdlt = abs(k1(5))/abs(k2(5)) - rone
    else
      perm = [1,0,2]
      sdlt = abs(k2(5))/abs(k1(5)) - rone
    end if
  else
    tr_t_top = .FALSE.
  end if

end subroutine tch_triangle_check

! ******************************************************************************
subroutine triangle_check_ir(mom1, mom2, unstable)
! ------------------------------------------------------------------------------
! Check for a soft and collinear triangle topology with two external massless legs
! ******************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_loop_parameters_decl_/**/REALKIND, only: mureg,zero,rone
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx
  implicit none
  integer, intent(in) :: mom1, mom2
  logical, intent(out) :: unstable
  logical :: zero_mass(3), soft_mass(3)
  complex(REALKIND) :: k1(5), k2(5), k12(5), p(0:3)
  real(REALKIND) :: coll, soft1, soft2, sdlt
  complex(REALKIND) :: k1_std(0:3), k2_std(0:3), k12_std(0:3)

  !! k1 and k2 are the external momenta flowing into the triangle
  k1(1:4)  = L(1:4,mom1)
  k1(5)    = L(5,mom1) + L(6,mom1)
  k2(1:4)  = L(1:4,mom2-mom1)
  k2(5)    = L(5,mom2-mom1) + L(6,mom2-mom1)
  k12(1:4) = L(1:4,mom2)
  k12(5)   = L(5,mom2) + L(6,mom2)

  zero_mass = .FALSE.
  zero_mass(1) = k1(5) .eq. zero
  zero_mass(2) = k2(5) .eq. zero
  zero_mass(3) = k12(5).eq. zero

  soft_mass(1) = abs(k1(5) )/mureg**2 < soft_thres .and. k1(5) .ne. zero
  soft_mass(2) = abs(k2(5) )/mureg**2 < soft_thres .and. k2(5) .ne. zero
  soft_mass(3) = abs(k12(5))/mureg**2 < soft_thres .and. k12(5) .ne. zero

  unstable = .false.

  call LC2Std_Rep_cmplx(k1(1:4), k1_std)
  call LC2Std_Rep_cmplx(k2(1:4), k2_std)
  call LC2Std_Rep_cmplx(k12(1:4), k12_std)

  ! check for small invariants
  if (zero_mass(1) .and. zero_mass(2)) then
    coll = abs(k12(5)/maxval(abs(k12(1:4)))**2)
    soft1 = abs(k1_std(0))/ max( maxval(abs(k12_std(0:3))), maxval(abs(k2_std(0:3))) )
    soft2 = abs(k2_std(0))/ max( maxval(abs(k12_std(0:3))), maxval(abs(k1_std(0:3))) )
    if (coll .lt. coll_thres) then
      unstable = .true.
    else if ( min(soft1,soft2) .lt. soft_thres) then
      unstable = .true.
    end if

  else if (zero_mass(1) .and. zero_mass(3)) then
    coll = abs(k2(5)/maxval(abs(k2(1:4)))**2)
    soft1 = abs(k1_std(0))/ max( maxval(abs(k12_std(0:3))), maxval(abs(k2_std(0:3))) )
    soft2 = abs(k12_std(0))/ max( maxval(abs(k1_std(0:3))), maxval(abs(k2_std(0:3))) )
    if (coll .lt. coll_thres) then
      unstable = .true.
    else if ( min(soft1,soft2) .lt. soft_thres) then
      unstable = .true.
    end if
  else if (zero_mass(2) .and. zero_mass(3)) then
    coll = abs(k1(5)/maxval(abs(k1(1:4)))**2)
    soft1 = abs(k2_std(0))/ max( maxval(abs(k1_std(0:3))), maxval(abs(k12_std(0:3))) )
    soft2 = abs(k12_std(0))/ max( maxval(abs(k1_std(0:3))), maxval(abs(k2_std(0:3))) )
    if (coll .lt. coll_thres) then
      unstable = .true.
    else if ( min(soft1,soft2) .lt. soft_thres) then
      unstable = .true.
    end if
  end if

  ! check for soft tch configuration
  sdlt = rone
  if(soft_mass(1) .AND. REAL(k12(5)*k2(5)) > 0) then
    if(abs(k2(5)) > abs(k12(5))) then
      sdlt = ABS(k2(5)/k12(5)) - rone
    else
      sdlt = ABS(k12(5)/k2(5)) - rone
    end if
    if (sdlt .ne. rone .and. sdlt .lt. delta_thres .and. (REAL(k2(5)) < 0)) unstable =.true.

  else if(soft_mass(2) .AND. REAL(k1(5)*k12(5)) > 0) then
    if(abs(k1(5)) > abs(k12(5))) then
      sdlt = ABS(k1(5)/k12(5)) - rone
    else
      sdlt = ABS(k12(5)/k1(5)) - rone
    end if
    if (sdlt .ne. rone .and. sdlt .lt. delta_thres .and. (REAL(k1(5)) < 0)) unstable =.true.

  else if(soft_mass(3) .AND. REAL(k1(5)*k2(5)) > 0) then
    if(abs(k1(5)) > abs(k2(5))) then
      sdlt = ABS(k1(5)/k2(5)) - rone
    else
      sdlt = ABS(k2(5)/k1(5)) - rone
    end if
    if (sdlt .ne. rone .and. sdlt .lt. delta_thres .and. (REAL(k1(5)) < 0)) unstable =.true.
  end if

end subroutine triangle_check_ir


! ******************************************************************************
subroutine tch_triangle_exact(r,G_in,msq,masses_ind,dlt,psq,p1,p2,&
& redcoeff)
! ******************************************************************************
! OpenLoops Reduction. Reduction of three point integrals with one
! external massless leg. Exact reduction formulae
! ------------------------------------------------------------------------------
! G_in     = input OL coefficient
! r        = tensor integral rank
! msq      = array of internal masses squared
! dlt      = \delta parameter
! p1, p2   = external momenta in the triangle
! redcoeff = (1) scalar triangle
!            (2) D0-pinch
!            (3) D1-pinch
!            (4) D2-pinch
!            (5) m0-tadpole
!            (6) m1-tadpole
!            (7) m2-tadpole
!            (8) rational part
! masses_ind = integer for the internal masses configuration
! ******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_in(:)
  complex(REALKIND), intent(in) :: msq(0:2), p1(1:5), p2(1:5)
  integer, intent(in) :: r, masses_ind
  real(REALKIND), intent(in) :: dlt, psq
  complex(REALKIND), intent(out) :: redcoeff(1:8)

  complex(REALKIND) :: p12(1:4), psq_delta, psq_delta2
  complex(REALKIND) :: contr_p12, contr_p1

  complex(REALKIND) :: p1p1(6:15), p1p2(6:15), p1p12(6:15)
  complex(REALKIND) :: p1p1p1(16:35), p1_p2(16:35), p1p12p12(16:35), p1p1p12(16:35)
  complex(REALKIND) :: P11, P22, T, U, Gmn, W0, P111, V1, Gp1, Gp2, Gp12, Gm, Gp
  complex(REALKIND) :: T3, T6, U3, V16, P1116, P111z3, P16, V
  complex(REALKIND) :: dlt2, dlt3, dlt4, dlt5, dlt6, zero, one, two, five
  complex(REALKIND) :: z, z0, z01, z1, ct1, ct2, z12, z2, zd, z10, z02, z11, z22
  integer :: i, j, k, n

  zero = 0._/**/REALKIND
  one  = 1._/**/REALKIND
  two  = 2._/**/REALKIND
  five = 5._/**/REALKIND

  redcoeff = zero
  psq_delta = psq*dlt
  psq_delta2 = psq_delta*dlt

  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-1 part ---------------------*
  !------------------------------------------------------------------------*
  p12(1:4) = p1(1:4) - p2(1:4)

  if(r == 1) then

    contr_p1  = SUM(G_in(2:5)*p1(1:4) )  !! N_\mu p_1^{\mu}
    contr_p12 = SUM(G_in(2:5)*p12(1:4))  !! N_\mu (p_1^{\mu} - p_2^{\mu})

    !! masses: (0,0,0)
    if (masses_ind == 0) then
      redcoeff(1) = -contr_p12/dlt - contr_p1
      redcoeff(2) = - contr_p12/psq_delta !scaleless bubble
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,0,0)
    else if (masses_ind == 1) then
      z0  = msq(0)/psq
      z01 = one + z0

      redcoeff(1) = -z01*contr_p12/dlt - contr_p1
      redcoeff(2) = - contr_p12/psq_delta !scaleless bubble
      redcoeff(3) = (2*contr_p12)/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -(2*contr_p12)/psq_delta2 - contr_p1/psq_delta

    !! masses: (0,m1,0)
    else if (masses_ind == 2) then
      z1  = msq(1)/psq
      dlt2 = dlt**2
      redcoeff(1) = - contr_p1 + 2*z1*contr_p12/dlt2+(z1*contr_p1-(one-z1)*contr_p12)/dlt
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -(contr_p1/psq_delta + 2*contr_p12/psq_delta2)
      redcoeff(5) = -contr_p12/psq_delta

    !! masses: (0,0,m2)
    else if (masses_ind == 3) then
      z2  = msq(2)/psq
      dlt2 = dlt**2
      redcoeff(1) = -contr_p1 - 2*z2*contr_p12/dlt2 - (z2*contr_p1 + contr_p12)/dlt
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -(contr_p1/psq_delta + 2*contr_p12/psq_delta2)
      redcoeff(5) = -contr_p12/psq_delta

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      redcoeff(1) = -contr_p12/dlt - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (0,m,m)
    else if (masses_ind == 5) then
      redcoeff(1) = contr_p12*(-1/dlt + msq(1)/psq_delta) - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -(2*contr_p12)/psq_delta2 - contr_p1/psq_delta

    !! masses: (m,0,m)
    else if (masses_ind == 6) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z01 = one + z0
      redcoeff(1) = -2*z0*contr_p12/dlt2 - (contr_p1*z0 + z01*contr_p12)/dlt - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m,m,0)
    else if (masses_ind == 7) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z01 = one + z0
      redcoeff(1) = 2*z0*contr_p12/dlt2 - (contr_p12- contr_p1*z0)/dlt - contr_p1
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta
      redcoeff(5) = -contr_p12/psq_delta

    !! masses: (m0,m1,m1)
    else if (masses_ind == 8) then
      z0 = msq(0)/psq
      z1 = msq(1)/psq
      redcoeff(1) = contr_p12*(z1-z0-one)/dlt - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,m0,m2)
    else if (masses_ind == 9) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z2  = msq(2)/psq

      redcoeff(1) = 2*(z0-z2)*contr_p12/dlt2 +(contr_p1*(z0-z2)-contr_p12)/dlt - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,m1,m0)
    else if (masses_ind == 10) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z01 = one + z0 - z1
      redcoeff(1) = 2*(z1-z0)*contr_p12/dlt2 + (contr_p1*(z1-z0) - z01*contr_p12)/dlt - &
      contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (0,m1,m2)
    else if (masses_ind == 11) then
      dlt2 = dlt**2
      z1  = msq(1)/psq
      z2  = msq(2)/psq
      z12 = z1 - z2

      redcoeff(1) = 2*z12*contr_p12/dlt2 + (contr_p1*z12+(z1-one)*contr_p12)/dlt-contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,0,m2)
    else if (masses_ind == 12) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z2  = msq(2)/psq
      z02 = one + z0
      redcoeff(1) = -2*z2*contr_p12/dlt2 - (z2*contr_p1 + z02*contr_p12)/dlt-contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,m1,0)
    else if (masses_ind == 13) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z01 = one + z0 - z1
      redcoeff(1) = 2*z1*contr_p12/dlt2 + (z1*contr_p1 - z01*contr_p12)/dlt - contr_p1
      redcoeff(2) = -contr_p12/psq_delta
      redcoeff(3) =  2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta

    !! masses: (m0,m1,m2)
    else if (masses_ind == 14) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z2  = msq(2)/psq
      z01 = one + z0 - z1
      z12 = z1-z2

      redcoeff(1) = 2*z12*contr_p12/dlt2 + (z12*contr_p1 - z01*contr_p12)/dlt - contr_p1
      redcoeff(3) = 2*contr_p12/psq_delta2 + (contr_p1 + contr_p12)/psq_delta
      redcoeff(4) = -2*contr_p12/psq_delta2 - contr_p1/psq_delta
      redcoeff(6) = (z1/z12)*(-contr_p12/psq_delta)
      redcoeff(7) = (z2/z12)*(contr_p12/psq_delta)

    end if

  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-2 part ---------------------*
  !------------------------------------------------------------------------*
  else if (r == 2) then

    dlt2 = dlt**2
    dlt3 = dlt**3
    dlt4 = dlt**4

    n = 6
    do i = 1, 4
      do j = i, 4
        p1p1(n)  = p1(i)*p1(j)
        p1p2(n)  = p12(i)*p12(j)
        p1p12(n) = p1(i)*p12(j) + p12(i)*p1(j)
        n = n + 1
      end do
    end do

    Gmn = (G_in(7)-G_in(14))/2      !! g^{\mu\nu}/4
    P11 = SUM(p1p1*G_in(6:15))/psq  !! N_{\mu\nu} p_1^{\mu} p_1^{\nu}
    T   = SUM(p1p2*G_in(6:15))/psq  !! N_{\mu\nu} (p_1^{\mu} - p_2^{\mu})(p_1^{\nu} - p_2^{\nu}) + (\mu <-> \nu)
    V   = SUM(p1p12*G_in(6:15))/psq !! N_{\mu\nu} (p_1^{\mu} - p_2^{\mu})(p_1^{\nu}) + (\mu <-> \nu)

    !! masses: (0,0,0)
    if (masses_ind == 0) then
      redcoeff(1) = (T/dlt2 + V/dlt + P11)*psq
      redcoeff(2) = (T/dlt2 + (V - T/2)/dlt)
      redcoeff(3) = (-3*T/dlt3 +(one+1/dlt)*Gmn-(T+2.5*V)/dlt2 - (1.5*P11 - 0.5*T+V)/dlt)
      redcoeff(4) = (3*T/dlt3 + 2.5*V/dlt2 - Gmn/dlt + 1.5*P11/dlt)
      redcoeff(8) = (-T/dlt2 - V/(2*dlt) + Gmn)

    !! masses: (m0,0,0)
    else if (masses_ind == 1) then
      z0  = msq(0)/psq
      z01 = one + z0
      redcoeff(1) = (T*(z01**2)/dlt2 + V*z01/dlt + P11)*psq
      redcoeff(2) = z01*T/dlt2 + (V - T/2)/dlt
      redcoeff(3) = (-3*z01*T/dlt3- ((z0 + five)*(V - T) + (two + 11*z01)*T)/(2*dlt2)  -&
      ((5*z0 + 8)*T + (z0 + 10*one)*P11 - (2*z0 + 7*one)*(P11 - V + T))/(2*dlt)        -&
      (V - T/2 + 3*P11/2) + Gmn*(z01/dlt + z01 + one + dlt))/(one + dlt)
      redcoeff(4) = 3*z01*T/dlt3+(z0+five)*V/(2*dlt2)+(3*one-z0)*P11/(2*dlt) -Gmn*z01/dlt
      redcoeff(5) = ((T/dlt2 + (T + V)/(2*dlt) + P11/2)*z0)/(one + dlt)
      redcoeff(8) = ((-T/dlt2 - (2*T+V)/(2*dlt) - V/2 + (one + dlt)*Gmn))/(one + dlt)

    !! masses: (0,m1,0)
    else if (masses_ind == 2) then
      z1 = msq(1)/psq
      redcoeff(1) = psq*((6*T*z1**2)/dlt4+P11+2*Gmn*z1-(3*z1*(2*T*(one-z1)-V*z1))/dlt3 -&
      (4*P11*z1 + 4*Gmn*(z1-one)*z1+2*V*(-one+2*z1))/(2*dlt)-(4*Gmn*z1**2-2*(-2*V*(two -&
      z1)*z1 + P11*z1**2 + T*(one - 4*z1 + z1**2)))/(2*dlt2))
      redcoeff(3) = Gmn + (T/2 -1.5*P11 - V +Gmn*(one-2*z1))/dlt+6*T*z1/dlt4+(-(V*(0.5 +&
      2*(one - z1)))-T*(one-z1)-2*Gmn*z1+P11*z1)/dlt2 - (3*(T*(one - 2*z1) - V*z1))/dlt3
      redcoeff(4) = (-6*T*z1)/dlt4 +(V*(five - z1) + 4*Gmn*z1 - 2*P11*z1)/(2*dlt2)     +&
      (3*(T*(one - z1) - V*z1))/dlt3 + (4*Gmn*(z1 - one) + 2*P11*(3*one + z1))/(4*dlt)
      redcoeff(5) =z1*(-3*T/dlt3 + (T*(one - z1) - 1.5*V*z1)/(dlt2*z1) + (Gmn +(-T+2*V -&
      P11*z1)/(2*z1))/dlt)
      redcoeff(8) = Gmn - T/dlt2 - (T + 2*V)/(4*dlt)

    !! masses: (0,0,m2)
    else if (masses_ind == 3) then
      dlt2 = dlt**2
      z2  = msq(2)/psq
      redcoeff(1) = psq*(P11 + (6*T*z2**2)/dlt4 + (V + 2*P11*z2 -2*Gmn*z2)/dlt         +&
      (T + 4*V*z2 - 2*Gmn*z2**2 + P11*z2**2)/dlt2 + (6*T*z2 + 3*V*z2**2)/dlt3)
      redcoeff(3) = Gmn - (6*T*z2)/dlt4 + ((P11 + T - V)*z2)/(2*(dlt + one))-(3*(V*z2  +&
      T*(one + z2)))/dlt3 + (4*Gmn*z2 - 2*P11*z2+2*T*(-one+z2)-V*(five+3*z2))/(2*dlt2) +&
      (4*Gmn*(one + z2) + 2*(P11*(-3*one - z2) +T*(one - z2) + V*(-two + z2)))/(4*dlt)
      redcoeff(4) = (1.5*P11 - Gmn)/dlt + 6*T*z2/dlt4+(V*5-4*Gmn*z2+2*P11*z2)/(2*dlt2) +&
      (3*(T + V*z2))/dlt3
      redcoeff(5) = 3*T*z2/dlt3 + ((V - P11 - T)*z2)/(2*(dlt + one)) +(-2*T*(one-z2)   -&
      2*V*(z2 - two) - 4*Gmn*z2 + 2*P11*z2)/(4*dlt) + (2*T*(one - z2) + 3*V*z2)/(2*dlt2)
      redcoeff(8) = Gmn - T/dlt2 + (T - 2*V)/(4*dlt)

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      z = msq(0)/psq
      redcoeff(1) = psq*((one-2*z)*T/dlt2 + (one-z)*V/dlt + 2*z*Gmn + P11)
      redcoeff(2) = T/dlt2 + (V - T/2)/dlt
      redcoeff(3) = (-3*T/dlt3 - (2.5*V +T)/dlt2 + Gmn/dlt - (3*P11-T+2*V)/(2*dlt) + Gmn)
      redcoeff(4) = (3*T/dlt3 + 2.5*V/dlt2 + (1.5*P11-Gmn)/dlt)
      redcoeff(8) = (-T/dlt2 - V/(2*dlt) + Gmn)

    !! masses: (0,m,m)
    else if(masses_ind == 5) then
      z1 = msq(1)/psq
      redcoeff(1) = psq*((one+z1**2 - 4*z1)*T/dlt2 + (one - 2*z1)*V/dlt + 2*z1*Gmn + P11)
      redcoeff(2) = (one - 2*z1)*T/dlt2 + (T + (two - z1)*(V - T))/(2*dlt) + z1*(V - T -&
      P11)/(2*(one + dlt))
      redcoeff(3) = -3*(one-z1)*T/dlt3 + ((5*z1-7*one)*T + (z1-five)*(V - T))/(2*dlt2) +&
      (one - z1)*Gmn/dlt - ((3*one + z1)*P11 - (z1 - two)*(V - T) + T)/(2*dlt) + Gmn   +&
      z1*(P11 - V + T)/(2*(one+dlt))
      redcoeff(4) = 3*(one - z1)*T/dlt3 + (five - z1)*V/(2*dlt2) - Gmn*(one - z1)/dlt  +&
      (3*one+z1)*P11/(2*dlt)
      redcoeff(8) = -(one + z1)*T/dlt2 - ((one + z1)*V - z1*T)/(2*dlt) + Gmn+0.5*z1*(V -&
      T - P11)/(one + dlt)

    !! masses: (m,0,m)
    else if (masses_ind == 6) then
      z0 = msq(0)/psq
      z01 = one + z0
      redcoeff(1) = psq*(P11 + (6*T*z0**2)/dlt4 + (3*z0*(2*T*z01+V*z0))/dlt3-(-2*V*z01 -&
      4*P11*z0 + 4*Gmn*z01*z0)/(2*dlt) - (4*Gmn*z0**2 -2*(T*((one + z0)**2 ) +2*V*(one +&
      z01)*z0 + P11*z0**2))/(2*dlt2))
      redcoeff(3) = Gmn + (-6*P11+2*T-4*V+4*Gmn*(-one+2*z01))/(4*dlt) -(6*T*z0)/dlt4   +&
      (-2*T*z01-V*(one +4*z01) + 4*Gmn*z0 - 2*P11*z0)/(2*dlt2) - (3*(T*(2*z01 - one)   +&
      V*z0))/dlt3
      redcoeff(4) = (0.5*P11*(3*one - z0) - Gmn*z01)/dlt +6*T*z0/dlt4+(V*(two+0.5*z01) -&
      2*Gmn*z0 + P11*z0)/dlt2 + (3*(T*z01 + V*z0))/dlt3
      redcoeff(5) = -(-3*z0*T/dlt3 + (z0*Gmn - 0.5*(z0*P11 + 2*V - T))/dlt - (2*T*z01  +&
      3*V*z0)/(2*dlt2))
      redcoeff(8) = Gmn - T/dlt2 + (0.25*T - 0.5*V)/dlt

    !! masses: (m,m,0)
    else if (masses_ind == 7) then
      z0 = msq(0)/psq
      redcoeff(1) = psq*((6*T*z0**2)/dlt4 + 2*(P11/2 + Gmn*z0) + (V*(one - z0) +2*(Gmn -&
      P11)*z0)/dlt + (3*z0*(V*z0-2*T))/dlt3 + (T*(one - 2*z0) - 4*V*z0 - 2*Gmn*z0**2   +&
      P11*z0**2)/dlt2)
      redcoeff(3) = Gmn + (6*T*z0)/dlt4 + ((V-P11 - T)*z0)/(2*(dlt +one))+(3*(-(T*(one -&
      z0)) + V*z0))/dlt3 + (4*Gmn*(one - z0) - 2*P11*(3*one - z0) + 2*T*(one + z0)     -&
      2*V*(two + z0))/(4*dlt) + (-2*Gmn*z0 + P11*z0 - T*(one+z0)+V*(1.5*z0-2.5*one))/dlt2
      redcoeff(4) = (1.5*P11-Gmn)/dlt - (6*T*z0)/dlt4 + (2.5*V + (2*Gmn -P11)*z0)/dlt2 +&
      (3*(T - V*z0))/dlt3
      redcoeff(5) = (-3*T*z0)/dlt3 + (-T/2 + V + Gmn*z0)/dlt +(T*(one+z0/2)-V*z0)/dlt2 -&
      (z0*(V - T + P11*dlt))/(2*dlt2*(dlt + one))
      redcoeff(8) = Gmn - T/dlt2 - (0.25*T + V/2)/dlt

    !! masses: (m0,m1,m1)
    else if(masses_ind == 8) then
      z0 = msq(0)/psq
      z1 = msq(1)/psq
      z01 = one + z0 - z1
      z10 = z0 - z1
      redcoeff(1) = psq*(P11 + (T*(one + 2*(z0 - 2*z1) + z10**2))/dlt2 + (V*(one + z0  -&
      2*z1))/dlt + 2*Gmn*z1)
      redcoeff(2) = ((one+z0-2*z1)*T/dlt2 + (T*(one/2 + z0 -1.5*z1)+V*(one-z1/2))/dlt + &
      0.5*(2*V - T - z1*P11))/(one+dlt)
      redcoeff(3) =(-1.5*P11 + T/2 - V - 3*T*z01/dlt3 + Gmn*(two + dlt +z0+z01/dlt-z1) -&
      ((V - T)*(7*one + 2*z10) + T*(8*one + 5*z10) + P11*(3*one - z10))/(2*dlt) + ((T  -&
      V)*(five + z10) + T*(-13*one - 11*z0 + 11*z1))/(2*dlt2))/(dlt + one)
      redcoeff(4) =((T*(11*one + 7*z10) + (V -T)*(five + z10))/(2*dlt2)+(3*T*z01)/dlt3 -&
      (one + one/dlt)*Gmn*z01 + (P11*(3*one - z0 + z1))/2 + (((P11 + V)*(5*one+z10))/2 -&
      P11*z01)/dlt)/(dlt + one)
      redcoeff(5) = z0*(P11/two + T/dlt2 + 0.5*(T + V)/dlt)/(dlt + one)
      redcoeff(8) = (dlt*Gmn -T*(one+z1)/dlt2 -((3*one+2*z1)*T+(V-T)*(one+z1))/(2*dlt) +&
      Gmn - 0.5*(V+z1*P11))/(dlt + one)

    !! masses: (m0,m0,m2)
    else if(masses_ind == 9) then
      z0 = msq(0)/psq
      z2 = msq(2)/psq
      z02 = z0 - z2
      redcoeff(1) = psq*(P11 + 2*Gmn*z0 + (6*T*z02**2)/dlt4 + (V*(one - z0) +2*Gmn*z02 -&
      2*P11*z02)/dlt + (3*z02*(-2*T + V*z02))/dlt3 +(T*(one-2*z0)-4*V*z02-2*Gmn*z02**2 +&
      P11*z02**2)/dlt2)
      redcoeff(2) = (-3*T/dlt3 + (P11 + T - V)/(2*(dlt + one)) + (Gmn - P11/2 -(T*(one +&
      z0/z02**2))/2 + V*(0.5*one + one/z02))/dlt + (-3*V+2*T*(one+one/z02))/(2*dlt2))*z02
      redcoeff(3) = Gmn + 6*T*z02/dlt4 + ((V-P11-T)*z02)/(2*(dlt+one))+(3*(T*(z02-one) +&
      V*z02))/dlt3 + (Gmn*(one - z02) - V*(one + z02/2) + (P11*(z02-3*one))/2 +(T*(one +&
      z02))/2)/dlt + (-4*Gmn*z02 + 2*P11*z02 - 2*T*(one + z02) +V*(3*z02-5*one))/(2*dlt2)
      redcoeff(4) = (-Gmn + 1.5*P11)/dlt -6*T*z02/dlt4+(2.5*V+2*Gmn*z02-P11*z02)/dlt2  +&
      (3*(T - V*z02))/dlt3
      redcoeff(6) = (z2/2)*T/(dlt*z02)
      redcoeff(8) = Gmn - T/dlt2 - (V/2 + (T*(z0 + z2))/(4*z02))/dlt

    !! masses: (m0,m1,m0)
    else if(masses_ind == 10) then
      z0 = msq(0)/psq
      z1 = msq(1)/psq
      z01 = one + z0 - z1
      z10 = z0 - z1
      redcoeff(1) = psq*((2*P11 + 4*Gmn*z1)/2 + (6*T*z10**2)/dlt4 + (3*z10*(2*T*z01 +   &
      V*z10))/dlt3-(2*V*(z1 -z01) - 4*P11*z10 + 4*Gmn*z01*z10)/(2*dlt) - (4*Gmn*z10**2 -&
      2*(T*((one + z0)**2 - 2*(two + z0)*z1 + z1**2) + 2*V*(one + z01)*z10 +            &
      P11*z10**2))/(2*dlt2))
      redcoeff(3) = Gmn + (-6*P11+2*T-4*V+4*Gmn*(-one+2*z01))/(4*dlt) -(6*T*z10)/dlt4  +&
      (-2*T*z01-V*(one +4*z01) + 4*Gmn*z10 - 2*P11*z10)/(2*dlt2) - (3*(T*(2*z01 - one) +&
      V*z10))/dlt3
      redcoeff(4) = (0.5*P11*(3*one - z10) - Gmn*z01)/dlt + (6*T*z10)/dlt4 + (V*(two   +&
      0.5*z01) - 2*Gmn*z10 + P11*z10)/dlt2 + (3*(T*z01 + V*z10))/dlt3
      redcoeff(6) = (-3*T/dlt3 + (Gmn - 0.5*(P11 +T*z1/z10**2+2*V/z10))/dlt - (2*T*z01 +&
      3*V*z10)/(2*dlt2*z10))*z1
      redcoeff(5) = -z0*(redcoeff(6)/z1 + (0.5*T/z10)/dlt)
      redcoeff(8) = Gmn - T/dlt2 + (T*(z0 + z1)/(4*z10) - 0.5*V)/dlt

    !! masses: (0,m1,m2)
    else if (masses_ind == 11) then
      z1 = msq(1)/psq
      z2 = msq(2)/psq
      z12 = z1 - z2

      redcoeff(1) = psq*(P11 + 2*Gmn*z1 + 6*T*z12**2/dlt4 + (3*z12*(2*T*(-one + z1)    +&
      V*z12))/dlt3 - (-2*V*(one - 2*z1) + 4*P11*z12 + 4*Gmn*(z1-one)*z12)/(2*dlt)      -&
      (4*Gmn*z12**2 - 2*(T*(one-4*z1 + z1**2) + 2*V*(z1-two)*z12 + P11*z12**2))/(2*dlt2))
      redcoeff(3) = Gmn + 6*T*z12/dlt4 + (3*(V*z12 + T*(z1 -one+z12)))/dlt3+(2*T*(2*z1 -&
      z12-one) - 4*Gmn*z12 + 2*P11*z12 + V*(z1 + 3*z12-five))/(2*dlt2) + ((P11 + T     -&
      V)*z2)/(2*(dlt + one)) + (Gmn*(one - z1 - z12) + (-(T*(z1 -one -z12))-V*(two-z2) -&
      P11*(3*one + z2))/2)/dlt
      redcoeff(4) = (Gmn*(z1-one)+(P11*(3*one+z1))/2)/dlt-6*T*z12/dlt4+(-(V*(z1-five)) +&
      (4*Gmn - 2*P11)*z12)/(2*dlt2) + (3*(T*(one - z1) - V*z12))/dlt3
      redcoeff(5) = z1*(-3*T/dlt3 + (Gmn + (-P11 - T*z1/z12**2 + 2*V/z12)/2)/dlt       -&
      (T*(z1-one) + 1.5*V*z12)/(dlt2*z12))
      redcoeff(6) = (3*T/dlt3 + (V-P11 - T)/(2*(dlt + one)) + (-(V*(one + z12/2))      -&
      Gmn*z12 + P11*z12/2 + (T*(one + z1/z12 + z12))/2)/(dlt*z12) + (3*V*z12 -2*T*(one -&
      z2))/(2*dlt2*z12))*z2
      redcoeff(8) = Gmn - T/dlt2 - (V/2 + (T*(z1 + z2))/(4*z12))/dlt

    !! masses: (m0,0,m2)
    else if (masses_ind == 12) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z2  = msq(2)/psq
      z02 = one + z0
      redcoeff(1) = psq*(P11 + (6*T*z2**2)/dlt4 + (V*z02 + 2*P11*z2 -2*Gmn*z02*z2)/dlt +&
      (T*z02**2 + 2*V*(2 + z0)*z2 - 2*Gmn*z2**2 + P11*z2**2)/dlt2 + (6*T*z02*z2        +&
      3*V*z2**2)/dlt3)
      redcoeff(3) = Gmn - (6*T*z2)/dlt4 + ((P11 + T - V)*(-z0 + z2))/(2*(dlt + one))   -&
      (3*(V*z2 + T*(z02 + z2)))/dlt3 + (4*Gmn*z2 - 2*P11*z2 + 2*T*(-one - 2*z0 + z2)   -&
      V*(five + z0 + 3*z2))/(2*dlt2) + (4*Gmn*(z02 + z2) + 2*(P11*(-3*one + z0 - z2)   +&
      T*(z02 - z2) + V*(-one - z02 + z2)))/(4*dlt)
      redcoeff(4) = (-2*P11*(z0 -3*one) - 4*Gmn*z02)/(4*dlt) + (6*T*z2)/dlt4           +&
      (V*(five + z0) - 4*Gmn*z2 + 2*P11*z2)/(2*dlt2) + (3*(T*z02 + V*z2))/dlt3
      redcoeff(5) = (T/dlt2 + (P11 + T - V)/(2*(dlt + one)) + (-T + V)/(2*dlt))*z0
      redcoeff(6) = (3*T*z2)/dlt3 + ((V - P11 - T)*z2)/(2*(dlt + one)) +(-2*T*(one-z2) -&
      2*V*(z2 - two) - 4*Gmn*z2 + 2*P11*z2)/(4*dlt) + (2*T*(z02 - z2) + 3*V*z2)/(2*dlt2)
      redcoeff(8) = Gmn - T/dlt2 + (T - 2*V)/(4*dlt)

    !! masses: (m0,m1,0)
    else if (masses_ind == 13) then
      dlt2 = dlt**2
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z01 = one + z0 - z1
      redcoeff(1) = psq*(P11 + (6*T*z1**2)/dlt4 + 2*Gmn*z1*(one + z01/dlt - z1/dlt2)   +&
      (V*(z01 - z1)- 2*P11*z1)/dlt +(3*z1*(-2*T*z01 +V*z1))/dlt3 +(-2*V*(one + z01)*z1 +&
      P11*z1**2 + T*((one + z0)**2 - 2*(two + z0)*z1 + z1**2))/dlt2)
      redcoeff(2) = (-3*T*z1)/dlt3 + (-T/2 + V + Gmn*z1)/dlt + (T*z01 - V*z1)/dlt2
      redcoeff(3) = (T/2 -1.5*P11 - V + (6*T*z1)/dlt4 +Gmn*(one+dlt+z01+(z01-3*z1)/dlt -&
      z1 - (2*z1)/dlt2) - (3*(T*(z01 - 3*z1) - V*z1))/dlt3 + (T*(-one - 3*z0 + 2*z1)   +&
      P11*(-3*one + z0 + 2*z1) + V*(-7*one-2*z0+4*z1))/(2*dlt)+(-(V*(five+z0-10*z1))/2 +&
      P11*z1 + T*(7*z1 - 5*z0 -4*one))/dlt2)/(dlt + one)
      redcoeff(4) = (P11*(two - z01/2))/dlt - Gmn*z01/dlt-6*T*z1/dlt4+(V*(two + z01/2) +&
      2*Gmn*z1 - P11*z1)/dlt2 + (3*(T*z01 - V*z1))/dlt3
      redcoeff(5) = ((P11 + 2*T/dlt2 + (T + V)/dlt)*z0)/(2*(dlt + one))
      redcoeff(6) = -((P11/dlt + V/dlt2)*z1)/2
      redcoeff(8) = Gmn - (one/dlt2 + one/(4*dlt))*T - V/(2*dlt)

    !! masses: (m0,m1,m2)
    else if (masses_ind == 14) then
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z2  = msq(2)/psq
      z12 = z1 - z2
      z01 = one + z0 - z1

      redcoeff(1) = psq*((P11*(dlt - z12)**2)/dlt2 - (2*Gmn*(- (dlt2*z1) - dlt*z01*z12 +&
      z12**2))/dlt2 + (V*(dlt2*(z01 - z1) - 2*dlt*(one + z01)*z12 + 3*z12**2))/dlt3    +&
      (T*(dlt2*(-2*(one + z0) + 2*z01 + z01**2) + 6*z12*(-(dlt*z01) + z12)))/dlt4)

      redcoeff(3) = (Gmn*(dlt*(one + dlt + z0) - 2*(dlt + one)*z1+(two +dlt)*z2))/dlt2 +&
      (P11*(dlt*(-3*(one+dlt) + z0) + 2*(dlt + one)*z1 - (two+3*dlt)*z2))/(2*dlt2*(dlt +&
      one)) - (V*(dlt*(five + 7*dlt +2*dlt2+z0+2*dlt*z0)-2*(3*one + 5*dlt + 2*dlt2)*z1 +&
      (6*one + 9*dlt + 2*dlt2)*z2))/(2*dlt3*(dlt +one))+(T*(dlt*(-6*one - 8*dlt - dlt2 +&
      dlt3 - (6*one + 10*dlt + 3*dlt2)*z0) + 2*(6*one + 12*dlt + 7*dlt2 + dlt3)*z1     +&
      (-12*one - 18*dlt - 4*dlt2 + dlt3)*z2))/(2*dlt4*(dlt + one))

      redcoeff(4) = (Gmn*(-(dlt*(one+z0))+(two+dlt)*z1-2*z2))/dlt2+(P11*(-(dlt*(-3*one +&
      z0))+(-two + dlt)*z1 + 2*z2))/(2*dlt2)+(3*T*(dlt+dlt*z0-(two+dlt)*z1+2*z2))/dlt4 +&
      (V*(dlt*(five + z0) - (6*one + dlt)*z1 + 6*z2))/(2*dlt3)

      redcoeff(5) = ((dlt2*P11 + (two + dlt)*T + dlt*V)*z0)/(2*dlt2*(dlt + one))
      redcoeff(6) = z1*(Gmn/dlt - P11/(2*dlt) + (V*(2*dlt - 3*z12))/(2*dlt2*z12)       +&
      (T*(-2*(3*one + dlt)*z1**2 - 2*z2*(dlt + dlt*z0 + 3*z2) + z1*(dlt*(two-dlt+2*z0) +&
      2*(6*one + dlt)*z2)))/(2*dlt3*z12**2))

      redcoeff(7) = z2*(-(Gmn/dlt) + P11/(2*dlt*(dlt + one)) + (V*(-2*dlt*(dlt + one)  +&
      (3*one + 2*dlt)*z1 - (3*one +2*dlt)*z2))/(2*dlt2*(dlt+one)*z12)+(T*((6*(one+dlt) +&
      dlt2)*z1**2 - 2*z1*(dlt*(dlt + one)*(one - dlt + z0) + (6*one + 5*dlt)*z2)       +&
      z2*(dlt*(dlt + one)*(two - dlt + 2*z0) + (6*one + 4*dlt-dlt2)*z2)))/(2*dlt3*(dlt +&
      one)*z12**2))

      redcoeff(8) = Gmn - V/(2*dlt) - (T*(2*dlt*z1 - (-4*one + dlt)*z12))/(4*dlt2*z12)

    end if

  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-3 part ---------------------*
  !------------------------------------------------------------------------*
  else if (r == 3) then
  dlt2 = dlt**2
  dlt3 = dlt**3
  dlt4 = dlt**4
  dlt5 = dlt**5
  dlt6 = dlt**6

  n = 16
  do i = 1, 4
    do j = i, 4
      do k = j, 4
        p1p1p1(n) = p1(i)*p1(j)*p1(k)
        p1_p2(n)  = p12(i)*p12(j)*p12(k)
        p1p12p12(n) = p12(i)*p12(j)*p1(k) + p12(i)*p12(k)*p1(j) + p12(k)*p12(j)*p1(i)
        p1p1p12(n)  = p1(i)*p1(j)*p12(k) + p1(i)*p1(k)*p12(j) + p1(k)*p1(j)*p12(i)
        n = n + 1
      end do
    end do
  end do

  P111 = SUM(p1p1p1(16:35)*G_in(16:35))/psq

  T  =   SUM(p1_p2(16:35)*G_in(16:35))/psq
  U  = - SUM(p1p12p12(16:35)*G_in(16:35))/psq
  V1 = - SUM(p1p1p12(16:35)*G_in(16:35))/psq
  W0 = U + T

  !! g^{\mu\nu} p^{\rho}_1 + g^{\rho\mu} p^{\nu}_1 + g^{\nu\rho} p^{\mu}_1
  Gp1 = 2*((2*G_in(17) - G_in(24))*p1(1) + (2*G_in(20) - G_in(30))*p1(2)    + &
        (G_in(21) - 2*G_in(33))*p1(3) + (G_in(22) - 2*G_in(34))*p1(4))

  !! g^{\mu\nu} p^{\rho}_2 + g^{\rho\mu} p^{\nu}_2 + g^{\nu\rho} p^{\mu}_2
  Gp2 = 2*((2*G_in(17) - G_in(24))*p2(1) + (2*G_in(20) - G_in(30))*p2(2)    + &
        (G_in(21) - 2*G_in(33))*p2(3) + (G_in(22) - 2*G_in(34))*p2(4))

  !! g^{\mu\nu} (p^{\rho}_1 - p^{\rho}_2)
  Gp12 = 2*((2*G_in(17) - G_in(24))*p12(1) + (2*G_in(20) - G_in(30))*p12(2) + &
        (G_in(21) - 2*G_in(33))*p12(3) + (G_in(22) - 2*G_in(34))*p12(4))

  Gp = (Gp1+Gp2)/12
  Gm = Gp12/12

  !! masses: (0,0,0)
    if (masses_ind == 0) then
      redcoeff(1) = psq*(-T/dlt3 + U/dlt2 + V1/dlt - P111)
      redcoeff(2) = (-(one - dlt/2 + dlt2/3)*T/dlt3 - (dlt-two)*U/(2*dlt2) + V1/dlt)
      redcoeff(3) = (11*T/3)/dlt4 + (T-10*U/3)/dlt3 - (Gm + T/2 + U + 17*V1/6)/dlt2    +&
      (2*T + 3*U - 6*V1 + 11*P111 - Gp1)/(6*dlt) - Gp
      redcoeff(4) = -(11*T/3)/dlt4+(10*U/3)/dlt3+(Gm+17*V1/6)/dlt2+(Gp1 -11*P111)/(6*dlt)
      redcoeff(8) = (5*T/3)/dlt3 - (T+4*U)/(3*dlt2) + (T/18 + (U - 5*V1)/6 -Gm)/dlt    -&
      (5*Gp1+2*Gp2)/36

    !! masses: (m0,0,0)
    else if (masses_ind == 1) then
      z   = msq(0)/psq
      z2  = z**2
      z1  = one + z
      zd  = z1/dlt
      z12 = z1**2
      T3  = T/3
      T6  = T/6
      U3  = U/3
      V16 = V1/6
      P16 = P111/6

      redcoeff(1) = (-T*zd**3 + U*zd**2 + V1*zd - P111)*psq
      redcoeff(2) = -T*z12/dlt3 + (T/2 + U)*z1/dlt2 - (T3 + U/2 - V1)/dlt
      redcoeff(3) =(-dlt2*Gp + (dlt*(11*P111 + 2*T +3*U -6*V1-3*Gp*(6*one+z)-3*Gm*(two +&
      3*z)))/6 + ((-(U*(10*one + z))+5*T*(5*one+6*z))*z1)/(3*dlt3)+(11*T*z12)/(3*dlt4) +&
      (22*P111 + T - 29*V1 - 10*P111*z + 5*(P111 -T+V1-U)*z-4*U*z-14*V1*z+3*Gp*(-6*one -&
      2*z + z2) - 3*Gm*(6*one + 10*z + 3*z2))/6 + (11*P111 + 2*T - 29*U-40*V1-5*P111*z +&
      7*T*z - 45*U*z - 13*V1*z + 2*P111*z2 + 11*T*z2-6*U*z2 + 3*V1*z2 + 3*Gp*(-two - z +&
      z2) - 3*Gm*(6*one + 11*z + 5*z2))/(6*dlt) + (-46*U -17*V1-58*U*z-4*V1*z-6*Gm*z12 -&
      6*U*z2 + V1*z2 + T*(31 + 78*z + 51*z2))/(6*dlt2))/(one + dlt)**2
      redcoeff(4) = ((Gm + Gp)*(2*one - z)*z1)/(2*dlt) + (U*(10*one + z)*z1)/(3*dlt3)  +&
      (Gm*z12)/dlt2 - (11*T*z12)/(3*dlt4) - (V1*(-17*one - 4*z + z2))/(6*dlt2)         -&
      (P111*(11*one - 5*z + 2*z2))/(6*dlt)
      redcoeff(5) = -5*z1*T3/dlt3 + ((z+4*one)*U3-(17*z + 13*one)*T6)/dlt2 + ((5.5*one +&
      2*z)*U3 - (z - five)*V16 + Gm*z1 - T6*(5*z + one))/dlt + (z*Gp)/2+P16*(4*z-five) +&
      Gm*(1.5*z + two) + T3 + U/2 + (five + z)*V16 + dlt*((z*Gp)/2 + (2*z - five)*P16  +&
      Gm*(z/2 +one))
      redcoeff(5) = z*redcoeff(5)/((one+ dlt)**2)
      redcoeff(8) = 5*z1*T3/dlt3 + ((6*z+4*one)*T3 -(z+4*one)*U3)/dlt2 +((z-five/3)*T6 -&
      (1.5*z + 3.5*one)*U3 + (z-five)*V16 - Gm*z1)/dlt + (T3 + U -7*Gp-z*P111-Gm*(6*z  +&
      9*one) - 5*V1)/6 - (5*Gp1 + 2*Gp2)*dlt/36
      redcoeff(8) = redcoeff(8)/(one + dlt)

    !! masses: (0,m1,0)
    else if (masses_ind == 2) then
      z1  = msq(1)/psq
      z11 = one - z1
      redcoeff(1) = psq*(-P111 - Gp1*z1/2 +(20*T*z1**3)/dlt6-(10*z1**2*(T*(3*one-4*z1) +&
      W0*z1))/dlt5+(Gp1*z1*(2*z1-one)+2*(V1*(one-3*z1)+3*P111*z1)-Gp12*z1*z11)/(2*dlt) +&
      (Gp1*(two - z1)*z1**2 - Gp12*z1*(z11**2 - 2*z1) - 2*(-(z1*(-6*V1*z11-3*P111*z1)) -&
      U*(-two + 3*z11**2)))/(2*dlt2) - (z1*(Gp12*z1**2-2*(6*T*(-z1+z11**2)-z1*(2*V1*z1 -&
      3*U*(one+2*z11)))))/dlt4+(3*Gp12*z11*z1**2-Gp1*z1**3-2*(T*(one - 8*z1+z1**2)*z11 -&
      z1*(z1*(3*V1*(3*one - z1) + P111*z1) - 9*U*((-2*z1**2)/3 + z11**2))))/(2*dlt3))

      redcoeff(3) = (Gp12-2*Gp1)/12+(20*T*z1**2)/dlt6-(10*z1*(T*(two-3*z1)+U*z1))/dlt5 +&
      (11*P111 + 2*T + 3*U - 6*V1 + 1.5*Gp12*z1 + Gp1*(-one + 4.5*z1))/(6*dlt)         +&
      (-3*Gp12*z1**2 + 3*z1*(U*(one + 12*z11) - 4*V1*z1) + T*(11*one - 63*z1           +&
      36*z1**2))/(3*dlt4) + (-3*Gp1*z1**2 - 3*Gp12*z1*(3*z1-two)+6*(V1*(7*one-3*z1)*z1 +&
      P111*z1**2 + (U*(36*z1 -10*one - 9*z1**2))/3 + T*(one - 3*z1 + z1**2)))/(6*dlt3) +&
      (-2*U*(6*one - 9*z1) - 30*P111*z1 - 3*Gp1*z1*(2*z1-3*one) + V1*(54*z1-34*one)    -&
      Gp12*(one - 15*z1 + 6*z1**2) - 6*T*z11)/(12*dlt2)

      redcoeff(4) = (-20*T*z1**2)/dlt6 - (Gp1*(5*z1 -two + z1**2) + 2*P111*(11 + 7*z1  +&
      2*z1**2))/(12*dlt) - (U*(13*one - 7*z1)*z1 - (Gp12 + 4*V1)*z1**2 + (T*(11*one    +&
      z1*(11*z1-38*one)))/3)/dlt4 + (3*Gp1*z1**2 + 2*(3*Gp12*(z1 - one)*z1             +&
      3*z1*(-(V1*(7*one - z1)) -P111*z1) + U*(10*one - 19*z1 + z1**2)))/(6*dlt3)       +&
      (3*Gp1*(z1 - 3*one)*z1 +Gp12*(one - 10*z1 + z1**2) + 2*(3*P111*z1*(5*one + z1)   -&
      V1*(8*z1 + z1**2 -17*one)))/(12*dlt2) + (10*z1*(U*z1 + 2*T*z11))/dlt5

      redcoeff(5) = z1*(-10*T*z1/dlt5 + (-2*U*(17*one - 8*z1)*z1 + 3*Gp12*z1**2        +&
      12*V1*z1**2-6*T*(one-3*z1+z1**2))/(6*dlt3*z1)+(-Gp12*z1/4+(Gp1*(z1-4*one)*z1)/12 +&
      (3*V1 - T -1.5*U+P111*z1*(3.5*one+z1))/3)/(dlt*z1)+(5*(3*U*z1+5*T*z11))/(3*dlt4) +&
      (-2*V1*(19*one - z1)*z1 + 6*T*z11 + 6*U*(2*z11 - z1) + z1*(3*Gp1*z1 - 6*P111*z1  -&
      5*Gp12*z11))/(12*dlt2*z1))

      redcoeff(8) = (2*Gp12 - 7*Gp1)/36 - (10*T*z1)/(3*dlt4) - (-2*Gp12*z1 - 8*V1*z1   +&
      T*(one + 3*z1) - U*(-16*one + 9*z1))/(12*dlt2) + (-8*T - 3*(U + 10*V1)+ 3*Gp1*z1 -&
      6*P111*z1 - 1.5*Gp12*(-z1 + 2*z11))/(36*dlt) + (5*(2*U*z1 + T*(2*z11-z1)))/(6*dlt3)

    !! masses: (0,0,m2)
    else if (masses_ind == 3) then
      z2 = msq(2)/psq
      redcoeff(1) = psq*(- P111 - (20*T*z2**3)/dlt6 - (10*z2**2*(3*T - U*z2))/dlt5     +&
      (Gp1*z2 - 2*(-V1+3*P111*z2))/(2*dlt)+(Gp12*z2+Gp1*2*z2**2 - 2*(-U + 3*z2*(-2*V1  +&
      P111*z2)))/(2*dlt2) + (z2*(Gp12*z2**2 -2*(6*T - z2*(9*U + 2*V1*z2))))/dlt4       +&
      (3*Gp12*z2**2 + Gp1*z2**3 - 2*(T + z2*(-9*U + z2*(-9*V1 + P111*z2))))/(2*dlt3))
      redcoeff(3) = (-2*Gp1 + Gp12)/12 + ((T - P111+U-V1)*z2**2)/(3*(one+dlt)**2)      +&
      (20*T*z2**2)/dlt6 + (10*z2*(-U*z2 + T*(2+z2)))/dlt5+(-Gp12*z2**2-4*V1*z2**2      -&
      U*z2*(13*one + 5*z2) + (T*(11*one + 5*(five-2*z2)*z2))/3)/dlt4 +(Gp1*(-two+z2**2 -&
      4*z2)-Gp12*(z2**2 - 2*z2) + 2*V1*(-6*one +7*z2 - 2*z2**2) + 2*U*(3*one - 7*z2    +&
      3*z2**2) + 2*T*(two - 7*z2 + 4*z2**2) + 2*P111*(11*one +7*z2 - z2**2 ))/(12*dlt) +&
      (-((Gp1 - Gp12)*z2**2) -2*V1*((7*one -2*z2)*z2) + 2*P111*((-7*one + z2)*z2)      +&
      2*U*((7*one-3*z2)*z2)-2*T*(z2*(4*z2-7*one)))/(12*(dlt + one))                    +&
      (-2*T*(one-3*z2)*(3*one - 2*z2) + 6*P111*z2*(five + z2) -3*Gp1*z2*(3*one+z2)     -&
      Gp12*(one + (five - 2*z2)*z2)-2*U*(6*one - 13*z2 + 5*z2**2) + 2*V1*(-17*one      -&
      19*z2+4*z2**2))/(12*dlt2)+(6*P111*z2**2-6*V1*z2*(7*one + 2*z2)+ 2*T*(3*one-10*z2 +&
      5*z2**2) - 3*z2*(Gp1*z2 + Gp12*(two + z2)) - 2*U*(10*one +17*z2-5*z2**2))/(6*dlt3)
      redcoeff(4) = (-11*P111 + Gp1)/(6*dlt) - 20*T*z2**2/dlt6 +(34*V1+Gp12-30*P111*z2 +&
      9*Gp1*z2)/(12*dlt2)-(10*(2*T*z2-U*z2**2))/dlt5 +(-11*T/3 + U*13*z2 + Gp12*z2**2  +&
      4*V1*z2**2)/dlt4 + (20*U+42*V1*z2 - 6*P111*z2**2 + 3*z2*(2*Gp12 + Gp1*z2))/(6*dlt3)
      redcoeff(5) = (-(T -P111 + U - V1)*z2**2)/(3*(one + dlt)**2) - (10*T*z2**2)/dlt5 +&
      (5*z2*(3*U*z2 +T*(2*z2-five)))/(3*dlt4)+(z2*((Gp1-Gp12)*z2-2*P111*(-7*one + z2)  -&
      2*V1*(-7*one + 2*z2)  + 2*U*(-7*one + 3*z2) + 2*T*(-7*one+4*z2)))/(12*(dlt+one)) +&
      (2*U*(17*one - 5*z2)*z2 + 3*Gp12*z2**2 + 12*V1*z2**2 - 2*T*(3*one - 10*z2        +&
      5*z2**2))/(6*dlt3) + (2*P111*z2*(z2 - 7*one) + 2*U*(-3*(one - z2)**2 + z2)       +&
      2*T*(-two + 7*z2 - 4*z2**2)+2*V1*(6*one - 7*z2 + 2*z2**2) +z2*(Gp1*(4*one-z2)    +&
      Gp12*(z2-two)))/(12*dlt) + (2*V1*(19*one -4*z2)*z2-6*P111*z2**2 + z2*(Gp12*(five -&
      2*z2)+3*Gp1*z2)+T*(6-22*z2 + 12*z2**2)+ 2*U*((z2 - two)*(-3*one + 5*z2)))/(12*dlt2)
      redcoeff(8) = (2*Gp12-7*Gp1)/36+10*T*z2/(3*dlt4)+(T-P111 +U-V1)*z2/(6*(dlt+one)) -&
      (5*(2*U*z2 + T*(z2-two)))/(6*dlt3)+(-6*Gp1*z2+12*P111*z2 + 12*V1*(z2-five)       +&
      3*Gp12*(z2 - two) - 6*U*(2*z2 - five) - 4*T*(3*z2-five))/(72*dlt) + (-2*Gp12*z2  -&
      8*V1*z2 +T*(-7*one +4*z2) + U*(-16*one + 5*z2))/(12*dlt2)

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      z = msq(0)/psq
      ct1 = one - 2*z
      ct2 = 4*z
      T6  = T/6
      U3  = U/3
      V16 = V1/6
      P1116 = (11*P111)/6
      P111z3 = (z*P111)/3

      redcoeff(1) = psq*((T*(6*z - one))/dlt3 - U*(ct2-one)/dlt2 +6*(ct1*V16-z*Gm)/dlt -&
      (z*Gp1)/2 - P111)

      redcoeff(3) = -(16*z-11*one)*T/(3*dlt4)+(4*(7*one-12*z)*T6+2*(ct2-five)*U3)/dlt3 +&
      ((12*z - 13*one)*U3 + (ct2 - 17*one)*V16 + Gm*(8*z - one) + 3*T6*(one-ct2))/dlt2 +&
      (P1116-ct1*Gp+P111z3+Gm*(14*z - two) - ct1*(T6 + 1.5*U3)+(6*z - 23*one)*V16)/dlt +&
      P1116 + 2*Gp*(z-one) - Gp*dlt + T/3 + U/2 - V1 + Gm*(6*z-one)
      redcoeff(3) = redcoeff(3)/(one + dlt)

      redcoeff(4) =(16*z-11*one)*T/(3*dlt4) + 2*(five-ct2)*U3/dlt3 +((one-8*z)*Gm-(ct2 -&
      17*one)*V16)/dlt2 + (ct1*(Gp+Gm)-P1116-P111z3)/dlt

      redcoeff(5) = ((8*z - 3*one)*T)/(3*dlt3) + ((3*one-ct2)*U3 +3*T6*(ct2-one))/dlt2 -&
      (ct2*Gm - ct1*(T6 + 1.5*U3) + 2*(z - 3*one)*V16)/dlt-ct2*Gm+P111z3-2*T6-1.5*U3 + V1
      redcoeff(5) = redcoeff(5)/(one + dlt)

      redcoeff(8) = 16*T6*(one -z)/dlt3 + ((11*one - 12*z)*T6 + (ct2 - 7*one)*U3)/dlt2 +&
      ((3*z - 4*one)*T/9 + (z - five/3)*U - (ct1 + 10*one)*V16 + Gm*(ct2 - one))/dlt   +&
      (Gm*(8*z - 3*one))/2 - 7*one*Gp/6 - P111z3 + ((-Gm - (7*Gp)/3)*dlt)/2 + (7*T6)/3 +&
      2*U3 - 11*V16
      redcoeff(8) = redcoeff(8)/(one + dlt)

    !! masses: (0,m,m)
    else if(masses_ind == 5) then
      z   = msq(1)/psq
      z2  = z**2
      ct1 = one - z
      T3  = T/3
      T6  = T/6
      U3  = U/3
      V16 = V1/6
      P16 = P111/6

      redcoeff(1) = psq*(-(ct1*(ct1**2-6*z)*T)/dlt3 + &
      (U*(3*z*(z-two)+one))/dlt2 + (-6*ct1*z*Gm+(ct1-2*z)*V1)/dlt -(z*Gp1)/2-P111)

      redcoeff(3) = ((11*(ct1**2) - 16*z)*T3)/dlt4 + (5*(6*z2-19*z+five)*T3+(19*z - z2- &
      10*one)*U3)/dlt3 + ((8*z-ct1**2)*Gm + (-138*z + 51*z2 + 31*one)*T6 +(49*z - 3*z2 -&
      23*one)*U3 + (8*z-17*one+z2)*V16)/dlt2+((Gp2*(5*z+z2-two))/12-2*Gm*(z2-15*z+two) +&
      P16*(7*z+11*one+2*z2) + (11*z2-17*z+two)*T6+(37.5*z-3*z2-14.5*one)*U3+(23*z+3*z2 -&
      40*one)*V16)/dlt + Gm*(25*z-1.5*z2-3*one)+Gp*(5*z+z2/2- 3*one)+ P16*(7*z+22*one) +&
      2.5*z*U+(15*z-29*one)*V16 + T6*(7*z + one)+ dlt*(20*P16+ Gp*(2.5*z - 3*one) + T3 -&
      (1.5*P111 - U/2 + V1) + Gm*(7.5*z-one))-Gp*dlt2
      redcoeff(3) = redcoeff(3)/((one + dlt)**2)

      redcoeff(4) = -((11*ct1**2 - 16*z)*T3)/dlt4+(z2-19*z + 10*one)*U3/dlt3 +((ct1**2 -&
      8*z)*Gm -(z2+8*z-17*one)*V16)/dlt2 - ((Gm+Gp)*(2.5*z+z2/2-one) + P16*(7*z + 2*z2 +&
      11*one))/dlt

      redcoeff(5) = (19*z-8*z2-3*one)*T3/dlt3 + ((62*z - 29*z2 - 9*one)*T6 + (z2 -11*z +&
      3*one)*U3)/dlt2 + (z*Gm*(z-five)+(17*z-11*z2-two)*T6 + (2*z2-18.5*z +4.5*one)*U3 +&
      (6*one-7*z-z2)*V16)/dlt + ((7*z + 4*z2)*P16 + z*((z*Gp)/2 + Gm*(1.5*z - 10*one)) -&
      2.5*z*U+(z2-7*z+12*one)*V16-T6*(7*z+one)) + dlt*(-5*z*Gm+(z2*Gp1)/12 +z*P16*(2*z +&
      7*one)-T3-U/2+V1)
      redcoeff(5) = redcoeff(5)/(one + dlt)**2

      redcoeff(8) = (((3*z2 - 19*z + 8*one)*T3)/dlt3 + ((6*z2 - 33*z + 11*one)*T6+(8*z -&
      7*one)*U3)/dlt2 + ((2*z-4*one/3)*T3 + (z-11*one)*V16 + 5*U3*(1.5*z-one)+Gm*(5*z  -&
      one))/dlt + ((-7*Gp)/6-z*P16+Gm*(5*z-1.5*one) +((-5*Gp1 -2*Gp2)*dlt)/36+(7*T6)/3 +&
      2*U3 - 11*V16))/(one + dlt)

    !! masses: (m,0,m)
    else if (masses_ind == 6) then
      z0 = msq(0)/psq
      z01 = one + z0

      redcoeff(1) = psq*(-P111 - 20*T*z0**3/dlt6-(10*z0**2*(T*(3*z01+z0)- W0*z0))/dlt5 +&
      (Gp1*(z0 + z0**2 ) + 2*(V1*z01 - 3*P111*z0))/(2*dlt) + (Gp12*(z0**3 + 2*z0**2    +&
      z0) + Gp1*z0**2*(two + z0) - 2*(-U*(one + z0**2 + 2*z0)+z0*(-2*V1*(3*one + 2*z0) +&
      3*P111*z0)))/(2*dlt2) + (z0*(Gp12*z0**2 - 2*(6*T*(one +z0**2+2*z0)+z0*(-U*(9*one +&
      6*z0)-2*V1*z0))))/dlt4 + (3*Gp12*z01*z0**2+Gp1*z0**3-2*(T*z01**3+z0*(-3*U*(3*one +&
      z0**2 + 4*z0) + z0*(P111*z0 - 3*V1*(3*one + z0)))))/(2*dlt3))

      redcoeff(3) = (-2*Gp1 + Gp12)/12 + ((W0 - P111 - V1)*z0)/(3*(dlt + one)) + (-2*T +&
      4*W0*(1.5*one - z0) + 4*V1*(z0 - 3*one) + 2*P111*(11*one + 2*z0) +Gp12*z0        -&
      Gp1*(two + 5*z0))/(12*dlt) + 20*T*z0**2/dlt6 +(10*z0*(-W0*z0+2*T*(z01+z0)))/dlt5 +&
      (-(Gp12*(one + 7*z0 + 6*z0**2)) -3*Gp1*(2*z0**2 + 3*z0)+V1*(-34*one+8*z0- 54*z0) +&
      30*P111*z0 - 2*W0*(6*one + z0) + 4*T*(1.5*one+z0))/(12*dlt2)+(T*(11*one+72*z0**2 +&
      86*z0)- 3*Gp12*z0**2 -3*z0*(W0*(one+12*z01)+4*V1*z0))/(3*dlt4)+(-3*Gp12*(3*z0**2 +&
      2*z0) - 3*Gp1*z0**2 + 2*(T*(13*one + 12*z0**2 +29*z0)-W0*(10*one+9*z0**2 +28*z0) +&
      3*P111*z0**2 - 3*V1*z0*(7*one + 3*z0)))/(6*dlt3)

      redcoeff(4) = -(Gp1*(z0**2 - two - z0)+2*P111*(11*one + 2*z0**2 -5*z0))/(12*dlt) +&
      (3*Gp1*(z0**2 +3*z0) + Gp12*z01**2 + 2*(-(V1*(z0**2 -17*one-4*z0))+3*P111*(z0**2 -&
      5*z0)))/(12*dlt2) - (20*T*z0**2)/dlt6 - (10*z0*(-W0*z0 + T*(two+3*z0)))/dlt5     -&
      (-((Gp12 + 4*V1)*z0**2) - W0*z0*(13*one + 7*z0) + (T*(11*one + 11*z0**2+22*z0    +&
      3*z0*(13*one + 7*z0)))/3)/dlt4 + (3*Gp1*z0**2 + 2*(3*z0*Gp12*z01 + U*(10*one     +&
      z0**2 + z0*11) - 3*z0*(P111*z0 -V1*(7*one + z0))))/(6*dlt3)

      redcoeff(5) = z0*((P111 + V1 - W0)/(3*(dlt + one))-(10*T*z0)/dlt5 - (5*(-3*W0*z0 +&
      T*(five + 8*z0)))/(3*dlt4) - (-2.5*Gp12*z01*z0 -1.5*Gp1*z0**2 + 3*P111*z0**2     +&
      T*(3*one + 2*z0) - V1*z0*(19*one + z0) - W0*(6*one + z0))/(6*dlt2*z0)            -&
      (-1.5*Gp12*z0 + T*(3/z0 + 18*one+11*z0)-(6*V1*z0 + W0*(17*one + 8*z0)))/(3*dlt3) +&
      (0.5*Gp1*(4*one + z0)*z0- 0.5*Gp12*z0 + P111*z0*(2*z0 - 7*one) + W0*(2*z0-3*one) +&
      2*V1*(3*one - z0) + T)/(6*dlt*z0))

      redcoeff(8) = (2*Gp12 - 7*Gp1)/36 - ((P111 + V1 - W0)*z0)/(3*(dlt + one))        +&
      (10*T*z0)/(3*dlt4) +(5*(-2*W0*z0 + T*(two + 3*z0)))/(6*dlt3)+(T*(4*z0 + 9*one)   -&
      2*Gp12*z0 + (W0*(z0-16*one) - 8*V1*z0))/(12*dlt2) + (6*P111*z0 - 3*Gp1*z0        -&
      1.5*Gp12*(two + z0) + 3*(W0*(five - 4*z0) + 2*V1*(2*z0-five)) - 5*T)/(36*dlt)

    !! masses: (m,m,0)
    else if (masses_ind == 7) then
      z0  = msq(0)/psq
      redcoeff(1) = psq*(-P111 - Gp1*z0/2 + 20*T*z0**3/dlt6-(10*z0**2*(3*T+U*z0))/dlt5 +&
      (2*V1*(one - 2*z0) - Gp12*z0 + 6*P111*z0 + Gp1*(z0 - one)*z0)/(2*dlt) + (U*(one  -&
      4*z0) + 2*V1*(z0-3*one)*z0 + Gp12*(z0 - 0.5*one)*z0+Gp1*z0**2-3*P111*z0**2)/dlt2 -&
      (-1.5*Gp12*z0**2 + Gp1*z0**3/2 - P111*z0**3 - T*(6*z0-one) - 3*z0*(3*V1*z0       +&
      U*(2*z0-3*one)))/dlt3-(z0*(Gp12*z0**2+2*(6*T*(z0-one)+ z0*(-9*U + 2*V1*z0))))/dlt4)

      redcoeff(3) = (Gp12-2*Gp1)/12 + 20*T*z0**2/dlt6 + ((T-P111 +U-V1)*z0**2)/(3*(one +&
      dlt)**2) + (10*z0*(T*(z0 - two) - U*z0))/dlt5+(z0*(-((Gp1-Gp12)*z0)+2*P111*(five +&
      z0) + 2*V1*(five + 2*z0) - 2*U*(five + 3*z0) - 2*T*(five +4*z0)))/(12*(dlt+one)) +&
      (-(Gp12*(z0 -two)*z0)-2*V1*(6*one+5*z0+2*z0**2)+T*(4*one+10*z0+8*z0**2)+U*(-2*z0 +&
      6*(one + z0)**2)-2*P111*(z0*(five+z0)-11*one)+Gp1*(z0*(8*one+ z0)-two))/(12*dlt) +&
      (V1*(7*one - 2*z0)*z0 - Gp12*(z0/2-one)*z0-Gp1*z0**2/2+P111*z0**2+(5*U*(z0*(five +&
      z0)-two))/3 + (T*(3*one + z0*(two + 5*z0)))/3)/dlt3 + (6*P111*(z0-five)*z0       -&
      3*Gp1*(z0-3*one)*z0 + Gp12*(-one + 13*z0 + 2*z0**2) + V1*(46*z0 -34*one+8*z0**2) -&
      2*U*(6*one + 5*z0*(one + z0)) - 2*T*(3*one + z0*(7*one + 6*z0)))/(12*dlt2)       +&
      (U*(13*one-5*z0)*z0-Gp12*z0**2-4*V1*z0**2+(T*(11*one-z0*(41*one + 10*z0)))/3)/dlt4

      redcoeff(4) = -20*T*z0**2/dlt6 + (Gp12*(one - 8*z0) + V1*(34*one -8*z0)-9*Gp1*z0 +&
      30*P111*z0)/(12*dlt2) + (10*z0*(2*T + U*z0))/dlt5 + (Gp1*(one-2*z0)-P111*(11*one +&
      2*z0))/(6*dlt)+(-Gp12*z0-7*V1*z0+Gp1*z0**2/2 -P111*z0**2-U*(8*z0-10*one)/3)/dlt3 +&
      (Gp12*z0**2 + (T*(16*z0-11*one))/3 + z0*(4*V1*z0-13*U))/dlt4

      redcoeff(5) = (-10*T*z0**2)/dlt5 + ((P111 - T - U + V1)*z0**2)/(3*(one +dlt)**2) +&
      (5*z0*(3*U*z0+T*(five+2*z0)))/(3*dlt4)+(Gp1*z0**2-Gp12*z0**2-2*P111*z0*(five+z0) -&
      2*V1*z0*(five + 2*z0) + 2*U*z0*(five + 3*z0) +2*T*z0*(five+4*z0))/(12*(dlt+one)) +&
      (Gp12*(z0-two)*z0 - Gp1*z0*(4*one + z0) +2*P111*z0*(7*one+z0)+2*V1*(6*one + 5*z0 +&
      2*z0**2) - 2*U*(3*one + 5*z0 + 3*z0**2) - 2*T*(2*one + 5*z0 + 4*z0**2))/(12*dlt) +&
      (2*U*(3*one - 2*z0 - z0**2) +2*(T+U)*(3*one+7*z0+6*z0**2)+z0*(3*Gp1*z0-6*P111*z0 -&
      Gp12*(five + 2*z0) - 2*V1*(19*one + 4*z0)))/(12*dlt2)+(-2*T*(3*one+2*z0+5*z0**2) +&
      z0*(3*Gp12*z0 + 12*V1*z0 - 2*U*(17*one + 5*z0)))/(6*dlt3)

      redcoeff(8) = (2*Gp12-7*Gp1)/36 - 10*T*z0/(3*dlt4) + ((T-P111 +U-V1)*z0)/(6*(dlt +&
      one)) - (5*(T*(z0-two) - 2*U*z0))/(6*dlt3) + (12*V1*(z0-five) + 3*Gp12*(z0-two)  +&
      6*Gp1*z0 - 12*P111*z0 - 6*U*(one + 2*z0) - 4*T*(4*one+3*z0))/(72*dlt)+(2*Gp12*z0 +&
      8*V1*z0 + T*(4*z0-one) + U*(5*z0-16*one))/(12*dlt2)

    !! masses: (m0,m1,m1)
    else if(masses_ind == 8) then
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z01 = one + z0 - z1

      redcoeff(1) = psq*(-P111 + (V1*(z01 - 2*z1))/dlt - (Gp1*z1)/2 -(6*Gm*z01*z1)/dlt -&
      (T*z01*((one + z0)**2 - 2*(4*one + z0)*z1 + z1**2))/dlt3 - ((T - W0)*(one +z0**2 +&
      z0*(two - 4*z1) - 6*z1 + 3*z1**2))/dlt2)

      redcoeff(2) = -((T*(one + 3*z0 - 8*z1))/2 + ((2*Gp1*(five + z0-z1) - Gp2*(10*one +&
      z0 - z1))*z1)/4 + 3*W0*(2.5*z1 - z0) -(P111*z1*(7*one-4*z0+4*z1))/2+ (V1*((7*one +&
      z0)*z1 - 12*one - z1**2))/2 + (T*(3*(one + z0**2) + z0*(6*one - 11*z1) - 19*z1   +&
      8*z1**2))/dlt3 + (W0*(z0*(z1 - 3*one) - 3*one - (z1 - 11*one)*z1) + (T*(15*one   +&
      3*z0*(9*one + 4*z0) - 84*z1 - 43*z0*z1 + 31*z1**2))/2)/dlt2 - dlt*(T/2 + 3*V1    -&
      1.5*W0 + 0.5*P111*z1*(7*one + 2*(z1 - z0)) + (z1*(5*Gp2 + Gp1*(z1 -five-z0)))/4) +&
      (6*Gm*(five + z0 - z1)*z1 + W0*(-9*one + 4*z0*(-3*one + z1) + 37*z1 - 4*z1**2)   +&
      T*(11*one + 6*z0*(3*one + z0) - 54*z1 -21*z0*z1+15*z1**2)+V1*(-6*one + z1*(7*one -&
      z0 + z1)))/(2*dlt))/(3*(one + dlt)**2)

      redcoeff(3) = -(3*dlt**2*Gp + dlt*(-5.5*P111 + T/2 + 3*V1 - 1.5*W0 + (3*Gm*(two  +&
      3*z0 - 15*z1))/2 + (3*Gp*(6*one + z0 - 5*z1))/2) + (V1*(29*one + 9*z0 -15*z1))/2 +&
      (P111*(-22*one + 5*z0 - 7*z1))/2 + (3*W0*(3*z0 - 5*z1))/2 + T*(-one/2-2*z0+4*z1) -&
      (11*T*(one + z0**2 + 2*z0*(one - z1) - 38*z1/11 + z1**2))/dlt4 + (-2*Gp2*(z0**2  +&
      2*z0*(one - z1) - (10*one - z1)*z1) + Gp1*(6*one+ z0**2 -2*z0*(z1 - 3*one)-30*z1 +&
      z1**2))/4 + (T*(-35*one - 66*z0-31*z0**2+114*z1+62*z0*z1- 31*z1**2) + W0*(10*one +&
      11*z0 + z0**2 - 19*z1 - 2*z0*z1 + z1**2))/dlt3 + (T*(-77*one - 136*z0 - 57*z0**2 +&
      236*z1 + 114*z0*z1 - 57*z1**2) + 6*Gm*(one + z0**2 + 2*z0*(one - z1) - 10*z1     +&
      z1**2) + 2*W0*(23*one + 29*z0 + 3*z0**2 - 49*z1 - 6*z0*z1 + 3*z1**2)-V1*(-17*one +&
      z0**2 + 8*z1 + z1**2 - 2*z0*(two + z1)))/(2*dlt2)+(T*(-31*one - 52*z0 - 17*z0**2 +&
      92*z1 + 34*z0*z1 - 17*z1**2) + Gp1*(two + 3*z0 + z0**2 - 15*z1 - 2*z0*z1 +z1**2) -&
      P111*(11*one - 5*z0 + 2*z0**2 + 7*z1 - 4*z0*z1 + 2*z1**2) - V1*(-40*one - 13*z0  +&
      3*z0**2 + 23*z1 - 6*z0*z1 + 3*z1**2) + W0*(29*one + 6*z0**2 + z0*(45*one-12*z1)  -&
      75*z1 + 6*z1**2) - Gp2*(one + 2.5*z0 - 12.5*z1 - 3*z0*z1 + 1.5*(z0**2            +&
      z1**2)))/(2*dlt))/(3*(one + dlt)**2)

      redcoeff(4) = -((11*T*(one + 2*z0 + (z0 - z1)**2 - 38*z1/11 ))/dlt4              +&
      P111*(11*one - 5*z0 + 2*z0**2 + 7*z1 - 4*z0*z1 + 2*z1**2) + (V1*(-17*one + z0**2 +&
      8*z1 + z1**2 - 2*z0*(two + z1)))/2 + (Gp2*(one + z0**2 + 2*z0*(one - z1) - 10*z1 +&
      z1**2) + Gp1*(-five + z0**2 + 20*z1 + z1**2 - 2*z0*(two + z1)))/4 + (T*(32*one   +&
      23*z0**2 + z0*(55*one - 46*z1) - 95*z1 + 23*z1**2) + W0*(-10*one - z0**2 + 19*z1 -&
      z1**2 + z0*(-11*one + 2*z1)))/dlt3 + (dlt*(Gp1*(-two +z0**2+5*z1+z1**2 - z0*(one +&
      2*z1)) + 2*P111*(11 + 2*z0**2 + 7*z1 + 2*z1**2 - z0*(5 + 4*z1))))/4 +(-3*Gm*(one +&
      z0**2 + 2*z0*(one - z1) - 10*z1 + z1**2) - 2*W0*(10*one + 11*z0 + z0**2 - 19*z1  -&
      2*z0*z1 + z1**2) + T*(31*one + 13*z0**2 + z0*(44*one -26*z1) - 76*z1 + 13*z1**2) +&
      V1*(-8.5*one + 4*z1 - z0*(two + z1) + (z0**2 + z1**2)/2))/dlt2 +(T*(10*one+z0**2 +&
      z0*(11*one - 2*z1) - 19*z1 + z1**2) + V1*(-17*one - 4*z0 + z0**2 + 8*z1 -2*z0*z1 +&
      z1**2) + (P111*(11*one - 5*z0 + 2*z0**2+7*z1-4*z0*z1 +2*z1**2))/2 + W0*(-10*one  -&
      z0**2 + 19*z1 - z1**2 + z0*(-11*one + 2*z1)) + (-(Gp1*(4*one+5*z0+z0**2 -25*z1   -&
      2*z0*z1 + z1**2)) + 2*Gp2*(one + 2*z0 + z0**2 - 10*z1 - 2*z0*z1                  +&
      z1**2))/4)/dlt)/(3*(one + dlt)**2)

      redcoeff(5) = -(z0*(T/2 - 1.5*W0  + (-2*Gp1*z01 + Gp2*(one + z01))/4 + (T*(one   +&
      19*z01/2) - W0*(3*one + z01))/dlt2 + (5*T*z01)/dlt3 + (V1*(-five - z0 + z1))/2   +&
      (P111*(five + 4*(z1 - z0)))/2 + dlt*((Gp2 - Gp1*z01)/4 + P111*(2.5*one - z0+z1)) -&
      ((-3*T*(4 + 3*(z0 - z1)))/2 + W0*(5.5*one + 2*(z0 - z1)) + 3*Gm*z01+(V1*(five-z0 +&
      z1))/2)/dlt))/(3*(one + dlt)**2)

      redcoeff(8) = -((dlt**2*(5*Gp1+2*Gp2))/12-(5*T*z01*(one+z1))/dlt3 + (2*W0*(3*one +&
      z01)*(one + z1) + T*(-26*one + 3*z1 + 19*z1**2 - z0*(24*one + 19*z1)))/(2*dlt2)  +&
      (36*Gm*z01*(one + z1) - 2*T*(64*one - 18*z1 - 27*z1**2 + 27*z0*(two + z1))       -&
      6*(V1*(-five + z0 - z1)*(one + z1) - W0*(15*one + 5*z0 + 6*z1 + 4*z0*z1          -&
      4*z1**2)))/(12*dlt) + (dlt*(Gp2*(one - 3*z0) + Gp1*(13*one - 3*z1**2 + 3*z0*(one +&
      z1)) + 2*(2*T + 15*V1 - 3*W0 + 3*P111*(z0 + 2*z0*z1 - 2*z1*(3*one + z1)))))/12   +&
      (Gp1*(11*one - 6*z1**2 + 6*z0*(one + z1)) - Gp2*(4*one - 3*z1**2 +3*z0*(two+z1)) +&
      2*(T*(-14*one - 12*z0 + 9*z1) + 3*P111*(z0 + 4*z0*z1 - 2*z1*(3*one + 2*z1))      +&
      3*(3*W0*(two + z0) + V1*(10*one - z0*(one-z1)+6*z1-z1**2))))/12)/(3*(one + dlt)**2)

    !! masses: (m0,m0,m2)
    else if(masses_ind == 9) then
      z0 = msq(0)/psq
      z2 = msq(2)/psq
      z02 = z0 - z2
      z22 = z02**2

      redcoeff(1) = psq*(P111*(-one + z02/dlt)**3 + (U*(dlt3*(one - 4*z0)+3*dlt2*(2*z0 -&
      3*one)*z02 + 18*dlt*z02**2 - 10*z02**3))/dlt5 + (T*(dlt3*(6*z0-one)+12*dlt2*(one -&
      z0)*z02 - 30*dlt*z02**2 + 20*z02**3))/dlt6 + (Gp1*(dlt - z02)*(z0**2+z2*(dlt+z2) -&
      z0*(dlt + dlt2 + 2*z2)))/(2*dlt3) - (Gp12*(-dlt + 2*z02)*(z0**2 + z2*(dlt + z2)  -&
      z0*(dlt + dlt2 + 2*z2)))/(2*dlt4) - (V1*(-dlt + z02)*(4*z0**2 + (dlt + z2)*(dlt  +&
      4*z2) - z0*(dlt*(5 + 2*dlt) + 8*z2)))/dlt4)

      redcoeff(3) = (Gp1*(-2*(one + dlt)*dlt2 +dlt*(9*one+8*dlt)*z0-((6*one+dlt*(9*one +&
      2*dlt))*z22)/(dlt+one)-dlt*(9*one+4*dlt)*z2))/(12*dlt3)+(P111*(11*(one+dlt)*dlt2 -&
      5*dlt*(3*one + 4*dlt)*z0+((6*one+dlt*(15*one+11*dlt))*z22)/(one+dlt)+dlt*(15*one +&
      22*dlt)*z2))/(6*dlt3*(dlt+one))+U*(((13*z02 -5*z22))/dlt4 + z22/(3*(one+dlt)**2) -&
      (10*z22)/dlt5 + (-10*one + 25*z0 + 5*z22 - 17*z2)/(3*dlt3) + (5*z0 + 3*(one+z22) -&
      7*z2)/(6*dlt) + (-5*z02 - 3*z22 + 2*z2)/(6*(dlt + one)) + (-6*one - 5*z0 - 5*z22 +&
      13*z2)/(6*dlt2))-V1*(-(((7*one-2*z02)*z02)/dlt3)+z22/(3*(one+dlt)**2)+4*z22/dlt4 +&
      (6*one + 5*z0 + 2*z22 - 7*z2)/(6*dlt) +(-5*z02-2*z22+2*z2)/(6*(dlt+one))+(17*one -&
      23*z0 -4*z22+19*z2)/(6*dlt2))+T*((10*(-2*one+z02)*z02)/dlt5+z22/(3*(one+dlt)**2) +&
      (20*z22)/dlt6 + (3*one + 2*z0 + 5*z22 - 10*z2)/(3*dlt3) + (2*one + 5*z02 + 4*z22 -&
      2*z2)/(6*dlt) -(3*one+7*z02+6*z22-4*z2)/(6*dlt2)-(5*z0+4*z22-7*z2)/(6*(dlt+one)) +&
      (11*one - 41*z0-10*z22+25*z2)/(3*dlt4))+Gp12*(one/12+ ((two - z02)*z02)/(2*dlt3) -&
      z22/dlt4 + z22/(12*(dlt + one)) + (-one + 13*z0 + 2*z0**2 - 5*z2 - 4*z0*z2       +&
      2*z2**2)/(12*dlt2) + (-z22 + 2*(z0 + z2))/(12*dlt))

      redcoeff(4) = (T*(-(dlt2*(11*one - 16*z0)) + 60*(dlt - z02)*z02))/(3*dlt6)       +&
      (Gp1*(dlt2*(two-4*z0)-9*dlt*z02+6*z02**2))/(12*dlt3)+Gp12*((one-8*z0)/(12*dlt2)  -&
      z02/dlt3 +z02**2/dlt4)-(P111*(dlt2*(11*one+2*z0)+3*z02*(-5*dlt+2*z02)))/(6*dlt3) +&
      (V1*(dlt2*(17*one - 4*z0) + 6*z02*(-7*dlt + 4*z02)))/(6*dlt4) +(U*(2*dlt2*(5*one -&
      4*z0) + 3*z02*(-13*dlt + 10*z02)))/(3*dlt5)

      redcoeff(5) = z0*((Gp12*(-3*(one + dlt)*dlt2*z0 + dlt*(-five - 4*dlt + dlt2)*z02 +&
      (6*one + 4*dlt - dlt2)*z22))/(12*dlt3*(dlt + one)*z02) +(V1*(6*dlt2-(dlt*(19*one +&
      33*dlt + 14*dlt2)*z02)/(one + dlt)**2 + (2*(6*one + 10*dlt + 3*dlt2)*z22)/(one   +&
      dlt)**2))/(6*dlt3*z02) + (Gp1*(-4*dlt*(one + dlt) + (3*one + 2*dlt)*z0 - (3*one  +&
      2*dlt)*z2))/(12*dlt2*(dlt + one)) + (P111*(dlt*(one + dlt)*(7*one +2*dlt)-(3*one +&
      5*dlt)*z0 + (3*one+5*dlt)*z2))/(6*(one+dlt)**2*dlt2)+U*((-17*one-5*z02)/(3*dlt3) -&
      z02/(3*(one + dlt)**2)+5*z02/dlt4+(five +3*z02)/(6*(dlt + one)) + (6*one + 5*z02 +&
      5*z22 - 3*z2)/(6*dlt2*z02) +(-five/6-z0/(2*z22)+(1.5*z0*z2)/z02 - (z0**2 + z0*z2 +&
      z2**2)/(2*z02))/dlt) +  T*( -z02/(3*(one + dlt)**2) - 10*z02/dlt5 + (10*(2.5*one +&
      z02))/(3*dlt4) + (five + 4*z02)/(6*(dlt + one)) - (2*z0**2 + z02**3*(five        +&
      4*z02))/(6*dlt*z02**3) + (-3*one - 2*z02 - 5*z22 + 3*z2)/(3*dlt3*z02)+((7*one)/6 +&
      z0/(2*z22) - (3*z0*z2)/z02 + (z0**2 + z0*z2 + z2**2)/z02)/dlt2))

      redcoeff(6) = z2*((Gp1*(4*dlt*(one + dlt) - (3*one + 2*dlt)*z02))/(12*dlt2*(dlt  +&
      one)) + (P111*(-7*dlt*(one + dlt) + (3*one + 5*dlt)*z02))/(6*(one +dlt)**2*dlt2) +&
      (Gp12*(3*(one + dlt)*dlt2*z0 + dlt*(5*one + (3*one - 2*dlt)*dlt)*z02 + (-6*one   +&
      (-4*one + dlt)*dlt)*z22))/(12*dlt3*(dlt + one)*z02)+ (V1*(-6*dlt2 + (dlt*(19*one +&
      12*dlt)*z02)/(dlt + one) - (2*(6*one + dlt*(10*one + 3*dlt))*z22)/(one           +&
      dlt)**2))/(6*dlt3*z02) + (U*(3*(one + dlt)**2*dlt3*z0 +3*(one+dlt)**2*dlt2*(-two +&
      dlt + z0)*z02 + (-30*one + dlt*(-50*one + 3*(dlt-five)*dlt))*z02**3 +dlt*(34*one +&
      dlt*(55*one + 3*(five - 2*dlt)*dlt))*z22))/(6*(one +dlt)**2*dlt4*z22)+T*((-7*one -&
      4*z02)/(6*(dlt + one)) + z02/(3*(one + dlt)**2) + 10*z02/dlt5 - (5*(five         +&
      2*z02))/(3*dlt4) + ((3*(one - z0))/z02 + 5*(two + z02))/(3*dlt3) +(-11*one-6*z02 -&
      (3*(z0 + z02))/z22)/(6*dlt2) + (2*(z0 + z02)**2 + z02*(-2*z0 + 4*z02**3          +&
      7*z22))/(6*dlt*z02**3)))

      redcoeff(8) = -(P111*(2*dlt*z0 + z02))/(6*dlt*(dlt + one)) + (Gp1*(-7*dlt        +&
      3*z02))/(36*dlt) + (Gp12*(2*dlt*(-3*one + 2*dlt) + 3*(4*one + dlt)*z0 +3*(-4*one +&
      dlt)*z2))/(72*dlt2) + (V1*(-5*dlt*(one + dlt) + (4*one + 5*dlt)*z0 - (4*one      +&
      3*dlt)*z2))/(6*dlt2*(dlt + one)) + (U*((20*z02)/dlt3 + (2*(z0 + z2))/(dlt + one) +&
      (-16*one + 5*(z0 + z2))/dlt2 - (z0 + 5*z2 + 2*z02*(z0 + z2))/(dlt*z02)))/12      +&
      (T*((-120*z02)/dlt4 - (30*(-two + z0 + z2))/dlt3 + (6*(z0 + z2))/(dlt + one)     +&
      (3*(z0*(-one + 4*z0) +(7*one-4*z2)*z2))/(dlt2*z02)+(-6*z0**3+2*(five-3*z2)*z2**2 +&
      2*z0*z2*(-7*one + 3*z2) + z0**2*(-8*one + 6*z2))/(dlt*z22)))/36

    !! masses: (m0,m1,m0)
    !! TODO: think about z10 at the denominators.
    else if(masses_ind == 10) then
      z0 = msq(0)/psq
      z1 = msq(1)/psq
      z01 = one + z0 - z1
      z10 = z0 - z1

      redcoeff(1) = psq*(-P111 - 0.5*Gp1*z1 - (20*T*z10**3)/dlt6 - (10*z10**2*(-W0*z10 +&
      T*(3*z01 + z10)))/dlt5 + (-Gp12*z01*z1 + Gp1*(z0 + z0**2 -3*z0*z1+z1*(2*z1-one)) +&
      2*(V1*(one + z0 - 3*z1) - 3*P111*z10))/(2*dlt) + (Gp12*(z0**3 + z0**2*(two-3*z1) -&
      z1*(one - 4*z1 + z1**2) + z0*(one - 6*z1 + 3*z1**2)) + Gp1*z10**2*(two + z10)    -&
      2*((T - W0)*(one + z0**2 + 2*z0*(one - 2*z1) - 6*z1 + 3*z1**2)+z10*(-2*V1*(3*one +&
      2*z0 - 3*z1) + 3*P111*z10)))/(2*dlt2) + (z10*(Gp12*z10**2 - 2*(6*T*(one + z0**2  -&
      2*z0*(z1 - one) - 3*z1 + z1**2) + z10*((T - W0)*(9*one +6*z10)-2*V1*z10))))/dlt4 +&
      (3*Gp12*z01*z10**2 + Gp1*z10**3 - 2*(T*(one + z0**3 - 3*z0**2*(z1 - one) - 9*z1  +&
      9*z1**2 - z1**3 + 3*z0*(one - 4*z1 + z1**2)) + z10*(3*(T - W0)*(3*one + z0**2    -&
      2*z0*(z1 - two) - 6*z1 + z1**2) + z10*(P111*z10 - 3*V1*(3*one + z10)))))/(2*dlt3))

      redcoeff(3) = (-2*Gp1 + Gp12)/12 + ((W0 - P111 - V1)*z0)/(3*(dlt + one)) + (-2*T +&
      2*W0*(3*one - 2*z0) + 4*V1*(z0 - 3*one) + 2*P111*(11*one + 2*z0) +Gp12*(z0+3*z1) +&
      Gp1*(-two - 5*z0 + 9*z1))/(12*dlt) +(20*T*z10**2)/dlt6+(10*z10*(-W0*z10+2*T*(z01 +&
      z10)))/dlt5 + (-(Gp12*(one + 7*z0 + 6*z0**2 - 15*z1 - 12*z0*z1 + 6*z1**2))       -&
      3*Gp1*(2*z0**2 + z0*(3*one - 4*z1) + z1*(-3*one+2*z1))+V1*(-34*one+8*z0- 54*z10) +&
      30*P111*z10 - 2*W0*(6*one - 8*z1 + z10) + 4*T*(1.5*one - 2*z1 + z10))/(12*dlt2)  +&
      (T*(11*one + 72*z0**2 + z0*(86*one - 144*z1) - 102*z1 + 72*z1**2)- 3*Gp12*z10**2 -&
      3*z10*(W0*(one + 12*z01) + 4*V1*z10))/(3*dlt4)+(-3*Gp12*(3*z0**2+2*z0*(one-3*z1) +&
      z1*(-two + 3*z1)) - 3*Gp1*z10**2 + 2*(T*(13*one + 12*z0**2 + z0*(29*one - 24*z1) -&
      45*z1 + 12*z1**2) + W0*(-10*one - 9*z0**2 + 36*z1-9*z1**2+2*z0*(-14*one + 9*z1)) +&
      3*P111*z10**2 - 3*V1*z10*(7*one + 3*z10)))/(6*dlt3)

      redcoeff(4) = -(Gp1*(z0**2 - two + 5*z1 + z1**2 - z0*(one +2*z1))+2*P111*(11*one +&
      2*z0**2 + 7*z1 + 2*z1**2 - z0*(five + 4*z1)))/(12*dlt) +(3*Gp1*(z0**2 +z0*(3*one -&
      2*z1) + (z1 - 3*one)*z1) + Gp12*(one   + z0**2 - 2*z0*(z1 - one)- 10*z1 + z1**2) +&
      2*(-(V1*(z0**2 - 17*one + 8*z1 + z1**2 - 2*z0*(two +z1)))+3*P111*(z0**2+z1*(five +&
      z1) - z0*(five + 2*z1))))/(12*dlt2) - (20*T*z10**2)/dlt6 - (10*z10*(-W0*z10      +&
      T*(two + 3*z10)))/dlt5 - (-((Gp12 + 4*V1)*z10**2) - W0*z10*(13*one + 7*z10)      +&
      (T*(11*one + 11*z0**2 - 22*z0*(z1 - one) - 38*z1 + 11*z1**2 + 3*z10*(13*one      +&
      7*z10)))/3)/dlt4 + (3*Gp1*z10**2 + 2*(3*Gp12*(z0 + z0**2 - 2*z0*z1 +(z1-one)*z1) -&
      (T - W0)*(10*one + z0**2 + z0*(11*one - 2*z1) - 19*z1 + z1**2) - 3*z10*(P111*z10 -&
      V1*(7*one + z10))))/(6*dlt3)

      redcoeff(5) = z0*((P111 + V1 - W0)/(3*(dlt + one))-(10*T*z10)/dlt5-(5*(-3*W0*z10 +&
      T*(five + 8*z10)))/(3*dlt4) - (-5*Gp12*z01*z10**2 -3*Gp1*z10**3+2*(3*P111*z10**3 +&
      T*(3*z0 + 2*z10**2) - z10*(V1*z10*(19*one + z10) + W0*(6*one - 3*z1              +&
      z10))))/(12*dlt2*z10**2) - (-3*Gp12*z10**2 + 2*(T*(3*(one -z1)+18*z10+11*z10**2) -&
      z10*(6*V1*z10 + W0*(17*one + 8*z10))))/(6*dlt3*z10) + (-(Gp12*(z0 +2*z1)*z10**2) +&
      Gp1*(3*one + z01)*z10**3 + 2*(P111*z10**3*(-7*one + 2*z10) + z10*(W0*(2*z0**2    +&
      2*z1*(3*one + z1) - z0*(3*one + 4*z1)) + 2*V1*(3*one - z10)*z10) + T*(2*z0**3    -&
      2*z0*z1*(3*one + z1) - 2*z0**2*(one+3*z1+ z10) + z0*(6*z1**2 + 3*z10 + z1*(6*one +&
      4*z10)))))/(12*dlt*z10**3))

      redcoeff(6) = z1*((10*T*z10)/dlt5 + (5*(-3*W0*z10 + T*(5*z01 + 3*z10)))/(3*dlt4) +&
      (-6*W0*(2*z01 - z1)*z10 - 2*V1*z10**2*(19*one+z10)+z10**2*(-5*Gp12*z01-3*Gp1*z10 +&
      6*P111*z10) + T*(6*z01*z1 + 6*(2*z01 - z1)*z10))/(12*dlt2*z10**2)-(3*Gp12*z10**2 +&
      12*V1*z10**2 + 2*W0*z10*(17*one + 8*z10) + T*(-6*one - 46*z0 + 52*z1             -&
      22*z10**2))/(6*dlt3*z10) - (-3*Gp12*z1*z10**2 + Gp1*z10**3*(4*one + z10)         +&
      2*(3*W0*z1*z10 + 6*V1*z10**2 + P111*z10**3*(-7*one + 2*z10) - T*z1*(2*z1         +&
      3*z10)))/(12*dlt*z10**3))

      redcoeff(8) = (-7*Gp1 + 2*Gp12)/36. - ((P111 + V1 - W0)*z0)/(3*(dlt + one))      +&
      (10*T*z10)/(3*dlt4) + (5*(-2*W0*z10 +T*(2*(one-z1)+3*z10)))/(6*dlt3)+(T*(4*z0**2 +&
      z0*(9*one - 16*z1) + 3*z1*(-five + 4*z1)) - 2*Gp12*z10**2 + z10*(W0*(-16*one +z0 +&
      9*z1) - 8*V1*z10))/(12*dlt2*z10) + (-6*Gp1*z10**3 - 3*Gp12*z10**2*(2*(one - z1)  +&
      z10) - 2*(-6*P111*z10**3 - 3*z10*(W0*(5*z0 - 4*z0**2 + z1 + 4*z0*z1)+2*V1*(-five +&
      2*z0)*z10) + T*(12*z0*z1 + 5*z10**2)))/(72*dlt*z10**2)

    !! masses: (0,m1,m2)
    else if (masses_ind == 11) then
      z1 = msq(1)/psq
      z2 = msq(2)/psq
      z12 = z1 - z2
      z22 = z12**2

      redcoeff(1) = (P111*psq*(z12-dlt)**3)/dlt3 - (Gp1*psq*(dlt - z12)**2*(dlt*z1     +&
      z12))/(2*dlt3)-(Gp12*psq*(z12-dlt)*(dlt*z1+z12)*((dlt*(z1-one))+2*z12))/(2*dlt4) -&
      (psq*V1*(dlt-z12)**2*(3*dlt*z1-dlt+4*z12))/dlt4+(psq*U*(dlt3*(one+3*(z1-two)*z1) -&
      3*dlt2*(3*one +(z1-6*one)*z1)*z12 - 6*dlt*(2*z1-3*one)*z12**2 - 10*z12**3))/dlt5 +&
      (psq*T*(dlt3*(z1-one)*(one + (z1 -8*one)*z1) + 12*dlt2*(one + (z1-3*one)*z1)*z12 +&
      30*dlt*(z1-one)*z12**2 + 20*z12**3))/dlt6

      redcoeff(3) = -(V1*((4*z12**2)/dlt4 + (z12*(3*z1 -7*one- 2*z2))/dlt3 + ((7*one   +&
      3*z1 - 2*z2)*z2)/(6*(dlt + one)) + z2**2/(3*(one + dlt)**2) + (17*one - 27*z1    +&
      (19*one + 3*z1 - 4*z2)*z2)/(6*dlt2) + (one + (z2*(-7*one - 3*z1 +2*z2))/6)/dlt)) +&
      U*(-((z12*(-13*one + 12*z1 - 5*z2))/dlt4) + ((7*one + 3*z12)*z2)/(6*(dlt+one))   +&
      z2**2/(3*(one + dlt)**2) + (-6*one + 9*z1 + 13*z2 + 3*z1*z2 - 5*z2**2)/(6*dlt2)  +&
      (-10*one + 36*z1 - 9*z1**2 - 17*z2+3*z1*z2+5*z2**2)/(3*dlt3)+(3*one + z2*(-7*one -&
      3*z12))/(6*dlt) - 10*z22/dlt5) + Gp1*(-one/6 + ((3*one - z1 - z12)*z12)/(4*dlt2) -&
      z2**2/(12*(dlt + one)) - z22/(2*dlt3)  + (-two + 5*z1 + z1**2 + 4*z12 - 2*z1*z12 +&
      z22)/(12*dlt)) + P111*((z12*(-five - z1 + z12))/(2*dlt2) + (11*one + 2*z1**2     -&
      z1*(-7*one + z12) - z12*(7*one + z12))/(6*dlt) - z2**2/(3*(one+dlt)**2)+z22/dlt3 +&
      (-2*z1**2 + z1*z12 - 7*z2 + z22)/(6*(dlt + one))) + Gp12*(one/12 +((2*(one - z1) -&
      z12)*z12)/(2*dlt3) - z12**2/dlt4 +(5*z1-z1**2-z12*(-2*(-one+z1) + z12))/(12*dlt) +&
      z2**2/(12*(dlt + one)) + (-one +(10*one-z1)*z1+(five-7*z1)*z12+2*z22)/(12*dlt2)) +&
      T*((10*z12*(-two + 2*z1 + z12))/dlt5+(-3*one + 14*z1 - 3*z1**2 + (-11*one + 9*z1 -&
      6*z12)*z12)/(6*dlt2) + (11*one-38*z1+11*z1**2+5*(-five+7*z1-2*z12)*z12)/(3*dlt4) +&
      (3*one - 19*z1 + 8*z1**2 + 5*z12*(2*(one - z1) + z12))/(3*dlt3) + z2**2/(3*(one  +&
      dlt)**2) + (-z1**2 + 5*z1*z12 + 7*z2 - 4*z22)/(6*(dlt + one)) + (20*z22)/dlt6    +&
      (two - 7*z1 + z1**2 + 7*z12 - 5*z1*z12 + 4*z22)/(6*dlt))

      redcoeff(4) = (Gp1*(-(dlt2*(z1*(five + z1)-two)) + 3*dlt*(z1-3*one)*z12          +&
      6*z22))/(12*dlt3) - (P111*(dlt2*(11*one + z1*(7*one + 2*z1))-3*dlt*(five+z1)*z12 +&
      6*z22))/(6*dlt3) + (Gp12*(dlt2*(one + (z1-10*one)*z1) + 12*dlt*(z1-one)*z12      +&
      12*z22))/(12*dlt4) + (V1*(-(dlt2*(-17*one + z1*(8*one+z1)))+6*dlt*(z1-7*one)*z12 +&
      24*z22))/(6*dlt4) + (U*(dlt2*(10*one + (z1-19*one)*z1) + 3*dlt*(7*z1-13*one)*z12 +&
      30*z22))/(3*dlt5) + (T*(dlt2*(-11*one + (38 - 11*z1)*z1) - 60*(dlt*(z1-one)*z12  +&
      z22)))/(3*dlt6)

      redcoeff(5) = z1*((V1*(6*dlt2 + dlt*(-19*one + z1)*z12 + 12*z22))/(6*dlt3*z12)   +&
      (Gp1*(-4*dlt + (3*one + dlt)*z1 -3*z2))/(12*dlt2)+(P111*(7*dlt+(-3*one+2*dlt)*z1 +&
      3*z2))/(6*dlt2) +(Gp12*(-3*dlt2*z1+5*dlt*(-one + z1)*z12 + 6*z22))/(12*dlt3*z12) +&
      (U*(-3*dlt3*z1 - 3*dlt2*(-2*one + 3*z1)*z12 + 30*z12**3 + 2*dlt*(-17*one         +&
      8*z1)*z22))/(6*dlt4*z22) - (T*(2*dlt4*z1**2 + 3*dlt3*(z1-one)*z1*z12+6*dlt2*(one +&
      (-3*one + z1)*z1)*z22 + 50*dlt*(-one + z1)*z12**3 + 60*z22**2))/(6*dlt5*z12**3))

      redcoeff(6) = z2*(-(V1*((-19*one+z1-4*z12)/(6*dlt2)-(7*one + z1 + 2*z12)/(6*(dlt +&
      one))+2*z12/dlt3+(6*one+z12*(7*one+z1+2*z12))/(6*dlt*z12)-z2/(3*(one+ dlt)**2))) +&
      T*((-7*one + z1-4*z12)/(6*(dlt+one))+10*z12/dlt5-(5*(5*(one-z1)+2*z12))/(3*dlt4) +&
      (3*one - 9*z1 + 3*z1**2 + 5*z12*(2*one-z2))/(3*dlt3*z12)+(3*z1**2+3*z1*(-one+z12 +&
      z22) - z12*(3*one + 11*z12 + 6*z22))/(6*dlt2*z22) + (2*z1**2-z1*z12*(-2*one+z22) +&
      z22*(two + z12*(7*one +4*z12)))/(6*dlt*z12*z22)-z2/(3*(one+dlt)**2))-U*((-17*one +&
      8*z1 - 5*z12)/(3*dlt3) + 5*z12/dlt4 + (7*one + 3*z12)/(6*(dlt +one))+(6*one-9*z1 +&
      z12*(13*one-2*z1+5*z12))/(6*dlt2*z12)-(3*z1+z12*(3*one+7*z12+3*z22))/(6*dlt*z22) +&
      z2/(3*(one + dlt)**2)) + (Gp1*(4*dlt*(one + dlt) - 3*(one + dlt)*z1 + (3*one     +&
      2*dlt)*z2))/(12*dlt2*(dlt + one)) - (P111*(7*dlt*(one+dlt)-3*(one+dlt)*z1+(3*one +&
      5*dlt)*z2))/(6*(one + dlt)**2*dlt2) + Gp12*(-z12/(2*dlt3) + (5*(one - z1)        +&
      2*z12)/(12*dlt2) +(3*z1 + z12*(-two + z2))/(12*dlt*z12) - z2/(12*(dlt + one))))

      redcoeff(8) = -(P111*(dlt*z1+z12))/(6*dlt*(dlt+one))+ Gp1*(3*z12-7*dlt)/(36*dlt) +&
      Gp12*(one/18 + (-two +4*z1-z12)/(24*dlt)+z12/(6*dlt2))-V1*(-2*z12/(3*dlt2)+(five -&
      z2)/(6*dlt) +z2/(6*(dlt+one)))+U*((-16*one+14*z1-5*z12)/(12*dlt2)+5*z12/(3*dlt3) +&
      z2/(6*(dlt +one))+(-6*z1+5*z12-2*z1*z12+2*z22)/(12*dlt*z12))+(T*((-30*(-two+4*z1 -&
      z12))/dlt3 - 120*z12/dlt4 + 6*z2/(dlt + one) - (3*(6*(z1-one)*z1 +7*(one-z1)*z12 +&
      4*z22))/(dlt2*z12) - (2*(6*z1**2-(five+3*z12)*z22 + 3*z1*(z12+z22)))/(dlt*z22)))/36

    !! masses: (m0,0,m2)
    else if(masses_ind == 12) then
      z0 = msq(0)/psq
      z2 = msq(2)/psq
      z02 = one + z0

      redcoeff(1) = psq*(-P111 - (20*T*z2**3)/dlt6 - (10*z2**2*(3*T*z02 - U*z2))/dlt5  +&
      (Gp1*z02*z2 - 2*(-V1*z02+3*P111*z2))/(2*dlt)+(Gp12*z02**2*z2+Gp1*(one+z02)*z2**2 -&
      2*(-U*z02**2 + z2*(-2*V1*(one + 2*z02) + 3*P111*z2)))/(2*dlt2) + (z2*(Gp12*z2**2 -&
      2*(6*T*z02**2 - z2*(3*U*(one + 2*z02) + 2*V1*z2))))/dlt4 + (3*Gp12*z02*z2**2     +&
      Gp1*z2**3 - 2*(T*z02**3 + z2*(-3*U*(3*one + 4*z0 + z0**2) + z2*(-3*V1*(two +z02) +&
      P111*z2))))/(2*dlt3))

      redcoeff(3) = (-2*Gp1 + Gp12)/12 + ((T - P111+U-V1)*(z0-z2)**2)/(3*(one+dlt)**2) +&
      (20*T*z2**2)/dlt6 + (10*z2*(-U*z2 + T*(2*z02 +z2)))/dlt5+(-Gp12*z2**2-4*V1*z2**2 -&
      U*z2*(13*one + 7*z0 + 5*z2) + (T*(11*z02**2 + 5*(five + 7*z0 -2*z2)*z2))/3)/dlt4 +&
      (Gp1*(-two - z0 + (z0 - z2)**2 -4*z2)-Gp12*(z0+(z0-z2)**2 - 2*z2) + 2*V1*(-6*one +&
      z0**2 + z0*(-five + z2) + 7*z2 - 2*z2**2) + 2*U*(3*one + z0*(five - 3*z2) - 7*z2 +&
      3*z2**2) + 2*T*(two + z0**2 -5*z0*(-one + z2) - 7*z2 + 4*z2**2) + 2*P111*(11*one +&
      2*z0**2 + 7*z2 - z2**2 - z0*(five+z2)))/(12*dlt) + (-((Gp1 - Gp12)*(z0 - z2)**2) -&
      2*V1*(z0**2+z0*(z2-five)+(7*one -2*z2)*z2) + 2*P111*(-2*z0**2 + (-7*one + z2)*z2 +&
      z0*(five+z2))+2*U*((7*one-3*z2)*z2+z0*(-five+3*z2))-2*T*(z0**2- 5*z0*(-one + z2) +&
      z2*(-7*one + 4*z2)))/(12*(dlt + one)) + (-2*T*(one+3*(z0-z2))*(two + z02 - 2*z2) +&
      6*P111*z2*(five - z0 + z2) -3*Gp1*z2*(3*one+z0+z2) - Gp12*(z02**2 + (five + 7*z0 -&
      2*z2)*z2)-2*U*(6*one+2*z0**2-7*z0*(-two + z2) - 13*z2 + 5*z2**2) + 2*V1*(-17*one +&
      z0**2-19*z2+4*z2**2-z0*(4*one+5*z2)))/(12*dlt2) + (6*P111*z2**2 - 6*V1*z2*(7*one +&
      z0 + 2*z2) + 2*T*(3*one+8*z0**2+z0*(11*one-10*z2)-10*z2+ 5*z2**2) - 3*z2*(Gp1*z2 +&
      Gp12*(2*z02 + z2)) - 2*U*(10*one + z0**2+17*z2-5*z2**2+z0*(11*one+13*z2)))/(6*dlt3)

      redcoeff(4) = (2*P111*(-11*one + (five-2*z0)*z0)+Gp1*(two+(one-z0)*z0))/(12*dlt) -&
      (20*T*z2**2)/dlt6 + (2*V1*(17*one+(4*one-z0)*z0)+Gp12*z02**2+6*P111*(z0-five)*z2 +&
      3*Gp1*(3*one+z0)*z2)/(12*dlt2)-(10*(2*T*z02*z2-U*z2**2))/dlt5 +((-11*T*z02**2)/3 +&
      U*(13*one + 7*z0)*z2 + Gp12*z2**2 + 4*V1*z2**2)/dlt4 + (2*U*(10*one + z0)*z02    +&
      6*V1*(7*one + z0)*z2 - 6*P111*z2**2 + 3*z2*(2*Gp12*z02 + Gp1*z2))/(6*dlt3)

      redcoeff(5) =z0*(((P111 - T - U+V1)*(z0-z2))/(3*(one+dlt)**2)-(10*T*z2)/(3*dlt4) -&
      (5*(T*(z02 - z2) - U*z2))/(3*dlt3) + (T*(7*one + 3*z0 -6*z2)+U*(8*one+2*z0-5*z2) +&
      Gp12*z2 + 4*V1*z2)/(6*dlt2) + (U*(10*one - 6*z2) + 2*T*(five + z0 - 4*z2)        +&
      (Gp1 - Gp12)*(z0-z2)+2*P111*(-five+2*z0+z2)+2*V1*(-five+z0+2*z2))/(12*(dlt+one)) +&
      (-2*T*(five + z0 - 4*z2) + Gp12*(z02 - z2)+Gp1*z2-2*P111*z2-2*V1*(-five+z0+2*z2) +&
      U*(-10*one + 6*z2))/(12*dlt))

      redcoeff(6) = ((T -P111 + U - V1)*(z0-z2)*z2)/(3*(one+dlt)**2)-(10*T*z2**2)/dlt5 +&
      (5*z2*(3*U*z2 +T*(-5*z02+2*z2)))/(3*dlt4)+(z2*((Gp1-Gp12)*(z2-z0)-2*P111*(-7*one +&
      2*z0 + z2) - 2*V1*(-7*one + z0 + 2*z2)  + 2*U*(-7*one + 3*z2) + 2*T*(-7*one - z0 +&
      4*z2)))/(12*(dlt+one))+(2*U*(17*one+8*z0 - 5*z2)*z2 + 3*Gp12*z2**2 + 12*V1*z2**2 -&
      2*T*(3*z02**2-5*(one+z02)*z2+5*z2**2))/(6*dlt3) +(2*P111*z2*(-7*one + 2*z0 + z2) +&
      2*U*(-3*(one - z2)**2 + z2) + 2*T*(-two + (7*one + z0)*z2 - 4*z2**2)+2*V1*(6*one +&
      (-7*one + z0)*z2 + 2*z2**2) +z2*(Gp1*(4*one+z0-z2) + Gp12*(z2-two-z0)))/(12*dlt) +&
      (2*V1*(19*one + z0 - 4*z2)*z2 - 6*P111*z2**2 + z2*(Gp12*(5*z02 - 2*z2)+3*Gp1*z2) +&
      T*(6*(one + z0*(one - z2)) - 22*z2 + 12*z2**2) + 2*U*(-2*z0*(-3*one + z2)        +&
      (-two + z2)*(-3*one + 5*z2)))/(12*dlt2)

      redcoeff(8) = (2*Gp12 - 7*Gp1)/36 + (10*T*z2)/(3*dlt4) + ((T - P111 + U -V1)*(z0 +&
      z2))/(6*(dlt+one)) - (5*(2*U*z2 + T*(-2*z02+z2)))/(6*dlt3)+(-6*Gp1*z2+12*P111*z2 +&
      12*V1*(-five + z0 + z2) + 3*Gp12*(-2*z02 + z2) - 6*U*(-five + 2*z0 + 2*z2)       -&
      4*T*(-5*one + 3*z0 + 3*z2))/(72*dlt) + (-2*Gp12*z2 - 8*V1*z2 +T*(-7*one+z0+4*z2) +&
      U*(-16*one - 4*z0 + 5*z2))/(12*dlt2)

    !! masses: (m0,m1,0)
    else if(masses_ind == 13) then
      z0 = msq(0)/psq
      z1 = msq(1)/psq
      z01 = one + z0 - z1
      z10 = z0 - z1
      redcoeff(1) = psq*(-P111 - Gp1*z1/2 + (20*T*z1**3)/dlt6 - (10*z1**2*(3*T*z01     +&
      U*z1))/dlt5+(2*V1*(z01-2*z1)+6*P111*z1 - Gp12*z01*z1 + Gp1*z1*(z1 -z01))/(2*dlt) +&
      ((Gp1*(one + z01)*z1**2)/2 - z1*(2*V1*(3*one + 2*z0 - 3*z1) + 3*P111*z1)         -&
      U*(-one - z0**2 - 2*z0*(one - 2*z1) - 3*(z1 - two)*z1) - (Gp12*z1*((one + z0)**2 -&
      2*(two + z0)*z1 + z1**2))/2)/dlt2 + (2*z1*(-(Gp12*z1**2)/2 + z1*(3*U*(one+2*z01) -&
      2*V1*z1) + 6*T*((one + z0)**2 - (3*one+2*z0)*z1 +z1**2)))/dlt4+(3*Gp12*z01*z1**2 -&
      Gp1*z1**3-2*(-(z1**2*(3*V1*(two+z01)+P111*z1))+3*U*z1*(3*one+z0**2-2*z0*(z1-two) -&
      6*z1 + z1**2) + T*z01*((one + z0)**2 - 2*(4*one + z0)*z1 + z1**2)))/(2*dlt3))

      redcoeff(3) = (-2*Gp1 + Gp12)/12 + ((-P111 + T + U - V1)*z0**2)/(3*(one+dlt)**2) +&
      (20*T*z1**2)/dlt6 + (10*z1*(-U*z1 + T*(-2*z01+z1)))/dlt5+(z0*(((-Gp1+Gp12)*z0)/2 -&
      V1*(z0 -five - 3*z1) - U*(five + 3*z1) + P111*(five - 2*z0 + 3*z1) - T*(five +z0 +&
      3*z1)))/(6*(dlt + one)) + (3*U*(13*one + 7*z0-12*z1)*z1-3*Gp12*z1**2-12*V1*z1**2 +&
      T*(11*(one + z0)**2 -3*(21*one+19*z0)*z1+36*z1**2))/(3*dlt4)+(-(Gp12*(z0 + z0**2 -&
      3*z1)) + Gp1*(-two-z0+z0**2+9*z1) + P111*(22*one + 4*z0**2 - 2*z0*(five + 3*z1)) +&
      2*V1*(-6*one+z0**2-z0*(five+3*z1))+2*T*(two + z0**2 + z0*(five+3*z1)) + U*(6*one +&
      2*z0*(five + 3*z1)))/(12*dlt) + (6*P111*(z0 - five)*z1 +3*Gp1*(3*one+z0-2*z1)*z1 +&
      2*V1*(-17*one - 4*z0 + z0**2 +3*(9*one+z0)*z1)-Gp12*(one+z0**2 + z0*(two - 9*z1) -&
      15*z1+6*z1**2)-2*T*(3*(one+z0**2-z1) +z0*(10*one + 3*z1)) - 2*U*(6*one + 2*z0**2 -&
      9*z1 + z0*(14*one + 3*z1)))/(12*dlt2) + (6*V1*(7*one +z0-3*z1)*z1 + 6*P111*z1**2 -&
      2*U*(10*one + z0**2 + z0*(11*one - 15*z1) - 9*(4*one - z1)*z1) - 3*z1*(Gp1*z1    +&
      Gp12*(-2*z01 + z1)) + 2*T*(3*one + 11*z0 + 5*z0**2 - 9*z1 + 3*z10**2))/(6*dlt3)

      redcoeff(4) = (-20*T*z1**2)/dlt6 + (10*z1*(2*T*z01 + U*z1))/dlt5 -(3*Gp12*z01*z1 -&
      1.5*Gp1*z1**2 + 3*z1*(V1*(6*one + z01) + P111*z1) - U*(10*one + 11*z0 - 19*z1    +&
      z10**2))/(3*dlt3) + (3*P111*(6*one - z01)*z1 -(3*Gp1*(two+z01)*z1)/2+(Gp12*((one +&
      z0)**2-2*(five+z0)*z1 + z1**2))/2 - V1*(-17*one -4*z0 + 8*z1 + z10**2))/(6*dlt2) -&
      (Gp1*(-two-z0+5*z1 + z10**2) +2*P111*(11*one - 5*z0 + 7*z1 + 2*z10**2))/(12*dlt) +&
      (Gp12*z1**2 - z1*(U*(6*one + 7*z01) - 4*V1*z1) - (T*(-38*z1 + 11*(one + 2*z0     +&
      z10**2)))/3)/dlt4

      redcoeff(5) = z0*(((P111 - T - U + V1)*z0)/(3*(one + dlt)**2) + 10*T*z1/(3*dlt4) -&
      (5*(T*(one + z0)+U*z1))/(3*dlt3)+(Gp12*(one+z0)-Gp1*z1+2*P111*z1-2*U*(five+3*z1) +&
      2*V1*(five - z0 + 3*z1) - 2*T*(five + z0 + 3*z1))/(12*dlt) + ((Gp1 - Gp12)*z0    +&
      P111*(-10*one + 4*z0 - 6*z1) + 2*V1*(-five + z0 - 3*z1) + 2*T*(five + z0 + 3*z1) +&
      U*(10*one + 6*z1))/(12*(dlt + one)) + (-Gp12*z1 - 4*V1*z1 + U*(8*one +2*z0+3*z1) +&
      T*(7*one + 3*(z0 + z1)))/(6*dlt2))

      redcoeff(6) = (-10*T*z1**2)/dlt5 + (5*z1*((5*T*z01)/3 +U*z1))/dlt4+(3*Gp12*z1**2 +&
      12*V1*z1**2 - 6*T*((one + z0)**2 - (3*one + 2*z0)*z1 + z1**2) + 2*U*z1*(-17*one  -&
      8*z10))/(6*dlt3) + (-4*T - 6*U + 12*V1 + z1*(-3*Gp12 - Gp1*(4*one + z10))        +&
      2*P111*z1*(7*one - 2*z10))/(12*dlt) + (z1*(-5*Gp12*z01 + 3*Gp1*z1) + 6*(T*z01    +&
      U*(2*z01 - z1) - P111*z1**2) - 2*V1*z1*(19 + z10))/(12*dlt2)

      redcoeff(8) = (2*Gp12 -7*Gp1)/36 + ((T - P111 + U - V1)*z0)/(6*(dlt + one))      -&
      (10*T*z1)/(3*dlt4) + (12*V1*(z0 - five) - 6*U*(one + 2*z0) - 4*T*(4*one + 3*z0)  -&
      3*Gp12*(2*z01 - z1) + 6*Gp1*z1 - 12*P111*z1)/(72*dlt) - (5*(-2*U*z1 + T*(-2*z01  +&
      z1)))/(6*dlt3) + (T*(-one + 7*z0 - 3*z1) + 2*Gp12*z1 + 8*V1*z1 + U*(-16*one-4*z0 +&
      9*z1))/(12*dlt2)

    !! masses: (m0,m1,m2)
    else if (masses_ind == 14) then
      z0  = msq(0)/psq
      z1  = msq(1)/psq
      z2  = msq(2)/psq
      z12 = z1 - z2
      z01 = one + z0 - z1
      z22 = z12**2

      redcoeff(1) = psq*(-P111-Gp1*z1/2+ 20*T*z12**3/dlt6+(V1*(z01-2*z1)-Gp12*z01*z1/2 +&
      3*P111*z12 - (Gp1*(z01 -z1)*z12)/2)/dlt+(U*(z01**2+2*z1*(z1-two-z0))-2*V1*(3*one +&
      2*z0-3*z1)*z12-(Gp12*(z01**2-2*z1)*z12)/2-3*P111*z22+(Gp1*(one+z01)*z22)/2)/dlt2 +&
      (-(T*z01*(z01**2-6*z1))-3*U*(z01**2 +2*(z01-z1))*z12+1.5*Gp12*z01*z22 +3*V1*(two +&
      z01)*z22 - Gp1*z12**3/2 + P111*z12**3)/dlt3 + (-30*T*z01*z22 - 10*U*z12**3)/dlt5 +&
      (12*T*(z01**2 - z1)*z12 + 6*U*(one + 2*z01)*z22 - Gp12*z12**3 - 4*V1*z12**3)/dlt4)

      redcoeff(3) = 20*T*z22/dlt6 + (-10*U*z22 + 10*T*z12*(z12-2*z01))/dlt5+(Gp1*(-two -&
      (z0 - z2)**2/(dlt + one)))/12 +(Gp12*(one+(z0-z2)**2/(dlt+one)))/12+(U*(2*z0**2  +&
      z2*(7*(one+dlt) + 3*(one + dlt)*z1 - z2-3*dlt*z2)-z0*(five+5*dlt+3*(one +dlt)*z1 +&
      z2-3*dlt*z2)))/(6*(one+dlt)**2)+(V1*(-((3*one+dlt)*z0**2)+z0*(five+ 5*dlt+3*(one +&
      dlt)*z1+3*z2-dlt*z2)+z2*(2*dlt*z2-7*(one+dlt)-3*(one+dlt)*z1)))/(6*(one+dlt)**2) +&
      (T*(-((-one + dlt)*z0**2) + z2*(7*(one+dlt) + 3*(one + dlt)*z1 - 2*z2 -4*dlt*z2) +&
      z0*(-5*(one+dlt)- 3*(one +dlt)*z1+z2+5*dlt*z2)))/(6*(one+dlt)**2)+(P111*(-2*(two +&
      dlt)*z0**2 + z2*(-7*(one + dlt) - 3*(one + dlt)*z1+(dlt-one)*z2)+z0*(5*(one+dlt) +&
      3*(one + dlt)*z1 + (five + dlt)*z2)))/(6*(one + dlt)**2) + (-Gp1*z22/2 +P111*z22 -&
      (Gp12*z12*(-2*(one + z0) + 3*z1 - z2))/2 + V1*z12*(7*one + z0 - 3*z1 + 2*z2)     +&
      (T*(3 + 8*z0**2 - 9*z1 + 3*z1**2 + z0*(11*one - 6*z1 - 10*z2) -10*z2+5*z2**2))/3 +&
      (U*(-10*one - z0**2 - 9*z1**2 + z0*(-11*one + 15*z1 - 13*z2) - 17*z2 + 5*z2**2   +&
      3*z1*(12*one + z2)))/3)/dlt3+((U*(3*one+z0*(5+3*z12)-(7*one+3*z1)*z2+3*z2**2))/6 +&
      (T*(two+z0**2+z0*(five+3*z1-5*z2)-(7*one+3*z1)*z2+4*z2**2))/6+(V1*(z0**2 - 6*one +&
      (7*one + 3*z1)*z2 - 2*z2**2 + z0*(-5*one - 3*z1 + z2)))/6 +(P111*(11*one+2*z0**2 +&
      (7*one + 3*z1)*z2 - z2**2 - z0*(five+3*z1+z2)))/6+(Gp12*(-z0**2+3*z1-(z2-two)*z2 +&
      z0*(-one + 2*z2)))/12 + (Gp1*(-two+z0**2+9*z1-4*z2+z2**2-z0*(one+2*z2)))/12)/dlt +&
      (-(Gp1*z12*(-3*one -z0+2*z1-z2))/4-(P111*z12*(5*one-z0+z2))/2+(V1*(-17*one+z0**2 +&
      z0*(-4*one + 3*z1 - 5*z2) - 3*z1*(-9*one + z2) -19*z2+4*z2**2))/6+(U*(-6-2*z0**2 +&
      z0*(-3*z1 + 7*(-2*one+z2))+13*z2-5*z2**2+3*z1*(3*one+z2)))/6+(Gp12*(-one - z0**2 -&
      6*z1**2 +z0*(-two+9*z1-7*z2)-5*z2+2*z2**2+3*z1*(five+z2)))/12+(T*(-3*(one+z0**2) +&
      11*z2 - 6*z2**2 + 3*z1*(one + z2) +z0*(-10*one-3*z1+9*z2)))/6)/dlt2+(-(Gp12*z22) -&
      4*V1*z22 + U*z12*(13*one + 7*z0 - 12*z1 +5*z2)+(T*(11*(one+z0**2)+36*z1**2+25*z2 -&
      10*z2**2 - 3*z1*(21*one + 5*z2) + z0*(22*one - 57*z1 + 35*z2)))/3)/dlt4

      redcoeff(4) = ((Gp1*(2 + z0 - (z0 - z1)**2-5*z1))/12+(P111*(-11*one-2*z0**2-7*z1 -&
      2*z1**2 + z0*(five + 4*z1)))/6)/dlt - (20*T*z22)/dlt6 + ((Gp12*(z01**2-8*z1))/12 +&
      (V1*(17*one + 4*z0 - (z0 - z1)**2 - 8*z1))/6 - (Gp1*(two+z01)*z12)/4+(P111*(five -&
      z0 + z1)*z12)/2)/dlt2 + ((U*(10*one + 11*z0 + (z0 -z1)**2-19*z1))/3-Gp12*z01*z12 +&
      V1*(-7*one - z0 + z1)*z12 +Gp1*z22/2-P111*z22)/dlt3+(20*T*z01*z12+10*U*z22)/dlt5 +&
      ((T*(38*z1-11*(one+2*z0+(z0-z1)**2)))/3-U*(6*one+7*z01)*z12+Gp12*z22+4*V1*z22)/dlt4

      redcoeff(5) = z0*(10*T*z12/(3*dlt4) + (-(Gp12*z12) - 4*V1*z12+T*(7*one+3*(z0+z1) -&
      6*z2) + U*(8*one + 2*z0 + 3*z1 - 5*z2))/(6*dlt2) + ((-5*U*z12)/3 -(5*T*(one + z0 -&
      z2))/3)/dlt3 + ((Gp1 - Gp12)*(z0 - z2))/(12*(dlt + one)) + (T*(5*(one+dlt)+(-one +&
      dlt)*z0 + 3*(one + dlt)*z1 -2*z2-4*dlt*z2))/(6*(dlt+one)**2)+(P111*(-5*(one+dlt) +&
      2*(two + dlt)*z0 - 3*(one + dlt)*z1 - z2 + dlt*z2))/(6*(dlt+one)**2)+(V1*((3*one +&
      dlt)*z0 - (one + dlt)*(five + 3*z1) + 2*dlt*z2))/(6*(dlt +one)**2)-(U*(2*z0-(one +&
      dlt)*(five + 3*z1) + z2 + 3*dlt*z2))/(6*(dlt + one)**2) + (-Gp1*z12/2 + P111*z12 -&
      T*(five + z0 + 3*z1 - 4*z2) - U*(five + 3*z12) + (Gp12*(one +z0-z2))/2-V1*(-five +&
      z0 - 3*z1 + 2*z2))/(6*dlt))

      redcoeff(6) = z1*(((P111*(7*one - 2*(z0 - z1)))/6 + (Gp1*(-4*one - z0 + z1))/12  -&
      T*z1**2/(3*z12**3) - U*z1/(2*z22) + V1/z12 - Gp12*z1/(4*z12))/dlt -10*T*z12/dlt5 +&
      ((-5*Gp12*z01)/12+(V1*(-19*one-z0+z1))/6+T*z01*z1/(2*z22)+(U*(2*z01-z1))/(2*z12) +&
      (Gp1/4 -P111/2)*z12)/dlt2+((25*T*z01)/3+5*U*z12)/dlt4+((U*(-17*one-8*(z0-z1)))/3 -&
      (T*(one + 2*z0 + z0**2 - 3*z1 - 2*z0*z1 + z1**2))/z12 + Gp12*z12/2 + 2*V1*z12)/dlt3)

      redcoeff(7) = z2*((10*T*z12)/dlt5 + (-Gp1*z12/4 +P111*z12/2-(U*(2*z0*(3*one+z12) -&
      z1*(9*one + 2*z12) + (two+z12)*(3*one+5*z12)))/(6*z12)+(T*(3*z1**2+3*z1*(-one-z0 +&
      z12 + z22) - z12*(3*one + 11*z12 +6*z22+3*z0*(one+z12))))/(6*z22)+(V1*(19*one+z0 +&
      3*z1 - 4*z2))/6 + (Gp12*(5*(one + z0) - 3*z1 - 2*z2))/12)/dlt2 - ((Gp1-Gp12)*(z0 -&
      z2))/(12*(dlt + one)) + (V1*(-((3*one + dlt)*z0) + (one + dlt)*(7*one + 3*z1)    -&
      2*dlt*z2))/(6*(dlt + one)**2) + (P111*(7*(one+dlt)-2*(two+dlt)*z0+3*(one+dlt)*z1 +&
      z2 - dlt*z2))/(6*(dlt + one)**2) + (U*(-7*(one+dlt) + 2*z0-3*(one + dlt)*z1 + z2 +&
      3*dlt*z2))/(6*(dlt + one)**2) + ((-5*T*(5*(one+z0)-3*z1-2*z2))/3 - 5*U*z12)/dlt4 +&
      (T*(z0 - (one + dlt)*(7*one + 3*z1) + 2*z2 + dlt*(-z0 + 4*z2)))/(6*(dlt+one)**2) +&
      (-(V1*(6*one + z12*(7*one -z0+z1+2*z12)))/(6*z12)+(U*(3*z1+z12*(3*one+z12*(7*one +&
      3*z12))))/(6*z22) + (T*(2*z1**2 - z1*z12*(-two + z22) + z22*(two + z12*(7*one+z0 +&
      4*z12))))/(6*z12**3) + (Gp1*(4*one + z0 - z2))/12+(P111*(-7*one+2*z0-3*z1+z2))/6 +&
      (Gp12*((two + z0 - z2)*z2 + z1*(one-z0+z2)))/(12*z12))/dlt+(-Gp12*z12/2-2*V1*z12 +&
      (U*(17*one + 8*z0 - 3*z1 - 5*z2))/3 + (T*(3*one+3*z0**2+z1+3*z1**2-10*z2-5*z1*z2 +&
      5*z2**2 - z0*(-6*one + z1 + 5*z2)))/(3*z12))/dlt3)

      redcoeff(8) = (-7*Gp1)/36 + Gp12/18 - 10*T*z12/(3*dlt4) - ((P111 - T - U+V1)*(z0 +&
      z2))/(6*(dlt + one)) +((5*U*z12)/3-(5*T*(-two-2*z0+3*z1+z2))/6)/dlt3+(Gp12*z12/6 +&
      (2*V1*z12)/3 + (U*(-16*one - 4*z0 +9*z1+5*z2))/12+(T*(-3*z1**2+z1*(-one+7*z0+z2) -&
      z2*(-7*one+z0+4*z2)))/(12*z12))/dlt2+(Gp1*z12/12-P111*z12/6+(V1*(-five+z0+z2))/6 +&
      (Gp12*(-two - 2*z0 + 3*z1 + z2))/24 + (U*(z2*(-five + 2*z0 + 2*z2)-z1*(one+2*(z0 +&
      z2))))/(12*z12)+(T*((five-3*(z0+z2))*z2**2-z1**2*(4*one+3*(z0+z2))+z1*z2*(-7*one +&
      6*(z0 + z2))))/(18*z22))/dlt

    end if

  end if

end subroutine tch_triangle_exact

! ******************************************************************************
subroutine tch_triangle_expand(r,G_in,msq,masses_ind,dlt,psq,p1,p2,&
& redcoeff)
! ******************************************************************************
! OpenLoops Reduction. Reduction of three point integrals with one
! external massless leg. Reduction formulae expanded in \delta
! ------------------------------------------------------------------------------
! G_in       = input coefficient
! r          = tensor integral rank
! msq        = array of squared mass
! dlt        = \delta parameter
! p1, p2     = external momenta in the triangle
! redcoeff   = (1-7) zero
!              (8) contracted coefficients
! masses_ind = integer for the internal masses configuration
! ******************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_error
  use ol_generic, only: to_string
  use trred,only: get_errflag_trred,reset_errflag_trred,               &
                  set_muUV2,set_muIR2,set_duv,set_dir,                 &
                  C_m00_p1,C_m00_P12,C_0mm_p1,C_0mm_P12,               &
                  C_m00_P12P12,C_m00_p1P12,C_m00_p1p1,C_m00_g,         &
                  C_0mm_P12P12,C_0mm_p1P12,C_0mm_p1p1,C_0mm_g,         &
                  C_000_p1,C_000_P12,C_000_g,C_000_p1p1,C_000_p1P12,   &
                  C_000_P12P12,C_mmm_P12,C_mmm_p1,C_mmm_g,C_mmm_p1p1,  &
                  C_mmm_p1P12,C_mmm_P12P12,C_000_P12P12P12,            &
                  C_000_p1P12P12,C_000_p1p1P12,C_000_p1p1p1,C_000_gp1, &
                  C_000_gP12,C_m00_P12P12P12,C_m00_p1P12P12,           &
                  C_m00_p1p1P12,C_m00_p1p1p1,C_m00_gp1,C_m00_gP12,     &
                  C_0mm_P12P12P12,C_0mm_p1P12P12,C_0mm_p1p1P12,        &
                  C_0mm_p1p1p1,C_0mm_gp1,C_0mm_gP12,C_mmm_P12P12P12,   &
                  C_mmm_p1P12P12,C_mmm_p1p1P12,C_mmm_p1p1p1,C_mmm_gp1, &
                  C_mmm_gP12,C_m0m1m1_P12,C_m0m1m1_p1,C_m0m1m1_P12P12, &
                  C_m0m1m1_p1p1,C_m0m1m1_p1P12,C_m0m1m1_g,             &
                  C_m0m1m1_P12P12P12,C_m0m1m1_p1P12P12,C_m0m1m1_gP12,  &
                  C_m0m1m1_p1p1P12,C_m0m1m1_p1p1p1,C_m0m1m1_gp1

  use ol_loop_parameters_decl_/**/REALKIND, only: mu2_UV, mu2_IR
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_UV, de1_IR
  implicit none
  complex(REALKIND), intent(in)  :: G_in(:)
  complex(REALKIND), intent(in)  :: msq(0:2),p1(1:5),p2(1:5)
  integer,           intent(in)  :: r,masses_ind
  real(REALKIND),    intent(in)  :: dlt,psq
  complex(REALKIND), intent(out) :: redcoeff(1:8)

  complex(REALKIND) :: psq_delta,psq_delta2,dlt2,dlt3,dlt4
  complex(REALKIND) :: contr_p1,contr_p2,contr_P12
  complex(REALKIND) :: contr_P12P12,contr_p1p1,contr_p1P12,contr_g,contr_gp1,contr_gP12
  complex(REALKIND) :: contr_P12P12P12,contr_p1p1p1,contr_p1P12P12,contr_p1p1P12
  complex(REALKIND) :: p1p1(6:15),p2p2(6:15),p1p2(6:15),p1P12(6:15),P12P12(6:15),P12_vec(1:4)
  complex(REALKIND) :: p1p1p1(16:35),p1p1P12(16:35),p1P12P12(16:35),P12P12P12(16:35)
  integer :: i,j,n,k,errflag

  redcoeff = 0._/**/REALKIND
  psq_delta = psq*dlt
  psq_delta2 = psq_delta*dlt
  dlt2 = dlt**2
  dlt3 = dlt**3
  dlt4 = dlt**4

  call set_muUV2(mu2_UV)
  call set_muIR2(mu2_IR)
#ifdef PRECISION_dp
  call set_duv(de1_UV)
  call set_dir(de1_IR)
#else
  call set_duv(real(de1_UV,kind=REALKIND))
  call set_dir(real(de1_IR,kind=REALKIND))
#endif

  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-1 part ---------------------*
  !------------------------------------------------------------------------*
  if(r == 1) then
    contr_p1  = SUM(G_in(2:5)*p1(1:4))
    contr_p2  = SUM(G_in(2:5)*p2(1:4))
    contr_P12 = contr_p1 - contr_p2

    call reset_errflag_trred
    !! masses: (0,0,0)
    if (masses_ind == 0) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = contr_P12*C_000_P12(psq,dlt) &
                  + contr_p1*C_000_p1(psq,dlt)
    !! masses: (m0,0,0)
    else if (masses_ind == 1) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = contr_P12*C_m00_P12(psq,msq(0),dlt) &
                  + contr_p1*C_m00_p1(psq,msq(0),dlt)

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = contr_P12*C_mmm_P12(psq,msq(0),dlt)  &
                  + contr_p1*C_mmm_p1(psq,msq(0),dlt)

    !! masses: (0,m,m)
    else if (masses_ind == 5) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = contr_P12*C_0mm_P12(psq,msq(1),dlt)  &
                  + contr_p1*C_0mm_p1(psq,msq(1),dlt)

    !! masses: (m0,m1,m1)
    else if (masses_ind == 8) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = contr_P12*C_m0m1m1_P12(psq,msq(0),msq(1),dlt)  &
                  + contr_p1*C_m0m1m1_p1(psq,msq(0),msq(1),dlt)
    end if
    errflag = get_errflag_trred()
    if (errflag .ne. 0) then
      call ol_error("In TrRed for rank 1 and masses id " // &
                    trim(to_string(masses_ind)) // ". Trred flag: " // &
                    to_string(errflag))
      contr_p1 = 0._/**/REALKIND
      redcoeff = redcoeff/contr_p1
    end if
  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-2 part ---------------------*
  !------------------------------------------------------------------------*
  else if(r == 2) then
    n = 6
    do i = 1, 4
      do j = i, 4
        p1p1(n) = p1(i)*p1(j)
        p1P12(n) = p1(i)*(p1(j)-p2(j)) + p1(j)*(p1(i)-p2(i))
        P12P12(n) = (p1(i)-p2(i))*(p1(j)-p2(j))
        n = n + 1
      end do
    end do

    contr_p1p1 = sum(p1p1*G_in(6:15))
    contr_p1P12 = sum(p1P12*G_in(6:15))
    contr_P12P12 = sum(P12P12*G_in(6:15))
    contr_g = 2*(G_in(7)-G_in(14))

    call reset_errflag_trred
    !! masses: (0,0,0)
    if(masses_ind == 0) then
        redcoeff(:7) = 0._/**/REALKIND
        redcoeff(8) = + contr_P12P12*C_000_P12P12(psq,dlt) &
                      + contr_p1p1*C_000_p1p1(psq,dlt)     &
                      + contr_p1P12*C_000_p1P12(psq,dlt)   &
                      + contr_g*C_000_g(psq,dlt)

    !! masses: (m0,0,0)
    else if(masses_ind == 1) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12*C_m00_P12P12(psq,msq(0),dlt) &
                    + contr_p1p1*C_m00_p1p1(psq,msq(0),dlt)     &
                    + contr_p1P12*C_m00_p1P12(psq,msq(0),dlt)   &
                    + contr_g*C_m00_g(psq,msq(0),dlt)

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12*C_mmm_P12P12(psq,msq(0),dlt) &
                    + contr_p1p1*C_mmm_p1p1(psq,msq(0),dlt)     &
                    + contr_p1P12*C_mmm_p1P12(psq,msq(0),dlt)   &
                    + contr_g*C_mmm_g(psq,msq(0),dlt)

    !! masses: (0,m,m)
    else if(masses_ind == 5) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12*C_0mm_P12P12(psq,msq(1),dlt) &
                    + contr_p1p1*C_0mm_p1p1(psq,msq(1),dlt)     &
                    + contr_p1P12*C_0mm_p1P12(psq,msq(1),dlt)   &
                    + contr_g*C_0mm_g(psq,msq(1),dlt)

    !! masses: (m0,m1,m1)
    else if(masses_ind == 8) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12*C_m0m1m1_P12P12(psq,msq(0),msq(1),dlt) &
                    + contr_p1p1*C_m0m1m1_p1p1(psq,msq(0),msq(1),dlt)     &
                    + contr_p1P12*C_m0m1m1_p1P12(psq,msq(0),msq(1),dlt)   &
                    + contr_g*C_m0m1m1_g(psq,msq(0),msq(1),dlt)
    end if

    errflag = get_errflag_trred()
    if (errflag .ne. 0) then
      call ol_error("In TrRed for rank 2 and masses id " // &
                    trim(to_string(masses_ind)) // ". Trred flag: " // &
                    to_string(errflag))
      contr_p1 = 0._/**/REALKIND
      redcoeff = redcoeff/contr_p1
    end if

  !------------------------------------------------------------------------*
  !--------------------- Formulae for the rank-3 part ---------------------*
  !------------------------------------------------------------------------*
  else if(r == 3) then
    n = 16

    P12_vec(1:4) = p1(1:4) - p2(1:4)

    do i = 1, 4
      do j = i, 4
        do k = j, 4
          p1p1p1(n) = p1(i)*p1(j)*p1(k)
          P12P12P12(n) = P12_vec(i)*P12_vec(j)*P12_vec(k)
          p1p1P12(n) = p1(i)*p1(j)*P12_vec(k) + &
                       p1(i)*p1(k)*P12_vec(j) + &
                       p1(k)*p1(j)*P12_vec(i)
          p1P12P12(n) = P12_vec(i)*P12_vec(j)*p1(k) + &
                        P12_vec(i)*P12_vec(k)*p1(j) + &
                        P12_vec(k)*P12_vec(j)*p1(i)
          n = n + 1
        end do
      end do
    end do

    contr_p1p1p1 = sum(p1p1p1(16:35)*G_in(16:35))
    contr_p1p1P12 = sum(p1p1P12(16:35)*G_in(16:35))
    contr_p1P12P12 = sum(p1P12P12(16:35)*G_in(16:35))
    contr_P12P12P12 = sum(P12P12P12(16:35)*G_in(16:35))

    contr_gp1 = 2*((2*G_in(17) - G_in(24))*p1(1) + &
               (2*G_in(20) - G_in(30))*p1(2) +     &
               (G_in(21) - 2*G_in(33))*p1(3) +     &
               (G_in(22) - 2*G_in(34))*p1(4))
    contr_gP12 = 2*((2*G_in(17) - G_in(24))*P12_vec(1) + &
                (2*G_in(20) - G_in(30))*P12_vec(2) +     &
                (G_in(21) - 2*G_in(33))*P12_vec(3) +     &
                (G_in(22) - 2*G_in(34))*P12_vec(4))

    call reset_errflag_trred
    !! masses: (0,0,0)
    if(masses_ind == 0) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12P12*C_000_P12P12P12(psq,dlt) &
                    + contr_p1P12P12*C_000_p1P12P12(psq,dlt)   &
                    + contr_p1p1P12*C_000_p1p1P12(psq,dlt)     &
                    + contr_p1p1p1*C_000_p1p1p1(psq,dlt)       &
                    + contr_gp1*C_000_gp1(psq,dlt)             &
                    + contr_gP12*C_000_gP12(psq,dlt)

    !! masses: (m0,0,0)
    else if(masses_ind == 1) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12P12*C_m00_P12P12P12(psq,msq(0),dlt) &
                    + contr_p1P12P12*C_m00_p1P12P12(psq,msq(0),dlt)   &
                    + contr_p1p1P12*C_m00_p1p1P12(psq,msq(0),dlt)     &
                    + contr_p1p1p1*C_m00_p1p1p1(psq,msq(0),dlt)       &
                    + contr_gp1*C_m00_gp1(psq,msq(0),dlt)             &
                    + contr_gP12*C_m00_gP12(psq,msq(0),dlt)

    !! masses: (m,m,m)
    else if(masses_ind == 4) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12P12*C_mmm_P12P12P12(psq,msq(0),dlt) &
                    + contr_p1P12P12*C_mmm_p1P12P12(psq,msq(0),dlt)   &
                    + contr_p1p1P12*C_mmm_p1p1P12(psq,msq(0),dlt)     &
                    + contr_p1p1p1*C_mmm_p1p1p1(psq,msq(0),dlt)       &
                    + contr_gp1*C_mmm_gp1(psq,msq(0),dlt)             &
                    + contr_gP12*C_mmm_gP12(psq,msq(0),dlt)
    !! masses: (0,m,m)
    else if(masses_ind == 5) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12P12*C_0mm_P12P12P12(psq,msq(1),dlt) &
                    + contr_p1P12P12*C_0mm_p1P12P12(psq,msq(1),dlt)   &
                    + contr_p1p1P12*C_0mm_p1p1P12(psq,msq(1),dlt)     &
                    + contr_p1p1p1*C_0mm_p1p1p1(psq,msq(1),dlt)       &
                    + contr_gp1*C_0mm_gp1(psq,msq(1),dlt)             &
                    + contr_gP12*C_0mm_gP12(psq,msq(1),dlt)
    !! masses: (m0,m1,m1)
    else if(masses_ind == 8) then
      redcoeff(:7) = 0._/**/REALKIND
      redcoeff(8) = + contr_P12P12P12*C_m0m1m1_P12P12P12(psq,msq(0),msq(1),dlt) &
                    + contr_p1P12P12*C_m0m1m1_p1P12P12(psq,msq(0),msq(1),dlt)   &
                    + contr_p1p1P12*C_m0m1m1_p1p1P12(psq,msq(0),msq(1),dlt)     &
                    + contr_p1p1p1*C_m0m1m1_p1p1p1(psq,msq(0),msq(1),dlt)       &
                    + contr_gp1*C_m0m1m1_gp1(psq,msq(0),msq(1),dlt)             &
                    + contr_gP12*C_m0m1m1_gP12(psq,msq(0),msq(1),dlt)
    end if

    errflag = get_errflag_trred()
    if (errflag .ne. 0) then
      call ol_error("In TrRed for rank 3 and masses id " // &
                    trim(to_string(masses_ind)) // ". Trred flag: " // &
                    to_string(errflag))
      contr_p1 = 0._/**/REALKIND
      redcoeff = redcoeff/contr_p1
    end if

  end if

end subroutine tch_triangle_expand

! ********************************************************************************
subroutine t_channel_triangle_reduction(perm,sdlt,mom,msq,Gin,Gout_A,Gout_A0,Gout_A1,&
Gout_A2, Gout_R1,A0msq,A0_0,A0_1,A0_2)
! --------------------------------------------------------------------------------
! Reduction of a t-channel triangle topology with one external massless leg
! --------------------------------------------------------------------------------
! perm: permutation. It gives the permutation order for the pinched subtopologies
! Gin: input openloops coefficient. It can be rank-1,2,3
! mom: external momenta of the triangle
! msq: array of squared internal masses
! Gout_A : coefficient of the scalar triangle
! Gout_Ai: coefficient of the scalar bubble, Di-pinch, i=0,1,2
! Gout_R1: rational terms
! ********************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  implicit none
  integer, intent(in) :: perm(3)
  real(REALKIND), intent(in) :: sdlt
  complex(REALKIND), intent(inout)  :: Gin(:)
  complex(REALKIND), intent(in) :: msq(0:2), mom(5,3)
  complex(REALKIND), intent(out) :: Gout_A(1), Gout_A0(1), Gout_A1(1), Gout_A2(1)
  complex(REALKIND), intent(out) :: Gout_R1
  complex(REALKIND), optional, intent(in) :: A0msq(:)
  complex(REALKIND), optional, intent(out) :: A0_0(1), A0_1(1), A0_2(1)
  complex(REALKIND) :: k1(5), k2(5), k12(5), masses(0:2), tadp0(1), tadp1(1), tadp2(1)

  k1  = mom(:,1)
  k2  = mom(:,2)
  k12 = mom(:,3)

  !! If needed a loop momentum shift is performed in order to set the D0 propagator in
  !! the opposite position to the external massless leg.
  if(perm(1) == 2) then

    call G_TensorShift_otf(Gin,k12(1:4))

    if(perm(2) == 0) then
      masses = (/msq(2),msq(0),msq(1)/)
      call triangle_zero_leg(k12,k2,sdlt,masses,Gin,Gout_A,Gout_A2,Gout_A0,Gout_A1,&
                             Gout_R1,tadp0,tadp1,tadp2)

    else
      masses = (/msq(2),msq(1),msq(0)/)
      call triangle_zero_leg(k2,k12,sdlt,masses,Gin,Gout_A,Gout_A2,Gout_A1,Gout_A0,&
                             Gout_R1,tadp0,tadp1,tadp2)

    end if

  else if(perm(1) == 0) then

    if(perm(2) == 2) then
      masses = (/msq(0),msq(2),msq(1)/)
      call triangle_zero_leg(k12,k1,sdlt,masses,Gin,Gout_A,Gout_A0,Gout_A2,Gout_A1,&
                             Gout_R1,tadp0,tadp1,tadp2)

    else
      masses = (/msq(0),msq(1),msq(2)/)
      call triangle_zero_leg(k1,k12,sdlt,masses,Gin,Gout_A,Gout_A0,Gout_A1,Gout_A2,&
                             Gout_R1,tadp0,tadp1,tadp2)

    end if

  else if(perm(1) == 1) then

    call G_TensorShift_otf(Gin,k1(1:4))

    if(perm(2) == 2) then
      masses = (/msq(1),msq(2),msq(0)/)
      call triangle_zero_leg(k2,k1,sdlt,masses,Gin,Gout_A,Gout_A1,Gout_A2,Gout_A0,&
                             Gout_R1,tadp0,tadp1,tadp2)

    else
      masses = (/msq(1),msq(0),msq(2)/)
      call triangle_zero_leg(k1,k2,sdlt,masses,Gin,Gout_A,Gout_A1,Gout_A0,Gout_A2,&
                             Gout_R1,tadp0,tadp1,tadp2)

    end if
  end if

  if(present(A0msq)) then
    if(present(A0_2)) then
      call tadpole_assignment(perm,A0msq,masses,tadp0,tadp1,tadp2,A0_0,A0_1,A0_2)
    else if(present(A0_1)) then
      call tadpole_assignment(perm,A0msq,masses,tadp0,tadp1,tadp2,A0_0,A0_1)
    else if(present(A0_0)) then
      call tadpole_assignment(perm,A0msq,masses,tadp0,tadp1,tadp2,A0_0)
    end if
  end if

end subroutine t_channel_triangle_reduction

!******************************************************************************
subroutine tadpole_assignment(perm,A0_msq,masses,tadp0,tadp1,tadp2,A0_0,A0_1,A0_2)
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  integer, intent(in) :: perm(3)
  complex(REALKIND), intent(in) :: A0_msq(:)
  complex(REALKIND), intent(in) :: masses(0:2)
  complex(REALKIND), intent(in) :: tadp0(1), tadp1(1), tadp2(1)
  complex(REALKIND), optional, intent(out) :: A0_0(1), A0_1(1), A0_2(1)
  complex(REALKIND) :: zero = 0._/**/REALKIND

  if(present(A0_0)) A0_0 = zero
  if(present(A0_1)) A0_1 = zero
  if(present(A0_2)) A0_2 = zero

  !! one tadpole
  if(size(A0_msq) == 1) then
    A0_0 = A0_0 + tadp0
  !! two tadpoles
  else if(size(A0_msq) == 2) then
    if(A0_msq(1) == A0_msq(2) .and. A0_msq(1) == zero) then
      A0_0 = zero
      A0_1 = zero
    else if(A0_msq(1) == A0_msq(2)) then
      A0_0 = A0_0 + tadp0
    else if(A0_msq(1) == zero) then
      A0_1 = A0_1 + tadp0
    else if(A0_msq(2) == zero) then
      A0_0 = A0_0 + tadp0
    else
      if(masses(0) == zero) then
        if(masses(1) == A0_msq(1)) then
          A0_0 = A0_0 + tadp0
          A0_1 = A0_1 + tadp1
        else
          A0_0 = A0_0 + tadp1
          A0_1 = A0_1 + tadp0
        end if
      else if(masses(0) == A0_msq(1)) then
        A0_0 = A0_0 + tadp0
        A0_1 = A0_1 + tadp1
      else
        A0_0 = A0_0 + tadp1
        A0_1 = A0_1 + tadp0
      end if
    end if
  !! three tadpoles
  else if(size(A0_msq) == 3) then
    if(A0_msq(1) == A0_msq(2) .and. A0_msq(2) == A0_msq(3) .and. &
      A0_msq(1) == zero) then
      A0_0 = zero
      A0_1 = zero
      A0_2 = zero
    else if(A0_msq(1) == A0_msq(2) .and. A0_msq(2) == A0_msq(3)) then
      A0_0 = A0_0 + tadp0
    else if(A0_msq(2) == A0_msq(3) .and. A0_msq(2) == zero) then
      A0_0 = A0_0 + tadp0
    else if(A0_msq(1) == A0_msq(3) .and. A0_msq(1) == zero) then
      A0_1 = A0_1 + tadp0
    else if(A0_msq(1) == A0_msq(2) .and. A0_msq(1) == zero) then
      A0_2 = A0_2 + tadp0
    else if(A0_msq(3) == zero) then
      if(A0_msq(1) == A0_msq(2)) then
        A0_0 = A0_0 + tadp0
      else
        if(masses(0) == zero) then
          if(masses(1) == A0_msq(1)) then
            A0_0 = A0_0 + tadp0
            A0_1 = A0_1 + tadp1
          else
            A0_0 = A0_0 + tadp1
            A0_1 = A0_1 + tadp0
          end if
        else
          if(masses(0) == A0_msq(1)) then
            A0_0 = A0_0 + tadp0
            A0_1 = A0_1 + tadp1
          else
            A0_0 = A0_0 + tadp1
            A0_1 = A0_1 + tadp0
          end if
        end if
      end if
    else if(A0_msq(2) == zero) then
      if(A0_msq(1) == A0_msq(3)) then
        A0_0 = A0_0 + tadp0
      else
        if(masses(0) == zero) then
          if(masses(1) == A0_msq(1)) then
            A0_0 = A0_0 + tadp0
            A0_2 = A0_2 + tadp1
          else
            A0_0 = A0_0 + tadp1
            A0_2 = A0_2 + tadp0
          end if
        else
          if(masses(0) == A0_msq(1)) then
            A0_0 = A0_0 + tadp0
            A0_2 = A0_2 + tadp1
          else
            A0_0 = A0_0 + tadp1
            A0_2 = A0_2 + tadp0
          end if
        end if
      end if
    else if(A0_msq(1) == zero) then
      if(A0_msq(2) == A0_msq(3)) then
        A0_1 = A0_1 + tadp0
      else
        if(masses(0) == zero) then
          if(masses(1) == A0_msq(2)) then
            A0_1 = A0_1 + tadp0
            A0_2 = A0_2 + tadp1
          else
            A0_1 = A0_1 + tadp1
            A0_2 = A0_2 + tadp0
          end if
        else
          if(masses(0) == A0_msq(2)) then
            A0_1 = A0_1 + tadp0
            A0_2 = A0_2 + tadp1
          else
            A0_1 = A0_1 + tadp1
            A0_2 = A0_2 + tadp0
          end if
        end if
      end if
    else
      if(perm(1) == 2) then
        A0_2 = A0_2 + tadp0
        if(perm(2) == 0) then
          A0_0 = A0_0 + tadp1
          A0_1 = A0_1 + tadp2
        else
          A0_1 = A0_1 + tadp1
          A0_0 = A0_0 + tadp2
        end if
      else if(perm(1) == 0) then
        A0_0 = A0_0 + tadp0
        if(perm(2) == 1) then
          A0_1 = A0_1 + tadp1
          A0_2 = A0_2 + tadp2
        else
          A0_2 = A0_2 + tadp1
          A0_1 = A0_1 + tadp2
        end if
      else if(perm(1) == 1) then
        A0_1 = A0_1 + tadp0
        if(perm(2) == 2) then
          A0_2 = A0_2 + tadp1
          A0_0 = A0_0 + tadp2
        else
          A0_0 = A0_0 + tadp1
          A0_2 = A0_2 + tadp2
        end if
      end if
    end if
  end if

end subroutine tadpole_assignment


!******************************************************************************
subroutine triangle_zero_leg(p1,p2,dlt,msq,Gin_A,Gout_A,Gout_A0,Gout_A1,Gout_A2,&
Gout_R1,A0_0,A0_1,A0_2)
!******************************************************************************
! Reduction of 3-point integrals with one massless external leg.
! ------------------------------------------------------------------------------
! p1, p2  = external momenta of the 3-point function
! dlt     = \delta parameter
! msq     = array of internal masses squared. (m0^2, m1^2, m2^2)
! Gin_A   = input  OpenLoops coefficient of the rank-2   3-point function
! Gout_A  = output OpenLoops coefficient of the scalar 3-point function
! Gout_Ai = output OpenLoops coefficient of the scalar bubble. i-th pinch
! Gout_R1 = rational part
! A0_i    = coefficients of the tadpoles
! ------------------------------------------------------------------------------
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_error
  use ol_loop_parameters_decl_/**/DREALKIND, only: polecheck_is
  use ol_parameters_decl_/**/REALKIND, only: zero
  implicit none
  complex(REALKIND), intent(in)  :: p1(1:5), p2(1:5), msq(0:2)
  real(REALKIND),    intent(in)  :: dlt
  complex(REALKIND), intent(in)  :: Gin_A(:)
  complex(REALKIND), intent(out) :: Gout_A(1), Gout_A0(1), Gout_A1(1), Gout_A2(1)
  complex(REALKIND), intent(out) :: Gout_R1
  complex(REALKIND), intent(out) :: A0_0(1), A0_1(1), A0_2(1)

  real(REALKIND) :: psq
  integer :: m_ind
  logical :: dlt_exp
  complex(REALKIND) :: redcoeff1(1:8), redcoeff2(1:8), redcoeff3(1:8)

  Gout_A  = zero
  Gout_A0 = zero
  Gout_A1 = zero
  Gout_A2 = zero
  A0_0 = zero
  A0_1 = zero
  A0_2 = zero

  redcoeff1 = zero
  redcoeff2 = zero
  redcoeff3 = zero

  if(REAL(p1(5)) < 0) then
    psq = ABS(p1(5))
  else
    psq = - ABS(p1(5))
  end if

  if (msq(1)==msq(2) .AND. msq(1)==msq(0) .AND. msq(0)==zero) then
    !! (0,0,0) masses configuration
    m_ind = 0
  else if (msq(1)==msq(2) .AND. msq(1)==zero) then
    !! (m0,0,0) masses configuration
    m_ind = 1
  else if (msq(0)==msq(2) .AND. msq(0)==zero) then
    !! (0,m1,0) masses configuration
    m_ind = 2
  else if (msq(0)==msq(1) .AND. msq(0)==zero) then
    !! (0,0,m2) masses configuration
    m_ind = 3
  else if (msq(1)==msq(2) .AND. msq(1)==msq(0)) then
    !! (m,m,m) masses configuration
    m_ind = 4
  else if (msq(1)==msq(2) .AND. msq(0)==zero) then
    !! (0,m,m) masses configuration
    m_ind = 5
  else if (msq(0)==msq(2) .AND. msq(1)==zero) then
    !! (m,0,m) masses configuration
    m_ind = 6
  else if (msq(0)==msq(1) .AND. msq(2)==zero) then
    !! (m,m,0) masses configuration
    m_ind = 7
  else if (msq(1)==msq(2)) then
    !! (m0,m1,m1) masses configuration
    m_ind = 8
  else if (msq(0)==msq(1)) then
    !! (m0,m0,m2) masses configuration
    m_ind = 9
  else if (msq(0)==msq(2)) then
    !! (m0,m1,m0) masses configuration
    m_ind = 10
  else if (msq(0)==zero) then
    !! (0,m1,m2) masses configuration
    m_ind = 11
  else if (msq(1)==zero) then
    !! (m0,0,m2) masses configuration
    m_ind = 12
  else if (msq(2)==zero) then
    !! (m0,m1,0) masses configuration
    m_ind = 13
  else
    !! (m0,m1,m2) masses configuration
    m_ind = 14
  end if

  ! local switcher to activate the expansions in the Gram-Determinant
  dlt_exp = (dlt < delta_thres) .AND. DeltaExp .AND. (REAL(p1(5)) < 0)

  ! TEMPORARY: expansions are switched off in case polecheck = 1
  !dlt_exp = dlt_exp .and. (polecheck_is == 0)

! TEMPORARY: expansions are switched off if they have not been implemented yet
  if((m_ind == 14) .or. (m_ind == 13) .or. (m_ind == 12) .or. (m_ind == 11) .or. &
     (m_ind == 10) .or. (m_ind ==  9) .or. (m_ind ==  7) .or. (m_ind ==  6) .or. &
     (m_ind ==  3) .or. (m_ind ==  2)) then
     dlt_exp = .false.
  end if

  if(.NOT. dlt_exp) then
    !! Exact reduction formulae
    if(size(Gin_A) == 5) then
      call tch_triangle_exact(1,Gin_A(1:5), msq,m_ind,dlt,psq,p1,p2,redcoeff1)

    else if(size(Gin_A) == 15) then
      call tch_triangle_exact(2,Gin_A(1:15),msq,m_ind,dlt,psq,p1,p2,redcoeff2)
      call tch_triangle_exact(1,Gin_A(1:5), msq,m_ind,dlt,psq,p1,p2,redcoeff1)

    else if(size(Gin_A) == 35) then
      call tch_triangle_exact(3,Gin_A(1:35),msq,m_ind,dlt,psq,p1,p2,redcoeff3)
      call tch_triangle_exact(2,Gin_A(1:15),msq,m_ind,dlt,psq,p1,p2,redcoeff2)
      call tch_triangle_exact(1,Gin_A(1:5), msq,m_ind,dlt,psq,p1,p2,redcoeff1)
    end if

  else
    !! Expansions, computed via the TrRed library
    if(size(Gin_A) == 5) then
      call tch_triangle_expand( 1,Gin_A(1:5) ,msq,m_ind,dlt,psq, p1,p2,redcoeff1)

    else if(size(Gin_A) == 15) then
      call tch_triangle_expand( 2,Gin_A(1:15),msq,m_ind,dlt,psq, p1,p2,redcoeff2)
      call tch_triangle_expand( 1,Gin_A(1:5) ,msq,m_ind,dlt,psq, p1,p2,redcoeff1)

    else if(size(Gin_A) == 35) then
      call tch_triangle_expand( 3,Gin_A(1:35),msq,m_ind,dlt,psq, p1,p2,redcoeff3)
      call tch_triangle_expand( 2,Gin_A(1:15),msq,m_ind,dlt,psq, p1,p2,redcoeff2)
      call tch_triangle_expand( 1,Gin_A(1:5) ,msq,m_ind,dlt,psq, p1,p2,redcoeff1)
    end if
  end if

  Gout_A(1)  = Gin_A(1)     + redcoeff1(1) + redcoeff2(1) + redcoeff3(1)
  Gout_A0(1) = Gout_A0(1)   + redcoeff1(2) + redcoeff2(2) + redcoeff3(2)
  Gout_A1(1) = Gout_A1(1)   + redcoeff1(3) + redcoeff2(3) + redcoeff3(3)
  Gout_A2(1) = Gout_A2(1)   + redcoeff1(4) + redcoeff2(4) + redcoeff3(4)
  Gout_R1    = redcoeff1(8) + redcoeff2(8) + redcoeff3(8)

  !! Tadpole contribution
  A0_0 = (redcoeff1(5) + redcoeff2(5) + redcoeff3(5))
  A0_1 = (redcoeff1(6) + redcoeff2(6) + redcoeff3(6))
  A0_2 = (redcoeff1(7) + redcoeff2(7) + redcoeff3(7))

end subroutine triangle_zero_leg

!************************************************************************************
subroutine scalar_MIs(momenta, masses2, Gsum, M2)
!************************************************************************************
  use KIND_TYPES, only: REALKIND,QREALKIND
  use ol_data_types_/**/REALKIND, only: hcl, met
  use ol_kinematics_/**/REALKIND, only: get_mass2
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres,coli_cache_use
#ifdef PRECISION_dp
  use ol_parameters_decl_/**/DREALKIND, only: a_switch, coli_cache_use
  use ol_kinematics_/**/QREALKIND, only: get_mass2_qp=>get_mass2
  use ol_loop_handling_/**/REALKIND, only: req_qp_cmp, upgrade_qp, downgrade_dp
  use ol_parameters_decl_/**/REALKIND, only: hybrid_zero_mode,hp_step_thres, &
                                             hp_err_thres,hp_switch,hybrid_dp_mode
  use ol_loop_reduction_/**/QREALKIND, only: avh_olo_interface_qp=>avh_olo_interface
#endif
  implicit none
  integer,           intent(in)    :: momenta(:),masses2(:)
  type(hcl),         intent(inout) :: Gsum
  type(met),         intent(inout) :: M2
  complex(REALKIND) :: M2add
#ifdef PRECISION_dp
  complex(QREALKIND) :: M2add_qp
#endif
  real(REALKIND)    :: error

#ifdef PRECISION_dp
  if(a_switch==1 .or. a_switch==7) then

    if (iand(Gsum%mode, hybrid_dp_mode) .ne. 0 .or. coli_cache_use .eq. 1) then
      call collier_scalars_interface(momenta, get_mass2(masses2), Gsum%cmp, M2add, error)
      if (hp_switch .eq. 1 .and. error > hp_step_thres .and. Gsum%error + error > hp_err_thres) then
        call upgrade_qp(Gsum)
      else
        M2%cmp = M2%cmp + real(M2add,kind=REALKIND)
        M2%sicount = M2%sicount + 1
      end if
    end if

    if (req_qp_cmp(Gsum)) then
      call avh_olo_interface_qp(momenta, get_mass2_qp(masses2), Gsum%cmp_qp, M2add_qp)
      M2%cmp_qp = M2%cmp_qp + real(M2add_qp,kind=QREALKIND)
      M2%sicount_qp = M2%sicount_qp + 1
    end if

  else if(a_switch==5) then

    if (iand(Gsum%mode, hybrid_dp_mode) .ne. 0) then
      call avh_olo_interface(momenta, get_mass2(masses2), Gsum%cmp, M2add)
      M2%cmp = M2%cmp + real(M2add,kind=REALKIND)
      M2%sicount = M2%sicount + 1
    end if

    if (req_qp_cmp(Gsum)) then
      call avh_olo_interface_qp(momenta, get_mass2_qp(masses2), Gsum%cmp_qp, M2add_qp)
      M2%cmp_qp = M2%cmp_qp + real(M2add_qp,kind=QREALKIND)
      M2%sicount_qp = M2%sicount_qp + 1
    end if

  end if
#else
  call avh_olo_interface(momenta, get_mass2(masses2), Gsum%cmp, M2add)
  M2%cmp = M2%cmp + real(M2add,kind=REALKIND)
#endif
end subroutine scalar_MIs

#ifdef PRECISION_dp
!************************************************************************************
subroutine collier_scalars_interface(momenta, masses2, Gsum, M2add, loss)
!************************************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx, collier_invariants
#ifdef USE_COLLIER
  use collier, only: A_cll,B_cll,C_cll,D_cll
#endif
  implicit none
  integer,                  intent(in)  :: momenta(:)
  complex(REALKIND),        intent(in)  :: masses2(:), Gsum(:)
  complex(REALKIND),        intent(out) :: M2add
  real(REALKIND), optional, intent(out) :: loss
  complex(REALKIND) :: TI(size(Gsum)), p(0:3,1:size(momenta)-1)
  complex(REALKIND) :: momenta_TI(0:3,size(momenta)-1)
  complex(REALKIND) :: sc_int(1),UV_part(1)
  real(REALKIND)    :: sc_err(1)
  integer :: np,i, k, int_mom(1:size(momenta)-1)
#ifdef USE_COLLIER
  i = momenta(1)
  int_mom(1) = i
  p(0:3,1) = L(1:4,i)
  do k = 2, size(momenta) - 1
    i = i + momenta(k)
    int_mom(k) = i
  end do

  np = size(masses2)
  select case (np)
  case (1)
    call A_cll(sc_int, UV_part, masses2(1), 0, Aerr=sc_err)
  case (2)
    call B_cll(sc_int, UV_part, collier_invariants(int_mom), masses2, 0, Berr=sc_err)
  case (3)
    call C_cll(sc_int, UV_part, collier_invariants(int_mom), masses2, 0, Cerr=sc_err)
  case (4)
    call D_cll(sc_int, UV_part, collier_invariants(int_mom), masses2, 0, Derr=sc_err)
  end select
  if (present(loss)) loss = max(0.,15.+log10(abs(sc_err(1)/real(sc_int(1)))))
  M2add = Gsum(1)*sc_int(1)
#else
  call ol_error("COLLIER library not compiled")
  M2add = 0
#endif
end subroutine collier_scalars_interface
#endif

!************************************************************************************
subroutine avh_olo_interface(momenta, masses2, Gsum, M2add)
!************************************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use avh_olo, only: olo_scale_prec, olo
  use ol_loop_parameters_decl_/**/REALKIND, only: mureg
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_UV, de1_IR, de2_i_IR
  implicit none
  complex(REALKIND), intent(in)  :: masses2(:), Gsum(:)
  integer, intent(in) :: momenta(:)
  complex(REALKIND), intent(out) :: M2add
  logical           :: tadpole_check
  complex(REALKIND) :: tadpole_mass, p1_2, p2_2, p3_2, p4_2, p12_2, p23_2
  complex(REALKIND) :: rslt(0:2), zero = 0._/**/REALKIND
  integer :: i, j, k

  tadpole_check = .FALSE.
  call olo_scale_prec(mureg)

  ! Check if the bubble given in input is a representation of a Tadpole.
  ! A Tadpole can be represented through a bubble with zero momentum:
  ! A(m) = m^2 B(p=0,m0=0,m1=m)
  if(size(momenta)==2) then
    tadpole_check = (momenta(1) == 0) .AND. (momenta(2) == 0)
  end if

  ! -- TADPOLES
  if(tadpole_check) then
    if(masses2(1)/=zero .AND. masses2(2)==zero) then
      tadpole_mass = masses2(1)
    else if(masses2(1)==zero .AND. masses2(2)/=zero) then
      tadpole_mass = masses2(2)
    else if(masses2(1)==zero .AND. masses2(2)==zero) then
      ! Tadpole with zero mass. An exact 0 is returned
      M2add = zero
      return
    end if

    call olo(rslt, tadpole_mass)
    M2add = (Gsum(1)/tadpole_mass)*(rslt(0) + rslt(1)*de1_UV)
    return
  end if

  ! -- BUBBLES
  if(size(masses2) == 2) then
    i = momenta(1)
    p1_2 = L(5,i) + L(6,i)
    call olo(rslt,p1_2,masses2(1),masses2(2))

  ! -- TRIANGLES
  else if(size(masses2) == 3) then

    i = momenta(1)
    j = momenta(2)
    p1_2 = L(5,i)   + L(6,i)
    p2_2 = L(5,j)   + L(6,j)
    p3_2 = L(5,i+j) + L(6,i+j)
    call olo(rslt,p1_2,p2_2,p3_2,masses2(1),masses2(2),masses2(3))

  ! -- BOXES
  else if(size(masses2) == 4) then
    i = momenta(1)
    j = momenta(2)
    k = momenta(3)
    p1_2  = L(5,i)     + L(6,i)
    p2_2  = L(5,j)     + L(6,j)
    p3_2  = L(5,k)     + L(6,k)
    p4_2  = L(5,i+j+k) + L(6,i+j+k)
    p12_2 = L(5,i+j)   + L(6,i+j)
    p23_2 = L(5,j+k)   + L(6,j+k)
    call olo(rslt, p1_2,p2_2,p3_2,p4_2,p12_2,p23_2,&
      masses2(1),masses2(2),masses2(3),masses2(4))
  else
    call ol_error('avh_olo_interface: integration called for a non-MI')
  end if

  if (size(masses2) .gt. 2) then
    M2add = Gsum(1)*(rslt(0) + rslt(1)*de1_IR + rslt(2)*de2_i_IR)
  else if (size(masses2) .eq. 2 .and. p1_2 .eq. 0 .and. sum(masses2) .eq. 0) then
    M2add = Gsum(1)*(rslt(0) + de1_UV - de1_IR)
  else
    M2add = Gsum(1)*(rslt(0) + rslt(1)*de1_UV)
  end if

end subroutine avh_olo_interface


#ifdef PRECISION_dp
!************************************************************************************
subroutine collier_scalar_box(momenta, masses2, rslt, err)
!************************************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx, collier_invariants
#ifdef USE_COLLIER
  use collier, only: D_cll
#endif
  implicit none
  integer,           intent(in)  :: momenta(3)
  complex(REALKIND), intent(in)  :: masses2(4)
  complex(REALKIND), intent(out) :: rslt(0:2)
  real(REALKIND),    intent(out) :: err
  complex(REALKIND) :: p(0:3,3), momenta_TI(0:3,3), sc_box(1), UV_part(1)
  real(REALKIND) :: Derr(1)
  integer :: k

  rslt = 0._/**/REALKIND
#ifdef USE_COLLIER
  p(0:3,1) = L(1:4,momenta(1))
  p(0:3,2) = L(1:4,momenta(2))
  p(0:3,3) = L(1:4,momenta(3))

  do k = 1, 3
    call LC2Std_Rep_cmplx(p(:,k), momenta_TI(:,k))
  end do

  call D_cll(sc_box, UV_part, collier_invariants(momenta), masses2, 0, Derr=Derr)
  rslt(0) = sc_box(1)
  err = abs(Derr(1)/real(rslt(0)))
#else
  call ol_error("COLLIER library not compiled")
#endif
end subroutine collier_scalar_box
#endif


!************************************************************************************
subroutine avh_olo_box(momenta, masses2, rslt)
!************************************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use avh_olo, only: olo_scale_prec, olo
  use ol_loop_parameters_decl_/**/REALKIND, only: mureg
  implicit none
  complex(REALKIND), intent(in)  :: masses2(4)
  integer,           intent(in)  :: momenta(3)
  complex(REALKIND), intent(out) :: rslt(0:2)
  complex(REALKIND) :: p1_2, p2_2, p3_2, p4_2, p12_2, p23_2
  integer :: i, j, k

  call olo_scale_prec(mureg)

  ! -- BOXES
  i = momenta(1)
  j = momenta(2)
  k = momenta(3)
  p1_2  = L(5,i)     + L(6,i)
  p2_2  = L(5,j)     + L(6,j)
  p3_2  = L(5,k)     + L(6,k)
  p4_2  = L(5,i+j+k) + L(6,i+j+k)
  p12_2 = L(5,i+j)   + L(6,i+j)
  p23_2 = L(5,j+k)   + L(6,j+k)
  call olo(rslt, p1_2,p2_2,p3_2,p4_2,p12_2,p23_2,&
  masses2(1),masses2(2),masses2(3),masses2(4))

end subroutine avh_olo_box


!******************************************************************************
subroutine compute_scalar_box(mom_ind, masses2in, RedSet, box)
!******************************************************************************
! Evaluation of a scalar box needed for the OPP reduction
! ------------------------------------------------------------------------------
! mom_ind   = indices of internal momenta
! masses2in = array of internal masses squared. (m0^2, m1^2, m2^2, m3^2)
! RedSet    = reduction set used to find the cut solution
! box       = scalar box data type. It contains the value of the box and the cut
!             solutions
! ------------------------------------------------------------------------------
  use KIND_TYPES, only: REALKIND,QREALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_data_types_/**/REALKIND, only: redset4, scalarbox
  use ol_kinematics_/**/REALKIND, only: get_LC_5,get_mass2
  use collier, only: tnten_cll
  use ol_parameters_decl_/**/DREALKIND, only: a_switch
#ifdef PRECISION_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_onshell_cut_qp=>box_onshell_cut
  use ol_kinematics_/**/QREALKIND, only: get_mass2_qp=>get_mass2
  use ol_kinematics_/**/DREALKIND, only: get_LC_5_qp
#endif
  implicit none
  integer,         intent(in)    :: mom_ind(3)
  integer,         intent(in)    :: masses2in(0:3)
  type(redset4),   intent(inout) :: RedSet
  type(scalarbox), intent(out)   :: box
  complex(REALKIND) :: masses2(0:3),p(1:5,3),q0_p(5),q0_m(5),q0_cuts(2,5),sc_box(0:2)
  integer           :: mom_box(3)
  real(REALKIND)    :: box_error,cut_error
#ifdef PRECISION_dp
  complex(QREALKIND) :: p_qp(1:5,3),q0_p_qp(5),q0_m_qp(5)
  real(QREALKIND)    :: box_error_qp,cut_error_qp
#endif

  masses2 = get_mass2(masses2in)

  p(1:5,1) = get_LC_5(mom_ind(1))
  p(1:5,2) = get_LC_5(mom_ind(2))
  p(1:5,3) = get_LC_5(mom_ind(3))

#ifdef PRECISION_dp
  if (RedSet%qp_computed) then
    p_qp(1:5,1) = get_LC_5_qp(mom_ind(1))
    p_qp(1:5,2) = get_LC_5_qp(mom_ind(2))
    p_qp(1:5,3) = get_LC_5_qp(mom_ind(3))
    call box_onshell_cut_qp(p_qp,get_mass2_qp(masses2in),RedSet%rsqp, &
                            q0_p_qp,q0_m_qp,cut_error_qp)
    q0_cuts(1,:) = q0_p_qp(:)
    q0_cuts(2,:) = q0_m_qp(:)
  else
    call box_onshell_cut(p,masses2,RedSet,q0_p,q0_m,cut_error)
    q0_cuts(1,:) = q0_p(:)
    q0_cuts(2,:) = q0_m(:)
  end if
#else
  call box_onshell_cut(p,masses2,RedSet,q0_p,q0_m,cut_error)
  q0_cuts(1,:) = q0_p(:)
  q0_cuts(2,:) = q0_m(:)
#endif


#ifdef PRECISION_dp
  if (a_switch == 5) then
  !!! OneLoop used for the scalar boxes
    mom_box(1) = mom_ind(1)
    mom_box(2) = mom_ind(2)-mom_ind(1)
    mom_box(3) = mom_ind(3)-mom_ind(2)
    call avh_olo_box(mom_box,masses2,sc_box)
    ! OneLoop  does not provide error estimator for integrals
    box_error = 10**(-17)
  else if (a_switch == 1 .or. a_switch == 7) then
  !!! Collier/DD call for scalar boxes
    call collier_scalar_box(mom_ind,masses2,sc_box,box_error)
  end if

  box_error = 15+log10(box_error)
  box = scalarbox(sc_box,          &
        q0_cuts,                   &
        cut_error,                 &
        box_error,                 &
        mom_ind=mom_ind,           &
        mom1=RedSet%redbasis%mom1, &
        mom2=RedSet%redbasis%mom2, &
        mom3=RedSet%mom3,          &
        perm=RedSet%perm,          &
        gd3=RedSet%gd3,            &
        masses2=masses2in,         &
        qp_computed=.false.,       &
        poles_qp=0,                &
        onshell_cuts_qp=0)

  if (a_switch == 1 .or. a_switch == 7) call check_box

#else
  !!! OneLoop used for the scalar boxes
  mom_box(1) = mom_ind(1)
  mom_box(2) = mom_ind(2)-mom_ind(1)
  mom_box(3) = mom_ind(3)-mom_ind(2)
  call avh_olo_box(mom_box,masses2,sc_box)
  box = scalarbox(sc_box,q0_cuts,cut_error,box_error)
#endif

#ifdef PRECISION_dp
contains
  subroutine check_box()
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_check_box,hp_err_thres,hp_qp_kinematics_init_mode
  use ofred_basis_construction_/**/QREALKIND, only: &
      construct_RedBasis_qp=>construct_RedBasis, &
      construct_p3scalars_qp=>construct_p3scalars
  use ol_data_types_/**/QREALKIND, only: basis_qp=>basis, redset4_qp=>redset4
  use ol_data_types_/**/QREALKIND, only: scalarbox_qp=>scalarbox
  use ol_loop_reduction_/**/QREALKIND, only: compute_scalar_box_qp=>compute_scalar_box
  use ol_external_decl_/**/DREALKIND, only: init_qp
  use ol_kinematics_/**/DREALKIND, only: init_qp_kinematics
  complex(QREALKIND) :: scalars(0:4),cttmp(0:4)
  type(basis_qp)     :: Redbasis_qp
  real(QREALKIND)    :: gd2,gd3
  type(scalarbox_qp) :: box_qp

  if (hp_check_box .eq. 1 .and. hp_switch .eq. 1 &
      !.and. errbox .gt. 10**(-real(15-hp_err_thres,kind=REALKIND))) then
      .and. box_error .gt. hp_err_thres) then
    if (.not. redset%qp_computed) then
      if (hp_qp_kinematics_init_mode .gt. 0 .and. .not. init_qp) call init_qp_kinematics
      call construct_RedBasis_qp(box%mom1,box%mom2,Redbasis_qp)
      call construct_p3scalars_qp(box%mom3,Redbasis_qp,scalars,gd2,gd3)
      redset%qp_computed = .true.
      redset%rsqp = redset4_qp(redbasis=Redbasis_qp, &
                               p3scalars=scalars,    &
                               perm=RedSet%perm,     &
                               mom3=RedSet%mom3,     &
                               gd2=gd2,gd3=gd3)
    end if

    call compute_scalar_box_qp(box%mom_ind,box%masses2,redset%rsqp,box_qp)
    box%qp_computed = .true.
    box%poles = box_qp%poles
    box%poles_qp = box_qp%poles
    box%onshell_cuts = box_qp%onshell_cuts
    box%onshell_cuts_qp = box_qp%onshell_cuts

  end if

  end subroutine check_box
#endif

end subroutine compute_scalar_box


! ********************************************************************
subroutine box_onshell_cut(p, m2, RedSet, q0_p, q0_m, error)
!---------------------------------------------------------------------
! Calculation of q0+, q0- for the quadruple cut
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_data_types_/**/REALKIND, only: basis, redset4
  use ol_parameters_decl_/**/REALKIND, only: hp_gamma_trig
  implicit none
  complex(REALKIND), intent(in)  :: p(1:5,1:3), m2(0:3)
  complex(REALKIND), intent(out) :: q0_p(1:5), q0_m(1:5)
  real(REALKIND), intent(out) :: error
  type(redset4), intent(in) :: RedSet

  complex(REALKIND) :: l(1:4,4), k1(1:5), k2(1:5), k3(1:5)
  complex(REALKIND) :: a1, a2
  complex(REALKIND) :: x1, x2, x3p, x3m, x4p, x4m, k3scalars(1:4)
  complex(REALKIND) :: dd0, dd1, dd2, dd3
  complex(REALKIND) :: c0, b0, A, B, R
  integer :: i, perm(3)
  type(basis) :: bas

  bas = RedSet%redbasis
  perm = RedSet%perm

  l = bas%li
  a1 = bas%alpha(1)
  a2 = bas%alpha(2)

  if(perm(1)+perm(2)==3) then ! [1,2,3]
    do i = 1, 4
      k3scalars(i) = cont_LC_cntrv(p(1:4,3),l(:,i))
    end do
    dd0 = m2(0)
    dd1 = m2(1) - p(5,1)
    dd2 = m2(2) - p(5,2)
    dd3 = m2(3) - p(5,3)
  else if(perm(1)+perm(2)==4) then ! [1,3,2]
    do i = 1, 4
      k3scalars(i) = cont_LC_cntrv(p(1:4,2),l(:,i))
    end do
    dd0 = m2(0)
    dd1 = m2(1) - p(5,1)
    dd2 = m2(3) - p(5,3)
    dd3 = m2(2) - p(5,2)
  else if(perm(1)+perm(2)==5) then ! [2,3,1]
    do i = 1, 4
      k3scalars(i) = cont_LC_cntrv(p(1:4,1),l(:,i))
    end do
    dd0 = m2(0)
    dd1 = m2(2) - p(5,2)
    dd2 = m2(3) - p(5,3)
    dd3 = m2(1) - p(5,1)
  end if

  x1 = (dd2 - dd0*(1._/**/REALKIND-a2) - dd1*a2)/bas%gamma
  x2 = (dd1 - dd0*(1._/**/REALKIND-a1) - dd2*a1)/bas%gamma

  A = k3scalars(3)/k3scalars(4)

  B = (dd0-dd3)/2 + k3scalars(1)*x1 + k3scalars(2)*x2
  B = B/k3scalars(4)

  R = (x1*x2 - dd0/bas%gamma)
  c0 = SQRT(B**2 - A*R)

  x4m = (- B + c0)/2
  x4p = (- B - c0)/2
  x3m = x4p/A
  x3p = x4m/A

  if (hp_gamma_trig) then
    error = max(log10(max(1._/**/REALKIND, ABS(x1))), log10(max(1._/**/REALKIND, ABS(x2))))
    error = max(error, max(log10(max(1._/**/REALKIND, ABS(x3m))), log10(max(1._/**/REALKIND, ABS(x3p)))))
  else
    error = max(log10(max(1._/**/REALKIND, ABS(x3m))), log10(max(1._/**/REALKIND, ABS(x3p))))
  end if

  ! The assumption is that the p_0 momentum is zero
  q0_p(1:4) = x1*l(1:4,1) + x2*l(1:4,2)
  q0_m(1:4) = q0_p(1:4)
  q0_p(1:4) = q0_p(1:4) + x3p*l(1:4,3) + x4p*l(1:4,4)
  q0_m(1:4) = q0_m(1:4) + x3m*l(1:4,3) + x4m*l(1:4,4)
  q0_p(5) = m2(0)
  q0_m(5) = m2(0)

end subroutine box_onshell_cut

! --- Calculation of the coefficients of the OPP reduction ---
! ********************************************************************
subroutine opp_numerator(Gtensor,q0,N0)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: Gtensor(:), q0(4)
  complex(REALKIND), intent(out) :: N0

  if(size(Gtensor)==1) then
    N0 = Gtensor(1)
  else if(size(Gtensor)==5) then
    N0 = Gtensor(1) + SUM(Gtensor(2:5)*q0(1:4))
  else
    call ol_error("opp_numerator: rank > 1 ")
  end if

end subroutine opp_numerator


! ********************************************************************
subroutine box_coefficient(p_offshell, m_offshell, q0_pm, Gtensor, box_coeff)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  implicit none
  complex(REALKIND), intent(in)  :: Gtensor(:)
  complex(REALKIND), intent(in)  :: p_offshell(:,:), m_offshell(:), q0_pm(2,5)
  complex(REALKIND), intent(out) :: box_coeff
  complex(REALKIND) :: q0_p(1:4), q0_m(1:4)
  complex(REALKIND) :: N0p, N0m, Gsum(size(Gtensor))
  complex(REALKIND) :: mom(1:4), prop_m, prop_p, Di_p, Di_m
  integer :: i

  !q0+, q0- for the quadruple cut
  q0_p(1:4) = q0_pm(1,1:4)
  q0_m(1:4) = q0_pm(2,1:4)

  ! Calculation of the numerator at the multiple cut
  call opp_numerator(Gtensor,q0_p,N0p)
  call opp_numerator(Gtensor,q0_m,N0m)

  Di_p = 1._/**/REALKIND
  Di_m = 1._/**/REALKIND

  do i = 1, size(m_offshell)
    prop_p = p_offshell(5,i) + q0_pm(1,5) + &
             2*cont_LC_cntrv(p_offshell(1:4,i),q0_p(1:4))
    prop_p = prop_p - m_offshell(i)

    prop_m = p_offshell(5,i) + q0_pm(2,5) + &
             2*cont_LC_cntrv(p_offshell(1:4,i),q0_m(1:4))
    prop_m = prop_m - m_offshell(i)

    Di_p = Di_p*prop_p
    Di_m = Di_m*prop_m
  end do

  box_coeff = (N0p/Di_p + N0m/Di_m)/2

end subroutine box_coefficient



!*************************************************************************************
subroutine OPP_reduction(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
!-------------------------------------------------------------------------------------
!*************************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_parameters_decl_/**/DREALKIND, only: hp_max_err
  implicit none
  integer,         intent(in) :: rank
  type(hcl),       intent(inout) :: Gtensor
  integer,         intent(in)  :: momenta(:), masses(:)
  type(met),       intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)

  if(size(masses) == 5) then
    call reduction_5points_ofr(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
  else if(size(masses) == 6) then
    call reduction_6points_ofr(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
  else if(size(masses) == 7) then
    call reduction_7points_ofr(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
  else if(size(masses) == 8) then
    call reduction_8points_ofr(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
  else if(size(masses) == 9) then
    call reduction_9points_ofr(rank, momenta, masses, Gtensor, M2, scboxes, all_scboxes)
  end if
  if (Gtensor%error > hp_max_err) hp_max_err = Gtensor%error

end subroutine OPP_reduction


!*************************************************************************************
subroutine TI_reduction_1(rank, momenta, masses2, Gtensor, M2add, scboxes, all_scboxes)
!-------------------------------------------------------------------------------------
! Tensor integral manager. Momenta passed as 4-dimensional vectors in light-cone
!*************************************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl
  use ol_parameters_decl_/**/DREALKIND, only: hp_max_err
  implicit none
  integer,           intent(in)  :: rank
  type(hcl), intent(in) :: Gtensor
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(:)
  complex(REALKIND),    intent(out) :: M2add
  integer,         intent(in), optional :: scboxes(:)
  type(scalarbox), intent(in), optional :: all_scboxes(:)
  real(REALKIND) :: box_error, error_OPP

  if(size(masses2) == 5) then
    call reduction_5points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  else if(size(masses2) == 6) then
    call reduction_6points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  else if(size(masses2) == 7) then
    call reduction_7points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  else if(size(masses2) == 8) then
    call reduction_8points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  else if(size(masses2) == 9) then
    call reduction_9points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  else if(size(masses2) == 10) then
    call reduction_10points(rank, momenta, masses2, Gtensor%cmp, M2add, scboxes, all_scboxes, box_error)
  end if
  error_OPP = Gtensor%error + box_error
  if (error_OPP > hp_max_err) hp_max_err = error_OPP

end subroutine TI_reduction_1


#ifdef PRECISION_dp
subroutine upgrade_scalar_box(scbox)
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_data_types_/**/QREALKIND, only: scalarbox_qp=>scalarbox
  use ol_data_types_/**/QREALKIND, only: basis_qp=>basis, redset4_qp=>redset4
  use ofred_basis_construction_/**/QREALKIND, only: construct_RedBasis_qp => construct_RedBasis, &
                                                    construct_p3scalars_qp=>construct_p3scalars
  use ol_loop_reduction_/**/QREALKIND, only: compute_scalar_box_qp=>compute_scalar_box
  implicit none
  type(scalarbox), intent(inout) :: scbox
  type(scalarbox_qp) :: scbox_qp
  type(basis_qp)     :: Redbasis_qp
  type(redset4_qp)   :: RedSet_qp
  real(QREALKIND)    :: gd2,gd3
  complex(QREALKIND) :: scalars(0:4)
  if (.not. scbox%qp_computed) then
    call construct_RedBasis_qp(scbox%mom1, &
                               scbox%mom2, &
                               Redbasis_qp)
    call construct_p3scalars_qp(scbox%mom3, &
                                Redbasis_qp,scalars,gd2,gd3)

    RedSet_qp = redset4_qp(redbasis=Redbasis_qp,&
                           p3scalars=scalars,   &
                           perm=scbox%perm,     &
                           mom3=scbox%mom3,     &
                           gd2=gd2,gd3=gd3)
    call compute_scalar_box_qp(scbox%mom_ind, &
                               scbox%masses2, &
                               RedSet_qp,     &
                               scbox_qp)
    scbox%qp_computed = .TRUE.
    scbox%poles_qp = scbox_qp%poles
    scbox%onshell_cuts_qp = scbox_qp%onshell_cuts
  end if

end subroutine upgrade_scalar_box
#endif

! --- PENTAGONS ---
! ********************************************************************
subroutine reduction_5points_ofr(rank, momenta, masses, Gtensor, M2, &
 scboxes, all_scboxes)
! ********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_kinematics_/**/REALKIND, only: get_LC_5,get_mass2
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: G_TensorShift_otf_qp=>G_TensorShift_otf
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp, req_qp_cmp, upgrade_qp_hcl, downgrade_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_coefficient_qp=>box_coefficient
  use ol_kinematics_/**/QREALKIND, only: get_LC_5_qp=>get_LC_5,get_mass2_qp=>get_mass2
#endif
  implicit none
  integer,   intent(in) :: rank
  type(hcl), intent(inout)  :: Gtensor
  integer,   intent(in)  :: momenta(:), masses(:)
  type(met), intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)
  real(REALKIND) :: error,error_step,max_error_opp

  complex(REALKIND) :: offshell_p(1:5,1), offshell_m2(1)
  complex(REALKIND) :: Gsum(size(Gtensor%cmp))
  complex(REALKIND) :: box_coeff, box, scalar_box(0:2), q0p_q0m(2,5)
#ifdef PRECISION_dp
  complex(QREALKIND) :: offshell_p_qp(1:5,1), offshell_m2_qp(1)
  complex(QREALKIND) :: Gsum_qp(size(Gtensor%cmp))
  complex(QREALKIND) :: box_coeff_qp, box_qp, scalar_box_qp(0:2), q0p_q0m_qp(2,5)
#endif

  integer :: k,momshift

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 pentagon")

  max_error_opp = 0
  do k = 1, 5
#ifdef PRECISION_dp
    error_step = all_scboxes(scboxes(k))%box_error + all_scboxes(scboxes(k))%cut_error
    max_error_opp = max(error_step, max_error_opp)
    error = Gtensor%error + error_step

    if (hp_switch .eq. 1 .and. error_step > hp_step_thres .and. error > hp_err_thres) call upgrade_qp_hcl(Gtensor)
    if (req_dp_cmp(Gtensor)) then
#endif
    q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
    scalar_box = all_scboxes(scboxes(k))%poles
    Gsum = Gtensor%cmp
    if(k == 5) then
      momshift = -momenta(1)
      call G_TensorShift_otf(Gsum,get_LC_5(momshift))
      offshell_p(1:5,1) = get_LC_5(momshift)
      offshell_m2(1)  = get_mass2(masses(1))
    else
      offshell_p(1:5,1) = get_LC_5(sum(momenta(:k)))
      offshell_m2(1)  = get_mass2(masses(k+1))
    end if

    call box_coefficient(offshell_p(:,:),offshell_m2(:),q0p_q0m,Gsum,box_coeff)
    box = box_coeff*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2%cmp = M2%cmp + box
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gtensor)) then
      call upgrade_scalar_box(all_scboxes(scboxes(k)))
      q0p_q0m_qp = all_scboxes(scboxes(k))%onshell_cuts_qp
      scalar_box_qp = all_scboxes(scboxes(k))%poles_qp

      Gsum_qp = Gtensor%cmp_qp
      if (k == 5) then
        momshift = -momenta(1)
        call G_TensorShift_otf_qp(Gsum_qp,get_LC_5_qp(momshift))
        offshell_p_qp(1:5,1) = get_LC_5_qp(momshift)
        offshell_m2_qp(1)  = get_mass2_qp(masses(1))
      else
        offshell_p_qp(1:5,1) = get_LC_5_qp(sum(momenta(:k)))
        offshell_m2_qp(1)  = get_mass2_qp(masses(k+1))
      end if

      call box_coefficient_qp(offshell_p_qp(:,:),offshell_m2_qp(:),q0p_q0m_qp,Gsum_qp,box_coeff_qp)
      box_qp = box_coeff_qp*(scalar_box_qp(0) + scalar_box_qp(1)*de1_IR + scalar_box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + box_qp
    end if
#endif
  end do
  Gtensor%error = Gtensor%error + max_error_opp

end subroutine reduction_5points_ofr


! --- HEXAGONS ---
! ********************************************************************
subroutine reduction_6points_ofr(rank, momenta, masses, Gtensor, M2, &
 scboxes, all_scboxes)
! ********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_kinematics_/**/REALKIND, only: get_LC_4,get_LC_5,get_mass2
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: G_TensorShift_otf_qp=>G_TensorShift_otf
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp, req_qp_cmp, upgrade_qp_hcl, downgrade_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_coefficient_qp=>box_coefficient
  use ol_kinematics_/**/QREALKIND, only: get_LC_4_qp=>get_LC_4,get_LC_5_qp=>get_LC_5,get_mass2_qp=>get_mass2
#endif
  implicit none
  integer,   intent(in) :: rank
  type(hcl), intent(inout)  :: Gtensor
  integer,   intent(in)  :: momenta(:), masses(:)
  type(met), intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)
  real(REALKIND) :: error,error_step,max_error_opp

  complex(REALKIND) :: offshell_p(1:5,2), offshell_m2(2)
  complex(REALKIND) :: Gsum(size(Gtensor%cmp))
  complex(REALKIND) :: box_coeff, box, scalar_box(0:2), q0p_q0m(2,5)
#ifdef PRECISION_dp
  complex(QREALKIND) :: offshell_p_qp(1:5,2), offshell_m2_qp(2)
  complex(QREALKIND) :: Gsum_qp(size(Gtensor%cmp))
  complex(QREALKIND) :: box_coeff_qp, box_qp, scalar_box_qp(0:2), q0p_q0m_qp(2,5)
#endif
  integer :: i,j,k,momshift

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 hexagon")

  max_error_opp = 0
  k = 0
  do i = 1, 5
    do j = i, 5
      k = k + 1
#ifdef PRECISION_dp
    error_step = all_scboxes(scboxes(k))%box_error + all_scboxes(scboxes(k))%cut_error
    max_error_opp = max(error_step, max_error_opp)
    error = Gtensor%error + error_step
    if (hp_switch .eq. 1 .and. error_step > hp_step_thres .and. error > hp_err_thres) call upgrade_qp_hcl(Gtensor)
    if (req_dp_cmp(Gtensor)) then
#endif
    q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
    scalar_box = all_scboxes(scboxes(k))%poles
    Gsum = Gtensor%cmp
    if(j == 5) then
      if (i == 1) then
        momshift = -sum(momenta(:2))
      else
        momshift = -sum(momenta(:1))
      end if
      call G_TensorShift_otf(Gsum,get_LC_4(momshift))
      offshell_p(1:5,1) = get_LC_5(momshift)
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i)) + momshift)

      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      ! in the case when D0 is pinched
      offshell_m2(1) = get_mass2(masses(1))
      offshell_m2(2) = get_mass2(masses(i+1))

    else

      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      offshell_p(1:5,1) = get_LC_5(sum(momenta(:i)))
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:j+1)))

      offshell_m2(1) = get_mass2(masses(i+1))
      offshell_m2(2) = get_mass2(masses(j+2))

    end if

    call box_coefficient(offshell_p(:,:),offshell_m2(:),q0p_q0m,Gsum,box_coeff)

    box = box_coeff*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2%cmp = M2%cmp + box
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gtensor)) then
      Gsum_qp = Gtensor%cmp_qp
      call upgrade_scalar_box(all_scboxes(scboxes(k)))
      q0p_q0m_qp = all_scboxes(scboxes(k))%onshell_cuts_qp
      scalar_box_qp = all_scboxes(scboxes(k))%poles_qp
      if(j == 5) then
        if (i == 1) then
          momshift = -sum(momenta(:2))
        else
          momshift = -sum(momenta(:1))
        end if
        call G_TensorShift_otf_qp(Gsum_qp,get_LC_4_qp(momshift))
        offshell_p_qp(1:5,1) = get_LC_5_qp(momshift)
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i)) + momshift)

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        ! in the case when D0 is pinched
        offshell_m2_qp(1) = get_mass2_qp(masses(1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i+1))

      else

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        offshell_p_qp(1:5,1) = get_LC_5_qp(sum(momenta(:i)))
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:j+1)))

        offshell_m2_qp(1) = get_mass2_qp(masses(i+1))
        offshell_m2_qp(2) = get_mass2_qp(masses(j+2))

      end if

      call box_coefficient_qp(offshell_p_qp(:,:),offshell_m2_qp(:),q0p_q0m_qp,Gsum_qp,box_coeff_qp)
      box_qp = box_coeff_qp*(scalar_box_qp(0) + scalar_box_qp(1)*de1_IR + scalar_box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + box_qp
    end if
#endif

    end do
  end do
  Gtensor%error = Gtensor%error + max_error_opp

end subroutine reduction_6points_ofr

! --- HEPTAGONS ---
! ********************************************************************
subroutine reduction_7points_ofr(rank, momenta, masses, Gtensor, M2, &
 scboxes, all_scboxes)
! ********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_kinematics_/**/REALKIND, only: get_LC_4,get_LC_5,get_mass2
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: G_TensorShift_otf_qp=>G_TensorShift_otf
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp,req_qp_cmp, &
                                           upgrade_qp_hcl,downgrade_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_coefficient_qp=>box_coefficient
  use ol_kinematics_/**/QREALKIND, only: get_LC_4_qp=>get_LC_4, &
                                         get_LC_5_qp=>get_LC_5, &
                                         get_mass2_qp=>get_mass2
#endif
  implicit none
  integer,   intent(in) :: rank
  type(hcl), intent(inout)  :: Gtensor
  integer,   intent(in)  :: momenta(:), masses(:)
  type(met), intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)
  real(REALKIND) :: error,error_step,max_error_opp

  complex(REALKIND) :: offshell_p(1:5,3), offshell_m2(3)
  complex(REALKIND) :: Gsum(size(Gtensor%cmp))
  complex(REALKIND) :: box_coeff, box, scalar_box(0:2), q0p_q0m(2,5)
#ifdef PRECISION_dp
  complex(QREALKIND) :: offshell_p_qp(1:5,3), offshell_m2_qp(3)
  complex(QREALKIND) :: Gsum_qp(size(Gtensor%cmp))
  complex(QREALKIND) :: box_coeff_qp, box_qp, scalar_box_qp(0:2), q0p_q0m_qp(2,5)
#endif
  integer :: i1,i2,i3,k,momshift

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  max_error_opp = 0
  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
    do i3 = i2, 5
      k = k + 1
#ifdef PRECISION_dp
    error_step = all_scboxes(scboxes(k))%box_error + all_scboxes(scboxes(k))%cut_error
    max_error_opp = max(error_step, max_error_opp)
    error = Gtensor%error + error_step
    if (hp_switch .eq. 1 .and. error_step > hp_step_thres .and. error > hp_err_thres) call upgrade_qp_hcl(Gtensor)
    if (req_dp_cmp(Gtensor)) then
#endif
    q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
    scalar_box = all_scboxes(scboxes(k))%poles
    Gsum = Gtensor%cmp
    if(i3 == 5) then
      if (i1 == 1) then
        if (i2 == 1) then
          momshift = -sum(momenta(:3))
        else
          momshift = -sum(momenta(:2))
        end if
      else
        momshift = -sum(momenta(:1))
      end if
      call G_TensorShift_otf(Gsum,get_LC_4(momshift))
      offshell_p(1:5,1) = get_LC_5(momshift)
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i1)) + momshift)
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i2+1)) + momshift)

      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      ! in the case when D0 is pinched
      offshell_m2(1) = get_mass2(masses(1))
      offshell_m2(2) = get_mass2(masses(i1+1))
      offshell_m2(3) = get_mass2(masses(i2+2))
    else
      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      offshell_p(1:5,1) = get_LC_5(sum(momenta(:i1)))
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i2+1)))
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i3+2)))

      offshell_m2(1) = get_mass2(masses(i1+1))
      offshell_m2(2) = get_mass2(masses(i2+2))
      offshell_m2(3) = get_mass2(masses(i3+3))
    end if

    call box_coefficient(offshell_p(:,:),offshell_m2(:),q0p_q0m,Gsum,box_coeff)
    box = box_coeff*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2%cmp = M2%cmp + box
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gtensor)) then
      Gsum_qp = Gtensor%cmp_qp
      call upgrade_scalar_box(all_scboxes(scboxes(k)))
      q0p_q0m_qp = all_scboxes(scboxes(k))%onshell_cuts_qp
      scalar_box_qp = all_scboxes(scboxes(k))%poles_qp
      if(i3 == 5) then
        if (i1 == 1) then
          if (i2 == 1) then
            momshift = -sum(momenta(:3))
          else
            momshift = -sum(momenta(:2))
          end if
        else
          momshift = -sum(momenta(:1))
        end if
        call G_TensorShift_otf_qp(Gsum_qp,get_LC_4_qp(momshift))
        offshell_p_qp(1:5,1) = get_LC_5_qp(momshift)
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i1)) + momshift)
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i2+1)) + momshift)

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        ! in the case when D0 is pinched
        offshell_m2_qp(1) = get_mass2_qp(masses(1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(3) = get_mass2_qp(masses(i2+2))
      else
        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        offshell_p_qp(1:5,1) = get_LC_5_qp(sum(momenta(:i1)))
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i2+1)))
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i3+2)))

        offshell_m2_qp(1) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i2+2))
        offshell_m2_qp(3) = get_mass2_qp(masses(i3+3))
      end if

      call box_coefficient_qp(offshell_p_qp(:,:),offshell_m2_qp(:),q0p_q0m_qp,Gsum_qp,box_coeff_qp)
      box_qp = box_coeff_qp*(scalar_box_qp(0) + scalar_box_qp(1)*de1_IR + scalar_box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + box_qp
    end if
#endif

    end do
    end do
  end do
  Gtensor%error = Gtensor%error + max_error_opp

end subroutine reduction_7points_ofr


! --- OCTAGON ---
! ********************************************************************
subroutine reduction_8points_ofr(rank, momenta, masses, Gtensor, M2, &
 scboxes, all_scboxes)
! ********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_kinematics_/**/REALKIND, only: get_LC_4,get_LC_5,get_mass2
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: G_TensorShift_otf_qp=>G_TensorShift_otf
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp,req_qp_cmp, &
                                           upgrade_qp_hcl,downgrade_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_coefficient_qp=>box_coefficient
  use ol_kinematics_/**/QREALKIND, only: get_LC_4_qp=>get_LC_4, &
                                         get_LC_5_qp=>get_LC_5, &
                                         get_mass2_qp=>get_mass2
#endif
  implicit none
  integer,   intent(in) :: rank
  type(hcl), intent(inout)  :: Gtensor
  integer,   intent(in)  :: momenta(:), masses(:)
  type(met), intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)
  real(REALKIND) :: error,error_step,max_error_opp

  complex(REALKIND) :: offshell_p(1:5,4), offshell_m2(4)
  complex(REALKIND) :: Gsum(size(Gtensor%cmp))
  complex(REALKIND) :: box_coeff, box, scalar_box(0:2), q0p_q0m(2,5)
#ifdef PRECISION_dp
  complex(QREALKIND) :: offshell_p_qp(1:5,4), offshell_m2_qp(4)
  complex(QREALKIND) :: Gsum_qp(size(Gtensor%cmp))
  complex(QREALKIND) :: box_coeff_qp, box_qp, scalar_box_qp(0:2), q0p_q0m_qp(2,5)
#endif
  integer :: i1,i2,i3,i4,k,momshift

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  max_error_opp = 0
  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
    do i3 = i2, 5
    do i4 = i3, 5
      k = k + 1

#ifdef PRECISION_dp
    error_step = all_scboxes(scboxes(k))%box_error + all_scboxes(scboxes(k))%cut_error
    max_error_opp = max(error_step, max_error_opp)
    error = Gtensor%error + error_step
    if (hp_switch .eq. 1 .and. error_step > hp_step_thres .and. error > hp_err_thres) call upgrade_qp_hcl(Gtensor)
    if (req_dp_cmp(Gtensor)) then
#endif
    q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
    scalar_box = all_scboxes(scboxes(k))%poles
    Gsum = Gtensor%cmp
    if(i4 == 5) then
      if (i1 == 1) then
        if (i2 == 1) then
          if (i3 == 1) then
            momshift = -sum(momenta(:4))
          else
            momshift = -sum(momenta(:3))
          end if
        else
          momshift = -sum(momenta(:2))
        end if
      else
        momshift = -sum(momenta(:1))
      end if
      call G_TensorShift_otf(Gsum,get_LC_4(momshift))
      offshell_p(1:5,1) = get_LC_5(momshift)
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i1)) + momshift)
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i2+1)) + momshift)
      offshell_p(1:5,4) = get_LC_5(sum(momenta(:i3+2)) + momshift)

      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      ! in the case when D0 is pinched
      offshell_m2(1) = get_mass2(masses(1))
      offshell_m2(2) = get_mass2(masses(i1+1))
      offshell_m2(3) = get_mass2(masses(i2+2))
      offshell_m2(4) = get_mass2(masses(i3+3))
    else
      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      offshell_p(1:5,1) = get_LC_5(sum(momenta(:i1)))
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i2+1)))
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i3+2)))
      offshell_p(1:5,4) = get_LC_5(sum(momenta(:i4+3)))

      offshell_m2(1) = get_mass2(masses(i1+1))
      offshell_m2(2) = get_mass2(masses(i2+2))
      offshell_m2(3) = get_mass2(masses(i3+3))
      offshell_m2(4) = get_mass2(masses(i4+4))
    end if

    call box_coefficient(offshell_p(:,:),offshell_m2(:),q0p_q0m,Gsum,box_coeff)
    box = box_coeff*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2%cmp = M2%cmp + box
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gtensor)) then
      Gsum_qp = Gtensor%cmp_qp
      call upgrade_scalar_box(all_scboxes(scboxes(k)))
      q0p_q0m_qp = all_scboxes(scboxes(k))%onshell_cuts_qp
      scalar_box_qp = all_scboxes(scboxes(k))%poles_qp
      if(i4 == 5) then
        if (i1 == 1) then
          if (i2 == 1) then
            if (i3 == 1) then
              momshift = -sum(momenta(:4))
            else
              momshift = -sum(momenta(:3))
            end if
          else
            momshift = -sum(momenta(:2))
          end if
        else
          momshift = -sum(momenta(:1))
        end if
        call G_TensorShift_otf_qp(Gsum_qp,get_LC_4_qp(momshift))
        offshell_p_qp(1:5,1) = get_LC_5_qp(momshift)
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i1)) + momshift)
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i2+1)) + momshift)
        offshell_p_qp(1:5,4) = get_LC_5_qp(sum(momenta(:i3+2)) + momshift)

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        ! in the case when D0 is pinched
        offshell_m2_qp(1) = get_mass2_qp(masses(1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(3) = get_mass2_qp(masses(i2+2))
        offshell_m2_qp(4) = get_mass2_qp(masses(i3+3))
      else
        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        offshell_p_qp(1:5,1) = get_LC_5_qp(sum(momenta(:i1)))
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i2+1)))
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i3+2)))
        offshell_p_qp(1:5,4) = get_LC_5_qp(sum(momenta(:i4+3)))

        offshell_m2_qp(1) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i2+2))
        offshell_m2_qp(3) = get_mass2_qp(masses(i3+3))
        offshell_m2_qp(4) = get_mass2_qp(masses(i4+4))
      end if

      call box_coefficient_qp(offshell_p_qp(:,:),offshell_m2_qp(:),q0p_q0m_qp,Gsum_qp,box_coeff_qp)
      box_qp = box_coeff_qp*(scalar_box_qp(0) + scalar_box_qp(1)*de1_IR + scalar_box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + box_qp
    end if
#endif

    end do
    end do
    end do
  end do
  Gtensor%error = Gtensor%error + max_error_opp

end subroutine reduction_8points_ofr


! --- NONAGONS ---
! NOT TESTED
! ********************************************************************
subroutine reduction_9points_ofr(rank, momenta, masses, Gtensor, M2, &
 scboxes, all_scboxes)
! ********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox, hcl, met
  use ol_kinematics_/**/REALKIND, only: get_LC_4,get_LC_5,get_mass2
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: G_TensorShift_otf_qp=>G_TensorShift_otf
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,hp_err_thres,hp_step_thres
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp,req_qp_cmp, &
                                           upgrade_qp_hcl,downgrade_dp
  use ol_loop_reduction_/**/QREALKIND, only: box_coefficient_qp=>box_coefficient
  use ol_kinematics_/**/QREALKIND, only: get_LC_4_qp=>get_LC_4, &
                                         get_LC_5_qp=>get_LC_5, &
                                         get_mass2_qp=>get_mass2
#endif
  implicit none
  integer,   intent(in) :: rank
  type(hcl), intent(inout)  :: Gtensor
  integer,   intent(in)  :: momenta(:), masses(:)
  type(met), intent(inout) :: M2
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(inout) :: all_scboxes(:)
  real(REALKIND) :: error,error_step,max_error_opp

  complex(REALKIND) :: offshell_p(1:5,5), offshell_m2(5)
  complex(REALKIND) :: Gsum(size(Gtensor%cmp))
  complex(REALKIND) :: box_coeff, box, scalar_box(0:2), q0p_q0m(2,5)
#ifdef PRECISION_dp
  complex(QREALKIND) :: offshell_p_qp(1:5,5), offshell_m2_qp(5)
  complex(QREALKIND) :: Gsum_qp(size(Gtensor%cmp))
  complex(QREALKIND) :: box_coeff_qp, box_qp, scalar_box_qp(0:2), q0p_q0m_qp(2,5)
#endif
  integer :: i1,i2,i3,i4,i5,k,momshift

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  max_error_opp = 0
  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
    do i3 = i2, 5
    do i4 = i3, 5
    do i5 = i4, 5
      k = k + 1

#ifdef PRECISION_dp
    error_step = all_scboxes(scboxes(k))%box_error + all_scboxes(scboxes(k))%cut_error
    max_error_opp = max(error_step, max_error_opp)
    error = Gtensor%error + error_step
    if (hp_switch .eq. 1 .and. error_step > hp_step_thres .and. error > hp_err_thres) call upgrade_qp_hcl(Gtensor)
    if (req_dp_cmp(Gtensor)) then
#endif
    q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
    scalar_box = all_scboxes(scboxes(k))%poles
    Gsum = Gtensor%cmp
    if(i5 == 5) then
      if (i1 == 1) then
        if (i2 == 1) then
          if (i3 == 1) then
            if (i4 == 1) then
              momshift = -sum(momenta(:5))
            else
              momshift = -sum(momenta(:4))
            end if
          else
            momshift = -sum(momenta(:3))
          end if
        else
          momshift = -sum(momenta(:2))
        end if
      else
        momshift = -sum(momenta(:1))
      end if
      call G_TensorShift_otf(Gsum,get_LC_4(momshift))
      offshell_p(1:5,1) = get_LC_5(momshift)
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i1)) + momshift)
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i2+1)) + momshift)
      offshell_p(1:5,4) = get_LC_5(sum(momenta(:i3+2)) + momshift)
      offshell_p(1:5,5) = get_LC_5(sum(momenta(:i4+3)) + momshift)

      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      ! in the case when D0 is pinched
      offshell_m2(1) = get_mass2(masses(1))
      offshell_m2(2) = get_mass2(masses(i1+1))
      offshell_m2(3) = get_mass2(masses(i2+2))
      offshell_m2(4) = get_mass2(masses(i3+3))
      offshell_m2(5) = get_mass2(masses(i4+4))
    else
      ! Assignment of the off-shell momenta and masses for all the sub-boxes
      offshell_p(1:5,1) = get_LC_5(sum(momenta(:i1)))
      offshell_p(1:5,2) = get_LC_5(sum(momenta(:i2+1)))
      offshell_p(1:5,3) = get_LC_5(sum(momenta(:i3+2)))
      offshell_p(1:5,4) = get_LC_5(sum(momenta(:i4+3)))
      offshell_p(1:5,5) = get_LC_5(sum(momenta(:i5+4)))

      offshell_m2(1) = get_mass2(masses(i1+1))
      offshell_m2(2) = get_mass2(masses(i2+2))
      offshell_m2(3) = get_mass2(masses(i3+3))
      offshell_m2(4) = get_mass2(masses(i4+4))
      offshell_m2(5) = get_mass2(masses(i5+5))
    end if

    call box_coefficient(offshell_p(:,:),offshell_m2(:),q0p_q0m,Gsum,box_coeff)
    box = box_coeff*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2%cmp = M2%cmp + box
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gtensor)) then
      Gsum_qp = Gtensor%cmp_qp
      call upgrade_scalar_box(all_scboxes(scboxes(k)))
      q0p_q0m_qp = all_scboxes(scboxes(k))%onshell_cuts_qp
      scalar_box_qp = all_scboxes(scboxes(k))%poles_qp
      if(i5 == 5) then
        if (i1 == 1) then
          if (i2 == 1) then
            if (i3 == 1) then
              if (i4 == 1) then
                momshift = -sum(momenta(:5))
              else
                momshift = -sum(momenta(:4))
              end if
            else
              momshift = -sum(momenta(:3))
            end if
          else
            momshift = -sum(momenta(:2))
          end if
        else
          momshift = -sum(momenta(:1))
        end if
        call G_TensorShift_otf_qp(Gsum_qp,get_LC_4_qp(momshift))
        offshell_p_qp(1:5,1) = get_LC_5_qp(momshift)
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i1)) + momshift)
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i2+1)) + momshift)
        offshell_p_qp(1:5,4) = get_LC_5_qp(sum(momenta(:i3+2)) + momshift)
        offshell_p_qp(1:5,5) = get_LC_5_qp(sum(momenta(:i4+3)) + momshift)

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        ! in the case when D0 is pinched
        offshell_m2_qp(1) = get_mass2_qp(masses(1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(3) = get_mass2_qp(masses(i2+2))
        offshell_m2_qp(4) = get_mass2_qp(masses(i3+3))
        offshell_m2_qp(5) = get_mass2_qp(masses(i4+4))
      else
        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        offshell_p_qp(1:5,1) = get_LC_5_qp(sum(momenta(:i1)))
        offshell_p_qp(1:5,2) = get_LC_5_qp(sum(momenta(:i2+1)))
        offshell_p_qp(1:5,3) = get_LC_5_qp(sum(momenta(:i3+2)))
        offshell_p_qp(1:5,4) = get_LC_5_qp(sum(momenta(:i4+3)))
        offshell_p_qp(1:5,5) = get_LC_5_qp(sum(momenta(:i5+4)))

        offshell_m2_qp(1) = get_mass2_qp(masses(i1+1))
        offshell_m2_qp(2) = get_mass2_qp(masses(i2+2))
        offshell_m2_qp(3) = get_mass2_qp(masses(i3+3))
        offshell_m2_qp(4) = get_mass2_qp(masses(i4+4))
        offshell_m2_qp(5) = get_mass2_qp(masses(i5+5))
      end if

      call box_coefficient_qp(offshell_p_qp(:,:),offshell_m2_qp(:),q0p_q0m_qp,Gsum_qp,box_coeff_qp)
      box_qp = box_coeff_qp*(scalar_box_qp(0) + scalar_box_qp(1)*de1_IR + scalar_box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + box_qp
    end if
#endif

    end do
    end do
    end do
    end do
  end do
  Gtensor%error = Gtensor%error + max_error_opp

end subroutine reduction_9points_ofr



! --- PENTAGONS ---
! ********************************************************************
subroutine reduction_5points(rank, momenta, masses2, Gtensor, M2add, &
 scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer, intent(in) :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:4), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: box_p(1:4,0:3,1:5), box_m2(0:3,5)
  complex(REALKIND) :: offshell_p(1:5,1,5), offshell_m2(1,5)
  complex(REALKIND) :: Gsum(size(Gtensor))

  complex(REALKIND) :: box_coeff(5), box(5), scalar_box(0:2), q0p_q0m(2,5)
  integer :: i, j, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 pentagon")

  M2add = 0._/**/REALKIND
  scalar_box = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  ! Assignment of the off-shell momenta for all the sub-boxes
  offshell_p(1:4,1,1:4) = momenta(1:4,1:4)
  offshell_p(5,1,1:4)   = momenta(5,1:4)
  offshell_m2(1,1:4)    = masses2(1:4)

  offshell_p(1:4,1,5) = - momenta(1:4,1)
  offshell_p(5,1,5)   = momenta(5,1)
  offshell_m2(1,5)    = masses2(0)

  do i = 1, 5

    q0p_q0m = all_scboxes(scboxes(i))%onshell_cuts
    scalar_box = all_scboxes(scboxes(i))%poles
    Gsum = Gtensor

    if (all_scboxes(scboxes(i))%box_error > box_error) then
      box_error = all_scboxes(scboxes(i))%box_error
    end if

    if(i == 5) call G_TensorShift_otf(Gsum,-momenta(1:4,1))

    call box_coefficient(offshell_p(:,:,i),&
    offshell_m2(:,i),q0p_q0m,Gsum,box_coeff(i))

    box(i) = box_coeff(i)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
    M2add = M2add + box(i)
  end do

end subroutine reduction_5points


! --- HEXAGONS ---
! ********************************************************************
subroutine reduction_6points(rank, momenta, masses2, Gtensor, M2add, &
 scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_loop_handling_/**/REALKIND, only: G_TensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:5), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: offshell_p(1:5,2), offshell_m2(2)

  complex(REALKIND) :: box_coeff(size(scboxes)), box(size(scboxes)), scalar_box(0:2)
  complex(REALKIND) :: q0p_q0m(2,5), momshift(5), Gsum(size(Gtensor))
  integer :: i, j, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 hexagon")

  M2add = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  k = 0
  do i = 1, 5
    do j = i, 5
      k = k + 1

      q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
      scalar_box = all_scboxes(scboxes(k))%poles
      Gsum = Gtensor

      if (all_scboxes(scboxes(k))%box_error > box_error) then
        box_error = all_scboxes(scboxes(k))%box_error
      end if

      if(j == 5) then
        if (i == 1) then
          momshift(1:5) = momenta(1:5,2)
        else
          momshift(1:5) = momenta(1:5,1)
        end if

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        ! in the case when D0 is pinched
        offshell_p(1:4,1) = - momshift(1:4)
        offshell_p(5,1)   = momshift(5)
        offshell_p(1:4,2) = momenta(1:4,i) - momshift(1:4)
        offshell_p(5,2)   = momenta(5,i) + momshift(5) - &
                            2*cont_LC_cntrv(momenta(1:4,i),momshift(1:4))

        offshell_m2(1) = masses2(0)
        offshell_m2(2) = masses2(i)

        call G_TensorShift_otf(Gsum,-momshift(1:4))
      else

        ! Assignment of the off-shell momenta and masses for all the sub-boxes
        offshell_p(1:5,1) = momenta(1:5,i)
        offshell_p(1:5,2) = momenta(1:5,j+1)

        offshell_m2(1) = masses2(i)
        offshell_m2(2) = masses2(j+1)
      end if

      call box_coefficient(offshell_p(:,:),&
      offshell_m2(:),q0p_q0m,Gsum,box_coeff(k))

      box(k) = box_coeff(k)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
      M2add = M2add + box(k)
    end do
  end do

end subroutine reduction_6points


! --- HEPTAGONS ---
! ********************************************************************
subroutine reduction_7points(rank, momenta, masses2, Gtensor, M2add, &
scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_loop_handling_/**/REALKIND, only: G_tensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:6), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: offshell_p(1:5,3), offshell_m2(3)

  complex(REALKIND) :: box_coeff(size(scboxes)), box(size(scboxes)), scalar_box(0:2)
  complex(REALKIND) :: q0p_q0m(2,5), momshift(5), Gsum(size(Gtensor))
  integer :: i1, i2, i3, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  M2add = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
      do i3 = i2, 5
        k = k + 1

        q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
        scalar_box = all_scboxes(scboxes(k))%poles
        Gsum = Gtensor

        if (all_scboxes(scboxes(k))%box_error > box_error) then
          box_error = all_scboxes(scboxes(k))%box_error
        end if

        if(i3 == 5) then
          if(i1 == 1) then
            if(i2 == 1) then
              momshift(1:5) = momenta(1:5,3)
            else
              momshift(1:5) = momenta(1:5,2)
            end if
          else
            momshift(1:5) = momenta(1:5,1)
          end if

          ! Assignment of the off-shell momenta and masses for all the sub-boxes
          ! in the case when D0 is pinched
          offshell_p(1:4,1) = - momshift(1:4)
          offshell_p(5,1) = momshift(5)

          offshell_p(1:4,2) = momenta(1:4,i1) - momshift(1:4)
          offshell_p(5,2) = momenta(5,i1) + momshift(5) - &
                            2*cont_LC_cntrv(momenta(1:4,i1),momshift(1:4))

          offshell_p(1:4,3) = momenta(1:4,i2+1) - momshift(1:4)
          offshell_p(5,3) = momenta(5,i2+1) + momshift(5) - &
                            2*cont_LC_cntrv(momenta(1:4,i2+1),momshift(1:4))

          offshell_m2(1) = masses2(0)
          offshell_m2(2) = masses2(i1)
          offshell_m2(3) = masses2(i2+1)

          call G_TensorShift_otf(Gsum,-momshift(1:4))
        else

          ! Assignment of the off-shell momenta and masses for all the sub-boxes
          offshell_p(1:5,1) = momenta(1:5,i1)
          offshell_p(1:5,2) = momenta(1:5,i2+1)
          offshell_p(1:5,3) = momenta(1:5,i3+2)

          offshell_m2(1) = masses2(i1)
          offshell_m2(2) = masses2(i2+1)
          offshell_m2(3) = masses2(i3+2)
        end if

      call box_coefficient(offshell_p(:,:),&
      offshell_m2(:),q0p_q0m,Gsum,box_coeff(k))

      box(k) = box_coeff(k)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
      M2add = M2add + box(k)
      end do
    end do
  end do

end subroutine reduction_7points

! --- OCTAGON ---
! ********************************************************************
subroutine reduction_8points(rank, momenta, masses2, Gtensor, M2add, &
scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_loop_handling_/**/REALKIND, only: G_tensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:7), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: offshell_p(1:5,4), offshell_m2(4)

  complex(REALKIND) :: box_coeff(size(scboxes)), box(size(scboxes)), scalar_box(0:2)
  complex(REALKIND) :: q0p_q0m(2,5), momshift(5), Gsum(size(Gtensor))
  integer :: i1, i2, i3, i4, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  M2add = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
      do i3 = i2, 5
        do i4 = i3, 5
          k = k + 1

          q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
          scalar_box = all_scboxes(scboxes(k))%poles
          Gsum = Gtensor

          if (all_scboxes(scboxes(k))%box_error > box_error) then
            box_error = all_scboxes(scboxes(k))%box_error
          end if

          if(i4 == 5) then
            if(i1 == 1) then
              if(i2 == 1) then
                if(i3 == 1) then
                  momshift(1:5) = momenta(1:5,4)
                else
                  momshift(1:5) = momenta(1:5,3)
                end if
              else
                momshift(1:5) = momenta(1:5,2)
              end if
            else
              momshift(1:5) = momenta(1:5,1)
            end if

            ! Assignment of the off-shell momenta and masses for all the sub-boxes
            ! in the case when D0 is pinched
            offshell_p(1:4,1) = - momshift(1:4)
            offshell_p(5,1) = momshift(5)

            offshell_p(1:4,2) = momenta(1:4,i1) - momshift(1:4)
            offshell_p(5,2) = momenta(5,i1) + momshift(5) - &
                              2*cont_LC_cntrv(momenta(1:4,i1),momshift(1:4))

            offshell_p(1:4,3) = momenta(1:4,i2+1) - momshift(1:4)
            offshell_p(5,3) = momenta(5,i2+1) + momshift(5) - &
                              2*cont_LC_cntrv(momenta(1:4,i2+1),momshift(1:4))

            offshell_p(1:4,4) = momenta(1:4,i3+2) - momshift(1:4)
            offshell_p(5,4) = momenta(5,i3+2) + momshift(5) - &
                              2*cont_LC_cntrv(momenta(1:4,i3+2),momshift(1:4))

            offshell_m2(1) = masses2(0)
            offshell_m2(2) = masses2(i1)
            offshell_m2(3) = masses2(i2+1)
            offshell_m2(4) = masses2(i3+2) !?

            call G_TensorShift_otf(Gsum,-momshift(1:4))
          else

            ! Assignment of the off-shell momenta and masses for all the sub-boxes
            offshell_p(1:5,1) = momenta(1:5,i1)
            offshell_p(1:5,2) = momenta(1:5,i2+1)
            offshell_p(1:5,3) = momenta(1:5,i3+2)
            offshell_p(1:5,4) = momenta(1:5,i4+3)

            offshell_m2(1) = masses2(i1)
            offshell_m2(2) = masses2(i2+1)
            offshell_m2(3) = masses2(i3+2)
            offshell_m2(4) = masses2(i4+3)
          end if

        call box_coefficient(offshell_p(:,:),&
        offshell_m2(:),q0p_q0m,Gsum,box_coeff(k))

        box(k) = box_coeff(k)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
        M2add = M2add + box(k)
        end do
      end do
    end do
  end do

end subroutine reduction_8points

! --- NONAGONS ---
! TODO:  <27-03-19, J.-N. Lang> !
! NOT TESTED
! ********************************************************************
subroutine reduction_9points(rank, momenta, masses2, Gtensor, M2add, &
scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_loop_handling_/**/REALKIND, only: G_tensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:8), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: offshell_p(1:5,5), offshell_m2(5)

  complex(REALKIND) :: box_coeff(size(scboxes)), box(size(scboxes)), scalar_box(0:2)
  complex(REALKIND) :: q0p_q0m(2,5), momshift(5), Gsum(size(Gtensor))
  integer :: i1, i2, i3, i4, i5, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  M2add = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
      do i3 = i2, 5
        do i4 = i3, 5
          do i5 = i4, 5
            k = k + 1

            q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
            scalar_box = all_scboxes(scboxes(k))%poles
            Gsum = Gtensor

            if (all_scboxes(scboxes(k))%box_error > box_error) then
              box_error = all_scboxes(scboxes(k))%box_error
            end if

            if(i5 == 5) then
              if(i1 == 1) then
                if(i2 == 1) then
                  if(i3 == 1) then
                    if(i4 == 1) then
                      momshift(1:5) = momenta(1:5,5)
                    else
                      momshift(1:5) = momenta(1:5,4)
                    end if
                  else
                    momshift(1:5) = momenta(1:5,3)
                  end if
                else
                  momshift(1:5) = momenta(1:5,2)
                end if
              else
                momshift(1:5) = momenta(1:5,1)
              end if

              ! Assignment of the off-shell momenta and masses for all the sub-boxes
              ! in the case when D0 is pinched
              offshell_p(1:4,1) = - momshift(1:4)
              offshell_p(5,1) = momshift(5)

              offshell_p(1:4,2) = momenta(1:4,i1) - momshift(1:4)
              offshell_p(5,2) = momenta(5,i1) + momshift(5) - &
                                2*cont_LC_cntrv(momenta(1:4,i1),momshift(1:4))

              offshell_p(1:4,3) = momenta(1:4,i2+1) - momshift(1:4)
              offshell_p(5,3) = momenta(5,i2+1) + momshift(5) - &
                                2*cont_LC_cntrv(momenta(1:4,i2+1),momshift(1:4))

              offshell_p(1:4,4) = momenta(1:4,i3+2) - momshift(1:4)
              offshell_p(5,4) = momenta(5,i3+2) + momshift(5) - &
                                2*cont_LC_cntrv(momenta(1:4,i3+2),momshift(1:4))

              offshell_p(1:4,5) = momenta(1:4,i4+3) - momshift(1:4)
              offshell_p(5,5) = momenta(5,i4+3) + momshift(5) - &
                                2*cont_LC_cntrv(momenta(1:4,i4+3),momshift(1:4))

              offshell_m2(1) = masses2(0)
              offshell_m2(2) = masses2(i1)
              offshell_m2(3) = masses2(i2+1)
              offshell_m2(4) = masses2(i3+2)
              offshell_m2(5) = masses2(i4+3)

              call G_TensorShift_otf(Gsum,-momshift(1:4))
            else

              ! Assignment of the off-shell momenta and masses for all the sub-boxes
              offshell_p(1:5,1) = momenta(1:5,i1)
              offshell_p(1:5,2) = momenta(1:5,i2+1)
              offshell_p(1:5,3) = momenta(1:5,i3+2)
              offshell_p(1:5,4) = momenta(1:5,i4+3)
              offshell_p(1:5,5) = momenta(1:5,i5+4)

              offshell_m2(1) = masses2(i1)
              offshell_m2(2) = masses2(i2+1)
              offshell_m2(3) = masses2(i3+2)
              offshell_m2(4) = masses2(i4+3)
              offshell_m2(5) = masses2(i5+4)
            end if

          call box_coefficient(offshell_p(:,:),&
          offshell_m2(:),q0p_q0m,Gsum,box_coeff(k))

          box(k) = box_coeff(k)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
          M2add = M2add + box(k)
          end do
        end do
      end do
    end do
  end do

end subroutine reduction_9points

! --- DECAGON ---
! TODO:  <27-03-19, J.-N. Lang> !
! NOT TESTED
! ********************************************************************
subroutine reduction_10points(rank, momenta, masses2, Gtensor, M2add, &
scboxes, all_scboxes, box_error)
! ********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_kinematics_/**/REALKIND, only: cont_LC_cntrv
  use ol_loop_handling_/**/REALKIND, only: G_tensorShift_otf
  use ol_data_types_/**/REALKIND, only: scalarbox
  use ol_loop_parameters_decl_/**/DREALKIND, only: de1_IR, de2_i_IR
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(0:9), Gtensor(:)
  complex(REALKIND), intent(out) :: M2add
  integer,         intent(in) :: scboxes(:)
  type(scalarbox), intent(in) :: all_scboxes(:)
  real(REALKIND), intent(out) :: box_error

  complex(REALKIND) :: offshell_p(1:5,6), offshell_m2(6)

  complex(REALKIND) :: box_coeff(size(scboxes)), box(size(scboxes)), scalar_box(0:2)
  complex(REALKIND) :: q0p_q0m(2,5), momshift(5), Gsum(size(Gtensor))
  integer :: i1, i2, i3, i4, i5, i6, k

  if(rank > 1) call ol_error("TI_reduction: reduction of a rank > 1 heptagon")

  M2add = 0._/**/REALKIND
  box_error = 0._/**/REALKIND

  k = 0
  do i1 = 1, 5
    do i2 = i1, 5
      do i3 = i2, 5
        do i4 = i3, 5
          do i5 = i4, 5
            do i6 = i5, 5
              k = k + 1

              q0p_q0m = all_scboxes(scboxes(k))%onshell_cuts
              scalar_box = all_scboxes(scboxes(k))%poles
              Gsum = Gtensor

              if (all_scboxes(scboxes(k))%box_error > box_error) then
                box_error = all_scboxes(scboxes(k))%box_error
              end if

              if(i6 == 5) then
                if(i1 == 1) then
                  if(i2 == 1) then
                    if(i3 == 1) then
                      if(i4 == 1) then
                        if(i5 == 1) then
                          momshift(1:5) = momenta(1:5,6)
                        else
                          momshift(1:5) = momenta(1:5,5)
                        end if
                      else
                        momshift(1:5) = momenta(1:5,4)
                      end if
                    else
                      momshift(1:5) = momenta(1:5,3)
                    end if
                  else
                    momshift(1:5) = momenta(1:5,2)
                  end if
                else
                  momshift(1:5) = momenta(1:5,1)
                end if

                ! Assignment of the off-shell momenta and masses for all the sub-boxes
                ! in the case when D0 is pinched
                offshell_p(1:4,1) = - momshift(1:4)
                offshell_p(5,1) = momshift(5)

                offshell_p(1:4,2) = momenta(1:4,i1) - momshift(1:4)
                offshell_p(5,2) = momenta(5,i1) + momshift(5) - &
                                  2*cont_LC_cntrv(momenta(1:4,i1),momshift(1:4))

                offshell_p(1:4,3) = momenta(1:4,i2+1) - momshift(1:4)
                offshell_p(5,3) = momenta(5,i2+1) + momshift(5) - &
                                  2*cont_LC_cntrv(momenta(1:4,i2+1),momshift(1:4))

                offshell_p(1:4,4) = momenta(1:4,i3+2) - momshift(1:4)
                offshell_p(5,4) = momenta(5,i3+2) + momshift(5) - &
                                  2*cont_LC_cntrv(momenta(1:4,i3+2),momshift(1:4))

                offshell_p(1:4,5) = momenta(1:4,i4+3) - momshift(1:4)
                offshell_p(5,5) = momenta(5,i4+3) + momshift(5) - &
                                  2*cont_LC_cntrv(momenta(1:4,i4+3),momshift(1:4))

                offshell_p(1:4,6) = momenta(1:4,i5+4) - momshift(1:4)
                offshell_p(5,6) = momenta(5,i5+4) + momshift(5) - &
                                  2*cont_LC_cntrv(momenta(1:4,i5+4),momshift(1:4))

                offshell_m2(1) = masses2(0)
                offshell_m2(2) = masses2(i1)
                offshell_m2(3) = masses2(i2+1)
                offshell_m2(4) = masses2(i3+2)
                offshell_m2(5) = masses2(i4+3)
                offshell_m2(6) = masses2(i5+4)

                call G_TensorShift_otf(Gsum,-momshift(1:4))
              else

                ! Assignment of the off-shell momenta and masses for all the sub-boxes
                offshell_p(1:5,1) = momenta(1:5,i1)
                offshell_p(1:5,2) = momenta(1:5,i2+1)
                offshell_p(1:5,3) = momenta(1:5,i3+2)
                offshell_p(1:5,4) = momenta(1:5,i4+3)
                offshell_p(1:5,5) = momenta(1:5,i5+4)
                offshell_p(1:5,6) = momenta(1:5,i6+5)

                offshell_m2(1) = masses2(i1)
                offshell_m2(2) = masses2(i2+1)
                offshell_m2(3) = masses2(i3+2)
                offshell_m2(4) = masses2(i4+3)
                offshell_m2(5) = masses2(i5+4)
                offshell_m2(6) = masses2(i6+5)
              end if

            call box_coefficient(offshell_p(:,:),&
            offshell_m2(:),q0p_q0m,Gsum,box_coeff(k))

            box(k) = box_coeff(k)*(scalar_box(0) + scalar_box(1)*de1_IR + scalar_box(2)*de2_i_IR)
            M2add = M2add + box(k)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine reduction_10points

end module ol_loop_reduction_/**/REALKIND
