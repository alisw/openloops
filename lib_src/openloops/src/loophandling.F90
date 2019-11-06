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


module ol_loop_handling_/**/REALKIND
  use ol_data_types_/**/REALKIND, only: hol, hcl
  use ol_parameters_decl_/**/DREALKIND, only: hp_switch,        &
                                              hybrid_zero_mode, &
                                              hybrid_dp_mode,   &
                                              hybrid_qp_mode,   &
                                              hybrid_dp_qp_mode,&
                                              hp_alloc_mode
  implicit none

#ifdef PRECISION_dp
  interface req_dp_cmp
    module procedure req_dp_cmp_hcl, req_dp_cmp_hol
  end interface

  interface req_qp_cmp
    module procedure req_qp_cmp_hcl, req_qp_cmp_hol
  end interface

  interface upgrade_qp
    module procedure upgrade_qp_hcl, upgrade_qp_hol
  end interface

  interface downgrade_dp
    module procedure downgrade_dp_hcl, downgrade_dp_hol
  end interface

#endif

  contains


!!!!!!!!!!!!!!!!!!!!!!!!!
!  Hybrid mode helpers  !
!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef PRECISION_dp
  subroutine hol_dealloc_hybrid(ol_coeff)
    use ol_data_types_/**/REALKIND, only: hol
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    type(hol), intent(inout) :: ol_coeff

    if (hp_switch .eq. 1) then
      if (hp_alloc_mode .eq. 0) then
        ol_coeff%j_qp(:,:,:,:) = ZERO
      else if (hp_alloc_mode .eq. 2) then
        if (allocated(ol_coeff%j_qp)) deallocate(ol_coeff%j_qp)
      end if
    end if

  end subroutine hol_dealloc_hybrid

  subroutine hcl_dealloc_hybrid(ol_coeff)
    use ol_data_types_/**/REALKIND, only: hol
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    type(hcl), intent(inout) :: ol_coeff

    if (hp_switch .eq. 1) then
      if (hp_alloc_mode .eq. 0) then
        ol_coeff%cmp_qp(:) = ZERO
      else if (hp_alloc_mode .eq. 2) then
        if (allocated(ol_coeff%cmp_qp)) deallocate(ol_coeff%cmp_qp)
      end if
    end if

  end subroutine hcl_dealloc_hybrid

  function req_dp_cmp_hol(Gin)
    type(hol), intent(in) :: Gin
    logical :: req_dp_cmp_hol

    if (iand(Gin%mode, hybrid_dp_mode) .ne. 0) then
      req_dp_cmp_hol = .true.
    else
      req_dp_cmp_hol = .false.
    end if
  end function req_dp_cmp_hol

  function req_dp_cmp_hcl(Gin)
    type(hcl), intent(in) :: Gin
    logical :: req_dp_cmp_hcl

    if (iand(Gin%mode, hybrid_dp_mode) .ne. 0) then
      req_dp_cmp_hcl = .true.
    else
      req_dp_cmp_hcl = .false.
    end if
  end function req_dp_cmp_hcl

  function req_qp_cmp_hol(Gin)
    type(hol), intent(in) :: Gin
    logical :: req_qp_cmp_hol

    if (Gin%mode .gt. hybrid_dp_mode .and. (hp_switch .eq. 1)) then
      req_qp_cmp_hol = .true.
    else
      req_qp_cmp_hol = .false.
    end if
  end function req_qp_cmp_hol

  function req_qp_cmp_hcl(Gin)
    type(hcl), intent(in) :: Gin
    logical :: req_qp_cmp_hcl

    if (Gin%mode .gt. hybrid_dp_mode .and. (hp_switch .eq. 1)) then
      req_qp_cmp_hcl = .true.
    else
      req_qp_cmp_hcl = .false.
    end if
  end function req_qp_cmp_hcl


  subroutine upgrade_qp_hcl(Gin)
    use KIND_TYPES, only:  QREALKIND
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    use ol_external_decl_/**/DREALKIND, only: init_qp
    use ol_kinematics_/**/DREALKIND, only: init_qp_kinematics
    use ol_parameters_decl_/**/REALKIND, only: hp_qp_kinematics_init_mode
    type(hcl), intent(inout) :: Gin

    if (hp_qp_kinematics_init_mode .gt. 0 .and. .not. init_qp) call init_qp_kinematics

    if (Gin%mode .eq. hybrid_dp_mode) then
      Gin%mode = hybrid_qp_mode
      if (hp_alloc_mode .gt. 1) call hcl_alloc_hybrid(Gin)
      Gin%cmp_qp(:) = cmplx(Gin%cmp(:),kind=QREALKIND)
      Gin%cmp(:) = ZERO
    else if (Gin%mode .eq. hybrid_dp_qp_mode) then
      Gin%mode = hybrid_qp_mode
      Gin%cmp_qp(:) = Gin%cmp_qp(:) + cmplx(Gin%cmp(:),kind=QREALKIND)
      Gin%cmp(:) = ZERO
    end if

  end subroutine upgrade_qp_hcl

  subroutine hcl_alloc_hybrid(ol_coeff)
    use KIND_TYPES, only: REALKIND
    use ol_data_types_/**/REALKIND, only: hcl

    type(hcl), intent(inout) :: ol_coeff
    integer :: rank

    if (.not. allocated(ol_coeff%cmp_qp)) then
      rank = size(ol_coeff%cmp)
      allocate(ol_coeff%cmp_qp(rank))
    end if
    ol_coeff%cmp_qp = 0
  end subroutine hcl_alloc_hybrid

  subroutine upgrade_qp_hol(Gin)
    use KIND_TYPES, only:  QREALKIND
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    use ol_external_decl_/**/DREALKIND, only: init_qp
    use ol_kinematics_/**/DREALKIND, only: init_qp_kinematics
    use ol_parameters_decl_/**/REALKIND, only: hp_qp_kinematics_init_mode
    type(hol), intent(inout) :: Gin

    if (hp_qp_kinematics_init_mode .gt. 0 .and. .not. init_qp) call init_qp_kinematics

    ! pure dp mode -> pure qp
    if (Gin%mode .eq. hybrid_dp_mode) then
      Gin%mode = 2
      if (hp_alloc_mode .gt. 1) call hol_alloc_hybrid(Gin)
      Gin%j_qp(:,:,:,:) = cmplx(Gin%j(:,:,:,:),kind=QREALKIND)
      Gin%j(:,:,:,:) = ZERO
    ! mixed dp.qp mode -> pure qp
    else if (Gin%mode .eq. hybrid_dp_qp_mode) then
      Gin%mode = 2
      Gin%j_qp(:,:,:,:) = Gin%j_qp(:,:,:,:) + cmplx(Gin%j(:,:,:,:),kind=QREALKIND)
      Gin%j(:,:,:,:) = ZERO
    end if

  end subroutine upgrade_qp_hol

  subroutine hol_alloc_hybrid(ol_coeff)
    use KIND_TYPES, only: REALKIND
    use ol_data_types_/**/REALKIND, only: hol

    type(hol), intent(inout) :: ol_coeff
    integer :: alpha,rank,beta,hel_states

    if (.not. allocated(ol_coeff%j_qp)) then
      hel_states = size(ol_coeff%hf)
      alpha = size(ol_coeff%j,1)
      rank = size(ol_coeff%j,2)
      beta = size(ol_coeff%j,3)
      allocate(ol_coeff%j_qp(alpha, rank, beta, hel_states))
    end if
    ol_coeff%j_qp = 0

  end subroutine hol_alloc_hybrid

  subroutine downgrade_dp_hcl(Gin)
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    type(hcl), intent(inout) :: Gin

    if (Gin%mode .eq. hybrid_qp_mode) then
      Gin%mode = hybrid_dp_mode
      Gin%cmp(:) = cmplx(Gin%cmp_qp(:),kind=DREALKIND)
    else if (Gin%mode .eq. hybrid_dp_qp_mode) then
      Gin%mode = hybrid_dp_mode
      Gin%cmp(:) = Gin%cmp(:) + cmplx(Gin%cmp_qp(:),kind=DREALKIND)
    end if
    call hcl_dealloc_hybrid(Gin)

  end subroutine downgrade_dp_hcl

  subroutine downgrade_dp_hol(Gin)
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/REALKIND, only: ZERO
    type(hol), intent(inout) :: Gin

    ! pure qp mode -> pure dp
    if (Gin%mode .eq. hybrid_qp_mode) then
      Gin%mode = hybrid_dp_mode
      Gin%j(:,:,:,:) = cmplx(Gin%j_qp(:,:,:,:),kind=DREALKIND)
    ! mixed dp.qp mode -> pure dp
    else if (Gin%mode .eq. hybrid_dp_qp_mode) then
      Gin%mode = hybrid_dp_mode
      Gin%j(:,:,:,:) = Gin%j(:,:,:,:) + cmplx(Gin%j_qp(:,:,:,:),kind=DREALKIND)
    end if
    call hol_dealloc_hybrid(Gin)

  end subroutine downgrade_dp_hol


#endif

  function merge_mode(mode1,mode2)
    integer, intent(in) :: mode1,mode2
    integer :: merge_mode

    select case (mode1)
    case (hybrid_zero_mode)
      select case (mode2)
      case (hybrid_zero_mode)
        merge_mode = hybrid_zero_mode
      case (hybrid_dp_mode)
        merge_mode = hybrid_dp_mode
      case (hybrid_qp_mode)
        merge_mode = hybrid_qp_mode
      case (hybrid_dp_qp_mode)
        merge_mode = hybrid_dp_qp_mode
      end select

    case (hybrid_dp_mode)
      select case (mode2)
      case (hybrid_zero_mode)
        merge_mode = hybrid_dp_mode
      case (hybrid_dp_mode)
        merge_mode = hybrid_dp_mode
      case (hybrid_qp_mode,hybrid_dp_qp_mode)
        merge_mode = hybrid_dp_qp_mode
      end select

    case (hybrid_qp_mode)
      select case (mode2)
      case (hybrid_zero_mode)
        merge_mode = hybrid_qp_mode
      case (hybrid_dp_mode,hybrid_dp_qp_mode)
        merge_mode = hybrid_dp_qp_mode
      case (hybrid_qp_mode)
        merge_mode = hybrid_qp_mode
      end select

    case (hybrid_dp_qp_mode)
      merge_mode = hybrid_dp_qp_mode
    end select

  end function merge_mode

!******************************************************************************
subroutine signflip_OLR(Ginout)
!------------------------------------------------------------------------------
! changing sign of colour factor because of inverted building
! direction for vertex of 3 coloured objects
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  type(hol), intent(inout) :: Ginout

  Ginout%j = -Ginout%j
#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. Ginout%mode .gt. hybrid_dp_mode) then
    Ginout%j_qp = -Ginout%j_qp
  end if
#endif

end subroutine signflip_OLR

!!*****************************************************************************
!! Functions for transposing an open loop (change of dressing direction)     !!
!!*****************************************************************************

!******************************************************************************
subroutine HGT_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: al ,be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do al=1,4
      do be=1,4
        Gout(al,r1:r2,be,h)=Gin%j(be,r1:r2,al,h)
      end do
    end do
  end do

  Gin%j(:,r1:r2,:,:)=Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do al=1,4
        do be=1,4
          Gout_qp(al,r1:r2,be,h)=Gin%j_qp(be,r1:r2,al,h)
        end do
      end do
    end do

  Gin%j_qp(:,r1:r2,:,:)=Gout_qp

  end if
#endif

end subroutine HGT_OLR


!******************************************************************************
subroutine HGT_invQ_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! change sign because of inverted direction of q
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: al, be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do al=1,4
      do be=1,4
        Gout(al,r1:r2,be,h)=-Gin%j(be,r1:r2,al,h)
      end do
    end do
  end do

  Gin%j(:,r1:r2,:,:)=Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do al=1,4
        do be=1,4
          Gout_qp(al,r1:r2,be,h)=-Gin%j_qp(be,r1:r2,al,h)
        end do
      end do
    end do

    Gin%j_qp(:,r1:r2,:,:)=Gout_qp

  end if
#endif

end subroutine HGT_invQ_OLR


!******************************************************************************
subroutine HGT_w2_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop
! NOTE: The factor 2 stemming from the metric when raising alpha before is
!       multiplied here
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: al ,be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do al=1,4
      do be=1,4
        Gout(al,r1:r2,be,h)=Gin%j(be,r1:r2,al,h)
      end do
    end do
  end do

  !! Factor 2 from raising alpha before
  Gin%j(:,r1:r2,:,:)=Gout+Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do al=1,4
        do be=1,4
          Gout_qp(al,r1:r2,be,h)=Gin%j_qp(be,r1:r2,al,h)
        end do
      end do
    end do

    Gin%j_qp(:,r1:r2,:,:)=Gout_qp+Gout_qp

  end if
#endif

end subroutine HGT_w2_OLR


!******************************************************************************
subroutine HGT_w2_invQ_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! change sign because of inverted direction of q
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop
! NOTE: The factor 2 stemming from the metric when raising alpha before is
!       multiplied here
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: al ,be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do al=1,4
      do be=1,4
        Gout(al,r1:r2,be,h)=-Gin%j(be,r1:r2,al,h)
      end do
    end do
  end do

  !! Factor 2 from raising alpha before
  Gin%j(:,r1:r2,:,:) = Gout + Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do al=1,4
        do be=1,4
          Gout_qp(al,r1:r2,be,h)=-Gin%j_qp(be,r1:r2,al,h)
        end do
      end do
    end do

    !! Factor 2 from raising alpha before
    Gin%j_qp(:,r1:r2,:,:) = Gout_qp + Gout_qp
  end if
#endif

end subroutine HGT_w2_invQ_OLR


!******************************************************************************
subroutine HGT_raise_alpha_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop is raised to become a contravariant
!         (light-cone) "active" index contracted with vertices/props to build
!         the loop against the usual direction
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop in usual direction (dir=0)
! NOTE: A factor 2 stemming from the metric is suppressed (to be multiplied
!       later)
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: be, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do be=1,4
    Gout(2,r1:r2,be,1:phys_hel)= Gin%j(be,r1:r2,1,1:phys_hel)
    Gout(1,r1:r2,be,1:phys_hel)= Gin%j(be,r1:r2,2,1:phys_hel)
    Gout(4,r1:r2,be,1:phys_hel)=-Gin%j(be,r1:r2,3,1:phys_hel)
    Gout(3,r1:r2,be,1:phys_hel)=-Gin%j(be,r1:r2,4,1:phys_hel)
  end do

  Gin%j(:,r1:r2,:,:) = Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do be=1,4
      Gout_qp(2,r1:r2,be,1:phys_hel)= Gin%j_qp(be,r1:r2,1,1:phys_hel)
      Gout_qp(1,r1:r2,be,1:phys_hel)= Gin%j_qp(be,r1:r2,2,1:phys_hel)
      Gout_qp(4,r1:r2,be,1:phys_hel)=-Gin%j_qp(be,r1:r2,3,1:phys_hel)
      Gout_qp(3,r1:r2,be,1:phys_hel)=-Gin%j_qp(be,r1:r2,4,1:phys_hel)
    end do

    Gin%j_qp(:,r1:r2,:,:) = Gout_qp
  end if
#endif

end subroutine HGT_raise_alpha_OLR


!******************************************************************************
subroutine HGT_raise_alpha_invQ_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(beta,l,alpha) for fixed l = tensor index
! change sign because of inverted direction of q
! alpha = covariant (light-cone) "frozen" open index, is untouched till the
!         last contraction in the loop is raised to become a contravariant
!         (light-cone) "active" index contracted with vertices/props to build
!         the loop against the usual direction
! beta  = contravariant (light-cone) "active" index contracted with vertices/
!         props to build the loop in usual direction
! NOTE: A factor 2 stemming from the metric is suppressed (to be multiplied
!       later)
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  integer :: be, h, phys_hel
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do be=1,4
      Gout(2,r1:r2,be,h)= -Gin%j(be,r1:r2,1,h)
      Gout(1,r1:r2,be,h)= -Gin%j(be,r1:r2,2,h)
      Gout(4,r1:r2,be,h)=  Gin%j(be,r1:r2,3,h)
      Gout(3,r1:r2,be,h)=  Gin%j(be,r1:r2,4,h)
    end do
  end do

  Gin%j(:,r1:r2,:,:) = Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do be=1,4
        Gout_qp(2,r1:r2,be,h)= -Gin%j_qp(be,r1:r2,1,h)
        Gout_qp(1,r1:r2,be,h)= -Gin%j_qp(be,r1:r2,2,h)
        Gout_qp(4,r1:r2,be,h)=  Gin%j_qp(be,r1:r2,3,h)
        Gout_qp(3,r1:r2,be,h)=  Gin%j_qp(be,r1:r2,4,h)
      end do
    end do

    Gin%j_qp(:,r1:r2,:,:) = Gout_qp

  end if
#endif

end subroutine HGT_raise_alpha_invQ_OLR


!******************************************************************************
subroutine HGT_lower_alpha_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(alpha,l,beta) for fixed l = tensor index lowering
! alpha (contravariant in lightcone)
! NOTE: A factor 1/2 stemming from the metric is suppressed
!       (cancelling a factor 2 from raising alpha before)
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do be=1,4
      Gout(be,r1:r2,2,h)= Gin%j(1,r1:r2,be,h)
      Gout(be,r1:r2,1,h)= Gin%j(2,r1:r2,be,h)
      Gout(be,r1:r2,4,h)=-Gin%j(3,r1:r2,be,h)
      Gout(be,r1:r2,3,h)=-Gin%j(4,r1:r2,be,h)
    end do
  end do

 Gin%j(:,r1:r2,:,:) = Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do be=1,4
        Gout_qp(be,r1:r2,2,h)= Gin%j_qp(1,r1:r2,be,h)
        Gout_qp(be,r1:r2,1,h)= Gin%j_qp(2,r1:r2,be,h)
        Gout_qp(be,r1:r2,4,h)=-Gin%j_qp(3,r1:r2,be,h)
        Gout_qp(be,r1:r2,3,h)=-Gin%j_qp(4,r1:r2,be,h)
      end do
    end do

   Gin%j_qp(:,r1:r2,:,:) = Gout_qp
  end if
#endif

end subroutine HGT_lower_alpha_OLR


!******************************************************************************
subroutine HGT_lower_alpha_invQ_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(alpha,l,beta) for fixed l = tensor index lowering
! change sign because of inverted direction of q
! alpha (contravariant in lightcone)
! NOTE: A factor 1/2 stemming from the metric is suppressed
!       (cancelling a factor 2 from raising alpha before)
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do be=1,4
      Gout(be,r1:r2,2,h)= -Gin%j(1,r1:r2,be,h)
      Gout(be,r1:r2,1,h)= -Gin%j(2,r1:r2,be,h)
      Gout(be,r1:r2,4,h)=  Gin%j(3,r1:r2,be,h)
      Gout(be,r1:r2,3,h)=  Gin%j(4,r1:r2,be,h)
    end do
  end do

  Gin%j(:,r1:r2,:,:) = Gout

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do be=1,4
        Gout_qp(be,r1:r2,2,h)= -Gin%j_qp(1,r1:r2,be,h)
        Gout_qp(be,r1:r2,1,h)= -Gin%j_qp(2,r1:r2,be,h)
        Gout_qp(be,r1:r2,4,h)=  Gin%j_qp(3,r1:r2,be,h)
        Gout_qp(be,r1:r2,3,h)=  Gin%j_qp(4,r1:r2,be,h)
      end do
    end do

    Gin%j_qp(:,r1:r2,:,:) = Gout_qp

  end if
#endif

end subroutine HGT_lower_alpha_invQ_OLR


!******************************************************************************
subroutine HGT_lower_alpha_w2_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(alpha,l,beta) for fixed l = tensor index lowering
! alpha (contravariant in lightcone)
! NOTE: The factor 1/2 stemming from the metric is multiplied here
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do be=1,4
      Gout(be,r1:r2,2,h)= Gin%j(1,r1:r2,be,h)
      Gout(be,r1:r2,1,h)= Gin%j(2,r1:r2,be,h)
      Gout(be,r1:r2,4,h)=-Gin%j(3,r1:r2,be,h)
      Gout(be,r1:r2,3,h)=-Gin%j(4,r1:r2,be,h)
    end do
  end do

  Gin%j(:,r1:r2,:,:) = Gout/2

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do be=1,4
        Gout_qp(be,r1:r2,2,h)= Gin%j_qp(1,r1:r2,be,h)
        Gout_qp(be,r1:r2,1,h)= Gin%j_qp(2,r1:r2,be,h)
        Gout_qp(be,r1:r2,4,h)=-Gin%j_qp(3,r1:r2,be,h)
        Gout_qp(be,r1:r2,3,h)=-Gin%j_qp(4,r1:r2,be,h)
      end do
    end do

    Gin%j_qp(:,r1:r2,:,:) = Gout_qp/2
  end if
#endif

end subroutine HGT_lower_alpha_w2_OLR


!******************************************************************************
subroutine HGT_lower_alpha_w2_invQ_OLR(Gin,r1,r2,hel)
!------------------------------------------------------------------------------
! Transpose G(alpha,l,beta) for fixed l = tensor index lowering
! change sign because of inverted direction of q
! alpha (contravariant in lightcone)
! NOTE: The factor 1/2 stemming from the metric is multiplied here
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  implicit none
  integer, intent(in) :: r1, r2, hel
  type(hol), intent(inout) :: Gin
  complex(REALKIND) :: Gout(4,r1:r2,4,size(Gin%hf))
  integer :: be, h, phys_hel
#ifdef PRECISION_dp
  complex(QREALKIND) :: Gout_qp(4,r1:r2,4,size(Gin%hf))
#endif

  phys_hel = size(Gin%hf)
  do h=1,phys_hel
    do be=1,4
      Gout(be,r1:r2,2,h)= -Gin%j(1,r1:r2,be,h)
      Gout(be,r1:r2,1,h)= -Gin%j(2,r1:r2,be,h)
      Gout(be,r1:r2,4,h)=  Gin%j(3,r1:r2,be,h)
      Gout(be,r1:r2,3,h)=  Gin%j(4,r1:r2,be,h)
    end do
  end do
 Gin%j(:,r1:r2,:,:) = Gout/2

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    do h=1,phys_hel
      do be=1,4
        Gout_qp(be,r1:r2,2,h)= -Gin%j_qp(1,r1:r2,be,h)
        Gout_qp(be,r1:r2,1,h)= -Gin%j_qp(2,r1:r2,be,h)
        Gout_qp(be,r1:r2,4,h)=  Gin%j_qp(3,r1:r2,be,h)
        Gout_qp(be,r1:r2,3,h)=  Gin%j_qp(4,r1:r2,be,h)
      end do
    end do
   Gin%j_qp(:,r1:r2,:,:) = Gout_qp/2
  end if
#endif

end subroutine HGT_lower_alpha_w2_invQ_OLR


!!*****************************************************************************
!!            Functions for shifting the loop momentum                       !!
!!*****************************************************************************

!******************************************************************************
subroutine HG1shiftOLR(HG1in,dQid,htot)
!------------------------------------------------------------------------------
! Shift the momentum q -> q + dQ for the tensor G1(beta,l,alpha), i.e.
! the scalar part of G1 is shifted by the vector part contracted with dQ
! (all in light-cone rep)
!******************************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, QREALKIND
  use ol_kinematics_/**/REALKIND, only: get_LC_4
  use ol_data_types_/**/REALKIND, only: hol
#ifdef PRECISION_dp
  use ol_kinematics_/**/DREALKIND, only: get_LC_4_qp
#endif
  integer,   intent(in)    :: htot
  integer,   intent(in)    :: dQid
  type(hol), intent(inout) :: HG1in
  complex(REALKIND) :: dQ(4)
  complex(REALKIND) :: G1shift(4,1,4,size(HG1in%hf))
#ifdef PRECISION_dp
  complex(QREALKIND) :: G1shift_qp(4,1,4,size(HG1in%hf))
  complex(QREALKIND) :: dQ_qp(4)
#endif

  dQ = get_LC_4(dQid)

  G1shift(:,1,:,:) = HG1in%j(:,2,:,:)*dQ(1) + HG1in%j(:,3,:,:)*dQ(2) &
                   + HG1in%j(:,4,:,:)*dQ(3) + HG1in%j(:,5,:,:)*dQ(4)

  HG1in%j(:,1,:,:) = HG1in%j(:,1,:,:) + G1shift(:,1,:,:)

#ifdef PRECISION_dp
  if (req_qp_cmp(HG1in)) then
    dQ_qp = get_LC_4_qp(dQid)
    G1shift_qp(:,1,:,:) = HG1in%j_qp(:,2,:,:)*dQ_qp(1) + HG1in%j_qp(:,3,:,:)*dQ_qp(2) &
                     + HG1in%j_qp(:,4,:,:)*dQ_qp(3) + HG1in%j_qp(:,5,:,:)*dQ_qp(4)

    HG1in%j_qp(:,1,:,:) = HG1in%j_qp(:,1,:,:) + G1shift_qp(:,1,:,:)
  end if
#endif

end subroutine HG1shiftOLR


!******************************************************************************
subroutine G_TensorShift(Gin,dQid)
! -----------------------------------------------------------------------------
! It performs a loop-momentum shift, i.e. q -> q + dQ, for the tensor Gin(l).
! rank-1, rank-2 and rank-3 supported
!******************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_data_types_/**/REALKIND, only: hcl
  use ol_kinematics_/**/REALKIND, only: get_LC_4
#ifdef PRECISION_dp
  use ol_loop_handling_/**/QREALKIND, only: &
              G1tensorshiftOLR_qp => G1tensorshiftOLR, &
              G2tensorshiftOLR_qp => G2tensorshiftOLR, &
              G3tensorshiftOLR_qp => G3tensorshiftOLR
  use ol_kinematics_/**/DREALKIND, only: get_LC_4_qp
#endif
  integer,   intent(in)    :: dQid
  type(hcl), intent(inout) :: Gin

  if (size(Gin%cmp) == 5) then
    call G1tensorshiftOLR(Gin%cmp,get_LC_4(dQid))
  else if (size(Gin%cmp) == 15) then
    call G2tensorshiftOLR(Gin%cmp,get_LC_4(dQid))
  else if (size(Gin%cmp) == 35) then
    call G3tensorshiftOLR(Gin%cmp,get_LC_4(dQid))
  end if

#ifdef PRECISION_dp
  if (req_qp_cmp(Gin)) then
    if (size(Gin%cmp) == 5) then
      call G1tensorshiftOLR_qp(Gin%cmp_qp,get_LC_4_qp(dQid))
    else if (size(Gin%cmp) == 15) then
      call G2tensorshiftOLR_qp(Gin%cmp_qp,get_LC_4_qp(dQid))
    else if (size(Gin%cmp) == 35) then
      call G3tensorshiftOLR_qp(Gin%cmp_qp,get_LC_4_qp(dQid))
    end if
  end if
#endif

end subroutine G_TensorShift

!******************************************************************************
subroutine G_TensorShift_otf(Gin,dQ)
! -----------------------------------------------------------------------------
! It performs a loop-momentum shift, i.e. q -> q + dQ, for the tensor Gin(l).
! rank-1, rank-2 and rank-3 supported
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: dQ(4)
  complex(REALKIND), intent(inout) :: Gin(:)

  if (size(Gin) == 5) then
    call G1tensorshiftOLR(Gin,dQ)
  else if (size(Gin) == 15) then
    call G2tensorshiftOLR(Gin,dQ)
  else if (size(Gin) == 35) then
    call G3tensorshiftOLR(Gin,dQ)
  end if

end subroutine G_TensorShift_otf

!******************************************************************************
subroutine G1tensorshiftOLR(G1in,dQ)
! -----------------------------------------------------------------------------
! Shift the momentum q -> q + dQ for the tensor G1tensor(l), i.e. the scalar
! part of G1 is shifted by the vector part contracted with dQ (all in
! light-cone rep)
! l = tensor index
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: dQ(4)
  complex(REALKIND), intent(inout) :: G1in(5)
  complex(REALKIND) :: G1shift

  G1shift = G1in(2)*dQ(1) + G1in(3)*dQ(2) + G1in(4)*dQ(3) + G1in(5)*dQ(4)
  G1in(1) = G1in(1) + G1shift
end subroutine G1tensorshiftOLR


!******************************************************************************
subroutine G2tensorshiftOLR(Gin,dQ)
! -----------------------------------------------------------------------------
! Shift the momentum q -> q + dQ for the tensor G2tensor(l), i.e. the vector
! and the scalar part  by the tensor and vector part contracted with dQ
! (all in light-cone rep)
! l = tensor index
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: dQ(4)
  complex(REALKIND), intent(inout) :: Gin(15)
  complex(REALKIND) :: g2q0, g3q1,  g4q2,  g5q3
  complex(REALKIND) :: g6q0, g7q1,  g8q2,  g9q3
  complex(REALKIND) :: g7q0, g10q1, g11q2, g12q3
  complex(REALKIND) :: g8q0, g11q1, g13q2, g14q3
  complex(REALKIND) :: g9q0, g12q1, g14q2, g15q3
  complex(REALKIND) :: g2, g3, g4

  g2q0=Gin(2)*dQ(1)
  g3q1=Gin(3)*dQ(2)
  g4q2=Gin(4)*dQ(3)
  g5q3=Gin(5)*dQ(4)

  g6q0=Gin(6)*dQ(1)
  g7q1=Gin(7)*dQ(2)
  g8q2=Gin(8)*dQ(3)
  g9q3=Gin(9)*dQ(4)

  g7q0 =Gin(7)*dQ(1)
  g10q1=Gin(10)*dQ(2)
  g11q2=Gin(11)*dQ(3)
  g12q3=Gin(12)*dQ(4)

  g8q0 =Gin(8)*dQ(1)
  g11q1=Gin(11)*dQ(2)
  g13q2=Gin(13)*dQ(3)
  g14q3=Gin(14)*dQ(4)

  g9q0 =Gin(9)*dQ(1)
  g12q1=Gin(12)*dQ(2)
  g14q2=Gin(14)*dQ(3)
  g15q3=Gin(15)*dQ(4)

  g2= g6q0  + g7q1  + g8q2  + g9q3
  g3= g10q1 + g11q2 + g12q3
  g4= g13q2 + g14q3

  Gin(1) = Gin(1) + (g2q0 + g3q1 + g4q2 + g5q3) + dQ(1)*g2 + dQ(2)*g3 + &
           dQ(3)*g4 + dQ(4)*g15q3
  Gin(2) = Gin(2) + g2 + g6q0
  Gin(3) = Gin(3) + g3 + g10q1 + g7q0
  Gin(4) = Gin(4) + g4 + g13q2 + g11q1 + g8q0
  Gin(5) = Gin(5) + g15q3 + g15q3 + g14q2 + g12q1 + g9q0

end subroutine G2tensorshiftOLR

!******************************************************************************
subroutine G3tensorshiftOLR(Gin,dQ)
! -----------------------------------------------------------------------------
! Shift the momentum q -> q + dQ for the tensor G3tensor(l), i.e. the vector
! and the scalar part by the tensor and vector part contracted with dQ
! (all in light-cone rep)
! l = tensor index
!******************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: dQ(4)
  complex(REALKIND), intent(inout) :: Gin(35)
  complex(REALKIND) :: g16q1, g17q2, g18q3, g19q4
  complex(REALKIND) :: g17q1, g20q2, g21q3, g22q4
  complex(REALKIND) :: g18q1, g21q2, g23q3, g24q4
  complex(REALKIND) :: g19q1, g22q2, g24q3, g25q4
  complex(REALKIND) :: g20q1, g26q2, g27q3, g28q4
  complex(REALKIND) :: g21q1, g27q2, g29q3, g30q4
  complex(REALKIND) :: g22q1, g28q2, g30q3, g31q4
  complex(REALKIND) :: g23q1, g29q2, g32q3, g33q4
  complex(REALKIND) :: g24q1, g30q2, g33q3, g34q4
  complex(REALKIND) :: g25q1, g31q2, g34q3, g35q4
  complex(REALKIND) :: g6, g7, g8, g9, g10, g11, g12, g13, g14, g15
  complex(REALKIND) :: Gtmp_r2(1:15)

  g16q1 =Gin(16)*dQ(1)
  g17q2 =Gin(17)*dQ(2)
  g18q3 =Gin(18)*dQ(3)
  g19q4 =Gin(19)*dQ(4)

  g17q1 =Gin(17)*dQ(1)
  g20q2 =Gin(20)*dQ(2)
  g21q3 =Gin(21)*dQ(3)
  g22q4 =Gin(22)*dQ(4)

  g18q1 =Gin(18)*dQ(1)
  g21q2 =Gin(21)*dQ(2)
  g23q3 =Gin(23)*dQ(3)
  g24q4 =Gin(24)*dQ(4)

  g19q1 =Gin(19)*dQ(1)
  g22q2 =Gin(22)*dQ(2)
  g24q3 =Gin(24)*dQ(3)
  g25q4 =Gin(25)*dQ(4)

  g20q1 =Gin(20)*dQ(1)
  g26q2 =Gin(26)*dQ(2)
  g27q3 =Gin(27)*dQ(3)
  g28q4 =Gin(28)*dQ(4)

  g21q1 =Gin(21)*dQ(1)
  g27q2 =Gin(27)*dQ(2)
  g29q3 =Gin(29)*dQ(3)
  g30q4 =Gin(30)*dQ(4)

  g22q1 =Gin(22)*dQ(1)
  g28q2 =Gin(28)*dQ(2)
  g30q3 =Gin(30)*dQ(3)
  g31q4 =Gin(31)*dQ(4)

  g23q1 =Gin(23)*dQ(1)
  g29q2 =Gin(29)*dQ(2)
  g32q3 =Gin(32)*dQ(3)
  g33q4 =Gin(33)*dQ(4)

  g24q1 =Gin(24)*dQ(1)
  g30q2 =Gin(30)*dQ(2)
  g33q3 =Gin(33)*dQ(3)
  g34q4 =Gin(34)*dQ(4)

  g25q1 =Gin(25)*dQ(1)
  g31q2 =Gin(31)*dQ(2)
  g34q3 =Gin(34)*dQ(3)
  g35q4 =Gin(35)*dQ(4)

  g6  = g16q1 + g17q2 + g18q3 + g19q4
  g7  = g17q1 + g20q2 + g21q3 + g22q4
  g8  = g18q1 + g21q2 + g23q3 + g24q4
  g9  = g19q1 + g22q2 + g24q3 + g25q4
  g10 = g20q1 + g26q2 + g27q3 + g28q4
  g11 = g21q1 + g27q2 + g29q3 + g30q4
  g12 = g22q1 + g28q2 + g30q3 + g31q4
  g13 = g23q1 + g29q2 + g32q3 + g33q4
  g14 = g24q1 + g30q2 + g33q3 + g34q4
  g15 = g25q1 + g31q2 + g34q3 + g35q4

  Gtmp_r2(6)  = 2*g16q1       + g6
  Gtmp_r2(7)  = g17q1 + g20q2 + g7
  Gtmp_r2(8)  = g18q1 + g23q3 + g8
  Gtmp_r2(9)  = g19q1 + g25q4 + g9
  Gtmp_r2(10) = 2*g26q2       + g10
  Gtmp_r2(11) = g27q2 + g29q3 + g11
  Gtmp_r2(12) = g28q2 + g31q4 + g12
  Gtmp_r2(13) = 2*g32q3       + g13
  Gtmp_r2(14) = g33q3 + g34q4 + g14
  Gtmp_r2(15) = 2*g35q4       + g15

  Gtmp_r2(2) = g6*dQ(1) + g7*dQ(2)  + g8*dQ(3)  + g9*dQ(4)  +&
  2*g16q1*dQ(1) - (g21q3+g22q4)*dQ(2) - g24q3*dQ(4)

  Gtmp_r2(3) = g7*dQ(1) + g10*dQ(2) + g11*dQ(3) + g12*dQ(4) +&
  2*g26q2*dQ(2) - (g21q3+g22q4)*dQ(1) - g30q3*dQ(4)

  Gtmp_r2(4) = g8*dQ(1) + g11*dQ(2) + g13*dQ(3) + g14*dQ(4) +&
  2*g32q3*dQ(3) - (g21q2+g24q4)*dQ(1) - g30q2*dQ(4)

  Gtmp_r2(5) = g9*dQ(1) + g12*dQ(2) + g14*dQ(3) + g15*dQ(4) +&
  2*g35q4*dQ(4) - (g22q2+g24q3)*dQ(1) - g30q2*dQ(3)

  Gtmp_r2(1) = g6*(dQ(1)**2) + g10*(dQ(2)**2) + g13*(dQ(3)**2) + g15*(dQ(4)**2) + &
  (g24q3*dQ(1) + g30q3*dQ(2))*dQ(4) + (g21q1*dQ(3) + g22q1*dQ(4))*dQ(2)

  call G2tensorshiftOLR(Gin(1:15),dQ)

  Gin(1)    = Gin(1)    + Gtmp_r2(1)
  Gin(2:15) = Gin(2:15) + Gtmp_r2(2:15)

end subroutine G3tensorshiftOLR

end module ol_loop_handling_/**/REALKIND
