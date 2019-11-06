!******************************************************************************!
!                                                                              !
!    ttred.f90                                                                 !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

!Last Modified: January 10, 2019
!Version 0.9.3

module trred
  use triangle_expansion_dp, only: &
  get_errflag_trred=>get_errflag, &
  reset_errflag_trred=>reset_errflag

  use triangle_aux_dp, only: &
  set_muUV2_dp=>set_muUV2, &
  set_muIR2_dp=>set_muIR2, &
  set_duv_dp=>set_duv, &
  set_dir_dp=>set_dir

  use triangle_aux_qp, only: &
  set_muUV2_qp=>set_muUV2, &
  set_muIR2_qp=>set_muIR2, &
  set_duv_qp=>set_duv, &
  set_dir_qp=>set_dir

  use triangle_reduction_dp, only: &
  C_0mm_P12_DP, C_0mm_p1_DP, C_0mm_p1p1_DP, C_0mm_p1P12_DP, C_0mm_P12P12_DP, &
  C_0mm_g_DP, C_0mm_p1p1p1_DP, C_0mm_p1p1P12_DP, C_0mm_p1P12P12_DP, &
  C_0mm_P12P12P12_DP, C_0mm_gp1_DP, C_0mm_gP12_DP, &
  C_m00_P12_DP, C_m00_p1_DP, C_m00_p1p1_DP, C_m00_p1P12_DP, C_m00_P12P12_DP, &
  C_m00_g_DP, C_m00_p1p1p1_DP, C_m00_p1p1P12_DP, C_m00_p1P12P12_DP, &
  C_m00_P12P12P12_DP, C_m00_gp1_DP, C_m00_gP12_DP, &
  C_mmm_P12_DP, C_mmm_p1_DP, C_mmm_p1p1_DP, C_mmm_p1P12_DP, C_mmm_P12P12_DP, &
  C_mmm_g_DP, C_mmm_p1p1p1_DP, C_mmm_p1p1P12_DP, C_mmm_p1P12P12_DP, &
  C_mmm_P12P12P12_DP, C_mmm_gp1_DP, C_mmm_gP12_DP, &
  C_000_P12_DP, C_000_p1_DP, C_000_p1p1_DP, C_000_p1P12_DP, C_000_P12P12_DP, &
  C_000_g_DP, C_000_p1p1p1_DP, C_000_p1p1P12_DP, C_000_p1P12P12_DP, &
  C_000_P12P12P12_DP, C_000_gp1_DP, C_000_gP12_DP, &
  C_m0m1m1_p1_DP, C_m0m1m1_P12_DP, &
  C_m0m1m1_p1p1_DP, C_m0m1m1_p1P12_DP, C_m0m1m1_P12P12_DP, C_m0m1m1_g_DP, &
  C_m0m1m1_p1p1p1_DP, C_m0m1m1_p1p1P12_DP, C_m0m1m1_p1P12P12_DP, &
  C_m0m1m1_P12P12P12_DP, C_m0m1m1_gp1_DP, C_m0m1m1_gP12_DP

  use triangle_reduction_qp, only: &
  C_0mm_P12_QP, C_0mm_p1_QP, C_0mm_p1p1_QP, C_0mm_p1P12_QP, C_0mm_P12P12_QP, &
  C_0mm_g_QP, C_0mm_p1p1p1_QP, C_0mm_p1p1P12_QP, C_0mm_p1P12P12_QP, &
  C_0mm_P12P12P12_QP, C_0mm_gp1_QP, C_0mm_gP12_QP, &
  C_m00_P12_QP, C_m00_p1_QP, C_m00_p1p1_QP, C_m00_p1P12_QP, C_m00_P12P12_QP, &
  C_m00_g_QP, C_m00_p1p1p1_QP, C_m00_p1p1P12_QP, C_m00_p1P12P12_QP, &
  C_m00_P12P12P12_QP, C_m00_gp1_QP, C_m00_gP12_QP, &
  C_mmm_P12_QP, C_mmm_p1_QP, C_mmm_p1p1_QP, C_mmm_p1P12_QP, C_mmm_P12P12_QP, &
  C_mmm_g_QP, C_mmm_p1p1p1_QP, C_mmm_p1p1P12_QP, C_mmm_p1P12P12_QP, &
  C_mmm_P12P12P12_QP, C_mmm_gp1_QP, C_mmm_gP12_QP, &
  C_000_P12_QP, C_000_p1_QP, C_000_p1p1_QP, C_000_p1P12_QP, C_000_P12P12_QP, &
  C_000_g_QP, C_000_p1p1p1_QP, C_000_p1p1P12_QP, C_000_p1P12P12_QP, &
  C_000_P12P12P12_QP, C_000_gp1_QP, C_000_gP12_QP, &
  C_m0m1m1_p1_QP, C_m0m1m1_P12_QP, &
  C_m0m1m1_p1p1_QP, C_m0m1m1_p1P12_QP, C_m0m1m1_P12P12_QP, C_m0m1m1_g_QP, &
  C_m0m1m1_p1p1p1_QP, C_m0m1m1_p1p1P12_QP, C_m0m1m1_p1P12P12_QP, &
  C_m0m1m1_P12P12P12_QP, C_m0m1m1_gp1_QP, C_m0m1m1_gP12_QP
  implicit none

  interface set_muUV2
    module procedure set_muUV2_dp, set_muUV2_qp
  end interface
  interface set_muIR2
    module procedure set_muIR2_dp, set_muIR2_qp
  end interface
  interface set_duv
    module procedure set_duv_dp, set_duv_qp
  end interface
  interface set_dir
    module procedure set_dir_dp, set_dir_qp
  end interface

  interface C_0mm_P12
    module procedure C_0mm_P12_DP, C_0mm_P12_QP
  end interface
  interface C_0mm_p1
    module procedure C_0mm_p1_DP, C_0mm_p1_QP
  end interface
  interface C_0mm_p1p1
    module procedure C_0mm_p1p1_DP, C_0mm_p1p1_QP
  end interface
  interface C_0mm_p1P12
    module procedure C_0mm_p1P12_DP, C_0mm_p1P12_QP
  end interface
  interface C_0mm_P12P12
    module procedure C_0mm_P12P12_DP, C_0mm_P12P12_QP
  end interface
  interface C_0mm_g
    module procedure C_0mm_g_DP, C_0mm_g_QP
  end interface
  interface C_0mm_p1p1p1
    module procedure C_0mm_p1p1p1_DP, C_0mm_p1p1p1_QP
  end interface
  interface C_0mm_p1p1P12
    module procedure C_0mm_p1p1P12_DP, C_0mm_p1p1P12_QP
  end interface
  interface C_0mm_p1P12P12
    module procedure C_0mm_p1P12P12_DP, C_0mm_p1P12P12_QP
  end interface
  interface C_0mm_P12P12P12
    module procedure C_0mm_P12P12P12_DP, C_0mm_P12P12P12_QP
  end interface
  interface C_0mm_gp1
    module procedure C_0mm_gp1_DP, C_0mm_gp1_QP
  end interface
  interface C_0mm_gP12
    module procedure C_0mm_gP12_DP, C_0mm_gP12_QP
  end interface
  interface C_m00_P12
    module procedure C_m00_P12_DP, C_m00_P12_QP
  end interface
  interface C_m00_p1
    module procedure C_m00_p1_DP, C_m00_p1_QP
  end interface
  interface C_m00_p1p1
    module procedure C_m00_p1p1_DP, C_m00_p1p1_QP
  end interface
  interface C_m00_p1P12
    module procedure C_m00_p1P12_DP, C_m00_p1P12_QP
  end interface
  interface C_m00_P12P12
    module procedure C_m00_P12P12_DP, C_m00_P12P12_QP
  end interface
  interface C_m00_g
    module procedure C_m00_g_DP, C_m00_g_QP
  end interface
  interface C_m00_p1p1p1
    module procedure C_m00_p1p1p1_DP, C_m00_p1p1p1_QP
  end interface
  interface C_m00_p1p1P12
    module procedure C_m00_p1p1P12_DP, C_m00_p1p1P12_QP
  end interface
  interface C_m00_p1P12P12
    module procedure C_m00_p1P12P12_DP, C_m00_p1P12P12_QP
  end interface
  interface C_m00_P12P12P12
    module procedure C_m00_P12P12P12_DP, C_m00_P12P12P12_QP
  end interface
  interface C_m00_gp1
    module procedure C_m00_gp1_DP, C_m00_gp1_QP
  end interface
  interface C_m00_gP12
    module procedure C_m00_gP12_DP, C_m00_gP12_QP
  end interface
  interface C_mmm_P12
    module procedure C_mmm_P12_DP, C_mmm_P12_QP
  end interface
  interface C_mmm_p1
    module procedure C_mmm_p1_DP, C_mmm_p1_QP
  end interface
  interface C_mmm_p1p1
    module procedure C_mmm_p1p1_DP, C_mmm_p1p1_QP
  end interface
  interface C_mmm_p1P12
    module procedure C_mmm_p1P12_DP, C_mmm_p1P12_QP
  end interface
  interface C_mmm_P12P12
    module procedure C_mmm_P12P12_DP, C_mmm_P12P12_QP
  end interface
  interface C_mmm_g
    module procedure C_mmm_g_DP, C_mmm_g_QP
  end interface
  interface C_mmm_p1p1p1
    module procedure C_mmm_p1p1p1_DP, C_mmm_p1p1p1_QP
  end interface
  interface C_mmm_p1p1P12
    module procedure C_mmm_p1p1P12_DP, C_mmm_p1p1P12_QP
  end interface
  interface C_mmm_p1P12P12
    module procedure C_mmm_p1P12P12_DP, C_mmm_p1P12P12_QP
  end interface
  interface C_mmm_P12P12P12
    module procedure C_mmm_P12P12P12_DP, C_mmm_P12P12P12_QP
  end interface
  interface C_mmm_gp1
    module procedure C_mmm_gp1_DP, C_mmm_gp1_QP
  end interface
  interface C_mmm_gP12
    module procedure C_mmm_gP12_DP, C_mmm_gP12_QP
  end interface
  interface C_000_P12
    module procedure C_000_P12_DP, C_000_P12_QP
  end interface
  interface C_000_p1
    module procedure C_000_p1_DP, C_000_p1_QP
  end interface
  interface C_000_p1p1
    module procedure C_000_p1p1_DP, C_000_p1p1_QP
  end interface
  interface C_000_p1P12
    module procedure C_000_p1P12_DP, C_000_p1P12_QP
  end interface
  interface C_000_P12P12
    module procedure C_000_P12P12_DP, C_000_P12P12_QP
  end interface
  interface C_000_g
    module procedure C_000_g_DP, C_000_g_QP
  end interface
  interface C_000_p1p1p1
    module procedure C_000_p1p1p1_DP, C_000_p1p1p1_QP
  end interface
  interface C_000_p1p1P12
    module procedure C_000_p1p1P12_DP, C_000_p1p1P12_QP
  end interface
  interface C_000_p1P12P12
    module procedure C_000_p1P12P12_DP, C_000_p1P12P12_QP
  end interface
  interface C_000_P12P12P12
    module procedure C_000_P12P12P12_DP, C_000_P12P12P12_QP
  end interface
  interface C_000_gp1
    module procedure C_000_gp1_DP, C_000_gp1_QP
  end interface
  interface C_000_gP12
    module procedure C_000_gP12_DP, C_000_gP12_QP
  end interface
  interface C_m0m1m1_p1
    module procedure C_m0m1m1_p1_DP, C_m0m1m1_p1_QP
  end interface
  interface C_m0m1m1_P12
    module procedure C_m0m1m1_P12_DP, C_m0m1m1_P12_QP
  end interface
  interface C_m0m1m1_p1p1
    module procedure C_m0m1m1_p1p1_DP, C_m0m1m1_p1p1_QP
  end interface
  interface C_m0m1m1_p1P12
    module procedure C_m0m1m1_p1P12_DP, C_m0m1m1_p1P12_QP
  end interface
  interface C_m0m1m1_P12P12
    module procedure C_m0m1m1_P12P12_DP, C_m0m1m1_P12P12_QP
  end interface
  interface C_m0m1m1_g
    module procedure C_m0m1m1_g_DP, C_m0m1m1_g_QP
  end interface
  interface C_m0m1m1_p1p1p1
    module procedure C_m0m1m1_p1p1p1_DP, C_m0m1m1_p1p1p1_QP
  end interface
  interface C_m0m1m1_p1p1P12
    module procedure C_m0m1m1_p1p1P12_DP, C_m0m1m1_p1p1P12_QP
  end interface
  interface C_m0m1m1_p1P12P12
    module procedure C_m0m1m1_p1P12P12_DP, C_m0m1m1_p1P12P12_QP
  end interface
  interface C_m0m1m1_P12P12P12
    module procedure C_m0m1m1_P12P12P12_DP, C_m0m1m1_P12P12P12_QP
  end interface
  interface C_m0m1m1_gp1
    module procedure C_m0m1m1_gp1_DP, C_m0m1m1_gp1_QP
  end interface
  interface C_m0m1m1_gP12
    module procedure C_m0m1m1_gP12_DP, C_m0m1m1_gP12_QP
  end interface

end module trred
