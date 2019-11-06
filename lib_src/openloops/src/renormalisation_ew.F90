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


module ol_ew_renormalisation_/**/REALKIND
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), save      ::   IRrational = -1 ! Rational contribution from IR divergent fermion self energies
  contains

! **********************************************************************
subroutine ew_renormalisation
! **********************************************************************
! EW renormalisation and R2 constants (in physical GeV units);
! conventions for UV/IR div. see subroutine loop_parameters_init.
! This subroutine is automatically (once) called by loop_parameters_init.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_msg, ol_error, ol_fatal
  use ol_generic, only: to_string
  use ol_parameters_decl_/**/REALKIND
#ifndef PRECISION_dp
  use ol_parameters_decl_/**/DREALKIND, only: LeadingColour, ew_renorm_scheme, cms_on, delta_alphamz_dimreg
#endif
  use ol_loop_parameters_decl_/**/REALKIND
  use ol_self_energy_integrals_/**/REALKIND
#ifndef PRECISION_dp
  use ol_loop_parameters_decl_/**/DREALKIND, only: &
    & nc, nf, N_lf, N_lu, N_ld, N_ll, nq_nondecoupl, CT_is_on, R2_is_on, TP_is_on, SwF, SwB
#endif
  implicit none

  real(REALKIND) :: deB_UV, deB_IR, deT_UV, deT_IR
  real(REALKIND) :: eps=1.e-17
  logical :: zeromasses(6)  ! maximal allowed number of quarks
  integer :: debug_ew_renorm = 0

!for debugging with ABC
  real(REALKIND), save         :: debug_norm = 1._/**/REALKIND

! self energies
  complex(REALKIND), save      ::   Tadpole   = 0
  complex(REALKIND), save      ::   SiW    = 0
  complex(REALKIND), save      ::   SiW0  = 0
  complex(REALKIND), save      ::   dSiW   = 0
  complex(REALKIND), save      ::   SiZZ   = 0
  complex(REALKIND), save      ::   dSiZZ  = 0
  complex(REALKIND), save      ::   SiAZZ  = 0
  complex(REALKIND), save      ::   dSiAZZ  = 0
  complex(REALKIND), save      ::   SiAZ0  = 0
  complex(REALKIND), save      ::   dSiAAheavy0  = 0
  complex(REALKIND), save      ::   dSiAA0  = 0
  complex(REALKIND), save      ::   PiAAlightZ  = 0
  complex(REALKIND), save      ::   dAlphaQED_MZ = 0
  complex(REALKIND), save      ::   dSiAAZ  = 0
  complex(REALKIND), save      ::   SiH    = 0
  complex(REALKIND), save      ::   dSiH   = 0
  complex(REALKIND), save      ::   SitL   = 0
  complex(REALKIND), save      ::   SitR   = 0
  complex(REALKIND), save      ::   SitS   = 0
  complex(REALKIND), save      ::   dSitL  = 0
  complex(REALKIND), save      ::   dSitR  = 0
  complex(REALKIND), save      ::   dSitS  = 0
  complex(REALKIND), save      ::   SibL   = 0
  complex(REALKIND), save      ::   SibR   = 0
  complex(REALKIND), save      ::   SibS   = 0
  complex(REALKIND), save      ::   dSibL  = 0
  complex(REALKIND), save      ::   dSibR  = 0
  complex(REALKIND), save      ::   dSibS  = 0
  complex(REALKIND), save      ::   SiuL   = 0
  complex(REALKIND), save      ::   SiuR   = 0
  complex(REALKIND), save      ::   SidL   = 0
  complex(REALKIND), save      ::   SidR   = 0
  complex(REALKIND), save      ::   SieL   = 0
  complex(REALKIND), save      ::   SieR   = 0
  complex(REALKIND), save      ::   SinL   = 0
  complex(REALKIND), save      ::   SinlL  = 0
  complex(REALKIND), save      ::   SilL   = 0
  complex(REALKIND), save      ::   SilR   = 0
  complex(REALKIND), save      ::   SilS   = 0
  complex(REALKIND), save      ::   dSilL  = 0
  complex(REALKIND), save      ::   dSilR  = 0
  complex(REALKIND), save      ::   dSilS  = 0

  complex(REALKIND), save      ::   dSitLRS = 0
  complex(REALKIND), save      ::   dSibLRS = 0
  complex(REALKIND), save      ::   dSilLRS = 0

! one-point functions
  complex(REALKIND), save      ::   A0W = 0
  complex(REALKIND), save      ::   A0Z = 0
  complex(REALKIND), save      ::   A0H = 0
  complex(REALKIND), save      ::   A0T = 0
  complex(REALKIND), save      ::   A0B = 0
  complex(REALKIND), save      ::   A0L = 0

! two-point
  complex(REALKIND), save      ::  B0000   = 0
  complex(REALKIND), save      ::  B00WW   = 0
  complex(REALKIND), save      ::  dB00WW  = 0
  complex(REALKIND), save      ::  B00W0   = 0
  complex(REALKIND), save      ::  dB00W0  = 0
  complex(REALKIND), save      ::  B0WW0   = 0
  complex(REALKIND), save      ::  dB0WW0  = 0
  complex(REALKIND), save      ::  B0ZWW   = 0
  complex(REALKIND), save      ::  dB0ZWW  = 0
  complex(REALKIND), save      ::  B00HH   = 0
  complex(REALKIND), save      ::  B00ZZ   = 0
  complex(REALKIND), save      ::  B0ZZH   = 0
  complex(REALKIND), save      ::  dB0ZZH  = 0
  complex(REALKIND), save      ::  B00ZH   = 0
  complex(REALKIND), save      ::  B00WH   = 0
  complex(REALKIND), save      ::  dB00WH  = 0
  complex(REALKIND), save      ::  B0WWH   = 0
  complex(REALKIND), save      ::  dB0WWH  = 0
  complex(REALKIND), save      ::  B00WZ   = 0
  complex(REALKIND), save      ::  dB00WZ  = 0
  complex(REALKIND), save      ::  B0WWZ   = 0
  complex(REALKIND), save      ::  dB0WWZ  = 0
  complex(REALKIND), save      ::  B0HHH   = 0
  complex(REALKIND), save      ::  dB0HHH  = 0
  complex(REALKIND), save      ::  B0HWW   = 0
  complex(REALKIND), save      ::  dB0HWW  = 0
  complex(REALKIND), save      ::  B0HZZ   = 0
  complex(REALKIND), save      ::  dB0HZZ  = 0
  complex(REALKIND), save      ::  B0BTW   = 0
  complex(REALKIND), save      ::  dB0BTW  = 0
  complex(REALKIND), save      ::  B0TBW   = 0
  complex(REALKIND), save      ::  dB0TBW  = 0
  complex(REALKIND), save      ::  B0WTB   = 0
  complex(REALKIND), save      ::  dB0WTB  = 0
  complex(REALKIND), save      ::  B1TBW   = 0
  complex(REALKIND), save      ::  dB1TBW  = 0
  complex(REALKIND), save      ::  B1BTW   = 0
  complex(REALKIND), save      ::  dB1BTW  = 0
  complex(REALKIND), save      ::  B0TT0   = 0
  complex(REALKIND), save      ::  dB0TT0  = 0
  complex(REALKIND), save      ::  B00TT   = 0
  complex(REALKIND), save      ::  dB00TT  = 0
  complex(REALKIND), save      ::  B1TT0   = 0
  complex(REALKIND), save      ::  dB1TT0  = 0
  complex(REALKIND), save      ::  B0TTZ   = 0
  complex(REALKIND), save      ::  dB0TTZ  = 0
  complex(REALKIND), save      ::  B1TTZ   = 0
  complex(REALKIND), save      ::  dB1TTZ  = 0
  complex(REALKIND), save      ::  B0ZTT  = 0
  complex(REALKIND), save      ::  dB0ZTT  = 0
  complex(REALKIND), save      ::  B0TTH   = 0
  complex(REALKIND), save      ::  dB0TTH  = 0
  complex(REALKIND), save      ::  B0HTT   = 0
  complex(REALKIND), save      ::  dB0HTT  = 0
  complex(REALKIND), save      ::  B1TTH   = 0
  complex(REALKIND), save      ::  dB1TTH  = 0
  complex(REALKIND), save      ::  B0BB0   = 0
  complex(REALKIND), save      ::  dB0BB0  = 0
  complex(REALKIND), save      ::  B00BB   = 0
  complex(REALKIND), save      ::  dB00BB  = 0
  complex(REALKIND), save      ::  B1BB0   = 0
  complex(REALKIND), save      ::  dB1BB0  = 0
  complex(REALKIND), save      ::  B0BBZ   = 0
  complex(REALKIND), save      ::  dB0BBZ  = 0
  complex(REALKIND), save      ::  B1BBZ   = 0
  complex(REALKIND), save      ::  dB1BBZ  = 0
  complex(REALKIND), save      ::  B0BBH   = 0
  complex(REALKIND), save      ::  dB0BBH  = 0
  complex(REALKIND), save      ::  B0HBB   = 0
  complex(REALKIND), save      ::  dB0HBB  = 0
  complex(REALKIND), save      ::  B1BBH   = 0
  complex(REALKIND), save      ::  dB1BBH  = 0
  complex(REALKIND), save      ::  B00TB   = 0
  complex(REALKIND), save      ::  dB00TB  = 0
  complex(REALKIND), save      ::  B100W   = 0
  complex(REALKIND), save      ::  B1000   = 0
  complex(REALKIND), save      ::  B100Z   = 0
  complex(REALKIND), save      ::  B0Z00   = 0
  complex(REALKIND), save      ::  dB0Z00  = 0
  complex(REALKIND), save      ::  B0W00   = 0
  complex(REALKIND), save      ::  dB0W00  = 0
  complex(REALKIND), save      ::  B0ZBB   = 0
  complex(REALKIND), save      ::  dB0ZBB  = 0

  complex(REALKIND), save      ::  B0L0W   = 0
  complex(REALKIND), save      ::  dB0L0W  = 0
  complex(REALKIND), save      ::  B0W0L  = 0
  complex(REALKIND), save      ::  dB0W0L  = 0
  complex(REALKIND), save      ::  B10LW   = 0
  complex(REALKIND), save      ::  dB10LW  = 0
  complex(REALKIND), save      ::  B1L0W   = 0
  complex(REALKIND), save      ::  dB1L0W  = 0
  complex(REALKIND), save      ::  B0LL0   = 0
  complex(REALKIND), save      ::  dB0LL0  = 0
  complex(REALKIND), save      ::  B00LL   = 0
  complex(REALKIND), save      ::  dB00LL  = 0
  complex(REALKIND), save      ::  B1LL0   = 0
  complex(REALKIND), save      ::  dB1LL0  = 0
  complex(REALKIND), save      ::  B0LLZ   = 0
  complex(REALKIND), save      ::  dB0LLZ  = 0
  complex(REALKIND), save      ::  B1LLZ   = 0
  complex(REALKIND), save      ::  dB1LLZ  = 0
  complex(REALKIND), save      ::  B0LLH   = 0
  complex(REALKIND), save      ::  dB0LLH  = 0
  complex(REALKIND), save      ::  B0HLL   = 0
  complex(REALKIND), save      ::  dB0HLL  = 0
  complex(REALKIND), save      ::  B1LLH   = 0
  complex(REALKIND), save      ::  dB1LLH  = 0
  complex(REALKIND), save      ::  B000L   = 0
  complex(REALKIND), save      ::  dB000L  = 0
  complex(REALKIND), save      ::  B0ZLL   = 0
  complex(REALKIND), save      ::  dB0ZLL  = 0


  complex(REALKIND), save      ::   cTW = 0
  complex(REALKIND), save      ::   cTZ = 0
  complex(REALKIND), save      ::   cTAZ = 0
  complex(REALKIND), save      ::   cTH = 0
  complex(REALKIND), save      ::   cTT = 0
  complex(REALKIND), save      ::   cTB = 0
  complex(REALKIND), save      ::   cTL = 0

!  complex(REALKIND), save      ::   lambda  = 0   ! fictous photon mass. Not used atm.

  complex(REALKIND), save      ::   dZuLPole   = 0.
  complex(REALKIND), save      ::   dZuRPole   = 0.
  complex(REALKIND), save      ::   dZtLPole   = 0.
  complex(REALKIND), save      ::   dZtRPole   = 0.

  complex(REALKIND), save      ::   dZeQEDEWPole, dcwEWPole, dswEWPole, dZAAEWPole, dZZAEWPole, dZAZEWPole
  complex(REALKIND), save      ::   dZMZ2EWPole, dZMW2EWPole, dZZZEWPole, dZWEWPole


  if (LeadingColour == 0) then
    cf = (nc**2-1)/(2.*nc)
    ca = nc
    tf = 0.5
  else
    cf = 0.5*nc
    ca = nc
    tf = 0
  end if

  zeromasses = [MU == 0, MD == 0, MS == 0, MC == 0, MB == 0, MT == 0]
  N_lf = count(zeromasses(:nf))

  zeromasses(1:3) = [MU == 0, MC == 0, MT == 0]
  N_lu = count(zeromasses(:3))

  zeromasses(1:3) = [MD == 0, MS == 0, MB == 0]
  N_ld = count(zeromasses(:3))

  zeromasses(1:3) = [ME == 0, MM == 0, ML == 0]
  N_ll = count(zeromasses(:3))

  Qf2sum = nc*N_lu*Qu2+nc*N_ld*Qd2+N_ll*Ql2

  if (N_lf /= 4 .and. N_lf /= 5) then
    call ol_error(2, 'in ew_renormalisation:')
    call ol_msg( 'N_lf = ' // to_string(N_lf) // 'is not supported.')
    call ol_fatal()
  end if

  Tadpole   = 0.
  SiW    = 0.
  SiW0   = 0.
  dSiW   = 0.
  SiZZ   = 0.
  dSiZZ  = 0.
  SiAZZ  = 0.
  dSiAZZ  = 0.
  SiAZ0  = 0.
  dSiAAheavy0  = 0.
  dSiAA0 = 0.
  PiAAlightZ  = 0.
  dSiAAZ = 0.
  SiH    = 0.
  dSiH   = 0.
  SitL   = 0.
  SitR   = 0.
  SitS   = 0.
  dSitL  = 0.
  dSitR  = 0.
  dSitS  = 0.
  SibL   = 0.
  SibR   = 0.
  SibS   = 0.
  dSibL  = 0.
  dSibR  = 0.
  dSibS  = 0.
  SiuL   = 0.
  SiuR   = 0.
  SidL   = 0.
  SidR   = 0.
  SieL   = 0.
  SieR   = 0.
  SinL   = 0.
  SinlL  = 0.
  SilL   = 0.
  SilR   = 0.
  SilS   = 0.
  dSilL  = 0.
  dSilR  = 0.
  dSilS  = 0.


! Do not use complex masses if CMS is not used --> only for debug & model=sm_vaux
  if ( cms_on == 0 ) then
    call masspowers(rME, 0._/**/REALKIND, ME, ME2, rME2)
    call masspowers(rMM, 0._/**/REALKIND, MM, MM2, rMM2)
    call masspowers(rML, 0._/**/REALKIND, ML, ML2, rML2)
    call masspowers(rMU, 0._/**/REALKIND, MU, MU2, rMU2)
    call masspowers(rMD, 0._/**/REALKIND, MD, MD2, rMD2)
    call masspowers(rMS, 0._/**/REALKIND, MS, MS2, rMS2)
    call masspowers(rMC, 0._/**/REALKIND, MC, MC2, rMC2)
    call masspowers(rMB, 0._/**/REALKIND, MB, MB2, rMB2)
    call masspowers(rMT, 0._/**/REALKIND, MT, MT2, rMT2)
    call masspowers(rMW, 0._/**/REALKIND, MW, MW2, rMW2)
    call masspowers(rMZ, 0._/**/REALKIND, MZ, MZ2, rMZ2)
    call masspowers(rMH, 0._/**/REALKIND, MH, MH2, rMH2)
  end if


  ! calculate one- and two-point functions
  call init_ol_self_energy_integrals(.true.)
  A0W = calcA0(MW2)
  A0Z = calcA0(MZ2)
  A0H = calcA0(MH2)
  A0T = calcA0(MT2)
  A0B = calcA0(MB2)
  A0L = calcA0(ML2)

  B0000   = calcB0(ZERO,ZERO,ZERO)
  B00WW   = calcB0(ZERO,MW2,MW2)
  dB00WW  = calcdB0(ZERO,MW2,MW2)
  B00W0   = calcB0(ZERO,MW2,ZERO)
  dB00W0  = calcdB0(ZERO,MW2,ZERO)
  B00HH   = calcB0(ZERO,MH2,MH2)
  B00ZZ   = calcB0(ZERO,MZ2,MZ2)
  B00ZH   = calcB0(ZERO,MZ2,MH2)
  B00WH   = calcB0(ZERO,MW2,MH2)
  dB00WH  = calcdB0(ZERO,MW2,MH2)
  B00WZ   = calcB0(ZERO,MW2,MZ2)
  dB00WZ  = calcdB0(ZERO,MW2,MZ2)
  B00TT   = calcB0(ZERO,MT2,MT2)
  dB00TT  = calcdB0(ZERO,MT2,MT2)
  B00BB   = calcB0(ZERO,MB2,MB2)
  dB00BB  = calcdB0(ZERO,MB2,MB2)
  B00TB   = calcB0(ZERO,MT2,MB2)
  dB00TB  = calcdB0(ZERO,MT2,MB2)
  B100W   = calcB1(ZERO,ZERO,MW2)
  B1000   = calcB1(ZERO,ZERO,ZERO)
  B100Z   = calcB1(ZERO,ZERO,MZ2)
  B10LW   = calcB1(ZERO,ML2,MW2)
  dB10LW  = calcdB1(ZERO,ML2,MW2)
  B000L   = calcB0(ZERO,ZERO,ML2)
  dB000L  = calcdB0(ZERO,ZERO,ML2)
  B00LL   = calcB0(ZERO,ML2,ML2)
  dB00LL  = calcdB0(ZERO,ML2,ML2)

  B0WW0   = calcB0(MW2,MW2,ZERO)
  dB0WW0  = calcdB0(MW2,MW2,ZERO)
  B0ZZH   = calcB0(MZ2,MZ2,MH2)
  dB0ZZH  = calcdB0(MZ2,MZ2,MH2)
  B0WWH   = calcB0(MW2,MW2,MH2)
  dB0WWH  = calcdB0(MW2,MW2,MH2)
  B0WWZ   = calcB0(MW2,MW2,MZ2)
  dB0WWZ  = calcdB0(MW2,MW2,MZ2)
  B0HHH   = calcB0(MH2,MH2,MH2)
  dB0HHH  = calcdB0(MH2,MH2,MH2)
  B0TT0   = calcB0(MT2,MT2,ZERO)
  dB0TT0  = calcdB0(MT2,MT2,ZERO)
  B1TT0   = calcB1(MT2,MT2,ZERO)
  dB1TT0  = calcdB1(MT2,MT2,ZERO)
  B0TTZ   = calcB0(MT2,MT2,MZ2)
  dB0TTZ  = calcdB0(MT2,MT2,MZ2)
  B1TTZ   = calcB1(MT2,MT2,MZ2)
  dB1TTZ  = calcdB1(MT2,MT2,MZ2)
  B0TTH   = calcB0(MT2,MT2,MH2)
  dB0TTH  = calcdB0(MT2,MT2,MH2)
  B1TTH   = calcB1(MT2,MT2,MH2)
  dB1TTH  = calcdB1(MT2,MT2,MH2)
  B0BB0   = calcB0(MB2,MB2,ZERO)
  dB0BB0  = calcdB0(MB2,MB2,ZERO)
  B1BB0   = calcB1(MB2,MB2,ZERO)
  dB1BB0  = calcdB1(MB2,MB2,ZERO)
  B0BBZ   = calcB0(MB2,MB2,MZ2)
  dB0BBZ  = calcdB0(MB2,MB2,MZ2)
  B1BBZ   = calcB1(MB2,MB2,MZ2)
  dB1BBZ  = calcdB1(MB2,MB2,MZ2)
  B0BBH   = calcB0(MB2,MB2,MH2)
  dB0BBH  = calcdB0(MB2,MB2,MH2)
  B1BBH   = calcB1(MB2,MB2,MH2)
  dB1BBH  = calcdB1(MB2,MB2,MH2)
  B0LL0   = calcB0(ML2,ML2,ZERO)
  dB0LL0  = calcdB0(ML2,ML2,ZERO)
  B1LL0   = calcB1(ML2,ML2,ZERO)
  dB1LL0  = calcdB1(ML2,ML2,ZERO)
  B0LLZ   = calcB0(ML2,ML2,MZ2)
  dB0LLZ  = calcdB0(ML2,ML2,MZ2)
  B1LLZ   = calcB1(ML2,ML2,MZ2)
  dB1LLZ  = calcdB1(ML2,ML2,MZ2)
  B0LLH   = calcB0(ML2,ML2,MH2)
  dB0LLH  = calcdB0(ML2,ML2,MH2)
  B1LLH   = calcB1(ML2,ML2,MH2)
  dB1LLH  = calcdB1(ML2,ML2,MH2)

  B0ZWW   = calcRB0(MZ2,MW2,MW2)
  dB0ZWW  = calcRdB0(MZ2,MW2,MW2)
  B0HWW   = calcRB0(MH2,MW2,MW2)
  dB0HWW  = calcRdB0(MH2,MW2,MW2)
  B0HZZ   = calcRB0(MH2,MZ2,MZ2)
  dB0HZZ  = calcRdB0(MH2,MZ2,MZ2)
  B0BTW   = calcRB0(MB2,MT2,MW2)
  dB0BTW  = calcRdB0(MB2,MT2,MW2)
  B0TBW   = calcRB0(MT2,MB2,MW2)
  dB0TBW  = calcRdB0(MT2,MB2,MW2)
  B0WTB   = calcRB0(MW2,MT2,MB2)
  dB0WTB  = calcRdB0(MW2,MT2,MB2)
  B1TBW   = calcRB1(MT2,MB2,MW2)
  dB1TBW  = calcRdB1(MT2,MB2,MW2)
  B1BTW   = calcRB1(MB2,MT2,MW2)
  dB1BTW  = calcRdB1(MB2,MT2,MW2)
  B0ZTT   = calcRB0(MZ2,MT2,MT2)
  dB0ZTT  = calcRdB0(MZ2,MT2,MT2)
  B0HTT   = calcRB0(MH2,MT2,MT2)
  dB0HTT  = calcRdB0(MH2,MT2,MT2)
  B0HBB   = calcRB0(MH2,MB2,MB2)
  dB0HBB  = calcRdB0(MH2,MB2,MB2)
  B0Z00   = calcRB0(MZ2,ZERO,ZERO)
  dB0Z00  = calcRdB0(MZ2,ZERO,ZERO)
  B0W00   = calcRB0(MW2,ZERO,ZERO)
  dB0W00  = calcRdB0(MW2,ZERO,ZERO)
  B0ZBB   = calcRB0(MZ2,MB2,MB2)
  dB0ZBB  = calcRdB0(MZ2,MB2,MB2)
  B0L0W   = calcRB0(ML2,ZERO,MW2)
  dB0L0W  = calcRdB0(ML2,ZERO,MW2)
  B0W0L   = calcRB0(MW2,ZERO,ML2)
  dB0W0L  = calcRdB0(MW2,ZERO,ML2)
  B1L0W   = calcRB1(ML2,ZERO,MW2)
  dB1L0W  = calcRdB1(ML2,ZERO,MW2)
  B0HLL   = calcRB0(MH2,ML2,ML2)
  dB0HLL  = calcRdB0(MH2,ML2,ML2)
  B0ZLL   = calcRB0(MZ2,ML2,ML2)
  dB0ZLL  = calcRdB0(MZ2,ML2,ML2)
  call init_ol_self_energy_integrals(.false.)


  ! calculate renormalisation constants
    if (SwB /= 0 .and. CT_is_on /= 0) then
    ! non-fermionic!

    !Bosonic Tadpole. Agrees with eq. 12 in NP B 383 (1992) 73
    Tadpole   = 0.75*MH2/MW2/sw*A0H                   &
              + (0.5*MH2/MW2/sw+3./sw)*A0W            &
              + (0.25*MH2/MW2/sw + 1.5/sw/cw2)*A0Z    &
              - (2*MW2)/sw - MZ2/(cw2*sw) ! rational contribution

    !This is the bosonic part of page 105-107 from Denner92 evaluated in selfenergies.nb

    dSiAAheavy0=dSiAAheavy0-3*B00WW-4*MW2*dB00WW

    dSiAA0 = dSiAA0-3*B00WW-4*MW2*dB00WW

    dSiAAZ=dSiAAZ-3*B0ZWW-(4*MW2+3*rMZ2)*dB0ZWW

    SiAZ0=SiAZ0+(2*MW2*B00WW)/(cw*sw)

    SiAZZ=SiAZZ+(rMZ2/3. - (-2 + 12*cw2)*MW2*B00WW &
        + ((4 + 12*cw2)*MW2 + (0.5 + 9*cw2)*rMZ2)*B0ZWW)/(3.*cw*sw)

    dSiAZZ=dSiAZZ+(2 + (3 + 54*cw2)*B0ZWW  &
        +  3*(8*MW2 + 24*cw2*MW2 + rMZ2 + 18*cw2*rMZ2)*dB0ZWW)/(18.*cw*sw)

    SiZZ=SiZZ-(((-1 + 4*cw2)*rMZ2)/3. - (2 - 8*cw2 + 24*cw4)*MW2*B00WW  &
        + ((-10 + 16*cw2 + 24*cw4)*MW2 + (-0.5 + 2*cw2 + 18*cw4)*rMZ2)*B0ZWW)/(6.*cw2*sw2) &
        - ((-2*rMZ2)/3. - 2*MH2*B00HH - 2*MZ2*B00ZZ + (2*MH2 - 10*MZ2 - rMZ2)*B0ZZH &
        - ((-MH2 + MZ2)**2*(-B00ZH + B0ZZH))/rMZ2)/(12.*cw2*sw2)

    dSiZZ=dSiZZ+(2._/**/REALKIND/3. &
        + ((MH2 - MZ2)**2*(B00ZH - B0ZZH))/rMZ2**2 + B0ZZH &
        + (2 - 8*cw2 - 3*(-1 + 4*cw2 + 36*cw4)*B0ZWW &
        - 3*(4*(-5 + 8*cw2 + 12*cw4)*MW2 + (-1 + 4*cw2 + 36*cw4)*rMZ2)*dB0ZWW)/3.  &
        - (2*MH2 - 10*MZ2 - rMZ2)*dB0ZZH &
        + ((MH2 - MZ2)**2*dB0ZZH)/rMZ2)/(12.*cw2*sw2)

    SiW=SiW+(-8*(rMW2/3. - 2*MW2*B00WW + (MW2**2*(B00W0 - B0WW0))/rMW2 &
        + (2*MW2 + 5*rMW2)*B0WW0) - ((-2*rMW2)/3. - 2*MH2*B00HH - 2*MW2*B00WW  &
        + ((MH2 - MW2)**2*(B00WH - B0WWH))/rMW2 + (2*MH2 - 10*MW2 - rMW2)*B0WWH)/sw2  &
        - ((2*(-1 + 4*cw2)*rMW2)/3. - 2*(1 + 8*cw2)*(MW2*B00WW + MZ2*B00ZZ) &
        + ((1 + 8*cw2)*(MW2 - MZ2)**2*(B00WZ - B0WWZ))/rMW2 &
        + ((54 - 10/cw2 + 16*cw2)*MW2 + (-1 + 40*cw2)*rMW2)*B0WWZ)/sw2)/12.

    dSiW=dSiW+(-2*(1._/**/REALKIND/3. + 5*B0WW0 + (MW2**2*(-B00W0 + B0WW0))/rMW2**2  &
        - (MW2**2*dB0WW0)/rMW2 + (2*MW2 + 5*rMW2)*dB0WW0))/3. &
        - (-2._/**/REALKIND/3. - B0WWH + ((-MH2 + MW2)**2*(-B00WH + B0WWH))/rMW2**2 &
        + (2*MH2 - 10*MW2 - rMW2)*dB0WWH - ((-MH2 + MW2)**2*dB0WWH)/rMW2)/(12.*sw2) &
        - ((2*(-1 + 4*cw2))/3. + (-1 + 40*cw2)*B0WWZ &
        + ((1 + 8*cw2)*(MW2 - MZ2)**2*(-B00WZ + B0WWZ))/rMW2**2 &
        - ((1 + 8*cw2)*(MW2 - MZ2)**2*dB0WWZ)/rMW2 + ((54 - 10/cw2 + 16*cw2)*MW2 &
        + (-1 + 40*cw2)*rMW2)*dB0WWZ)/(12.*sw2)

    SiW0=SiW0+(-8*(2*MW2*B00W0 - 2*MW2*B00WW - MW2**2*dB00W0) &
          - (-2*MH2*B00HH + (2*MH2 - 10*MW2)*B00WH - 2*MW2*B00WW - (MH2 - MW2)**2*dB00WH)/sw2 &
          - ((54 - 10/cw2 + 16*cw2)*MW2*B00WZ - 2*(1 + 8*cw2)*(MW2*B00WW + MZ2*B00ZZ) &
          - (1 + 8*cw2)*(MW2 - MZ2)**2*dB00WZ)/sw2)/12.

    SiH=SiH+((3*MH2*(A0H + 3*MH2*B0HHH))/MW2 &
        + (2*(-12*MW2**2 + (MH2 + 6*MW2)*A0W  &
        + (MH2**2 + 12*MW2**2 - 4*MW2*rMH2)*B0HWW))/MW2 &
        + (-12*MZ2**2 + (MH2 + 6*MZ2)*A0Z  &
        + (MH2**2 + 12*MZ2**2 - 4*MZ2*rMH2)*B0HZZ)/(cw2*MZ2))/(8.*sw2)

    dSiH=dSiH+((9*MH2**2*lambdaHHH**2*dB0HHH)/MW2 + 4*(-2*B0HWW &
        + (MH2**2/(2.*MW2) + 6*MW2 - 2*rMH2)*dB0HWW) +(2*(-2*B0HZZ &
        + (MH2**2/(2.*MZ2) + 6*MZ2 - 2*rMH2)*dB0HZZ))/cw2)/(8.*sw2)

    SitL=SitL-(1 + (2 + MB2/MW2)*B1TBW)/(2.*sw2) - (4*(1 + 2*B1TT0))/9. &
        - (MT2*(B1TTH + B1TTZ))/(4.*MW2*sw2) - (1 + 2*B1TTZ)*gZu(2)**2

    SitR=SitR-(MT2*B1TBW)/(2.*MW2*sw2) - (4*(1 + 2*B1TT0))/9. &
        - (MT2*(B1TTH + B1TTZ))/(4.*MW2*sw2) - (1 + 2*B1TTZ)*gZu(1)**2

    SitS=SitS-(MB2*B0TBW)/(2.*MW2*sw2) - (4*(-2 + 4*B0TT0))/9. &
        - (MT2*(-B0TTH + B0TTZ))/(4.*MW2*sw2) - (-2 + 4*B0TTZ)*gZu(1)*gZu(2)

    dSitL=dSitL-((2 + MB2/MW2)*dB1TBW)/(2.*sw2) - (8*dB1TT0)/9. &
        - (MT2*(dB1TTH + dB1TTZ))/(4.*MW2*sw2) - 2*dB1TTZ*gZu(2)**2

    dSitR=dSitR-(MT2*dB1TBW)/(2.*MW2*sw2) - (8*dB1TT0)/9. &
        - (MT2*(dB1TTH + dB1TTZ))/(4.*MW2*sw2) - 2*dB1TTZ*gZu(1)**2

    dSitS=dSitS-(MB2*dB0TBW)/(2.*MW2*sw2) - (16*dB0TT0)/9. &
        - (MT2*(-dB0TTH + dB0TTZ))/(4.*MW2*sw2) - 4*dB0TTZ*gZu(1)*gZu(2)

    if (MB /= 0._/**/REALKIND) then
      SibL=SibL+(-1 - 2*B1BB0)/9. &          ! Photon
          - (1 + 2*B1BBZ)*gZd(2)**2 &                     ! Z
          - (1 + (2 + MT2/MW2)*B1BTW)/(2.*sw2) &          ! W
          - (MB2*(B1BBH + B1BBZ))/(4.*MW2*sw2)            ! H & X

      SibR=SibR+(-1 - 2*B1BB0)/9. &           ! Photon
          - (1 + 2*B1BBZ)*gZd(1)**2 &                     ! Z
          - (MB2*B1BTW)/(2.*MW2*sw2) &                    ! W
          - (MB2*(B1BBH + B1BBZ))/(4.*MW2*sw2)            ! H & X

      SibS=SibS+(2 - 4*B0BB0)/9. & ! Photon
          - (MB2*(-B0BBH + B0BBZ))/(4.*MW2*sw2) & ! H & X
          - (MT2*B0BTW)/(2.*MW2*sw2) & ! W
          - (-2 + 4*B0BBZ)*gZd(1)*gZd(2) ! Z

    else
      SibL=SibL+(-1 - IRrational - 2*B1BB0)/9. &          ! Photon
          - (1 + 2*B1BBZ)*gZd(2)**2 &                     ! Z
          - (1 + (2 + MT2/MW2)*B1BTW)/(2.*sw2)            ! W

      SibR=SibR+(-1 - IRrational - 2*B1BB0)/9. &          ! Photon
          - (1 + 2*B1BBZ)*gZd(1)**2 &                     ! Z
          - (MB2*B1BTW)/(2.*MW2*sw2)                      ! W

    end if

    dSibL=dSibL+(-2*dB1BB0)/9. & ! Photon
        - (MB2*(dB1BBH + dB1BBZ))/(4.*MW2*sw2) & ! H & X
        - ((2 + MT2/MW2)*dB1BTW)/(2.*sw2) & ! W
        - 2*dB1BBZ*gZd(2)**2 ! Z

    dSibR=dSibR+(-2*dB1BB0)/9. & ! Photon
        - (MB2*(dB1BBH + dB1BBZ))/(4.*MW2*sw2) & ! H & X
        - (MB2*dB1BTW)/(2.*MW2*sw2) & ! W
        - 2*dB1BBZ*gZd(1)**2  ! Z

    dSibS=dSibS+(-4*dB0BB0)/9. & ! Photon
        - (MB2*(-dB0BBH + dB0BBZ))/(4.*MW2*sw2) & ! H & X
        - (MT2*dB0BTW)/(2.*MW2*sw2) & ! W
        - 4*dB0BBZ*gZd(1)*gZd(2) ! Z

    SiuL=SiuL+(-4.*(1. + IRrational + 2.*B1000))/9. &  ! Photon
             - (1 + 2*B100W)/(2.*sw2) &                ! W
             - (1 + 2*B100Z)*gZu(2)**2                 ! Z

    SiuR=SiuR+(-4.*(1. + IRrational + 2.*B1000))/9. &  ! Photon
               - (1 + 2*B100Z)*gZu(1)**2               ! Z

    SidL=SidL+(-1. - IRrational - 2.*B1000)/9. &       ! Photon
             - (1 + 2*B100W)/(2.*sw2) &                ! W
             - (1 + 2*B100Z)*gZd(2)**2                 ! Z

    SidR=SidR+(-1. - IRrational - 2.*B1000)/9. &       ! Photon
             - (1 + 2*B100Z)*gZd(1)**2                 ! Z

    SieL=SieL-1. - IRrational - 2.*B1000 &             ! Photon
             - (1 + 2*B100W)/(2.*sw2) &                ! W
             - (1 + 2*B100Z)*gZl(2)**2                 ! Z

    SieR=SieR-1. - IRrational - 2.*B1000 &             ! Photon
             - (1 + 2*B100Z)*gZl(1)**2                 ! Z

    SinL=SinL-(1. + 2.*B100W)/(2.*sw2) &               ! W
             - (1. + 2*B100Z)*gZn(2)**2                 ! Z


    if (ML /= 0._/**/REALKIND) then
      SilL=SilL-1. - 2.*B1LL0 &  ! Photon
          - (1. + 2*B1LLZ)*gZl(2)**2 & ! Z
          - (1. + 2*B1L0W)/(2.*sw2) &  ! W
          - (ML2*(B1LLH + B1LLZ))/(4.*MW2*sw2) ! H & X

      SilR=SilR-1. - 2.*B1LL0 &  ! Photon
          - (1 + 2*B1LLZ)*gZl(1)**2 & ! Z
          - (ML2*B1L0W)/(2.*MW2*sw2) & ! W
          - (ML2*(B1LLH + B1LLZ))/(4.*MW2*sw2) ! H & X

      SilS=SilS+(2. - 4.*B0LL0) & ! Photon
          - (ML2*(-B0LLH + B0LLZ))/(4.*MW2*sw2)  & ! H & X
          - (-2. + 4*B0LLZ)*gZl(1)*gZl(2) ! Z

      dSilL=dSilL-2.*dB1LL0 & ! Photon
          - (ML2*(dB1LLH + dB1LLZ))/(4.*MW2*sw2) & ! H & X
          - dB1L0W/sw2 & ! W
          - 2.*dB1LLZ*gZl(2)**2 ! Z

      dSilR=dSilR-2.*dB1LL0 & ! Photon
          - (ML2*(dB1LLH + dB1LLZ))/(4.*MW2*sw2) & ! H & X
          - (ML2*dB1L0W)/(2.*MW2*sw2) & ! W
          - 2.*dB1LLZ*gZl(1)**2        ! Z

      dSilS=dSilS-4.*dB0LL0 & ! Photon
          - (ML2*(-dB0LLH + dB0LLZ))/(4.*MW2*sw2) & ! H & X
          - (ML2*dB0L0W)/(2.*MW2*sw2) & ! W
          - 4.*dB0LLZ*gZl(1)*gZl(2)    ! Z

      SinlL=SinlL-(1. + (2.+ML2/MW2)*B10LW)/(2.*sw2) &  ! W
             - (1 + 2*B100Z)*gZn(2)**2                  ! Z

    else
      SilL=SieL
      SilR=SieR
      SinlL=SinL
    end if

  end if

  if (SwF /= 0 .and. CT_is_on /= 0) then
    ! fermionic

    !Fermionic Tadpole. Agrees with eq. 12 in NP B 383 (1992) 73
    Tadpole   = Tadpole - nc*(2*MT2/MW2/sw*A0T + 2*MB2/MW2/sw*A0B) - 2.*ML2/MW2/sw*A0L

!This is the fermionic part based on page 105-107 from Denner92

    ! contributions from light fermions incl. bottom
    PiAAlightZ=PiAAlightZ  &
      -     2.*4.*rMZ2*(1._/**/REALKIND/3. - B0Z00)/3. &      ! light-leptons
      -        4.*(rMZ2*(1._/**/REALKIND/3. - B0ZLL) - 2.*ML2*(B0ZLL-B00LL))/3. & ! tau-lepton
      -   2*20*nc*rMZ2*(1._/**/REALKIND/3. - B0Z00)/27. &       ! light-quarks
      -      4*nc*(rMZ2*(1._/**/REALKIND/3. - B0ZBB) - 2.*MB2*(B0ZBB-B00BB))/27. ! b-quark

    ! contributions from heavy fermions = top-quark
    dSiAAheavy0=dSiAAheavy0-(16*nc*(1-3*B00TT-6*MT2*dB00TT))/81.

    dSiAAZ=dSiAAZ &
         - 2.*4.*(1._/**/REALKIND/3.-B0Z00-rMZ2*dB0Z00)/3. &               ! light-leptons
         -        4.*(1._/**/REALKIND/3.-B0ZLL-(2.*ML2+rMZ2)*dB0ZLL)/3. &  ! tau-lepton
         - 2.*20.*nc*(1._/**/REALKIND/3.-B0Z00-rMZ2*dB0Z00)/27. &          ! light-quarks
         -     4.*nc*(1._/**/REALKIND/3.-B0ZBB-(2.*MB2+rMZ2)*dB0ZBB)/27. & ! b-quark
         -    16.*nc*(1._/**/REALKIND/3.-B0ZTT-(2.*MT2+rMZ2)*dB0ZTT)/27.   ! t-quark

    ! contribution from heavy + light fermions (in dimreg)
    dSiAA0=dSiAA0+4._/**/REALKIND/3.*( &
                                 (  2.*B0000) & ! light-leptons
          +                      (     B00LL) & ! tau-lepton
          + nc*1._/**/REALKIND/9*( 2.*B0000) & ! down-quarks
          + nc*4._/**/REALKIND/9*( 2.*B0000) & ! up-quarks
          + nc*1._/**/REALKIND/9*( B00BB ) &    ! b_quark
          + nc*4._/**/REALKIND/9*( B00TT ) &    ! t-quark
          )

!      SiAZ0=SiAZ0+0 !no fermionic contribution

    SiAZZ=SiAZZ-2._/**/REALKIND/3.*( &
             2.*(rMZ2/3. - rMZ2*B0Z00)*(gZl(1) + gZl(2))   & ! light-lepton
        +       (rMZ2/3. + 2*ML2*B00LL + (-2*ML2 - rMZ2)*B0ZLL)*(gZl(1) + gZl(2)) & ! tau-lepton
        + 2.*nc*(rMZ2/3. - rMZ2*B0Z00)*(gZd(1) + gZd(2))/3.  & ! down-quark
        - 2.*nc*(rMZ2/3. - rMZ2*B0Z00)*(gZu(1) + gZu(2))*2._/**/REALKIND/3. & ! up-quark
        +   nc*(rMZ2/3. + 2*MB2*B00BB + (-2*MB2 - rMZ2)*B0ZBB)*(gZd(1) + gZd(2))/3. &     ! b-quark
        -   nc*(rMZ2/3. + 2*MT2*B00TT + (-2*MT2 - rMZ2)*B0ZTT)*(gZu(1) + gZu(2))*2._/**/REALKIND/3. &  ! top-quark
        )

    dSiAZZ=dSiAZZ-2._/**/REALKIND/9.*( &
        -    2.*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZl(1) + gZl(2)) & ! light-lepton
        -       (-1 + 3*B0ZLL + 3*(2*ML2 + rMZ2)*dB0ZLL)*(gZl(1) + gZl(2)) & ! tau-lepton
        - 2.*nc*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZd(1) + gZd(2))/3. & ! down-quark
        + 2.*nc*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZu(1) + gZu(2))*2._/**/REALKIND/3. & ! up-quark
        -    nc*(-1 + 3*B0ZBB + 3*(2*MB2 + rMZ2)*dB0ZBB)*(gZd(1) + gZd(2))/3. & ! b-quark
        +    nc*(-1 + 3*B0ZTT + 3*(2*MT2 + rMZ2)*dB0ZTT)*(gZu(1) + gZu(2))*2._/**/REALKIND/3. & ! top-quark
       )

    SiZZ=SiZZ-2._/**/REALKIND/3.*( &
        -      3*rMZ2*(-1 + 3*B0Z00)*(gZn(1)**2 + gZn(2)**2)/3. &  ! neutrinos
        -      2*rMZ2*(-1 + 3*B0Z00)*(gZl(1)**2 + gZl(2)**2)/3. &  ! light-leptons
        + (3*ML2*B0ZLL)/(4.*cw2*sw2) + (rMZ2/3. + 2*ML2*B00LL - (2*ML2 + rMZ2)*B0ZLL)*(gZl(1)**2 + gZl(2)**2) &  ! tau-lepton
        -   2*nc*rMZ2*(-1 + 3*B0Z00)*(gZd(1)**2 + gZd(2)**2)/3. & ! down-quark
        -   2*nc*rMZ2*(-1 + 3*B0Z00)*(gZu(1)**2 + gZu(2)**2)/3. & ! up-quark
        + nc*((3*MB2*B0ZBB)/(4.*cw2*sw2) + (rMZ2/3. + 2*MB2*B00BB - (2*MB2 + rMZ2)*B0ZBB)*(gZd(1)**2 + gZd(2)**2)) &  ! b-quark
        + nc*((3*MT2*B0ZTT)/(4.*cw2*sw2) + (rMZ2/3. + 2*MT2*B00TT - (2*MT2 + rMZ2)*B0ZTT)*(gZu(1)**2 + gZu(2)**2)) & ! t-quark
        )

    dSiZZ=dSiZZ-2._/**/REALKIND/3.*( &
        - 3.*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZn(1)**2 + gZn(2)**2)/3. & ! neutrinos
        - 2.*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZl(1)**2 + gZl(2)**2)/3. & ! light-leptons
        +    (3*ML2*dB0ZLL)/(4.*cw2*sw2) - ((-1 + 3*B0ZLL + 3*(2*ML2 + rMZ2)*dB0ZLL)*(gZl(1)**2 + gZl(2)**2))/3. & ! tau-lepton
        - 2.*nc*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZd(1)**2 + gZd(2)**2)/3. &  ! down-quark
        - 2.*nc*(-1 + 3*B0Z00 + 3*rMZ2*dB0Z00)*(gZu(1)**2 + gZu(2)**2)/3. &  ! up-quark
        + nc*((3*MB2*dB0ZBB)/(4.*cw2*sw2) - ((-1 + 3*B0ZBB + 3*(2*MB2 + rMZ2)*dB0ZBB)*(gZd(1)**2 + gZd(2)**2))/3.) & ! b-quark
        + nc*((3*MT2*dB0ZTT)/(4.*cw2*sw2) - ((-1 + 3*B0ZTT + 3*(2*MT2 + rMZ2)*dB0ZTT)*(gZu(1)**2 + gZu(2)**2))/3.) & ! t-quark
        )

    SiW=SiW-( &
           2*rMW2*(1._/**/REALKIND/3. - B0W00) & ! light-leptons
         + (rMW2/3. + ML2*B00LL + ((ML2 - 2*rMW2)*B0W0L)/2. + (ML2**2*(-B000L + B0W0L))/(2.*rMW2)) & ! heavy-leptons
         + 2*nc*rMW2*(1._/**/REALKIND/3.-B0W00) & ! light-quarks
         + nc*(rMW2/3. + MB2*B00BB + MT2*B00TT & ! heavy-quarks
           + ((MB2 + MT2 - 2*rMW2)*B0WTB)/2.  &
           + ((MB2 - MT2)**2*(-B00TB + B0WTB))/(2.*rMW2)) &
         )/(3.*sw2)

    dSiW=dSiW+( &
         - 2*(1._/**/REALKIND/3.- B0W00 - rMW2*dB0W00) & ! light-leptons
         - 1._/**/REALKIND/2.*(2._/**/REALKIND/3.+ ML2**2*(B000L - B0W0L)/rMW2**2 - 2*B0W0L & ! heavy-leptons
                  +  ML2**2*dB0W0L/rMW2 + (ML2 - 2*rMW2)*dB0W0L) &
         - 2*nc*(1._/**/REALKIND/3. - B0W00 - rMW2*dB0W00) & ! light-quark
         -   nc/2.*(2._/**/REALKIND/3. + ((MB2 - MT2)**2* (B00TB - B0WTB))/rMW2**2 - 2*B0WTB & ! heavy-quark
                  +  ((MB2 - MT2)**2*dB0WTB)/rMW2 + (MB2 + MT2 - 2*rMW2)*dB0WTB) &
          )/(3.*sw2)

    SiW0=SiW0 &
        - 1/(3.*sw2)*(ML2*B00LL + ML2*B000L/2. + ML2**2*dB000L/2.) &
        - nc/(3.*sw2)*(MB2*B00BB + ((MB2 + MT2)*B00TB)/2. +  &
          MT2*B00TT +  ((MB2 - MT2)**2*dB00TB)/2.)

    SiH=SiH &
       -   (ML2*(2*A0L + (4*ML2 - rMH2)*B0HLL))/(2.*MW2*sw2) &  ! tau
       -(MB2*nc*(2*A0B + (4*MB2 - rMH2)*B0HBB))/(2.*MW2*sw2) &  ! b-quark
       -(MT2*nc*(2*A0T + (4*MT2 - rMH2)*B0HTT))/(2.*MW2*sw2)    ! t-quark


    dSiH=dSiH &
       -   (ML2*(-B0HLL + (4*ML2 - rMH2)*dB0HLL))/(2.*MW2*sw2) & ! tau
       -(MB2*nc*(-B0HBB + (4*MB2 - rMH2)*dB0HBB))/(2.*MW2*sw2) & ! b-quark
       -(MT2*nc*(-B0HTT + (4*MT2 - rMH2)*dB0HTT))/(2.*MW2*sw2)   ! t-quark

  end if


! set RCs

!gauge bosons
  dAlphaQED_MZ=4*pi/alpha_QED*(1-alpha_QED_0/alpha_QED_MZ)  ! running of alpha from 0 to MZ

  dZAAEWnreg = -(dSiAAheavy0 + PiAAlightZ/rMZ2 + dAlphaQED_MZ)
  dZAAEWdreg = -dSiAA0
  if (delta_alphamz_dimreg) then
    dZAAEW = dZAAEWdreg
  else
    dZAAEW = dZAAEWnreg
  end if

  dZWEW   = -(dSiW)
  if (cms_on == 0) then
    dZMW2EW = real(SiW)
  else if (imag(MW2) /= 0 .and. cms_on > 0) then ! CMS
    cTW     = (MW2-rMW2)*dSiW+4.*(rMW2-MW2)  ! 1st order expansion & C^w term
    dZMW2EW = SiW + cTW
  else if (imag(MW2) == 0  .and. cms_on > 0) then ! CMS I & II
    dZMW2EW = SiW
  end if

  dZZAEW  = 2.*SiAZ0/MZ2
  dZAZEW  = -2.*(SiAZZ/rMZ2)
  dZZZEW  = -dSiZZ
  if (cms_on == 0) then ! on-shell
    dZZAEW = real(dZZAEW)
    dZAZEW = real(dZAZEW)
    dZMZ2EW = real(SiZZ)
  else if (imag(MZ2) /= 0 .and. cms_on > 0) then ! CMS
    cTZ = -(rMZ2-MZ2)*dSiZZ ! 1st order expansion
    dZMZ2EW = SiZZ + cTZ
    cTAZ = -2.*(MZ2-rMZ2)/rMZ2*dSiAZZ ! 1st order expansion
    dZAZEW = dZAZEW + cTAZ
  else if (imag(MZ2) == 0  .and. cms_on > 0) then ! CMS-I & CMS-II
    dZMZ2EW = SiZZ
  end if


!Higgs
  dtEW    = -Tadpole
  dZHEW   = -(dSiH)
  if (imag(MH2) == 0 .and. cms_on == 0) then ! on-shell
    dZMH2EW = real(SiH)
  else if (imag(MH) /= 0 .and. cms_on > 0) then ! CMS
    cTH = -(rMH2-MH2)*dSiH ! 1st order expansion
    dZMH2EW = SiH + cTH
  else if (imag(MH2) == 0 .and. cms_on > 0) then ! CMS-I & CMS-II
    dZMH2EW = SiH
  else
    call ol_fatal('on-shell EW renormalization with finite width not supported!')
  end if

!Fermions
  !light leptons
  dZeLEW = -SieL
  dZeREW = -SieR
  dZnLEW = -SinL
  dZnlLEW = -SinlL

  !tau
  if (ML /= 0._/**/REALKIND) then
    dSilLRS = dSilL+dSilR+2.*dSilS
    dZlLEW = -(SilL+rML2*dSilLRS)
    dZlREW = -(SilR+rML2*dSilLRS)
    if (cms_on == 0) then ! on-shell
      dZlLEW = -real(dZlLEW)
      dZlREW = -real(dZlREW)
      dZMLEW  = real(dZMLEW)
    else if (imag(ML2) /= 0 .and. cms_on > 0) then ! CMS
      cTL     = -(rML2-ML2)*dSilLRS+(rML2-ML2)/rML2*1.*4. ! 1st order expansion & C^l term
      dZMlEW  = ML*0.5*(SilL+SilR+2.*SilS + cTL)
    else if (imag(ML2) == 0 .and. cms_on > 0) then ! CMS-I & CMS-II
      dZMlEW = ML*0.5*(SilL+SilR+2.*SilS)
    else
      call ol_fatal('on-shell EW renormalization with finite width not supported!')
    end if
  else
    dZMLEW = 0.
    dZlLEW = -SilL
    dZlREW = -SilR
    if (cms_on /= 2) then ! CMS-I / on-shell
      dZlLEW = real(dZlLEW)
      dZlREW = real(dZlREW)
    end if
  end if

  !heavy quarks
  !Top
  dSitLRS = dSitL+dSitR+2.*dSitS
  dZtLEW = -(SitL+rMT2*dSitLRS)
  dZtREW = -(SitR+rMT2*dSitLRS)
  if (cms_on == 0) then ! on-shell
    dZtLEW = real(dZtLEW)
    dZtREW = real(dZtREW)
    dZMTEW = real(MT*0.5*(SitL+SitR+2.*SitS))
  else if (imag(MT2) /= 0 .and. cms_on > 0) then ! CMS
    cTT    = -(rMT2-MT2)*dSitLRS+(rMT2-MT2)/rMT2*4._/**/REALKIND/9.*4. ! 1st order expansion & C^t term
    dZMTEW = MT*0.5*(SitL+SitR+2.*SitS + cTT)
  else if (imag(MT2) == 0 .and. cms_on > 0) then ! CMS-I & CMS-II
    dZMTEW = MT*0.5*(SitL+SitR+2.*SitS)
  else
    call ol_fatal('on-shell EW renormalization with finite width not supported!')
  end if

  !Bottom
  if (MB /= 0._/**/REALKIND) then
    dSibLRS = dSibL+dSibR+2.*dSibS
    dZbLEW = -(SibL+rMB2*dSibLRS)
    dZbREW = -(SibR+rMB2*dSibLRS)
    if (cms_on == 0) then ! on-shell
      dZbLEW = real(dZbLEW)
      dZbREW = real(dZbREW)
      dZMBEW  = real(MB*0.5*(SibL+SibR+2.*SibS))
    else if (imag(MB2) /= 0 .and. cms_on > 0) then ! CMS
      cTB     = -(rMB2-MB2)*dSibLRS+(rMB2-MB2)/rMB2*1._/**/REALKIND/9.*4 ! 1st order expansion & C^b term
      dZMBEW  = MB*0.5*(SibL+SibR+2.*SibS + cTB)
    else if (imag(MB2) == 0 .and. cms_on > 0) then ! CMS-I & CMS-II
      dZMBEW  = MB*0.5*(SibL+SibR+2.*SibS)
    else
      call ol_fatal('on-shell EW renormalization with finite width not supported!')
    end if
  else
    dZMBEW = 0.
    dZbLEW = -SibL
    dZbREW = -SibR
    if (cms_on /= 2) then ! CMS-I / on-shell
      dZbLEW = real(dZbLEW)
      dZbREW = real(dZbREW)
    end if
  end if

  !light quarks
  dZuLEW = -SiuL
  dZuREW = -SiuR
  dZdLEW = -SidL
  dZdREW = -SidR


! weak mixing angle
  dcwEW     = cw/2.*(dZMW2EW/MW2-dZMZ2EW/MZ2)
  dswEW     =  -1.*cw/sw*dcwEW

! charge renormalization
  ! a(0)
  dZe0QEDEWnreg = -0.5*dZAAEWnreg - sw/cw*SiAZ0/MZ2
  dZe0QEDEWdreg = -0.5*dZAAEWdreg - sw/cw*SiAZ0/MZ2

  ! Gmu
  dZeGmuQEDEW = dswEW/sw - 1./sw/cw*SiAZ0/MZ2
  dZeGmuQEDEW = dZeGmuQEDEW - 0.5/sw2*(6.+(7.-4.*sw2)/(2.*sw2)*log(cw2))
  dZeGmuQEDEW = dZeGmuQEDEW +  0.5*(dZMW2EW-SiW0)/MW2

  ! a(MZ)
  dZeZQEDEW = -0.5*(dZAAEWnreg+dAlphaQED_MZ) - sw/cw*SiAZ0/MZ2   !NB: dAlphaQED_MZ drops out -> reg independent

  if (ew_renorm_scheme == 0) then ! on-shell scheme = alpha(0) scheme
    if (delta_alphamz_dimreg) then
      dZeQEDEW = dZe0QEDEWdreg ! dim-reg
    else
      dZeQEDEW = dZe0QEDEWnreg ! n-reg
    end if
  else if (ew_renorm_scheme == 1) then ! Gmu scheme
    dZeQEDEW = dZeGmuQEDEW
  else if (ew_renorm_scheme == 2) then ! alpha(mZ) scheme
    dZeQEDEW = dZeZQEDEW
  else
    call ol_error(2, 'in ew_renormalisation: ew_renorm_scheme= ' // to_string(ew_renorm_scheme) // ' not supported!')
    call ol_fatal()
  end if



  if (cms_on /= 2) then ! CMS-I / on-shell: truncate imaginary parts
   ! light fermions
   dZeLEW = real(dZeLEW)
   dZeREW = real(dZeREW)
   dZnLEW = real(dZnLEW)
   dZnlLEW = real(dZnlLEW)
   dZuLEW = real(dZuLEW)
   dZuREW = real(dZuREW)
   dZdLEW = real(dZdLEW)
   dZdREW = real(dZdREW)
   ! diagonal WFRCs
   dZAAEW = real(dZAAEW)
   dZAAEWdreg = real(dZAAEWdreg)
   dZZZEW  = real(dZZZEW)
   dZWEW = real(dZWEW)
   dZHEW  = real(dZHEW)
   ! charge renormalisation
   dZeQEDEW = real(dZeQEDEW)
   dZe0QEDEWnreg = real(dZe0QEDEWnreg)
   dZeGmuQEDEW = real(dZeGmuQEDEW)
   dZeZQEDEW = real(dZeZQEDEW)
  end if


  dZtLEWcc = dZtLEW
  dZtREWcc = dZtREW
  dZbLEWcc = dZbLEW
  dZbREWcc = dZbREW
  dZuLEWcc = dZuLEW
  dZuREWcc = dZuREW
  dZdLEWcc = dZdLEW
  dZdREWcc = dZdREW
  dZeLEWcc = dZeLEW
  dZeREWcc = dZeREW
  dZlLEWcc = dZlLEW
  dZlREWcc = dZlREW
  dZnLEWcc = dZnLEW
  dZnlLEWcc = dZnlLEW



!pure pole contributions for debugging
!     dZuRPole=-gZu(1)**2*de1_UV - 4.*cONE*(de1_UV-de1_IR)/9.
!     dZuLPole=-(gZu(2)**2+1/(2*sw2))*de1_UV - 4.*cONE*(de1_UV-de1_IR)/9.
!
!     dZtRPole=- gZu(1)**2*de1_UV - 4.*cONE*de1_UV/9. - 2d0*(MT2/4d0/sw2/MW2)*de1_UV-8d0/9d0*de1_IR
!     dZtLPole=- (gZu(2)**2+1/(2*sw2))*de1_UV - 4*cONE*de1_UV/9. - (MT2/4d0/sw2/MW2)*de1_UV-8d0/9d0*de1_IR
!
!     dZAAEWPole = -32d0/3d0*(de1_UV-de1_IR)+3d0*de1_UV-16d0/9d0*de1_IR
!     dZAAEWPole = -32d0/3d0*(de1_UV)+3d0*de1_UV   !pure UV
!     dZZAEWPole = 4d0*cw/sw*de1_UV
!     dZAZEWPole = -1d0*(30d0*cw2+1)/3d0/sw/cw*de1_UV & !Bosonic
!                   -8d0/3d0*de1_UV*(8d0*sw2-3d0)/sw/cw   !Fermionic
!     dZZZEWPole = de1_UV/6d0/sw2/cw2*(18d0*cw4+2d0*cw2-1d0) & !& !Bosonic
!                   -de1_UV*2d0/3d0*(6d0-12d0*sw2+16d0*sw**4)/cw2/sw2  !Fermionic
!     dZWEWPole = 10d0/3d0*de1_UV+1d0/12d0/sw2*(40*cw2-2)*de1_UV-2d0*de1_IR - 4d0/sw2*de1_UV
!
!     dZeQEDEWPole = -0.5d0*dZAAEWPole-0.5*sw/cw*dZZAEWPole     ! on-shell scheme
!
!     dZMZ2EWPole = -1d0/3d0/sw2/cw2*de1_UV*(MZ2*(9d0*cw**4+cw2-7d0/2d0)+MW2*(12*cw2-6d0)) & !Bosonic
!                  + 2d0/3d0*de1_UV*(MZ2*(6d0-12d0*sw2+16d0*sw**4)/cw2/sw2-9d0/4d0*MT2/sw2/cw2)  !Fermionic
!     dZMW2EWPole = -10d0/3*MW2*de1_UV - 1d0/12d0/sw2*((40d0*cw2+38d0)*MW2-(16d0*cw2+12d0)*MZ2)*de1_UV & !Bosonic
!                   + 1d0/3d0/sw2*de1_UV*(12*MW2 - 9d0/2d0*MT2) !Fermionic
!     dcwEWPole     = cw/2.*(dZMW2EWPole/MW2-dZMZ2EWPole/MZ2)

! check explicit poles for debugging
!      print*, de1_UV, de1_IR, dZuREW-dZuRPole, dZuLEW-dZuLPole
!      print*, de1_UV, de1_IR, dZtREW-dZtRPole, dZtLEW-dZtLPole
!      print*, de1_UV, de1_IR, dZeQEDEW-dZeQEDEWPole
!      print*, de1_UV, de1_IR, dcwEW-dcwEWPole
!      print*, de1_UV, de1_IR, dZAAEW-dZAAEWPole
!      print*, de1_UV, de1_IR, dZZAEW-dZZAEWPole
!      print*, de1_UV, de1_IR, dZAZEW-dZAZEWPole
!      print*, de1_UV, de1_IR, dZZZEW-dZZZEWPole
!      print*, de1_UV, de1_IR, dZWEW-dZWEWPole
!      print*, de1_UV, de1_IR, dZMZ2EW-dZMZ2EWPole
!      print*, de1_UV, de1_IR, dZMW2EW-dZMW2EWPole


! print rs constans
  if (debug_ew_renorm .ge. 1) then
    debug_norm = alpha_QED/4/pi
    print*, "==================================="
    print*, "==             DEBUG             =="
    print*, "==  EW renormalization constants =="
    print*, "==================================="
    print*, "loop_parameters_status", loop_parameters_status
    print*, "de1_UV", de1_UV
    print*, "de1_IR", de1_IR
    print*, "ew_renorm_scheme: ", ew_renorm_scheme
    print*, "Tadpole ", real(dtEW*debug_norm/eQED), aimag(dtEW*debug_norm/eQED)
    print*, "WF_V"
    print*, "dZAAEW  " , real(dZAAEW*debug_norm), aimag(dZAAEW*debug_norm)
    print*, "dZAZEW  " , real(dZAZEW*debug_norm), aimag(dZAZEW*debug_norm)
    print*, "dZZAEW  " , real(dZZAEW*debug_norm), aimag(dZZAEW*debug_norm)
    print*, "dZZZEW  " , real(dZZZEW*debug_norm), aimag(dZZZEW*debug_norm)
    print*, "dZWEW   " , real(dZWEW*debug_norm), aimag(dZWEW*debug_norm)
    print*, "dZHEW   " , real(dZHEW*debug_norm), aimag(dZHEW*debug_norm)
    print*, "M_V"
    print*, "dZMW2EW " , real(dZMW2EW*debug_norm), aimag(dZMW2EW*debug_norm)
    print*, "dZMZ2EW " , real(dZMZ2EW*debug_norm), aimag(dZMZ2EW*debug_norm)
    print*, "dZMH2EW " , real(dZMH2EW*debug_norm), aimag(dZMH2EW*debug_norm)
    print*, "WF_F"
    print*, "dZuREW  " , real(dZuREW*debug_norm), aimag(dZuREW*debug_norm)
    print*, "dZuLEW  " , real(dZuLEW*debug_norm), aimag(dZuLEW*debug_norm)
    print*, "dZdREW  " , real(dZdREW*debug_norm), aimag(dZdREW*debug_norm)
    print*, "dZdLEW  " , real(dZdLEW*debug_norm), aimag(dZdLEW*debug_norm)
    print*, "dZbREW  " , real(dZbREW*debug_norm), aimag(dZbREW*debug_norm)
    print*, "dZbLEW  " , real(dZbLEW*debug_norm), aimag(dZbLEW*debug_norm)
    print*, "dZtREW  " , real(dZtREW*debug_norm), aimag(dZtREW*debug_norm)
    print*, "dZtLEW  " , real(dZtLEW*debug_norm), aimag(dZtLEW*debug_norm)
    print*, "dZeREW  " , real(dZeREW*debug_norm), aimag(dZeREW*debug_norm)
    print*, "dZeLEW  " , real(dZeLEW*debug_norm), aimag(dZeLEW*debug_norm)
    print*, "dZlREW  " , real(dZlREW*debug_norm), aimag(dZlREW*debug_norm)
    print*, "dZlLEW  " , real(dZlLEW*debug_norm), aimag(dZlLEW*debug_norm)
    print*, "dZnLEW  " , real(dZnLEW*debug_norm), aimag(dZnLEW*debug_norm)
    print*, "dZnlLEW " , real(dZnlLEW*debug_norm), aimag(dZnlLEW*debug_norm)
    print*, "M_F"
    print*, "dZMBEW  " , real(dZMBEW*debug_norm), aimag(dZMBEW*debug_norm)
    print*, "dZMTEW  " , real(dZMTEW*debug_norm), aimag(dZMTEW*debug_norm)
    print*, "dZMLEW  " , real(dZMLEW*debug_norm), aimag(dZMLEW*debug_norm)
    print*, "C"
    print*, "dZeQEDEW" , real(dZeQEDEW*debug_norm), aimag(dZeQEDEW*debug_norm)
    print*, "dcwEW   " , real(dcwEW*debug_norm), aimag(dcwEW*debug_norm)
    print*, "dcwEW/cw", real(dcwEW/cw*debug_norm), aimag(dcwEW/cw*debug_norm)
    print*, "dswEW   " , real(dswEW*debug_norm), aimag(dswEW*debug_norm)
    print*, "dswEW/sw", real( dswEW/sw*debug_norm), aimag(dswEW/sw*debug_norm)
    print*, "==================================="
   end if
!   if (debug_ew_renorm .ge. 2) then
!     debug_norm = alpha_QED/4/pi
!     print*, "==================================="
!     print*, "==             DEBUG             =="
!     print*, "==    EW renorm. self energies   =="
!     print*, "==================================="
!     print*, "SitL" , SitL*debug_norm
!     print*, "SitR" , SitR*debug_norm
!     print*, "SitS" , SitS*debug_norm
!     print*, "dSitL" , dSitL*debug_norm
!     print*, "dSitR" , dSitR*debug_norm
!     print*, "dSitS" , dSitS*debug_norm
!     print*, "==================================="
!   end if

! calculate counterterms & R2

  ! some abbreviations
  ! Sum of squared quark masses
  sumMQ2 = MU2 + MD2 + MS2 + MC2
  sumMQ2Q2 = 4._/**/REALKIND/9.*MU2 + 1._/**/REALKIND/9.*MD2 + 4._/**/REALKIND/9.*MS2 + 1._/**/REALKIND/9.*MC2
  sumMQ2QI = 2._/**/REALKIND/6.*MU2 + 1._/**/REALKIND/6.*MD2 + 2._/**/REALKIND/6.*MS2 + 1._/**/REALKIND/6.*MC2
  sumMQ2QUD = 2._/**/REALKIND/3.*MD2+1._/**/REALKIND/3.*MU2 + 2._/**/REALKIND/3.*MS2+1._/**/REALKIND/3.*MC2   ! Sum(Qu*mD2 - Qd*mU2)
  sumMQ4 = MU2**2 + MD2**2 + MS2**2 + MC2**2
  sumMUD = MU*MD + MC*MS
  sumMUD2 = (MU2+MD2)**2 + (MC2+MS2)**2
  sumMU2 = MU2 + MC2
  sumMD2 = MD2 + MD2
  sumMU4 = MU2**2 + MC2**2
  sumMD4 = MD2**2 + MD2**2
  if (nf > 4) then
    sumMQ2 = sumMQ2 + MB2
    sumMQ2Q2 = sumMQ2Q2 + 1._/**/REALKIND/9.*MB2
    sumMQ2QI = sumMQ2QI + 1._/**/REALKIND/6.*MB2
    sumMQ2QUD = sumMQ2QUD + 2._/**/REALKIND/3.*MB2
    sumMQ4 = sumMQ4 + MB2**2
    sumMD2 = sumMD2 + MB2
    sumMD4 = sumMD4 + MB2**2
  end if
  if (nf > 5) then
    sumMQ2 = sumMQ2 + MT2
    sumMQ2Q2 = sumMQ2Q2 + 4._/**/REALKIND/9.*MT2
    sumMQ2QI = sumMQ2QI + 2._/**/REALKIND/6.*MT2
    sumMQ2QUD = sumMQ2QUD + 1._/**/REALKIND/3.*MT2
    sumMQ4 = sumMQ4 + MT2**2
    sumMUD = sumMUD + MT*MB
    sumMUD2 = sumMUD2 + (MT2+MB2)**2
    sumMU2 = sumMU2 + MT2
    sumMU4 = sumMU4 + MT2**2
  end if
!  sumQ = 1.
!  sumQI = 3._/**/REALKIND/2.
!  sumQ2 = 5._/**/REALKIND/3.
!  sumI2 = 3._/**/REALKIND/2.

  sumML2 = ME2 + MM2 + ML2
  sumML4 = ME2**2 + MM2**2 + ML2**2
  sumML2Q2 = sumML2
  sumML2QI = sumML2/2.


  ! firt set all counterterms to zero
    EWctWW = 0
    EWctZZ = 0
    EWctAZ = 0.
    EWctAA = 0.
    EWctHH = 0.
    EWctXX = 0.
    EWctPP = 0.
    EWctXA = 0.
    EWctXZ = 0.
    EWctPW = 0.
    EWctuu   = 0.
    EWctdd   = 0.
    EWcttt   = 0.
    EWctbb   = 0.
    EWctee   = 0.
    EWctLL   = 0.
    EWctnn   = 0.
    EWctWWWW = 0.
    EWctWWZZ = 0.
    EWctWWAZ = 0.
    EWctWWAA = 0.
    EWctR2AAAA = 0.
    EWctR2AAAZ = 0.
    EWctR2AAZZ = 0.
    EWctR2AZZZ = 0.
    EWctR2ZZZZ = 0.
    EWctAWW = 0.
    EWctZWW = 0.
    EWctSSSS1 = 0.
    EWctSSSS2 = 0.
    EWctSSSS3 = 0.
    EWctHHH = 0.
    EWctHXX = 0.
    EWctHPP = 0.
    EWctWWXX = 0.
    EWctWWHH = 0.
    EWctZZPP = 0.
    EWctZAPP = 0.
    EWctAAPP = 0.
    EWctZZHH = 0.
    EWctZZXX = 0.
    EWctZAHH = 0.
    EWctWZPH = 0.
    EWctWAPH = 0.
    EWctWZPX = 0.
    EWctWAPX = 0.
    EWctAXH = 0.
    EWctZXH = 0.
    EWctAPP = 0.
    EWctZPP = 0.
    EWctWPH = 0.
    EWctWPX = 0.
    EWctHWW = 0.
    EWctHZZ = 0.
    EWctHZA = 0.
    EWctPWZ = 0.
    EWctPWA = 0.
    EWctHAA = 0.
    EWctAuu  = 0.
    EWctAdd  =  0.
    EWctAtt  =  0.
    EWctAbb  = 0.
    EWctAee  = 0.
    EWctALL  = 0.
    EWctAnn  = 0.
    dgZu = 0.
    dgZd = 0.
    dgZl = 0.
    dgZn = 0.
    EWctVuu  = 0.
    EWctVdd  = 0.
    EWctVtt  = 0.
    EWctVbb  = 0.
    EWctVee  = 0.
    EWctVll  = 0.
    EWctVnn  = 0.
    EWctVnlnl  = 0.
    EWctVdu  = 0.
    EWctVbt  = 0.
    EWctVen  = 0.
    EWctVLn  = 0.
    EWctVud  = 0.
    EWctVtb  = 0.
    EWctVne  = 0.
    EWctVnL  = 0.
    EWctGuu  = 0.
    EWctGdd  = 0.
    EWctGtt  = 0.
    EWctGbb  = 0.
    EWctHtt = 0.
    EWctHbb = 0.
    EWctHLL = 0.
    EWctXtt = 0.
    EWctXbb = 0.
    EWctXLL = 0.
    EWctPdu = 0.
    EWctPud = 0.
    EWctPtb = 0.
    EWctPbt = 0.
    EWctPnL = 0.
    EWctPLn = 0.
    EWctAUWUW = 0.
    EWctZUWUW = 0.
    EWctWUWUZ = 0.
    EWctWUZUW = 0.
    EWctWUWUA = 0.
    EWctWUAUW = 0.
    EWctHUZUZ = 0.
    EWctHUWUW   = 0.
    EWctXUWUW   = 0.
    EWctPUZUW  = 0.
    EWctPUWUZ  = 0.
    EWctPUWUA  = 0.



  ! set UV counterterms
  if (CT_is_on /= 0) then
    ! Only UV counterterms

    ! VV Vector propagators
    EWctWW = [ dZWEW  , -1.*MW2*dZWEW  - dZMW2EW , ZERO ]
    EWctZZ = [ dZZZEW , -1.*MZ2*dZZZEW - dZMZ2EW , ZERO ]
    EWctAZ = [ 0.5*dZAZEW+0.5*dZZAEW , -0.5*MZ2*dZZAEW  , ZERO ]
    EWctAA = [ dZAAEW , ZERO , ZERO ]

    ! SS scalar propagators
    EWctHH = [ dZHEW , MH2*dZHEW + dZMH2EW ]
    EWctXX = [ ZERO, - dtEW/(2.*sw) + dZMZ2EW ]
    EWctPP = [ ZERO, - dtEW/(2.*sw) + dZMW2EW ]

!    EWctXA = MZ/4.*dZZAEW
!    EWctXZ = MZ/4.*(dZZZEW + dZMZ2EW/MZ2)
!    EWctPW = (dZWEW + dZMW2EW/MW2)

    ! FF fermionic propagators
    EWctuu   = [ 0.5*(dZuREW+dZuREWcc), 0.5*(dZuLEW+dZuLEWcc), ZERO , ZERO ]
    EWctdd   = [ 0.5*(dZdREW+dZdREWcc), 0.5*(dZdLEW+dZdLEWcc), ZERO , ZERO ]
    EWcttt   = [ 0.5*(dZtREW+dZtREWcc), 0.5*(dZtLEW+dZtLEWcc), MT/2.*(dZtREW+dZtLEWcc)+dZMTEW , MT/2.*(dZtLEW+dZtREWcc)+dZMTEW ]
    EWctbb   = [ 0.5*(dZbREW+dZbREWcc), 0.5*(dZbLEW+dZbLEWcc), MB/2.*(dZbREW+dZbLEWcc)+dZMBEW , MB/2.*(dZbLEW+dZbREWcc)+dZMBEW ]
    EWctee   = [ 0.5*(dZeREW+dZeREWcc), 0.5*(dZeLEW+dZeLEWcc), ZERO , ZERO ]
    EWctLL   = [ 0.5*(dZLREW+dZLREWcc), 0.5*(dZLLEW+dZLLEWcc), ML/2.*(dZLREW+dZLLEWcc)+dZMLEW , ML/2.*(dZLLEW+dZLREWcc)+dZMLEW ]
    EWctnn   = [ ZERO , 0.5*(dZnLEW+dZnLEWcc) , ZERO , ZERO ]
    EWctnlnl = [ ZERO , 0.5*(dZnlLEW+dZnlLEWcc) , ZERO , ZERO ]

    !VVVV
    EWctWWWW = 1./sw2*(2*dZeQEDEW - 2.*dswEW/sw + 2.*dZWEW) * [2,-1]
    EWctWWZZ = (-cw2/sw2*(2*dZeQEDEW - 2./cw2*dswEW/sw + dZWEW + dZZZEW) + cw/sw*dZAZEW)*[2,-1]
    EWctWWAZ = (cw/sw*(2*dZeQEDEW - 1./cw2*dswEW/sw + dZWEW + 0.5*dZZZEW + 0.5*dZAAEW) &
    - 0.5*dZAZEW- 0.5*cw2/sw2*dZZAEW) * [2,-1]
    EWctWWAA = (-1.*(2*dZeQEDEW + dZWEW + dZAAEW) + cw/sw*dZZAEW) * [2,-1]

    !VVV
    EWctAWW = (dZeQEDEW + dZWEW + 0.5*dZAAEW - 0.5*cw/sw*dZZAEW)
    EWctZWW = -cw/sw * (dZeQEDEW - dswEW/cw2/sw + dZWEW + 0.5*dZZZEW) + 0.5*dZAZEW

    !SSSS
    EWctSSSS1 = 2.*dZeQEDEW - 2.*dswEW/sw + dZMH2EW/MH2 + 1./(2.*sw) * dtEW/MH2 - dZMW2EW/MW2
    EWctSSSS2 = EWctSSSS1 + 2.*dZHEW
    EWctSSSS3 = EWctSSSS1 + dZHEW

    EWctHHHH = 3*EWctSSSS2
    EWctHHXX = EWctSSSS3
    EWctHHPP = EWctSSSS3
    EWctXXXX = 3*EWctSSSS1
    EWctXXPP = EWctSSSS1
    EWctPPPP = 2*EWctSSSS1

    !SSS
    EWctHHH = dZeQEDEW - dswEW/sw + dZMH2EW/MH2 + 0.5/sw * dtEW/MH2 - 0.5*dZMW2EW/MW2 + 1.5*dZHEW
    EWctHXX = EWctHHH - dZHEW
    EWctHPP = EWctHHH - dZHEW

    !VVSS
    EWctWWXX = 0.5/sw2*(2.*dZeQEDEW - 2.*dswEW/sw + dZWEW)
    EWctWWHH = EWctWWXX + 0.5/sw2*dZHEW
    EWctZZPP = (sw2-cw2)**2/2./sw2/cw2*(2*dZeQEDEW + 2./(sw2-cw2)/cw2*dswEW/sw + dZZZEW) + (sw2-cw2)/sw/cw*dZAZEW
    EWctZAPP = (sw2-cw2)/sw/cw*(2*dZeQEDEW + 1./(sw2-cw2)/cw2*dswEW/sw + 0.5*dZZZEW + 0.5*dZAAEW) + &
                (sw2-cw2)**2/4./sw2/cw2*dZZAEW + dZAZEW
    EWctAAPP = 2.*(2.*dZeQEDEW + dZAAEW) + (sw2-cw2)/sw/cw*dZZAEW
    EWctZZXX = 0.5/sw2/cw2*(2.*dZeQEDEW + 2.*(sw2-cw2)/cw2*dswEW/sw + dZZZEW )
    EWctZZHH = EWctZZXX + 0.5/sw2/cw2*dZHEW
    EWctZAHH = 0.25/sw2/cw2*dZZAEW
    EWctWZPH = -0.5/cw*(2.*dZeQEDEW - dcwEW/cw + 0.5*dZWEW + 0.5*dZHEW + 0.5*dZZZEW) - 0.25/sw*dZAZEW
    EWctWAPH = -0.5/sw*(2.*dZeQEDEW - dswEW/sw + 0.5*dZWEW + 0.5*dZHEW + 0.5*dZAAEW) - 0.25/cw*dZZAEW
    EWctWZPX = -1*CI/2./cw* ( 2.*dZeQEDEW - dcwEW/cw + 0.5*dZWEW + 0.5*dZZZEW ) - CI*0.25/sw*dZAZEW
    EWctWAPX = -1*CI/2./sw* ( 2.*dZeQEDEW - dswEW/sw + 0.5*dZWEW + 0.5*dZAAEW ) - CI*0.25/cw*dZZAEW

    !VSS
    EWctAXH = -0.25*CI/cw/sw*dZZAEW
    EWctZXH = -0.5*CI/cw/sw* ( dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw + 0.5*dZHEW + 0.5*dZZZEW )
    EWctAPP = -1.*( dZeQEDEW + 0.5*dZAAEW + 0.25*(sw2-cw2)/sw/cw*dZZAEW )
    EWctZPP = -0.5*(sw2-cw2)/sw/cw*( dZeQEDEW + dswEW/sw/cw2/(sw2-cw2) + 0.5*dZZZEW) - 0.5*dZAZEW
    EWctWPH = -0.5/sw*( dZeQEDEW - dswEW/sw + 0.5*dZWEW + 0.5*dZHEW )
    EWctWPX = -0.5*CI/sw*( dZeQEDEW - dswEW/sw + 0.5*dZWEW )

    !SVV
    EWctHWW = MW/sw * ( dZeQEDEW - dswEW/sw + 0.5*dZMW2EW/MW2 + 0.5*dZHEW + dZWEW )
    EWctHZZ = MW/sw/cw2 * ( dZeQEDEW + (2.*sw2-cw2)/cw2*dswEW/sw + 0.5*dZMW2EW/MW2 + 0.5*dZHEW + &
              dZZZEW )
    EWctHZA = 0.5*MW/sw/cw2 * dZZAEW
    EWctPWZ = -1.*MW*sw/cw * ( dZeQEDEW + dswEW/cw2/sw + 0.5*dZMW2EW/MW2 + 0.5*dZWEW + 0.5*dZZZEW )-&
              0.5*MW*dZAZEW
    EWctPWA = -1.*MW * ( dZeQEDEW + 0.5*dZMW2EW/MW2 + 0.5*dZWEW + 0.5*dZAAEW ) - 0.5*MW*sw/cw*dZZAEW

    !VFF
    ! Gff
    EWctGuu  = [ dZuREW , dZuLEW ]
    EWctGdd  = [ dZdREW , dZdLEW ]
    EWctGtt  = [ dZtREW , dZtLEW ]
    EWctGbb  = [ dZbREW , dZbLEW ]
    ! Aff
    EWctAuu  = [ -2._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZuREW+dZuREWcc))  + gZu(1)/2.*dZZAEW , &
                 -2._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZuLEW+dZuLEWcc))  + gZu(2)/2.*dZZAEW ]
    EWctAdd  = [  1._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZdREW+dZdREWcc))  + gZd(1)/2.*dZZAEW , &
                  1._/**/REALKIND/3. * (dZeQEDEW+  0.5*dZAAEW + 0.5*(dZdLEW+dZdLEWcc))  + gZd(2)/2.*dZZAEW ]
    EWctAtt  = [ -2._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZtREW+dZtREWcc))  + gZu(1)/2.*dZZAEW , &
                 -2._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZtLEW+dZtLEWcc))  + gZu(2)/2.*dZZAEW ]
    EWctAbb  = [  1._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZbREW+dZbREWcc))  + gZd(1)/2.*dZZAEW , &
                  1._/**/REALKIND/3. * (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZbLEW+dZbLEWcc))  + gZd(2)/2.*dZZAEW ]
    EWctAee  = [  (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZeREW+dZeREWcc))  + gZl(1)/2.*dZZAEW , &
                  (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZeLEW+dZeLEWcc))  + gZl(2)/2.*dZZAEW ]
    EWctALL  = [  (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZlREW+dZlREWcc))  + gZl(1)/2.*dZZAEW ,  &
                  (dZeQEDEW + 0.5*dZAAEW + 0.5*(dZlLEW+dZlLEWcc))  + gZl(2)/2.*dZZAEW ]

    EWctAnn  = [ ZERO , gZn(2)*0.5*dZZAEW ]

    ! Zff
    dgZu =  [ gZu(1) * (dZeQEDEW + dswEW/cw2/sw) , 0.5*gZLH * (dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw) + &
                gZu(1) * (dZeQEDEW + dswEW/cw2/sw) ]
    dgZd =  [ gZd(1) * (dZeQEDEW + dswEW/cw2/sw) , -0.5*gZLH * (dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw) + &
                gZd(1) * (dZeQEDEW + dswEW/cw2/sw) ]
    dgZl =  [ gZl(1) * (dZeQEDEW + dswEW/cw2/sw) , -0.5*gZLH * (dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw) + &
                gZl(1) * (dZeQEDEW + dswEW/cw2/sw) ]
    dgZn =  [ ZERO , 0.5*gZLH * (dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw) ]
    EWctVuu  = [ dgZu(1) + gZu(1) * (0.5*dZZZEW + 0.5*(dZuREW+dZuREWcc)) - 2._/**/REALKIND/6.*dZAZEW , &
                 dgZu(2) + gZu(2) * (0.5*dZZZEW + 0.5*(dZuLEW+dZuLEWcc)) - 2._/**/REALKIND/6.*dZAZEW ]
    EWctVdd  = [ dgZd(1) + gZd(1) * (0.5*dZZZEW + 0.5*(dZdREW+dZdREWcc)) + 1._/**/REALKIND/6.*dZAZEW , &
                 dgZd(2) + gZd(2) * (0.5*dZZZEW + 0.5*(dZdLEW+dZdLEWcc)) + 1._/**/REALKIND/6.*dZAZEW ]
    EWctVtt  = [ dgZu(1) + gZu(1) * (0.5*dZZZEW + 0.5*(dZtREW+dZtREWcc)) - 2._/**/REALKIND/6.*dZAZEW , &
                 dgZu(2) + gZu(2) * (0.5*dZZZEW + 0.5*(dZtLEW+dZtLEWcc)) - 2._/**/REALKIND/6.*dZAZEW ]
    EWctVbb  = [ dgZd(1) + gZd(1) * (0.5*dZZZEW + 0.5*(dZbREW+dZbREWcc)) + 1._/**/REALKIND/6.*dZAZEW , &
                 dgZd(2) + gZd(2) * (0.5*dZZZEW + 0.5*(dZbLEW+dZbLEWcc)) + 1._/**/REALKIND/6.*dZAZEW ]
    EWctVee  = [ dgZl(1) + gZl(1) * (0.5*dZZZEW + 0.5*(dZeREW+dZeREWcc)) + 0.5*dZAZEW , &
                 dgZl(2) + gZl(2) * (0.5*dZZZEW + 0.5*(dZeLEW+dZeLEWcc)) + 0.5*dZAZEW ]
    EWctVll  = [ dgZl(1) + gZl(1) * (0.5*dZZZEW + 0.5*(dZlREW+dZlREWcc)) + 0.5*dZAZEW , &
                 dgZl(2) + gZl(2) * (0.5*dZZZEW + 0.5*(dZlLEW+dZlLEWcc)) + 0.5*dZAZEW ]
    EWctVnn  = [ ZERO , dgZn(2) + gZn(2)* (0.5*dZZZEW + 0.5*(dZnLEW+dZnLEWcc)) ]
    EWctVnlnl = [ ZERO , dgZn(2) + gZn(2)* (0.5*dZZZEW + 0.5*(dZnlLEW+dZnlLEWcc)) ]


    ! Wff
    EWctVdu  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZuLEW + dZdLEWcc)
    EWctVud  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZuLEWcc + dZdLEW)
    EWctVbt  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZtLEW + dZbLEWcc)
    EWctVtb  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZtLEWcc + dZbLEW)
    EWctVen  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZeLEW + dZnLEWcc)
    EWctVne  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZeLEWcc + dZnLEW)

    EWctVln  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZlLEW + dZnlLEWcc)
    EWctVnl  = dZeQEDEW-dswEW/sw+0.5*dZWEW+0.5*(dZlLEWcc + dZnlLEW)

    !SFF
    EWctHtt = -0.5/sw*( dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 + 0.5*dZHEW) &
              -0.5/sw*0.5*[dZtREW+dZtLEWcc,dZtREWcc+dZtLEW ]
    if (MB /= 0._/**/REALKIND) then
      EWctHbb = -0.5/sw*( dZeQEDEW - dswEW/sw + dZMBEW/MB - 0.5*dZMW2EW/MW2 + 0.5*dZHEW) &
                -0.5/sw*0.5*[dZbREW+dZbLEWcc,dZbREWcc+dZbLEW ]
    else
      EWctHbb = 0.
    end if

    if (ML /= 0._/**/REALKIND) then
      EWctHll = -0.5/sw*( dZeQEDEW - dswEW/sw + dZMLEW/ML - 0.5*dZMW2EW/MW2 + 0.5*dZHEW) &
                -0.5/sw*0.5*[dZlREW+dZlLEWcc,dZlREWcc+dZlLEW ]
    else
      EWctHll = 0.
    end if

    EWctXtt = CI/sw*0.5*( dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 )*[1,-1] + &
              CI/sw*0.25*[(dZtREW+dZtLEWcc),-1.*(dZtREWcc+dZtLEW)]
    if (MB /= 0._/**/REALKIND) then
      EWctXbb = -CI/sw*0.5*( dZeQEDEW - dswEW/sw + dZMBEW/MB - 0.5*dZMW2EW/MW2 )*[1,-1] &
                -CI/sw*0.25*[(dZbREW+dZbLEWcc),-1.*(dZbREWcc+dZbLEW)]
    else
      EWctXbb = 0.
    end if

    if (ML /= 0._/**/REALKIND) then
      EWctXll = -CI/sw*0.5*( dZeQEDEW - dswEW/sw + dZMLEW/ML - 0.5*dZMW2EW/MW2 )*[1,-1] &
                -CI/sw*0.25*[(dZlREW+dZlLEWcc),-1.*(dZlREWcc+dZlLEW)]
    else
      EWctXll = 0.
    end if


    if (MB /= 0._/**/REALKIND) then
      EWctPtb = 1./sqrt2/sw*[ -1*MB*(dZeQEDEW - dswEW/sw + dZMBEW/MB - 0.5*dZMW2EW/MW2 + 0.5*(dZtLEWcc+dZbREW)), &
                MT*(dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 + 0.5*(dZtREWcc+dZbLEW)) ]
      EWctPbt = 1./sqrt2/sw*[  MT*(dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 + 0.5*(dZtREW+dZbLEWcc)), &
                -1*MB*(dZeQEDEW - dswEW/sw + dZMBEW/MB - 0.5*dZMW2EW/MW2 + 0.5*(dZtLEW+dZbREWcc)) ]
    else
      EWctPtb = 1./sqrt2/sw*[ ZERO , &
                MT*(dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 + 0.5*(dZtREWcc+dZbLEW)) ]
      EWctPbt = 1./sqrt2/sw*[  MT*(dZeQEDEW - dswEW/sw + dZMTEW/MT - 0.5*dZMW2EW/MW2 + 0.5*(dZtREW+dZbLEWcc)), &
                ZERO ]
    end if

    if (ML /= 0._/**/REALKIND) then
      EWctPnl = 1./sqrt2/sw*[ -1*ML*(dZeQEDEW - dswEW/sw + dZMLEW/ML - 0.5*dZMW2EW/MW2 + 0.5*(dZnlLEWcc+dZlREW)), ZERO]
      EWctPbt = 1./sqrt2/sw*[  ZERO, &
                -1*ML*(dZeQEDEW - dswEW/sw + dZMLEW/ML - 0.5*dZMW2EW/MW2 + 0.5*(dZnlLEW+dZlREWcc)) ]
    end if

    ! VUU
    EWctAUWUW = ( dZeQEDEW + 0.5*dZAAEW ) - 0.5*cw/sw*dZZAEW
    EWctZUWUW = -1*cw/sw*( dZeQEDEW - dswEW/cw2/sw + 0.5*dZZZEW ) + 0.5*dZAZEW
    EWctWUWUZ = cw/sw * ( dZeQEDEW - dswEW/cw2/sw + 0.5*dZWEW )
    EWctWUZUW = -1.*EWctWUWUZ
    EWctWUWUA = -1.*(dZeQEDEW + 0.5*dZWEW )
    EWctWUAUW = -1.*EWctWUWUA

    ! SUU
    EWctHUZUZ  =  ( dZeQEDEW + (2.*sw2-cw2)/cw2*dswEW/sw + 0.5*dZMW2EW/MW2 + 0.5*dZHEW )
    EWctHUWUW  =  ( dZeQEDEW - dswEW/sw + 0.5*dZMW2EW/MW2 + 0.5*dZHEW )
    EWctXUWUW  =  (dZeQEDEW - dswEW/sw + 0.5*dZMW2EW/MW2)
    EWctPUZUW  =  (dZeQEDEW + (sw2-cw2)/cw2*dswEW/sw + 0.5*dZMW2EW )
    EWctPUWUZ  =  ( dZeQEDEW + dswEW/(sw2-cw2)/cw2/sw + 0.5*dZMW2EW/MW2 )
    EWctPUWUA  =  ( dZeQEDEW + 0.5*dZMW2EW/MW2 )

  end if


  ! add R2
  if (R2_is_on /= 0) then

! twopoint
        EWctAA=EWctAA+[cONE*(-3. - (10._/**/REALKIND*nc)/9.), &
               2.*MW2 + 4*nc*sumMQ2Q2 + 4*sumML2Q2,cONE*(2._/**/REALKIND/3.)]

        EWctAZ=EWctAZ+[1./(2.*cw*sw) + cw/sw + nc/(2.*cw*sw) - (2.*sw)/cw -  &
                (10.*nc*sw)/(9.*cw),(-2.*cw*MW2)/sw &
                - (2.*nc*sumMQ2QI)/(cw*sw) - (2.*sumML2QI)/(cw*sw)  &
                + (4.*nc*sumMQ2Q2*sw)/cw + (4.*sumML2Q2*sw)/cw ,(-2.*cw)/(3.*sw)]

        EWctZZ=EWctZZ+[1./cw2 + nc/cw2 - 1/(2.*cw2*sw2) - cw2/sw2 -  &
                nc/(2.*cw2*sw2) - (2.*sw2)/cw2 - (10.*nc*sw2)/(9.*cw2), &
               (-4.*nc*sumMQ2QI)/cw2 + (-4.*sumML2QI)/cw2 + (2.*cw2*MW2)/sw2 &
               + (nc*sumMQ2)/(2.*cw2*sw2) + sumML2/(2.*cw2*sw2) &
               + (4.*nc*sumMQ2Q2*sw2)/cw2 + (4.*sumML2Q2*sw2)/cw2 , &
               (2.*cw2)/(3.*sw2)]

        EWctWW=EWctWW+[-(3 + nc)/(2*sw2), &
                      (2*MW2)/sw2 + (nc*sumMQ2)/(2.*sw2) + sumML2/(2.*sw2),2/(3.*sw2)]

        EWctHH=EWctHH+[-1/(12.*sw2) - 1/(24.*cw2*sw2)  &
               - (nc*sumMQ2)/(6.*MW2*sw2) - sumML2/(6.*MW2*sw2), &
               (5.*MW2)/(2.*sw2) + (11.*MW2)/(8.*cw4*sw2) - MZ2/(8.*cw2*sw2)  &
               - (nc*sumMQ4)/(MW2*sw2) - sumML4/(MW2*sw2)]

        EWctXX=EWctXX+[-1./(12.*sw2) - 1./(24.*cw2*sw2)  &
               - (nc*sumMQ2)/(6.*MW2*sw2) - sumML2/(6.*MW2*sw2), &
               - MH2/(8.*cw2*sw2) + MW2/(2.*sw2) + (3.*MW2)/(8.*cw4*sw2)  &
               - (nc*sumMQ4)/(MW2*sw2) - sumML4/(MW2*sw2)]

        EWctPP=EWctPP+[-1/(12*sw2) - 1/(24*cw2*sw2) &
               - (sumMQ2*nc)/(6*MW2*sw2) - sumML2/(6*MW2*sw2), &
               - MH2/(8*sw2) + MW2/(4*sw2) + (3*MW2)/(8*cw4*sw2) + (3*MW2)/(8*cw2*sw2) &
               - MZ2/(8*sw2) - (nc*sumMUD2)/(2*MW2*sw2) - sumML4/(2*MW2*sw2)]

        EWctdd=EWctdd+[1/(9.*cw2),(9. + 18.*cw2 - 8.*sw2)/(36.*cw2*sw2),ZERO,ZERO]*(-1.)

        EWctbb=EWctbb+[1/(9.*cw2),(9 + 18*cw2 - 8*sw2)/(36.*cw2*sw2), &
               -MB/(9.*cw2),-MB/(9.*cw2)]*(-1)

        EWctuu=EWctuu+[4./(9.*cw2),(9. + 18.*cw2 - 8.*sw2)/(36.*cw2*sw2),ZERO,ZERO]*(-1.)

        EWcttt=EWcttt+[4./(9.*cw2),(9. + 18.*cw2 - 8.*sw2)/(36.*cw2*sw2), &
               (2.*MT)/(9.*cw2),(2.*MT)/(9.*cw2)]*(-1.)

        EWctee=EWctee+[1./cw2,(1. + 2*cw2)/(4.*cw2*sw2),ZERO,ZERO]*(-1.)

        EWctll=EWctll+[1/cw2,(1 + 2*cw2)/(4.*cw2*sw2),(-4 + 5/cw2)*ML, &
                (-4 + 5/cw2)*ML]*(-1.)

        EWctnn=EWctnn+[ZERO,(1. + 2.*cw2)/(4.*cw2*sw2),ZERO,ZERO]*(-1.)

        EWctnlnl=EWctnlnl+[ZERO,(1. + 2.*cw2)/(4.*cw2*sw2),ZERO,ZERO]*(-1.)


! three-point
        EWctHtt=EWctHtt+[(2.*cw2*(9.*MB2 + MW2*(9. + 64.*sw2)) +  &
                    MW2*(9. - 96.*sw2 + 128.*sw4))/(144.*cw2*MW2*sw3), &
                    (2.*cw2*(9.*MB2 + MW2*(9. + 64.*sw2)) +  &
                    MW2*(9. - 96.*sw2 + 128.*sw4))/(144.*cw2*MW2*sw3)]

        EWctHbb=EWctHbb+[(2.*cw2*(9.*MT2 + MW2*(9. + 16.*sw2)) +  &
                    MW2*(9. - 48.*sw2 + 32.*sw4))/(144.*cw2*MW2*sw3), &
                    (2.*cw2*(9.*MT2 + MW2*(9. + 16.*sw2)) +  &
                    MW2*(9. - 48.*sw2 + 32.*sw4))/(144.*cw2*MW2*sw3)]

        EWctHll=EWctHll+[(ML*(1 - 16*sw2 + cw2*(2 + 32*sw2) + 32*sw4))/(16.*cw2*sw3), &
                 (ML*(1 - 16*sw2 + cw2*(2 + 32*sw2) + 32*sw4))/(16.*cw2*sw3)]

        EWctXtt=EWctXtt+[-(CI*(2*cw2*(9*MB2 + MW2*(9 + 64*sw2)) +  &
                    MW2*(9 - 96*sw2 + 128*sw4)))/(144*cw2*MW2*sw3), &
               (CI*(2*cw2*(9*MB2 + MW2*(9 + 64*sw2)) +  &
                     MW2*(9 - 96*sw2 + 128*sw4)))/(144*cw2*MW2*sw3)]

        EWctXbb=EWctXbb+[(CI* &
                   (2*cw2*(9*MT2 + MW2*(9 + 16*sw2)) + MW2*(9 - 48*sw2 + 32*sw4)) &
                   )/(144.*cw2*MW2*sw3), &
               -(CI*(2*cw2*(9*MT2 + MW2*(9 + 16*sw2)) +  &
                    MW2*(9 - 48*sw2 + 32*sw4)))/(144.*cw2*MW2*sw3)]

        EWctXll=EWctXll+[-(CI*ML*(1 - 16*sw2 + cw2*(2 + 32*sw2) + 32*sw4))/ &
                 (16.*cw2*sw3),(CI*ML*(1 - 16*sw2 + cw2*(2 + 32*sw2) + 32*sw4))/ &
                 (16.*cw2*sw3)]

        EWctPtb=EWctPtb+[(MB*(cw2*(18*MT2 + MW2*(27 - 46*sw2)) +  &
                    MW2*(39 - 46*sw2)*sw2))/(72.*cw2*MW2*sqrt2*sw3), &
               (MT*(MW2*sw2*(-87 + 46*sw2) +  &
                    cw2*(-18*MB2 + MW2*(-27 + 46*sw2))))/(72.*cw2*MW2*sqrt2*sw3)]

        EWctPbt=EWctPbt+[(23*MT)/(36.*sqrt2*sw) -  &
                (29*MT)/(24.*cw2*sqrt2*sw) + (23*MT*sw)/(36.*cw2*sqrt2) -  &
                (3*MT)/(8.*sqrt2*sw3) - (MB2*MT)/(4.*MW2*sqrt2*sw3), &
               (-23*MB)/(36.*sqrt2*sw) + (13*MB)/(24.*cw2*sqrt2*sw) -  &
                (23*MB*sw)/(36.*cw2*sqrt2) + (3*MB)/(8.*sqrt2*sw3) +  &
                (MB*MT2)/(4.*MW2*sqrt2*sw3)]

        EWctPnl=EWctPnl+[(ML*(cw2*(3 + 2*sw2) + sw2*(15 + 2*sw2)))/ &
                 (8.*cw2*sqrt2*sw3),ZERO]

        EWctPln=EWctPln+[ZERO,ML/(4.*sqrt2*sw) + (15*ML)/(8.*cw2*sqrt2*sw) +  &
                 (ML*sw)/(4.*cw2*sqrt2) + (3*ML)/(8.*sqrt2*sw3)]

        EWctAuu=EWctAuu+[16/(27.*cw2),-8/(27.*cw2) + 2/(3.*sw2) + 1/(3.*cw2*sw2)]

        EWctAdd=EWctAdd+[-2/(27.*cw2),4/(27.*cw2) - 1/(3.*sw2) - 1/(6.*cw2*sw2)]

        EWctAtt=EWctAtt+[16/(27.*cw2) + MT2/(12.*MW2*sw2), &
               -8/(27.*cw2) + 2/(3.*sw2) + 1/(3.*cw2*sw2) - MB2/(12.*MW2*sw2) +  &
                MT2/(6.*MW2*sw2)]

        EWctAbb=EWctAbb+[-2/(27.*cw2) + MB2/(12.*MW2*sw2), &
               4/(27.*cw2) - 1/(3.*sw2) - 1/(6.*cw2*sw2) - MB2/(12.*MW2*sw2) +  &
                MT2/(6.*MW2*sw2)]

        EWctAee=EWctAee+[-2/cw2,-(1/sw2) - 1/(2.*cw2*sw2)]

        EWctAll=EWctAll+[-2/cw2 - ML2/(4.*MW2*sw2), &
                -(1/sw2) - 1/(2.*cw2*sw2) - ML2/(4.*MW2*sw2)]

        EWctVuu=EWctVuu+[(16*sw)/(27.*cw3), &
               -7/(9.*cw*sw) + 1/(cw3*sw) + (16*sw)/(27.*cw) - (4*sw)/(3.*cw3) +  &
                1/(2.*cw*sw3) - cw/sw3 - 1/(4.*cw3*sw3) + (16*sw3)/(27.*cw3)]

        EWctVdd=EWctVdd+[(-2*sw)/(27.*cw3), &
               7/(9.*cw*sw) - 1/(2.*cw3*sw) - (2*sw)/(27.*cw) + sw/(3.*cw3) -  &
                1/(2.*cw*sw3) + cw/sw3 + 1/(4.*cw3*sw3) - (2*sw3)/(27.*cw3)]

        EWctVtt=EWctVtt+[MT2/(12.*cw*MW2*sw) + (16*sw)/(27.*cw) +  &
                (16*sw3)/(27.*cw3),-7/(9.*cw*sw) + 1/(cw3*sw) -  &
                MB2/(12.*cw*MW2*sw) + MT2/(6.*cw*MW2*sw) + (16*sw)/(27.*cw) -  &
                (4*sw)/(3.*cw3) + 1/(2.*cw*sw3) - cw/sw3 - 1/(4.*cw3*sw3) +  &
                (16*sw3)/(27.*cw3)]

        EWctVbb=EWctVbb+[MB2/(12.*cw*MW2*sw) - (2*sw)/(27.*cw) -  &
                (2*sw3)/(27.*cw3),7/(9.*cw*sw) - 1/(2.*cw3*sw) -  &
                MB2/(12.*cw*MW2*sw) + MT2/(6.*cw*MW2*sw) - (2*sw)/(27.*cw) +  &
                sw/(3.*cw3) - 1/(2.*cw*sw3) + cw/sw3 + 1/(4.*cw3*sw3) -  &
                (2*sw3)/(27.*cw3)]

        EWctVee=EWctVee+[(-2*sw)/cw3, &
               1/(cw*sw) - 3/(2.*cw3*sw) - (2*sw)/cw + (3*sw)/cw3 -  &
                1/(2.*cw*sw3) + cw/sw3 + 1/(4.*cw3*sw3) - (2*sw3)/cw3]

        EWctVll=EWctVll+[(-2*sw)/cw3, &
               1/(cw*sw) - 3/(2.*cw3*sw) - ML2/(4*cw*MW2*sw) - (2*sw)/cw &
               + (3*sw)/cw3 - 1/(2.*cw*sw3) + cw/sw3 + 1/(4.*cw3*sw3) - (2*sw3)/cw3]

        EWctVnn=EWctVnn+[ZERO,-(1/(cw*sw)) + 1/(2.*cw*sw3) - cw/sw3 - 1/(4.*cw3*sw3)]

        EWctVnlnl=EWctVnlnl+[ZERO, &
            -(1/(cw*sw)) - ML2/(4.*cw*MW2*sw) + 1/(2.*cw*sw3) - cw/sw3 -1/(4.*cw3*sw3)]

        EWctVdu=EWctVdu-5/(9.*cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVbt=EWctVbt-5/(9.*cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVen=EWctVen-(1/cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVud=EWctVud-5/(9.*cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVtb=EWctVtb-5/(9.*cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVne=EWctVne-(1/cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctVnl=EWctVnl-(1/cw2) - 2/sw2 + 1/(2.*cw2*sw2)

        EWctHHH=EWctHHH-1/(4.*sw2) - 1/(8.*cw2*sw2) + (3*MW2)/(2.*MH2*sw2) +  &
               (3*MW2)/(4.*cw4*MH2*sw2) - (nc*sumMQ4)/(MH2*MW2*sw2) - sumML4/(MH2*MW2*sw2)

        EWctHXX=EWctHXX-1/(4.*sw2) - 1/(8.*cw2*sw2) + (3*MW2)/(2.*MH2*sw2) +  &
               (3*MW2)/(4.*cw4*MH2*sw2) - (nc*sumMQ4)/(MH2*MW2*sw2) - sumML4/(MH2*MW2*sw2)

        EWctHPP=EWctHPP-0.125 + (7*MW2)/(2.*MH2) - (2*MW2)/(cw2*MH2) -  &
               3/(8.*sw2) + (9*MW2)/(4.*MH2*sw2) - (nc*sumMQ4)/(MH2*MW2*sw2) - sumML4/(MH2*MW2*sw2) &
              - sw2/(8.*cw2) + (7*MW2*sw2)/(2.*cw2*MH2) + (3*MW2*sw2)/(4.*cw4*MH2)

        EWctAXH=EWctAXH+(-5*CI)/(12.*sw2)

        EWctZXH=EWctZXH+(CI*(MW2 + 22*cw4*MW2 &
              + 2*cw2*(4*nc*sumMQ2 + 4*sumML2 + MW2*sw2)))/(48.*cw3*MW2*sw3)

        EWctAPP=EWctAPP+1/(24.*cw2) + 13/(24.*sw2) + (nc*sumMQ2)/(3.*MW2*sw2) + sumML2/(3.*MW2*sw2)

        EWctZPP=EWctZPP+(nc*sumMQ2)/(3.*cw*MW2*sw) + sumML2/(3.*cw*MW2*sw) &
              + sw/(48.*cw) + 1/(24.*cw*sw3) - (25*cw)/(48.*sw3) + sw3/(48.*cw3)  &
              - (nc*sumMQ2)/(6.*cw*MW2*sw3) - sumML2/(6.*cw*MW2*sw3)

        EWctWPH=EWctWPH+1/(48.*cw2*sw) + 23/(48.*sw3) + (nc*sumMQ2)/(6.*MW2*sw3) + sumML2/(6.*MW2*sw3)

        EWctWPX=EWctWPX+CI/(48.*cw2*sw) + (23*CI)/(48.*sw3)  &
              + (CI*nc*sumMQ2)/(6.*MW2*sw3) + (CI*sumML2)/(6.*MW2*sw3)

        EWctHAA=EWctHAA-(MW/sw) - (2*nc*sumMQ2Q2)/(MW*sw) - (2*sumML2Q2)/(MW*sw)

        EWctHZA=EWctHZA+MW/(2.*cw) - (2*nc*sumMQ2Q2)/(cw*MW) - (2*sumML2Q2)/(cw*MW) &
               + (3*cw*MW)/(2.*sw2) + (nc*sumMQ2QI)/(cw*MW*sw2) + sumML2QI/(cw*MW*sw2)

        EWctHZZ=EWctHZZ+MW/sw + (2*nc*sumMQ2QI)/(cw2*MW*sw) + (2*sumML2QI)/(cw2*MW*sw)  &
               - (2*nc*sumMQ2Q2*sw)/(cw2*MW) - (2*sumML2Q2*sw)/(cw2*MW) - (2*MW)/sw3   &
               - (nc*sumMQ2)/(2.*cw2*MW*sw3) - sumML2/(2.*cw2*MW*sw3)

        EWctHWW=EWctHWW+(-2*MW)/sw3 - (nc*sumMQ2)/(2.*MW*sw3) - sumML2/(2.*MW*sw3)

        EWctPWA=EWctPWA+MW/(2.*sw2) + (nc*sumMQ2QUD)/(2.*MW*sw2) + sumML2/(2.*MW*sw2)

        EWctPWZ=EWctPWZ+MW/(2.*cw*sw) + (nc*sumMQ2QUD)/(2.*cw*MW*sw) + sumML2/(2.*cw*MW*sw)

        EWctAWW=EWctAWW - (17+6*nc)/(6.*sw2)

        EWctZWW=EWctZWW + cw*(17+6*nc)/(6.*sw3)

! EWctGff see V2
        EWctGuu=EWctGuu + [8/(9.*cw2),-4/(9.*cw2) + 1/sw2 + 1/(2.*cw2*sw2)]*(-1.)

        EWctGdd=EWctGdd + [2/(9.*cw2),-4/(9.*cw2) + 1/sw2 + 1/(2.*cw2*sw2)]*(-1.)

        EWctGtt=EWctGtt + [8/(9.*cw2) + MT2/(2.*MW2*sw2), &
                  -4/(9.*cw2) + 1/sw2 + 1/(2.*cw2*sw2) + MB2/(4.*MW2*sw2) +  &
                  MT2/(4.*MW2*sw2)]*(-1.)

        EWctGbb=EWctGbb + [2/(9.*cw2) + MB2/(2.*MW2*sw2), &
                  -4/(9.*cw2) + 1/sw2 + 1/(2.*cw2*sw2) + MB2/(4.*MW2*sw2) +  &
                  MT2/(4.*MW2*sw2)]*(-1.)

! four-point

        EWctHHHH=EWctHHHH-3/(2.*sw2) - 3/(4.*cw2*sw2) + (11*MW2)/(2.*MH2*sw2) +  &
               (11*MW2)/(4.*cw4*MH2*sw2) - (5*nc*sumMQ4)/(MH2*MW2*sw2) - (5*sumML4)/(MH2*MW2*sw2)

        EWctXXXX=EWctXXXX-3/(2.*sw2) - 3/(4.*cw2*sw2) + (11*MW2)/(2.*MH2*sw2) +  &
               (11*MW2)/(4.*cw4*MH2*sw2) - (5*nc*sumMQ4)/(MH2*MW2*sw2) - (5*sumML4)/(MH2*MW2*sw2)

        EWctHHXX=EWctHHXX-3/(2.*sw2) - 3/(4.*cw2*sw2) + (11*MW2)/(2.*MH2*sw2) +  &
               (11*MW2)/(4.*cw4*MH2*sw2) - (5*nc*sumMQ4)/(MH2*MW2*sw2) - (5*sumML4)/(MH2*MW2*sw2)

        EWctHHPP=EWctHHPP-3/(2.*sw2) - 3/(4.*cw2*sw2) + (11*MW2)/(2.*MH2*sw2) +  &
               (11*MW2)/(4.*cw4*MH2*sw2) - (5*nc*sumMQ4)/(MH2*MW2*sw2) - (5*sumML4)/(MH2*MW2*sw2)

        EWctXXPP=EWctXXPP-0.125 + (41*MW2)/(12.*MH2) - (2*MW2)/(cw2*MH2) -  &
               5/(8.*sw2) - 1/(8.*cw2*sw2) + (7*MW2)/(3.*MH2*sw2) + (5*MW2)/(12.*cw2*MH2*sw2) &
              - (5*nc*sumMQ4)/(3.*MH2*MW2*sw2) - (5*sumML4)/(3.*MH2*MW2*sw2) &
              - sw2/(8.*cw2) + (41*MW2*sw2)/(12.*cw2*MH2) +  &
               (11*MW2*sw2)/(12.*cw4*MH2)

        EWctPPPP=EWctPPPP-0.5 + (11*MW2)/(3.*MH2) - 3/(2.*sw2) +  &
               (11*MW2)/(2.*MH2*sw2) - (10*nc*sumMQ4)/(3.*MH2*MW2*sw2) - (10*sumML4)/(3.*MH2*MW2*sw2) &
              -  sw2/(2.*cw2) + (11*MW2*sw2)/(2.*MH2) +  (22*MW2*sw4)/(3.*cw2*MH2) + (11*MW2*sw6)/(6.*cw4*MH2)

        EWctAAHH=EWctAAHH+1/(12.*sw2) - (nc*sumMQ2Q2)/(MW2*sw2) - sumML2Q2/(MW2*sw2)

        EWctAAXX=EWctAAXX+1/(12.*sw2) - (nc*sumMQ2Q2)/(MW2*sw2) - sumML2Q2/(MW2*sw2)

        EWctZAHH=EWctZAHH+2./(3.*cw*sw) - (nc*sumMQ2Q2)/(cw*MW2*sw) - sumML2Q2/(cw*MW2*sw)&
              - 1/(2.*sw2) +  cw/(12.*sw3) + (nc*sumMQ2QI)/(2.*cw*MW2*sw3) + sumML2QI/(2.*cw*MW2*sw3)

        EWctZAXX=2/(3.*cw*sw) - (nc*sumMQ2Q2)/(cw*MW2*sw) - sumML2Q2/(cw*MW2*sw)  &
              - 1/(2.*sw2) + cw/(12.*sw3) + (nc*sumMQ2QI)/(2.*cw*MW2*sw3) + sumML2QI/(2.*cw*MW2*sw3)

        EWctZZHH=EWctZZHH-1/(24.*cw2) - (nc*sumMQ2Q2)/(cw2*MW2) - sumML2Q2/(cw2*MW2) &
               - 1/(8.*sw2) - 19/(24.*sw4) - 1/(48.*cw4*sw4)  &
               + (nc*sumMQ2QI)/(cw2*MW2*sw2) + sumML2QI/(cw2*MW2*sw2) &
               - (nc*sumMQ2)/(3.*cw2*MW2*sw4) - sumML2/(3.*cw2*MW2*sw4)

        EWctZZXX=EWctZZXX-1/(24.*cw2) - 1/(8.*sw2) &
               - (nc*sumMQ2Q2)/(cw2*MW2) - sumML2Q2/(cw2*MW2)  &
               + (nc*sumMQ2QI)/(cw2*MW2*sw2) + sumML2QI/(cw2*MW2*sw2) &
               - 19/(24.*sw4) - 1/(48.*cw4*sw4) &
               -  (nc*sumMQ2)/(3.*cw2*MW2*sw4) - sumML2/(3.*cw2*MW2*sw4)

        EWctWWHH=EWctWWHH-19/(24.*sw4) - 1/(48.*cw2*sw4) &
               - (nc*sumMQ2)/(3.*MW2*sw4) - sumML2/(3.*MW2*sw4)

        EWctWWXX=EWctWWXX-19/(24.*sw4) - 1/(48.*cw2*sw4) &
               - (nc*sumMQ2)/(3.*MW2*sw4) - sumML2/(3.*MW2*sw4)

        EWctWAPH=EWctWAPH+1/(2.*sw) + 1/(48.*cw2*sw) - cw/(2.*sw2)  &
               + 11/(48.*sw3) + (nc*sumMD2)/(4.*MW2*sw3) + sumML2/(12.*MW2*sw3) &
               + (nc*sumMU2)/(6.*MW2*sw3)

        EWctWAPX=EWctWAPX-CI/(2.*sw) - CI/(48.*cw2*sw) + (CI*cw)/(2.*sw2)  &
               - (11*CI)/(48.*sw3) - (CI*nc*sumMD2)/(4.*MW2*sw3) - (CI*sumML2)/(12.*MW2*sw3)  &
               - (CI*nc*sumMU2)/(6.*MW2*sw3)

        EWctWZPH=EWctWZPH+1/(24.*cw3) + 1/(2.*sw) + 5/(48.*cw*sw2) + cw/(2.*sw2) - 1/(48.*cw3*sw2) &
               + (nc*sumMD2)/(4.*cw*MW2*sw2) + sumML2/(12.*cw*MW2*sw2)   &
               + (nc*sumMU2)/(6.*cw*MW2*sw2) + 7/(48.*cw*sw4) - (7*cw)/(48.*sw4)

        EWctWZPX=EWctWZPX-CI/(24.*cw3) - CI/(2.*sw) - (5*CI)/(48.*cw*sw2)  &
              - (CI*cw)/(2.*sw2) + CI/(48.*cw3*sw2)  &
              - (CI*nc*sumMD2)/(4.*cw*MW2*sw2) - (CI*nc*sumML2)/(12.*cw*MW2*sw2) &
              - (CI*nc*sumMU2)/(6.*cw*MW2*sw2)  &
              - (7*CI)/(48.*cw*sw4) + (7*CI*cw)/(48.*sw4)

        EWctAAPP=EWctAAPP-1._/**/REALKIND/12. - 11/(6.*sw2) - sw2/(12.*cw2) &
              - (10*nc*sumMQ2)/(9.*MW2*sw2) - (4.*nc*sumML2)/(3.*MW2*sw2)

        EWctZAPP=EWctZAPP-0.25 - 1/(2.*cw*sw) &
              - (10*nc*sumMD2)/(9.*cw*MW2*sw) -  (4.*sumML2)/(3.*cw*MW2*sw) &
              - (10*nc*sumMU2)/(9.*cw*MW2*sw) + sw/(12.*cw) - 1/(4.*sw2)  &
              - sw2/(4.*cw2) + (7*cw)/(6.*sw3) &
              + (7*nc*sumMD2)/(12.*cw*MW2*sw3) +  (5*sumML2)/(12.*cw*MW2*sw3)&
              + (nc*sumMU2)/(2.*cw*MW2*sw3) + sw3/(12.*cw3)

        EWctZZPP=EWctZZPP+1._/**/REALKIND/48.  &
              - (10*nc*sumMD2)/(9.*cw**2*MW**2) - (4.*nc*sumMD2)/(3.*cw**2*MW**2) &
              - (10*nc*sumMU2)/(9.*cw**2*MW**2) - 37/(48.*sw**4)  &
              - 1/(24.*cw**2*sw**4) - (nc*sumMD2)/(3.*cw**2*MW**2*sw**4)  &
              - (nc*sumMU2)/(3.*cw**2*MW**2*sw**4) + 43/(24.*sw**2)   &
              + (7*nc*sumMD2)/(6.*cw**2*MW2*sw**2) + (5.*sumML2)/(6.*cw**2*MW2*sw**2)  &
              + (nc*sumMU2)/(cw**2*MW**2*sw**2) - sw**4/(48.*cw**4)

        EWctWWPP=EWctWWPP-4 - 1/(48.*cw2) + (4*cw)/sw + 119/(48.*sw2) - &
                15/(16.*sw4) - (nc*sumMQ2)/(3.*MW2*sw4) - sumML2/(3.*MW2*sw4)

        EWctR2AAAA=4._/**/REALKIND/3*(2._/**/REALKIND + (17*nc)/27._/**/REALKIND)

        EWctR2AAAZ=-1/(cw*sw) + (4*cw)/(3.*sw) - nc/(3.*cw*sw) + &
               (4*sw)/cw + (68*nc*sw)/(81.*cw)

        EWctR2AAZZ=4._/**/REALKIND/3 - 2._/**/REALKIND/cw2 - (2*nc)/(3.*cw2) - &
               4._/**/REALKIND/(3.*sw2) + 1._/**/REALKIND/(2.*cw2*sw2) &
                + (5*nc)/(18.*cw2*sw2) + (4*sw2)/cw2 + (68*nc*sw2)/(81.*cw2)

        EWctR2AZZZ=(-4*cw)/(3.*sw) + 3._/**/REALKIND/(2.*cw3*sw) + (5*nc)/(6.*cw3*sw) - &
               (3*sw)/cw3 - (nc*sw)/cw3 + (4*cw)/(3.*sw3) - 1._/**/REALKIND/(4.*cw3*sw3) - &
               nc/(4.*cw3*sw3) + (4*sw3)/cw3 + (68*nc*sw3)/(81.*cw3)

        EWctR2ZZZZ=-4._/**/REALKIND/3. + 3/cw4 + (5*nc)/(3.*cw4) + 8/(3.*sw2) - &
                   1/(cw4*sw2) - nc/(cw4*sw2) - (4*sw2)/cw4 - &
                   (4*nc*sw2)/(3.*cw4) - 4/(3.*sw4) + 1/(4.*cw4*sw4) + &
                   nc/(4.*cw4*sw4) + (4*sw4)/cw4 + (68*nc*sw4)/(81.*cw4)

        EWctWWAA=EWctWWAA+[23/(3*sw2) + (25*nc)/(9*sw2), &
               -4/sw2 - 11*nc/(9*sw2)]

        EWctWWAZ=EWctWWAZ+[3/(cw*sw) + (25*nc)/(9.*cw*sw) - 11/(4.*cw*sw3) - &
                 (14*cw)/(3.*sw3) - (11*nc)/(4.*cw*sw3), &
                  -(1/(cw*sw)) - (11*nc)/(9.*cw*sw) + 5/(4.*cw*sw3) + (3*cw)/sw3 &
                  + (5*nc)/(4.*cw*sw3)]

        EWctWWZZ=EWctWWZZ+[3/cw2 + (25*nc)/(9*cw2) + 14/(3*sw4) + 11*(1+nc)/(4*cw2*sw4) &
                           - 14/(3*sw2) - 11*(1+nc)/(2*cw2*sw2), &
                           -1/cw2 - 11*nc/(9*cw2) - 3/(sw4) - 5*(1+nc)/(4*cw2*sw4) &
                           + 3/sw2 + 5*(1+nc)/(2*cw2*sw2)]

        EWctWWWW=EWctWWWW+[-17/(2.*sw4) - (5*nc)/(2.*sw4), &
                            19/(6.*sw4) + (3*nc)/(2.*sw4)]

  end if



! add fourpoint-Tadpoles (calculated in D=4)

  if (TP_is_on /= 0 .and. SwB /= 0) then

    EWctAA = EWctAA + [ ZERO, 8.*A0W, ZERO ]

    EWctZZ = EWctZZ + [ ZERO, 0.25/(sw2*cw2)*(A0Z+A0H) + &
                      ((sw2-cw2)**2/(2.*sw2*cw2)+6.*cw2/sw2)*A0W, ZERO]

    EWctAZ = EWctAZ + [ ZERO, ((sw2-cw2)/(sw*cw)-6.*cw/sw)*A0W, ZERO]

    EWctWW = EWctWW + [ ZERO, (0.25/sw2+3.*cw2/sw2)*A0Z + 3.5/sw2*A0W + &
                               0.25/sw2*A0H, ZERO]


    EWctHH = EWctHH+[ZERO, -1./(2.*sw2)*(4.*A0W+2./cw2*A0Z+ &
                    (3*MH2)/(4*MW2)*A0H+MH2/(4*MW2)*A0Z+MH2/(2*MW2)*A0W)]

    EWctXX = EWctXX+[ZERO, -1./(2.*sw2)*(4.*A0W+2./cw2*A0Z+ &
                    MH2/(4*MW2)*A0H+(3.*MH2)/(4*MW2)*A0Z+MH2/(2*MW2)*A0W)]

    EWctPP = EWctPP+[ZERO, -1./(2.*sw2)*(4.*A0W+2.*(sw2-cw2)**2/cw2*A0Z+ &
                    MH2/(4*MW2)*A0H+MH2/(4*MW2)*A0Z+MH2/MW2*A0W)]


  end if


        if (debug_ew_renorm .eq. 2) then
          print*, "DEBUG"
          print*, "DEBUG EWctHH=", EWctHH
          print*, "DEBUG EWctXX=", EWctXX
          print*, "DEBUG EWctPP=", EWctPP
          print*, "DEBUG EWctAA=", EWctAA
          print*, "DEBUG EWctAZ=", EWctAZ
          print*, "DEBUG EWctZZ=", EWctZZ
          print*, "DEBUG EWctWW=", EWctWW
          print*, "DEBUG EWctuu=", EWctuu
          print*, "DEBUG EWctdd=", EWctdd
          print*, "DEBUG EWcttt=", EWcttt
          print*, "DEBUG EWctbb=", EWctbb
          print*, "DEBUG EWctee=", EWctee
          print*, "DEBUG EWctnn=", EWctnn
          print*, "DEBUG EWctVbt=", EWctVbt
          print*, "DEBUG EWctVdu=", EWctVdu
        end if


! set back possible complex masses
  if ( cms_on == 0 ) then
    call masspowers(rME, wME, ME, ME2, rME2)
    call masspowers(rMM, wMM, MM, MM2, rMM2)
    call masspowers(rML, wML, ML, ML2, rML2)
    call masspowers(rMU, wMU, MU, MU2, rMU2)
    call masspowers(rMD, wMD, MD, MD2, rMD2)
    call masspowers(rMS, wMS, MS, MS2, rMS2)
    call masspowers(rMC, wMC, MC, MC2, rMC2)
    call masspowers(rMB, wMB, MB, MB2, rMB2)
    call masspowers(rMT, wMT, MT, MT2, rMT2)
    call masspowers(rMW, wMW, MW, MW2, rMW2)
    call masspowers(rMZ, wMZ, MZ, MZ2, rMZ2)
    call masspowers(rMX, wMX, MX, MX2, rMX2)
    call masspowers(rMY, wMY, MY, MY2, rMY2)
    call masspowers(rMH, wMH, MH, MH2, rMH2)
  end if


   contains

  subroutine masspowers(rM, Ga, M, M2, rM2)
    use KIND_TYPES, only: REALKIND
    use ol_parameters_decl_/**/REALKIND, only: CI
    implicit none
    real(REALKIND),    intent(in)  :: rM, Ga
    complex(REALKIND), intent(out) :: M,  M2
    real(REALKIND),    intent(out) :: rM2
    M2  = rM*rM - CI*rM*Ga
    M   = sqrt(M2)
    rM2 = real(M2)
  end subroutine masspowers

end subroutine ew_renormalisation


subroutine photon_factors(photonid, ew_renorm, bornfactor, loopfactor)
! ****************************************************************************************
! EW photon factors.
! bornfactor: alpha(0)/alpha for each on-shell photo
! loopfactor: dZe-dZe(0) for each on-shell photon
!               + dZAA(epsilon)-dZAA(DeltaAlpha) for each initial-state off-shell photon
! ****************************************************************************************
  use ol_parameters_decl_/**/DREALKIND, only: ew_scheme, ew_renorm_scheme, &
                     &  onshell_photons_lsz, offshell_photons_lsz, delta_alphamz_dimreg
  use ol_parameters_decl_/**/REALKIND, only: alpha_QED, alpha_QED_0, alpha_QED_Gmu, pi
  use ol_loop_parameters_decl_/**/REALKIND, only: countertermnorm, &
                                              & dZeQEDEW, dZe0QEDEWnreg, dZe0QEDEWdreg, dZeGmuQEDEW, &
                                              & dZAAEW, dZAAEWdreg
  implicit none
  integer, intent(in) :: photonid(:)
  integer, intent(in) :: ew_renorm
  real(REALKIND),  intent(out) :: bornfactor
  real(REALKIND),  intent(out), optional :: loopfactor
  integer :: n_gamma=0, n_onshell_gamma=0, n_offshell_gamma=0

  !n_gamma=sum(photonid/photonid,photonid.ne.0)
  n_onshell_gamma=sum(photonid/photonid,photonid.gt.0)
  n_offshell_gamma=sum(photonid/photonid,photonid.lt.0)

  bornfactor=1.
  if (onshell_photons_lsz .and. (ew_scheme > 0)) then
    bornfactor=(alpha_QED_0/alpha_QED)**n_onshell_gamma
  end if

  if (ew_scheme == 0 .and. offshell_photons_lsz) then
    bornfactor=(alpha_QED_Gmu/alpha_QED_0)**n_offshell_gamma
  end if


  if (present(loopfactor) .and. ew_renorm > 0) then
    loopfactor=0.

    ! on-shell photons: dZe(Gmu/MZ) -> dZe(0)
    if (onshell_photons_lsz .and. (ew_renorm_scheme > 0)) then
      loopfactor=(dZe0QEDEWnreg-dZeQEDEW)*alpha_QED*countertermnorm*4.*pi*n_onshell_gamma
    end if

    ! off-photons: dZe(0,n-reg) -> dZe(Gmu)
    if (offshell_photons_lsz .and. ew_renorm_scheme == 0) then
        loopfactor=loopfactor+(dZeGmuQEDEW-dZeQEDEW)*alpha_QED*countertermnorm*4.*pi*n_offshell_gamma
    end if

    ! #off-shell photons: dZAA(n-reg) -> dZAA(dimreg)
    if (.not. delta_alphamz_dimreg .and. offshell_photons_lsz) then
      loopfactor=loopfactor+(dZAAEWdreg-dZAAEW)*alpha_QED*countertermnorm*2*pi*n_offshell_gamma
    end if

  end if
end subroutine photon_factors


end module ol_ew_renormalisation_/**/REALKIND
