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


module ol_qcd_renormalisation_/**/REALKIND
  implicit none
  contains

! **********************************************************************
subroutine qcd_renormalisation
! **********************************************************************
! QCD renormalisation constants (in physical GeV units);
! conventions for UV/IR div. see subroutine loop_parameters_init.
! This subroutine is automatically called by qcd_parameters_init.
! **********************************************************************
! dZMC  = charm-quark mass RC        : MC_bare = MC*(1+dZMC)
! dZMB  = bottom-quark mass RC       : MB_bare = MB*(1+dZMB)
! dZMT  = top-quark mass RC          : MT_bare = MT*(1+dZMT)
! dZYC  = charm-quark yukawa RC      : YC_bare = YC*(1+dZYC)
! dZYB  = bottom-quark yukawa RC     : YB_bare = YB*(1+dZYB)
! dZYT  = top-quark yukawa RC        : YT_bare = YT*(1+dZYT)
! dZg   = gluon-field RC             : A_bare  = (1+1/2*dZg)*A_ren
! dZq   = massless-quark field RC    : Q_bare  = (1+1/2*dZq)*Q_ren
! dZb   = bottom-quark field RC      : idem
! dZt   = top-quark field RC         : idem
! dgQCD = strong coupling RC         : g_bare  = (1+dgQCD)*g_ren
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_msg, ol_error, ol_fatal
  use ol_generic, only: to_string
  use ol_parameters_decl_/**/REALKIND
  use ol_loop_parameters_decl_/**/REALKIND
#ifndef PRECISION_dp
  use ol_parameters_decl_/**/DREALKIND, only: LeadingColour, model
  use ol_loop_parameters_decl_/**/DREALKIND, only: &
    & nc, nf, nf_up, nf_down, N_lf, nq_nondecoupl, CT_is_on, R2_is_on, SwF, SwB, &
    bubble_vertex
#endif
  implicit none

  real(REALKIND)    :: deC_UV, deC_IR, deB_UV, deB_IR, deT_UV, deT_IR
  complex(REALKIND) :: dummy_complex
#ifdef PRECISION_dp
  logical :: zeromasses(6) ! maximal possible number of quarks; only to determine N_lf

  zeromasses = [MU==0,MD==0,MS==0,MC==0,MB==0,MT==0]
  N_lf = count(zeromasses(:nf))

  if (N_lf /= 3 .and. N_lf /= 4 .and. N_lf /= 5) then
    call ol_error(2, 'in qcd_renormalisation:')
    call ol_msg( 'N_lf = ' // to_string(N_lf) // 'is not supported.')
    call ol_fatal()
  end if
#endif

  if (LeadingColour == 0) then
    cf = (nc*nc-1)/(2._/**/REALKIND*nc)
    ca = nc
    tf = 0.5_/**/REALKIND
  else
    cf = 0.5_/**/REALKIND*nc
    ca = nc
    tf = 0
  end if

  ! modified poles ((mu^2/MT^2)^eps)/eps
  if (MC /= 0) then
    deC_UV = de1_UV + real(log(mu2_UV/MC2))
    deC_IR = de1_IR + real(log(mu2_IR/MC2))
  end if
  if (MB /= 0) then
    deB_UV = de1_UV + real(log(mu2_UV/MB2))
    deB_IR = de1_IR + real(log(mu2_IR/MB2))
  end if
  deT_UV   = de1_UV + real(log(mu2_UV/MT2))
  deT_IR   = de1_IR + real(log(mu2_IR/MT2))

  dZq   = 0
  dZc   = 0
  dZb   = 0
  dZt   = 0
  dZMC  = 0
  dZMB  = 0
  dZMT  = 0
  dZYC  = 0
  dZYB  = 0
  dZYT  = 0
  dZg   = 0
  dgQCD = 0

  if (SwB /= 0) then
    ! non-fermionic
    ! On-shell renormalisation constant for gluon wave functions
    dZg   = (5*ca)/3 * (deT_UV - deT_IR)
    ! On-shell renormalisation constants for quark wave functions
    dZq   = -cf * (deT_UV - deT_IR)       ! massless quarks
    dZc   = dZq
    dZb   = dZq
    dZt   = dZq
    if (MC /= 0) then
      dZc  = -cf * (deC_UV + 2*deC_IR + 4) ! massive charm-quark
      if (LambdaMC2 == 0) then
        dZMC = -cf * (4 + 3*(de1_UV+log(mu2_UV/MC2))) ! on-shell
      else
        dZMC = -cf * 3*(de1_UV+log(mu2_UV/LambdaMC2)) ! MSbar
      end if
    end if
    if (YC /=0) then
      if (LambdaYC2 == 0) then
        dZYC = -cf * (4 + 3*(de1_UV+log(mu2_UV/YC2))) ! on-shell
      else
        dZYC = -cf * 3*(de1_UV+log(mu2_UV/LambdaYC2)) ! MSbar
      end if
    end if
    if (MB /= 0) then
      dZb  = -cf * (deB_UV + 2*deB_IR + 4) ! massive bottom-quark
      if (LambdaMB2 == 0) then
        dZMB = -cf * (4 + 3*(de1_UV+log(mu2_UV/MB2))) ! on-shell
      else
        dZMB = -cf * 3*(de1_UV+log(mu2_UV/LambdaMB2)) ! MSbar
      end if
    end if
    if (YB /=0) then
      if (LambdaYB2 == 0) then
        dZYB = -cf * (4 + 3*(de1_UV+log(mu2_UV/YB2))) ! on-shell
      else
        dZYB = -cf * 3*(de1_UV+log(mu2_UV/LambdaYB2)) ! MSbar
      end if
    end if
    ! top-mass renormalisation
    ! on-shell at complex pole p^2 = MT^2: dMT = MT * (-cf * (4 + 3*(de1_UV+log(mu2_UV/MT2))) + dZt)
    ! MSbar:                               dMT = MT * (-cf * (3*(de1_UV+log(mu2_UV/LambdaMT2))) + dZt)
    if (MT /= 0) then
      dZt   = -cf * (deT_UV + 2*deT_IR + 4) ! massive top-quark
      if (LambdaMT2 == 0) then
        dZMT = -cf * (4 + 3*(de1_UV+log(mu2_UV/MT2))) ! on-shell
      else
        dZMT = -cf * 3*(de1_UV+log(mu2_UV/LambdaMT2)) ! MSbar
      end if
    end if
    if (LambdaYT2 == 0) then
      dZYT = -cf * (4 + 3*(de1_UV+log(mu2_UV/YT2))) ! on-shell
    else
      dZYT = -cf * 3*(de1_UV+log(mu2_UV/LambdaYT2)) ! MSbar
    end if
    ! MS-bar renormalization constant for gQCD, YM-contribution
    dgQCD = -(11*ca)/6 * (de1_UV + log(mu2_UV/muren2))
  end if
  if (SwF /= 0) then
    ! fermionic
    dZg   = dZg - (4*tf)/3 * N_lf * (deT_UV - deT_IR)
    ! MS-bar renormalization constant for gQCD; contribution for nf quarks
    dgQCD = dgQCD + (2*tf*nf)/3 * (de1_UV + log(mu2_UV/muren2))
    if (nf > 3) then
      if (MC /= 0) then
        dZg = dZg - (4*tf)/3 * deC_UV
        if (nq_nondecoupl < 4 .or. muren <= rMC) then
          dgQCD = dgQCD + (2*tf)*real(log(muren2/MC2))/3
        end if
      end if
    end if
    if (nf > 4) then
      if (MB /= 0) then
        dZg = dZg - (4*tf)/3 * deB_UV
        if (nq_nondecoupl < 5 .or. muren <= rMB) then
          dgQCD = dgQCD + (2*tf)*real(log(muren2/MB2))/3
        end if
      end if
    end if
    if (nf > 5) then
      if (MT /= 0) then
        dZg = dZg - (4*tf)/3 * deT_UV
        ! top-quark decoupling term
        if (nq_nondecoupl < 6 .or. muren <= rMT) then
          dgQCD = dgQCD + (2*tf)*real(log(muren2/MT2))/3
        end if
      end if
    end if
  end if

!   ! MS-bar renormalisation constants
!   if (SwB /= 0) then
!     dZq  = -cf * de1_UV
!     dZb  = dZq
!     dZt  = dZq
!     dZg  = (5*ca)/3 * de1_UV
!     dZMB = -3*cf * de1_UV
!     dZMT = dZMB
!     dgQCD = -(11*ca)/6 * de1_UV
!   end if
!   if (SwF /= 0) then
!     ! only fermionic contribution
!     dgQCD = dgQCD + nf * (de1_UV/3)
!     dZg   = dZg - 2*nf * (de1_UV/3)
!   end if

  ! Sum of squared quark masses
  MQ2sum = MU2 + MD2 + MS2
  MQ2sumpairs = YU2 + YD2
  YQD2sum = YD2 + YS2
  YQU2sum = YU2
  YQD2sumpairs = YD2
  YQU2sumpairs = YU2
  nf_up = 1
  nf_down = 2
  YC2pair = 0
  YB2pair = 0
  YT2pair = 0
  if (nf > 3) then
    MQ2sum = MQ2sum + MC2
    MQ2sumpairs = MQ2sumpairs + YS2 + YC2
    YQU2sum = YQU2sum + YC2
    YQD2sumpairs = YQD2sumpairs + YS2
    YQU2sumpairs = YQU2sumpairs + YC2
    nf_up = nf_up + 1
    YC2pair = YC2
  end if
  if (nf > 4) then
    MQ2sum = MQ2sum + MB2
    YQD2sum = YQD2sum + YB2
    nf_down = nf_down + 1
  end if
  if (nf > 5) then
    MQ2sum = MQ2sum + MT2
    MQ2sumpairs = MQ2sumpairs + YB2 + YT2
    YQU2sum = YQU2sum + YT2
    YQD2sumpairs = YQD2sumpairs + YB2
    YQU2sumpairs = YQU2sumpairs + YT2
    nf_up = nf_up + 1
    YB2pair = YB2
    YT2pair = YT2
  end if

  ! set all counterterms to zero
  ctqq   = 0
  ctcc   = 0
  ctbb   = 0
  cttt   = 0
  ctGG   = 0
  ctGqq  = 0
  ctGcc  = 0
  ctGbb  = 0
  ctGtt  = 0
  ctVVV  = 0
  ctVVVV = 0
  ctVdt  = 0
  ctVsc  = 0
  ctVst  = 0
  ctVbu  = 0
  ctVbc  = 0
  ctVbt  = 0
  ctVtt  = 0
  ctVcc  = 0
  ctVbb  = 0
  ctVqq  = 0
  ctSud  = 0
  ctSus  = 0
  ctSub  = 0
  ctScd  = 0
  ctScs  = 0
  ctScb  = 0
  ctStd  = 0
  ctSts  = 0
  ctStb  = 0
  ctSdu  = 0
  ctSdc  = 0
  ctSdt  = 0
  ctSsu  = 0
  ctSsc  = 0
  ctSst  = 0
  ctSbu  = 0
  ctSbc  = 0
  ctSbt  = 0
  ctSqq  = 0
  ctScc  = 0
  ctSbb  = 0
  ctStt  = 0
  ! set pure R2 terms to zero
  ctZGG  = 0
  ctHGG  = 0
  ctAAGG = 0
  ctAZGG = 0
  ctZZGG = 0
  ctWWGG = 0
  ctHHGG = 0
  ctHXGG = 0
  ctXXGG = 0
  ctPPGG = 0
  ctAGGG = [ 0, 0 ]
  ctZGGG = [ 0, 0 ]
  R2GGGG = 0
  ctHEFTggh   = [0,0,0,0,0]
  ctHEFTgggh  = 0
  ctHEFTggggh = 0
  R2HEFTggggh = 0
  R2HEFThqq   = 0
  R2HEFTghqq  = 0


  if (CT_is_on /= 0) then
    ! Only UV counterterms

    if (bubble_vertex .eq. 1) then
      ctqq   = [ ZERO, ZERO ]
      ctcc   = [ ZERO, ZERO ]
      ctbb   = [ ZERO, ZERO ]
      cttt   = [ ZERO, ZERO ]
      ctGG   = [ ZERO, ZERO , ZERO ]
    else
      dummy_complex = dZq
      ctqq   = [ dummy_complex, ZERO ]
      dummy_complex = dZc
      ctcc   = [ dummy_complex, MC * (dZMC + dZc) ]
      dummy_complex = dZb
      ctbb   = [ dummy_complex, MB * (dZMB + dZb) ]
      dummy_complex = dZt
      cttt   = [ dummy_complex, MT * (dZMT + dZt) ]
      dummy_complex = dZg
      ctGG   = [ dummy_complex, ZERO , ZERO ]
    end if
    ctGqq = (dgQCD + dZq + dZg/2)
    ctGcc = (dgQCD + dZc + dZg/2)
    ctGbb = (dgQCD + dZb + dZg/2)
    ctGtt = (dgQCD + dZt + dZg/2)
    ctVVV  = (dgQCD + 1.5_/**/REALKIND * dZg)
    ctVVVV = (2 * (dgQCD + dZg))
    ctVdt  = (dZt/2 + dZq/2)
    ctVsc  = (dZc/2 + dZq/2)
    ctVst  = (dZt/2 + dZq/2)
    ctVbu  = (dZq/2 + dZb/2)
    ctVbc  = (dZc/2 + dZb/2)
    ctVbt  = (dZt/2 + dZb/2)
    ctVtt  = (dZt)
    ctVcc  = (dZc)
    ctVbb  = (dZb)
    ctVqq  = (dZq)
    ctSud  = [ -YD * (dZq/2 + dZq/2       ), YU * (dZq/2 + dZq/2) ]
    ctSus  = [ -YS * (dZq/2 + dZq/2       ), YU * (dZq/2 + dZq/2) ]
    ctSub  = [ -YB * (dZq/2 + dZb/2 + dZYB), YU * (dZq/2 + dZb/2) ]
    ctScd  = [ -YD * (dZc/2 + dZq/2       ), YC * (dZc/2 + dZq/2 + dZYC) ]
    ctScs  = [ -YS * (dZc/2 + dZq/2       ), YC * (dZc/2 + dZq/2 + dZYC) ]
    ctScb  = [ -YB * (dZc/2 + dZb/2 + dZYB), YC * (dZc/2 + dZb/2 + dZYC) ]
    ctStd  = [ -YD * (dZt/2 + dZq/2       ), YT * (dZt/2 + dZq/2 + dZYT) ]
    ctSts  = [ -YS * (dZt/2 + dZq/2       ), YT * (dZt/2 + dZq/2 + dZYT) ]
    ctStb  = [ -YB * (dZt/2 + dZb/2 + dZYB), YT * (dZt/2 + dZb/2 + dZYT) ]
    ctSdu  = [ -YU * (dZq/2 + dZq/2       ), YD * (dZq/2 + dZq/2) ]
    ctSdc  = [ -YC * (dZq/2 + dZc/2 + dZYC), YD * (dZq/2 + dZc/2) ]
    ctSdt  = [ -YT * (dZq/2 + dZt/2 + dZYT), YD * (dZq/2 + dZt/2) ]
    ctSsu  = [ -YU * (dZq/2 + dZq/2       ), YS * (dZq/2 + dZq/2) ]
    ctSsc  = [ -YC * (dZq/2 + dZc/2 + dZYC), YS * (dZq/2 + dZc/2) ]
    ctSst  = [ -YT * (dZq/2 + dZt/2 + dZYT), YS * (dZq/2 + dZt/2) ]
    ctSbu  = [ -YU * (dZb/2 + dZq/2       ), YB * (dZb/2 + dZq/2 + dZYB) ]
    ctSbc  = [ -YC * (dZb/2 + dZq/2 + dZYB), YB * (dZb/2 + dZq/2 + dZYB) ]
    ctSbt  = [ -YT * (dZb/2 + dZt/2 + dZYT), YB * (dZb/2 + dZt/2 + dZYB) ]
    ctSqq  = (dZq)
    ctScc  = (dZc + dZYC)
    ctSbb  = (dZb + dZYB)
    ctStt  = (dZt + dZYT)

    ! HEFT
    ! finite two-loop contribution: gGGH*11/4/pi*as -> 11
    ctHEFTggh   = (11+2*dgQCD + dZg) * [ rZERO, -rONE, rZERO, rZERO, rONE]
    ctHEFTgggh  = (11+3*dgQCD + 1.5_/**/REALKIND*dZg)
    ctHEFTggggh = (11+4*dgQCD + 2*dZg)

    ! 2HDM
    if (trim(model) == "2hdm") then
      thdmctHpcs = [ -YS*(-thdmYuk3) * (dZc/2 + dZq/2       ), YC/thdmTB * (dZc/2 + dZq/2 + dZYC) ]
      thdmctHpsc = [ -YC/thdmTB * (dZc/2 + dZq/2 + dZYC), YS*(-thdmYuk3) * (dZc/2 + dZq/2       ) ]
      thdmctHptb = [ -YB*(-thdmYuk3) * (dZt/2 + dZb/2 + dZYB), YT/thdmTB * (dZt/2 + dZb/2 + dZYT) ]
      thdmctHpbt = [ -YT/thdmTB * (dZt/2 + dZb/2 + dZYT), YB*(-thdmYuk3) * (dZt/2 + dZb/2 + dZYB) ]
    end if
  end if

  if (R2_is_on /= 0) then
    ! Add R2 contribution
    ! ctff = [ dZf - cf , mf * (dZmf + dZf - 2*cf) ]
    ! ctVV = [ dZV - (ca/3*(1/2+LHV) + 2*tf*nf/3) , - m2*(dZm2+dZV) + 4*tf*MQ2sum , LHV ]
    if (SwB /= 0) then
      ! non-fermionic
      dummy_complex = -cf
      if (bubble_vertex .eq. 0) then
        ctqq   = ctqq + [ dummy_complex, ZERO ]
        ctcc   = ctcc + [ dummy_complex, -2*MC*cf ]
        ctbb   = ctbb + [ dummy_complex, -2*MB*cf ]
        cttt   = cttt + [ dummy_complex, -2*MT*cf ]
      end if
      dummy_complex = -0.5_/**/REALKIND*ca
      if (bubble_vertex .eq. 0) then
        ctGG   = ctGG + [ dummy_complex, ZERO , cONE ]
      end if
      ctGqq = ctGqq - 2*cf
      ctGcc = ctGcc - 2*cf
      ctGbb = ctGbb - 2*cf
      ctGtt = ctGtt - 2*cf
      ctVVV  = ctVVV - (11*ca)/12
      ! ctVVVV = ctVVVV; handled in the Feynman rule
      ctVdt  = ctVdt - 2*cf
      ctVsc  = ctVsc - 2*cf
      ctVst  = ctVst - 2*cf
      ctVbu  = ctVbu - 2*cf
      ctVbc  = ctVbc - 2*cf
      ctVbt  = ctVbt - 2*cf
      ctVtt  = ctVtt - 2*cf
      ctVcc  = ctVcc - 2*cf
      ctVbb  = ctVbb - 2*cf
      ctVqq  = ctVqq - 2*cf
      ctSud  = ctSud - 4*cf*gPud
      ctSus  = ctSus - 4*cf*gPus
      ctSub  = ctSub - 4*cf*gPub
      ctScd  = ctScd - 4*cf*gPcd
      ctScs  = ctScs - 4*cf*gPcs
      ctScb  = ctScb - 4*cf*gPcb
      ctStd  = ctStd - 4*cf*gPtd
      ctSts  = ctSts - 4*cf*gPts
      ctStb  = ctStb - 4*cf*gPtb
      ctSdu  = ctSdu - 4*cf*gPdu
      ctSdc  = ctSdc - 4*cf*gPdc
      ctSdt  = ctSdt - 4*cf*gPdt
      ctSsu  = ctSsu - 4*cf*gPsu
      ctSsc  = ctSsc - 4*cf*gPsc
      ctSst  = ctSst - 4*cf*gPst
      ctSbu  = ctSbu - 4*cf*gPbu
      ctSbc  = ctSbc - 4*cf*gPbc
      ctSbt  = ctSbt - 4*cf*gPbt
      ctSqq  = ctSqq - 4*cf
      ctScc  = ctScc - 4*cf
      ctSbb  = ctSbb - 4*cf
      ctStt  = ctStt - 4*cf

      ! HEFT
      ctHEFTggh   = ctHEFTggh + [ 1, 89, 14, -17, -93 ] * (nc/24._/**/REALKIND)
      ctHEFTgggh  = ctHEFTgggh - 15*nc/8._/**/REALKIND
      R2HEFTggggh = 0.125_/**/REALKIND
      R2HEFThqq   = 0.25_/**/REALKIND*(nc-1._/**/REALKIND/nc)
      R2HEFTghqq  = 0.25_/**/REALKIND*(3._/**/REALKIND/nc-5*nc)

      ! 2HDM
     if (trim(model) == "2hdm") then
       thdmctHpcs = thdmctHpcs - 4*cf*thdmHpcs
       thdmctHpsc = thdmctHpsc - 4*cf*thdmHpsc
       thdmctHptb = thdmctHptb - 4*cf*thdmHptb
       thdmctHpbt = thdmctHpbt - 4*cf*thdmHpbt
     end if
    end if
    if (SwF /= 0) then
      ! fermionic
      dummy_complex = -(2*tf*nf)/3
      if (bubble_vertex .eq. 0) then
        ctGG   = ctGG + [ dummy_complex, 4*tf*MQ2sum , ZERO ]
      end if
      ctVVV  = ctVVV - (4*tf*nf)/3
      ! pure R2 terms
      ! ZGG R2 coupling: 4/3*sum_q(a_q)
      ctZGG = (nf_down-nf_up)*(-2*tf)/(3*cw*sw)
      ! HGG R2 coupling: 2*sum_q(m_q*v_q)
      ctHGG = -2*tf*MQ2sum/(sw*MW) ! *MQ*YQ/MQ2sum in Feynman rule
      ! VVGG R2 coupling: 2/3*sum(v1*v2+a1*a2)
      ctAAGG = (nf_up*4+nf_down)*(tf*4)/27
      ctAZGG = (nf_up*(-6+16*sw2)+nf_down*(-3+4*sw2)) * tf/(27*cw*sw)
      ctZZGG = (nf_up*(9-24*sw2+32*sw2**2)+nf_down*(9-12*sw2+ 8*sw2**2)) * tf/(54*cw2*sw2)
      ctWWGG = int(nf/2)*tf/(3*sw2) ! tf/(3*sw2) per SU(2) doublet
      ! SSGG R2 coupling: 2*sum(v1*v2-a1*a2)
      ctHHGG = tf*MQ2sum/(sw2*MW2) ! *YQ2/MQ2sum in Feynman rule
      ctHXGG = 0
      ctXXGG = tf*MQ2sum/(sw2*MW2) ! *YQ2/MQ2sum in Feynman rule
      ctPPGG = tf*MQ2sumpairs/(sw2*MW2)
      ! VGGG R2 coupling: [ 4/3*tf*sum_q(v_q), -12*tf*CI*sum_q(a_q) ]
      dummy_complex = -(4*tf)/9*(2*nf_up-nf_down)
      ctAGGG = [ dummy_complex, ZERO ]
      ctZGGG = [ (nf_up*(3-8*sw2)+nf_down*(-3+4*sw2)) * tf/(9*cw*sw) , (nf_up-nf_down)*3*tf*CI/(sw*cw) ]
      R2GGGG = int(2*tf) ! switch on the 4-gluon R2
      ! 2HDM
      if (trim(model) == "2hdm") then
        ! SGG
        ! thdmctA0GG = 0
        thdmctHGG = -(thdmYuk1*YQD2sum + thdmCA/thdmSB*YQU2sum)/(MW*sw)
        thdmctHh0GG = -(thdmYuk2*YQD2sum + thdmSA/thdmSB*YQU2sum)/(MW*sw)
        ! SSGG
        ! thdmctHA0GG = 0
        ! thdmctXHhGG = 0
        ! thdmctHhA0GG = 0
        ! thdmctHXGG = 0 (effectively unchanged wrt. SM)
        ! XXGG unchanged wrt. SM
        ! PPGG unchanged wrt. SM
        thdmctHHGG = (thdmYuk1**2*YQD2sum + thdmCA**2/thdmSB**2*YQU2sum)/(2*MW2*sw2)
        thdmctHHhGG = (thdmYuk1*thdmYuk2*YQD2sum + thdmCA*thdmSA/thdmSB**2*YQU2sum)/(2*MW2*sw2)
        thdmctXA0GG = (-thdmYuk3*YQD2sum + YQU2sum/thdmTB)/(2*MW2*sw2)
        thdmctHhHhGG = (thdmYuk2**2*YQD2sum + thdmSA**2/thdmSB**2*YQU2sum)/(2*MW2*sw2)
        thdmctA0A0GG = (thdmYuk3**2*YQD2sum + YQU2sum/thdmTB**2)/(2*MW2*sw2)
        thdmctPHpGG = (-thdmYuk3*YQD2sumpairs + YQU2sumpairs/thdmTB)/(2*MW2*sw2)
        thdmctHpHpGG = (thdmYuk3**2*YQD2sumpairs + YQU2sumpairs/thdmTB**2)/(2*MW2*sw2)
      end if
    end if

  end if
end subroutine qcd_renormalisation

end module ol_qcd_renormalisation_/**/REALKIND


module ol_qcd_offshell_selfenergies/**/REALKIND
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_loop_parameters_decl_/**/REALKIND, only: de1_IR, de1_UV, de2_i_IR
  use ol_debug, only: ol_msg, ol_error, ol_fatal
  use ol_generic, only: to_string
  use ol_self_energy_integrals_/**/REALKIND
  implicit none

  contains

    ! Gluon self-energy
function gluon_ofsse(p2,pid)
  use KIND_TYPES, only: DREALKIND, REALKIND
  use ol_parameters_decl_/**/REALKIND
  use ol_loop_parameters_decl_/**/REALKIND
#ifndef PRECISION_dp
  use ol_loop_parameters_decl_/**/DREALKIND, only: &
    & nc, nf, N_lf, bubble_vertex, CT_is_on, R2_is_on
#endif
  complex(REALKIND), intent(in) :: p2
  integer,           intent(in) :: pid
  complex(REALKIND) :: cc1,gluon_ofsse(3),B01,B02,B03,B04
  complex(REALKIND) :: cc1R1,dZ

  !complex(REALKIND) :: B0_1,B0_2,B1_1,B1_2,B11_1,B11_2,B00_1,B00_2
  real(REALKIND) :: ncc,nlf
  integer :: i

  if (pid .ne. 21) then
    call ol_fatal('Cannot use gluon_ofsse for pid other than 21. pid=' //  &
                  trim(to_string(pid)) // '.')
  end if

  ncc = real(nc,kind=REALKIND)
  nlf = real(N_lf,kind=REALKIND)
  ! gs**2/(16*pi**2) stripped, assert(Nf == 6)

  call init_ol_self_energy_integrals(.true.)
  B01 = calcB0(p2,ZERO,ZERO)
  cc1 = +(5*ncc-2*nlf)*(B01)/real(3,kind=REALKIND)

  do i = N_lf, nf-1
    select case (i)
    case (5)
      B02 = calcB0(ZERO,MT2,MT2)
      B03 = calcB0(p2,MT2,MT2)
      cc1 = cc1 + 4*MT2*(B02-B03)/(3*p2) &
                - 2*B03/3
    case (4)
      B02 = calcB0(ZERO,MB2,MB2)
      B03 = calcB0(p2,MB2,MB2)
      cc1 = cc1 + 4*MB2*(B02-B03)/(3*p2) &
                - 2*B03/3
    case (3)
      B02 = calcB0(ZERO,MC2,MC2)
      B03 = calcB0(p2,MC2,MC2)
      cc1 = cc1 + 4*MC2*(B02-B03)/(3*p2) &
                - 2*B03/3
    case default
      call ol_fatal('Flavor scheme N_lf=' // trim(to_string(N_lf)) // ' not implemented for bubble_vertex.')
  end select
  end do
  call init_ol_self_energy_integrals(.false.)


  if (CT_is_on .eq. 0) then
    dZ = 0
  else
    dZ = dZg
  end if

  ! only R1
  cc1R1 = (nf*2 + ncc)/real(9,kind=REALKIND)
  gluon_ofsse(1) =  (dZ - cc1 - cc1R1)  ! w^\mu_out = p^2 w^\mu_in
  gluon_ofsse(2) =  0                   ! w^\mu_out = w^\mu_in
  gluon_ofsse(3) = -(-cc1 + cc1R1)      ! w^\mu_out = (w_in.p) p^\mu

  if (R2_is_on .eq. 0) then
    gluon_ofsse(1) = gluon_ofsse(1) + ca/2 + (2*tf*nf)/3
    gluon_ofsse(2) = gluon_ofsse(2) - 4*tf*MQ2sum
    gluon_ofsse(3) = gluon_ofsse(3) - cONE
  end if

  end function

    ! Quark self-energy
function quark_ofsse(p2,pid)
  use KIND_TYPES, only: DREALKIND, REALKIND
  use ol_parameters_decl_/**/REALKIND
  use ol_loop_parameters_decl_/**/REALKIND
  use ol_self_energy_integrals_/**/REALKIND
#ifndef PRECISION_dp
  use ol_loop_parameters_decl_/**/DREALKIND, only: &
    & nc, N_lf, bubble_vertex, CT_is_on, R2_is_on
#endif
  complex(REALKIND), intent(in) :: p2
  integer,           intent(in) :: pid
  complex(REALKIND) :: cc1,cc2,cc3,quark_ofsse(2),B01,B02,B03
  complex(REALKIND) :: cc1R1,cc2R1,cc3R1,fac
  complex(REALKIND) :: cc1R2,cc2R2,cc3R2

  complex(REALKIND) :: B0_1,A0_1,dZ,dM,M,M2


  select case (pid)
  case (1)
    dZ = dZq
    dM = 0
    M = MD
    M2 = MD2
  case (2)
    dZ = dZq
    dM = 0
    M = MU
    M2 = MU2
  case (3)
    dZ = dZq
    dM = 0
    M = MS
    M2 = MS2
  case (4)
    dZ = dZc
    dM = dZMC
    M = MC
    M2 = MC2
  case (5)
    dZ = dZb
    dM = dZMB
    M = MB
    M2 = MB2
  case (6)
    dZ = dZt
    dM = dZMT
    M = MT
    M2 = MT2
  case default
    call ol_fatal('Cannot use quark_ofsse for pidse other than 1,2,3,4,5,6. ' // &
                  'pid=' // trim(to_string(pid)) // '.')
  end select
  if (CT_is_on .eq. 0) then
    dZ = 0
    dM = 0
  end if

  call init_ol_self_energy_integrals(.true.)
  B0_1 = calcB0(p2,ZERO,M2)
  A0_1 = calcA0(M2)
  call init_ol_self_energy_integrals(.false.)

  fac = rONE*cf
  quark_ofsse(1) = dZ + fac*((M2+p2) * B0_1 - A0_1 -p2)/p2  ! w^i_out = pslash^{ij} w^j_in
  quark_ofsse(2) = M*(dM + dZ) + fac*(4*M*B0_1-2*M)         ! w^i_out = w^i_in

  if (R2_is_on .eq. 0) then
    quark_ofsse(1) = quark_ofsse(1) + cf
    quark_ofsse(2) = quark_ofsse(2) + 2*M*cf
  end if

  end function

end module ol_qcd_offshell_selfenergies/**/REALKIND
