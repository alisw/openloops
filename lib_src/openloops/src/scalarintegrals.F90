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



module ol_self_energy_integrals_/**/REALKIND
! !interfaces for scalar one-point & two-point functions
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND, only: ZERO
  use ol_parameters_decl_/**/DREALKIND, only: cms_on
  use ol_loop_parameters_decl_/**/REALKIND, only: de1_UV, de1_IR
  use ol_loop_parameters_decl_/**/DREALKIND, only: se_integral_switch, coli_cache_use, &
    & a_switch
#ifdef USE_ONELOOP
    use avh_olo_/**/REALKIND
#endif
  implicit none
#ifndef PRECISION_dp
  integer, parameter, private :: dp = selected_real_kind(15)
  real(REALKIND), private :: eps=1.e-33
#else
  real(REALKIND), private :: eps=1.e-17
#endif

  contains

  subroutine init_ol_self_energy_integrals(init)
#ifdef USE_COLLIER
  use collier, only: setmode_cll
  use cache, only: SwitchOffCacheSystem_cll, SwitchOnCacheSystem_cll
#endif
  implicit none
  logical :: init
  if(init) then
    if ((se_integral_switch == 1 .or. se_integral_switch == 7) &
      .and. coli_cache_use == 1) then
        call SwitchOffCacheSystem_cll
    end if
    if (se_integral_switch == 1 .and. a_switch /= 1) then
      call setmode_cll(1)
    end if
    if (se_integral_switch == 7 .and. a_switch /= 7) then
      call setmode_cll(2)
    end if
  else
    if (se_integral_switch == 1 .and. a_switch == 1) then
      call setmode_cll(1)
    end if
    if (se_integral_switch == 7 .and. a_switch /= 1) then
      call setmode_cll(2)
    end if
    if ((se_integral_switch == 1 .or. se_integral_switch == 7) &
      .and. coli_cache_use == 1) then
        call SwitchOnCacheSystem_cll
    end if
  end if
  end subroutine


  function calcA0(m12_in)
#ifdef USE_COLLIER
    use collier_coefs, only: A0_cll
#endif
#ifdef USE_ONELOOP
    use avh_olo_/**/REALKIND
#endif
    complex(REALKIND) calcA0
    complex(REALKIND), intent(in) :: m12_in
    complex(DREALKIND) m12
    complex(DREALKIND) A0_coli
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt(0:2)
#endif

    calcA0 = 0

    m12 = m12_in

#if defined(USE_COLLIER)
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      call A0_cll(A0_coli,m12)
      calcA0 = A0_coli
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      call olo_a0(rslt,m12_in)
      calcA0 = rslt(0) + rslt(1)*de1_UV
    end if
#endif

      return
  end function calcA0

  function calcB0(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: B0_cll
#endif
#ifdef USE_ONELOOP
    use avh_olo_/**/REALKIND
#endif
    complex(REALKIND) calcB0
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    complex(REALKIND) p2q
    complex(DREALKIND) p2
    complex(DREALKIND) m12
    complex(DREALKIND) m22
    complex(DREALKIND) B0_coli
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt(0:2)
#endif

    calcB0 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#ifdef USE_COLLIER
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B0_cll(B0_coli,p2,m12,m22)
      calcB0 = B0_coli
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      call olo_b0(rslt,real(p2q),m12_in,m22_in)
      if (p2q .eq. 0 .and. m12_in .eq. 0 .and. m22_in .eq. 0) then
        calcB0 = de1_UV - de1_IR
      else
        calcB0 = rslt(0) + rslt(1)*de1_UV
      end if
    end if
#endif

      return
  end function calcB0


  function calcdB0(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: DB0_cll
#endif
    implicit none
    complex(REALKIND) calcdB0
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    real(REALKIND) :: p2q
    complex(DREALKIND) :: p2
    complex(DREALKIND) :: m12
    complex(DREALKIND) :: m22
    complex(DREALKIND) DB0_coli
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt(0:2)
#endif

    calcdB0 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#ifdef USE_COLLIER
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call DB0_cll(DB0_coli,p2,m12,m22)
      calcdB0 = DB0_coli
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      if (p2q == 0. .and. m12 == 0. .and. m22 == 0. ) then
        calcdB0 = 0.
      else
        call olo_db0(rslt,p2q,m12_in,m22_in)
        calcdB0 = rslt(0) + rslt(1)*de1_IR
      end if
    end if
#endif

    return
  end function calcdB0


  function calcB1(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: B_cll
#endif
    complex(REALKIND) calcB1
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    complex(DREALKIND) :: p2
    complex(DREALKIND) :: m12
    complex(DREALKIND) :: m22
    complex(DREALKIND) B1_coli
#ifdef USE_COLLIER
    complex(DREALKIND) B(0:1,0:1), Buv(0:1,0:1)
#endif
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)
#endif

    calcB1 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#ifdef USE_COLLIER
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,1)
      calcB1 = B(1,0)
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      if (p2 .eq. 0 .and. m12_in .eq. 0 .and. m22_in .eq. 0) then
        calcB1 = rslt_b1(0) + (de1_IR - de1_UV)/2
      else
        calcB1 = rslt_b1(0) + rslt_b1(1)*de1_UV
      end if
    end if
#endif

    return
  end function calcB1

  function calcB00(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: B_cll
#endif
    complex(REALKIND) calcB00
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    complex(DREALKIND) :: p2
    complex(DREALKIND) :: m12
    complex(DREALKIND) :: m22
    complex(DREALKIND) B1_coli
#ifdef USE_COLLIER
    complex(DREALKIND) B(0:2,0:1), Buv(0:2,0:1)
#endif
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)
#endif

    calcB00 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#if defined(USE_COLLIER)
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,2)
      calcB00 = B(1,0)
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      calcB00 = rslt_b00(0) + rslt_b00(1)*de1_UV
    end if
#endif

    return
  end function calcB00

  function calcB11(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: B_cll
#endif
    complex(REALKIND) calcB11
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    complex(DREALKIND) :: p2
    complex(DREALKIND) :: m12
    complex(DREALKIND) :: m22
    complex(DREALKIND) B1_coli
#ifdef USE_COLLIER
    complex(DREALKIND) B(0:1,0:2), Buv(0:1,0:2)
#endif
#ifdef USE_ONELOOP
    complex(REALKIND) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)
#endif

    calcB11 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#if defined(USE_COLLIER)
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,2)
      calcB11 = B(0,2)
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      calcB11 = rslt_b11(0) + rslt_b11(1)*de1_UV
    end if
#endif

    return
  end function calcB11

  function calcdB1(p2_in,m12_in,m22_in)
#ifdef USE_COLLIER
    use collier_coefs, only: DB1_cll
#endif
    implicit none
    complex(REALKIND) calcdB1
    complex(REALKIND), intent(in) :: p2_in
    complex(REALKIND), intent(in) :: m12_in
    complex(REALKIND), intent(in) :: m22_in
    real(REALKIND) :: p2q
    complex(DREALKIND) :: p2, m12, m22
    complex(DREALKIND) DB1_coli

    calcdB1 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in

#ifdef USE_COLLIER
    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call DB1_cll(DB1_coli,p2,m12,m22)
      calcdB1 = DB1_coli
    end if
#endif

#ifdef USE_ONELOOP
    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      if (abs(p2q) > eps) then
        calcdB1 = - (m22_in-m12_in)/2/p2q**2*(calcB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
       &          + (m22_in-m12_in-p2q)/2/p2q*calcdB0(p2_in,m12_in,m22_in)
      else
        if (abs(m12_in) < eps .and. abs(m22_in) < eps) then
          calcdB1 = 0
        else if (abs(m12_in) < eps .and. abs(m22_in) > eps) then
          calcdB1 = -1/m22_in/6
        else if (abs(m12_in) > eps .and. abs(m22_in) < eps) then
          calcdB1 = -1/m12_in/6
        else if (m12_in == m22_in) then
          calcdB1 = -1/m12_in/12
        else
          calcdB1 = -(2.*m12_in**3+3.*m12_in**2*m22_in-6.*m12_in*m22_in**2+m22_in**3 &
       &              +6.*m12_in**2*m22_in*log(m22_in/m12_in))/6./(m12_in-m22_in)**4
        end if
      end if
    end if
#endif

    return
  end function calcdB1


  function calcRB0(p2_in,m12_in,m22_in)
    complex(REALKIND) calcRB0
    complex(REALKIND), intent(in) :: p2_in, m12_in, m22_in

    calcRB0 = calcB0(p2_in,m12_in,m22_in)
    if (imag(p2_in) == 0 .and. cms_on == 1) then
      if (real(p2_in) .gt. real(m12_in+m22_in)) then
        calcRB0 = real(calcRB0)
      end if
    end if
  end function calcRB0

  function calcRdB0(p2_in,m12_in,m22_in)
    complex(REALKIND) calcRdB0
    complex(REALKIND), intent(in) :: p2_in, m12_in, m22_in

    calcRdB0 = calcdB0(p2_in,m12_in,m22_in)
    if (imag(p2_in) == 0 .and. cms_on == 1) then
      if (real(p2_in) .gt. real(m12_in+m22_in)) then
        calcRdB0 = real(calcRdB0)
      end if
    end if
  end function calcRdB0

  function calcRB1(p2_in,m12_in,m22_in)
    complex(REALKIND) calcRB1
    complex(REALKIND), intent(in) :: p2_in, m12_in, m22_in

    if (imag(p2_in) == 0 .and. cms_on == 1 &
      .and. (real(p2_in) > real(m12_in+m22_in)) &
      .and. (p2_in /= 0)) then
!(m12-m02)/2/p2*(B0 - B0(0,m0,m1)) - B0/2
          calcRB1 = (m22_in-m12_in)/2./p2_in*(calcRB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
                  - 1./2.*calcRB0(p2_in,m12_in,m22_in)
    else
      calcRB1 = calcB1(p2_in,m12_in,m22_in)
    end if
  end function calcRB1

  function calcRdB1(p2_in,m12_in,m22_in)
    complex(REALKIND) calcRdB1
    complex(REALKIND), intent(in) :: p2_in, m12_in, m22_in

    if (imag(p2_in) == 0 .and. cms_on == 1 &
      .and. (real(p2_in) > real(m12_in+m22_in)) &
      .and. (p2_in /= 0)) then
!-(m1^2-m0^2)/2/Q4*(real(B0(Q2,m0,m1)) - B0(0,m0,m1)) + (m1^2-m0^2-Q2)/2/Q2*real(dB0(Q2,m0,m1))
          calcRdB1 = -(m22_in-m12_in)/2./p2_in**2*(calcRB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
                  + (m22_in-m12_in-p2_in)/2./p2_in*calcRdB0(p2_in,m12_in,m22_in)
    else
      calcRdB1 = calcdB1(p2_in,m12_in,m22_in)
    end if

  end function calcRdB1


end module ol_self_energy_integrals_/**/REALKIND

