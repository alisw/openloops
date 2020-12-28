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


module ol_kinematics_/**/REALKIND
  use KIND_TYPES
  implicit none

  interface internal_momenta
    module procedure internal_momenta_six, internal_momenta_std
  end interface

  interface get_mass
    module procedure get_mass_id,get_mass_idarr
  end interface

  interface get_mass2
    module procedure get_mass2_id,get_mass2_idarr
  end interface

  interface get_rmass2
    module procedure get_rmass2_id,get_rmass2_idarr
  end interface

  interface conv_mom_scatt2in
    module procedure conv_mom_scatt2in_mexpl, conv_mom_scatt2in_mids
  end interface

  interface init_kinematics
    module procedure init_kinematics_mexpl, init_kinematics_mids
  end interface


  real(REALKIND) :: collthres = 1.E-4
  real(REALKIND) :: softthres = 1.E-4

contains

! **********************************************************************
subroutine Std2LC_Rep(P,L)
! **********************************************************************
! Lorentz -> light-cone representation
! P(0:3)   = Lorentz momentum P^mu (contravariant)
! L(1:4)   = light-cone representation L^A (contravariant)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND, only: CI
  implicit none
  real(REALKIND),    intent(in)  :: P(0:3)
  complex(REALKIND), intent(out) :: L(1:4)
  L(1) =  P(0) - P(3)
  L(2) =  P(0) + P(3)
  L(3) = -P(1) - CI*P(2)
  L(4) = -P(1) + CI*P(2)
end subroutine Std2LC_Rep


! **********************************************************************
subroutine Std2LC_cmplx(P,L)
! **********************************************************************
! complex version needed for OPP Reduction
! Lorentz -> light-cone representation
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND, only: CI
  implicit none
  complex(REALKIND), intent(in)  :: P(0:3)
  complex(REALKIND), intent(out) :: L(1:4)
  L(1) =  P(0) -    P(3)
  L(2) =  P(0) +    P(3)
  L(3) = -P(1) - CI*P(2)
  L(4) = -P(1) + CI*P(2)
end subroutine Std2LC_cmplx


! **********************************************************************
function get_LC_5(mom) result(P)
! **********************************************************************
! P(1:4) = light-cone representation
! P(5)   = mass
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  implicit none
  integer, intent(in) :: mom
  complex(REALKIND) :: P(1:5)
  if (mom .gt. 0) then
    P(1:4) = L(1:4,mom)
    P(5) = L(5,mom) + L(6,mom)
  else
    P(1:4) = -L(1:4,-mom)
    P(5) = L(5,-mom) + L(6,-mom)
  end if

end function get_LC_5

#ifdef PRECISION_dp
function get_LC_5_qp(mom) result(P)
  use KIND_TYPES, only: QREALKIND
  use ol_momenta_decl_/**/QREALKIND, only: L_qp=>L
  use ol_external_decl_/**/REALKIND, only: init_qp
  implicit none
  integer, intent(in) :: mom
  complex(QREALKIND) :: P(1:5)
  if (.not. init_qp) call init_qp_kinematics
  if (mom .gt. 0) then
    P(1:4) = L_qp(1:4,mom)
    P(5) = L_qp(5,mom) + L_qp(6,mom)
  else
    P(1:4) = -L_qp(1:4,-mom)
    P(5) = L_qp(5,-mom) + L_qp(6,-mom)
  end if

end function get_LC_5_qp
#endif



! **********************************************************************
function get_LC_4(mom) result(P)
! **********************************************************************
! P(1:4) = light-cone representation
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  implicit none
  integer, intent(in) :: mom
  complex(REALKIND) :: P(1:4)
  if (mom .gt. 0) then
    P(1:4) = L(1:4,mom)
  else
    P(1:4) = -L(1:4,-mom)
  end if

end function get_LC_4

#ifdef PRECISION_dp
function get_LC_4_qp(mom) result(P)
  use KIND_TYPES, only: QREALKIND
  use ol_momenta_decl_/**/QREALKIND, only: L_qp=>L
  use ol_external_decl_/**/REALKIND, only: init_qp
  implicit none
  integer, intent(in) :: mom
  complex(QREALKIND) :: P(1:4)
  if (.not. init_qp) call init_qp_kinematics
  if (mom .gt. 0) then
    P(1:4) = L_qp(1:4,mom)
  else
    P(1:4) = -L_qp(1:4,-mom)
  end if

end function get_LC_4_qp
#endif


! **********************************************************************
function get_LC_mass2(mom) result(m2)
! **********************************************************************
! m2   = mass squared
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  implicit none
  integer, intent(in) :: mom
  complex(REALKIND) :: m2
  if (mom .gt. 0) then
    m2 = L(5,mom) + L(6,mom)
  else
    m2 = L(5,-mom) + L(6,-mom)
  end if
end function get_LC_mass2


! **********************************************************************
function get_mass_id(mid) result(m)
! **********************************************************************
! m   = complex mass
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  complex(REALKIND) :: m

  select case (mid)
  case (0)
    m = ZERO
  case (nME)
    m = ME
  case (nMM)
    m = MM
  case (nML)
    m = ML
  case (nMU)
    m = MU
  case (nMD)
    m = MD
  case (nMS)
    m = MS
  case (nMC)
    m = MC
  case (nMB)
    m = MB
  case (nMT)
    m = MT
  case (nMW)
    m = MW
  case (nMZ)
    m = MZ
  case (nMH)
    m = MH
  case (nMX)
    m = MX
  case (nMY)
    m = MY
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_mass_id


! **********************************************************************
function get_mass_idarr(mids) result(m)
! **********************************************************************
! m   = complex mass array
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  implicit none
  integer, dimension(:), intent(in) :: mids
  complex(REALKIND), dimension(size(mids)) :: m
  integer :: i

  do i = 1, size(mids)
    m(i) = get_mass_id(mids(i))
  end do
end function get_mass_idarr


! **********************************************************************
function get_mass2_id(mid) result(m2)
! **********************************************************************
! m2   = complex mass squared
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  complex(REALKIND) :: m2

  select case (mid)
  case (0)
    m2 = ZERO
  case (nME)
    m2 = ME2
  case (nMM)
    m2 = MM2
  case (nML)
    m2 = ML2
  case (nMU)
    m2 = MU2
  case (nMD)
    m2 = MD2
  case (nMS)
    m2 = MS2
  case (nMC)
    m2 = MC2
  case (nMB)
    m2 = MB2
  case (nMT)
    m2 = MT2
  case (nMW)
    m2 = MW2
  case (nMZ)
    m2 = MZ2
  case (nMH)
    m2 = MH2
  case (nMX)
    m2 = MX2
  case (nMY)
    m2 = MY2
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_mass2_id


! **********************************************************************
function get_rmass2_id(mid) result(m2)
! **********************************************************************
! m2   = real mass squared
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  real(REALKIND) :: m2

  select case (mid)
  case (0)
    m2 = ZERO
  case (nME)
    m2 = rME2
  case (nMM)
    m2 = rMM2
  case (nML)
    m2 = rML2
  case (nMU)
    m2 = rMU2
  case (nMD)
    m2 = rMD2
  case (nMS)
    m2 = rMS2
  case (nMC)
    m2 = rMC2
  case (nMB)
    m2 = rMB2
  case (nMT)
    m2 = rMT2
  case (nMW)
    m2 = rMW2
  case (nMZ)
    m2 = rMZ2
  case (nMH)
    m2 = rMH2
  case (nMX)
    m2 = rMX2
  case (nMY)
    m2 = rMY2
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_rmass2_id


! **********************************************************************
function get_mass2_idarr(mids) result(m2)
! **********************************************************************
! m2   = complex mass squared array
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  integer, dimension(:),        intent(in) :: mids
  complex(REALKIND), dimension(size(mids)) :: m2
  integer :: i

  do i = 1, size(mids)
    m2(i) = get_mass2_id(mids(i))
  end do
end function get_mass2_idarr


! **********************************************************************
function get_rmass2_idarr(mids) result(m2)
! **********************************************************************
! m2   = real mass squared array
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND
  integer, dimension(:),     intent(in) :: mids
  real(REALKIND), dimension(size(mids)) :: m2
  integer :: i

  do i = 1, size(mids)
    m2(i) = get_rmass2_id(mids(i))
  end do
end function get_rmass2_idarr


! **********************************************************************
subroutine LC2Std_Rep(L,P)
! **********************************************************************
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: L(1:4)
  real(REALKIND),    intent(out) :: P(0:3)
  P(0) =  real(L(1)+L(2))*0.5_/**/REALKIND
  P(1) = -real(L(3)+L(4))*0.5_/**/REALKIND
  P(2) = aimag(L(4)-L(3))*0.5_/**/REALKIND
  P(3) =  real(L(2)-L(1))*0.5_/**/REALKIND
end subroutine LC2Std_Rep


! **********************************************************************
subroutine LC2Std_Rep_D(L,P)
! **********************************************************************
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! P(0:3)   here is a double-realkind vector
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: L(1:4)
  real(DREALKIND),    intent(out) :: P(0:3)
  P(0) =  real(L(1)+L(2))*0.5_/**/REALKIND
  P(1) = -real(L(3)+L(4))*0.5_/**/REALKIND
  P(2) = aimag(L(4)-L(3))*0.5_/**/REALKIND
  P(3) =  real(L(2)-L(1))*0.5_/**/REALKIND
end subroutine LC2Std_Rep_D


! **********************************************************************
subroutine LC2Std_Rep_cmplx(L,P)
! **********************************************************************
! complex version
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND, only: CI
  implicit none
  complex(REALKIND), intent(in)  :: L(1:4)
  complex(REALKIND), intent(out) :: P(0:3)
  P(0) =  (L(1)+L(2))*0.5_/**/REALKIND
  P(1) = -(L(3)+L(4))*0.5_/**/REALKIND
  P(2) = -CI*(L(4)-L(3))*0.5_/**/REALKIND
  P(3) =  (L(2)-L(1))*0.5_/**/REALKIND
end subroutine LC2Std_Rep_cmplx


! **********************************************************************
function cont_L_cmplx(A)
! Contraction of a complex Lorentz vector in standard representation with itself
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_L_cmplx
  complex(REALKIND), intent(in) :: A(0:3)
  cont_L_cmplx = A(0)*A(0) - A(1)*A(1) - A(2)*A(2) - A(3)*A(3)
end function cont_L_cmplx


! **************************************************************************
function cont_LC_cntrv(V1,V2)
! Contraction of contravariant Lorentz vectors in LightCone representation.
! Contraction of V1 with itself or with a second vector V2
! **************************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND) :: cont_LC_cntrv
  complex(REALKIND), intent(in) :: V1(1:4)
  complex(REALKIND), optional, intent(in) :: V2(1:4)
  if(present(V2)) then
    cont_LC_cntrv = (V1(1)*V2(2) + V2(1)*V1(2) - V1(3)*V2(4) - V2(3)*V1(4))/2
  else
    cont_LC_cntrv = V1(1)*V1(2) - V1(3)*V1(4)
  end if
end function cont_LC_cntrv


#ifdef PRECISION_dp
#ifdef USE_RAMBO
! *********************************************************************
subroutine rambo(sqrt_s, m_ex, p_rambo)
! *********************************************************************
! Calls Rambo for 2 -> n-2 or 1 -> n-1 PS
! *********************************************************************
! sqrt_s         = total cms energy
! m_ex(n)        = external particle masses
! p_rambo(0:3,n) = momenta, n = 1,2 incoming; n = 3,..,n outgoing
! *********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_external_decl_/**/DREALKIND, only: n_scatt
  use ol_generic, only: to_string
  use ol_debug, only: ol_fatal
  implicit none
  real(REALKIND), intent(in)  :: sqrt_s, m_ex(:)
  real(REALKIND), intent(out) :: p_rambo(0:3,size(m_ex))

  if (n_scatt == 2) then
    call rambo_2scatt(sqrt_s, m_ex, p_rambo)
  else if (n_scatt == 1) then
    call rambo_decay(sqrt_s, m_ex, p_rambo)
  else
    call ol_fatal("phase-space not available for scattering of: " // to_string(n_scatt) // &
            &     " -> " // to_string(size(m_ex)-n_scatt))
  end if

end subroutine rambo


subroutine rand_sphere(r, v)
  ! Random point, uniformly distributed on the surface of a sphere with radius r:
  ! z/r in [-1,1], phi in [0,2*pi),
  ! x=sqrt(r^2-z^2)*cos(phi), y=sqrt(r^2-z^2)*sin(phi)
  use ol_parameters_decl_/**/REALKIND, only: pi
  use ol_ramboX, only: rans
  implicit none
  real(REALKIND), intent(in) :: r
  real(REALKIND), intent(out) :: v(3)
  real(REALKIND) :: u, phi, rho
  call rans(u)
  v(3) = r*(2*u-1)
  rho = sqrt(r**2-v(3)**2)
  call rans(phi)
  phi = 2*pi*phi
  v(1) = rho*cos(phi)
  v(2) = rho*sin(phi)
end subroutine rand_sphere


subroutine decay3(E_in, m, psp)
  ! Random phase space point for a 1->2 decay with energy E_in.
  ! Set up the decay in the centre of mass system, then apply a Lorentz boost.
  ! Note that the points will not be uniformly distributed if E_in /= m(1).
  use ol_parameters_decl_/**/DREALKIND, only: psp_tolerance
  use ol_debug, only: ol_error, ol_msg, ol_fatal
  implicit none
  real(REALKIND), intent(in) :: E_in, m(3)
  real(REALKIND), intent(out) :: psp(0:3,3)
  real(REALKIND) :: E1, m12, m22, m32, E2a, E3a, p1_abs, p2a_abs, p2a(3), gamma, beta(3), beta_p2a
  if (m(1) <= m(2) + m(3)) then
    call ol_error(2,"3-particle interaction:")
    call ol_msg("mass condition m1+m2>m3 (production) or m1>m2+m3 (decay) not satisfied.")
    call ol_fatal()
  end if
  if (abs(E_in/m(1)-1) < psp_tolerance) then
    E1 = m(1)
  else if (E_in < m(1)) then
    call ol_fatal("3-particle interaction energy too low.")
    return
  else
    E1 = E_in
  end if
  m12 = m(1)**2
  m22 = m(2)**2
  m32 = m(3)**2
  ! centre of mass system: (m1,0) = (E2,p2a) + (E3,-p2a)
  E2a = (m12 + m22 - m32)/(2*m(1))
  E3a = (m12 - m22 + m32)/(2*m(1))
!   p2a_abs = sqrt((m12**2 + m22**2 + m32**2 - 2*m12*m22 - 2*m12*m32 - 2*m22*m32)/(4*m12))
  p2a_abs = sqrt(E2a**2-m22)
  call rand_sphere(p2a_abs, p2a)
  ! Lorentz boost (Ep,pp) = L.(E,p) such that L.(m1,0) = (sqrt,p1)
  ! L_00 = gamma, L_0i = L_i0 = -gamma*beta_i, L_ij = delta_ij + (gamma-1)*beta_i*beta_j/beta^2
  ! Ep = gamma*(E-beta.p)
  ! pp = p + ((gamma-1)*beta.p/beta^2 - E*gamma)*beta
  ! p1 and boost parameters
  if (E1 == m(1)) then
    psp(0,1) = m(1)
    psp(1:3,1) = 0
    psp(0,2) = E2a
    psp(0,3) = E3a
    psp(1:3,2) =  p2a
    psp(1:3,3) = -p2a
  else
    p1_abs = sqrt(E1**2-m12)
    psp(0,1) = E1
    call rand_sphere(p1_abs, psp(1:3,1))
    gamma = E1/m(1)
    beta = -psp(1:3,1)/E1
    ! p2 and p3
    beta_p2a = sum(beta*p2a)
    psp(0,2) = gamma*(E2a-beta_p2a)
    psp(0,3) = gamma*(E3a+beta_p2a)
    p2a = p2a + ((gamma-1)*beta_p2a/sum(beta**2)) * beta
    psp(1:3,2) =  p2a - (gamma*E2a)*beta
    psp(1:3,3) = -p2a - (gamma*E3a)*beta
  end if
end subroutine decay3


! *********************************************************************
subroutine rambo_2scatt(sqrt_s, m_ex, p_rambo)
! *********************************************************************
! Generate 2 -> n-2 PS-point with P(1) + P(2) = P(3) + ... + P(n)
! Apply cleaning procedure to get full numerical precision
! *********************************************************************
! sqrt_s         = total cms energy
! m_ex(n)        = external particle masses
! p_rambo(0:3,n) = momenta, n = 1,2 incoming; n = 3,..,n outgoing
! beam momenta:
!   p_rambo(0:3,1) = (EA,0,0,kA)
!   p_rambo(0:3,2) = (EB,0,0,kB)
!   EA + EB        = sqrt_s
!   kA + kB        = 0
! *********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal
  use ol_ramboX, only: rambo0 => rambo
  implicit none
  real(REALKIND), intent(in)  :: sqrt_s, m_ex(:)
  real(REALKIND), intent(out) :: p_rambo(0:3,size(m_ex))
  real(REALKIND) :: E, MA2, MB2, dEAB
  real(REALKIND) :: p_scatt(4,size(m_ex)-2), wgt
  integer :: n
  n = size(m_ex)
  if (n >= 4) then
    E = sqrt_s*0.5_/**/REALKIND
    ! beam momenta
    if((m_ex(1) == 0) .and. (m_ex(2) == 0)) then
      p_rambo(0,1) =  E
      p_rambo(1,1) =  0
      p_rambo(2,1) =  0
      p_rambo(3,1) =  E
      p_rambo(0,2) =  E
      p_rambo(1,2) =  0
      p_rambo(2,2) =  0
      p_rambo(3,2) = -E
    else
      MA2  = m_ex(1)*m_ex(1)
      MB2  = m_ex(2)*m_ex(2)
      dEAB = (MA2 - MB2) / (2*sqrt_s)
      p_rambo(0,1) =  E + dEAB
      p_rambo(1,1) =  0
      p_rambo(2,1) =  0
      p_rambo(3,1) =  sqrt(p_rambo(0,1)**2-MA2)
      p_rambo(0,2) =  E - dEAB
      p_rambo(1,2) =  0
      p_rambo(2,2) =  0
      p_rambo(3,2) = -p_rambo(3,1)
    end if
    call rambo0(n-2, sqrt_s, m_ex(3:n), p_scatt, wgt)
    p_rambo(  0,3:n) = p_scatt(  4,1:n-2)
    p_rambo(1:3,3:n) = p_scatt(1:3,1:n-2)
  else if (n == 3) then
    ! reverse 1->2 decay
    call decay3(sqrt_s, [m_ex(3),m_ex(1),m_ex(2)], p_rambo)
    p_scatt(:,1) = p_rambo(:,1)
    p_rambo(:,1) = p_rambo(:,2)
    p_rambo(:,2) = p_rambo(:,3)
    p_rambo(:,3) = p_scatt(:,1)
  else
    call ol_fatal("2->0 scattering not possible -- use decay instead.")
  end if

end subroutine rambo_2scatt


! *********************************************************************
subroutine rambo_decay(sqrt_s, m_ex, p_rambo)
! *********************************************************************
! Generate 1 -> n-1 PS-point with P(1) = P(2) + ... + P(n)
! Apply cleaning procedure to get full numerical precision
! *********************************************************************
! sqrt_s         = total cms energy
! m_ex(n)        = external particle masses
! p_rambo(0:3,n) = momenta, n = 2,..,n outgoing
! *********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal, ol_msg
  use ol_parameters_decl_/**/DREALKIND, only: psp_tolerance
  use ol_ramboX, only: rambo0 => rambo
  implicit none
  real(REALKIND), intent(in)  :: sqrt_s, m_ex(:)
  real(REALKIND), intent(out) :: p_rambo(0:3,size(m_ex))
  real(REALKIND) :: p_scatt(4,size(m_ex)-1), wgt, p1
  integer :: n, k
  n = size(m_ex)
  if (sqrt_s < m_ex(1)) then
    call ol_fatal("energy in decay lower than mass.")
  end if
  if (n >= 3) then
    if( m_ex(1) == 0 ) then
      call ol_msg("Warning: decay of massless particle!")
    else
      p_rambo(0,1) =  sqrt_s
      p_rambo(1,1) =  0
      p_rambo(2,1) =  0
      p_rambo(3,1) =  sqrt(sqrt_s*sqrt_s-m_ex(1)*m_ex(1))
    end if
    call rambo0(n-1, sqrt_s, m_ex(2:n), p_scatt, wgt)
    p_rambo(  0,2:n) = p_scatt(  4,1:n-1)
    p_rambo(1:3,2:n) = p_scatt(1:3,1:n-1)
  else
    ! particle with energy sqrt_s coming from a random direction
    if (abs(m_ex(1)-m_ex(2))/sqrt_s > psp_tolerance) then
      call ol_fatal("two particle processes require external particles with equal mass.")
    end if
    p1 = sqrt(sqrt_s**2-m_ex(1)**2)
    p_rambo(0,1) = sqrt_s
    call rand_sphere(p1, p_rambo(1:3,1))
    p_rambo(:,2) = p_rambo(:,1)
  end if
end subroutine rambo_decay


subroutine rambo_c(sqrt_s, m_ex, n, p_rambo) bind(c,name="ol_rambo")
  use KIND_TYPES, only: REALKIND
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), intent(in)  :: n
  real(c_double), intent(out) :: p_rambo(0:3,n)
  real(c_double), intent(in)  :: sqrt_s, m_ex(n)
  real(REALKIND) :: f_p_rambo(0:3,n)
  real(REALKIND) :: f_sqrt_s, f_m_ex(n)
  f_m_ex = m_ex
  f_sqrt_s = sqrt_s
  call rambo(f_sqrt_s, f_m_ex, f_p_rambo)
  p_rambo = f_p_rambo
end subroutine rambo_c

! #ifdef USE_RAMBO
#else

subroutine rambo(sqrt_s, m_ex, p_rambo)
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal
  implicit none
  real(REALKIND), intent(in)  :: sqrt_s, m_ex(:)
  real(REALKIND), intent(out) :: p_rambo(0:3,size(m_ex))
  p_rambo = 0 ! prevent compiler warning
  call ol_fatal('Rambo is not available.')
end subroutine rambo

subroutine rambo_c(sqrt_s, m_ex, n, p_rambo) bind(c,name="ol_rambo")
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), intent(in)  :: n
  real(c_double), intent(out) :: p_rambo(0:3,n)
  real(c_double), intent(in)  :: sqrt_s, m_ex(n)
  p_rambo = 0 ! prevent compiler warning
  call ol_fatal('Rambo is not available.')
end subroutine rambo_c

! #ifdef USE_RAMBO
#endif
! #ifdef PRECISION_dp
#endif


! **********************************************************************
subroutine clean_mom_in(P_in, m_ext2, P, n)
! remove numerical inaccuracies in n -> 0 point (momenta all incoming)
! P_in(0:3,n) = original PS point in double precision
! m_ext2      = squared external masses
! P(0:3,n)    = improved PS point in working precision
! 3-momentum of n-th particle recomputed with momentum conservation;
! all energies recomputed with on-shell condition;
! 3-momenta (and energies) of incoming particles shifted by +/- eps*p1,
! so that energy conservation is fulfilled up to terms of O(eps^3)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_msg
  use ol_generic, only: to_string
  use ol_parameters_decl_/**/DREALKIND, only: psp_tolerance
  implicit none
  real(REALKIND), intent(in)  :: P_in(0:3,n)
  integer,        intent(in)  :: n
  real(REALKIND), intent(in)  :: m_ext2(n)
  real(REALKIND), intent(out) :: P(0:3,n)
  real(REALKIND)  :: E_ref, P0(n), P2(n)
  real(REALKIND)  :: E0(2), E1(2), E2(2)
  real(REALKIND)  :: E0_tot, E1_tot, E2_tot
  real(REALKIND)  :: eps1, eps
  integer         :: nex, i, pout_max_pos, pout_max_pos_arr(1)
  real            :: prec

  P = P_in

!   call dirty_mom(P_in, P, n, 9)

  ! print momena before cleaning
!   do i = 1, n
!     write(*,*) P(:,i), (P(0,i)**2-P(1,i)**2-P(2,i)**2-P(3,i)**2)-m_ext2(i)
!   end do
!   write(*,*) sum(P(0,:)), sum(P(1,:)), sum(P(2,:)), sum(P(3,:))
!   write(*,*)

  E_ref = 0.5_/**/REALKIND * sum(abs(P(0,:)))
  ! check momentum conservation
  do i = 0, 3
    prec = abs(sum(P(i,:)))/E_ref
    if (prec > psp_tolerance) then
      call ol_msg("=== WARNING ===")
      call ol_msg("OpenLoops subroutine clean_mom: inconsistent phase space point.")
      call ol_msg("Momentum conservation is only satisfied to " // to_string(-log10(prec)) // "digits.")
      call ol_msg("===============")
    end if
  end do

  ! check on-shell conditions
  do nex = 1, n
    P2(nex) = P(1,nex)*P(1,nex) + P(2,nex)*P(2,nex) + P(3,nex)*P(3,nex)
    P0(nex) = sign(sqrt(P2(nex) + m_ext2(nex)), P(0,nex))
    prec = abs(P(0,nex)-P0(nex))/E_ref
    if(prec > psp_tolerance) then
      call ol_msg("=== WARNING ===")
      call ol_msg("OpenLoops subroutine clean_mom: inconsistent phase space point.")
      call ol_msg("On-shell condition is only satisfied to " // to_string(-log10(prec)) // "digits.")
      call ol_msg("===============")
    end if
  end do

  ! position of the outgoing momentum with the largest energy
  pout_max_pos_arr = maxloc(abs(P(0,3:)))
  pout_max_pos = 2 + pout_max_pos_arr(1)
  ! fix 3-momentum by momentum conservation
  do i = 1, 3
    P(i,pout_max_pos) = P(i,pout_max_pos) - sum(P(i,:))
  end do

  ! enforce on-shell conditions
  P0(pout_max_pos) = sign(sqrt(P(1,pout_max_pos)**2 + P(2,pout_max_pos)**2 + P(3,pout_max_pos)**2 + m_ext2(pout_max_pos)), &
                        & P(0,pout_max_pos))
  P(0,:) = P0

  E0_tot = sum(P(0,:))

  ! lowest-order energy terms
  E0(1)  = P(0,1)
  E0(2)  = P(0,2)

  ! 1st order energy coefficients
  E1(1)  = P2(1)/E0(1)
  E1(2)  = -(P(1,1)*P(1,2)+P(2,1)*P(2,2)+P(3,1)*P(3,2))/E0(2)
  E1_tot = E1(1)+E1(2)

  ! 2nd order energy coefficients
  ! for quad-precision applications it is recommended to use the formulas w.o.
  ! beam-alignement instabilities both for E2(1) and E2(2)
!   E2(1)  = (P2(1)-E1(1)**2)/(2*E0(1))   ! direct formula
  E2(1)  = P2(1)*m_ext2(1)/(2*E0(1)**3) ! default = equivalent fomula w.o. beam-alignment instabilities
  E2(2)  = (P2(1)-E1(2)**2)/(2*E0(2))   ! default = direct formula
!   E2(2)  = 0.5_/**/REALKIND/E0(2)**3*(P2(1)*(m_ext2(2) + P(1,2)**2 + P(2,2)**2) & ! equivalent formula w.o. beam-alignment instabilities
!          + (P(1,1)**2+P(2,1)**2)*P(3,2)**2 &
!          - (P(1,1)*P(1,2)+P(2,1)*P(2,2))*(P(1,1)*P(1,2)+P(2,1)*P(2,2)+2*P(3,1)*P(3,2)))
  E2_tot = E2(1) + E2(2)

  ! 1st order shift
  eps1 = -E0_tot/E1_tot
  ! 2nd order shift
  eps  = eps1*(1-eps1*E2_tot/E1_tot)

  ! shifted IS momenta
  do i = 1, 3
    P(i,2) = P(i,2) - P(i,1)*eps
    P(i,1) = P(i,1) + P(i,1)*eps
  end do

  ! shifted IS energies
  do nex = 1, 2
    if(m_ext2(nex) == 0 .and. P(1,nex) == 0 .and. P(2,nex) == 0) then
      ! exact formula for m = 0 along beam-axis
      P(0,nex) = sign(abs(P(3,nex)), real(P_in(0,nex), REALKIND))
    else
      P(0,nex) = E0(nex) + E1(nex)*eps + E2(nex)*eps**2
    end if
  end do

  ! print momenta after cleaning
!    do i = 1, n
!      write(*,*) P(:,i), (P(0,i)**2-P(1,i)**2-P(2,i)**2-P(3,i)**2)-m_ext2(i)
!    end do
!    write(*,*) sum(P(0,:)), sum(P(1,:)), sum(P(2,:)), sum(P(3,:))
!   write(*,*)

end subroutine clean_mom_in


! **********************************************************************
subroutine clean_mom_scatt(P_in, m_ext2, P, n)
! same as clean_mom_in but for 2-> n-2 PS point
! This routine is not used internally.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  real(REALKIND), intent(in)  :: P_in(0:3,n)
  integer,        intent(in)  :: n
  real(REALKIND), intent(in)  :: m_ext2(n)
  real(REALKIND), intent(out) :: P(0:3,n)
  real(REALKIND) :: Q_in(0:3,n)
  real(REALKIND) :: Q(0:3,n)

  Q_in(:,1:2) = P_in(:,1:2)
  Q_in(:,3:) = -P_in(:,3:)

  call clean_mom_in(Q_in, m_ext2, Q, n)

  P(:,1:2) = Q(:,1:2)
  P(:,3:) = -Q(:,3:)

end subroutine clean_mom_scatt


#ifdef USE_RAMBO
! **********************************************************************
subroutine dirty_mom(P_in, P, n, DIG)
! introduces random noise in every PS component after DIG digits
! (to test PS-point cleaning).
! Needs RAMBO subroutine rans().
! **********************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_ramboX, only: rans
  implicit none
  integer,        intent(in)  :: n, DIG
  real(REALKIND), intent(in)  :: P_in(0:3,n)
  real(REALKIND), intent(out) :: P(0:3,n)
  real(DREALKIND) :: x, shift
  integer        :: nex, i

  shift = 10._/**/DREALKIND**(-DIG)

  do nex = 1, n
    do i = 0, 3
      call rans(x)
      P(i,nex) = P_in(i,nex)*(1+(x-0.5_/**/REALKIND)*shift)
    end do
  end do

end subroutine dirty_mom
! #ifdef USE_RAMBO
#else
subroutine dirty_mom(P_in,P, n, DIG)
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal
  implicit none
  integer,        intent(in)  :: n, DIG
  real(REALKIND), intent(in)  :: P_in(0:3,n)
  real(REALKIND), intent(out) :: P(0:3,n)
  call ol_fatal('dirty_mom() requires Rambo.')
end subroutine dirty_mom
! #ifdef USE_RAMBO
#endif

subroutine init_kinematics_mexpl(P_scatt, m_ext2, P_in_clean, perm_inv, n)
  use ol_external_decl_/**/DREALKIND, only: init_qp
  integer,           intent(in)  :: n
  real(DREALKIND),   intent(in)  :: P_scatt(0:3,n)
  real(REALKIND),    intent(in)  :: m_ext2(n)
  real(REALKIND),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)
  real(REALKIND) :: P(0:3,n)

#ifdef PRECISION_dp
  init_qp = .false.
#else
  init_qp = .true.
#endif

  call conv_mom_scatt2in(P_scatt, m_ext2, P, perm_inv, n)
  call internal_momenta(P, n)

end subroutine init_kinematics_mexpl

subroutine init_kinematics_mids(P_scatt, m_ext2, P_in_clean, perm_inv, n, qp_kinematics)
#ifdef PRECISION_dp
  use ol_parameters_decl_/**/REALKIND, only: hp_switch,hp_qp_kinematics_init_mode,&
                                             hp_nsi,hp_nsi_qp,hp_ndrs,hp_ndrs_qp, &
                                             hp_nred,hp_nred_qp,hp_max_err
  use ol_kinematics_/**/QREALKIND, only: init_kinematics_mids_qp=>init_kinematics_mids
#endif
  use ol_external_decl_/**/DREALKIND, only: init_qp
  integer,         intent(in)  :: n
  real(DREALKIND), intent(in)  :: P_scatt(0:3,n)
  integer,         intent(in)  :: m_ext2(n)
  real(REALKIND),  intent(inout) :: P_in_clean(0:3,n)
  integer,         intent(in)  :: perm_inv(n)
  logical,         intent(in)  :: qp_kinematics
#ifdef PRECISION_dp
  real(QREALKIND) :: P_qp(0:3,n)
  integer :: i1

  init_qp = .false.
#else
  init_qp = .true.
#endif

  call conv_mom_scatt2in_mids(P_scatt, m_ext2, P_in_clean, perm_inv, n)
#ifdef PRECISION_dp
  if (hp_switch .eq. 1) then
    hp_nsi = 0
    hp_nsi_qp = 0
    hp_ndrs = 0
    hp_ndrs_qp = 0
    hp_nred = 0
    hp_nred_qp = 0
    hp_max_err = 0._/**/REALKIND
  end if
  ! initialize qp kinematics for hybrid mode
  if ((hp_qp_kinematics_init_mode .eq. 0 .and. qp_kinematics) .and. &
      (hp_switch .eq. 1)) then
    call init_kinematics_mids_qp(P_scatt, m_ext2, P_qp, perm_inv, n, &
      (hp_qp_kinematics_init_mode .eq. 0 .and. qp_kinematics))
  end if
  call internal_momenta_six(P_in_clean, n, m_ext2, &
       (hp_qp_kinematics_init_mode .eq. 0 .and. qp_kinematics))
#else
  call internal_momenta_six(P_in_clean, n, m_ext2, .false.)
#endif

end subroutine init_kinematics_mids

#ifdef PRECISION_dp
subroutine init_qp_kinematics()
  use ol_external_decl_/**/REALKIND, only: nParticles,M_ex,init_qp
  use KIND_TYPES, only: QREALKIND
  use ol_kinematics_/**/QREALKIND, only: conv_mom_scatt2in_cache, &
    internal_momenta_six_qp=>internal_momenta_six
  real(QREALKIND) :: P_in_clean_qp(0:3,nParticles)
  real(QREALKIND) :: P_qp(0:3,nParticles)
  call conv_mom_scatt2in_cache(P_in_clean_qp, nParticles)
  call internal_momenta_six_qp(P_in_clean_qp, nParticles, M_ex, .false.)
  init_qp = .true.
  call internal_momenta_six(real(P_in_clean_qp,kind=REALKIND), nParticles, M_ex, .true.)
end subroutine init_qp_kinematics
#endif

! **********************************************************************
subroutine conv_mom_scatt2in_mexpl(P_scatt, m_ext2, P_in_clean, perm_inv, n)
! Keep incoming momenta and reverse outgoing momenta.
! Apply phase space cleaning and crossing.
! Note: calls init_kin_arrays -> allocation of kinematic arrays
! ToDo: cleaning for n_scatt /= 2
! **********************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_external_decl_/**/REALKIND, only: nParticles, P_ex, inverse_crossing
  use ol_external_decl_/**/DREALKIND, only: n_scatt
  use ol_parameters_decl_/**/REALKIND, only: scalefactor
  use ol_parameters_init_/**/REALKIND, only: init_kin_arrays
!#ifdef PRECISION_dp
!  use ol_parameters_init_qp, only: init_kin_arrays_qq => init_kin_arrays
!#endif
  integer,           intent(in)  :: n
  real(DREALKIND),   intent(in)  :: P_scatt(0:3,n)
  real(REALKIND),    intent(in)  :: m_ext2(n)
  real(REALKIND),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)
  real(REALKIND) :: P_in(0:3,n)
  real(REALKIND) :: P_clean(0:3,n), m_ext2_perm(n)
  integer        :: k
  nParticles = n
  call init_kin_arrays(nParticles)
  P_ex(:,1:n) = P_scatt
  inverse_crossing(1:n) = perm_inv
  do k = 1, n
    m_ext2_perm(perm_inv(k)) = m_ext2(k)
  end do
  P_in(:,1:n_scatt) =   scalefactor * P_scatt(:,1:n_scatt)
  P_in(:,n_scatt+1:)  = - scalefactor * P_scatt(:,n_scatt+1:)
  if (n_scatt == 2 .and. n > 2) then
    ! Clean momenta to get full numerical precision.
    ! Do the cleaning in the original permutation where the first two momenta are incoming.
    ! Otherwise the beam alignment (zero components) might be spoiled by the cleaning.
    call clean_mom_in(P_in, m_ext2_perm, P_clean, n)
  else
    P_clean = P_in
  end if
  do k = 1, n
    P_in_clean(:,k) = P_clean(:,perm_inv(k))
  end do
end subroutine conv_mom_scatt2in_mexpl

subroutine conv_mom_scatt2in_mids(P_scatt, m_ext, P_in_clean, perm_inv, n)
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_external_decl_/**/REALKIND, only: M_ex
  use ol_external_decl_/**/DREALKIND, only: n_scatt
  use ol_parameters_decl_/**/REALKIND, only: scalefactor
  use ol_parameters_init_/**/REALKIND, only: init_kin_arrays
!#ifdef PRECISION_dp
!  use ol_parameters_init_qp, only: init_kin_arrays_qq => init_kin_arrays
!#endif
  implicit none
  integer,           intent(in)  :: n
  real(DREALKIND),   intent(in)  :: P_scatt(0:3,n)
  integer,           intent(in)  :: m_ext(n)
  real(REALKIND),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)

  call conv_mom_scatt2in_mexpl(P_scatt,get_rmass2(m_ext),P_in_clean,perm_inv,n)
  M_ex(1:n) = m_ext
end subroutine conv_mom_scatt2in_mids

#ifndef PRECISION_dp
! **********************************************************************
subroutine conv_mom_scatt2in_cache(P_in_clean,n)
! ToDo: Initializes the momenta in QP from the DP cache
! **********************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_external_decl_/**/DREALKIND, only: nParticles, P_ex_dp=>P_ex, M_ex_dp=>M_ex, inverse_crossing_dp=>inverse_crossing
  use ol_external_decl_/**/REALKIND, only: P_ex
  use ol_external_decl_/**/DREALKIND, only: n_scatt
  use ol_parameters_decl_/**/REALKIND, only: scalefactor
  use ol_parameters_init_/**/REALKIND, only: init_kin_arrays
  integer,           intent(in)  :: n
  real(REALKIND),    intent(out) :: P_in_clean(0:3,n)
  real(REALKIND) :: P_in(0:3,n)
  real(REALKIND) :: P_clean(0:3,n), m_ext2_perm(n)
  integer        :: k
  call init_kin_arrays(n)
  P_ex(:,1:n) = P_ex_dp(:,1:n)
  do k = 1, n
    m_ext2_perm(inverse_crossing_dp(k)) = get_rmass2(M_ex_dp(k))
  end do
  P_in(:,1:n_scatt) =   scalefactor * P_ex(:,1:n_scatt)
  P_in(:,n_scatt+1:)  = - scalefactor * P_ex(:,n_scatt+1:)
  if (n_scatt == 2 .and. n > 2) then
    ! Clean momenta to get full numerical precision.
    ! Do the cleaning in the original permutation where the first two momenta are incoming.
    ! Otherwise the beam alignment (zero components) might be spoiled by the cleaning.
    call clean_mom_in(P_in, m_ext2_perm, P_clean, n)
  else
    P_clean = P_in
  end if
  do k = 1, n
    P_in_clean(:,k) = P_clean(:,inverse_crossing_dp(k))
  end do
end subroutine conv_mom_scatt2in_cache

#endif

! **********************************************************************
subroutine conv_mom_os(P_decay, P_in, n)
! reverse decay products
! ToDo: Apply phase space cleaning
! **********************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_parameters_decl_/**/DREALKIND, only: scalefactor
  implicit none
  integer,         intent(in)  :: n
  real(DREALKIND), intent(in)  :: P_decay(0:3,n)
  real(REALKIND),  intent(out) :: P_in(0:3,n)
  integer         :: k

  P_in  = - scalefactor * P_decay

end subroutine conv_mom_os



#ifdef PRECISION_dp

! **********************************************************************
subroutine write_INmom(P_ex, n, unit)
! **********************************************************************
! Write the information on the four momenta.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none

  real(REALKIND), intent(in) :: P_ex(0:3,n)
  integer,        intent(in) :: n, unit
  integer        :: i, k
  real(REALKIND) :: Ptot(0:3), Pabs(0:3), PR(0:3)

  do i = 0, 3
    Ptot(i) = P_ex(i,1)
    Pabs(i) = abs(P_ex(i,1))
    do k = 2, n
      Ptot(i) = Ptot(i) + P_ex(i,K)
      Pabs(i) = Pabs(i)+abs(P_ex(i,K))
    end do
    PR(i) = Ptot(i)/Pabs(i)
  end do

  write (unit,*) "------------------------------------", "-----------------------------------------"
  write (unit,*) " ", n, " -> 0  Phase space point:"
  write (unit,*) "------------------------------------", "-----------------------------------------"
  write (unit,*) "n        E             px             ", "py              pz               m "
  do i = 1, n
    write (unit,'(i2,1x,5e15.7)') i, P_ex(0,i), P_ex(1,i), P_ex(2,i), P_ex(3,i), sqrt(abs(cont_LL(P_ex(0,i),P_ex(0,i))))
  end do
  write (unit,*) "------------------------------------", "-----------------------------------------"
  write (unit,'(A2,1x,4e15.7)') 'T', PR(0), PR(1), PR(2), PR(3)
  write (unit,*) "------------------------------------", "-----------------------------------------"
  write (unit,*)
  contains
  function cont_LL(A, B)
    ! Contraction of (real) Lorentz vectors in standard representation.
    ! Do not use the contractions module to avoid dependency cycles.
    use KIND_TYPES, only: REALKIND
    implicit none
    real(REALKIND) :: cont_LL
    real(REALKIND), intent(in) :: A(0:3), B(0:3)
    cont_LL = A(0)*B(0) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
  end function cont_LL
end subroutine write_INmom

! #ifdef PRECISION_dp
#endif



! **********************************************************************
subroutine internal_momenta_std(P, Npart)
! **********************************************************************
! P(0:3,Npart) = external real-valued four-momenta (standard representation)
! Npart        = total (in & out) external particle number
! Q(1:5,1:2^Npart-1) = internal four-momenta in light-cone representation;
!                      the fifth component is the complex-valued squared momentum.
! Numbering of internal momenta:
!   Sum_i s(i)*P(i) => Q(Sum_i s(i)*2^(i-1)), s(i) = 0, 1
!   so that Q(J1) + Q(J2) = Q(J1+J2)
! QInvariantsMatrix(i,j) = (p_i+p_j)^2 for i /= j, otherwise undefined.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: Q, QInvariantsMatrix, &
                                          collconf, softconf
  implicit none

  real(REALKIND), intent(in) :: P(0:3,Npart)
  integer,        intent(in)  :: Npart
  integer :: i, j
  integer :: Jmax
  integer :: i1, i2 ! conventional particle numbers
  integer :: l1, l2 ! individual 2^(i-1) particles numbers
  integer :: s1, s2 ! sums of 2^(i-1) particle numbers
  integer :: r1, r2 ! inverse of s1, s2, ...

  collconf = 0
  softconf = 0
  i = 2**Npart - 2

  Q(:,i+1) = 0

  if (Npart == 2) then
    call Std2LC_Rep(P(:,1),Q(1:4,1))
    Q(1:4,2) = - Q(1:4,1)
    Q(5,1) = Q(1,1)*Q(2,1) - Q(3,1)*Q(4,1)
    Q(5,2) = Q(5,1)
  else if (Npart == 3) then
    call Std2LC_Rep(P(:,1),Q(1:4,1))
    call Std2LC_Rep(P(:,2),Q(1:4,2))
    Q(1:4,3) = Q(1:4,1) + Q(1:4,2)
    Q(1:4,4) = - Q(1:4,3)
    Q(1:4,5) = - Q(1:4,2)
    Q(1:4,6) = - Q(1:4,1)
    Q(5,1) = Q(1,1)*Q(2,1) - Q(3,1)*Q(4,1)
    Q(5,2) = Q(1,2)*Q(2,2) - Q(3,2)*Q(4,2)
    Q(5,3) = Q(1,3)*Q(2,3) - Q(3,3)*Q(4,3)
    Q(5,4) = Q(5,3)
    Q(5,5) = Q(5,2)
    Q(5,6) = Q(5,1)
  else
    call intmom(P, Npart, i)
  end if

  do i = 1, Npart
    ! QInvariantsMatrix(i,i) = m_ex2(i)
    do j = i + 1, Npart
      QInvariantsMatrix(i,j) = Q(5,2**(i-1)+2**(j-1))
      QInvariantsMatrix(j,i) = QInvariantsMatrix(i,j)
    end do
  end do

end subroutine internal_momenta_std


! **********************************************************************
subroutine internal_momenta_six(P, Npart, ext_masses, use_qp_kinematics)
! **********************************************************************
! P(0:3,Npart)       = external real-valued four-momenta
!                      (standard-representation)
! Npart              = number of external particle (in & out)
! L(1:6,1:2^Npart-1) = internal four-momenta in light-cone representation;
!                      the fifth component is the real-valued sum of external
!                      masses while the
!                      the sixth component is the sum of the scalar-products
! Numbering of internal momenta:
!   Sum_i s(i)*P(i) => Q(Sum_i s(i)*2^(i-1)), s(i) = 0, 1
!   so that Q(J1) + Q(J2) = Q(J1+J2)
! QInvariantsMatrix(i,j) = (p_i+p_j)^2 for i /= j, otherwise undefined.
! invariants_mode    = flag used to decide in which way the invariants
!                      from the internal momenta are computed
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1
  use ol_momenta_decl_/**/REALKIND, only: L, collconf, softconf
#ifdef PRECISION_dp
  use ol_parameters_decl_/**/REALKIND, only: hp_switch, hp_mode, sync_qp_kinematics
  use ol_momenta_decl_/**/QREALKIND, only: L_qp=>L
#endif
  implicit none

  integer,        intent(in) :: Npart
  real(REALKIND), intent(in) :: P(0:3,Npart)
  integer,        intent(in) :: ext_masses(Npart)
  logical,        intent(in) :: use_qp_kinematics
  complex(REALKIND) :: zero = 0._/**/REALKIND
  real(REALKIND)    :: maxinv
  integer :: i, j
  integer :: Jmax
  integer :: i1, i2 ! conventional particle numbers
  integer :: l1, l2 ! individual 2^(i-1) particles numbers
  integer :: s1, s2 ! sums of 2^(i-1) particle numbers
  integer :: r1, r2 ! inverse of s1, s2, ...
  integer :: t1

  i = 2**Npart - 2
  L(:,i+1) = 0
  L(5:6,:) = zero
  collconf = 0
  softconf = 0

  if (Npart == 2) then
#ifdef PRECISION_dp
    if (use_qp_kinematics .and. hp_switch .eq. 1 .and. sync_qp_kinematics .eq. 1) then
      L(1:4,1) = L_qp(1:4,1)
    else
      call Std2LC_Rep(P(:,1),L(1:4,1))
    end if
#else
    call Std2LC_Rep(P(:,1),L(1:4,1))
#endif
    L(5,1) = get_rmass2(ext_masses(1))
    L(6,1) = zero
    L(5:6,2) = L(5:6,1)
  else if (Npart == 3) then
#ifdef PRECISION_dp
    if (use_qp_kinematics .and. hp_switch .eq. 1 .and. sync_qp_kinematics .eq. 1) then
      L(1:4,1:6) = L_qp(1:4,1:6)
    else
      call Std2LC_Rep(P(:,1),L(1:4,1))
      call Std2LC_Rep(P(:,2),L(1:4,2))
      L(1:4,3) = L(1:4,1) + L(1:4,2)
      L(1:4,4) = - L(1:4,3)
      L(1:4,5) = - L(1:4,2)
      L(1:4,6) = - L(1:4,1)
    end if
#else
    call Std2LC_Rep(P(:,1),L(1:4,1))
    call Std2LC_Rep(P(:,2),L(1:4,2))
    L(1:4,3) = L(1:4,1) + L(1:4,2)
    L(1:4,4) = - L(1:4,3)
    L(1:4,5) = - L(1:4,2)
    L(1:4,6) = - L(1:4,1)
#endif
    L(5,1) = get_rmass2(ext_masses(1))
    L(5,2) = get_rmass2(ext_masses(2))
    L(6,1) = zero
    L(6,2) = zero
    L(5,3) = get_rmass2(ext_masses(1)) + get_rmass2(ext_masses(2))
#ifdef PRECISION_dp
    ! overwrite dp invariants by qp ones
    if (use_qp_kinematics .and. hp_switch .eq. 1 .and. sync_qp_kinematics .eq. 1) then
      L(6,3) = L_qp(6,3)
    else
      L(6,3) = 2*cont_LC_cntrv(L(1:4,1),L(1:4,2))
    end if
#else
    L(6,3) = 2*cont_LC_cntrv(L(1:4,1),L(1:4,2))
#endif
    L(5:6,4) = L(5:6,3)
    L(5:6,5) = L(5:6,2)
    L(5:6,6) = L(5:6,1)
  else
    ! compute masses and scalar products
    call intmom_six(P, Npart, i, get_rmass2(ext_masses), use_qp_kinematics)

  end if

#ifdef PRECISION_dp
  maxinv = 0
  do i = 1, Npart
    do j = i + 1, Npart
      t1 = 2**(i-1)+2**(j-1)
      maxinv = max(abs(L(5,t1) + L(6,t1)), maxinv)
    end do
  end do
  if (hp_mode .ge. 2) then
    do i = 1, Npart
      if (maxval(abs(L(1:4,2**(i-1))))**2/maxinv .lt. softthres) softconf  = 2**(i-1)
    end do
  end if
#endif

end subroutine internal_momenta_six


! **********************************************************************
subroutine intmom(P_ex,Npart,Jmax)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: Q
  implicit none

  integer,        intent(in) :: Npart, Jmax
  real(REALKIND), intent(in) :: P_ex(0:3,Npart)
  integer :: A
  integer :: i1 ! conventional particle numbers
  integer :: l1 ! individual 2^(i-1) particles numbers
  integer :: s1 ! sums of 2^(i-1) particle numbers
  integer :: r1 ! inverse of s1, s2, ...

  l1 = 1 ! ext mom 1 <= i1 <= Npart
  do i1 = 1, Npart
    s1 = l1
    r1 = Jmax + 1 - s1
    call Std2LC_Rep(P_ex(0,i1),Q(1,l1))
    l1 = 2*l1
    do A = 1, 4
      Q(A,r1) = -Q(A,s1)
    end do
    call intmom_rec(Npart, Jmax, i1, s1, 1)
  end do !i1

  do s1 = 1, Jmax/2 ! squared momenta
    Q(5,s1) = Q(1,s1)*Q(2,s1) - Q(3,s1)*Q(4,s1)
    Q(5,Jmax+1-s1) = Q(5,s1)
  end do

end subroutine intmom


! **********************************************************************
subroutine intmom_six(P_ex,Npart,Jmax,ext_masses,use_qp_kinematics)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_momenta_decl_/**/REALKIND, only: L
  use ol_loop_parameters_decl_/**/REALKIND, only: zero
  implicit none

  integer,        intent(in) :: Npart, Jmax
  real(REALKIND), intent(in) :: ext_masses(Npart)
  real(REALKIND), intent(in) :: P_ex(0:3,Npart)
  logical, intent(in)        :: use_qp_kinematics
  integer :: i1 ! conventional particle numbers
  integer :: l1 ! individual 2^(i-1) particles numbers

  l1 = 1
  do i1 = 1, Npart
    L(5,l1) = ext_masses(i1)
    L(6,l1) = zero
    L(5:6,Jmax + 1 - l1) = L(5:6,l1)
    call Std2LC_Rep(P_ex(0,i1),L(1,l1))
    L(1:4,Jmax + 1 - l1) = -L(1:4,l1)
    call intmom_rec_six(Npart, Jmax, i1, l1, 1, use_qp_kinematics)
    l1 = 2*l1
  end do !i1

end subroutine intmom_six


! **********************************************************************
recursive subroutine intmom_rec(Npart, Jmax, i1, s1, x)
! **********************************************************************
  use ol_momenta_decl_/**/REALKIND, only: Q
  implicit none
  integer,        intent(in) :: Npart, Jmax, i1, s1, x
  integer :: A
  integer :: ix ! conventional particle numbers
  integer :: lx ! individual 2^(i-1) particles numbers
  integer :: sx ! sums of 2^(i-1) particle numbers
  integer :: rx ! inverse of sx
  logical :: last

  last = .false.
  if (2*x+2 == Npart .or. 2*x+3 == Npart) then
    last = .true.
  end if

  lx = 1 ! adding ext mom 1 <= ix < i1
  do ix = 1, i1 - 1
    sx = s1 + lx
    rx = Jmax + 1 - sx
    if ( (last .eqv. .false.) .or. (mod(Npart,2) == 1 .or. (sx < rx)) ) then  ! avoid double determination for even Npart
      do A = 1, 4
        Q(A,sx) = Q(A,s1) + Q(A,lx)
        Q(A,rx) = -Q(A,sx)
      end do
    end if
    lx = 2*lx
    if ( last .eqv. .false. ) then ! recursion
      call intmom_rec(Npart, Jmax, ix, sx, x+1)
    end if
  end do  !ix

end subroutine intmom_rec


! **********************************************************************
recursive subroutine intmom_rec_six(Npart, Jmax, i1, s1, x, use_qp_kinematics)
! **********************************************************************
  use ol_momenta_decl_/**/REALKIND, only: L,collconf,softconf
#ifdef PRECISION_dp
  use ol_parameters_decl_/**/REALKIND, only: hp_mode, hp_switch, sync_qp_kinematics
  use ol_momenta_decl_/**/QREALKIND, only: L_qp=>L
  use ol_kinematics_/**/QREALKIND, only: cont_LC_cntrv_qp=>cont_LC_cntrv
#endif
  implicit none
  integer,    intent(in) :: Npart, Jmax, i1, s1, x
  logical,    intent(in) :: use_qp_kinematics
  real(REALKIND) :: sp
  integer :: A
  integer :: ix ! conventional particle numbers
  integer :: lx ! individual 2^(i-1) particles numbers
  integer :: sx ! sums of 2^(i-1) particle numbers
  integer :: rx ! inverse of sx
  logical :: last

  last = .false.
  if (2*x+2 == Npart .or. 2*x+3 == Npart) then
    last = .true.
  end if
  lx = 1 ! adding ext mom 1 <= ix < i1
  do ix = 1, i1 - 1
    sx = s1 + lx
    rx = Jmax + 1 - sx
    if ( (last .eqv. .false.) .or. (mod(Npart,2) == 1 .or. (sx < rx)) ) then
    ! avoid double determination for even Npart
#ifdef PRECISION_dp
      ! overwrite dp invariants by qp ones (not the masses squared terms though)
      if (use_qp_kinematics .and. hp_switch .eq. 1 .and. sync_qp_kinematics .eq. 1) then
        L(1:5,sx) = L_qp(1:5,sx)
        L(1:4,rx) = L_qp(1:4,rx)
        sp = 2*cont_LC_cntrv_qp(L_qp(1:4,s1),L_qp(1:4,lx))
        L(6,sx) = sp + L(6,s1)
        L(6,3) = L_qp(6,3)
      else
        L(1:5,sx) = L(1:5,s1) + L(1:5,lx)
        L(1:4,rx) = -L(1:4,sx)
        sp = 2*cont_LC_cntrv(L(1:4,s1),L(1:4,lx))
        L(6,sx) = sp + L(6,s1)
      end if
#else
      L(1:5,sx) = L(1:5,s1) + L(1:5,lx)
      L(1:4,rx) = -L(1:4,sx)
      sp = 2*cont_LC_cntrv(L(1:4,s1),L(1:4,lx))
      L(6,sx) = sp + L(6,s1)
#endif
      L(5:6,rx) = L(5:6,sx)
#ifdef PRECISION_dp
      if(hp_mode .ge. 2) then
        if (abs(sp/(maxval(abs(L(1:4,s1)))*maxval(abs(L(1:4,lx))))) .lt. collthres) then
          collconf = ior(collconf,sx)
        end if
      end if
#endif
    end if
    lx = 2*lx
    if ( last .eqv. .false. ) then ! recursion
      call intmom_rec_six(Npart, Jmax, ix, sx, x+1, use_qp_kinematics)
    end if
  end do  !ix

end subroutine intmom_rec_six


! **********************************************************************
function squeeze_onshell(pinv, masses)
! **********************************************************************
! If 'pinv' is "close" to an element of 'masses', return the mass (positive or negative).
! Otherwise return pinv.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_parameters_decl_/**/DREALKIND, only: ti_os_thresh, mureg
  implicit none
  complex(REALKIND) :: squeeze_onshell
  complex(REALKIND), intent(inout) :: pinv
  real(REALKIND) :: masses(:), mass
  integer :: k
  squeeze_onshell = pinv
  do k = 1, size(masses)
    mass = masses(k)
    if (k /= 1 .and. mass == 0) cycle
    if (abs(abs(pinv)-mass**2)/mureg**2 < ti_os_thresh) then
      squeeze_onshell = sign(mass*mass, real(pinv))
    end if
  end do
end function squeeze_onshell


! **********************************************************************
function momenta_invariants(moms) result(invs)
! **********************************************************************
! Calculate the list of invariants from the momenta 'moms' (complex standard rep.)
! as used by Collier. Apply 'squeeze_onshell' to each invariant with the masses in the theory.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_external_decl_/**/REALKIND, only: binom2
  use ol_parameters_decl_/**/DREALKIND, only: model
  use ol_parameters_decl_/**/REALKIND, only: &
    & wMW, rMW, wMZ, rMZ, wMH, rMH, &
    & wMC, rMC, wMB, rMB, wMT, rMT, &
    & rME, wME, rMM, wMM, rML, wML, &
    & wMA0, rMA0, wMHH, rMHH, wMHp, rMHp
  implicit none
  complex(REALKIND), intent(in) :: moms(:,:)
  complex(REALKIND) :: invs(binom2(size(moms,2)+1))
  complex(REALKIND) :: moms0(0:3,0:size(moms,2))
  real(REALKIND) :: masses(0:12)
  integer :: n, k, a, b
  n = size(moms,2) + 1
  moms0(:,0) = 0
  moms0(:,1:n-1) = moms
  do k = 1, size(invs)
    invs(k) = cont_L_cmplx(moms0(:,mod(k-1,n)) - moms0(:,mod(k+((k-1)/n),n)))
  end do
  masses = 0
  n = 9
  if (wMW == 0) masses(1) = rMW
  if (wMZ == 0) masses(2) = rMZ
  if (wMH == 0) masses(3) = rMH
  if (wMC == 0) masses(4) = rMC
  if (wMB == 0) masses(5) = rMB
  if (wMT == 0) masses(6) = rMT
  if (wME == 0) masses(7) = rME
  if (wMM == 0) masses(8) = rMM
  if (wML == 0) masses(9) = rML
  if (trim(model) == "2hdm") then
    n = 12
    if (wMA0 == 0) masses(10) = rMA0
    if (wMHH == 0) masses(11) = rMHH
    if (wMHp == 0) masses(12) = rMHp
  end if
  do k = 1, size(invs)
    invs(k) = squeeze_onshell(invs(k), masses(0:n))
  end do
end function momenta_invariants

! **********************************************************************
function collier_invariants(moms) result(invs)
! **********************************************************************
! Calculate the list of invariants from the momenta-indices 'moms'
! as used by Collier.
! This function is meant to be used only to compute invariants for
! scalar integrals, namely one-loop bubbles, triangles and boxes.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_error
  use ol_external_decl_/**/REALKIND, only: binom2
  use ol_momenta_decl_/**/REALKIND, only: L
  implicit none
  integer, intent(in) :: moms(:)
  complex(REALKIND) :: invs(binom2(size(moms)+1))

  if(size(moms)==1) then
    invs(1) = L(5,moms(1)) + L(6,moms(1))
  else if(size(moms)==2) then
    invs(1) = L(5,moms(1))         + L(6,moms(1))
    invs(2) = L(5,moms(2)-moms(1)) + L(6,moms(2)-moms(1))
    invs(3) = L(5,moms(2))         + L(6,moms(2))
  else if(size(moms)==3) then
    invs(1) = L(5,moms(1))         + L(6,moms(1))
    invs(2) = L(5,moms(2)-moms(1)) + L(6,moms(2)-moms(1))
    invs(3) = L(5,moms(3)-moms(2)) + L(6,moms(3)-moms(2))
    invs(4) = L(5,moms(3))         + L(6,moms(3))
    invs(5) = L(5,moms(2))         + L(6,moms(2))
    invs(6) = L(5,moms(3)-moms(1)) + L(6,moms(3)-moms(1))
  else
    call ol_error('Collier invariants computed for a non-MI')
    invs(:) = 0
  end if

end function collier_invariants


end module ol_kinematics_/**/REALKIND
