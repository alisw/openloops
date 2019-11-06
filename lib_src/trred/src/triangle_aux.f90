!******************************************************************************!
!                                                                              !
!    triangle_aux.f90                                                          !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module triangle_aux_DP
  implicit none
  integer, parameter ::                            &
    sp = kind(1.0),                                &
  ! double precision
    dp = selected_real_kind(2*precision(1.0_sp)),  &
  ! quad precision
    qp = selected_real_kind(2*precision(1.0_dp))

  integer, parameter :: i8 = selected_int_kind(2*8)

  integer, parameter :: recursion_threshold = 1000

  real(dp), parameter :: pi = 3.14159265358979323846264338327950288

  real(dp), parameter :: target_precision = 10d0**(-15._dp)

  complex(dp), parameter :: cone = cmplx(1,0,kind=dp)
  complex(dp), parameter :: cnul = cmplx(0,0,kind=dp)
  complex(dp), parameter :: cima = cmplx(0,1,kind=dp)
  real(dp), parameter :: rone = real(1,kind=dp)
  real(dp), parameter :: rnul = real(0,kind=dp)

  real(dp), save :: muUV2 = 0._dp, muIR2 = 0._dp
  real(dp), save :: duv = 0._dp, dir = 0._dp

  contains

  subroutine set_muUV2(m)
    real(dp), intent(in) :: m
    muUV2 = m
  end subroutine set_muUV2

  subroutine set_muIR2(m)
    real(dp), intent(in) :: m
    muIR2 = m
  end subroutine set_muIR2

  subroutine set_duv(d)
    real(dp), intent(in) :: d
    duv = d
  end subroutine set_duv

  subroutine set_dir(d)
    real(dp), intent(in) :: d
    dir = d
  end subroutine set_dir

  function gamma_int(n) result(res)
    integer, intent(in) :: n
    real(dp) :: res
    res = factorial (n-1)
  end function gamma_int

  function acoth(x)
    real(dp), intent(in) :: x
    real(dp) :: acoth
    acoth = cone/2._dp*log( (x+cone)/(x-cone) )
  end function acoth

  function A0(m2)
    complex(dp), intent(in) :: m2
    real(dp) :: A0
    A0 = m2*(cone-log(m2/muUV2) + duv)
  end function A0

  function B0_zero()
    real(dp) :: B0_zero
    B0_zero = log(muUV2/muIR2) + duv - dir
  end function B0_zero

  function factorial (n) result (res)

    implicit none
    integer, intent (in) :: n
    real(dp) :: res
    integer :: i

    res = product ((/(real(i,kind=dp), i = 1, n)/))

  end function factorial


  function choose(n,k) result(res)
    integer, intent(in) :: n,k
    integer :: i
    real(dp)  :: res

    if (k .lt. 0_i8 .or. k .gt. n) then
      res = 0._dp
    else if (k .eq. 0_i8 .or. k .eq. n) then
      res = 1._dp
    else
      res = 1._dp
      do i = 0, min(k,n-k)-1
        res = res * real(n - i,kind=dp) / real(i + 1,kind=dp)
      end do
    end if
  end function choose

  function log1pX(x,n) result(Cnm)
    ! evaluates the log(1+x) for |x| < 0.9,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.

    complex(dp), intent(in) :: x
    integer,    intent(in) :: n
    complex(dp) :: Cnm,estim
    integer    :: steps

    if (abs(x) .gt. 0.9) then
      write(*,*)  "log1pX should not be used for x>0.9."
      stop
    end if

    ! this estim is a higher bound for x<=0.9
    estim = 300*log(target_precision)/(-16._dp)
    steps = ceiling(abs((log(target_precision*estim)-log(x))/log(x)))
    ! refine step
    steps = ceiling(abs((log(target_precision*steps)-log(x))/log(x)))+1+n
    Cnm = log1pXRec(x,n,0,steps)

  end function log1pX

  function log1pXdX(x,n) result(Cm)
    ! Evaluates the log(1+x)/x for |x| < 1,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.
    complex(dp), intent(in) :: x
    integer,    intent(in) :: n
    complex(dp) :: Cm,Cmp1,m

    m = n
    Cm = (-x)**m/cmplx(m+1,kind=dp)

    m = m+1
    Cmp1 = (-x)**m/cmplx(m+1,kind=dp)

    do while (abs(Cmp1/Cm) .gt. target_precision)
      Cm = Cm + Cmp1
      m = m + 1
      Cmp1 = (-x)**m/cmplx(m+1,kind=dp)
    end do
    Cm = Cm + Cmp1

  end function log1pXdX

  recursive function log1pXRec(x,n,m1,m2) result(Cnm)
    ! evaluates the log(1+x) for |x| < 0,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.

    complex(dp), intent(in) :: x
    integer,    intent(in) :: n,m1,m2
    complex(dp) :: Cnm

    ! check
    if (m1 .lt. m2)  then
      if (n .gt. m1) then
        Cnm = x*(-log1pXRec(x,n,m1+1,m2))
      else
        Cnm = x*(cone/cmplx(m1+1,kind=dp)-log1pXRec(x,n,m1+1,m2))
      end if
    else
      Cnm = x/cmplx(m1,kind=dp)
    end if

  end function log1pXRec

  function HarmNum(n) result(Hn)
    ! Computes the harmonic number Hn = Sum_k=1^n 1/k
    integer, intent(in) :: n
    complex(dp) :: Hn
    integer    :: i
    Hn = 0
    do i = 1, n
      Hn = Hn + cone/i
    end do
  end function HarmNum

  function Lphi(v,n,offset) result(lp)
    ! Computes LerchPhi(v,n) = Sum_k=0^\infty v^k/(k+n+2)
    complex(dp), intent(in)           :: v
    integer,    intent(in)           :: n
    integer,    intent(in), optional :: offset
    complex(dp) :: lp,lpn
    integer    :: k

    if (present(offset)) then
      k = offset
      lpn = v**k/cmplx(k+n,kind=dp)
    else
      k = 0
      lp = cone/cmplx(n,kind=dp)
    end if

    k = k+1
    lpn = v**k/cmplx(k+n,kind=dp)
    do while (abs(lpn/lp) .gt. target_precision)
      lp = lp + lpn
      k = k + 1
      lpn = v**k/cmplx(k+n,kind=dp)
    end do
    lp = lp + lpn
  end function Lphi

  function LphiLog(v,n) result(lp)
    ! Computes LerchPhi(v,n) + Log(1-v) = Sum_k=1^\infty v^k/k/(k+n-1)
    complex(dp), intent(in)           :: v
    integer,    intent(in)           :: n
    complex(dp) :: lp,lpn
    integer    :: k

    if (n .lt. 1) then
      write (*,*) 'ERROR: called LphiLog with n<1'
      stop
    end if

    k = 0
    lp = v/cmplx(n,kind=dp)

    k = k+1
    lpn = v**k/(cmplx(k,kind=dp)*cmplx(k+n-1,kind=dp))
    do while (abs(lpn/lp) .gt. target_precision)
      lp = lp + lpn
      k = k + 1
      lpn = v**k/(cmplx(k,kind=dp)*cmplx(k+n-1,kind=dp))
      write(*,*) "abs(lpn/lp):", abs(lpn/lp)
    end do
    lp = lp + lpn
  end function LphiLog

  function zlogzf(w) result(res)
    ! Computes z Log(z) - (1 + z) Log(1 + z)
    complex(dp), intent(in) :: w
    complex(dp) :: z,res,fac,sum,sumn,n,ffac

    if (abs(w) .gt. 1/100._dp .and. abs(w) .lt. 100._dp) then
    ! stable region, no expansion performed
      res = w*log(w) - (1+w)*log(1+w)
    else
      ! formula fulfills: zlogzf(1/w) = w*zlogzf(w)
      if (abs(w) .ge. 2._dp) then
        z = 1/w
        ffac = w
      else
        z = w
        ffac = cone
      end if

      n = 1
      fac = cone
      sum = -cone/fac

      n = n + 1
      fac = fac + cmplx(n,kind=dp)
      sumn = z/fac

      do while (abs(sumn/sum) .gt. target_precision)
        sum = sum + sumn
        n = n + 1
        sumn = fac*sumn
        fac = fac + cmplx(n,kind=dp)
        sumn = -1*sumn*z/fac
      end do
      sum = sum + sumn

      res = z*(log(z)-1+z*sum/2._dp)*ffac
    end if
  end function zlogzf

  function SvM1(z,n) result(Sv)
    ! Computes (1/(1+z))^n - 1 for small z
    complex(dp), intent(in) :: z
    integer,    intent(in) :: n
    integer    :: q
    complex(dp) :: Sv,Svnq

    if (abs(z) .gt. 1) then
      Sv = cone/(1+z)**n - cone
    else if (abs(z*n) .gt. 4) then
      Sv = exp(-n*log(1+z))-1
    else
      q = 1
      if (n .eq. 0) then
        Sv = -z
      else
        Sv = -z*cmplx(q+n-1,kind=dp)/cmplx(q,kind=dp)
      end if

      q = q + 1
      if (n .eq. 0) then
        Svnq = -z*Sv
      else
        Svnq = (-z)*cmplx(q+n-1,kind=dp)/cmplx(q,kind=dp)*Sv
      end if

      do while (abs(Svnq/Sv) .gt. target_precision)
        Sv = Sv + Svnq
        q = q + 1
        if (n .eq. 0) then
          Svnq = -z*Svnq
        else
          Svnq = (-z)*cmplx(q+n-1,kind=dp)/cmplx(q,kind=dp)*Svnq
        end if
      end do
      Sv = Sv + Svnq
    end if

  end function SvM1

  function Sv1(z,n) result(Sv)
    ! Computes ((1+z)^n - 1
    complex(dp), intent(in) :: z
    integer,    intent(in) :: n
    integer    :: q
    complex(dp) :: Sv,Svnq

    if (n .eq. 0) then
      Sv = z
    else
      Sv = (cone+z)**n - cone
      !Sv = exp(n*log(cone+z))-cone

      if (abs(Sv) .lt. 0.1) then
        q = 1
        Sv = z*cmplx(n-q+1,kind=dp)/cmplx(q,kind=dp)

        q = q + 1
        Svnq = z*cmplx(n-q+1,kind=dp)/cmplx(q,kind=dp)*Sv
        do while (abs(Svnq/Sv) .gt. target_precision)
          Sv = Sv + Svnq
          q = q + 1
          Svnq = z*cmplx(n-q+1,kind=dp)/cmplx(q,kind=dp)*Svnq
        end do
        Sv = Sv + Svnq
      end if
    end if

  end function Sv1

  function A0mB0(p2,m2) result(AB)
    ! Computes A0(m2)-m2 B0(p2,0,m2) + (p2+m2) B0^1(p2,0,m2) for small p2/m2

    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2
    complex(dp) :: w,AB,ABk
    integer :: k

    w = p2/m2
    if (abs(w) .gt. 0.5_dp) then
      write (*,*) 'ERROR: called A0mB0 with p2/m2>=0.5'
      stop
    end if

    k = 1
    AB = (-w)**k*2/cmplx((k+1)*(k+2),kind=dp)

    k = k + 1
    ABk = (-w)**k*2/cmplx((k+1)*(k+2),kind=dp)

    do while (abs(ABk/AB) .gt. target_precision)
      AB = AB + ABk
      k = k + 1
      ABk = (-w)**k*2/cmplx((k+1)*(k+2),kind=dp)
    end do
    AB = (AB + ABk)*p2

  end function A0mB0

  function A0mB0_p1p1p1(p2,m2,d) result(AB)
    ! Computes: A0 - m2*(B0^{0,(0)} + B0^{0,(1)} + B0^{0,(2)}) + (a*(B0^{1,(0)} + B0^{1,(1)}) + r)/b
    ! with:
    ! a = ((2*m2^2-5*m2*p2+11*p2^2)*(1+d)^2)
    ! r = -m2 p2 (1 + d)
    ! b = (-5*p2*(1+d)+2*m2*(2+d))

    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/((1+w)*(5*w*(1 + d)-2*(2+d)))
    c1 = w**3*(30+78*d+38*d**2+2*w**2*(1+d)*(-6+11*d**2)+w*(18+d*(12-d*(21+11*d))))/6
    c2 = (1+w)*(-6+w*(6+(3-5*d)*d+w*(1+d)*(-6+11*d**2)))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_p1p1p1

  function A0mB0_p1p1P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(3+d-5*(1+d)*w)
    c1 = -w**3*(-d*(11+9*d)+12*(-1+w)+d*(15+d*(5+6*d))*w-2*(1+d)*(6+d*(17+6*d))*w**2)/(6*(1+w))
    c2 = (3-3*d+(-3+d*(4+9*d))*w+(1+d)*(6+d*(17+6*d))*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)
  end function A0mB0_p1p1P12

  function A0mB0_p1P12P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(-2+5*(1+d)*w)
    c1 = w**3*(-16-6*d*(4+3*d)+30*w+d*(46+3*d*(5+d))*w+2*(14+d*(20-3*(-1+d)*d))*w**2)/(6*(1+w))
    c2 = -(-6*d-(22+33*d+9*d**2)*w+(1+d)*(-14+3*(-2+d)*d)*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_p1P12P12

  function A0mB0_P12P12P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(-1+d+5*(1+d)*w)
    c1 = (w**3*(-4-d*(5+11*d)+28*w+d*(39+d*(9+2*d))*w+2*(10+d*(13+d-2*d**2))*w**2))/(6*(1+w))
    c2 = -(-3-9*d-(19+d*(26+5*d))*w+(1+d)*(-10+d*(-3+2*d))*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_P12P12P12

  function A0mB0_gP12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = (w**2*(-1+d**2*(1-2*w)*w+w*(2+w-2*w**2)-d*(1+w)*(1+w+2*w**2)))/(6*(1+w))
    c2 = -(-1+2*d+d*(2+d)*w+(1+d)*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_gP12

  function A0mB0_gp1(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = (w**2*(4+4*d+10*w-(-4+d)*d*w+2*w**2+4*(-1+d**2)*w**3))/(6*(1 + w))
    c2 = (2-d-2*w**2+d**2*w*(1+2*w))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_gp1

  function A0mB0_0mm_p1p1(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(-8-9*d+w*(8+6*d))/6
    c2 = (2+w*(4+3*d))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1

  function A0mB0_0mm_p1P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(7+6*d-2*w*(1+2*d))/6
    c2 = -(-2+w+2*w*d)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1P12

  function A0mB0_0mm_P12P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(4+(3-2*w)*d)/6
    c2 = (2-w*d)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_P12P12

  function A0mB0_0mm_p1p1p1(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(7*w*(1+d)+2*(2+d))
    c1 = w**3*(-18*(1+w)*(3+w*(-1+3*w))-9*(12+w*(3+w+6*w**2))*d + &
                (-52+w*(15+w*(10+33*w)))*d**2+11*w*(2+w*(-1+3*w))*d**3)/(12*(1+w))
    c2 = (6+w*(18+(9-7*d)*d-w*(1+d)*(-18+11*d**2)))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1p1

  function A0mB0_0mm_p1p1P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(3+d+7*w*(1+d))
    c1 = w**3*(63+127*d+66*d**2-(36+d*(69+d*(61+12*d)))*w+(-9+d*(59+68*d+6*d**2))*w**2 &
               -3*(1+d)*(18+d*(17+6*d))*w**3)/(12*(1+w))

    c2 = (3-3*d+w*(9-d*(8+15*d)+w*(1+d)*(18+d*(17+6*d))))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1P12

  function A0mB0_0mm_p1P12P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(2+7*w*(1+d))
    c1 = w**3*(-16-6*d*(5+2*d)+42*w+d*(85+3*(7-2*d)*d)*w &
          + (-116+d*(-161+3*(-16+d)*d))*w**2+3*(2+d*(8-3*(-1+d)*d))*w**3)/(12*(1+w))
    c2 = (-6*d+w*(38+3*d*(17+5*d)+(1+d)*(-2+3*(-2+d)*d)*w))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1P12P12

  function A0mB0_0mm_P12P12P12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(1-d+7*(1+d)*w)

    c1 = w**3*(-7-d*(11+2*d)+34*w+d*(63+(9-4*d)*d)*w                       &
               +(-85+d*(-103+2*(-11+d)*d))*w**2+3*(-2+d+d**2-2*d**3)*w**3) &
               /(12*(1+w))

    c2 = (-3-9*d+w*(29+2*w+d*(34+7*d+(-1+d)*(1+2*d)*w)))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_P12P12P12

  function A0mB0_0mm_gp1(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2

    c1 = -w**2*(28+2*w*(-9+2*w)+d*(29+4*(-4+w)*w))/6
    c2 = (2+w*(6+5*d-2*(1+d)*w))*log1pXdX(w,3)

    AB = prefac*(c1+c2)

  end function A0mB0_0mm_gp1

  function A0mB0_0mm_gP12(p2,m2,d) result(AB)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2

    c1 = w**2*(12+11*d-2*(4+5*d)*w)/6
    c2 =  (2-(4+5*d)*w)*log1pXdX(w,3)

    AB = prefac*(c1+c2)

  end function A0mB0_0mm_gP12

end module triangle_aux_DP

module triangle_aux_QP
  implicit none
  integer, parameter ::                            &
    sp = kind(1.0),                                &
  ! double precision
    dp = selected_real_kind(2*precision(1.0_sp)),  &
  ! quad precision
    qp = selected_real_kind(2*precision(1.0_dp))

  integer, parameter :: i8 = selected_int_kind(2*8)

  integer, parameter :: recursion_threshold = 1000

  real(qp), parameter :: target_precision = 10d0**(-34._qp)

  real(qp), parameter :: pi = 3.14159265358979323846264338327950288

  complex(qp), parameter :: cone = cmplx(1,0,kind=qp)
  complex(qp), parameter :: cnul = cmplx(0,0,kind=qp)
  complex(qp), parameter :: cima = cmplx(0,1,kind=qp)
  real(qp), parameter :: rone = real(1,kind=qp)
  real(qp), parameter :: rnul = real(0,kind=qp)

  real(qp), save :: muUV2 = 0._qp, muIR2 = 0._qp
  real(qp), save :: duv = 0._qp, dir = 0._qp

  contains

  subroutine set_muUV2(m)
    real(qp), intent(in) :: m
    muUV2 = m
  end subroutine set_muUV2

  subroutine set_muIR2(m)
    real(qp), intent(in) :: m
    muIR2 = m
  end subroutine set_muIR2

  subroutine set_duv(d)
    real(qp), intent(in) :: d
    duv = d
  end subroutine set_duv

  subroutine set_dir(d)
    real(qp), intent(in) :: d
    dir = d
  end subroutine set_dir

  function gamma_int(n) result(res)
    integer, intent(in) :: n
    real(qp) :: res
    res = factorial (n-1)
  end function gamma_int

  function acoth(x)
    real(qp), intent(in) :: x
    real(qp) :: acoth
    acoth = cone/2._qp*log( (x+cone)/(x-cone) )
  end function acoth

  function A0(m2)
    complex(qp), intent(in) :: m2
    real(qp) :: A0
    A0 = m2*(cone-log(m2/muUV2) + duv)
  end function A0

  function B0_zero()
    real(qp) :: B0_zero
    B0_zero = log(muUV2/muIR2) + duv - dir
  end function B0_zero

  function factorial (n) result (res)

    implicit none
    integer, intent (in) :: n
    real(qp) :: res
    integer :: i

    res = product ((/(real(i,kind=qp), i = 1, n)/))

  end function factorial

  function choose(n,k) result(res)
    integer, intent(in) :: n,k
    integer :: i
    real(qp)  :: res

    if (k .lt. 0_i8 .or. k .gt. n) then
      res = 0._qp
    else if (k .eq. 0_i8 .or. k .eq. n) then
      res = 1._qp
    else
      res = 1._qp
      do i = 0, min(k,n-k)-1
        res = res * real(n - i,kind=qp) / real(i + 1,kind=qp)
      end do
    end if
  end function choose

  function log1pX(x,n) result(Cnm)
    ! evaluates the log(1+x) for |x| < 0.9,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.

    complex(qp), intent(in) :: x
    integer,    intent(in) :: n
    complex(qp) :: Cnm,estim
    integer    :: steps

    if (abs(x) .gt. 0.9) then
      write(*,*)  "log1pX should not be used for x>0.9."
      stop
    end if

    ! this estim is a higher bound for x<=0.9
    estim = 300*log(target_precision)/(-16._qp)
    steps = ceiling(abs((log(target_precision*estim)-log(x))/log(x)))
    ! refine step
    steps = ceiling(abs((log(target_precision*steps)-log(x))/log(x)))+1+n
    Cnm = log1pXRec(x,n,0,steps)

  end function log1pX

  function log1pXdX(x,n) result(Cm)
    ! Evaluates the log(1+x)/x for |x| < 1,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.
    complex(qp), intent(in) :: x
    integer,    intent(in) :: n
    complex(qp) :: Cm,Cmp1,m

    m = n
    Cm = (-x)**m/cmplx(m+1,kind=qp)

    m = m+1
    Cmp1 = (-x)**m/cmplx(m+1,kind=qp)

    do while (abs(Cmp1/Cm) .gt. target_precision)
      Cm = Cm + Cmp1
      m = m + 1
      Cmp1 = (-x)**m/cmplx(m+1,kind=qp)
    end do
    Cm = Cm + Cmp1

  end function log1pXdX

  recursive function log1pXRec(x,n,m1,m2) result(Cnm)
    ! evaluates the log(1+x) for |x| < 0,
    ! For n>0 the first n orders in x^1, x^2, .. x^(n-1) are ommited.

    complex(qp), intent(in) :: x
    integer,    intent(in) :: n,m1,m2
    complex(qp) :: Cnm

    ! check
    if (m1 .lt. m2)  then
      if (n .gt. m1) then
        Cnm = x*(-log1pXRec(x,n,m1+1,m2))
      else
        Cnm = x*(cone/cmplx(m1+1,kind=qp)-log1pXRec(x,n,m1+1,m2))
      end if
    else
      Cnm = x/cmplx(m1,kind=qp)
    end if

  end function log1pXRec

  function HarmNum(n) result(Hn)
    ! Computes the harmonic number Hn = Sum_k=1^n 1/k
    integer, intent(in) :: n
    complex(qp) :: Hn
    integer    :: i
    Hn = 0
    do i = 1, n
      Hn = Hn + cone/i
    end do
  end function HarmNum

  function Lphi(v,n,offset) result(lp)
    ! Computes LerchPhi(v,n) = Sum_k=0^\infty v^k/(k+n+2)
    complex(qp), intent(in)           :: v
    integer,    intent(in)           :: n
    integer,    intent(in), optional :: offset
    complex(qp) :: lp,lpn
    integer    :: k

    if (present(offset)) then
      k = offset
      lpn = v**k/cmplx(k+n,kind=qp)
    else
      k = 0
      lp = cone/cmplx(n,kind=qp)
    end if

    k = k+1
    lpn = v**k/cmplx(k+n,kind=qp)
    do while (abs(lpn/lp) .gt. target_precision)
      lp = lp + lpn
      k = k + 1
      lpn = v**k/cmplx(k+n,kind=qp)
    end do
    lp = lp + lpn
  end function Lphi

  function LphiLog(v,n) result(lp)
    ! Computes LerchPhi(v,n) + Log(1-v) = Sum_k=1^\infty v^k/k/(k+n-1)
    complex(qp), intent(in)           :: v
    integer,    intent(in)           :: n
    complex(qp) :: lp,lpn
    integer    :: k

    if (n .lt. 1) then
      write (*,*) 'ERROR: called LphiLog with n<1'
      stop
    end if

    k = 0
    lp = v/cmplx(n,kind=qp)

    k = k+1
    lpn = v**k/(cmplx(k,kind=qp)*cmplx(k+n-1,kind=qp))
    do while (abs(lpn/lp) .gt. target_precision)
      lp = lp + lpn
      k = k + 1
      lpn = v**k/(cmplx(k,kind=qp)*cmplx(k+n-1,kind=qp))
      write(*,*) "abs(lpn/lp):", abs(lpn/lp)
    end do
    lp = lp + lpn
  end function LphiLog

  function zlogzf(w) result(res)
    ! Computes z Log(z) - (1 + z) Log(1 + z)
    complex(qp), intent(in) :: w
    complex(qp) :: z,res,fac,sum,sumn,n,ffac

    if (abs(w) .gt. 1/100._qp .and. abs(w) .lt. 100._qp) then
    ! stable region, no expansion performed
      res = w*log(w) - (1+w)*log(1+w)
    else
      ! formula fulfills: zlogzf(1/w) = w*zlogzf(w)
      if (abs(w) .ge. 2._qp) then
        z = 1/w
        ffac = w
      else
        z = w
        ffac = cone
      end if

      n = 1
      fac = cone
      sum = -cone/fac

      n = n + 1
      fac = fac + cmplx(n,kind=qp)
      sumn = z/fac

      do while (abs(sumn/sum) .gt. target_precision)
        sum = sum + sumn
        n = n + 1
        sumn = fac*sumn
        fac = fac + cmplx(n,kind=qp)
        sumn = -1*sumn*z/fac
      end do
      sum = sum + sumn

      res = z*(log(z)-1+z*sum/2._qp)*ffac
    end if
  end function zlogzf

  function SvM1(z,n) result(Sv)
    ! Computes (1/(1+z))^n - 1 for small z
    complex(qp), intent(in) :: z
    integer,    intent(in) :: n
    integer    :: q
    complex(qp) :: Sv,Svnq

    if (abs(z) .gt. 1) then
      Sv = cone/(1+z)**n - cone
    else if (abs(z*n) .gt. 4) then
      Sv = exp(-n*log(1+z))-1
    else
      q = 1
      if (n .eq. 0) then
        Sv = -z
      else
        Sv = -z*cmplx(q+n-1,kind=qp)/cmplx(q,kind=qp)
      end if

      q = q + 1
      if (n .eq. 0) then
        Svnq = -z*Sv
      else
        Svnq = (-z)*cmplx(q+n-1,kind=qp)/cmplx(q,kind=qp)*Sv
      end if

      do while (abs(Svnq/Sv) .gt. target_precision)
        Sv = Sv + Svnq
        q = q + 1
        if (n .eq. 0) then
          Svnq = -z*Svnq
        else
          Svnq = (-z)*cmplx(q+n-1,kind=qp)/cmplx(q,kind=qp)*Svnq
        end if
      end do
      Sv = Sv + Svnq
    end if

  end function SvM1

  function Sv1(z,n) result(Sv)
    ! Computes ((1+z)^n - 1
    complex(qp), intent(in) :: z
    integer,    intent(in) :: n
    integer    :: q
    complex(qp) :: Sv,Svnq

    if (n .eq. 0) then
      Sv = z
    else
      Sv = (cone+z)**n - cone
      !Sv = exp(n*log(cone+z))-cone

      if (abs(Sv) .lt. 0.1) then
        q = 1
        Sv = z*cmplx(n-q+1,kind=qp)/cmplx(q,kind=qp)

        q = q + 1
        Svnq = z*cmplx(n-q+1,kind=qp)/cmplx(q,kind=qp)*Sv
        do while (abs(Svnq/Sv) .gt. target_precision)
          Sv = Sv + Svnq
          q = q + 1
          Svnq = z*cmplx(n-q+1,kind=qp)/cmplx(q,kind=qp)*Svnq
        end do
        Sv = Sv + Svnq
      end if
    end if

  end function Sv1

  function A0mB0(p2,m2) result(AB)
    ! Computes A0(m2)-m2 B0(p2,0,m2) + (p2+m2) B0^1(p2,0,m2) for small p2/m2

    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2
    complex(qp) :: w,AB,ABk
    integer :: k

    w = p2/m2
    if (abs(w) .gt. 0.5_qp) then
      write (*,*) 'ERROR: called A0mB0 with p2/m2>=0.5'
      stop
    end if

    k = 1
    AB = (-w)**k*2/cmplx((k+1)*(k+2),kind=qp)

    k = k + 1
    ABk = (-w)**k*2/cmplx((k+1)*(k+2),kind=qp)

    do while (abs(ABk/AB) .gt. target_precision)
      AB = AB + ABk
      k = k + 1
      ABk = (-w)**k*2/cmplx((k+1)*(k+2),kind=qp)
    end do
    AB = (AB + ABk)*p2

  end function A0mB0

  function A0mB0_p1p1p1(p2,m2,d) result(AB)
    ! Computes: A0 - m2*(B0^{0,(0)} + B0^{0,(1)} + B0^{0,(2)}) + (a*(B0^{1,(0)} + B0^{1,(1)}) + r)/b
    ! with:
    ! a = ((2*m2^2-5*m2*p2+11*p2^2)*(1+d)^2)
    ! r = -m2 p2 (1 + d)
    ! b = (-5*p2*(1+d)+2*m2*(2+d))

    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/((1+w)*(5*w*(1 + d)-2*(2+d)))
    c1 = w**3*(30+78*d+38*d**2+2*w**2*(1+d)*(-6+11*d**2)+w*(18+d*(12-d*(21+11*d))))/6
    c2 = (1+w)*(-6+w*(6+(3-5*d)*d+w*(1+d)*(-6+11*d**2)))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_p1p1p1

  function A0mB0_p1p1P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(3+d-5*(1+d)*w)
    c1 = -w**3*(-d*(11+9*d)+12*(-1+w)+d*(15+d*(5+6*d))*w-2*(1+d)*(6+d*(17+6*d))*w**2)/(6*(1+w))
    c2 = (3-3*d+(-3+d*(4+9*d))*w+(1+d)*(6+d*(17+6*d))*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)
  end function A0mB0_p1p1P12

  function A0mB0_p1P12P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(-2+5*(1+d)*w)
    c1 = w**3*(-16-6*d*(4+3*d)+30*w+d*(46+3*d*(5+d))*w+2*(14+d*(20-3*(-1+d)*d))*w**2)/(6*(1+w))
    c2 = -(-6*d-(22+33*d+9*d**2)*w+(1+d)*(-14+3*(-2+d)*d)*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_p1P12P12

  function A0mB0_P12P12P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(-1+d+5*(1+d)*w)
    c1 = (w**3*(-4-d*(5+11*d)+28*w+d*(39+d*(9+2*d))*w+2*(10+d*(13+d-2*d**2))*w**2))/(6*(1+w))
    c2 = -(-3-9*d-(19+d*(26+5*d))*w+(1+d)*(-10+d*(-3+2*d))*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_P12P12P12

  function A0mB0_gP12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = (w**2*(-1+d**2*(1-2*w)*w+w*(2+w-2*w**2)-d*(1+w)*(1+w+2*w**2)))/(6*(1+w))
    c2 = -(-1+2*d+d*(2+d)*w+(1+d)*w**2)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_gP12

  function A0mB0_gp1(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = (w**2*(4+4*d+10*w-(-4+d)*d*w+2*w**2+4*(-1+d**2)*w**3))/(6*(1 + w))
    c2 = (2-d-2*w**2+d**2*w*(1+2*w))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_gp1

  function A0mB0_0mm_p1p1(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(-8-9*d+w*(8+6*d))/6
    c2 = (2+w*(4+3*d))*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1

  function A0mB0_0mm_p1P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(7+6*d-2*w*(1+2*d))/6
    c2 = -(-2+w+2*w*d)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1P12

  function A0mB0_0mm_P12P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2
    c1 = w**2*(4+(3-2*w)*d)/6
    c2 = (2-w*d)*log1pXdX(w,3)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_P12P12

  function A0mB0_0mm_p1p1p1(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(7*w*(1+d)+2*(2+d))
    c1 = w**3*(-18*(1+w)*(3+w*(-1+3*w))-9*(12+w*(3+w+6*w**2))*d + &
                (-52+w*(15+w*(10+33*w)))*d**2+11*w*(2+w*(-1+3*w))*d**3)/(12*(1+w))
    c2 = (6+w*(18+(9-7*d)*d-w*(1+d)*(-18+11*d**2)))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1p1

  function A0mB0_0mm_p1p1P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(3+d+7*w*(1+d))
    c1 = w**3*(63+127*d+66*d**2-(36+d*(69+d*(61+12*d)))*w+(-9+d*(59+68*d+6*d**2))*w**2 &
               -3*(1+d)*(18+d*(17+6*d))*w**3)/(12*(1+w))

    c2 = (3-3*d+w*(9-d*(8+15*d)+w*(1+d)*(18+d*(17+6*d))))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1p1P12

  function A0mB0_0mm_p1P12P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(2+7*w*(1+d))
    c1 = w**3*(-16-6*d*(5+2*d)+42*w+d*(85+3*(7-2*d)*d)*w &
          + (-116+d*(-161+3*(-16+d)*d))*w**2+3*(2+d*(8-3*(-1+d)*d))*w**3)/(12*(1+w))
    c2 = (-6*d+w*(38+3*d*(17+5*d)+(1+d)*(-2+3*(-2+d)*d)*w))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_p1P12P12

  function A0mB0_0mm_P12P12P12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2/(1-d+7*(1+d)*w)

    c1 = w**3*(-7-d*(11+2*d)+34*w+d*(63+(9-4*d)*d)*w                       &
               +(-85+d*(-103+2*(-11+d)*d))*w**2+3*(-2+d+d**2-2*d**3)*w**3) &
               /(12*(1+w))

    c2 = (-3-9*d+w*(29+2*w+d*(34+7*d+(-1+d)*(1+2*d)*w)))*log1pXdX(w,4)
    AB = prefac*(c1+c2)

  end function A0mB0_0mm_P12P12P12

  function A0mB0_0mm_gp1(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2

    c1 = -w**2*(28+2*w*(-9+2*w)+d*(29+4*(-4+w)*w))/6
    c2 = (2+w*(6+5*d-2*(1+d)*w))*log1pXdX(w,3)

    AB = prefac*(c1+c2)

  end function A0mB0_0mm_gp1

  function A0mB0_0mm_gP12(p2,m2,d) result(AB)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp) :: w,AB,c1,c2,prefac

    w = p2/m2
    prefac = m2

    c1 = w**2*(12+11*d-2*(4+5*d)*w)/6
    c2 =  (2-(4+5*d)*w)*log1pXdX(w,3)

    AB = prefac*(c1+c2)

  end function A0mB0_0mm_gP12

end module triangle_aux_QP
