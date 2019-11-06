!******************************************************************************!
!                                                                              !
!    b0.f90                                                                    !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module b0_DP
  use triangle_aux_DP, only: target_precision,dp,cone,cnul,choose,Lphi,zlogzf,duv,muUV2
  implicit none

  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(dp), parameter :: B_exp_thr = real(3,kind=dp)

  contains

  function B0_0_n(p2,m2,n) result(B0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,    intent(in)    :: n
    complex(dp)             :: B0

    if (n==0) then
      B0 = 2._dp + log(muUV2/p2) + duv
    else
      B0 = (-1)**n/ cmplx(n,kind=dp)
    end if

  end function B0_0_n

  function B0_n(p2,m2,n) result(B0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: z,B0

    z = m2(1)/p2
    if (n .eq. 0) then
      B0 = B0_0(z,p2)
    else
      if (abs(z) .gt. B_exp_thr) then
        B0 = B0_n_large_z(z,n)
      else
        B0 = B0_n_small_z(z,n)
      end if
    end if

  end function B0_n

  function B0_0(z,p2) result(B0)
    complex(dp), intent(in) :: z,p2
    complex(dp) :: B0

    ! log terms become numerically instable for large/small z
    !B0 = 2._dp + log(muUV2/p2) + z*log(z/(z + cone))-log(z + cone)
    B0 = 2._dp + log(muUV2/p2) + zlogzf(z) + duv

  end function B0_0

  function B0_n_small_z(z,n) result(B0)
    complex(dp), intent(in) :: z
    integer      , intent(in) :: n
    complex(dp) :: B0,sum
    integer       :: k

    if (n .lt. 1) then
      write (*,*) 'ERROR: called B0_n with n<1'
      stop
    end if

    sum = cnul
    do k = 1, n
      sum = sum + (cone/(z + cone))**k/cmplx(k,kind=dp)
    end do

    B0 = (-cone)**n*(z*log(z/(z + cone)) + &
                   (cone/(z + cone))**n/cmplx(n,kind=dp) +  z*sum)

  end function B0_n_small_z

  function B0_n_large_z(z,n) result(B0)
    complex(dp), intent(in) :: z
    integer      , intent(in) :: n
    complex(dp) :: v,B0

    v = 1/(cone+z)

    if (n .lt. 0) then
      write (*,*) 'ERROR: called B0_n_exp with n<0'
      stop
    else if (n .eq. 0) then
      B0 = v*(1+ z*(-Lphi(v,1,1))/cmplx(n,kind=dp))
    else
      B0 = (-1)**n*v**(n+1)*(1+ z*(cone-n*Lphi(v,n+1)))/cmplx(n,kind=dp)
    end if

  end function B0_n_large_z

  function B0_0_n_EP1(p2,m2,n) result(B0)
    complex(dp), intent(in) :: p2,m2
    integer,    intent(in)    :: n
    complex(dp)             :: B0

    if (n==0) then
      B0 = 1._dp
    else
      B0 = 0._dp
    end if

  end function B0_0_n_EP1

end module b0_DP

module b0_QP
  use triangle_aux_QP, only: target_precision,qp,cone,cnul,choose,Lphi,zlogzf,duv,muUV2
  implicit none

  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(qp), parameter :: B_exp_thr = real(3,kind=qp)

  contains

  function B0_0_n(p2,m2,n) result(B0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,    intent(in)    :: n
    complex(qp)             :: B0

    if (n==0) then
      B0 = 2._qp + log(muUV2/p2) + duv
    else
      B0 = (-1)**n/ cmplx(n,kind=qp)
    end if

  end function B0_0_n

  function B0_n(p2,m2,n) result(B0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: z,B0

    z = m2(1)/p2
    if (n .eq. 0) then
      B0 = B0_0(z,p2)
    else
      if (abs(z) .gt. B_exp_thr) then
        B0 = B0_n_large_z(z,n)
      else
        B0 = B0_n_small_z(z,n)
      end if
    end if

  end function B0_n

  function B0_0(z,p2) result(B0)
    complex(qp), intent(in) :: z,p2
    complex(qp) :: B0

    ! log terms become numerically instable for large/small z
    !B0 = 2._qp + log(muUV2/p2) + z*log(z/(z + cone))-log(z + cone)
    B0 = 2._qp + log(muUV2/p2) + zlogzf(z) + duv

  end function B0_0

  function B0_n_small_z(z,n) result(B0)
    complex(qp), intent(in) :: z
    integer      , intent(in) :: n
    complex(qp) :: B0,sum
    integer       :: k

    if (n .lt. 1) then
      write (*,*) 'ERROR: called B0_n with n<1'
      stop
    end if

    sum = cnul
    do k = 1, n
      sum = sum + (cone/(z + cone))**k/cmplx(k,kind=qp)
    end do

    B0 = (-cone)**n*(z*log(z/(z + cone)) + &
                   (cone/(z + cone))**n/cmplx(n,kind=qp) +  z*sum)

  end function B0_n_small_z

  function B0_n_large_z(z,n) result(B0)
    complex(qp), intent(in) :: z
    integer      , intent(in) :: n
    complex(qp) :: v,B0

    v = 1/(cone+z)

    if (n .lt. 0) then
      write (*,*) 'ERROR: called B0_n_exp with n<0'
      stop
    else if (n .eq. 0) then
      B0 = v*(1+ z*(-Lphi(v,1,1))/cmplx(n,kind=qp))
    else
      B0 = (-1)**n*v**(n+1)*(1+ z*(cone-n*Lphi(v,n+1)))/cmplx(n,kind=qp)
    end if

  end function B0_n_large_z

  function B0_0_n_EP1(p2,m2,n) result(B0)
    complex(qp), intent(in) :: p2,m2
    integer,    intent(in)    :: n
    complex(qp)             :: B0

    if (n==0) then
      B0 = 1._qp
    else
      B0 = 0._qp
    end if

  end function B0_0_n_EP1

end module b0_QP
