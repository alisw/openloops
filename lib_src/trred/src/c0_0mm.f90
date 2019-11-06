!******************************************************************************!
!                                                                              !
!    c0_0mm.f90                                                                !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module c0_0mm_DP
  use triangle_aux_DP, only: target_precision,dp,cone,cnul,choose
  implicit none
  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(dp), parameter :: C_0mm_exp_thr = real(2,kind=dp)

  contains

  function C0_n_0mm(p2,m2,n) result(C0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: z,C0,loss

    z = m2(1)/p2
    if (abs(z) .gt. C_0mm_exp_thr) then
      C0 = C0_n_exp_0mm(z,p2,n,loss)
    else
      C0 = C0_n_0mm_small_z(z,p2,n)
    end if

  end function C0_n_0mm

  function C0_n_0mm_small_z(z,p2,n) result(Cn)
    ! Implements the closed formula which is stable for small z <= 2.
    ! z : m^2/p^2
    complex(dp), intent(in) :: z,p2
    integer,       intent(in) :: n
    complex(dp) :: Cn,sum
    integer       :: k

    if (n .lt. 0) then
      write (*,*) 'ERROR: called C0_n with n<0'
      stop
    end if

    sum = cnul
    if (n .gt. 0) then
      do k = 0, n-1
        sum = sum - (cone/(cone + z)**(n - k)/cmplx(k - n,kind=dp))
      end do
      sum = sum*(-1)**n
    end if

    Cn = (-(-1)**n*(log(cone + 1/z)) + sum)/(p2*(1+n))

  end function C0_n_0mm_small_z

  function C0_n_exp_0mm(z,p2,n,loss) result(Cn)
    complex(dp), intent(in) :: z,p2
    complex(dp), intent(out), optional :: loss
    integer,       intent(in) :: n
    complex(dp) :: w,Cn,Cnq
    integer       :: q

    if (n .lt. 0) then
      write (*,*) 'ERROR: called C0_n_exp_0mm with n<0'
      stop
    end if

    w = 1/z

    q = n+1
    Cn = C0_nq_coeff_0mm(w,n,q)
    Cnq = cnul
    if (present(loss)) then
      loss = Cn
    end if

    q = q+1
    Cnq = (-w)*(q-1)**2/cmplx(q,kind=dp)/cmplx(q-1-n,kind=dp)*Cn
    do while (abs(Cnq/Cn) .gt. target_precision)
      Cn = Cn + Cnq
      q = q + 1
      Cnq = (-w)*(q-1)**2/cmplx(q,kind=dp)/cmplx(q-1-n,kind=dp)*Cnq
    end do
    Cn = (Cn + Cnq)/cmplx(n+1,kind=dp)/p2

    if (present(loss)) then
     loss = log(abs(loss/Cn))/log(10._dp)
    end if

  end function C0_n_exp_0mm

  function C0_nq_coeff_0mm(w,n,q) result(C0)
    complex(dp), intent(in) :: w
    integer,       intent(in) :: n,q
    complex(dp) :: C0

    C0 = (-w)**cmplx(q,kind=dp)*cmplx(choose(q-1,n)/q,kind=dp)

  end function C0_nq_coeff_0mm

end module c0_0mm_DP

module c0_0mm_QP
  use triangle_aux_QP, only: target_precision,qp,cone,cnul,choose
  implicit none
  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(qp), parameter :: C_0mm_exp_thr = real(2,kind=qp)

  contains

  function C0_n_0mm(p2,m2,n) result(C0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: z,C0,loss

    z = m2(1)/p2
    if (abs(z) .gt. C_0mm_exp_thr) then
      C0 = C0_n_exp_0mm(z,p2,n,loss)
    else
      C0 = C0_n_0mm_small_z(z,p2,n)
    end if

  end function C0_n_0mm

  function C0_n_0mm_small_z(z,p2,n) result(Cn)
    ! Implements the closed formula which is stable for small z <= 2.
    ! z : m^2/p^2
    complex(qp), intent(in) :: z,p2
    integer,       intent(in) :: n
    complex(qp) :: Cn,sum
    integer       :: k

    if (n .lt. 0) then
      write (*,*) 'ERROR: called C0_n with n<0'
      stop
    end if

    sum = cnul
    if (n .gt. 0) then
      do k = 0, n-1
        sum = sum - (cone/(cone + z)**(n - k)/cmplx(k - n,kind=qp))
      end do
      sum = sum*(-1)**n
    end if

    Cn = (-(-1)**n*(log(cone + 1/z)) + sum)/(p2*(1+n))

  end function C0_n_0mm_small_z

  function C0_n_exp_0mm(z,p2,n,loss) result(Cn)
    complex(qp), intent(in) :: z,p2
    complex(qp), intent(out), optional :: loss
    integer,       intent(in) :: n
    complex(qp) :: w,Cn,Cnq
    integer       :: q

    if (n .lt. 0) then
      write (*,*) 'ERROR: called C0_n_exp_0mm with n<0'
      stop
    end if

    w = 1/z

    q = n+1
    Cn = C0_nq_coeff_0mm(w,n,q)
    Cnq = cnul
    if (present(loss)) then
      loss = Cn
    end if

    q = q+1
    Cnq = (-w)*(q-1)**2/cmplx(q,kind=qp)/cmplx(q-1-n,kind=qp)*Cn
    do while (abs(Cnq/Cn) .gt. target_precision)
      Cn = Cn + Cnq
      q = q + 1
      Cnq = (-w)*(q-1)**2/cmplx(q,kind=qp)/cmplx(q-1-n,kind=qp)*Cnq
    end do
    Cn = (Cn + Cnq)/cmplx(n+1,kind=qp)/p2

    if (present(loss)) then
     loss = log(abs(loss/Cn))/log(10._qp)
    end if

  end function C0_n_exp_0mm

  function C0_nq_coeff_0mm(w,n,q) result(C0)
    complex(qp), intent(in) :: w
    integer,       intent(in) :: n,q
    complex(qp) :: C0

    C0 = (-w)**cmplx(q,kind=qp)*cmplx(choose(q-1,n)/q,kind=qp)

  end function C0_nq_coeff_0mm

end module c0_0mm_QP
