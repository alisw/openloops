!******************************************************************************!
!                                                                              !
!    c0_m0m1m1.f90                                                             !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module c0_m0m1m1_DP
  use triangle_aux_DP, only: target_precision,dp,cone,cnul,choose
  implicit none
  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(dp), parameter :: C_0mm_exp_thr = real(2,kind=dp)

  contains

  function Crec(yk,yr,n)
    complex(dp), intent(in) :: yk,yr
    integer,       intent(in) :: n
    complex(dp)             :: Crec(n+1)
    integer                   :: i,l

    Crec(n+1) = yk**n * (cone-yk)**(n+1) / (yk-yr)**(n+1)
    do l = n, 1, -1
      Crec(l) = cnul
      do i = 1, n+1-l
        Crec(l) = Crec(l) + Crec(l+i) * lambda(yk,yr,i,n)
      end do
      Crec(l) = Crec(l) / (1+n-l)

    end do

    contains

    function lambda(yk,yr,i,n) result(lam)
      complex(dp), intent(in) :: yk,yr
      integer,       intent(in) :: i,n
      complex(dp)             :: lam

      lam = (1+n)*(cone/(yr-yk)**i - cone/(cone-yk)**i)-n/(-yk)**i

    end function lambda

  end function Crec


  function C0_n_m0m1m1(p2,m2,n) result(C0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: z0,z1,yr1,yr2,C0,loss,Cc1(n+1),Cc2(n+1)
    integer                   :: j

    z0 = m2(1)/p2
    z1 = m2(2)/p2

    yr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    yr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    Cc1 = Crec(yr1,yr2,n)
    Cc2 = Crec(yr2,yr1,n)

    C0 = log((yr1-cone)/yr1) * Cc1(1) + log((yr2-cone)/yr2) * Cc2(1)

    do j = 2, n+1
      C0 = C0 + ((-1)**(1+j) * yr1**(1-j) - (1-yr1)**(1-j))/(cone*(j-1)) * Cc1(j) &
              + ((-1)**(1+j) * yr2**(1-j) - (1-yr2)**(1-j))/(cone*(j-1)) * Cc2(j)
    end do

    C0 = C0 / (p2 * (1+n))

  end function C0_n_m0m1m1

end module c0_m0m1m1_DP

module c0_m0m1m1_QP
  use triangle_aux_QP, only: target_precision,qp,cone,cnul,choose
  implicit none
  ! for values z=m2/p2 > THRESHOLD expansion formulas are used
  real(qp), parameter :: C_0mm_exp_thr = real(2,kind=qp)

  contains

  function Crec(yk,yr,n)
    complex(qp), intent(in) :: yk,yr
    integer,       intent(in) :: n
    complex(qp)             :: Crec(n+1)
    integer                   :: i,l

    Crec(n+1) = yk**n * (cone-yk)**(n+1) / (yk-yr)**(n+1)
    do l = n, 1, -1
      Crec(l) = cnul
      do i = 1, n+1-l
        Crec(l) = Crec(l) + Crec(l+i) * lambda(yk,yr,i,n)
      end do
      Crec(l) = Crec(l) / (1+n-l)

    end do

    contains

    function lambda(yk,yr,i,n) result(lam)
      complex(qp), intent(in) :: yk,yr
      integer,       intent(in) :: i,n
      complex(qp)             :: lam

      lam = (1+n)*(cone/(yr-yk)**i - cone/(cone-yk)**i)-n/(-yk)**i

    end function lambda

  end function Crec


  function C0_n_m0m1m1(p2,m2,n) result(C0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: z0,z1,yr1,yr2,C0,loss,Cc1(n+1),Cc2(n+1)
    integer                   :: j

    z0 = m2(1)/p2
    z1 = m2(2)/p2

    yr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    yr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    Cc1 = Crec(yr1,yr2,n)
    Cc2 = Crec(yr2,yr1,n)

    C0 = log((yr1-cone)/yr1) * Cc1(1) + log((yr2-cone)/yr2) * Cc2(1)

    do j = 2, n+1
      C0 = C0 + ((-1)**(1+j) * yr1**(1-j) - (1-yr1)**(1-j))/(cone*(j-1)) * Cc1(j) &
              + ((-1)**(1+j) * yr2**(1-j) - (1-yr2)**(1-j))/(cone*(j-1)) * Cc2(j)
    end do

    C0 = C0 / (p2 * (1+n))

  end function C0_n_m0m1m1

end module c0_m0m1m1_QP
