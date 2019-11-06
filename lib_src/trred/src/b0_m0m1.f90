!******************************************************************************!
!                                                                              !
!    b0_m0m1.f90                                                               !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module b0_m0m1_DP
  use triangle_aux_DP, only: target_precision,dp,cone,cnul,choose,pi,cima,duv,muUV2
  implicit none

  contains

  function Brec(yk,yr,n)
    complex(dp), intent(in) :: yk,yr
    integer,       intent(in) :: n
    complex(dp)             :: Brec(n+1)
    integer                   :: i,l

    Brec(n) = yk**n * (cone-yk)**(n) / (yk-yr)**(n)
    do l = n-1, 1, -1
      Brec(l) = cnul
      do i = 1, n-l
        Brec(l) = Brec(l) + Brec(l+i) * lambda(yk,yr,i,n)
      end do
      Brec(l) = Brec(l) / (n-l)

    end do

    contains

    function lambda(yk,yr,i,n) result(lam)
      complex(dp), intent(in) :: yk,yr
      integer,       intent(in) :: i,n
      complex(dp)             :: lam

      lam = n*(cone/(yr-yk)**i - cone/(cone-yk)**i-cone/(-yk)**i)

    end function lambda

  end function Brec

  function B0_0_m0m1(p2,m2) result(B0)
    complex(dp), intent(in) :: p2,m2(0:)
    complex(dp)             :: m0,m1,B0,z0,z1,xr1,xr2
    integer                   :: j

    m0 = sqrt(m2(0))
    m1 = sqrt(m2(1))
    z0 = m2(0)/p2
    z1 = m2(1)/p2

    xr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    xr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    !B0 = 2*cone - log(p2/muUV2) &
    !   + (xr1 - 1)*log(1 - xr1) + (xr2 - 1)*log(1 - xr2) &
    !   - xr1*log(-xr1) - xr2*log(-xr2) + pi*cima
    ! stable imaginary part
    B0 = 2*cone - log(p2/muUV2) + duv &
       + (xr1 - 1)*log((1 - xr1)/(-xr1)) + (xr2 - 1)*log((1 - xr2)/(-xr2)) &
       - log(-xr1*xr2)


  end function B0_0_m0m1


  function B0_n_m0m1(p2,m2,n) result(B0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: z0,z1,xr1,xr2,B0,loss,Bc1(n),Bc2(n)
    integer                   :: j

    if (n .eq. 0) then
      B0 = B0_0_m0m1(p2,m2)
      return
    end if

    z0 = m2(1)/p2
    z1 = m2(2)/p2

    xr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    xr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    Bc1 = Brec(xr1,xr2,n)
    Bc2 = Brec(xr2,xr1,n)

    B0 = (-cone)**n + log((xr1-cone)/xr1) * Bc1(1) + log((xr2-cone)/xr2) * Bc2(1)

    do j = 2, n
      B0 = B0 + ((-1)**(1+j) * xr1**(1-j) - (1-xr1)**(1-j))/(cone*(j-1)) * Bc1(j) &
              + ((-1)**(1+j) * xr2**(1-j) - (1-xr2)**(1-j))/(cone*(j-1)) * Bc2(j)
    end do

    B0 = B0 / (n)

  end function B0_n_m0m1

end module b0_m0m1_DP

module b0_m0m1_QP
  use triangle_aux_QP, only: target_precision,qp,cone,cnul,choose,pi,cima,duv,muUV2
  implicit none

  contains

  function Brec(yk,yr,n)
    complex(qp), intent(in) :: yk,yr
    integer,       intent(in) :: n
    complex(qp)             :: Brec(n+1)
    integer                   :: i,l

    Brec(n) = yk**n * (cone-yk)**(n) / (yk-yr)**(n)
    do l = n-1, 1, -1
      Brec(l) = cnul
      do i = 1, n-l
        Brec(l) = Brec(l) + Brec(l+i) * lambda(yk,yr,i,n)
      end do
      Brec(l) = Brec(l) / (n-l)

    end do

    contains

    function lambda(yk,yr,i,n) result(lam)
      complex(qp), intent(in) :: yk,yr
      integer,       intent(in) :: i,n
      complex(qp)             :: lam

      lam = n*(cone/(yr-yk)**i - cone/(cone-yk)**i-cone/(-yk)**i)

    end function lambda

  end function Brec

  function B0_0_m0m1(p2,m2) result(B0)
    complex(qp), intent(in) :: p2,m2(0:)
    complex(qp)             :: m0,m1,B0,z0,z1,xr1,xr2
    integer                   :: j

    m0 = sqrt(m2(0))
    m1 = sqrt(m2(1))
    z0 = m2(0)/p2
    z1 = m2(1)/p2

    xr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    xr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    !B0 = 2*cone - log(p2/muUV2) &
    !   + (xr1 - 1)*log(1 - xr1) + (xr2 - 1)*log(1 - xr2) &
    !   - xr1*log(-xr1) - xr2*log(-xr2) + pi*cima
    ! stable imaginary part
    B0 = 2*cone - log(p2/muUV2) + duv &
       + (xr1 - 1)*log((1 - xr1)/(-xr1)) + (xr2 - 1)*log((1 - xr2)/(-xr2)) &
       - log(-xr1*xr2)


  end function B0_0_m0m1


  function B0_n_m0m1(p2,m2,n) result(B0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: z0,z1,xr1,xr2,B0,loss,Bc1(n),Bc2(n)
    integer                   :: j

    if (n .eq. 0) then
      B0 = B0_0_m0m1(p2,m2)
      return
    end if

    z0 = m2(1)/p2
    z1 = m2(2)/p2

    xr1 = (cone + z0 - z1 - sqrt(4*z1 + (cone + z0 - z1)**2))/2
    xr2 = (cone + z0 - z1 + sqrt(4*z1 + (cone + z0 - z1)**2))/2

    Bc1 = Brec(xr1,xr2,n)
    Bc2 = Brec(xr2,xr1,n)

    B0 = (-cone)**n + log((xr1-cone)/xr1) * Bc1(1) + log((xr2-cone)/xr2) * Bc2(1)

    do j = 2, n
      B0 = B0 + ((-1)**(1+j) * xr1**(1-j) - (1-xr1)**(1-j))/(cone*(j-1)) * Bc1(j) &
              + ((-1)**(1+j) * xr2**(1-j) - (1-xr2)**(1-j))/(cone*(j-1)) * Bc2(j)
    end do

    B0 = B0 / (n)

  end function B0_n_m0m1

end module b0_m0m1_QP
