!******************************************************************************!
!                                                                              !
!    b0_mm.f90                                                                 !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module b0_mm_DP
  use triangle_aux_DP, only: target_precision,dp,i8,cone,cnul,choose,Lphi,zlogzf, &
                          acoth,factorial,gamma_int,recursion_threshold,duv,muUV2
  implicit none

  complex(dp), dimension(0:recursion_threshold+1) :: HyperPn,HyperP_diff

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Optimized Recursive Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function B0_n_mm(p2,m2,n) result(Bn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp) :: z,x,y,Bn,output

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y)) -2._dp) - log(p2/muUV2) + duv
    else
      call b0_table_init3(z, y, n, output)
      Bn = - 2*output
    end if

  end function B0_n_mm

  subroutine b0_table_init3(z, y, r, output)
    complex(dp), intent(in) :: y, z
    integer,       intent(in) :: r
    complex(dp), intent(out) :: output
    complex(dp), dimension(0:r,0:2*r) :: g_old, f0_old, g_new, f0_new
    integer ::  n, l, k

    do l = 0, r
      k = l + r
      g_old(l,k) = -(z**k)*y/(y**(2*l))*log((-1+y)/(1+y))/2
      f0_old(l,k) = (z**k)/(y**(2*l))
    end do

    do n = 1, r
      do l = 0, (r-n)
        k = l + (r-n)
        if (n == 1) then
          f0_new(l,k+1) = f0_old(l,k+1)
        end if
        if ( n >= 2) then
          f0_new(l,k+1) = cone/cmplx(n-1,kind=dp)*( -(k+1)/z*f0_old(l,k+2) + (4*l)/z*f0_old(l+1,k+3) )
          f0_old(l,k+1) = f0_new(l,k+1)
        end if
        g_new(l,k) = cone/cmplx(n,kind=dp)*( (-k)/z*g_old(l,k+1) + (4*l-2)/z*g_old(l+1,k+2) + 1/(2*z)*f0_new(l,k+1) )
        g_old(l,k) = g_new(l,k)
      end do
    end do
    output = g_new(0,0)
  end subroutine b0_table_init3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Optimized Explicit Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function B0_n_mm_init(p2,m2,n) result(Bn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp) :: Bn
    complex(dp) :: z,x,y,sum, Pdiff, Pn
    integer       :: k, a, j, i

    ! print *, "n=", n
    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y))-2._dp) - log(p2/muUV2) + duv
    else
      do j = 0, n-1
        HyperPn(j) = HyperGeo( cone/2, -j, 2, 4/(x*y**2) )
        HyperP_diff(j) = HyperGeo( cone/2, -j, 1, -4/x )
      end do
      sum = cnul
      do j = 0, n-1
        if (j == 0) then
          Pn = y
          Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        else
          Pn = 2*(-1)**j * x**(-1-j) / y * HyperPn(-(1-j))
          Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        end if
        sum = sum + Pdiff*Pn
      end do
      Pn = 2*(-1)**n * x**(-1-n)/y*HyperPn(-(1-n))

      Bn = - ( Pn*log( (1+y)/(-1+y) ) + sum ) / (z**n)
    end if

  end function B0_n_mm_init

  function B0_n_mm_update(p2,m2,n) result(Bn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp) :: Bn
    complex(dp) :: z,x,y,sum, Pdiff, Pn
    integer       :: k, a, j, i

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n > 0) then
      j = n-1
      HyperPn(j) = HyperGeo( cone/2, -j, 2, 4/(x*y**2) )
      HyperP_diff(j) = HyperGeo( cone/2, -j, 1, -4/x )
      sum = cnul
      do j = 0, n-1
        if (j == 0) then
          Pn = y
        else
          Pn = 2*(-1)**j * x**(-1-j) / y * HyperPn(-(1-j))
        end if
        Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        sum = sum + Pdiff*Pn
      end do
      Pn = 2*(-1)**n * x**(-1-n)/ y * HyperPn(-(1-n))

      Bn = - ( Pn*log( (1+y)/(-1+y) ) + sum ) / (z**n )
    end if

  end function B0_n_mm_update

  function HyperGeo(a,b,c,value) result(res)
    complex(dp), intent(in) :: a,value
    integer,       intent(in) :: b,c
    complex(dp) :: res
    complex(dp) :: a_tmp,b_tmp,c_tmp
    integer       :: k

    res = cone
    if (b<0) then
      a_tmp = a
      b_tmp = cmplx(b,kind=dp)
      c_tmp = cmplx(c,kind=dp)
      do k = 1, -b
        if (k >=2) then
          a_tmp = a_tmp*(a+k-1)
          b_tmp = b_tmp*(b+k-1)
          c_tmp = c_tmp*(c+k-1)
          res = res + a_tmp*b_tmp/c_tmp * (value)**k/factorial (k)
        else if (k == 1) then
          res = res + a_tmp*b_tmp/c_tmp * (value)
        end if
      end do
    end if
  end function HyperGeo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Naive Explicit Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function B0_n_mm_explict(p2,m2,n) result(Bn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp) :: z,x,y,Bn,g,f, part1, part2, prod, sum
    integer       :: k, a, j, i

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y)) -2._dp) - log(p2/muUV2) + duv
    else if (n ==1) then
      Bn = -(-P_n(1,x)*log((-1+y)/(1+y))+cone/x)/z
    else
      sum = cnul
      do j = 0, n-1
        sum = sum + P_diff(n,j,x)*P_n(j,x)
      end do
      Bn = - ( P_n(n,x)*log( (1+y)/(-1+y) ) + sum ) / (z**n * factorial (n))
    end if

  end function B0_n_mm_explict

    function P_n(n,x) result(Pn)
    integer,       intent(in) :: n
    complex(dp), intent(in) :: x
    complex(dp) :: b(0:n), c(0:n)
    complex(dp) :: a(0:n), Pn, sum, y
    integer       :: k

    ! x = cone/z
    y = sqrt(4._dp/x+1._dp)

    if (n==0) then
      Pn = y
    else
      a(0) = 1._dp
      b(0) = 1
      c(0) = 1
      a(1) = cone/2
      b(1) = 1 - n
      c(1) = 2
      sum = cnul
      do k = 0, n-1
        if (k >=2) then
          a(k) = a(k-1)*( a(1)+k-1 )
          b(k) = b(k-1)*( b(1)+k-1 )
          c(k) = c(k-1)*( c(1)+k-1 )
        end if
        sum = sum + a(k)*b(k)/c(k) * (4._dp/(x*y**2))**k/factorial (k)
      end do

      Pn = 2*(-1)**n * x**(-1-n) * gamma_int(n+1) * sum / y
    end if

  end function P_n

  function P_diff(n,j,x) result(P_diff_nj)
    integer,       intent(in) :: j, n
    complex(dp), intent(in) :: x
    complex(dp) :: b(0:n), c(0:n), n_i8, j_i8
    complex(dp) :: a(0:n), P_diff_nj, sum, y
    integer       :: k

    !n_i8 = n
    !j_i8 = j
    ! x = cone/z
    y = sqrt(4._dp/x+1._dp)

    a(0) = 1._dp
    b(0) = 1
    c(0) = 1
    a(1) = cone/2
    b(1) = 1 + j - n
    c(1) = 1
    sum = cnul
    do k = 0, n -j -1
      if (k >=2) then
        a(k) = a(k-1)*( a(1)+k-1 )
        b(k) = b(k-1)*( b(1)+k-1 )
        c(k) = c(k-1)*( c(1)+k-1 )
      end if
      sum = sum + a(k)*b(k)/c(k) * (-4._dp/x)**k/factorial (k)
    end do

    P_diff_nj = (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)*gamma_int(n+1)*sum/((j-n)*gamma_int(1+j))

  end function P_diff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Naive Recursive Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine b0_table_init(z, y, r, edge_g, edge_f, edge_f0)
    complex(dp), intent(in) :: y, z
    integer,       intent(in) :: r
    complex(dp), intent(out) :: edge_g(0:r), edge_f(0:r), edge_f0(0:r)
    complex(dp) :: g(-r:r,-r:r),f(-r:r,-r:r), f0(-r:r,-r:r)
    integer       ::  n, k

    do k = -r,r
      !g(0,k) = y/(y**(2*k))*acoth(y)
      g(0,k) = -y/(y**(2*k))*log((-1+y)/(1+y))/2
      f(0,k) = cnul
      f0(0,k) = cone/(y**(2*k))
    end do

    edge_g(0) = g(0,r)
    edge_f(0) = f(0,r)
    edge_f0(0) = f0(0,r)

    do n = 1, r
      do k = -(r-n), (r-n)
        if ( n >= 2) then
          f0(n-1,k) = cone/cmplx(n-1,kind=dp)* k/(4*z) * ( f0(n-2,k-1) - 2*f0(n-2,k) + f0(n-2,k+1) )
          f0(n-1,k-1) = cone/cmplx(n-1,kind=dp)* (k-1)/(4*z) * ( f0(n-2,k-2) - 2*f0(n-2,k-1) + f0(n-2,k) )
        end if
        g(n,k) = cone/cmplx(n,kind=dp)*( (k-cone/2._dp)/(4*z)*(g(n-1,k-1) - 2*g(n-1,k) + g(n-1,k+1)) &
                                         + 1/(8*z)*(f0(n-1,k-1) - f0(n-1,k)) )
        ! f(n,k) = 1/cmplx(n,kind=dp)*( (k-cone/2._dp) * (f(n-1,k) - f(n-1,k+1)) + cone/2._dp*f0(n-1,k) )

        if (k == (r-n) ) then
          ! print *, "n,k", n,k
          edge_g(n) = g(n,k)
          ! edge_f(n) = f(n,k)
          ! edge_f0(n-1) = f0(n-1,k)
        end if
      end do
    end do
  end subroutine b0_table_init

  subroutine b0_table_init2(z, y, r, output)
    complex(dp), intent(in)  :: y, z
    integer,       intent(in)  :: r
    complex(dp), intent(out) :: output
    complex(dp), dimension(0:r,0:r,0:2*r) :: g, f0
    integer ::  n, l, k

    do l = 0, r
      k = l + r
      !g(0,l,k) = (z**k)*y/(y**(2*l))*acoth(y)
      g(0,l,k) = -(z**k)*y/(y**(2*l))*log((-1+y)/(1+y))/2
      f0(0,l,k) = (z**k)/(y**(2*l))
    end do

    do n = 1, r
      do l = 0, (r-n)
        k = l + (r-n)
        if ( n >= 2) then
          f0(n-1,l,k+1) = cone/cmplx(n-1,kind=dp)*( -(k+1)/z*f0(n-2,l,k+2) + (4*l)/z*f0(n-2,l+1,k+3) )
        end if
        g(n,l,k) = cone/cmplx(n,kind=dp)*( (-k)/z*g(n-1,l,k+1) + (4*l-2)/z*g(n-1,l+1,k+2) + 1/(2*z)*f0(n-1,l,k+1) )
      end do
    end do
    output = g(r,0,0)
  end subroutine b0_table_init2


end module b0_mm_DP

module b0_mm_QP
  use triangle_aux_QP, only: target_precision,qp,i8,cone,cnul,choose,Lphi,zlogzf, &
                          acoth,factorial,gamma_int,recursion_threshold,duv,muUV2
  implicit none

  complex(qp), dimension(0:recursion_threshold+1) :: HyperPn,HyperP_diff

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Optimized Recursive Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function B0_n_mm(p2,m2,n) result(Bn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp) :: z,x,y,Bn,output

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y)) -2._qp) - log(p2/muUV2) + duv
    else
      call b0_table_init3(z, y, n, output)
      Bn = - 2*output
    end if

  end function B0_n_mm

  subroutine b0_table_init3(z, y, r, output)
    complex(qp), intent(in) :: y, z
    integer,       intent(in) :: r
    complex(qp), intent(out) :: output
    complex(qp), dimension(0:r,0:2*r) :: g_old, f0_old, g_new, f0_new
    integer ::  n, l, k

    do l = 0, r
      k = l + r
      g_old(l,k) = -(z**k)*y/(y**(2*l))*log((-1+y)/(1+y))/2
      f0_old(l,k) = (z**k)/(y**(2*l))
    end do

    do n = 1, r
      do l = 0, (r-n)
        k = l + (r-n)
        if (n == 1) then
          f0_new(l,k+1) = f0_old(l,k+1)
        end if
        if ( n >= 2) then
          f0_new(l,k+1) = cone/cmplx(n-1,kind=qp)*( -(k+1)/z*f0_old(l,k+2) + (4*l)/z*f0_old(l+1,k+3) )
          f0_old(l,k+1) = f0_new(l,k+1)
        end if
        g_new(l,k) = cone/cmplx(n,kind=qp)*( (-k)/z*g_old(l,k+1) + (4*l-2)/z*g_old(l+1,k+2) + 1/(2*z)*f0_new(l,k+1) )
        g_old(l,k) = g_new(l,k)
      end do
    end do
    output = g_new(0,0)
  end subroutine b0_table_init3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Optimized Explicit Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function B0_n_mm_init(p2,m2,n) result(Bn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp) :: Bn
    complex(qp) :: z,x,y,sum, Pdiff, Pn
    integer       :: k, a, j, i

    ! print *, "n=", n
    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y))-2._qp) - log(p2/muUV2) + duv
    else
      do j = 0, n-1
        HyperPn(j) = HyperGeo( cone/2, -j, 2, 4/(x*y**2) )
        HyperP_diff(j) = HyperGeo( cone/2, -j, 1, -4/x )
      end do
      sum = cnul
      do j = 0, n-1
        if (j == 0) then
          Pn = y
          Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        else
          Pn = 2*(-1)**j * x**(-1-j) / y * HyperPn(-(1-j))
          Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        end if
        sum = sum + Pdiff*Pn
      end do
      Pn = 2*(-1)**n * x**(-1-n)/y*HyperPn(-(1-n))

      Bn = - ( Pn*log( (1+y)/(-1+y) ) + sum ) / (z**n)
    end if

  end function B0_n_mm_init

  function B0_n_mm_update(p2,m2,n) result(Bn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp) :: Bn
    complex(qp) :: z,x,y,sum, Pdiff, Pn
    integer       :: k, a, j, i

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n > 0) then
      j = n-1
      HyperPn(j) = HyperGeo( cone/2, -j, 2, 4/(x*y**2) )
      HyperP_diff(j) = HyperGeo( cone/2, -j, 1, -4/x )
      sum = cnul
      do j = 0, n-1
        if (j == 0) then
          Pn = y
        else
          Pn = 2*(-1)**j * x**(-1-j) / y * HyperPn(-(1-j))
        end if
        Pdiff =  (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)/((j-n)) * HyperP_diff(-(1+j-n))
        sum = sum + Pdiff*Pn
      end do
      Pn = 2*(-1)**n * x**(-1-n)/ y * HyperPn(-(1-n))

      Bn = - ( Pn*log( (1+y)/(-1+y) ) + sum ) / (z**n )
    end if

  end function B0_n_mm_update

  function HyperGeo(a,b,c,value) result(res)
    complex(qp), intent(in) :: a,value
    integer,       intent(in) :: b,c
    complex(qp) :: res
    complex(qp) :: a_tmp,b_tmp,c_tmp
    integer       :: k

    res = cone
    if (b<0) then
      a_tmp = a
      b_tmp = cmplx(b,kind=qp)
      c_tmp = cmplx(c,kind=qp)
      do k = 1, -b
        if (k >=2) then
          a_tmp = a_tmp*(a+k-1)
          b_tmp = b_tmp*(b+k-1)
          c_tmp = c_tmp*(c+k-1)
          res = res + a_tmp*b_tmp/c_tmp * (value)**k/factorial (k)
        else if (k == 1) then
          res = res + a_tmp*b_tmp/c_tmp * (value)
        end if
      end do
    end if
  end function HyperGeo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Naive Explicit Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function B0_n_mm_explict(p2,m2,n) result(Bn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp) :: z,x,y,Bn,g,f, part1, part2, prod, sum
    integer       :: k, a, j, i

    z = m2(1)/p2
    x = cone/z
    y = sqrt(4*z + cone)

    if (n==0) then
      Bn = -(log(z) - y*log((-1+y)/(1+y)) -2._qp) - log(p2/muUV2) + duv
    else if (n ==1) then
      Bn = -(-P_n(1,x)*log((-1+y)/(1+y))+cone/x)/z
    else
      sum = cnul
      do j = 0, n-1
        sum = sum + P_diff(n,j,x)*P_n(j,x)
      end do
      Bn = - ( P_n(n,x)*log( (1+y)/(-1+y) ) + sum ) / (z**n * factorial (n))
    end if

  end function B0_n_mm_explict

    function P_n(n,x) result(Pn)
    integer,       intent(in) :: n
    complex(qp), intent(in) :: x
    complex(qp) :: b(0:n), c(0:n)
    complex(qp) :: a(0:n), Pn, sum, y
    integer       :: k

    ! x = cone/z
    y = sqrt(4._qp/x+1._qp)

    if (n==0) then
      Pn = y
    else
      a(0) = 1._qp
      b(0) = 1
      c(0) = 1
      a(1) = cone/2
      b(1) = 1 - n
      c(1) = 2
      sum = cnul
      do k = 0, n-1
        if (k >=2) then
          a(k) = a(k-1)*( a(1)+k-1 )
          b(k) = b(k-1)*( b(1)+k-1 )
          c(k) = c(k-1)*( c(1)+k-1 )
        end if
        sum = sum + a(k)*b(k)/c(k) * (4._qp/(x*y**2))**k/factorial (k)
      end do

      Pn = 2*(-1)**n * x**(-1-n) * gamma_int(n+1) * sum / y
    end if

  end function P_n

  function P_diff(n,j,x) result(P_diff_nj)
    integer,       intent(in) :: j, n
    complex(qp), intent(in) :: x
    complex(qp) :: b(0:n), c(0:n), n_i8, j_i8
    complex(qp) :: a(0:n), P_diff_nj, sum, y
    integer       :: k

    !n_i8 = n
    !j_i8 = j
    ! x = cone/z
    y = sqrt(4._qp/x+1._qp)

    a(0) = 1._qp
    b(0) = 1
    c(0) = 1
    a(1) = cone/2
    b(1) = 1 + j - n
    c(1) = 1
    sum = cnul
    do k = 0, n -j -1
      if (k >=2) then
        a(k) = a(k-1)*( a(1)+k-1 )
        b(k) = b(k-1)*( b(1)+k-1 )
        c(k) = c(k-1)*( c(1)+k-1 )
      end if
      sum = sum + a(k)*b(k)/c(k) * (-4._qp/x)**k/factorial (k)
    end do

    P_diff_nj = (-1)**(-j+n) *x**(j-n) *y**(1+2*j-2*n)*gamma_int(n+1)*sum/((j-n)*gamma_int(1+j))

  end function P_diff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         Naive Recursive Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine b0_table_init(z, y, r, edge_g, edge_f, edge_f0)
    complex(qp), intent(in) :: y, z
    integer,       intent(in) :: r
    complex(qp), intent(out) :: edge_g(0:r), edge_f(0:r), edge_f0(0:r)
    complex(qp) :: g(-r:r,-r:r),f(-r:r,-r:r), f0(-r:r,-r:r)
    integer       ::  n, k

    do k = -r,r
      !g(0,k) = y/(y**(2*k))*acoth(y)
      g(0,k) = -y/(y**(2*k))*log((-1+y)/(1+y))/2
      f(0,k) = cnul
      f0(0,k) = cone/(y**(2*k))
    end do

    edge_g(0) = g(0,r)
    edge_f(0) = f(0,r)
    edge_f0(0) = f0(0,r)

    do n = 1, r
      do k = -(r-n), (r-n)
        if ( n >= 2) then
          f0(n-1,k) = cone/cmplx(n-1,kind=qp)* k/(4*z) * ( f0(n-2,k-1) - 2*f0(n-2,k) + f0(n-2,k+1) )
          f0(n-1,k-1) = cone/cmplx(n-1,kind=qp)* (k-1)/(4*z) * ( f0(n-2,k-2) - 2*f0(n-2,k-1) + f0(n-2,k) )
        end if
        g(n,k) = cone/cmplx(n,kind=qp)*( (k-cone/2._qp)/(4*z)*(g(n-1,k-1) - 2*g(n-1,k) + g(n-1,k+1)) &
                                         + 1/(8*z)*(f0(n-1,k-1) - f0(n-1,k)) )
        ! f(n,k) = 1/cmplx(n,kind=qp)*( (k-cone/2._qp) * (f(n-1,k) - f(n-1,k+1)) + cone/2._qp*f0(n-1,k) )

        if (k == (r-n) ) then
          ! print *, "n,k", n,k
          edge_g(n) = g(n,k)
          ! edge_f(n) = f(n,k)
          ! edge_f0(n-1) = f0(n-1,k)
        end if
      end do
    end do
  end subroutine b0_table_init

  subroutine b0_table_init2(z, y, r, output)
    complex(qp), intent(in)  :: y, z
    integer,       intent(in)  :: r
    complex(qp), intent(out) :: output
    complex(qp), dimension(0:r,0:r,0:2*r) :: g, f0
    integer ::  n, l, k

    do l = 0, r
      k = l + r
      !g(0,l,k) = (z**k)*y/(y**(2*l))*acoth(y)
      g(0,l,k) = -(z**k)*y/(y**(2*l))*log((-1+y)/(1+y))/2
      f0(0,l,k) = (z**k)/(y**(2*l))
    end do

    do n = 1, r
      do l = 0, (r-n)
        k = l + (r-n)
        if ( n >= 2) then
          f0(n-1,l,k+1) = cone/cmplx(n-1,kind=qp)*( -(k+1)/z*f0(n-2,l,k+2) + (4*l)/z*f0(n-2,l+1,k+3) )
        end if
        g(n,l,k) = cone/cmplx(n,kind=qp)*( (-k)/z*g(n-1,l,k+1) + (4*l-2)/z*g(n-1,l+1,k+2) + 1/(2*z)*f0(n-1,l,k+1) )
      end do
    end do
    output = g(r,0,0)
  end subroutine b0_table_init2


end module b0_mm_QP
