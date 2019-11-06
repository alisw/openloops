!******************************************************************************!
!                                                                              !
!    c0_mmm.f90                                                                !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module c0_mmm_DP
  use triangle_aux_DP, only: target_precision,dp,cone,cnul,recursion_threshold,&
                          muUV2,muIR2
  implicit none

  complex(dp) :: edge_g(0:recursion_threshold),    &
                   edge_f(0:recursion_threshold),    &
                   edge_f0(0:recursion_threshold),   &
                   NewEdge_g(0:recursion_threshold), &
                   NewEdge_f(0:recursion_threshold), &
                   NewEdge_f0(0:recursion_threshold)
  contains

  function C0_n_mmm(p2,m2,n) result(Cn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp) :: y,Cn,g,f
    integer       :: k

    y = sqrt(1+4*m2(1)/p2)

    call gfnm(y,n,0,g,f)
    Cn = (-g*log((cone+y)/(-cone+y))/y + f)/p2

  end function C0_n_mmm

  recursive subroutine gfnm(y,n,m,g,f)
    complex(dp), intent(in) :: y
    integer,       intent(in) :: n,m
    complex(dp), intent(out) :: g,f
    complex(dp) :: gnm,gtmp1,ftmp1,gtmp2,ftmp2,ftmp3

    if (n .eq. 0) then
      g = cone/y**(2*m)
      f = cnul
    else
      call gfnm(y,n-1,m,gtmp1,ftmp1)
      call gfnm(y,n-1,m+1,gtmp2,ftmp2)
      ftmp3 = fnm(y,n-1,m+1)
      g = -1/cmplx(n+1,kind=dp)*((n-m-cone/2)*gtmp1+(cone/2+m)*gtmp2)
      f = -1/cmplx(n+1,kind=dp)*((n-m-cone/2)*ftmp1+(cone/2+m)*ftmp2-ftmp3)
    end if

  end subroutine gfnm

  recursive function fnm(y,n,m) result(f)
    complex(dp), intent(in) :: y
    integer,       intent(in) :: n,m
    complex(dp) :: f

    if (n .eq. 0) then
      f = -1/y**(2*m)
    else
      f = -((n-m)*fnm(y,n-1,m)+m*fnm(y,n-1,m+1))/cmplx(n+1,kind=dp)
    end if

  end function fnm

  !subroutine table_init(y, r, edge_g, edge_f, edge_f0)
  subroutine table_init(y, r)
    complex(dp), intent(in) :: y
    integer,       intent(in) :: r
    complex(dp) :: g(0:r,0:r),f(0:r,0:r), f0(0:r,0:r)
    !complex(dp), intent(out) :: edge_g(0:r), edge_f(0:r), edge_f0(0:r)
    integer :: m, n

    do m = 0,r
      g(0,m) = cone/y**(2*m)
      f(0,m) = cnul
      f0(0,m) = -1/y**(2*m)
    end do

    edge_g(0) = g(0,r)
    edge_f(0) = f(0,r)
    edge_f0(0) = f0(0,r)

    do n = 1, r
      do m = 0, (r-n)
        if ( n >= 2) then
          f0(n-1,m+1) = -((n-m-2)*f0(n-2,m+1) + (m+1)*f0(n-2,m+2))/cmplx(n,kind=dp)
        end if
        g(n,m) = -1/cmplx(n+1,kind=dp)*( (n-m-cone/2)*g(n-1,m) + (cone/2+m)*g(n-1,m+1) )
        f(n,m) = -1/cmplx(n+1,kind=dp)*( (n-m-cone/2)*f(n-1,m) + (cone/2+m)*f(n-1,m+1) - f0(n-1,m+1) )

        if (m == (r-n) ) then
          edge_g(n) = g(n,m)
          edge_f(n) = f(n,m)
          edge_f0(n-1) = f0(n-1,m+1)
        end if
      end do
    end do
  end subroutine table_init

  !subroutine table_update(y, r, edge_g, edge_f, edge_f0, NewEdge_g, NewEdge_f, NewEdge_f0)
  subroutine table_update(y, r)
    complex(dp), intent(in) :: y
    integer,       intent(in) :: r
    integer       :: m, n
    complex(dp) :: edge_g_tmp


    edge_g_tmp = cone/y**(2*r)
    !NewEdge_g(0) = cone/y**(2*r)
    NewEdge_g(0) = edge_g_tmp
    NewEdge_f(0) = cnul
    NewEdge_f0(0) = -cone/y**(2*r)
    do n = 1, r
      m = r - n
      if ( n >= 2) then
        NewEdge_f0(n-1) = -((n-m-2)*edge_f0(n-2) + (m+1)*NewEdge_f0(n-2))/n
      end if
      edge_g_tmp = -cone/(n+1)*((n-m-cone/2)*edge_g(n-1)+(cone/2+m)*edge_g_tmp)
      NewEdge_g(n) = edge_g_tmp
      !write(*,*) "NewEdge_g(n):", NewEdge_g(n),n
      NewEdge_f(n) = -cone/(n+1)*((n-m-cone/2)*edge_f(n-1)+(cone/2+m)*NewEdge_f(n-1)-NewEdge_f0(n-1) )
    end do
    edge_g = NewEdge_g
    !write(*,*) "edge_g(25):", edge_g(25)
    !write(*,*) "edge_g(26):", edge_g(26)
    !write(*,*) "edge_g(27):", edge_g(27)
    !write(*,*) "edge_g(28):", edge_g(28)
    !write(*,*) "edge_g(29):", edge_g(29)
    edge_f = NewEdge_f
    edge_f0 = NewEdge_f0
  end subroutine table_update

  function C0_n_mmm_init(p2,m2,n) result(C0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: C0,y

    y = sqrt(1+4*m2(1)/p2)
    call table_init(y, n)
    C0 = ( -edge_g(n) *log((cone+y)/(-cone+y))/y + edge_f(n) )/p2

  end function C0_n_mmm_init


  function C0_n_mmm_update(p2,m2,n) result(C0)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(dp)             :: C0,y

    y = sqrt(1+4*m2(1)/p2)
    !call table_update(y, n, edge_g, edge_f, edge_f0, NewEdge_g, NewEdge_f, NewEdge_f0)
    call table_update(y, n)
    C0 = ( -NewEdge_g(n) *log((cone+y)/(-cone+y))/y + NewEdge_f(n) )/p2

  end function C0_n_mmm_update


end module c0_mmm_DP

module c0_mmm_QP
  use triangle_aux_QP, only: target_precision,qp,cone,cnul,recursion_threshold,&
                          muUV2,muIR2
  implicit none

  complex(qp) :: edge_g(0:recursion_threshold),    &
                   edge_f(0:recursion_threshold),    &
                   edge_f0(0:recursion_threshold),   &
                   NewEdge_g(0:recursion_threshold), &
                   NewEdge_f(0:recursion_threshold), &
                   NewEdge_f0(0:recursion_threshold)
  contains

  function C0_n_mmm(p2,m2,n) result(Cn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp) :: y,Cn,g,f
    integer       :: k

    y = sqrt(1+4*m2(1)/p2)

    call gfnm(y,n,0,g,f)
    Cn = (-g*log((cone+y)/(-cone+y))/y + f)/p2

  end function C0_n_mmm

  recursive subroutine gfnm(y,n,m,g,f)
    complex(qp), intent(in) :: y
    integer,       intent(in) :: n,m
    complex(qp), intent(out) :: g,f
    complex(qp) :: gnm,gtmp1,ftmp1,gtmp2,ftmp2,ftmp3

    if (n .eq. 0) then
      g = cone/y**(2*m)
      f = cnul
    else
      call gfnm(y,n-1,m,gtmp1,ftmp1)
      call gfnm(y,n-1,m+1,gtmp2,ftmp2)
      ftmp3 = fnm(y,n-1,m+1)
      g = -1/cmplx(n+1,kind=qp)*((n-m-cone/2)*gtmp1+(cone/2+m)*gtmp2)
      f = -1/cmplx(n+1,kind=qp)*((n-m-cone/2)*ftmp1+(cone/2+m)*ftmp2-ftmp3)
    end if

  end subroutine gfnm

  recursive function fnm(y,n,m) result(f)
    complex(qp), intent(in) :: y
    integer,       intent(in) :: n,m
    complex(qp) :: f

    if (n .eq. 0) then
      f = -1/y**(2*m)
    else
      f = -((n-m)*fnm(y,n-1,m)+m*fnm(y,n-1,m+1))/cmplx(n+1,kind=qp)
    end if

  end function fnm

  !subroutine table_init(y, r, edge_g, edge_f, edge_f0)
  subroutine table_init(y, r)
    complex(qp), intent(in) :: y
    integer,       intent(in) :: r
    complex(qp) :: g(0:r,0:r),f(0:r,0:r), f0(0:r,0:r)
    !complex(qp), intent(out) :: edge_g(0:r), edge_f(0:r), edge_f0(0:r)
    integer :: m, n

    do m = 0,r
      g(0,m) = cone/y**(2*m)
      f(0,m) = cnul
      f0(0,m) = -1/y**(2*m)
    end do

    edge_g(0) = g(0,r)
    edge_f(0) = f(0,r)
    edge_f0(0) = f0(0,r)

    do n = 1, r
      do m = 0, (r-n)
        if ( n >= 2) then
          f0(n-1,m+1) = -((n-m-2)*f0(n-2,m+1) + (m+1)*f0(n-2,m+2))/cmplx(n,kind=qp)
        end if
        g(n,m) = -1/cmplx(n+1,kind=qp)*( (n-m-cone/2)*g(n-1,m) + (cone/2+m)*g(n-1,m+1) )
        f(n,m) = -1/cmplx(n+1,kind=qp)*( (n-m-cone/2)*f(n-1,m) + (cone/2+m)*f(n-1,m+1) - f0(n-1,m+1) )

        if (m == (r-n) ) then
          edge_g(n) = g(n,m)
          edge_f(n) = f(n,m)
          edge_f0(n-1) = f0(n-1,m+1)
        end if
      end do
    end do
  end subroutine table_init

  !subroutine table_update(y, r, edge_g, edge_f, edge_f0, NewEdge_g, NewEdge_f, NewEdge_f0)
  subroutine table_update(y, r)
    complex(qp), intent(in) :: y
    integer,       intent(in) :: r
    integer       :: m, n
    complex(qp) :: edge_g_tmp


    edge_g_tmp = cone/y**(2*r)
    !NewEdge_g(0) = cone/y**(2*r)
    NewEdge_g(0) = edge_g_tmp
    NewEdge_f(0) = cnul
    NewEdge_f0(0) = -cone/y**(2*r)
    do n = 1, r
      m = r - n
      if ( n >= 2) then
        NewEdge_f0(n-1) = -((n-m-2)*edge_f0(n-2) + (m+1)*NewEdge_f0(n-2))/n
      end if
      edge_g_tmp = -cone/(n+1)*((n-m-cone/2)*edge_g(n-1)+(cone/2+m)*edge_g_tmp)
      NewEdge_g(n) = edge_g_tmp
      !write(*,*) "NewEdge_g(n):", NewEdge_g(n),n
      NewEdge_f(n) = -cone/(n+1)*((n-m-cone/2)*edge_f(n-1)+(cone/2+m)*NewEdge_f(n-1)-NewEdge_f0(n-1) )
    end do
    edge_g = NewEdge_g
    !write(*,*) "edge_g(25):", edge_g(25)
    !write(*,*) "edge_g(26):", edge_g(26)
    !write(*,*) "edge_g(27):", edge_g(27)
    !write(*,*) "edge_g(28):", edge_g(28)
    !write(*,*) "edge_g(29):", edge_g(29)
    edge_f = NewEdge_f
    edge_f0 = NewEdge_f0
  end subroutine table_update

  function C0_n_mmm_init(p2,m2,n) result(C0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: C0,y

    y = sqrt(1+4*m2(1)/p2)
    call table_init(y, n)
    C0 = ( -edge_g(n) *log((cone+y)/(-cone+y))/y + edge_f(n) )/p2

  end function C0_n_mmm_init


  function C0_n_mmm_update(p2,m2,n) result(C0)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    complex(qp)             :: C0,y

    y = sqrt(1+4*m2(1)/p2)
    !call table_update(y, n, edge_g, edge_f, edge_f0, NewEdge_g, NewEdge_f, NewEdge_f0)
    call table_update(y, n)
    C0 = ( -NewEdge_g(n) *log((cone+y)/(-cone+y))/y + NewEdge_f(n) )/p2

  end function C0_n_mmm_update


end module c0_mmm_QP
