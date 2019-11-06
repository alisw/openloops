!******************************************************************************!
!                                                                              !
!    c0_000.f90                                                                !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module c0_000_DP
  use triangle_aux_DP, only: dp,cone,cnul,dir,muUV2,muIR2
  implicit none

  contains

  function C0_n_000(p2,m2,n) result(Cn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    integer       :: k
    complex(dp) :: Cn,sum

    sum = cnul
    if (n .ge. 0) then
      do k = 0, n-1
        sum = sum + ( (-1)**k )/(cone+k)*(-1)**(n-k+1)/cmplx(n-k,kind=dp)
      end do
      sum = sum + 2*(-1)**n/(cone+n)*log(p2/muIR2)
    end if

    Cn = - sum/(2*p2)

    if (dir .ne. 0) then
      Cn = Cn + dir*C0_n_000_EP1(p2,m2,n)
    end if

  end function C0_n_000

  function C0_n_000_EP1(p2,m2,n) result(Cn)
    complex(dp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    integer       :: k
    complex(dp) :: Cn

    Cn = -(-1)**(1+n)/(cone+n)/(p2)

  end function C0_n_000_EP1

end module c0_000_DP

module c0_000_QP
  use triangle_aux_QP, only: qp,cone,cnul,dir,muUV2,muIR2
  implicit none

  contains

  function C0_n_000(p2,m2,n) result(Cn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    integer       :: k
    complex(qp) :: Cn,sum

    sum = cnul
    if (n .ge. 0) then
      do k = 0, n-1
        sum = sum + ( (-1)**k )/(cone+k)*(-1)**(n-k+1)/cmplx(n-k,kind=qp)
      end do
      sum = sum + 2*(-1)**n/(cone+n)*log(p2/muIR2)
    end if

    Cn = - sum/(2*p2)

    if (dir .ne. 0) then
      Cn = Cn + dir*C0_n_000_EP1(p2,m2,n)
    end if

  end function C0_n_000

  function C0_n_000_EP1(p2,m2,n) result(Cn)
    complex(qp), intent(in) :: p2,m2(:)
    integer,       intent(in) :: n
    integer       :: k
    complex(qp) :: Cn

    Cn = -(-1)**(1+n)/(cone+n)/(p2)

  end function C0_n_000_EP1

end module c0_000_QP
