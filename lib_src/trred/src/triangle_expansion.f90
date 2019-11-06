!******************************************************************************!
!                                                                              !
!    triangle_expansion.f90                                                    !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module triangle_expansion_DP
  use c0_0mm_DP, only: C0_n_0mm
  use c0_m00_DP, only: C0_n_m00, C0_n_m00_EP1
  use c0_mmm_DP, only: C0_n_mmm, C0_n_mmm_init, C0_n_mmm_update
  use c0_000_DP, only: C0_n_000, C0_n_000_EP1
  use b0_mm_DP, only: B0_n_mm, B0_n_mm_init, B0_n_mm_update
  use b0_DP, only: B0_0_n, B0_n
  use b0_m0m1_DP, only: B0_n_m0m1
  use triangle_aux_DP, only: dp,target_precision,cone,cnul,rnul,muUV2,muIR2

  implicit none
  integer, parameter :: max_terms = 300
  integer, save      :: errflag = 0

  ! error flags:
  ! 1: write (*,*) 'ERROR: called coeff_template with invalid d=', d
  ! 2: write (*,*) 'ERROR: called coeff_template with offset<0'
  ! 4: write (*,*) 'ERROR: called list_template with zero clist.'
  ! 8: write (*,*) 'ERROR: called matrix_template with zero cmatrix.'
  ! 16: write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'

  abstract interface
    function coeff_func(p2,m2,n)
      use triangle_aux_DP, only: dp
      implicit none
      complex(dp), intent(in) :: p2,m2(:)
      integer,       intent(in) :: n
      complex(dp)             :: coeff_func
    end function coeff_func
  end interface

  contains

  function get_errflag()
    integer :: get_errflag
    get_errflag = errflag
  end function get_errflag

  subroutine reset_errflag
    errflag = 0
  end subroutine reset_errflag

  function coeff_template(p2,m2,d,offset,coeff_init,coeff_update) result(coeff)
    ! Computes the coefficient function func(p2,d,m2,offset) :=
    ! \sum_{n=offset}^\infty d^{n-offset} func^(n)
    complex(dp),                  intent(in) :: m2(:)
    real(dp),                     intent(in) :: p2,d
    integer, optional,              intent(in) :: offset
    procedure(coeff_func), pointer, intent(in) :: coeff_init,coeff_update
    complex(dp) :: coeff,coeffn
    integer    :: n,m

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called coeff_template with invalid d=', d
      errflag = ior(1, errflag)
      return
    end if

    if (present(offset)) then
      if (offset .lt. 0) then
        !write (*,*) 'ERROR: called coeff_template with offset<0'
        errflag = ior(2, errflag)
        return
      end if
      m = offset
    else
      m = 0
    end if

    n = m
    coeff = coeff_init(cmplx(p2,kind=dp),m2,n)

    n = n + 1
    coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=dp),m2,n)

    do while (abs(coeffn/coeff) .gt. target_precision)
      coeff = coeff + coeffn
      n = n + 1
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16, errflag)
        coeffn = cnul
        exit
      end if
      coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=dp),m2,n)
    end do
    coeff = coeff + coeffn

  end function coeff_template

  function list_template(p2,m2,d,clist,list_init,list_update) result(list)
    ! Computes \sum_i clist(i)*func(p2,d,m2,i) in an efficient way,
    ! computing equal coefficients only once.
    complex(dp),                  intent(in) :: m2(:)
    real(dp),                     intent(in) :: p2,d
    complex(dp),                  intent(in) :: clist(0:)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(dp) :: list,listn
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called list_template with invalid d=',d
      errflag = ior(1, errflag)
      list = 0
      return
    end if

    do m = 1, size(clist)
      if (clist(m-1) .ne. 0._dp) then
        exit
      end if
    end do
    if (m .eq. size(clist) .and. clist(m-1) .eq. 0._dp  ) then
      !write (*,*) 'ERROR: called list_template with zero clist.'
      errflag = ior(4,errflag)
      list = 0
    end if
    m = m - 1

    n = m
    list = clist(m)*list_init(cmplx(p2,kind=dp),m2,n)

    n = n + 1
    listn = 0._dp
    do k = m, min(n,size(clist)-1)
      if (clist(k) .ne. 0._dp) then
        listn = listn + d**(n-k)*clist(k)
      end if
    end do
    listn = listn*list_update(cmplx(p2,kind=dp),m2,n)

    ! the condition n .lt. size(clist)-1 ensures that for d=0 all coefficients
    ! are taken into account. For instance for a list with a zero
    ! clist=[1.,1.,0.,1.] this would otherwise lead to a mistake.
    do while (abs(listn/list) .gt. target_precision .or. n .lt. size(clist)-1)
      list = list + listn
      n = n + 1
      listn = cnul
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16, errflag)
        exit
      end if
      do k = m, min(n,size(clist)-1)
        if (clist(k) .ne. 0._dp) then
          listn = listn + d**(n-k)*clist(k)
        end if
      end do
      listn = listn*list_update(cmplx(p2,kind=dp),m2,n)
    end do
    list = list + listn

  end function list_template

  function matrix_template(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero,list_init,list_update) result(list)
    ! Computes \sum_i cmatrix(k,i)*func(p2,d,m2,i) in an efficient way,
    complex(dp),                  intent(in) :: m2(:)
    real(dp),                     intent(in) :: p2,d
    integer,                        intent(in) :: mcoef,ncoef
    complex(dp),                  intent(in) :: cmatrix(mcoef,0:ncoef-1)
    logical,                        intent(in) :: cmatrixzero(0:ncoef-1)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(dp) :: list(mcoef),listn(mcoef)
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called matrix_template with invalid d=',d
      errflag = ior(1, errflag)
      list = 0
      return
    end if

    do m = 1, ncoef
      if (sum(abs(cmatrix(:,m-1))) .ne. 0._dp) then
        exit
      end if
    end do
    if (m .eq. ncoef .and. cmatrixzero(m-1) .eqv. .true.) then
      !write (*,*) 'ERROR: called matrix_template with zero cmatrix.'
      errflag = ior(8, errflag)
      list = 0
      return
    end if
    m = m - 1

    n = m
    list(:) = cmatrix(:,m)*list_init(cmplx(p2,kind=dp),m2,n)

    n = n + 1
    listn(:) = 0._dp
    do k = m, min(n,ncoef-1)
      if (cmatrixzero(k) .eqv. .false.) then
        listn(:) = listn(:) + d**(n-k)*cmatrix(:,k)
      end if
    end do
    listn(:) = listn(:)*list_update(cmplx(p2,kind=dp),m2,n)

    do while (n .lt. ncoef-1 .or. .not. convergent(listn,list,mcoef))
      list(:) = list(:) + listn(:)
      n = n + 1
      listn(:) = cnul
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16, errflag)
        exit
      end if
      do k = m, min(n,ncoef-1)
        if (cmatrixzero(k) .eqv. .false.) then
          listn(:) = listn(:) + d**(n-k)*cmatrix(:,k)
        end if
      end do
      listn(:) = listn(:)*list_update(cmplx(p2,kind=dp),m2,n)
    end do
    list(:) = list(:) + listn(:)

    contains

    function convergent(c_next,c,mcoef)
      integer,       intent(in) :: mcoef
      complex(dp), intent(in) :: c_next(mcoef),c(mcoef)
      complex(dp) :: sumc
      logical :: convergent
      integer :: i

      convergent = .true.
      do i = 1, mcoef
        if (abs(c_next(i)/c(i)) .gt. target_precision) then
          convergent = .false.
          exit
        end if
      end do

    end function convergent

  end function matrix_template

  !!!!!!!!
  !  B0  !
  !!!!!!!!

  function B0d_coeff(p2,m2,d,offset) result(B0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function B0d_coeff

  function B0d_list(p2,m2,d,clist) result(B0)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    complex(dp), intent(in) :: clist(0:)
    complex(dp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = list_template(p2,[m2],d,clist,coeff_init,coeff_init)

  end function B0d_list

  function B0d_matrix(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero) result(B0)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d
    integer,       intent(in) :: mcoef,ncoef
    complex(dp), intent(in) :: cmatrix(mcoef,ncoef)
    logical,       intent(in) :: cmatrixzero(ncoef)
    complex(dp) :: B0(mcoef)
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = matrix_template(p2,[m2],d,mcoef,ncoef, &
                         cmatrix,cmatrixzero,coeff_init,coeff_init)

  end function B0d_matrix

  function B0d_mm_opt(p2,m2,d,offset) result(B0)
    complex(dp), intent(in)        :: m2
    real(dp),    intent(in)        :: p2,d
    integer,    intent(in), optional :: offset
    complex(dp)                    :: B0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function B0d_mm_opt

  function B0d_mm_opt_list(p2,m2,d,clist) result(B0)
    complex(dp), intent(in)      :: m2
    real(dp),    intent(in)      :: p2,d
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = list_template(p2,[m2],d,clist,coeff_init,coeff_update)

  end function B0d_mm_opt_list

  function B0d_0_coeff(p2,d,offset) result(B0)
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp)                       :: B0
    procedure(coeff_func),      pointer :: coeff_init

    coeff_init => B0_0_n
    B0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)
  end function B0d_0_coeff

  function B0d_0_list(p2,d,clist) result(B0)
    real(dp),    intent(in)      :: p2,d
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_0_n
    B0 = list_template(p2,[cnul],d,clist,coeff_init,coeff_init)

  end function B0d_0_list

  function B0d_m0m1_coeff(p2,m02,m12,d,offset) result(B0)
    use b0_m0m1_DP, only: B0_n_m0m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(dp)                       :: B0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = coeff_template(p2,[m02,m12],d,offset,coeff_init,coeff_init)

  end function B0d_m0m1_coeff

  function B0d_m0m1_list(p2,m02,m12,d,clist) result(B0)
    use b0_m0m1_DP, only: B0_n_m0m1
    real(dp),    intent(in)      :: p2,d
    complex(dp), intent(in)      :: m02,m12
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = list_template(p2,[m02,m12],d,clist,coeff_init,coeff_init)

  end function b0d_m0m1_list

  !!!!!!!!
  !  C0  !
  !!!!!!!!

  function C0d_0mm_coeff(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func),     pointer  :: coeff_init

    coeff_init => C0_n_0mm
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_0mm_coeff

  function C0d_0mm_matrix(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero) result(C0)
    use c0_m0m1m1_DP, only: C0_n_m0m1m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m2
    integer,       intent(in)           :: mcoef,ncoef
    complex(dp), intent(in)           :: cmatrix(mcoef,ncoef)
    logical,       intent(in)           :: cmatrixzero(mcoef)
    complex(dp)                       :: C0(mcoef)
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_0mm
    C0 = matrix_template(p2,[m2],d,mcoef,ncoef,cmatrix,cmatrixzero,coeff_init,coeff_init)

  end function C0d_0mm_matrix


  function C0d_m00_coeff(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_m00_coeff

  function C0d_m00_EP1_coeff(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00_EP1
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_m00_EP1_coeff

  function C0d_mmm_table(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => C0_n_mmm_init
    coeff_update => C0_n_mmm_update
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_mmm_table

  function C0d_000_coeff(p2,d,offset) result(C0)
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000
    C0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)

  end function C0d_000_coeff

  function C0d_000_EP1_coeff(p2,d,offset) result(C0)
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000_EP1
    C0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)

  end function C0d_000_EP1_coeff

  function C0d_m0m1m1_coeff(p2,m02,m12,d,offset) result(C0)
    use c0_m0m1m1_DP, only: C0_n_m0m1m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = coeff_template(p2,[m02,m12],d,offset,coeff_init,coeff_init)

  end function C0d_m0m1m1_coeff

  function C0d_m0m1m1_list(p2,m02,m12,d,clist) result(C0)
    use c0_m0m1m1_DP, only: C0_n_m0m1m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m02,m12
    complex(dp), intent(in)           :: clist(0:)
    complex(dp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = list_template(p2,[m02,m12],d,clist,coeff_init,coeff_init)

  end function C0d_m0m1m1_list

end module triangle_expansion_DP

module triangle_expansion_QP
  use c0_0mm_QP, only: C0_n_0mm
  use c0_m00_QP, only: C0_n_m00, C0_n_m00_EP1
  use c0_mmm_QP, only: C0_n_mmm, C0_n_mmm_init, C0_n_mmm_update
  use c0_000_QP, only: C0_n_000, C0_n_000_EP1
  use b0_mm_QP, only: B0_n_mm, B0_n_mm_init, B0_n_mm_update
  use b0_QP, only: B0_0_n, B0_n
  use b0_m0m1_QP, only: B0_n_m0m1
  use triangle_aux_QP, only: qp,target_precision,cone,cnul,rnul,muUV2,muIR2

  use triangle_expansion_DP, only: max_terms, errflag

  implicit none

  abstract interface
    function coeff_func(p2,m2,n)
      use triangle_aux_QP, only: qp
      implicit none
      complex(qp), intent(in) :: p2,m2(:)
      integer,       intent(in) :: n
      complex(qp)             :: coeff_func
    end function coeff_func
  end interface

  contains

  function coeff_template(p2,m2,d,offset,coeff_init,coeff_update) result(coeff)
    ! Computes the coefficient function func(p2,d,m2,offset) :=
    ! \sum_{n=offset}^\infty d^{n-offset} func^(n)
    complex(qp),                  intent(in) :: m2(:)
    real(qp),                     intent(in) :: p2,d
    integer, optional,              intent(in) :: offset
    procedure(coeff_func), pointer, intent(in) :: coeff_init,coeff_update
    complex(qp) :: coeff,coeffn
    integer    :: n,m

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called coeff_template with invalid d=', d
      errflag = ior(1, errflag)
      coeff = 0
      return
    end if

    if (present(offset)) then
      if (offset .lt. 0) then
        !write (*,*) 'ERROR: called coeff_template with offset<0'
        errflag = ior(2, errflag)
        coeff = 0
        return
      end if
      m = offset
    else
      m = 0
    end if

    n = m
    coeff = coeff_init(cmplx(p2,kind=qp),m2,n)

    n = n + 1
    coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=qp),m2,n)

    do while (abs(coeffn/coeff) .gt. target_precision)
      coeff = coeff + coeffn
      n = n + 1
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16, errflag)
        coeffn = cnul
        exit
      end if
      coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=qp),m2,n)
    end do
    coeff = coeff + coeffn

  end function coeff_template

  function list_template(p2,m2,d,clist,list_init,list_update) result(list)
    ! Computes \sum_i clist(i)*func(p2,d,m2,i) in an efficient way,
    ! computing equal coefficients only once.
    complex(qp),                  intent(in) :: m2(:)
    real(qp),                     intent(in) :: p2,d
    complex(qp),                  intent(in) :: clist(0:)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(qp) :: list,listn
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called list_template with invalid d=',d
      errflag = ior(1, errflag)
      list = 0
      return
    end if

    do m = 1, size(clist)
      if (clist(m-1) .ne. 0._qp) then
        exit
      end if
    end do
    if (m .eq. size(clist) .and. clist(m-1) .eq. 0._qp  ) then
      !write (*,*) 'ERROR: called list_template with zero clist.'
      errflag = ior(4, errflag)
      list = 0
    end if
    m = m - 1

    n = m
    list = clist(m)*list_init(cmplx(p2,kind=qp),m2,n)

    n = n + 1
    listn = 0._qp
    do k = m, min(n,size(clist)-1)
      if (clist(k) .ne. 0._qp) then
        listn = listn + d**(n-k)*clist(k)
      end if
    end do
    listn = listn*list_update(cmplx(p2,kind=qp),m2,n)

    ! the condition n .lt. size(clist)-1 ensures that for d=0 all coefficients
    ! are taken into account. For instance for a list with a zero
    ! clist=[1.,1.,0.,1.] this would otherwise lead to a mistake.
    do while (abs(listn/list) .gt. target_precision .or. n .lt. size(clist)-1)
      list = list + listn
      n = n + 1
      listn = cnul
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16,errflag)
        exit
      end if
      do k = m, min(n,size(clist)-1)
        if (clist(k) .ne. 0._qp) then
          listn = listn + d**(n-k)*clist(k)
        end if
      end do
      listn = listn*list_update(cmplx(p2,kind=qp),m2,n)
    end do
    list = list + listn

  end function list_template

  function matrix_template(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero,list_init,list_update) result(list)
    ! Computes \sum_i cmatrix(k,i)*func(p2,d,m2,i) in an efficient way,
    complex(qp),                  intent(in) :: m2(:)
    real(qp),                     intent(in) :: p2,d
    integer,                        intent(in) :: mcoef,ncoef
    complex(qp),                  intent(in) :: cmatrix(mcoef,0:ncoef-1)
    logical,                        intent(in) :: cmatrixzero(0:ncoef-1)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(qp) :: list(mcoef),listn(mcoef)
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      !write (*,*) 'ERROR: called matrix_template with invalid d=',d
      errflag = ior(1, errflag)
      list = 0
      return
    end if

    do m = 1, ncoef
      if (sum(abs(cmatrix(:,m-1))) .ne. 0._qp) then
        exit
      end if
    end do
    if (m .eq. ncoef .and. cmatrixzero(m-1) .eqv. .true.) then
      !write (*,*) 'ERROR: called matrix_template with zero cmatrix.'
      errflag = ior(8, errflag)
      list = 0
      return
    end if
    m = m - 1

    n = m
    list(:) = cmatrix(:,m)*list_init(cmplx(p2,kind=qp),m2,n)

    n = n + 1
    listn(:) = 0._qp
    do k = m, min(n,ncoef-1)
      if (cmatrixzero(k) .eqv. .false.) then
        listn(:) = listn(:) + d**(n-k)*cmatrix(:,k)
      end if
    end do
    listn(:) = listn(:)*list_update(cmplx(p2,kind=qp),m2,n)

    do while (n .lt. ncoef-1 .or. .not. convergent(listn,list,mcoef))
      list(:) = list(:) + listn(:)
      n = n + 1
      listn(:) = cnul
      if (n .gt. max_terms) then
        !write (*,*) 'ERROR: series did not converge within the first `max_terms` coefficients'
        errflag = ior(16, errflag)
        exit
      end if
      do k = m, min(n,ncoef-1)
        if (cmatrixzero(k) .eqv. .false.) then
          listn(:) = listn(:) + d**(n-k)*cmatrix(:,k)
        end if
      end do
      listn(:) = listn(:)*list_update(cmplx(p2,kind=qp),m2,n)
    end do
    list(:) = list(:) + listn(:)

    contains

    function convergent(c_next,c,mcoef)
      integer,       intent(in) :: mcoef
      complex(qp), intent(in) :: c_next(mcoef),c(mcoef)
      complex(qp) :: sumc
      logical :: convergent
      integer :: i

      convergent = .true.
      do i = 1, mcoef
        if (abs(c_next(i)/c(i)) .gt. target_precision) then
          convergent = .false.
          exit
        end if
      end do

    end function convergent

  end function matrix_template

  !!!!!!!!
  !  B0  !
  !!!!!!!!

  function B0d_coeff(p2,m2,d,offset) result(B0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function B0d_coeff

  function B0d_list(p2,m2,d,clist) result(B0)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    complex(qp), intent(in) :: clist(0:)
    complex(qp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = list_template(p2,[m2],d,clist,coeff_init,coeff_init)

  end function B0d_list

  function B0d_matrix(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero) result(B0)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d
    integer,       intent(in) :: mcoef,ncoef
    complex(qp), intent(in) :: cmatrix(mcoef,ncoef)
    logical,       intent(in) :: cmatrixzero(ncoef)
    complex(qp) :: B0(mcoef)
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = matrix_template(p2,[m2],d,mcoef,ncoef, &
                         cmatrix,cmatrixzero,coeff_init,coeff_init)

  end function B0d_matrix

  function B0d_mm_opt(p2,m2,d,offset) result(B0)
    complex(qp), intent(in)        :: m2
    real(qp),    intent(in)        :: p2,d
    integer,    intent(in), optional :: offset
    complex(qp)                    :: B0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function B0d_mm_opt

  function B0d_mm_opt_list(p2,m2,d,clist) result(B0)
    complex(qp), intent(in)      :: m2
    real(qp),    intent(in)      :: p2,d
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = list_template(p2,[m2],d,clist,coeff_init,coeff_update)

  end function B0d_mm_opt_list

  function B0d_0_coeff(p2,d,offset) result(B0)
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp)                       :: B0
    procedure(coeff_func),      pointer :: coeff_init

    coeff_init => B0_0_n
    B0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)
  end function B0d_0_coeff

  function B0d_0_list(p2,d,clist) result(B0)
    real(qp),    intent(in)      :: p2,d
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_0_n
    B0 = list_template(p2,[cnul],d,clist,coeff_init,coeff_init)

  end function B0d_0_list

  function B0d_m0m1_coeff(p2,m02,m12,d,offset) result(B0)
    use b0_m0m1_QP, only: B0_n_m0m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(qp)                       :: B0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = coeff_template(p2,[m02,m12],d,offset,coeff_init,coeff_init)

  end function B0d_m0m1_coeff

  function B0d_m0m1_list(p2,m02,m12,d,clist) result(B0)
    use b0_m0m1_QP, only: B0_n_m0m1
    real(qp),    intent(in)      :: p2,d
    complex(qp), intent(in)      :: m02,m12
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = list_template(p2,[m02,m12],d,clist,coeff_init,coeff_init)

  end function b0d_m0m1_list

  !!!!!!!!
  !  C0  !
  !!!!!!!!

  function C0d_0mm_coeff(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func),     pointer  :: coeff_init

    coeff_init => C0_n_0mm
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_0mm_coeff

  function C0d_0mm_matrix(p2,m2,d,mcoef,ncoef,cmatrix,cmatrixzero) result(C0)
    use c0_m0m1m1_QP, only: C0_n_m0m1m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m2
    integer,       intent(in)           :: mcoef,ncoef
    complex(qp), intent(in)           :: cmatrix(mcoef,ncoef)
    logical,       intent(in)           :: cmatrixzero(mcoef)
    complex(qp)                       :: C0(mcoef)
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_0mm
    C0 = matrix_template(p2,[m2],d,mcoef,ncoef,cmatrix,cmatrixzero,coeff_init,coeff_init)

  end function C0d_0mm_matrix


  function C0d_m00_coeff(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_m00_coeff

  function C0d_m00_EP1_coeff(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00_EP1
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_m00_EP1_coeff

  function C0d_mmm_table(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => C0_n_mmm_init
    coeff_update => C0_n_mmm_update
    C0 = coeff_template(p2,[m2],d,offset,coeff_init,coeff_init)

  end function C0d_mmm_table

  function C0d_000_coeff(p2,d,offset) result(C0)
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000
    C0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)

  end function C0d_000_coeff

  function C0d_000_EP1_coeff(p2,d,offset) result(C0)
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000_EP1
    C0 = coeff_template(p2,[cnul],d,offset,coeff_init,coeff_init)

  end function C0d_000_EP1_coeff

  function C0d_m0m1m1_coeff(p2,m02,m12,d,offset) result(C0)
    use c0_m0m1m1_QP, only: C0_n_m0m1m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = coeff_template(p2,[m02,m12],d,offset,coeff_init,coeff_init)

  end function C0d_m0m1m1_coeff

  function C0d_m0m1m1_list(p2,m02,m12,d,clist) result(C0)
    use c0_m0m1m1_QP, only: C0_n_m0m1m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m02,m12
    complex(qp), intent(in)           :: clist(0:)
    complex(qp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = list_template(p2,[m02,m12],d,clist,coeff_init,coeff_init)

  end function C0d_m0m1m1_list

end module triangle_expansion_QP
