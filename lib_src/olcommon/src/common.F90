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


! ol_dilog configuration
! stop series expansion when the result does not change anymore
#define COMPARE
! check expansion depth required to achieve full precision;
! store the number of points which required n terms in bin(n)
!#CHECK_DEPTH


#ifdef PRECISION_dp

module ol_generic
  ! precision independent generic routines
  implicit none

  character(len=26), private, parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), private, parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface to_string
    module procedure integer_to_string, integer1_to_string, integer2_to_string, &
          & double_to_string, complex_to_string, single_to_string, &
          & integerlist_to_string, doublelist_to_string, complexlist_to_string, &
          logical_to_string
  end interface to_string

  interface to_int
    module procedure string_to_integer
  end interface to_int

  interface relative_deviation
    module procedure relative_deviation_dp, relative_deviation_qp
  end interface relative_deviation

  contains

  function integer_to_string(x)
    implicit none
    integer :: x
    character(12) :: integer_to_string
    write(integer_to_string,*) x
    integer_to_string = adjustl(integer_to_string)
  end function integer_to_string

  function integer1_to_string(x)
    use KIND_TYPES, only: intkind1
    implicit none
    integer(intkind1) :: x
    character(12) :: integer1_to_string
    write(integer1_to_string,*) x
    integer1_to_string = adjustl(integer1_to_string)
  end function integer1_to_string

  function integer2_to_string(x)
    use KIND_TYPES, only: intkind2
    implicit none
    integer(intkind2) :: x
    character(12) :: integer2_to_string
    write(integer2_to_string,*) x
    integer2_to_string = adjustl(integer2_to_string)
  end function integer2_to_string

  function integerlist_to_string(x,del,sep)
    implicit none
    integer :: x(:)
    character(13*size(x)+1) :: integerlist_to_string
    logical, optional, intent(in) :: del
    character(1), optional, intent(in) :: sep
    character(1) :: seperator
    integer :: k
    if (present(sep)) then
      seperator = sep
    else
      seperator = ","
    end if
    integerlist_to_string = ""
    if (present(del)) then
      if (del) integerlist_to_string = "["
    end if
    if (size(x) /= 0) integerlist_to_string = trim(integerlist_to_string) // trim(integer_to_string(x(1)))
    do k = 2, size(x)
      integerlist_to_string = trim(integerlist_to_string) // seperator // trim(integer_to_string(x(k)))
    end do
    if (present(del)) then
      if (del) integerlist_to_string = trim(integerlist_to_string) // "]"
    end if
  end function integerlist_to_string

  function logical_to_string(x)
    implicit none
    logical :: x
    character(5) :: logical_to_string
    !write(integer_to_string,*) x
    if (x .eqv. .true.) then
      logical_to_string = "True"
    else
      logical_to_string = "False"
    end if
  end function logical_to_string

  function doublelist_to_string(x,del,sep)
    use KIND_TYPES, only: DREALKIND
    implicit none
    real(DREALKIND) :: x(:)
    character(29*size(x)+1) :: doublelist_to_string
    logical, optional, intent(in) :: del
    character(1), optional, intent(in) :: sep
    character(1) :: seperator
    integer :: k
    if (present(sep)) then
      seperator = sep
    else
      seperator = ","
    end if
    doublelist_to_string = ""
    if (present(del)) then
      if (del) doublelist_to_string = "["
    end if
    if (size(x) /= 0) doublelist_to_string = trim(doublelist_to_string) // trim(double_to_string(x(1)))
    do k = 2, size(x)
      doublelist_to_string = trim(doublelist_to_string) // seperator // trim(double_to_string(x(k)))
    end do
    if (present(del)) then
      if (del) doublelist_to_string = trim(doublelist_to_string) // "]"
    end if
  end function doublelist_to_string

  function complexlist_to_string(x,del,sep)
    use KIND_TYPES, only: DREALKIND
    implicit none
    complex(DREALKIND) :: x(:)
    character(60*size(x)+1) :: complexlist_to_string
    logical, optional, intent(in) :: del
    character(1), optional, intent(in) :: sep
    character(1) :: seperator
    integer :: k
    if (present(sep)) then
      seperator = sep
    else
      seperator = ","
    end if
    complexlist_to_string = ""
    if (present(del)) then
      if (del) complexlist_to_string = "["
    end if
    if (size(x) /= 0) complexlist_to_string = trim(complexlist_to_string) // trim(complex_to_string(x(1)))
    do k = 2, size(x)
      complexlist_to_string = trim(complexlist_to_string) // seperator // trim(complex_to_string(x(k)))
    end do
    if (present(del)) then
      if (del) complexlist_to_string = trim(complexlist_to_string) // "]"
    end if
  end function complexlist_to_string

  function string_to_integerlist(c_in)
  ! convert a comma/space/slash separated string of numbers into an array of integers
    implicit none
    character(len=*), intent(in) :: c_in
    character(len(c_in)+1) :: c
    integer, allocatable :: string_to_integerlist(:)
    integer i, n, pos1
    logical last_seperator

    c = c_in // " "

    n=0
    pos1=0
    last_seperator =  .false.
    do i = 1, len(c)
      if (c(i:i) == "[" .or. c(i:i) == "]") c(i:i) = " "

      if (c(i:i) == ',' .or. c(i:i) == ' ' .or. c(i:i) == "/" ) then
        if (last_seperator)  then
          pos1 = i
          cycle
        end if
        n = n+1
        pos1 = i
        last_seperator = .true.
      else
        last_seperator = .false.
      end if
    end do

    allocate(string_to_integerlist(n))

    n=0
    pos1=0
    last_seperator =  .false.
    do i = 1, len(c)
      if (c(i:i) == ',' .or. c(i:i) == ' ' .or. c(i:i) == "/") then
        if (last_seperator)  then
          pos1 = i
          cycle
        end if
        n = n+1
        string_to_integerlist(n) = string_to_integer(c(pos1+1:i-1))
        pos1 = i
        last_seperator = .true.
      else
        last_seperator = .false.
      end if
    end do
  end function string_to_integerlist


  function double_to_string(x)
    use KIND_TYPES, only: DREALKIND
    implicit none
    real(DREALKIND) :: x
    character(28) :: double_to_string
    character(26) :: str
    integer :: k, epos, mantissaendpos
    logical :: leading0
    write(str,*) x
    str = adjustl(str)
    epos = index(to_lowercase(str), "e")
    if (epos == 0) then
      mantissaendpos = len(trim(str))
    else
      mantissaendpos = epos-1
    end if
    double_to_string = str(1:mantissaendpos)
    do k = mantissaendpos, 1, -1
      if (double_to_string(k:k) == "0") then
        double_to_string(k:k) = " "
      else
        exit
      end if
    end do
    if (epos /= 0) then
      double_to_string = trim(double_to_string) // "e"
      if (str(epos+1:epos+1) == "+" .or. str(epos+1:epos+1) == "-") then
        double_to_string = trim(double_to_string) // str(epos+1:epos+1)
        epos = epos + 1
      end if
      leading0 = .true.
      do k = epos+1, len(str)
        if (str(k:k) == "0" .and. leading0) then
          cycle
        else
          leading0 = .false.
          double_to_string = trim(double_to_string) // str(k:k)
        end if
      end do
    end if
    double_to_string = trim(double_to_string) // "_dp"
  end function double_to_string

  function single_to_string(x)
    implicit none
    real(selected_real_kind(6)) :: x
    character(20) :: single_to_string
    character(18) :: str
    integer :: k, epos, mantissaendpos
    logical :: leading0
    write(str,*) x
    str = adjustl(str)
    epos = index(to_lowercase(str), "e")
    if (epos == 0) then
      mantissaendpos = len(trim(str))
    else
      mantissaendpos = epos-1
    end if
    single_to_string = str(1:mantissaendpos)
    do k = mantissaendpos, 1, -1
      if (single_to_string(k:k) == "0") then
        single_to_string(k:k) = " "
      else
        exit
      end if
    end do
    if (epos /= 0) then
      single_to_string = trim(single_to_string) // "e"
      if (str(epos+1:epos+1) == "+" .or. str(epos+1:epos+1) == "-") then
        single_to_string = trim(single_to_string) // str(epos+1:epos+1)
        epos = epos + 1
      end if
      leading0 = .true.
      do k = epos+1, len(str)
        if (str(k:k) == "0" .and. leading0) then
          cycle
        else
          leading0 = .false.
          single_to_string = trim(single_to_string) // str(k:k)
        end if
      end do
    end if
    single_to_string = trim(single_to_string) // "_sp"
  end function single_to_string



  function complex_to_string(x)
    use KIND_TYPES, only: DREALKIND
    implicit none
    complex(DREALKIND) :: x
    character(59) :: complex_to_string
    complex_to_string = "(" // trim(double_to_string(real(x))) // &
                      & "," // trim(double_to_string(aimag(x))) // ")"
  end function complex_to_string


  function string_to_integer(c)
    implicit none
    character(len=*), intent(in) :: c
    integer :: string_to_integer
    integer :: stat
    read(c,*,iostat=stat) string_to_integer
    if (stat /= 0) then
!      print*, "[OpenLoops] Error: string to integer conversion not possible for string: ", c
      string_to_integer = -huge(string_to_integer)
    end if
  end function string_to_integer


  pure function factorial(k)
    implicit none
    integer, intent(in) :: k
    integer :: factorial, i
    factorial = 1
    do i = 1, k
      factorial = factorial * i
    end do
  end function factorial


  pure function binomial(n, k)
    implicit none
    integer :: binomial
    integer, intent(in) :: n, k
    integer :: i
    if (k < 0 .or. k > n) then
      binomial = 0
    else
      binomial = 1
      do i = 1, min(k, n-k)
        binomial = (binomial * (n-i+1)) / (i)
      end do
    end if
  end function binomial


  function nth_permutation(arr, n)
    ! return the n-th permutation of arr, defined by the set of canonically
    ! order permutations of [1,..,size(arr)] (same ordering as Mathematicas Permutation[]).
    implicit none
    integer, intent(in) :: arr(:), n
    integer :: nth_permutation(size(arr))
    integer :: arrcp(size(arr)), sz, nn, pos, ppos, fac
    sz = size(arr)
    arrcp = arr
    ppos = n - 1
    do nn = 1, size(arr)
      fac = factorial(sz-nn)
      pos = ppos/fac
      ppos = ppos - pos*fac
      nth_permutation(nn) = arrcp(pos+1)
      arrcp(pos+1:sz-nn) = arrcp(pos+2:sz-nn+1)
    end do
  end function nth_permutation


  function perm_pos(perm)
    ! Unique mapping of a permutation perm(1:n) -> integer in [1..n!].
    ! In a canonically ordered list of all permutations
    ! (like Mathematicas Permutations[Range[n]]),
    ! perm is at position perm_pos
    implicit none
    integer :: perm_pos
    integer, intent(in) :: perm(:)
    integer :: n, k, pos1, perm2(size(perm))
    perm_pos = 1
    perm2 = perm
    do n = size(perm)-1, 1, -1
      pos1 = perm2(1)
      do k = 1, n
        if (perm2(k+1) > pos1) then
          perm2(k) = perm2(k+1)-1
        else
          perm2(k) = perm2(k+1)
        end if
      end do
      perm_pos = perm_pos + (pos1-1)*factorial(n)
    end do
  end function perm_pos


  recursive subroutine compositions2(compos, n, maxpart)
    ! Store all unordered compositions of the integer n into integers >= 2
    ! in the allocatable array compos. Each composition compos(:,k) is ordered
    ! by constituents, larger first. The compositions array is ordered by
    ! constituents at the beginning of the composition, smaller first,
    ! independent of the composition length (i.e. not sorted by length).
    ! The optional parameter maxpart sets the maximal allowed size of the constituents.
    implicit none
    integer, allocatable, intent(out) :: compos(:,:)
    integer, intent(in) :: n
    integer, intent(in), optional :: maxpart
    integer :: maxp, j, k, next
    integer, allocatable :: thisc(:,:), lowerc(:,:)
    if (n == 0) then
      allocate(compos(0,1))
      return
    else if (n == 1) then
      allocate(thisc(1,10))
    else
      allocate(thisc(n/2,10))
    end if
    if (present(maxpart)) then
      maxp = maxpart
    else
      maxp = n
    end if
    next = 1
    do k = 2, min(n-2, maxp)
      call compositions2(lowerc, n-k, k)
      do j = 1, size(lowerc,2)
        if (size(thisc,2) == next) then
          allocate(compos(n/2,next))
          compos = thisc
          deallocate(thisc)
          allocate(thisc(n/2,10*next))
          thisc(:,1:next) = compos
          deallocate(compos)
        end if
        thisc(:,next) = 0
        thisc(1,next) = k
        thisc(2:size(lowerc,1)+1,next) = lowerc(:,j)
        next = next + 1
      end do
      deallocate(lowerc)
    end do
    next = next - 1
    if (n <= maxp) then
      next = next + 1
      thisc(:,next) = 0
      thisc(1,next) = n
    end if
    allocate(compos(size(thisc,1),next))
    compos = thisc(:,1:next)
    deallocate(thisc)
  end subroutine compositions2


  function compositions(n, k)
    ! compositions of n into k parts.
    ! compositions(1:k,j) for j=1:binomial(n+k-1,k-1)
    ! are the canonically ordered compositions.
    implicit none
    integer, intent(in) :: n, k
    integer :: compositions(k,binomial(n+k-1,k-1))
    integer :: compos(k-1,binomial(n+k-1,k-1))
    integer :: lowersz, pos, a, b, c, comp(k-1)
    lowersz = 1
    do a = 1, k-1
      pos = 0
      do b = 1, lowersz
        comp(1:a-1) = compositions(1:a-1,b)
        do c = 0, n-sum(comp(1:a-1))
          pos = pos + 1
          compos(1:a-1,pos) = comp(1:a-1)
          compos(a,pos) = c
        end do
      end do
      lowersz = pos
      compositions(1:a,1:lowersz) = compos(1:a,1:lowersz)
    end do
    if (k > 0) then
      do pos = 1, lowersz
        compositions(k,pos) = n - sum(compositions(1:k-1,pos))
      end do
    end if
  end function compositions


  function relative_deviation_dp(a, b) result (relative_deviation)
    use KIND_TYPES, only: DREALKIND
    implicit none
    real(DREALKIND), intent(in) :: a, b
    real(DREALKIND) :: relative_deviation
    if (a == b) then
      relative_deviation = 0
    else if ( a == 0 .or. b == 0) then
      relative_deviation = huge(a)
    else
      relative_deviation = max(abs(a/b-1), abs(b/a-1))
    end if
  end function relative_deviation_dp


  function relative_deviation_qp(a, b) result (relative_deviation)
    use KIND_TYPES, only: QREALKIND
    implicit none
    real(QREALKIND), intent(in) :: a, b
    real(QREALKIND) :: relative_deviation
    if (a == b) then
      relative_deviation = 0
    else if ( a == 0 .or. b == 0) then
      relative_deviation = huge(a)
    else
      relative_deviation = max(abs(a/b-1), abs(b/a-1))
    end if
  end function relative_deviation_qp


  function digit_agreement(a, b)
    use KIND_TYPES, only: DREALKIND
    implicit none
    real(DREALKIND), intent(in) :: a, b
    real(DREALKIND) :: digit_agreement
    if (a == b) then
      digit_agreement = 16
    else if (a == 0 .or. b == 0) then
      digit_agreement = 0
    else
      digit_agreement = -log10(relative_deviation(a, b))
    end if
  end function digit_agreement


  function random_string(n)
    ! Return a printable random string of n 6-bit characters
    ! from the set [0-9A-Za-z$_]
    implicit none
    integer, intent(in) :: n
    integer :: j
    integer(1) :: ran(n), i
    character(len=n) :: random_string
    open(42, file='/dev/urandom', access='stream', form='unformatted')
    read(42) ran
    close(42)
    do j = 1, n
      i = ishft(ran(j),-2) + 43 ! 43:     +     37+ 8=45
      if (i > 43) i = i + 1     ! 45:     -     46+ 2=48
      if (i > 45) i = i + 2     ! 48- 57: 0-9   58+ 7=65
      if (i > 57) i = i + 7     ! 65- 90: A-Z   91+ 6=97
      if (i > 90) i = i + 6     ! 97-122: a-z
      random_string(j:j) = char(i)
    end do
  end function random_string


  function to_lowercase(instr)
    ! return instr with uppercase letters converted to lowercase
    implicit none
    character(*), intent(in) :: instr
    character(len(instr)) :: to_lowercase
    integer :: i, n
    to_lowercase = instr
    do i = 1, len(to_lowercase)
      n = index(upper_case, to_lowercase(i:i))
      if (n /= 0) to_lowercase(i:i) = lower_case(n:n)
    end do
  end function to_lowercase


  function count_substring(s1, s2) result(c)
  ! counts the occurance of string s2 in string s1
    character(*), intent(in) :: s1, s2
    integer :: c, p, posn

    c = 0
    if(len(s2) == 0) return
    p = 1
    do
      posn = index(s1(p:), s2)
      if(posn == 0) return
      c = c + 1
      p = p + posn + len(s2)
    end do
  end function count_substring

end module ol_generic



module ol_iso_c_utilities
  ! function c_f_string_ptr(cptr)
  !   type(c_ptr), intent(in) :: cptr
  !   character(kind=c_char), pointer :: c_f_string_ptr(:)
  !   - convert a null terminated C character array pointer to a Fortran string pointer;
  ! subroutine c_f_string_static(c_str, f_str, maxlen):
  !   character(kind=c_char), dimension(*), intent(in) :: c_str
  !   integer, intent(in) :: maxlen
  !   character(len=maxlen), intent(out) :: f_str
  !   - convert a null terminated C character array to a Fortran string;
  ! subroutine c_f_string_alloc(c_str, f_str):
  !   character(kind=c_char), dimension(*), intent(in) :: c_str
  !   character(len=:), allocatable, intent(out) :: f_str
  !   - convert a null terminated C character array to a Fortran allocatable string;
  use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_long, c_short
  implicit none

  character(kind=c_char), save, target, private :: dummy_string(1) = "?"

  interface
    function strlen(string) bind(c)
      ! int strlen(char *string)
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), value :: string
      integer(c_int) :: strlen
    end function strlen
  end interface

  interface c_f_string
    module procedure c_f_string_static!, c_f_string_alloc
  end interface c_f_string

  contains

  function c_f_string_ptr(cptr)
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_associated, c_f_pointer
    implicit none
    type(c_ptr), intent(in) :: cptr
    character(kind=c_char), pointer :: c_f_string_ptr(:)
    if (c_associated(cptr)) then
      call c_f_pointer(cptr, c_f_string_ptr, shape = [strlen(cptr)])
    else
      ! if cptr is a null pointer, return 'dummy_string'
      c_f_string_ptr => dummy_string
    end if
  end function c_f_string_ptr


  subroutine c_f_string_static(c_str, f_str, maxlen)
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_loc, c_f_pointer
    implicit none
    character(kind=c_char), dimension(*), intent(in), target :: c_str
    integer, intent(in) :: maxlen
    type(c_ptr) :: c_str_ptr
    character(len=maxlen), intent(out) :: f_str
    character(kind=c_char), pointer :: f_str_ptr(:)
    integer :: slen, i
    c_str_ptr = c_loc(c_str)
    slen = strlen(c_str_ptr)
    call c_f_pointer(c_str_ptr, f_str_ptr, shape = [slen])
    f_str = ""
    do i = 1, slen
      f_str(i:i) = f_str_ptr(i)
    end do
  end subroutine c_f_string_static


!   subroutine f_c_string_static(f_str, c_str, maxlen)
!     use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_loc, c_f_pointer
!     implicit none
!     character(len=*), intent(in) :: f_str
!     character(kind=c_char), dimension(size(f_str)), intent(out), target :: c_str
!     type(c_ptr) :: c_str_ptr
!     character(len=maxlen), intent(out) :: f_str
!     character(kind=c_char), pointer :: f_str_ptr(:)
!     integer :: slen, i
!     c_str_ptr = c_loc(c_str)
!     slen = strlen(c_str_ptr)
!     call c_f_pointer(c_str_ptr, f_str_ptr, shape = [slen])
!     f_str = ""
!     do i = 1, slen
!       f_str(i:i) = f_str_ptr(i)
!     end do
!   end subroutine c_f_string_static
!

! deactivate to circumvent a bug in certain gfortran versions.
! Previously used in register_process_c, olp_setparameter_c, olp_printparameter_c, olp_start_c
!   subroutine c_f_string_alloc(c_str, f_str)
!     use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_loc, c_f_pointer
!     implicit none
!     character(kind=c_char), dimension(*), intent(in), target :: c_str
!     type(c_ptr) :: c_str_ptr
!     character(len=:), allocatable, intent(out) :: f_str
!     character(kind=c_char), pointer :: f_str_ptr(:)
!     integer :: slen, i
!     c_str_ptr = c_loc(c_str)
!     slen = strlen(c_str_ptr)
!     call c_f_pointer(c_str_ptr, f_str_ptr, shape = [slen])
!     if (allocated(f_str)) deallocate(f_str)
!     allocate(character(len=size(f_str_ptr))::f_str) ! len=slen does not work
!     do i = 1, slen
!       f_str(i:i) = f_str_ptr(i)
!     end do
!   end subroutine c_f_string_alloc


end module ol_iso_c_utilities


module ol_cwrappers
  ! err = opendir(dirname)
  !   open directory; only one directory can be open at a time;
  !   err=0 if successful
  ! err = readdir(entryname)
  !   read next entry from directory;
  !   err=0 if successful; entryname="" if all entries were retrieved;
  ! closedir()
  !   close directory
  ! mkdir(dirname)
  !   create directory
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char
  use ol_iso_c_utilities, only: c_f_string
  implicit none
  private
  public :: opendir, readdir, closedir, mkdir, direntry_length
  public :: stdout_off, stdout_on
  integer, parameter :: direntry_length = 256
  interface
    function c_opendir(dirname) bind(c,name="ol_c_opendir")
      ! int ol_c_opendir(const char *dirname)
      use, intrinsic :: iso_c_binding, only: c_char, c_int
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: dirname
      integer(c_int) :: c_opendir
    end function c_opendir
    function c_readdir(entryname) bind(c,name="ol_c_readdir")
      ! int ol_c_readdir(char* entryname)
      use, intrinsic :: iso_c_binding, only: c_char, c_int
      implicit none
      character(kind=c_char), intent(out) :: entryname(256) ! =direntry_length
      integer(c_int) :: c_readdir
    end function c_readdir
    subroutine c_closedir() bind(c,name="ol_c_closedir")
      ! void ol_c_closedir()
      implicit none
    end subroutine c_closedir
    function c_mkdir(dirname) bind(c,name="ol_c_mkdir")
      ! int ol_c_mkdir(char* dirname)
      use, intrinsic :: iso_c_binding, only: c_char, c_int
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: dirname
      integer(c_int) :: c_mkdir
    end function c_mkdir
  end interface

  interface
    subroutine stdout_off() bind(c,name="ol_c_stdout_off")
      implicit none
    end subroutine stdout_off
    subroutine stdout_on() bind(c,name="ol_c_stdout_on")
      implicit none
    end subroutine stdout_on
  end interface

  contains

  function opendir(dirname)
    implicit none
    character(len=*), intent(in) :: dirname
    integer :: opendir
    opendir = c_opendir(trim(dirname) // c_null_char)
    if (opendir == 127) then
      print *, "[OpenLoops] opendir: a directory is already open."
    else if (opendir /= 0) then
      print *, "[OpenLoops] opendir: error", opendir
    end if
  end function opendir

  function readdir(entryname)
    implicit none
    character(*), intent(out) :: entryname
    integer :: readdir
    character(kind=c_char) :: c_entryname(direntry_length)
    if (len(entryname) < 256) then
      print *, "[OpenLoops] readdir argument length <256."
      readdir = 127
      return
    end if
    readdir = c_readdir(c_entryname)
    if (readdir /= 0) then
      print *, "[OpenLoops] readdir: error reading directory content."
    end if
    entryname = ""
    call c_f_string(c_entryname, entryname, direntry_length)
  end function readdir

  subroutine closedir()
    implicit none
    call c_closedir()
  end subroutine closedir

  function mkdir(dirname)
    implicit none
    character(len=*), intent(in) :: dirname
    integer :: mkdir
    mkdir = c_mkdir(trim(dirname) // c_null_char)
  end function mkdir

end module ol_cwrappers


module ol_dlfcn
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_ptr, c_funptr, &
    & c_null_char, c_associated, c_f_procpointer
  implicit none
  private
  public :: RTLD_LAZY, RTLD_NOW, RTLD_GLOBAL, RTLD_LOCAL
  public :: dlopen, dlsym, dlclose
  ! dlopen modes:
  integer(c_int), bind(c,name="ol_c_rtld_lazy") :: RTLD_LAZY
  integer(c_int), bind(c,name="ol_c_rtld_now") :: RTLD_NOW
  integer(c_int), bind(c,name="ol_c_rtld_global") :: RTLD_GLOBAL
  integer(c_int), bind(c,name="ol_c_rtld_local") :: RTLD_LOCAL

  interface
    function c_dlopen(file, mode) bind(c,name="dlopen")
      ! void *dlopen(const char *file, int mode);
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: file
      integer(c_int), value :: mode
      type(c_ptr) :: c_dlopen
    end function c_dlopen
    function c_dlsym(lib, sym) bind(c,name="dlsym")
      ! void *dlsym(void *lib, const char *sym);
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_funptr
      implicit none
      type(c_ptr), value :: lib
      character(kind=c_char), dimension(*), intent(in) :: sym
      type(c_funptr) :: c_dlsym
    end function c_dlsym
    function c_dlclose(lib) bind(c,name="dlclose")
      ! int dlclose(void *lib);
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int
      implicit none
      type(c_ptr), value :: lib
      integer(c_int) :: c_dlclose ! status
    end function c_dlclose
    function c_dlerror() bind(c,name="dlerror")
      ! char *dlerror(void);
      use, intrinsic :: iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr) :: c_dlerror
    end function c_dlerror
  end interface

  contains

  function dlerror()
    use ol_iso_c_utilities, only: c_f_string_ptr
    implicit none
    character(kind=c_char), dimension(:), pointer :: dlerror
    dlerror => c_f_string_ptr(c_dlerror())
  end function

  function dlopen(file, mode, fatal)
    ! fatal: 0=silent (default), 1=warning, 2=error
    implicit none
    character(len=*), intent(in) :: file
    integer(c_int), intent(in) :: mode
    integer, intent(in), optional :: fatal
    type(c_ptr) :: dlopen
    dlopen = c_dlopen(trim(file) // c_null_char, mode)
    if (present(fatal)) then
      if (fatal == 1 .and. .not. c_associated(dlopen)) then
        print *, "[OpenLoops] dlopen:", dlerror()
      else if (fatal == 2 .and. .not. c_associated(dlopen)) then
        print *, "[OpenLoops] error in dlopen:", dlerror()
        stop
      end if
    end if
  end function dlopen

  function dlsym(lib, sym, fatal) result(f_funp)
    ! fatal: 0=silent (default), 1=warning, 2=error
    implicit none
    type(c_ptr), intent(in) :: lib
    character(len=*), intent(in) :: sym
    integer, intent(in), optional :: fatal
    type(c_funptr) :: c_funp
    procedure(), pointer :: f_funp
    c_funp = c_dlsym(lib, trim(sym) // c_null_char)
    if (present(fatal)) then
      if (fatal == 1 .and. .not. c_associated(c_funp)) then
        print *, "[OpenLoops] dlsym:", dlerror()
      else if (fatal == 2 .and. .not. c_associated(c_funp)) then
        print *, "[OpenLoops] error in dlsym:", dlerror()
        stop
      end if
    end if
    if (c_associated(c_funp)) then
      call c_f_procpointer(c_funp, f_funp)
    else
      f_funp => null()
    end if
  end function dlsym

  subroutine dlclose(lib, fatal)
    ! fatal: 0=silent (default), 1=warning, 2=error
    implicit none
    type(c_ptr), intent(in) :: lib
    integer, intent(in), optional :: fatal
    integer(c_int) :: status
    status = c_dlclose(lib)
    if (present(fatal)) then
      if (fatal == 1 .and. status /= 0) then
        print *, "[OpenLoops] dlclose:", dlerror()
      else if (fatal == 2 .and. status /= 0) then
        print *, "[OpenLoops] error in dlclose:", dlerror()
        stop
      end if
    end if
  end subroutine dlclose

end module ol_dlfcn

#endif



module ol_dilog_/**/REALKIND
  use KIND_TYPES, only: REALKIND
  implicit none
  real(REALKIND), parameter :: pi2_6 = 8/3._/**/REALKIND*atan(1._/**/REALKIND)**2
#ifdef PRECISION_sp
  integer, parameter :: max_n = 15
#else
  integer, parameter :: max_n = 29
#endif
#ifdef CHECK_DEPTH
  integer, save :: bin(max_n) = 0
#endif
  real(REALKIND), parameter :: G4  = 6,         G6  = G4 * 4* 5, G8  = G6 * 6* 7, G10 = G8 * 8* 9, G12 = G10*10*11
  real(REALKIND), parameter :: G14 = G12*12*13, G16 = G14*14*15, G18 = G16*16*17, G20 = G18*18*19, G22 = G20*20*21
  real(REALKIND), parameter :: G24 = G22*22*23, G26 = G24*24*25, G28 = G26*26*27, G30 = G28*28*29, G32 = G30*30*31
#ifndef PRECISION_sp
  real(REALKIND), parameter :: G34 = G32*32*33, G36 = G34*34*35, G38 = G36*36*37, G40 = G38*38*39, G42 = G40*40*41
  real(REALKIND), parameter :: G44 = G42*42*43, G46 = G44*44*45, G48 = G46*46*47, G50 = G48*48*49, G52 = G50*50*51
  real(REALKIND), parameter :: G54 = G52*52*53, G56 = G54*54*55, G58 = G56*56*57, G60 = G58*58*59, G62 = G60*60*61
#endif
  real(REALKIND), parameter :: B2n(max_n) = [ & ! BernoulliB[2*n]/Gamma[2*n+2]
    &                                     1._/**/REALKIND / (      6 * G4 ) & !  1
    & ,                                  -1._/**/REALKIND / (     30 * G6 ) & !  2
    & ,                                   1._/**/REALKIND / (     42 * G8 ) & !  3
    & ,                                  -1._/**/REALKIND / (     30 * G10) & !  4
    & ,                                   5._/**/REALKIND / (     66 * G12) & !  5
    & ,                                -691._/**/REALKIND / (   2730 * G14) & !  6
    & ,                                   7._/**/REALKIND / (      6 * G16) & !  7
    & ,                               -3617._/**/REALKIND / (    510 * G18) & !  8
    & ,                               43867._/**/REALKIND / (    798 * G20) & !  9
    & ,                             -174611._/**/REALKIND / (    330 * G22) & ! 10
    & ,                              854513._/**/REALKIND / (    138 * G24) & ! 11
    & ,                          -236364091._/**/REALKIND / (   2730 * G26) & ! 12
    & ,                             8553103._/**/REALKIND / (      6 * G28) & ! 13
    & ,                        -23749461029._/**/REALKIND / (    870 * G30) & ! 14
    & ,                       8615841276005._/**/REALKIND / (  14322 * G32) & ! 15
#ifndef PRECISION_sp
    & ,                      -7709321041217._/**/REALKIND / (    510 * G34) & ! 16
    & ,                       2577687858367._/**/REALKIND / (      6 * G36) & ! 17
    & ,               -26315271553053477373._/**/REALKIND / (1919190 * G38) & ! 18
    & ,                    2929993913841559._/**/REALKIND / (      6 * G40) & ! 19
    & ,              -261082718496449122051._/**/REALKIND / (  13530 * G42) & ! 20
    & ,              1520097643918070802691._/**/REALKIND / (   1806 * G44) & ! 21
    & ,            -27833269579301024235023._/**/REALKIND / (    690 * G46) & ! 22
    & ,            596451111593912163277961._/**/REALKIND / (    282 * G48) & ! 23
    & ,       -5609403368997817686249127547._/**/REALKIND / (  46410 * G50) & ! 24
    & ,         495057205241079648212477525._/**/REALKIND / (     66 * G52) & ! 25
    & ,     -801165718135489957347924991853._/**/REALKIND / (   1590 * G54) & ! 26
    & ,    29149963634884862421418123812691._/**/REALKIND / (    798 * G56) & ! 27
    & , -2479392929313226753685415739663229._/**/REALKIND / (    870 * G58) & ! 28
    & , 84483613348880041862046775994036021._/**/REALKIND / (    354 * G60) & ! 29
#endif
  & ]

  contains

! **********************************************************************
function Li2conv(z)
! Complex dilogarithm Li2(z) calculated from a series expansion in terms of Bernoulli numbers.
! This is supposed to be used in the complex region |z| < 1 && Re(z) < 1/2.
! In this region the expansion (max_n = 29) is deep enough to achieve quadruple precision.
! Required expansion depths: sp:9, dp:13, ep:15, qp:24
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: z
  complex(REALKIND) :: Li2conv
  complex(REALKIND) :: Lz2, Lz2n1, Li2conv_next
  integer           :: n
  Lz2n1   = -log(1-z)
  Lz2     = Lz2n1*Lz2n1
  Li2conv = Lz2n1 - 0.25_/**/REALKIND*Lz2
  do n = 1, max_n
    Lz2n1 = Lz2n1 * Lz2
    Li2conv_next = Li2conv + Lz2n1 * B2n(n)
#ifdef COMPARE
    if (Li2conv_next == Li2conv) exit
#endif
    Li2conv = Li2conv_next
  end do
#ifdef CHECK_DEPTH
  bin(n) = bin(n) + 1
#endif
end function Li2conv


! **********************************************************************
function Li2(z)
! Complex dilogarithm Li2(z)
! Map z to the complex region |z| < 1 && Re(z) < 1/2 and call the series expansion.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: z
  complex(REALKIND) :: Li2
  if (real(z) <= 0.5) then
    if (abs(z) <= 1) then
      Li2 =   Li2conv(z)
    else
      Li2 = - Li2conv(1/z)     -   pi2_6                   - .5_/**/REALKIND * log( -z)**2
    end if
  else
    if (abs(1-z) <= 1) then
      Li2 = - Li2conv(1-z)     +   pi2_6 - log(z)*log(1-z)
    else
      Li2 =   Li2conv(1/(1-z)) + 2*pi2_6 - log(z)*log(1-z) + .5_/**/REALKIND * log(z-1)**2
    end if
  end if
end function Li2

end module ol_dilog_/**/REALKIND
