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


module ol_tensor_bookkeeping
  ! A tensor is stored in a 1-dimensional array for rank zero to maximal rank.
  ! rank_to_size(r) = binimial(r+3,3): number of components for rank r.
  ! tensor_size(r) = sum_rr=0..r rank_to_size(rr) = binimial(r+4,4):
  !   number of components for rank zero up to r.
  ! hr(inc,from): tensor position when the rank is increased by adding
  !   a Lorenz index with value inc, starting from component 'from'.
  ! qproductrules: rules to calculate the direct product of a momentum q
  !   (contravariant light cone) with itself up to a given rank:
  !   tensor(i) = qproductrules(1,i)*q(2,i)
  ! l2lc: see below
  use ol_generic, only: binomial, compositions
  use ol_data_types_/**/DREALKIND, only: l2lc_data
  implicit none
  integer, save :: initialised_rank = -1
  integer, allocatable, save :: rank_to_size(:), tensor_size(:), hr(:,:)
  integer, allocatable, save :: qproductrules(:,:)
  type(l2lc_data), allocatable :: l2lc(:)
  contains

  function tensor_pos(a, b, c, d)
    implicit none
    integer :: tensor_pos
    integer, intent(in) :: a, b, c, d
    tensor_pos = (24*(1+d) &
              & + 12*(c+d)*(1+c+d) &
              & + 4*(b+c+d)*(1+b+c+d)*(2+b+c+d) &
              & + (a+b+c+d)*(1+a+b+c+d)*(2+a+b+c+d)*(3+a+b+c+d))/24
  end function tensor_pos

  subroutine init_tensorbookkeeping(maxrank)
    use KIND_TYPES, only: DREALKIND
    use ol_data_types_/**/DREALKIND, only: carray2
    use ol_parameters_decl_/**/DREALKIND, only: CI
    implicit none
    integer, intent(in) :: maxrank
    integer :: pos, r, i, j, compos(4,binomial(maxrank+3,3)), compo(4), lri(4), lrj(4)
    integer :: pos_hr(4), lr(4,binomial(maxrank+4,4)), symm(binomial(maxrank+4,4))
    type(carray2) :: lor2lc_marr(0:maxrank)
    complex(DREALKIND) :: tmp
    integer :: rdat, cdat, len_rdata, len_cdata
    integer :: rdata(2,binomial(maxrank+3,3)), cdata(2,binomial(maxrank+3,3))

    if (allocated(rank_to_size)) then
      deallocate(rank_to_size)
      deallocate(tensor_size)
      deallocate(hr)
      deallocate(qproductrules)
    end if
    allocate(rank_to_size(0:maxrank))
    rank_to_size = [(binomial(r+3,3), r=0, maxrank)]
    allocate(tensor_size(-1:maxrank))
    tensor_size = [(binomial(r+4,4), r=-1, maxrank)]
    allocate(hr(4,tensor_size(maxrank-1)))
    allocate(qproductrules(2,tensor_size(maxrank)))
    initialised_rank = maxrank
    pos = 0
    lr = 0
    do r = 0, maxrank - 1
      compos(:,1:rank_to_size(r)) = compositions(r,4)
      do i = rank_to_size(r), 1, -1
        pos = pos + 1
        compo = compos(:,i)
        symm(pos) = count(compo /= 0)
        pos_hr(1) = tensor_pos(compo(1)+1,compo(2),compo(3),compo(4))
        pos_hr(2) = tensor_pos(compo(1),compo(2)+1,compo(3),compo(4))
        pos_hr(3) = tensor_pos(compo(1),compo(2),compo(3)+1,compo(4))
        pos_hr(4) = tensor_pos(compo(1),compo(2),compo(3),compo(4)+1)
        hr(:,pos) = pos_hr
        lr(1,pos_hr(1)) = pos
        lr(2,pos_hr(2)) = pos
        lr(3,pos_hr(3)) = pos
        lr(4,pos_hr(4)) = pos
      end do
    end do
    compos = compositions(maxrank,4)
    do i = rank_to_size(maxrank), 1, -1
      pos = pos + 1
      symm(pos) = count(compos(:,i) /= 0)
    end do
    qproductrules = 0
    do r = 1, tensor_size(maxrank-1)
      do i = 1, 4
        if (qproductrules(1,hr(i,r)) == 0) qproductrules(:,hr(i,r)) = [i,r]
      end do
    end do

    ! Conversion rules from a contravariant Lorentz tensor to a contravariant light-cone tensor.
    ! lightcone(rank_to_size(r-1)+i) = sum(lor2lc_marr(r)%arr(:,i) * lorentz(rank_to_size(r-1)+1:rank_to_size(r)))
    ! vector conversion: lc^mu = M^mu_nu lor^nu with M={{1,0,0,-1},{1,0,0,1},{0,-1,-i,0},{0,-1,i,0}}
    ! tensor conversion:
    ! lc^mu1..mun = M^mu1_nu1 * ... * M^mun_nun * lor^nu1..nun = M^mu1..mun_nu1..nun lor^nu1..nun
    ! M^mu1..mun_nu1..nun = M^mu1..mun-1_nu1..nun-1 * M^mun_nun
    ! with symmetrised tensor components i <- (mu1,..,mun), j <- (nu1,..,nun)
    ! and lr(i,a) gives the tensor position with one Lorentz index with value a less starting from component j;
    ! if lr(i,a) does not exist (because i has no indices with value a), omit the term in the sum
    ! M(i,j) = sum_a,b=1..4 M^a_b * M(lr(i,a),lr(j,b)) / symm
    ! with symm = number of non-zero index values
    do r = 0, maxrank
      allocate(lor2lc_marr(r)%arr(rank_to_size(r),rank_to_size(r)))
    end do
    lor2lc_marr(0)%arr(1,1) = 1
    do r = 1, maxrank
      do i = 1, rank_to_size(r)
        do j = 1, rank_to_size(r)
          lri(1) = lr(1,i+tensor_size(r-1)) - tensor_size(r-2)
          lri(2) = lr(2,i+tensor_size(r-1)) - tensor_size(r-2)
          lri(3) = lr(3,i+tensor_size(r-1)) - tensor_size(r-2)
          lri(4) = lr(4,i+tensor_size(r-1)) - tensor_size(r-2)
          lrj(1) = lr(1,j+tensor_size(r-1)) - tensor_size(r-2)
          lrj(2) = lr(2,j+tensor_size(r-1)) - tensor_size(r-2)
          lrj(3) = lr(3,j+tensor_size(r-1)) - tensor_size(r-2)
          lrj(4) = lr(4,j+tensor_size(r-1)) - tensor_size(r-2)
          tmp = 0
          if (lri(1) > 0 .and. lrj(1) > 0) tmp = tmp + lor2lc_marr(r-1)%arr(lrj(1),lri(1))
          if (lri(1) > 0 .and. lrj(4) > 0) tmp = tmp - lor2lc_marr(r-1)%arr(lrj(4),lri(1))
          if (lri(2) > 0 .and. lrj(1) > 0) tmp = tmp + lor2lc_marr(r-1)%arr(lrj(1),lri(2))
          if (lri(2) > 0 .and. lrj(4) > 0) tmp = tmp + lor2lc_marr(r-1)%arr(lrj(4),lri(2))
          if (lri(3) > 0 .and. lrj(2) > 0) tmp = tmp - lor2lc_marr(r-1)%arr(lrj(2),lri(3))
          if (lri(3) > 0 .and. lrj(3) > 0) tmp = tmp - CI*lor2lc_marr(r-1)%arr(lrj(3),lri(3))
          if (lri(4) > 0 .and. lrj(2) > 0) tmp = tmp - lor2lc_marr(r-1)%arr(lrj(2),lri(4))
          if (lri(4) > 0 .and. lrj(3) > 0) tmp = tmp + CI*lor2lc_marr(r-1)%arr(lrj(3),lri(4))
          lor2lc_marr(r)%arr(j,i) = tmp/symm(i+tensor_size(r-1))
        end do
      end do
    end do

    ! convert lor2lc_marr to a representation for minimal effort Lorentz->lightcone conversion
    ! lightcone(i+tensor_size(rank-1)) = sum for j, a in l2lc(rank)%arr(i)%r: a*lorentz(j)
    !                            + CI * sum for j, a in l2lc(rank)%arr(i)%c: a*lorentz(j)
    if (allocated(l2lc)) then
      do r = 1, size(l2lc)
        do i = 1, size(l2lc(r)%arr)
          deallocate(l2lc(r)%arr(i)%r)
          deallocate(l2lc(r)%arr(i)%c)
        end do
        deallocate(l2lc(r)%arr)
      end do
      deallocate(l2lc)
    end if
    allocate(l2lc(maxrank))

    do r = 1, maxrank
      allocate(l2lc(r)%arr(rank_to_size(r)))
      do i = 1, rank_to_size(r)
        len_rdata = 0
        len_cdata = 0
        do j = 1, rank_to_size(r)
          rdat = int(real(lor2lc_marr(r)%arr(j,i)))
          cdat = int(aimag(lor2lc_marr(r)%arr(j,i)))
          if (rdat /= 0) then
            len_rdata = len_rdata + 1
            rdata(1,len_rdata) = j + tensor_size(r-1)
            rdata(2,len_rdata) = rdat
          end if
          if (cdat /= 0) then
            len_cdata = len_cdata + 1
            cdata(1,len_cdata) = j + tensor_size(r-1)
            cdata(2,len_cdata) = cdat
          end if
        end do
        allocate(l2lc(r)%arr(i)%r(2,len_rdata))
        allocate(l2lc(r)%arr(i)%c(2,len_cdata))
        l2lc(r)%arr(i)%r(:,:) = rdata(:,1:len_rdata)
        l2lc(r)%arr(i)%c(:,:) = cdata(:,1:len_cdata)
      end do
    end do

    do r = 0, maxrank
      deallocate(lor2lc_marr(r)%arr)
    end do

  end subroutine init_tensorbookkeeping

end module ol_tensor_bookkeeping
