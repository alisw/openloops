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


module ind_bookkeeping_/**/REALKIND

contains

!******************************************************************************
integer function CalcSind(Sarr)
!------------------------------------------------------------------------------
! determines the integer label associated with the particle subset Sarr
! A subset S of the external particles n=1,...,Npart can be
! characterised by an array Sarr(1:Npart) with components
!    Sarr(n)=1 if particle n belongs to the subset
!    Sarr(n)=0 otherwise.
! This information can be encoded in an integer label
!    Sind = Sum_{n} Sarr(n)*2^{n-1}
!******************************************************************************
 implicit none
 integer, intent(in)  :: Sarr(:)
 integer :: n_part, n,shift

 n_part=size(Sarr)

 shift = 1
 CalcSind = 0
 do n=1,n_part
   CalcSind = CalcSind+Sarr(n)*shift
   shift=shift*2
 end do

end function CalcSind


!******************************************************************************
integer function AddPart(n,Sind)
!------------------------------------------------------------------------------
!  adds particle n to subset Sind
!******************************************************************************
  implicit none
  integer, intent(in)  :: n
  integer, intent(in)  :: Sind

  AddPart = IBSET(Sind, n-1)

end function AddPart


!******************************************************************************
integer function CalcHind(Harr)
!------------------------------------------------------------------------------
! Determines the helicity label associated with the helicity configuration Harr
! A global helicity configuration for the external particles n=1,...,Npart
! can be characterised by an array Harr(1:Npart) with components
!    Harr(n)=0  if  particle n is unpolarised (helicity summed)
!    Harr(n)=1,2,3  for physical helicty states
! This information can be encoded in an integer helicity label
!    Hind = Sum_{n} Harr(n)*4^{n-1}
!******************************************************************************
 implicit none
 integer,    intent(in)  :: Harr(:)
 integer  :: n_part,n,shift

  n_part=size(Harr)

  shift = 1
  CalcHind = 0
  do n=1,n_part
    CalcHind = CalcHind+Harr(n)*shift
    shift=shift*4
  end do

end function CalcHind


!******************************************************************************
integer function SetPartHel(n,h,Hind)
!------------------------------------------------------------------------------
!  sets helicity of particle(n) to h=0,1,2,3 in helicity configuration Hind
!******************************************************************************
  implicit none
  integer, intent(in)  :: n,h
  integer, intent(in) :: Hind
  integer  :: Hind2

  Hind2=Hind

  call MVBITS(h,0,2,Hind2,2*n-2)
  SetPartHel=Hind2

end function SetPartHel


!******************************************************************************
integer(intkind2) function ProjSind(Sind,n_part)
!------------------------------------------------------------------------------
! transformation of particle subset label that can be exploited
! to project out helicity configuration of a certain particle subset
! using the IAND (see function ProjHind) In binary notation
!  ProjSind(xyz...) = xxyyzz... for z,y,w = 0,1
! n= total number of external particles
!******************************************************************************
  use KIND_TYPES, only: intkind2
  implicit none
  integer(intkind2), intent(in)  :: Sind
  integer(intkind2), intent(in)  :: n_part
  integer :: n

  ProjSind=0

  do n=0, floor(LOG(1.*Sind)/LOG(2.))+1
  call MVBITS(Sind,n,1,ProjSind,2*n)
  call MVBITS(Sind,n,1,ProjSind,2*n+1)
  end do

end function ProjSind


!******************************************************************************
integer function ProjHind(Sind,Hind,n_part)
!------------------------------------------------------------------------------
! Hind     = label of global helicity configuration
! Sind     = label of external particle subset
! ProjHind = label of reduces helicity configuration
!            including only given particle subset
! n_part   = total number of ext particles
!******************************************************************************
  use KIND_TYPES, only: intkind2
  implicit none
  integer(intkind2),   intent(in)  :: Sind, Hind
  integer(intkind2), intent(in) :: n_part

  ProjHind=IAND(Hind,ProjSind(Sind,n_part))

end function ProjHind

end module ind_bookkeeping_/**/REALKIND



module ol_helicity_init_/**/REALKIND

contains

!******************************************************************************
subroutine heltable(n_in, n, tab)
!------------------------------------------------------------------------------
! helicity table for tree vertex with K-legs
!  K=size(n)
!  1,....,K-1 / K  = incoming / outgoing legs
!  n(1:K)          = # helicity states for each leg
!  t(1:K-1,n(K))   = helicity table
!  t(i_j,i_K)      = j-th wave function input helicity
!                    associated with output helicity i_K
!******************************************************************************
  implicit none
  integer, intent(in)  :: n_in(:)
  integer, intent(out) :: n(:), tab(:,:)
  integer :: i1, i2, i3, i4, i5, i6

  n=n_in

  if (size(n) == 3) then
    i3 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        tab(1,i3) = i1
        tab(2,i3) = i2
        i3 = i3 + 1 ! <=> i3 = i2 + n(2) * (i1-1)
      end do
    end do
  else if (size(n) == 4) then
    i4 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          tab(1,i4) = i1
          tab(2,i4) = i2
          tab(3,i4) = i3
          i4 = i4 + 1 ! <=> i4 = i3 + n(3) * (i2-1 + n(2) * (i1-1))
        end do
      end do
    end do
  else if (size(n) == 5) then
    i5 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          do i4 = 1, n(4)
            tab(1,i5) = i1
            tab(2,i5) = i2
            tab(3,i5) = i3
            tab(4,i5) = i4
            i5 = i5 + 1
          end do
        end do
      end do
    end do
  else if (size(n) == 6) then
    i6 = 1
    do i1 = 1, n(1)
      do i2 = 1, n(2)
        do i3 = 1, n(3)
          do i4 = 1, n(4)
            do i5 = 1, n(5)
              tab(1,i6) = i1
              tab(2,i6) = i2
              tab(3,i6) = i3
              tab(4,i6) = i4
              tab(5,i6) = i5
              i6 = i6 + 1
            end do
          end do
        end do
      end do
    end do
  end if

end subroutine heltable

end module ol_helicity_init_/**/REALKIND


module hel_bookkeeping_/**/REALKIND
contains

!******************************************************************************
subroutine helbookkeeping_prop(ntry, N_in, N_out, n)
! helicity mappings in OpenLoops step. Propagators
! -----------------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! N_in(1:n)  = input open-loop of type hol
! N_out(1:n) = output open loop of type hol
!******************************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  use ind_bookkeeping_/**/REALKIND, only : ProjHind
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  type(hol),     intent(in)    :: N_in
  type(hol),     intent(inout)   :: N_out
  integer(intkind2) :: nhel
  integer :: h

  nhel = 0
  do h = 1, size(N_in%hf)
    if (N_in%hf(h) /= -1_intkind2) nhel = nhel + 1
  end do
  n = nhel

  N_out%hf = N_in%hf

end subroutine helbookkeeping_prop

!******************************************************************************
subroutine helbookkeeping_last_prop(ntry, N_in, n)
! helicity mappings in OpenLoops step. Last Propagator step
! -----------------------------------------------------------------------------
! ntry       = 1 (2) for 1st (subsequent) PS points
! N_in(1:n)  = input open-loop of type hol
! N_out(1:n) = output open loop of type hol
!******************************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  type(hol),     intent(in)    :: N_in
  integer(intkind2) :: nhel
  integer :: h

  nhel = 0
  do h = 1, size(N_in%hf)
    if (N_in%hf(h) /= -1_intkind2) nhel = nhel + 1
  end do
  n = nhel

end subroutine helbookkeeping_last_prop

!******************************************************************************
subroutine helbookkeeping_ol_vert3(ntry, W, N_in, N_out, n, t)
! helicity mappings in OpenLoops step - Three point vertex
! -----------------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! W(1:n(1))     = external wave function
! N_in(1:n(3))  = input OpenLoop of type hol
! N_out(1:n(2)) = output OpenLoop of type hol
! helicity table t(1:2,n(3)) links helicity states of N_in(h3)
! with W(t(1,h3)) and N_out(t(2,h3))
! n_part        = total number of ext particles
!******************************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  use ol_debug, only: ol_msg
  use ind_bookkeeping_/**/REALKIND, only : ProjHind
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3)
  integer(intkind2), intent(inout) :: t(:,:)
  type(wfun),    intent(in)    :: W(:)
  type(hol),     intent(in)    :: N_in
  type(hol),     intent(inout) :: N_out
  integer(intkind2)  :: n_expart
  integer(intkind2)  :: h1, h2, h3, k3
  integer(intkind2)  :: nhel_wf, nhel_in
  integer  :: s(3), n2_hel, k
  integer(intkind2)  :: hel3, hel2(size(N_in%hf)), hel1, subset, helaux
  logical  :: nintest1(size(W)),nintest3(size(N_in%hf))
  integer(intkind2) :: hf_min, hf_aux, i, imin, iq, taux!, taux(2)

  s(1)=size(W)
  s(3)=size(N_in%hf)
  s(2)=size(N_out%hf)

  n_expart = W(1)%n_part

  if (s(1)>n(1)) then
    call ol_msg(5, "Helicity bookkeeping, 3-pt vertex: inconsistent size of input WF")
  end if

  if (s(3)>n(3)) then
    call ol_msg(5, "Helicity bookkeeping, 3-pt vertex: inconsistent size of input OL")
  end if

  subset = W(1)%t

  !! Setting n(1) equal to the number of non-vanishing helicity configurations
  !! of the external wavefunction
  nhel_wf = 0
  do h1 = 1, s(1)
    if (W(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
  end do
  n(1) = nhel_wf

  !! Setting n(3) equal to the number of non-vanishing helicity configurations
  !! of the input open-loop
  nhel_in = 0
  do h3 = 1, s(3)
    if (N_in%hf(h3) /= -1_intkind2) nhel_in = nhel_in + 1
  end do
  n(3) = nhel_in

  nintest1 = .true.
  nintest3 = .true.

  hel2 = -1_intkind2

  !! assignment of the first entry of the helicity table:
  !! wavefunction - input openloop
  do h1=1,n(1)
    hel1 = W(h1)%hf
    do h3=1,n(3)
      if (nintest3(h3)) then
        hel3=N_in%hf(h3)
        if (ProjHind(subset,hel3,n_expart)==hel1) then
          t(1,h3)=h1
          if(nintest3(h3)) nintest3(h3)=.false.
          if(nintest1(h1)) nintest1(h1)=.false.
          hel2(h3)=hel3-hel1 !used below to determine h2=t(2,h3)
        end if
      end if
    end do
  end do


  do h1=1,n(1)
    if (nintest1(h1)) then
      call ol_msg(5,"Helicity bookkeeping, 3-pt vertex: irrelevant helicities in external WF")
    end if
    exit
  end do

  do h3=1,n(3)
    if (nintest3(h3)) then
      call ol_msg(5,"Helicity bookkeeping, 3-pt vertex: type 1 inconsistency of helicity mappings")
    end if
  end do

  !! assignment of the second entry of the helicity table:
  !! output openloop - input openloop
  h2 = 0_intkind2
  nintest3 =.true.

  do h3=1,n(3)
    if (nintest3(h3)) then
      h2=h2+1
      do k3=h3,n(3)
        if (hel2(k3)==hel2(h3)) then
          t(2,k3)=h2
          nintest3(k3)=.false.
        end if
      end do
    end if
  end do

  n(2)=h2

  !! assignment of global helicity label hf for the output openloop
  N_out%hf(:) = -1_intkind2

  do h3 = 1, n(3)
    h2=t(2,h3)
    h1=t(1,h3)
    N_out%hf(h2) = N_in%hf(h3) - W(h1)%hf ! global helicity after WF-helicity sum
  end do

  !! sorting of the components of N_out(1:n(2)) and adapting the tables accordingly
  do i = 1, n(2)-1
    hf_min = N_out%hf(i)
    imin = i
    do iq = i+1, n(2)        ! look for component with helicity < hf_min
      if (N_out%hf(iq) >= hf_min) cycle
      hf_min = N_out%hf(iq)
      imin = iq
    end do
    if (imin > i) then
       hf_aux         = N_out%hf(i)
       N_out%hf(i)    = N_out%hf(imin)
       N_out%hf(imin) = hf_aux
    end if
  end do

  do i = 1, n(2)
    do iq = 1, n(3)
      if (hel2(iq) == N_out%hf(i)) then
        t(2,iq) = i
      end if
    end do
  end do

  !! Final debugging checks
  if (s(2)<n(2)) then
    call ol_msg(5,"subroutine helbookkeeping_ol_vert3: inconsistent size of output open-loop")
  end if

  do h3=1,n(3)
    if (nintest3(h3)) then
      call ol_msg(5,"In helicity bookkeeping for 3-point vertex: type 2 inconsistency of helicity mappings")
    end if
  end do

end subroutine helbookkeeping_ol_vert3


!******************************************************************************
subroutine helbookkeeping_ol_vert4(ntry, W1, W2, N_in, N_out, n, t)
! helicity mappings in OpenLoops step - Four point vertex
! -----------------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! W1(1:n(1))    = first  external wave function
! W2(1:n(2))    = second external wave function
! N_in(1:n(4))  = input open-loop of type hol
! N_out(1:n(3)) = output open loop of type hol
! npart_i      = total number of ext particles in the subtree i
! helicity table t(1:3,n(4)) links helicity states of N_in(h4)
!                with W1(t(1,h4)), W2(t(2,h4)) and N_out(t(3,h4))
!******************************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  use ol_debug, only: ol_msg
  use ind_bookkeeping_/**/REALKIND, only : ProjHind
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4)
  integer(intkind2), intent(inout) :: t(:,:)
  type(wfun),    intent(in)    :: W1(:), W2(:)
  type(hol),     intent(in)    :: N_in
  type(hol),     intent(inout) :: N_out
  !! n_part_i = total number of ext particles in the subtree i
  integer(intkind2)  :: npart_1, npart_2, subset1, subset2
  integer(intkind2)  :: h1, h2, h3, h4, k4
  integer(intkind2)  :: nhel_wf1, nhel_wf2, nhel_in
  integer  :: s(4)
  integer(intkind2)  :: hel3(n(4)), hel1, hel2, hel4
  logical  :: nintest1(n(1)), nintest2(n(2)), nintest4(n(4))
  integer(intkind2) :: hf_min, hf_aux, i, imin, iq, taux

  s(1)=size(W1)
  s(2)=size(W2)
  s(4)=size(N_in%hf)
  s(3)=size(N_out%hf)

  if (s(1)>n(1)) then
    call ol_msg(5, "In helicity bookkeeping for 4-point vertex: inconsistent size of 1st input wavefunction")
  end if

  if (s(2)>n(2)) then
    call ol_msg(5, "In helicity bookkeeping for 4-point vertex: inconsistent size of 2nd input wavefunction")
  end if

  if (s(3)>n(3)) then
    call ol_msg(5, "In helicity bookkeeping for 4-point vertex: inconsistent size of output open-loop")
  end if

  if (s(4)>n(4)) then
    call ol_msg(5, "In helicity bookkeeping for 4-point vertex: inconsistent size of input open-loop")
  end if

  npart_1 = W1(1)%n_part
  npart_2 = W2(1)%n_part

  subset1 = W1(1)%t
  subset2 = W2(1)%t

  !! Setting n(1) equal to the number of non-vanishing helicity configurations
  !! of the first external wavefunction
  nhel_wf1 = 0
  do h1 = 1, s(1)
    if (W1(h1)%hf /= -1_intkind2) nhel_wf1 = nhel_wf1 + 1
  end do
  n(1) = nhel_wf1

  !! Setting n(2) equal to the number of non-vanishing helicity configurations
  !! of the second external wavefunction
  nhel_wf2 = 0
  do h2 = 1, s(2)
    if (W2(h2)%hf /= -1_intkind2) nhel_wf2 = nhel_wf2 + 1
  end do
  n(2) = nhel_wf2

  !! Setting n(4) equal to the number of non-vanishing helicity configurations
  !! of the input openloop
  nhel_in = 0
  do h4 = 1, s(4)
    if (N_in%hf(h4) /= -1_intkind2) nhel_in = nhel_in + 1
  end do
  n(4) = nhel_in

  nintest1 = .true.
  nintest2 = .true.
  nintest4 = .true.

  !! assignment of the first entry of the helicity table: wavefunction_1 - input openloop
  do h1=1,n(1)
    hel1 = W1(h1)%hf
    do h4=1,n(4)
      if (nintest4(h4)) then
        hel4=N_in%hf(h4)
        if (ProjHind(subset1,hel4,npart_1)==hel1) then
          t(1,h4)=h1
          if(nintest4(h4)) nintest4(h4)=.false.
          if(nintest1(h1)) nintest1(h1)=.false.
          hel3(h4) = hel4 - hel1 !used below to determine h3=t(3,h4)
        end if
      end if
    end do
  end do

  do h1=1,n(1)
    if (nintest1(h1)) then
      call ol_msg(5, "In helicity bookkeeping for 4-point vertex: irrelevant helicities in the first external WF")
    end if
    exit
  end do

  do h4=1,n(4)
    if (nintest4(h4)) then
      call ol_msg(5,"In helicity bookkeeping for 4-point vertex: type 1.1 inconsistency of helicity mappings")
    end if
  end do

  !! assignment of the second entry of the helicity table: wavefunction_2 - input openloop
  nintest4 =.true.

  do h2=1,n(2)
    hel2 = W2(h2)%hf
    do h4=1,n(4)
      if (nintest4(h4)) then
        hel4=N_in%hf(h4)
        if (ProjHind(subset2,hel4,npart_2)==hel2) then
          t(2,h4)=h2
          if(nintest4(h4)) nintest4(h4)=.false.
          if(nintest2(h2)) nintest2(h2)=.false.
          hel3(h4) = hel3(h4) - hel2 !used below to determine h3=t(3,h4)
        end if
      end if
    end do
  end do

  do h2=1,n(2)
    if (nintest2(h2)) then
      call ol_msg(5, "In helicity bookkeeping for 4-point vertex: irrelevant helicities in the second external WF")
    end if
    exit
  end do

  do h4=1,n(4)
    if (nintest4(h4)) then
      call ol_msg(5,"In helicity bookkeeping for 4-point vertex: type 1.2 inconsistency of helicity mappings")
    end if
  end do

  h3 = 0_intkind2

  nintest4 =.true.

  do h4=1,n(4)
    if (nintest4(h4)) then
        h3=h3+1
        do k4=h4,n(4)
          if (hel3(k4)==hel3(h4)) then
            t(3,k4)=h3
            nintest4(k4)=.false.
          end if
        end do
    end if
  end do

  n(3)=h3

  !! Assignment of global helicity label hf for the output openloop
  N_out%hf = -1_intkind2
  do h4 = 1, n(4)
    h3=t(3,h4)
    h2=t(2,h4)
    h1=t(1,h4)
    ! global helicity after Wfun-helicity sum
    N_out%hf(h3) = N_in%hf(h4) - W1(h1)%hf - W2(h2)%hf
  end do

  !! sorting the components of N_out(1:n(3)) and adapting the tables accordingly
  do i = 1, n(3)-1
    hf_min = N_out%hf(i)
    imin = i
    do iq = i+1, n(3)             ! look for component with helicity < hf_min
      if (N_out%hf(iq) >= hf_min) cycle
      hf_min = N_out%hf(iq)
      imin = iq
    end do
    if (imin > i) then
       hf_aux         = N_out%hf(i)
       N_out%hf(i)    = N_out%hf(imin)
       N_out%hf(imin) = hf_aux
    end if
  end do

  do i = 1, n(3)
    do iq = 1, n(4)
      if (hel3(iq) == N_out%hf(i)) then
        t(3,iq) = i
      end if
    end do
  end do

  do h4=1,n(4)
    if (nintest4(h4)) then
      call ol_msg(5,"In helicity bookkeeping for 4-point vertex: type 2 inconsistency of helicity mappings")
    end if
  end do

end subroutine helbookkeeping_ol_vert4


!************************************************************************
subroutine helbookkeeping_ol_last_vert3(ntry, W, N_in, n, t)
! helicity mappings in OpenLoops check last step - Three point vertex
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! W(1:n(1))     = external wave function
! N_in(1:n(3))  = input open-loop of type hol
! N_out(1:n(2)) = output open loop of type hol
! helicity table t(1:2,n(3)) links helicity states of N_in(h3)
! with W(t(1,h3)) and N_out(t(2,h3))
! n_part        = total number of ext particles
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  use ol_debug, only: ol_msg
  use ind_bookkeeping_/**/REALKIND, only : ProjHind
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3)
  integer(intkind2), intent(inout) :: t(:,:)
  type(wfun),  intent(in)    :: W(:)
  type(hol),   intent(in)    :: N_in
  integer(intkind2)  :: n_expart
  integer(intkind2)  :: h1, h3
  integer(intkind2)  :: nhel_wf, nhel_in
  integer  :: s(3)
  integer(intkind2)  :: hel3, hel2(size(N_in%hf)), hel1, subset
  logical  :: nintest1(size(W)),nintest3(size(N_in%hf))

  s(1)=size(W)
  s(3)=size(N_in%hf)
  s(2)=1

  n_expart = W(1)%n_part

  if (s(1)>n(1)) then
     call ol_msg(5,"Helicity bookkeeping for last-step 3-point vertex: inconsistent size of input wavefunction")
  end if

  if (s(3)>n(3)) then
     call ol_msg(5,"Helicity bookkeeping for last-step 3-point vertex: inconsistent size of input open-loop")
  end if

  if (s(1)/=s(3)) then
     call ol_msg(5,"Helicity bookkeeping for last-step 3-point vertex: invalid number of helicity states")
  end if

  subset = W(1)%t

  !! Setting n(1) equal to the number of the non-vanishing helicity configurations
  !! of the external wavefunction
  nhel_wf = 0
  do h1 = 1, s(1)
    if (W(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
  end do
  n(1) = nhel_wf

  !! Setting n(3) equal to the number of the non-vanishing helicity configurations
  !! of the input openloop
  nhel_in = 0
  do h3 = 1, s(3)
    if (N_in%hf(h3) /= -1_intkind2) nhel_in = nhel_in + 1
  end do
  n(3) = nhel_in

  nintest1 = .true.
  nintest3 = .true.

  do h1=1,n(1)
    hel1 = W(h1)%hf
    do h3=1,n(3)
      if (nintest3(h3)) then
        hel3=N_in%hf(h3)
        if (ProjHind(subset,hel3,n_expart)==hel1) then
          t(1,h3)=h1
          if(nintest3(h3)) nintest3(h3)=.false.
          if(nintest1(h1)) nintest1(h1)=.false.
          hel2(h3)=hel3-hel1 !used below to determine h2=t(2,h3)
        end if
      end if
    end do
  end do

  do h1=1,n(1)
    if (nintest1(h1)) then
      call ol_msg(5,"In helicity bookkeeping for 3-point vertex - last step: irrelevant helicities in external WF")
    end if
    exit
  end do

  do h3=1,n(3)
    if (nintest3(h3)) then
      call ol_msg(5,"In helicity bookkeeping for 3-point vertex - last step: type 1 inconsistency of helicity mappings")
    end if
  end do

  t(2,1) = 1_intkind2

  !! Final debugging checks
  !! Check of global helicity label hf for the output openloop
  do h3 = 1, n(3)
    h1 = t(1,h3)
    if( N_in%hf(h3) - W(h1)%hf /= 0) then
      call ol_msg(5,"In helicity bookkeeping for 3-point vertex - last step: non-zero helicity of the last tensor integral")
    end if
  end do

  if (s(2)<n(2)) then
     call ol_msg(5,"In helicity bookkeeping for 3-point vertex - last step: inconsistent size s(2)")
  end if

  do h3=1,n(3)
    if (nintest3(h3)) then
      call ol_msg(5,"In helicity bookkeeping for 3-point vertex - last step: type 2 inconsistency of helicity mappings")
    end if
  end do

end subroutine helbookkeeping_ol_last_vert3

!************************************************************************
subroutine helbookkeeping_ol_last_vert4(ntry, W1, W2, N_in, n, t)
! helicity mappings in the last four point OpenLoops vertex
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! W1(1:n(1))    = first  external wave function
! W2(1:n(2))    = second external wave function
! N_in(1:n(4))  = input open-loop of type hol
! helicity table t(1:3,n(4)) links helicity states of N_in(h4)
!                with W1(t(1,h4)), W2(t(2,h4)) and N_out(t(3,h4))
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, hol
  use ol_debug, only: ol_msg
  use ind_bookkeeping_/**/REALKIND, only : ProjHind
  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4)
  integer(intkind2), intent(inout) :: t(:,:)
  type(wfun),  intent(in)    :: W1(:), W2(:)
  type(hol),   intent(in)    :: N_in
  integer(intkind2)  :: n_part1, n_part2, subset1, subset2  ! n_part_i = total number of ext particles in the subtree i
  integer(intkind2)  :: h1, h2, h3, h4
  integer(intkind2)  :: nhel_wf1, nhel_wf2, nhel_in
  integer  :: s(4)
  integer(intkind2)  :: hel3(n(4)), hel1, hel2, hel4
  logical  :: nintest1(n(1)), nintest2(n(2)), nintest4(n(4))

  if (ntry/=1) then
    if (ntry <= 3) call ol_msg(5, "Helicity bookkeeping for last-step 4-point vertex called for N_PSP > 1")
  end if

  s(1)=size(W1)
  s(2)=size(W2)
  s(4)=size(N_in%hf)
  s(3)=1

  n_part1 = W1(1)%n_part
  n_part2 = W2(1)%n_part

  if (s(1)>n(1)) then
     call ol_msg(5,"Helicity bookkeeping last-step 4-pt vertex: inconsistent size of 1st input wavefunction")
  end if

  if (s(2)>n(2)) then
     call ol_msg(5,"Helicity bookkeeping last-step 4-pt vertex: inconsistent size of 2nd input wavefunction")
  end if

  if (s(3)>n(3)) then
     call ol_msg(5,"Helicity bookkeeping last-step 4-pt vertex: inconsistent size of input open-loop")
  end if

  if (s(4)>n(4)) then
     call ol_msg(5,"Helicity bookkeeping last-step 4-pt vertex: invalid number of helicity states")
  end if

  subset1 = W1(1)%t
  subset2 = W2(1)%t

  !! Setting n(1) equal to the number of non-vanishing helicity configurations
  !! of the first external wavefunction
  nhel_wf1 = 0
  do h1 = 1, s(1)
    if (W1(h1)%hf /= -1_intkind2) nhel_wf1 = nhel_wf1 + 1
  end do
  n(1) = nhel_wf1

  ! Setting n(2) equal to the number of non-vanishing helicity configurations
  ! of the second external wavefunction
  nhel_wf2 = 0
  do h2 = 1, s(2)
    if (W2(h2)%hf /= -1_intkind2) nhel_wf2 = nhel_wf2 + 1
  end do
  n(2) = nhel_wf2

  !! Setting n(4) equal to the number of non-vanishing helicity configurations
  !! of the input openloop
  nhel_in = 0
  do h4 = 1, s(4)
    if (N_in%hf(h4) /= -1_intkind2) nhel_in = nhel_in + 1
  end do
  n(4) = nhel_in

  nintest1 = .true.
  nintest2 = .true.
  nintest4 = .true.

  do h1=1,n(1)
    hel1 = W1(h1)%hf
    do h4=1,n(4)
      if (nintest4(h4)) then
        hel4=N_in%hf(h4)
        if (ProjHind(subset1,hel4,n_part1)==hel1) then
          t(1,h4)=h1
          if(nintest4(h4)) nintest4(h4)=.false.
          if(nintest1(h1)) nintest1(h1)=.false.
        end if
      end if
    end do
  end do

  do h1=1,n(1)
    if (nintest1(h1)) then
      call ol_msg(5,"Helicity bookkeeping for last-step 4-point vertex: irrelevant helicities in 1st input wavefunction")
    end if
    exit
  end do

  do h4=1,n(4)
    if (nintest4(h4)) then
      call ol_msg(5,"Helicity bookkeeping for last-step 4-point vertex: type 1.1 inconsistency of helicity mappings")
    end if
  end do

  nintest4 =.true.

  do h2=1,n(2)
    hel2 = W2(h2)%hf
    do h4=1,n(4)
      if (nintest4(h4)) then
        hel4=N_in%hf(h4)
        if (ProjHind(subset2,hel4,n_part2)==hel2) then
          t(2,h4)=h2
          if(nintest4(h4)) nintest4(h4)=.false.
          if(nintest2(h2)) nintest2(h2)=.false.
        end if
      end if
    end do
  end do

  do h2=1,n(2)
    if (nintest2(h2)) then
      call ol_msg(5,"Helicity bookkeeping for last-step 4-point vertex: irrelevant helicities in 2nd input wavefunction")
    end if
    exit
  end do

  do h4=1,n(4)
    if (nintest4(h4)) then
      call ol_msg(5,"Helicity bookkeeping for last-step 4-point vertex: type 1.2 inconsistency of helicity mappings")
    end if
  end do

  nintest4 =.true.

  t(3,1) = 1_intkind2

  !! Check of global helicity label hf for the output openloop
  do h4 = 1, n(4)
    h2=t(2,h4)
    h1=t(1,h4)
    if(N_in%hf(h4) - W1(h1)%hf - W2(h2)%hf /= 0) then
      call ol_msg(5,"Helicity bookkeeping for last-step 4-point vertex: non-zero helicity of the last tensor integral")
    end if
  end do

end subroutine helbookkeeping_ol_last_vert4


end module hel_bookkeeping_/**/REALKIND



module ol_merging_/**/REALKIND
  use ol_loop_handling_/**/REALKIND, only: hp_switch,hybrid_zero_mode, &
                                           hybrid_dp_mode,hybrid_qp_mode, &
                                           hybrid_dp_qp_mode,merge_mode
  use ol_parameters_decl_/**/REALKIND, only: hp_automerge
  implicit none
  contains

function valid_merge(Ginout, Gin)
  use ol_data_types_/**/REALKIND, only: hol
  type(hol), intent(in)    :: Gin
  type(hol), intent(inout) :: Ginout
  logical :: valid_merge

  if (Gin%mode .eq. hybrid_zero_mode) then
    valid_merge = .false.
  else
    valid_merge = .true.
  end if

end function valid_merge

#ifdef PRECISION_dp
!************************************************************************
subroutine helicity_matching_check(mergestep, merge_mism, hf0, hf1)
! ----------------------------------------------------------------------
! Check if the Helicity configurations of two OpenLoops are equal.
! Done only for the first phase space point
! ----------------------------------------------------------------------
! - hf0 helicities of the absorber open-loop
! - hf1 helicities of the absorbed open-loop
! **********************************************************************
  use KIND_TYPES, only: intkind2
  implicit none
  integer(intkind2), intent(in)    :: hf0(:), hf1(:)
  integer, intent(inout) :: mergestep
  integer(intkind2), intent(inout) :: merge_mism(:)
  integer :: i

  do i = 1, size(hf0)
    if(hf0(i) /= hf1(i)) then
      merge_mism(1) = merge_mism(1) + 1_intkind2
      merge_mism(mergestep) = merge_mism(1)
      exit
    end if
  end do

end subroutine helicity_matching_check

#endif

!************************************************************************
subroutine merge_map(ntry, nth_mism, merge_tabs, merge_hels, Gin_0, Gin_1)
! ----------------------------------------------------------------------
! Merging of two OpenLoops coefficients, mapping the helicities
! accordingly
! ----------------------------------------------------------------------
! - ntry = number of phase space point
! - Gin_0 absorber open-loop
! - Gin_1 absorbed open-loop
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: hol
#ifdef PRECISION_dp
  use ol_parameters_decl_/**/QREALKIND, only: zero
#endif
  implicit none
  integer(intkind1), intent(in)  :: ntry
  integer(intkind2), intent(in)  :: nth_mism
  integer(intkind2), intent(inout)  :: merge_tabs(:,:,:), merge_hels(:)
  type(hol), intent(inout) :: Gin_0
  type(hol), intent(in)    :: Gin_1
  complex(REALKIND) :: G_aux(size(Gin_0%j,1),size(Gin_0%j,2),size(Gin_0%j,3),size(Gin_0%j,4))
#ifdef PRECISION_dp
  complex(QREALKIND) :: G_aux_qp(size(Gin_0%j,1),size(Gin_0%j,2),size(Gin_0%j,3),size(Gin_0%j,4))
#endif
  integer :: s0,s1,i1,i2,htot,nhel_1,nhel_2,h,k,h1,h2
  integer(intkind2) :: hf_aux(size(Gin_0%hf)), mhel_loc
  integer(intkind2) ::  merge_tabs_loc(size(merge_tabs,1),2)
  logical :: eq, assgn

  s0 = size(Gin_0%j,2)
  s1 = size(Gin_1%j,2)

  !! Counting the number of non-vanishing helicities of the absorber
  nhel_1 = 0
  do i1 = 1, size(Gin_0%hf)
    if(Gin_0%hf(i1) /= - 1_intkind2) nhel_1 = nhel_1 + 1
  end do

  !! Counting the number of non-vanishing helicities of the absorbed
  nhel_2 = 0
  do i2 = 1, size(Gin_1%hf)
    if(Gin_1%hf(i2) /= - 1_intkind2) nhel_2 = nhel_2 + 1
  end do

  !! If it is the first phase space spoint evaluation the helicities of the output
  !! openloops respect the initial order of the ones of the absorber.
  !! All the other helicities are appended. This is in order to guarantee that
  !! the helicity tables for the absorber will match in future steps.
  if(ntry == 1_intkind1) then
    hf_aux = Gin_0%hf
    htot = nhel_1
    do i2 = 1, nhel_2
      eq = .FALSE.
      do i1 = 1, nhel_1
        if(Gin_0%hf(i1) == Gin_1%hf(i2)) then
          eq = .TRUE.
          exit
        end if
      end do

      if(.not. eq) then
        htot = htot + 1
        hf_aux(htot) = Gin_1%hf(i2)
      end if

    end do
  else
    htot = nhel_1
    hf_aux = Gin_0%hf
  end if

  if(ntry == 1_intkind1) then
    h = 0
    do k = 1, htot
      assgn = .false.
      do i1 = 1, nhel_1
        if(hf_aux(k) == Gin_0%hf(i1)) then
          merge_tabs_loc(k,1) = i1
          assgn = .true.
          exit
        end if
      end do
      if(.not. assgn) then
        h = h + 1
        merge_tabs_loc(k,1) = nhel_1 + h
      end if
    end do

    h = 0
    do k = 1, htot
      assgn = .false.
      do i2 = 1, nhel_2
        if(hf_aux(k) == Gin_1%hf(i2)) then
          merge_tabs_loc(k,2) = i2
          assgn = .true.
          exit
        end if
      end do
      if(.not. assgn) then
        h = h + 1
        merge_tabs_loc(k,2) = nhel_2 + h
      end if
    end do
    mhel_loc = htot
  end if

  if (Gin_0%mode .eq. hybrid_dp_mode .or. Gin_0%mode .eq. hybrid_dp_qp_mode) then
    G_aux = Gin_0%j
  end if
#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. Gin_0%mode .ge. hybrid_qp_mode) then
    G_aux_qp = Gin_0%j_qp
  end if
#endif

  if(ntry == 1) then
    merge_tabs(:,:,nth_mism) = merge_tabs_loc
    merge_hels(nth_mism) = mhel_loc
  end if

  do k = 1, merge_hels(nth_mism)
    h1 = merge_tabs(k,1,nth_mism)
    h2 = merge_tabs(k,2,nth_mism)
    if (Gin_1%mode .eq. hybrid_dp_mode .or. Gin_1%mode .eq. hybrid_dp_qp_mode) then
      if (Gin_0%mode .eq. hybrid_dp_mode .or. Gin_0%mode .eq. hybrid_dp_qp_mode) then
        G_aux(:,1:s1,:,k) = G_aux(:,1:s1,:,h1) + Gin_1%j(:,1:s1,:,h2)
      else
        G_aux(:,1:s1,:,k) = Gin_1%j(:,1:s1,:,h2)
        if (s0 .gt. s1) then
          G_aux(:,s1+1:,:,k) = 0
        end if
      end if
    end if
#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. Gin_1%mode .ge. hybrid_qp_mode) then
    if (hp_switch .eq. 1 .and. Gin_0%mode .ge. hybrid_qp_mode) then
      G_aux_qp(:,1:s1,:,k) = G_aux_qp(:,1:s1,:,h1) + Gin_1%j_qp(:,1:s1,:,h2)
    else
      G_aux_qp(:,1:s1,:,k) = Gin_1%j_qp(:,1:s1,:,h2)
      if (s0 .gt. s1) then
        G_aux_qp(:,s1+1:,:,k) = 0
      end if
    end if
  end if
#endif
  end do

  if(ntry == 1) Gin_0%hf = hf_aux
  ! only compute the the result if at least one dp channel is used
  if ((Gin_0%mode .eq. hybrid_dp_mode .or. Gin_0%mode .eq. hybrid_dp_qp_mode) .or. &
      (Gin_1%mode .eq. hybrid_dp_mode .or. Gin_1%mode .eq. hybrid_dp_qp_mode))then
    Gin_0%j  = G_aux
  end if
#ifdef PRECISION_dp
  ! only compute the the result if at least one qp channel is used
  if (hp_switch .eq. 1 .and. ((Gin_0%mode .ge. hybrid_qp_mode) .or. (Gin_1%mode .ge. hybrid_qp_mode))) then
    Gin_0%j_qp  = G_aux_qp
  end if
#endif

end subroutine merge_map

!************************************************************************
subroutine ol_merge_2(ntry, mergestep, merge_mism, merge_tables, merge_hels, Gin_0, Gin_1)
! ----------------------------------------------------------------------
! Merging of two OpenLoops coefficients.
! In the first phase space point it is checked that the helicity
! states of both OL coefficients match and be safefly summed over.
! If the helicity configurations do not match the merging
! is performed accordingly checking the matching helicity states.
! ----------------------------------------------------------------------
! - ntry = number of phase space point
! - Gin_0 absorber open-loop
! - Gin_1 absorbed open-loop
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: hol
#ifdef PRECISION_qp
  use ol_merging_/**/DREALKIND, only: helicity_matching_check
#else
  use ol_loop_handling_/**/REALKIND, only: upgrade_qp,req_qp_cmp,   &
                                           hp_switch,hp_alloc_mode, &
                                           hol_alloc_hybrid
#endif
  implicit none
  integer(intkind1), intent(in)  :: ntry
  integer,           intent(inout) :: mergestep
  integer(intkind2), intent(inout) :: merge_mism(:), merge_tables(:,:,:), merge_hels(:)
  type(hol), intent(inout) :: Gin_0
  type(hol), intent(in)    :: Gin_1
  integer :: s0, s1, i
  logical :: HelMatch_loc

  s0 = size(Gin_0%j,2)
  s1 = size(Gin_1%j,2)

  if(ntry == 1) then
    call helicity_matching_check(mergestep, merge_mism, Gin_0%hf, Gin_1%hf)
  end if

  if (ntry .ne. 1_intkind1 .and. (.not. valid_merge(Gin_0, Gin_1))) return

  if (merge_mism(1) == 0_intkind2) then ! no helicity mismatch
    HelMatch_loc =  .true.
  else
    HelMatch_loc = (merge_mism(mergestep) == 0_intkind2)
  end if

  Gin_0%error = max( abs(Gin_0%error), abs(Gin_1%error) )

  if (Gin_1%mode .eq. hybrid_zero_mode .and. Gin_0%mode .eq. hybrid_zero_mode) then
    Gin_0%ndrs = 0
    Gin_0%nred = 0
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = 0
    Gin_0%nred_qp = 0
#endif
  else if (Gin_1%mode .ne. hybrid_zero_mode .and. Gin_0%mode .eq. hybrid_zero_mode) then
    Gin_0%ndrs = Gin_1%ndrs
    Gin_0%nred = Gin_1%nred
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = Gin_1%ndrs_qp
    Gin_0%nred_qp = Gin_1%nred_qp
#endif
  else if (Gin_1%mode .ne. hybrid_zero_mode .and. Gin_0%mode .ne. hybrid_zero_mode) then
    Gin_0%ndrs = Gin_0%ndrs + Gin_1%ndrs
    Gin_0%nred = Gin_0%nred + Gin_1%nred
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = Gin_0%ndrs_qp + Gin_1%ndrs_qp
    Gin_0%nred_qp = Gin_0%nred_qp + Gin_1%nred_qp
#endif
  end if

  if(HelMatch_loc) then
#ifdef PRECISION_dp
    if (hp_switch .eq. 1 .and. Gin_1%mode .ge. hybrid_qp_mode) then
      if (Gin_0%mode .ge. hybrid_qp_mode) then
        Gin_0%j_qp(:,1:s1,:,:) = Gin_0%j_qp(:,1:s1,:,:) + Gin_1%j_qp(:,1:s1,:,:)
      else
        if (hp_alloc_mode .gt. 1) call hol_alloc_hybrid(Gin_0)
        Gin_0%j_qp(:,1:s1,:,:) = Gin_1%j_qp(:,1:s1,:,:)
        if (s0 .gt. s1) then
          Gin_0%j_qp(:,s1+1:,:,:) = 0
        end if
      end if
    end if
#endif
    if (Gin_1%mode .eq. hybrid_dp_mode .or. Gin_1%mode .eq. hybrid_dp_qp_mode) then
      if (Gin_0%mode .eq. hybrid_dp_mode .or. Gin_0%mode .eq. hybrid_dp_qp_mode) then
        Gin_0%j(:,1:s1,:,:) = Gin_0%j(:,1:s1,:,:) + Gin_1%j(:,1:s1,:,:)
      else
        Gin_0%j(:,1:s1,:,:) = Gin_1%j(:,1:s1,:,:)
        if (s0 .gt. s1) then
          Gin_0%j(:,s1+1:,:,:) = 0
        end if
      end if
    end if
    Gin_0%mode = merge_mode(Gin_0%mode,Gin_1%mode)
  else
    call merge_map(ntry,merge_mism(mergestep),merge_tables,merge_hels,Gin_0,Gin_1)
    Gin_0%mode = merge_mode(Gin_0%mode,Gin_1%mode)
  end if

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. hp_automerge) then
    if ( (req_qp_cmp(Gin_0)) .and. (Gin_0%mode .ne. hybrid_qp_mode) ) then
      call upgrade_qp(Gin_0)
    end if
  end if
#endif

end subroutine ol_merge_2

!************************************************************************
subroutine ol_merge_last(Gin_0, Gin_1)
! ----------------------------------------------------------------------
! Merging of two closed OpenLoops coefficients.
! ----------------------------------------------------------------------
! - Gin_0 absorber open-loop
! - Gin_1 absorbed open-loop
! **********************************************************************
  use KIND_TYPES, only: intkind1
  use ol_data_types_/**/REALKIND, only: hol
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: upgrade_qp,req_qp_cmp,   &
                                           hp_switch,hp_alloc_mode, &
                                           hol_alloc_hybrid
#endif
  implicit none
  type(hol), intent(inout) :: Gin_0
  type(hol), intent(in)    :: Gin_1
  integer :: s0, s1

  s0 = size(Gin_0%j,2)
  s1 = size(Gin_1%j,2)

  Gin_0%error = max(abs(Gin_0%error), abs(Gin_1%error))

  if (Gin_1%mode .eq. hybrid_zero_mode .and. Gin_0%mode .eq. hybrid_zero_mode) then
    Gin_0%ndrs = 0
    Gin_0%nred = 0
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = 0
    Gin_0%nred_qp = 0
#endif
  else if (Gin_1%mode .ne. hybrid_zero_mode .and. Gin_0%mode .eq. hybrid_zero_mode) then
    Gin_0%ndrs = Gin_1%ndrs
    Gin_0%nred = Gin_1%nred
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = Gin_1%ndrs_qp
    Gin_0%nred_qp = Gin_1%nred_qp
#endif
  else if (Gin_1%mode .ne. hybrid_zero_mode .and. Gin_0%mode .ne. hybrid_zero_mode) then
    Gin_0%ndrs = Gin_0%ndrs + Gin_1%ndrs
    Gin_0%nred = Gin_0%nred + Gin_1%nred
#ifdef PRECISION_dp
    Gin_0%ndrs_qp = Gin_0%ndrs_qp + Gin_1%ndrs_qp
    Gin_0%nred_qp = Gin_0%nred_qp + Gin_1%nred_qp
#endif
  end if

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. Gin_1%mode .ge. hybrid_qp_mode) then
    if (Gin_0%mode .ge. hybrid_qp_mode) then
      Gin_0%j_qp(:,1:s1,:,:) = Gin_0%j_qp(:,1:s1,:,:) + Gin_1%j_qp(:,1:s1,:,:)
    else
      if (hp_alloc_mode .gt. 1) call hol_alloc_hybrid(Gin_0)
      Gin_0%j_qp(:,1:s1,:,:) = Gin_1%j_qp(:,1:s1,:,:)
      if (s0 .gt. s1) then
        Gin_0%j_qp(:,s1+1:,:,:) = 0
      end if
    end if
  end if
#endif
  if (Gin_1%mode .eq. hybrid_dp_mode .or. Gin_1%mode .eq. hybrid_dp_qp_mode) then
    if (Gin_0%mode .eq. hybrid_dp_mode .or. Gin_0%mode .eq. hybrid_dp_qp_mode) then
      Gin_0%j(:,1:s1,:,:) = Gin_0%j(:,1:s1,:,:) + Gin_1%j(:,1:s1,:,:)
    else
      Gin_0%j(:,1:s1,:,:) = Gin_1%j(:,1:s1,:,:)
      if (s0 .gt. s1) then
        Gin_0%j(:,s1+1:,:,:) = 0
      end if
    end if
  end if
  Gin_0%mode = merge_mode(Gin_0%mode,Gin_1%mode)

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. hp_automerge) then
    if ( (req_qp_cmp(Gin_0)) .and. (Gin_0%mode .ne. hybrid_qp_mode) ) then
      call upgrade_qp(Gin_0)
    end if
  end if
#endif

end subroutine ol_merge_last

!************************************************************************
subroutine ol_merge(ntry, mergestep, merge_mism, merge_tables, merge_hels, Gout, Gins)
! Merging of an arbitrary number of OpenLoops coefficients.
! ----------------------------------------------------------------------
! - ntry = number of phase space point
! - Gout   absorber open-loop
! - Gins_i absorbed open-loops, i > 0
! **********************************************************************
  use KIND_TYPES, only: intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: hol
  implicit none
  integer(intkind1),       intent(in)    :: ntry
  integer,                 intent(inout) :: mergestep
  integer(intkind2),       intent(inout) :: merge_mism(:), merge_tables(:,:,:), merge_hels(:)
  type(hol),               intent(inout) :: Gout
  type(hol), dimension(:), intent(in)    :: Gins

  integer :: i
  do i = 1, size(Gins)
    mergestep = mergestep + 1
    call ol_merge_2(ntry, mergestep, merge_mism, merge_tables, merge_hels, Gout, Gins(i))
  end do

end subroutine ol_merge

!************************************************************************
subroutine ol_merge_tensor2(tout,tin)
! Merging of two OpenLoops tensors.
! ----------------------------------------------------------------------
! - tout absorber tensor
! - tin incoming tensor
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_data_types_/**/REALKIND, only: hcl
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: upgrade_qp,req_qp_cmp,   &
                                           hp_switch,hp_alloc_mode, &
                                           hcl_alloc_hybrid, &
                                           hybrid_zero_mode
#endif

  type(hcl), intent(in)    :: tin
  type(hcl), intent(inout) :: tout
  integer :: rein,reout

  rein = size(tin%cmp)
  reout = size(tout%cmp)
  tout%error = max(tout%error, tin%error)

  if (tin%mode .eq. hybrid_zero_mode .and. tout%mode .eq. hybrid_zero_mode) then
    tout%ndrs = 0
    tout%nred = 0
#ifdef PRECISION_dp
    tout%ndrs_qp = 0
    tout%nred_qp = 0
#endif
  else if (tin%mode .ne. hybrid_zero_mode .and. tout%mode .eq. hybrid_zero_mode) then
    tout%ndrs = tin%ndrs
    tout%nred = tin%nred
#ifdef PRECISION_dp
    tout%ndrs_qp = tin%ndrs_qp
    tout%nred_qp = tin%nred_qp
#endif
  else if (tin%mode .ne. hybrid_zero_mode .and. tout%mode .ne. hybrid_zero_mode) then
    tout%ndrs = tout%ndrs + tin%ndrs
    tout%nred = tout%nred + tin%nred
#ifdef PRECISION_dp
    tout%ndrs_qp = tout%ndrs_qp + tin%ndrs_qp
    tout%nred_qp = tout%nred_qp + tin%nred_qp
#endif
  end if

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. tin%mode .ge. hybrid_qp_mode) then
    if (tout%mode .ge. hybrid_qp_mode) then
      tout%cmp_qp(1:rein) = tout%cmp_qp(1:rein) + tin%cmp_qp(1:rein)
    else
      if (hp_alloc_mode .gt. 1) call hcl_alloc_hybrid(tout)
      tout%cmp_qp(1:rein) = tin%cmp_qp(1:rein)
      if (reout .gt. rein) then
        tout%cmp_qp(rein+1:) = 0
      end if
    end if
  end if
#endif

  if (tin%mode .eq. hybrid_dp_mode .or. tin%mode .eq. hybrid_dp_qp_mode) then
    if (tout%mode .eq. hybrid_dp_mode .or. tout%mode .eq. hybrid_dp_qp_mode) then
      tout%cmp(1:rein) = tout%cmp(1:rein) + tin%cmp(1:rein)
    else
      tout%cmp(1:rein) = tin%cmp(1:rein)
      if (reout .gt. rein) then
        tout%cmp(rein+1:) = 0
      end if
    end if
  end if
  tout%mode = merge_mode(tout%mode,tin%mode)

#ifdef PRECISION_dp
  if (hp_switch .eq. 1 .and. hp_automerge) then
    if ( (req_qp_cmp(tout)).AND.(tout%mode.ne.hybrid_qp_mode) ) then
      call upgrade_qp(tout)
    end if
  end if
#endif

end subroutine ol_merge_tensor2

!************************************************************************
subroutine ol_merge_tensors(tout,tin)
! Merging of an abirtrary number of OpenLoops tensors.
! ----------------------------------------------------------------------
! - tout absorber tensor
! - tin array of incoming tensors
! **********************************************************************

  use ol_data_types_/**/REALKIND, only: hcl
  type(hcl),               intent(inout) :: tout
  type(hcl), dimension(:), intent(in)    :: tin
  integer :: i

  do i = 1, size(tin)
    call ol_merge_tensor2(tout, tin(i))
  end do
end subroutine ol_merge_tensors

end module ol_merging_/**/REALKIND
