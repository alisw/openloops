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


#ifdef USE_QCDLOOP
#ifdef PRECISION_dp
module ol_qcdloop_interface
  implicit none
  interface ol_qcdloop
    function qli1(m1, mu2, ep)
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli1
      real(DREALKIND), intent(in) :: m1, mu2
      integer, intent(in) :: ep
    end function qli1
    function qli1c(m1, mu2, ep)
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli1c
      complex(DREALKIND), intent(in) :: m1, mu2
      integer, intent(in) :: ep
    end function qli1c
    function qli2(p1, m1, m2, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli2
      real(DREALKIND), intent(in) :: p1
      real(DREALKIND), intent(in) :: m1, m2
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli2
    function qli2c(p1, m1, m2, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli2c
      real(DREALKIND), intent(in) :: p1
      complex(DREALKIND), intent(in) :: m1, m2
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli2c
    function qli3(p1, p2, p3, m1, m2, m3, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli3
      real(DREALKIND), intent(in) :: p1, p2, p3
      real(DREALKIND), intent(in) :: m1, m2, m3
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli3
    function qli3c(p1, p2, p3, m1, m2, m3, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli3c
      real(DREALKIND), intent(in) :: p1, p2, p3
      complex(DREALKIND), intent(in) :: m1, m2, m3
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli3c
    function qli4(p1, p2, p3, p4, s12, s23, m1, m2, m3, m4, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli4
      real(DREALKIND), intent(in) :: p1, p2, p3, p4, s12, s23
      real(DREALKIND), intent(in) :: m1, m2, m3, m4
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli4
    function qli4c(p1, p2, p3, p4, s12, s23, m1, m2, m3, m4, mu2, ep);
      use KIND_TYPES, only: DREALKIND
      implicit none
      complex(DREALKIND) :: qli4c
      real(DREALKIND), intent(in) :: p1, p2, p3, p4, s12, s23
      complex(DREALKIND), intent(in) :: m1, m2, m3, m4
      real(DREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli4c
    function qli1q(m1, mu2, ep)
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli1q
      real(QREALKIND), intent(in) :: m1, mu2
      integer, intent(in) :: ep
    end function qli1q
    function qli1qc(m1, mu2, ep)
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli1qc
      complex(QREALKIND), intent(in) :: m1, mu2
      integer, intent(in) :: ep
    end function qli1qc
    function qli2q(p1, m1, m2, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli2q
      real(QREALKIND), intent(in) :: p1
      real(QREALKIND), intent(in) :: m1, m2
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli2q
    function qli2qc(p1, m1, m2, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli2qc
      real(QREALKIND), intent(in) :: p1
      complex(QREALKIND), intent(in) :: m1, m2
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli2qc
    function qli3q(p1, p2, p3, m1, m2, m3, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli3q
      real(QREALKIND), intent(in) :: p1, p2, p3
      real(QREALKIND), intent(in) :: m1, m2, m3
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli3q
    function qli3qc(p1, p2, p3, m1, m2, m3, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli3qc
      real(QREALKIND), intent(in) :: p1, p2, p3
      complex(QREALKIND), intent(in) :: m1, m2, m3
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli3qc
    function qli4q(p1, p2, p3, p4, s12, s23, m1, m2, m3, m4, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli4q
      real(QREALKIND), intent(in) :: p1, p2, p3, p4, s12, s23
      real(QREALKIND), intent(in) :: m1, m2, m3, m4
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli4q
    function qli4qc(p1, p2, p3, p4, s12, s23, m1, m2, m3, m4, mu2, ep);
      use KIND_TYPES, only: QREALKIND
      implicit none
      complex(QREALKIND) :: qli4qc
      real(QREALKIND), intent(in) :: p1, p2, p3, p4, s12, s23
      complex(QREALKIND), intent(in) :: m1, m2, m3, m4
      real(QREALKIND), intent(in) :: mu2
      integer, intent(in) :: ep
    end function qli4qc
  end interface ol_qcdloop
end module ol_qcdloop_interface
! #ifdef PRECISION_dp
#endif
! #ifdef USE_QCDLOOP
#endif


module ol_loop_routines_/**/REALKIND
  use ol_debug, only: ol_fatal, ol_msg, ol_error
  implicit none
  contains

! ****************************************************
subroutine tensor_integral(rank, momenta, masses_2, TI)
! ****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_error, ol_msg
  use ol_generic, only: to_string
#ifdef USE_COLLIER
  use ol_parameters_decl_/**/DREALKIND, only: current_processname, rZERO
  use ol_loop_parameters_decl_/**/DREALKIND, only: tensor_reduction_error
  use ol_external_decl_/**/REALKIND, only: nParticles, P_ex, crossing, inverse_crossing
  use ol_tensor_bookkeeping, only: rank_to_size, tensor_size
  use ol_Std2LC_converter_/**/REALKIND, only: lorentz2lc_tensor
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx, momenta_invariants
  use collier, only: tnten_cll
  use collier_init, only: GetErrFlag_cll
#endif
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses_2(:)
  complex(REALKIND), intent(out) :: TI(:)
#ifdef USE_COLLIER
  complex(REALKIND), allocatable :: T2dim(:,:)
  complex(REALKIND) :: momenta_TI(0:3,size(momenta,2)), T_Lor(size(TI))
  integer           :: l, r
  integer, external :: coli_get_error_code

  do l = 1, size(momenta,2)
    call LC2Std_Rep_cmplx(momenta(:,l), momenta_TI(:,l))
  end do

#ifdef PRECISION_dp
  call tnten_cll(T_Lor, TI, momenta_TI, momenta_invariants(momenta_TI), masses_2, size(masses_2), rank)
  call GetErrFlag_cll(tensor_reduction_error) ! Get overall error code from COLLIER
  if (tensor_reduction_error < -9) then
    call ol_error("Tensor Integral Reduction with COLLIER yields error code: " // to_string(tensor_reduction_error))
    call ol_msg(1,"process: " // current_processname )
    call ol_msg(1,"crossing:" // to_string(crossing(1:nParticles)))
    call ol_msg(1,"phase space point:")
    do l = 1, nParticles
      print*, P_ex(:,l)
      crossing(inverse_crossing(l)) = l
    end do
    T_Lor = 0
    T_Lor = 1./T_Lor
  end if
#else
  T_Lor = 0
  call ol_fatal('in tensor_integral: Collier is not available in quad precision')
#endif
  call lorentz2lc_tensor(rank, T_Lor, TI)
#else
  TI = 0 ! prevent compiler warning
  call ol_fatal('in tensor_integral: Collier is not available')
#endif
end subroutine tensor_integral



! ****************************************************
subroutine tensor_integral_contract(rank, momenta, masses_2, Gtensor, M2)
! ****************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses_2(:), Gtensor(:)
  complex(REALKIND), intent(out) :: M2
  complex(REALKIND) :: TI(size(Gtensor))
  call tensor_integral(rank, momenta, masses_2, TI)
  M2 = tensor_contract(Gtensor, TI)
end subroutine tensor_integral_contract



! ****************************************************
subroutine scalar_integral(momenta, masses_2)
! ****************************************************
  use KIND_TYPES, only: REALKIND
#ifdef USE_COLLIER
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep, momenta_invariants
  use collier, only: tnten_cll
#endif
  implicit none
  complex(REALKIND), intent(in)  :: momenta(:,:), masses_2(:)
#ifdef USE_COLLIER
  complex(REALKIND) :: T2dim(1,0:0), tuv(0:0)
  complex(REALKIND) :: momenta_TI(0:3,size(momenta,2))
  real(REALKIND)    :: mom(0:3)
  integer           :: l
  do l = 1, size(momenta,2)
    call LC2Std_Rep(momenta(:,l), mom)
    momenta_TI(:,l) = mom
  end do

#ifdef PRECISION_dp
  call tnten_cll(T2dim(1,:), tuv, momenta_TI, momenta_invariants(momenta_TI), masses_2, size(masses_2), 0)
#endif

#else
  call ol_fatal('in scalar_integral: Collier is not available')
#endif
end subroutine scalar_integral



! ***************************************************
subroutine fake_tensor_integral(rank, momenta, masses_2, Gtensor, M2)
! Build a tensor as the direct product of the loop momentum with itself up to rank 'rank',
! divided by the denominator of the corresponding tensor integral.
! ***************************************************
! rank     = highest rank in the loop
! momenta  = list of the N-1 momenta flowing inside the loop without the loop momentum.
!            momenta(N) = 0 is not included.
! masses_2 = list of the N squared masses of the loop propagators
! TI       = array containing all independent tensor components up to highest rank (contravariant light-cone);
!            divided by the denominator 'den' corresponding to the loop with 'momenta' and 'masses_2'.
! den      = product((momenta(i)+Qloop)^2 - masses_2(i), i=1..N)
! The loop momentum 'pseudotree_momentum' is taken from the module pseudotree
! QloopLC  = loop momentum in light-cone representation
! ***************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR, tensor_size
  use ol_pseudotree_/**/REALKIND, only: pseudotree_momentum
  use ol_contractions_/**/REALKIND, only: cont_V
  use ol_kinematics_/**/REALKIND, only: Std2LC_Rep
  implicit none
  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses_2(:), Gtensor(:)
  complex(REALKIND), intent(out) :: M2
  complex(REALKIND) :: TI(size(Gtensor)), QloopLC(4), den
  integer           :: i, l

  call Std2LC_Rep(pseudotree_momentum, QloopLC)
  den = cont_V(QloopLC) - masses_2(1)
  do i = 1, size(masses_2)-1
    den = den * (cont_V(momenta(:,i) + QloopLC) - masses_2(i+1))
  end do
  TI(1) = 1 / den
  do l = 1, tensor_size(rank-1)
    do i = 1,4
      TI(HR(i,l)) = QloopLC(i)*TI(l)
    end do
  end do

  M2 = tensor_contract(Gtensor, TI)

end subroutine fake_tensor_integral



! ***************************************************
subroutine TI_call(rank, momenta, masses_2, Gsum, M2)
! ***************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/DREALKIND, only: a_switch, &
    & ti_monitor, pid_string, stability_logdir, max_parameter_length
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx, momenta_invariants
  use ol_generic, only: to_string
  implicit none
  integer,           intent(in)    :: rank
  complex(REALKIND), intent(in)    :: momenta(:,:), masses_2(:), Gsum(:)
  real(REALKIND),    intent(inout) :: M2
  complex(REALKIND) :: M2add
  integer :: outunit = 44, k
  character(len=max_parameter_length) :: outfile
  complex(REALKIND) :: momenta_TI(0:3,size(momenta,2))
  if (a_switch == 0) then
    call fake_tensor_integral(rank, momenta, masses_2, Gsum, M2add)
  else if (a_switch == 1 .or. a_switch == 7) then
    ! COLI/DD: full tensor integral
    if (rank > size(masses_2) .and. a_switch == 7) then
      call ol_fatal("rank > #loop-propagators not supported by TI library DD (part of COLLIER).")
    end if
    call tensor_integral_contract(rank, momenta, masses_2, Gsum, M2add)
  else if (a_switch == 3) then
    ! COLI/DD: only scalar integrals
    call scalar_integral(momenta, masses_2)
    M2add = sum(Gsum)
  else if (a_switch == 4) then
    ! Do nothing
    M2add = sum(Gsum)
  else if (a_switch == 5) then
    ! CutTools
    call cuttools_interface(rank, momenta, masses_2, Gsum, M2add)
  else
    call ol_fatal('in TI_call: amp_switch out of range: ' // to_string(a_switch))
  end if
  M2 = M2 + real(M2add)
  if (ti_monitor > 0) then
    outfile = trim(stability_logdir) // "/ti_monitor_" // trim(to_string(a_switch)) // ".log"
    open(unit=outunit, file=outfile, form='formatted', position='append')
    if (ti_monitor > 1) write(outunit,*) ' '
    write(outunit,*) 'm2add= ', real(M2add), M2
    if (ti_monitor > 1) then
      write(outunit,*) 'rank= ', rank
      write(outunit,*) 'masses2= ', masses_2
      do k = 1, size(momenta,2)
        call LC2Std_Rep_cmplx(momenta(:,k), momenta_TI(:,k))
        write(outunit,*) 'p= ', real(momenta_TI(:,k))
      end do
      write(outunit,*) 'mominv= ', real(momenta_invariants(momenta_TI))
    end if
    close(outunit)
  end if
end subroutine TI_call


!************************************************************************************
subroutine TI_call_OL(qt_pow, rank, momenta, masses, Gsum_hcl, M2, scboxes, all_scboxes)
!************************************************************************************
! Calculation of closed one-loop integrals in the on-the-fly reduction mode.
!------------------------------------------------------------------------------------
! qt_pow      : power of \tilde{q}^2. If not zero this is a R1 rational integral
! rank        : tensor integral rank
! momenta     : set of indices used to read out the momenta from the array
!               of internal momenta L
! masses      : internal masses ids
! Gsum        : coefficient of the tensor integral
! M2          : squared matrix element
! scboxes     : indices of the scalar boxes used for the final OPP-like reduction
! all_scboxes : values of already computed scalar boxes
!************************************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND
  use ol_data_types_/**/REALKIND, only: scalarbox,hcl,met
  use ol_parameters_decl_/**/DREALKIND, only: a_switch
  use ol_loop_parameters_decl_/**/REALKIND, only: de1_IR, de2_i_IR
  use ol_loop_reduction_/**/REALKIND, only: TI_reduction, OPP_reduction, scalar_MIs
#ifdef PRECISION_dp
  use ol_loop_handling_/**/REALKIND, only: req_dp_cmp, req_qp_cmp, upgrade_qp
  use ol_loop_reduction_/**/REALKIND, only: upgrade_scalar_box
  use ol_loop_reduction_/**/QREALKIND, only: TI_reduction_qp=>TI_reduction
  use ol_loop_routines_/**/QREALKIND, only: TI_call_qt2_qp=>TI_call_qt2, &
                                            TI_call_qp=>TI_call
#endif
  use ol_parameters_decl_/**/REALKIND, only: hybrid_zero_mode, hp_switch
  implicit none
  integer,                   intent(in)    :: qt_pow,rank,momenta(:),masses(:)
  type(met),                 intent(inout) :: M2
  type(hcl),                 intent(inout) :: Gsum_hcl
  integer, optional,         intent(in)    :: scboxes(:)
  type(scalarbox), optional, intent(inout) :: all_scboxes(:)
  complex(REALKIND) :: TI(size(Gsum_hcl%cmp)), p_std(0:3,1:size(momenta)-1)
  complex(REALKIND) :: M2add, box(0:2)
  logical :: single_box
#ifdef PRECISION_dp
  integer :: u
  complex(QREALKIND) :: M2add_qp,box_qp(0:2)
#endif

  if (Gsum_hcl%mode .eq. hybrid_zero_mode) return
  M2%ndrs = M2%ndrs + Gsum_hcl%ndrs
  M2%nred = M2%nred + Gsum_hcl%nred
#ifdef PRECISION_dp
  if (hp_switch .eq. 1) then
    M2%ndrs_qp = M2%ndrs_qp + Gsum_hcl%ndrs_qp
    M2%nred_qp = M2%nred_qp + Gsum_hcl%nred_qp
  end if
#endif

  !! R1 rational term integral
  if(qt_pow > 0) then
    call TI_call_qt2(qt_pow,rank,momenta,masses,Gsum_hcl%cmp,M2%cmp)
#ifdef PRECISION_dp
    if (req_qp_cmp(Gsum_hcl)) then
      call TI_call_qt2_qp(qt_pow,rank,momenta,masses,Gsum_hcl%cmp_qp,M2%cmp_qp)
    end if
#endif
    return
  end if

  if(size(masses) .ge. 10) then
    call external_library
    return
  end if

  !! Single scalar box that has already been computed
  single_box = present(scboxes)
  if (single_box) single_box = size(scboxes)==1
  if (single_box) then
#ifdef PRECISION_dp
    if (hp_switch .eq. 1 .and. (.not. req_qp_cmp(Gsum_hcl))) then
      if (all_scboxes(scboxes(1))%qp_computed) call upgrade_qp(Gsum_hcl)
    end if
#endif

#ifdef PRECISION_dp
    if (req_dp_cmp(Gsum_hcl)) then
#endif
    box = all_scboxes(scboxes(1))%poles
    if(a_switch == 1 .or. a_switch == 7) then
      M2add = Gsum_hcl%cmp(1)*box(0)
    else if(a_switch == 5) then
      M2add = Gsum_hcl%cmp(1)*(box(0) + box(1)*de1_IR + box(2)*de2_i_IR)
    end if
    M2%cmp = M2%cmp + real(M2add)
    M2%sicount = M2%sicount + 1
#ifdef PRECISION_dp
    end if
#endif

#ifdef PRECISION_dp
    if (req_qp_cmp(Gsum_hcl)) then
      call upgrade_scalar_box(all_scboxes(scboxes(1)))
      box_qp = all_scboxes(scboxes(1))%poles_qp
      M2add_qp = Gsum_hcl%cmp_qp(1)*(box_qp(0) + box_qp(1)*de1_IR + box_qp(2)*de2_i_IR)
      M2%cmp_qp = M2%cmp_qp + real(M2add_qp,kind=QREALKIND)
      M2%sicount_qp = M2%sicount_qp + 1
    end if
#endif
    return
  end if

  !! Integration of N-point integrals with N <= 4 and reduction of N-point with N > 4.
  if(size(masses) .le. 4) then
    call scalar_MIs(momenta,masses,Gsum_hcl,M2)
  else
    call OPP_reduction(rank,momenta,masses,Gsum_hcl,M2,scboxes,all_scboxes)
  end if

  contains

  subroutine external_library()
    use ol_kinematics_/**/REALKIND, only: get_mass2
    use ol_momenta_decl_/**/REALKIND, only: L
#ifdef PRECISION_dp
    use ol_parameters_decl_/**/DREALKIND, only: coli_cache_use
    use ol_kinematics_/**/QREALKIND, only: get_mass2_qp=>get_mass2
#endif
    integer :: i,k
    complex(QREALKIND) :: p(1:5,1:size(momenta)-1)
    !! Reduction of 10-point functions not available yet.
    !! The chosen external reduction library is called for this purpose
    call ol_msg(1,"Reduction of N-point integral with N >= 10 not available. Using external library")

    i = momenta(1)
    p(1:4,1) = L(1:4,i)
    p(5,1) = L(5,i) + L(6,i)
    do k = 2, size(momenta) - 1
      i = i + momenta(k)
      p(1:4,k) = L(1:4,i)
      p(5,k) = L(5,i) + L(6,i)
    end do

#ifdef PRECISION_dp
    if (req_dp_cmp(Gsum_hcl) .or. coli_cache_use .eq. 1) then
#endif
    call TI_call(rank, cmplx(p(1:4,1:size(momenta)-1),kind=REALKIND), &
                 get_mass2(masses), Gsum_hcl%cmp, M2%cmp)
#ifdef PRECISION_dp
    end if
    if (req_qp_cmp(Gsum_hcl)) then
      call TI_call_qp(rank,p(1:4,1:size(momenta)-1),get_mass2_qp(masses), &
                      Gsum_hcl%cmp_qp,M2%cmp_qp)
    end if
#endif
    return
  end subroutine external_library
end subroutine TI_call_OL


! *******************************************************************
subroutine TI_call_qt2(qt2power, rank, momenta, masses, Gsum, M2)
! *******************************************************************
! TI-call with integrand propto qtilde^2 (qtilde^2)^qt2power.
! R1 rational term integral
! *******************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_parameters_decl_/**/REALKIND, only: CI
  use ol_kinematics_/**/REALKIND, only: get_LC_mass2,get_LC_4,get_mass2
  use ol_debug, only: ol_fatal, ol_msg, ol_error
  implicit none
  integer,           intent(in)    :: qt2power, rank
  !complex(REALKIND), intent(in)    :: momenta(:,:), masses_2(:), Gsum(:)
  complex(REALKIND), intent(in)    :: Gsum(:)
  integer,           intent(in)    :: momenta(:), masses(:)
  real(REALKIND),    intent(inout) :: M2
  complex(REALKIND) :: M2add
  complex(REALKIND) :: p1p1, p12(4), zero

  zero = 0._/**/REALKIND

  if(rank==0) then
    if(qt2power==1) then
      if (size(momenta)==2) then
        p1p1 = get_LC_mass2(momenta(1))
        M2add = - (get_mass2(masses(1)) + get_mass2(masses(2)) - p1p1/3)*Gsum(1)/2
      else if(size(momenta)==3) then
        M2add = - Gsum(1)/2
      else
        call ol_error('in TI_call_qt2: rank=0, qt2power=1, number of propagators !=2,3')
        M2add = zero
      end if
    else if (qt2power==2) then
      if(size(momenta)==4) then
        M2add = - Gsum(1)/6
      else
        call ol_error('in TI_call_qt2: rank=0, qt2power=2, number of propagators !=4')
        M2add = zero
      end if
    end if
  else if (rank==1 ) then
    if (size(momenta)==3 .AND. qt2power==1) then
      p12 = get_LC_4(momenta(1)) + get_LC_4(momenta(1) + momenta(2))
      M2add = - Gsum(1)/2 + SUM(Gsum(2:5)*p12(1:4))/6
    else
      call ol_error('in TI_call_qt2: rank=1, qt2power!=1 OR number of propagators !=3')
      M2add = zero
    end if
  else
    call ol_error('in TI_call_qt2: R1 integral with rank > 1')
    M2add = zero
  end if

  M2 = M2 + real(M2add)

end subroutine TI_call_qt2


! ****************************************************
function TI2_call(rank, momenta, masses_2, Gsum, TI)
! TI_call with precalculated tensor integrals.
! Used for loop^2 processes.
! Returns the contribution to the amplitude
! ****************************************************
  use KIND_TYPES, only: REALKIND
  use ol_debug, only: ol_fatal, ol_msg, ol_error
  use ol_parameters_decl_/**/DREALKIND, only: a_switch
  use ol_generic, only: to_string
  implicit none
  integer,           intent(in)    :: rank
  complex(REALKIND), intent(in)    :: momenta(:,:), masses_2(:), Gsum(:), TI(:)
  complex(REALKIND) :: TI2_call
  if (a_switch == 0) then
    call fake_tensor_integral(rank, momenta, masses_2, Gsum, TI2_call)
  else if (a_switch == 1 .or. a_switch == 7) then
    ! COLI/DD: contract with precalculated tensor integral
    TI2_call = tensor_contract(Gsum, TI)
  else if (a_switch == 4) then
    ! Do nothing
    TI2_call = sum(Gsum)
  else if (a_switch == 5) then
    ! CutTools
    call cuttools_interface(rank, momenta, masses_2, Gsum, TI2_call)
  else
    call ol_error(2, 'in TI2_call: amp_switch out of range: ' // to_string(a_switch))
    call ol_msg('note that modes 2 and 3 are not supported in loop^2.')
    call ol_fatal()
  end if
end function TI2_call



! **********************************************************************
function tensor_contract(G, TI)
! Contract two tensors from rank 0 to the highest rank of G.
! TI must be of equal or higher rank as G.
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: G(:), TI(:)
  complex(REALKIND) :: tensor_contract
  tensor_contract = sum(G*TI(1:size(G)))
end function tensor_contract



! ****************************************
subroutine loop_trace(G_in, Gtensor)
! ****************************************
! Close spinor or vector loop line after all vertex insertions
! Gtensor(l) = Tr[G_in(beta,l,alpha)] = Sum(G_in(alpha,l,alpha), alpha=1..4)
! alpha = covariant (light-cone)
! beta  = contravariant (light-cone)
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_in(:,:,:)
  complex(REALKIND), intent(out) :: Gtensor(:)
  Gtensor = G_in(1,:,1) + G_in(2,:,2) + G_in(3,:,3) + G_in(4,:,4)
end subroutine loop_trace



! *********************************************
subroutine loop_cont_VV(G_in, J_V2, J_V1, G_out)
! *********************************************
! contraction of tensor coefficients with vector currents for pseudo-tree
! G_out(l) =  J_V1(alpha) * G_in(beta,l,alpha) * J_V2(beta)
! J_Vi  = CONTRAVARIANT current
! alpha = COVARIANT index
! beta  = CONTRAVARIANT index
  use KIND_TYPES, only: REALKIND
  use ol_contractions_/**/REALKIND, only: cont_VV
  implicit none
  complex(REALKIND), intent(in)  :: G_in(:,:,:), J_V1(4), J_V2(4)
  complex(REALKIND), intent(out) :: G_out(:)
  integer :: l

  do l = 1, size(G_in,2)
    G_out(l) = J_V1(1) * cont_VV(G_in(:,l,1), J_V2) + J_V1(2) * cont_VV(G_in(:,l,2), J_V2) &
             + J_V1(3) * cont_VV(G_in(:,l,3), J_V2) + J_V1(4) * cont_VV(G_in(:,l,4), J_V2)
  end do

end subroutine loop_cont_VV



! *************************************
subroutine loop_cont_QA(G_in, Q, A, G_out)
! *************************************
! Contraction of tensor coefficients with quark-antiquark currents (attach loop line)
! G_out(l) = Q(alpha) * G_in(beta,l,alpha) * A(beta)
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_in(:,:,:), Q(4), A(4)
  complex(REALKIND), intent(out) :: G_out(:)

  G_out = Q(1)*A(1)*G_in(1,:,1) + Q(2)*A(1)*G_in(2,:,1) + Q(3)*A(1)*G_in(3,:,1) + Q(4)*A(1)*G_in(4,:,1) &
        + Q(1)*A(2)*G_in(1,:,2) + Q(2)*A(2)*G_in(2,:,2) + Q(3)*A(2)*G_in(3,:,2) + Q(4)*A(2)*G_in(4,:,2) &
        + Q(1)*A(3)*G_in(1,:,3) + Q(2)*A(3)*G_in(2,:,3) + Q(3)*A(3)*G_in(3,:,3) + Q(4)*A(3)*G_in(4,:,3) &
        + Q(1)*A(4)*G_in(1,:,4) + Q(2)*A(4)*G_in(2,:,4) + Q(3)*A(4)*G_in(3,:,4) + Q(4)*A(4)*G_in(4,:,4)

end subroutine loop_cont_QA



! *****************************
subroutine G0initialisation(G0)
! *****************************
! Initialise the rank 0 tensor coefficient of the cut spinor or vector loop line;
! parameterised only by loop momentum
! G0(beta,1,alpha) = delta(alpha,beta)
! index notation: G0(beta,l,alpha)
! l = tensor index
! alpha = covariant (light-cone) "frozen" open index, is untouched till the last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/props to build the loop
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(out) :: G0(4,1,4)
  G0 = 0
  G0(1,1,1) = 1
  G0(2,1,2) = 1
  G0(3,1,3) = 1
  G0(4,1,4) = 1
end subroutine G0initialisation

! *****************************
subroutine G0initialisationOLR(preFactor,G0)
! *****************************
! Initialise the rank 0 tensor coefficient of the cut spinor or vector loop line;
! parameterised only by loop momentum
! G0(beta,1,alpha) = delta(alpha,beta)*preFactor
! index notation: G0(beta,l,alpha)
! l = tensor index
! alpha = covariant (light-cone) "frozen" open index, is untouched till the last contraction in the loop
! beta  = contravariant (light-cone) "active" index contracted with vertices/props to build the loop
  use KIND_TYPES, only: REALKIND
  implicit none
  complex(REALKIND), intent(in) :: preFactor
  complex(REALKIND), intent(out) :: G0(4,1,4)
  G0 = 0
  G0(1,1,1) = preFactor
  G0(2,1,2) = preFactor
  G0(3,1,3) = preFactor
  G0(4,1,4) = preFactor
end subroutine G0initialisationOLR


subroutine cts_numerator(q, amp)
  use KIND_TYPES, only: REALKIND
  use ol_loop_momentum_/**/REALKIND, only: loop_mom_tens
  use ol_tensor_storage_/**/REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: q(0:3)
  complex(REALKIND), intent(out) :: amp
  complex(REALKIND) :: Qtensor(array_length_stored), QloopLC(1:4)
  call loop_mom_tens(q, Qtensor)
  amp = tensor_contract(Qtensor, tensor_stored) ! contract up to the length of Qtensor
end subroutine cts_numerator


subroutine cuttools_interface(rank, momenta, masses2, Gtensor, M2)
  use KIND_TYPES, only: REALKIND, DREALKIND
  use ol_loop_parameters_decl_/**/DREALKIND, only: opprootsvalue, mureg, de1_UV, de1_IR, de2_i_IR
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep
  use ol_tensor_storage_/**/REALKIND
  use ol_tensor_bookkeeping, only: tensor_size
#ifdef USE_CUTTOOLS
  use cts_numdummies, only: dpnumdummy, mpnumdummy
#endif
  implicit none

  integer,           intent(in)  :: rank
  complex(REALKIND), intent(in)  :: momenta(:,:), masses2(:), Gtensor(:)
  complex(REALKIND), intent(out) :: M2
#ifdef USE_CUTTOOLS
  complex(DREALKIND) :: cts_amp_array(0:2), cts_ampcc, cts_ar1
  real(REALKIND)     :: mom(0:3)
  complex(DREALKIND) :: masses2_dp(size(masses2))
  real(DREALKIND)    :: cts_pp(0:3,0:size(masses2)-1)
  integer            :: l
  logical            :: cts_stable

  if (de1_UV /= de1_IR) then
    call ol_fatal('pole1_UV != pole1_IR is not allowed with CutTools.')
  end if

  tensor_stored(:size(Gtensor)) = Gtensor
  rank_stored = rank
  array_length_stored = tensor_size(rank)

  masses2_dp = masses2

  cts_pp(:,0) = 0
  do l = 1, size(momenta,2)
    call LC2Std_Rep(momenta(:,l), mom)
    cts_pp(:,l) = mom
  end do

#ifdef PRECISION_dp
  call ctsxcut(3, opprootsvalue, mureg, size(masses2), cts_numerator, mpnumdummy, rank, &
             & cts_pp, masses2_dp, cts_amp_array, cts_ampcc, cts_ar1, cts_stable)
#else
  call ctsxcut(6, opprootsvalue, mureg, size(masses2), dpnumdummy, cts_numerator, rank, &
             & cts_pp, masses2_dp, cts_amp_array, cts_ampcc, cts_ar1, cts_stable)
#endif

  M2 = cts_amp_array(0) + cts_amp_array(1)*de1_IR + cts_amp_array(2)*de2_i_IR
#else
  M2 = 0 ! prevent compiler warning
  call ol_fatal('cuttools_interface: CutTools is not available')
! #ifdef USE_CUTTOOLS
#endif
end subroutine cuttools_interface

end module ol_loop_routines_/**/REALKIND
