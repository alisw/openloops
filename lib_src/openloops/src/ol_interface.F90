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


module openloops
  use KIND_TYPES, only: DREALKIND
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_char, c_int, c_double, c_null_char
  use ol_init, only: set_init_error_fatal, set_parameter, get_parameter, parameters_flush, &
      & tree_parameters_flush, cleanup, set_if_modified
  use ol_version, only: welcome, openloops_version_string
  use ol_parameters_decl_/**/DREALKIND,  only: procname_length, max_parameter_length
  use ol_external_decl_/**/DREALKIND, only: n_scatt
  use ol_debug, only: get_error, ol_msg, ol_error, ol_fatal, ol_fatal_update, error
  use ol_parameters_init_/**/DREALKIND, only: parameters_write
  implicit none
  private
  ! from module init_iu
  public :: set_init_error_fatal, get_error
  public :: set_parameter, get_parameter, parameters_flush, tree_parameters_flush
  ! from module use ol_version
  public :: welcome, openloops_version_string
  ! process interface
  public :: n_external, amplitudetype, phase_space_point, start, finish
  public :: tree_colbasis_dim, tree_colbasis, tree_colourflow
  public :: register_process, register_process_id
  public :: evaluate_tree
  public :: evaluate_cc, evaluate_loopcc, evaluate_cc2
  public :: evaluate_ccmatrix, evaluate_loopccmatrix, evaluate_ccmatrix2
  public :: evaluate_ccewmatrix, evaluate_ccewmatrix2
  public :: evaluate_sc, evaluate_loopsc, evaluate_sc2
  public :: evaluate_sctensor, evaluate_sctensor2, evaluate_loopsctensor
  public :: evaluate_scpowheg, evaluate_loopscpowheg, evaluate_scpowheg2
  public :: evaluate_stensor, evaluate_loopstensor, evaluate_stensor2
  public :: evaluate_tree_colvect, evaluate_tree_colvect2
  public :: evaluate_full, evaluate_loop, evaluate_loop2
  public :: evaluate_loopbare, evaluate_ct, evaluate_loopct, evaluate_r2
  public :: evaluate_iop, evaluate_iop2
  public :: evaluate_pt, evaluate_schsf
  public :: evaluate_associated
  public :: evaluate_poles
  ! Print parameters
  public :: ol_printparameter
  ! used in BLHA interface
  public :: rval_size, stop_invalid_id

  interface register_process
    module procedure register_process_string, register_process_id
  end interface register_process

  interface printparameter
    module procedure ol_printparameter
  end interface printparameter

  interface evaluate_stensor
    module procedure evaluate_scpowheg
  end interface evaluate_stensor

  interface evaluate_loopstensor
    module procedure evaluate_loopscpowheg
  end interface evaluate_loopstensor

  interface evaluate_stensor2
    module procedure evaluate_scpowheg2
  end interface evaluate_stensor2


  type process_handle
    integer :: n_particles = 0
    integer :: max_point = -1
    integer :: tensor_rank = -1
    integer :: qcd_powers(2)
    ! allocatable length character members are not supported in gfortran (tested with 4.8.1)
    character(len=procname_length) :: process_name
    character(len=max_parameter_length) :: library_name
    integer, allocatable :: permutation(:)
    integer, allocatable :: pol(:)
    integer, allocatable :: extid(:)
    integer, allocatable :: photon_id(:)
    type(c_ptr) :: library_handle = c_null_ptr
    integer :: amplitude_type ! 1=Tree, 11=Loop, 12=LoopInduced
    integer :: content = 0 ! bitwise: 2^0=tree, 2^1=loop, 2^2=loop2, 2^3=pt
    integer :: n_in = 2 ! Phase-space for n_in -> n-n_in
    integer :: associated_ew = 0
    integer :: associated_born_1 = 0
    integer :: associated_born_2 = 0
    integer :: replace_loop = 0
    integer :: stability_mode = -1
    integer :: a_switch = -1
    integer :: a_switch_rescue = -1
    integer :: redlib_qp = -1
    real(DREALKIND), allocatable :: masses(:)
    real(DREALKIND), allocatable :: last_psp(:,:)
    real(DREALKIND) :: last_muren
    real(DREALKIND) :: last_alpha_QCD
    integer :: loop_parameters_status
    integer, allocatable :: last_perm(:)
    integer, allocatable :: last_pol(:)
    logical :: last_zero = .true.
    procedure(), pointer, nopass :: set_permutation => null()
    procedure(), pointer, nopass :: pol_init => null()
    procedure(), pointer, nopass :: set_photons => null()
    procedure(), pointer, nopass :: tree => null()
    procedure(), pointer, nopass :: loop => null()
    procedure(), pointer, nopass :: loopcr => null()
    procedure(), pointer, nopass :: ct => null()
    procedure(), pointer, nopass :: pt => null()
    procedure(), pointer, nopass :: schsf => null()
    procedure(), pointer, nopass :: rambo => null()
    procedure(), pointer, nopass :: tree_colbasis_dim => null()
    procedure(), pointer, nopass :: tree_colbasis => null()
    procedure(), pointer, nopass :: tree_colvect => null()
    procedure(), pointer, nopass :: iop => null()
    procedure(), pointer, nopass :: loopcc => null()
  end type process_handle

  ! process handle array
  integer, save :: last_process_id = 0
  type(process_handle), save, allocatable :: process_handles(:)

  type processinfos
    integer :: EWorder(0:1)
    integer :: QCDorder(0:1)
    integer :: OLMode
    integer :: API
    integer :: LeadingColour
    integer :: NF
    integer :: NC
    integer :: CKMORDER
    integer :: POLSEL
    character :: ME, MM, ML, MU, MD, MS, MC, MB
    integer :: YE, YM, YL, YU, YD, YS, YC, YB, YT
    character :: CC
    character(len=max_parameter_length) :: LIBNAME
    character(len=127) :: MODEL
    character(len=127) :: PROC
    character(len=127) :: MAP
    character(len=127) :: MAPPERM
    character(len=127) :: APPROX
    character(len=4) :: ID
    character(len=4) :: TYPE
    character(len=4) :: LTYPE
  end type processinfos
  type(processinfos), save, allocatable :: process_infos(:)
  type(processinfos), save, allocatable :: loaded_libs(:)

  type extparticle
    integer :: id
    integer :: pol
    logical :: is_initial
  end type extparticle

  ! array for shopping list
  character(len=max_parameter_length), save, allocatable :: shopped_processes(:)
  logical, save :: shopping_list_open = .false.
  integer, parameter :: fh_shopping = 998
  ! PDG IDs of particles which carry adjoint colour
  integer :: pdgadjoint(1) = [21]

#if __APPLE__
  character(len=5), parameter :: dynlib_extension='dylib'
#else
  character(len=2), parameter :: dynlib_extension='so'
#endif

  character(len=4) :: loops_flags = "tlsp" ! used in this order to set content bits


  contains


  pure function rval_size(n_part, amp_type)
    implicit none
    integer :: rval_size
    integer, intent(in) :: n_part, amp_type
    select case (amp_type)
      case  (1, 12) ! Tree, LoopInduced
        rval_size = 1
      case (11) ! Loop
        rval_size = 4
      case (2) ! ccTreeo
        rval_size = (n_part*(n_part-1))/2
      case (3) ! scTree
        rval_size = 2*n_part*n_part
      case (4) ! scTree_polvect
        rval_size = n_part
      case (0)
        rval_size = 0
      case default
        ! unknown amp_type
        rval_size = 0
    end select
  end function rval_size


  function get_process_handle(lib, libname, proc, content, amptype, n_in, perm, pol, extid, photon_id, qcd_powers)
    ! [in] lib: a shared library handle
    ! [in] proc: a full process name, '<lib>_<subproc>_<id>'
    ! [in] perm: integer array with the crossing
    ! [in] content: integer with binary tags for tree, loop, loop2, pt
    ! [in] amptype: integer to specify BLHA matrix element type
    ! return process handle of type process_handle
    ! note: error handling is done in dlsym
    use KIND_TYPES, only: DREALKIND
    use ol_dlfcn, only: dlsym
    use ol_loop_parameters_decl_/**/DREALKIND, only: &
         & stability_mode,a_switch,a_switch_rescue,redlib_qp,hel_mem_opt_switch
    implicit none
    type(c_ptr), intent(in) :: lib
    character(len=*), intent(in) :: libname
    character(len=*), intent(in) :: proc
    integer, intent(in) :: content, amptype, n_in
    integer, intent(in), optional :: perm(:)
    integer, intent(in), optional :: pol(:)
    integer, intent(in), optional :: extid(:)
    integer, intent(in), optional :: photon_id(:)
    integer, intent(in), optional :: qcd_powers(2)
    type(process_handle) :: get_process_handle
    integer :: k
    procedure(), pointer :: tmp_fun

    ! number of external particles
    tmp_fun => dlsym(lib, "ol_f_n_external_" // trim(proc))

    call tmp_fun(get_process_handle%n_particles)
    get_process_handle%library_name = trim(libname)
    get_process_handle%process_name = trim(proc)
    allocate(get_process_handle%permutation(get_process_handle%n_particles))
    if (present(perm)) then
      ! check correct size of the permutation
      if (get_process_handle%n_particles /= size(perm)) then
        call ol_fatal('error: registered process with wrong size of particle permutation')
        return
      end if
      get_process_handle%permutation = perm
    else
      get_process_handle%permutation = [(k, k=1, get_process_handle%n_particles)]
    end if
    get_process_handle%library_handle = lib
    get_process_handle%set_permutation => dlsym(lib, "ol_f_set_permutation_" // trim(proc))
    get_process_handle%rambo => dlsym(lib, "ol_f_rambo_" // trim(proc))
    get_process_handle%amplitude_type = amptype
    get_process_handle%tree => dlsym(lib, "ol_f_amp2_" // trim(proc))
    get_process_handle%loop => dlsym(lib, "ol_f_vamp2_" // trim(proc))
    get_process_handle%ct => dlsym(lib, "ol_f_ctamp2_" // trim(proc))
    get_process_handle%pt => dlsym(lib, "ol_f_ptamp2_" // trim(proc))
    get_process_handle%schsf => dlsym(lib, "ol_f_schsfamp2_" // trim(proc))
    get_process_handle%content = content
    get_process_handle%n_in = n_in
    get_process_handle%qcd_powers = qcd_powers
    ! loop correlators
    get_process_handle%loopcr => dlsym(lib, "ol_f_vampcr2_" // trim(proc))
    allocate(get_process_handle%last_psp(4,get_process_handle%n_particles))
    get_process_handle%last_psp=0
    allocate(get_process_handle%last_perm(get_process_handle%n_particles))
    get_process_handle%last_perm=0
    allocate(get_process_handle%last_pol(get_process_handle%n_particles))
    get_process_handle%last_pol=0
    get_process_handle%loop_parameters_status = -1
    ! external masses and highest tensor rank
    tmp_fun => dlsym(lib, "ol_f_get_masses_" // trim(proc))
    allocate(get_process_handle%masses(get_process_handle%n_particles))
    call tmp_fun(get_process_handle%masses)
    allocate(get_process_handle%pol(get_process_handle%n_particles))
    get_process_handle%pol_init => dlsym(lib, "ol_f_pol_init_" // trim(proc))
    if (present(pol)) then
      ! check correct size of the polarization vector
      if (get_process_handle%n_particles /= size(pol)) then
        call ol_fatal('error: registered process with wrong size of polarization vector')
        return
      end if
      get_process_handle%pol = pol
    else
      get_process_handle%pol = 0
    end if
    allocate(get_process_handle%extid(get_process_handle%n_particles))
    if (present(extid)) then
      ! check correct size of the extid vector
      if (get_process_handle%n_particles /= size(extid)) then
        call ol_fatal('error: registered process with wrong size of extid')
        return
      end if
      get_process_handle%extid = extid
    else
      get_process_handle%extid = 0
    end if
    allocate(get_process_handle%photon_id(get_process_handle%n_particles))
    if (present(photon_id)) then
      ! check correct size of the photon_id vector
      if (get_process_handle%n_particles /= size(photon_id)) then
        call ol_fatal('error: registered process with wrong size of extid')
        return
      end if
      get_process_handle%photon_id = photon_id
    else
      get_process_handle%photon_id = 0
    end if
    get_process_handle%set_photons => dlsym(lib, "ol_f_set_photons_" // trim(proc))
    if (btest(content, 1)) then
      tmp_fun => dlsym(lib, "ol_f_max_point_" // trim(proc))
      call tmp_fun(get_process_handle%max_point)
      tmp_fun => dlsym(lib, "ol_f_tensor_rank_" // trim(proc))
      call tmp_fun(get_process_handle%tensor_rank)
    end if
    ! colour basis
    get_process_handle%tree_colbasis_dim => dlsym(lib, "ol_tree_colbasis_dim_" // trim(proc))
    get_process_handle%tree_colbasis => dlsym(lib, "ol_tree_colbasis_" // trim(proc))
    get_process_handle%tree_colvect => dlsym(lib, "ol_tree_colvect_" // trim(proc))
    ! stability settings
    get_process_handle%a_switch=a_switch
    get_process_handle%a_switch_rescue=a_switch_rescue
    get_process_handle%redlib_qp=redlib_qp
    get_process_handle%stability_mode=stability_mode
    ! optimized helicity bookkeeping
    tmp_fun => dlsym(lib, "ol_hel_mem_opt_" // trim(proc))
    if (associated(tmp_fun)) then
      call tmp_fun(hel_mem_opt_switch)
    else
!      call ol_fatal_update()
    end if
    tmp_fun => dlsym(lib, "ol_f_iopamp2_" // trim(proc))
    if (associated(tmp_fun)) then
      get_process_handle%iop => tmp_fun
    end if
    get_process_handle%loopcc => dlsym(lib, "ol_loopcc_" // trim(proc))
  end function get_process_handle

  function register_process_lib(libname, proc, content, amptype, n_in, pol, perm, extid, photon_id, qcd_powers)
    ! [in] libname: name of the process library
    ! [in] proc: a full process name, '<lib>_<subproc>_<id>'
    ! [in] perm: integer array with the crossing
    ! [in] content: integer with binary tags for tree, loop, loop2, pt
    ! [in] amptype: integer to specify BLHA matrix element type
    ! return (integer) process id to be used in OLP_EvalSubProcess
    use KIND_TYPES, only: DREALKIND
    use ol_dlfcn, only: dlopen, RTLD_LAZY
    use ol_loop_parameters_decl_/**/DREALKIND, only: maxpoint, maxrank, do_pole_checks
    implicit none
    character(len=*), intent(in) :: libname
    character(len=*), intent(in) :: proc
    integer, intent(in) :: content, amptype, n_in
    integer, intent(in), optional :: pol(:)
    integer, intent(in), optional :: perm(:)
    integer, intent(in), optional :: extid(:)
    integer, intent(in), optional :: photon_id(:)
    integer, intent(in), optional :: qcd_powers(2)
    type(c_ptr) :: lib
    logical :: same_perm, same_pol, same_photon_id
    integer :: register_process_lib
    integer :: j, k
    type(process_handle) :: prochandle
    type(process_handle), allocatable :: process_handles_bak(:)

    lib = dlopen(libname, RTLD_LAZY, 2)
    prochandle = get_process_handle(lib, libname, proc, content, amptype, n_in, perm=perm, pol=pol, &
                                                   & extid=extid, photon_id=photon_id, qcd_powers=qcd_powers)

    if (error > 1) return
    ! Check if the process was registered before with the same permutation, polarization and amptype.
    ! If yes, return the previously assigned id
    do k = 1, last_process_id
      if ((trim(proc) /= trim(process_handles(k)%process_name)) .or. &
        & (trim(libname) /= trim(process_handles(k)%library_name)) ) then
        cycle
      end if
      if (present(perm)) then
        same_perm = all(perm == process_handles(k)%permutation)
      else
        ! perm not present means 1,2,..,n
        same_perm = all(process_handles(k)%permutation == [(j, j=1, process_handles(k)%n_particles)])
      end if
      if (present(pol)) then
        same_pol = all(pol == process_handles(k)%pol)
      else
        same_pol = all(process_handles(k)%pol == 0)
      end if
      if (present(photon_id)) then
        same_photon_id = all(photon_id == process_handles(k)%photon_id)
      else
        same_photon_id = all(process_handles(k)%photon_id == 0)
      end if
      if (same_perm .and. &
        & same_pol .and. &
        & same_photon_id .and. &
        & (amptype == process_handles(k)%amplitude_type) ) then
        register_process_lib = k
        return
      end if
    end do
    if (.not. allocated(process_handles)) then
      allocate(process_handles(1))
    end if
    if (last_process_id == size(process_handles)) then
      allocate(process_handles_bak(last_process_id))
      process_handles_bak = process_handles
      deallocate(process_handles)
      allocate(process_handles(2*last_process_id))
      process_handles(1:last_process_id) = process_handles_bak
      deallocate(process_handles_bak)
    end if
    last_process_id = last_process_id + 1
    process_handles(last_process_id) = prochandle
    if (maxpoint < process_handles(last_process_id)%max_point) then
      call set_parameter("maxpoint",process_handles(last_process_id)%max_point)
    end if
    if (maxrank < process_handles(last_process_id)%tensor_rank) then
      call set_parameter("maxrank",process_handles(last_process_id)%tensor_rank)
    end if
    register_process_lib = last_process_id

    if (do_pole_checks) call check_poles(register_process_lib)
  end function register_process_lib


  subroutine unregister_processes()
    ! Close all process libraries and nullify process handles.
    use ol_dlfcn, only: dlclose
    implicit none
    integer :: id
    do id = 1, last_process_id
      call dlclose(process_handles(id)%library_handle)
      process_handles(id)%n_particles = 0
      process_handles(id)%content = 0
      if (allocated(process_handles(id)%permutation)) deallocate(process_handles(id)%permutation)
      if (allocated(process_handles(id)%pol)) deallocate(process_handles(id)%pol)
      if (allocated(process_handles(id)%extid)) deallocate(process_handles(id)%extid)
      if (allocated(process_handles(id)%masses)) deallocate(process_handles(id)%masses)
      process_handles(id)%library_handle = c_null_ptr
      process_handles(id)%set_permutation => null()
      process_handles(id)%tree => null()
      process_handles(id)%loop => null()
      process_handles(id)%ct => null()
      process_handles(id)%pt => null()
      process_handles(id)%schsf => null()
    end do
    if (allocated(process_handles)) deallocate(process_handles)
    last_process_id = 0
  end subroutine unregister_processes


  function register_process_string(process_in, amptype)
    ! process: string with format 2->n-2
    ! amptype: integer 1,2,3,4,11,12
    ! return (integer) process id to be used in evaluate_process
    use KIND_TYPES, only: DREALKIND
    use ol_generic, only: to_int, string_to_integerlist, count_substring, to_string, to_lowercase
    use ol_parameters_decl_/**/DREALKIND, only: &
      & install_path,flavour_mapping_on, coupling_qcd, coupling_ew, &
      & write_shopping_list, add_associated_ew, loop_order_ew, loop_order_qcd, &
      & check_collection, partial_normal_order
    implicit none
    character(len=*), intent(in) :: process_in
    integer, intent(in) :: amptype
    integer :: register_process_string
    character(len=max_parameter_length) :: tmp
    character(len=max_parameter_length) :: inp, outp
    character(len=max_parameter_length) :: process, proc, libhandle, permstring
    integer :: librarytype
    integer:: check
    integer :: n_in, n_out
    type(extparticle), allocatable :: ext(:)
    integer, allocatable :: perm(:)
    integer, allocatable :: pol(:)
    integer, allocatable :: extid(:)
    integer, allocatable :: photon_id(:)
    integer :: coupling_ew_bak(0:1), coupling_qcd_bak(0:1)
    integer :: loop_order_ew_bak, loop_order_qcd_bak
    logical :: check_collection_bak
    integer :: associated_ew, replace_loop = 0
    integer :: i
    logical :: decay = .false.
    character(len=max_parameter_length) :: outstring

    call parameters_flush() ! make sure that pid_string is set
    register_process_string = -1

    call get_environment_variable("OpenLoopsPath", tmp)
    if (len_trim(tmp) /= 0) then
      call set_parameter("install_path", tmp, error)
    end if

    call ol_msg(3,"registering process: " // trim(process_in) )

    ! process: in -> out
    if (index(process_in, ">") > 0) then

      if (index(process_in, "->") > 0) then
        inp = adjustl(trim(process_in(1:index(process_in, "->")-1)))
        outp = adjustl(process_in(index(process_in, "->")+2:len(process_in)))
      else
        inp = adjustl(trim(process_in(1:index(process_in, ">")-1)))
        outp = adjustl(process_in(index(process_in, ">")+1:len(process_in)))
      end if

      n_in = size(process_to_extparticlelist(inp))
      n_out = size(process_to_extparticlelist(outp))
      if (error > 0 .or. n_in == 0 .or. n_out == 0 .or. n_in+n_out < 3) then
        call ol_error("register_process: invalid argument: " // trim(process_in) )
      end if
      allocate (ext(n_in+n_out))
      ext(1:n_in) = process_to_extparticlelist(inp, .true.)
      ext(n_in+1:) = process_to_extparticlelist(outp, .false.)

      ! check for unresolved/on-/off-shell photons
      allocate (photon_id(size(ext)))
      call check_photon_id(ext, photon_id)

      ! charge conjugate final state particles
      call charge_conj(ext)

      ! flavour mapping
      if (flavour_mapping_on > 0) then
        call flavour_mapping(ext, flavour_mapping_on)
      end if

      ! determine normal ordering
      allocate (perm(size(ext)))
      allocate (extid(size(ext)))
      proc=""
      if(partial_normal_order) then
        call normal_order(ext(:n_in), perm(:n_in), extid(:n_in), proc)
        call normal_order(ext(n_in+1:), perm(n_in+1:), extid(n_in+1:), proc)
        perm(n_in+1:)=perm(n_in+1:)+n_in
      else
        call normal_order(ext, perm, extid, proc)
      end if

      if (proc == "") then
        call ol_error("register_process: invalid argument: " // trim(process_in))
        return
      end if

      ! permute polarization states
      allocate (pol(size(ext)))
      do i=1, size(ext)
        pol(perm(i)) = ext(i)%pol
      end do

      call ol_msg(3,"check process library for: " // trim(proc) // ", " // to_string([(ext(i)%id,i=1,size(ext))]))

      if (amptype == 99 .or. write_shopping_list ) then ! write shopping list
        ! charge conjugate back final state particles to write shopping list
        call charge_conj(ext)
        register_process_string = write_shop_list(ext, proc)
      else

        ! find process
        register_process_string = loop_over_libraries(proc, amptype, n_in, perm, pol, extid, process_in, photon_id)

        coupling_ew_bak  = coupling_ew
        coupling_qcd_bak  = coupling_qcd
        loop_order_ew_bak = loop_order_ew
        loop_order_qcd_bak = loop_order_qcd

        ! for >4q processes: if not found directly try to load tree and loop
        ! amplitudes from different libraries.
        if (register_process_string < 0 .and.  nextpid(ext,[1,2,3,4,5,6]) .ge. 4) then
         if (amptype > 10 .and. coupling_ew(0) /= -1 &
                          .and. (coupling_ew(1) == -1 .or. coupling_ew(1) == 0) &
                          .and. (coupling_qcd(1) == -1  .or. coupling_qcd(1) == 1)) then
            replace_loop = loop_over_libraries(proc, 1, n_in, perm, pol, extid, process_in, photon_id)
            if (replace_loop > 0) then
              loop_order_ew = coupling_ew(0)
              coupling_ew = -1
              coupling_qcd = -1
            end if
          else if (amptype > 10 .and. coupling_qcd(0) /= -1 &
                                .and. (coupling_qcd(1) == -1 .or. coupling_qcd(1) == 0) &
                                .and. (coupling_ew(1) == -1  .or. coupling_ew(1) == 1)) then
            replace_loop = loop_over_libraries(proc, 1, n_in, perm, pol, extid, process_in, photon_id)
            if (replace_loop > 0) then
              loop_order_qcd = coupling_qcd(0)
              coupling_ew = -1
              coupling_qcd = -1
            end if
          end if

          ! interference born -> squared born (for 4q processes)
          if (replace_loop > 0) then
            register_process_string = loop_over_libraries(proc, amptype, n_in, perm, pol, extid, process_in, photon_id)
            if (register_process_string /= -1) then
              process_handles(replace_loop)%replace_loop=register_process_string
              register_process_string = replace_loop
            end if
          end if
        end if

        ! register associate ew one-loop amplitude
        if (register_process_string > 0 .and. abs(add_associated_ew) > 0 .and. coupling_ew(1) == 0 ) then
          coupling_ew(1) = 1
          coupling_qcd(1) = 0
          associated_ew = loop_over_libraries(proc, amptype, n_in, perm, pol, extid, process_in, photon_id)
          process_handles(register_process_string)%associated_ew = associated_ew
        end if

        ! register associate qcd-ew born amplitudes
        if (register_process_string > 0 .and. abs(add_associated_ew) > 1) then
          check_collection_bak = check_collection
          check_collection=.false.
          if (coupling_ew_bak(0) > -1) then
            coupling_ew(0) = coupling_ew_bak(0)+1
            coupling_qcd(0) = -1
          else if (coupling_qcd_bak(0) > -1) then
            coupling_ew(0) = -1
            coupling_qcd(0) = coupling_qcd_bak(0)-1
          else
            call ol_error("associate ew order selection not valid.")
          end if
          associated_ew = loop_over_libraries(proc, 1, n_in, perm, pol, extid, process_in, photon_id)
          process_handles(register_process_string)%associated_born_1= associated_ew
        end if

        if (register_process_string > 0 .and. abs(add_associated_ew) > 2) then
          if (coupling_ew_bak(0) > -1) then
            coupling_ew(0) = coupling_ew_bak(0)+2
            coupling_qcd(0) = -1
          else if (coupling_qcd_bak(0) > -1) then
            coupling_ew(0) = -1
            coupling_qcd(0) = coupling_qcd_bak(0)-2
          else
            call ol_error("associate ew order selection not valid.")
          end if
          associated_ew = loop_over_libraries(proc, 1, n_in, perm, pol, extid, process_in, photon_id)
          process_handles(register_process_string)%associated_born_2= associated_ew
        end if

        loop_order_ew = loop_order_ew_bak
        loop_order_qcd = loop_order_qcd_bak
        coupling_ew  = coupling_ew_bak
        coupling_qcd  = coupling_qcd_bak
        check_collection=check_collection_bak

     end if

     if (register_process_string < 1) then ! process not found
       outstring = adjustl("register_process: process " // trim(process_in))
       outstring = trim(outstring) // " @" // " tree="
       if (coupling_ew(0) > -1) then
         outstring = trim(outstring) // trim(to_string(coupling_ew(0))) // ","
       else
         outstring = trim(outstring) // "?,"
       end if
       if (coupling_qcd(0) > -1) then
         outstring = trim(outstring) // trim(to_string(coupling_qcd(0)))
       else
         outstring = trim(outstring) // "?"
       end if
       if (amptype > 10) then
          outstring = trim(outstring) // " loop="
          if (coupling_ew(1) > -1) then
            outstring = trim(outstring) // trim(to_string(coupling_ew(0)+coupling_ew(1))) // ","
          else if (loop_order_ew > -1) then
            outstring = trim(outstring) // trim(to_string(loop_order_ew)) // ","
          else
            outstring = trim(outstring) // "?,"
          end if
          if (coupling_qcd(1) > -1) then
            outstring = trim(outstring) // trim(to_string(coupling_qcd(0)+coupling_qcd(1)))
          else if (loop_order_qcd > -1) then
            outstring = trim(outstring) // trim(to_string(loop_order_qcd))
          else
            outstring = trim(outstring) // "?"
          end if
       end if
       outstring = trim(outstring) // " (EW,QCD) not found!"
       call ol_msg(outstring)
     end if

    else ! direct library loader

      ! read permutation
      permstring =  process_in(index(process_in, '[')+1:index(process_in, ']')-1)
      if (len_trim(permstring) /= 0) then
        allocate(perm(size(string_to_integerlist(permstring))))
        perm = string_to_integerlist(permstring)
        libhandle = process_in(:index(process_in,"[")-1)
      else
        libhandle = trim(process_in)
      end if

      !register
      librarytype = 0
      do
        if (allocated(perm)) then
          register_process_string =  check_process(libhandle, amptype, librarytype, 2, perm_in=perm)
          if (error > 1) return
        else
          register_process_string =  check_process(libhandle, amptype, librarytype, 2)
          if (error > 1) return
        end if
        if (register_process_string > 0) then ! found & registered
          exit
        else if (register_process_string == 0) then ! look in next library type
          librarytype = librarytype + 1
        else
          call ol_msg("register_process: library " // trim(libhandle) // " not found!")
          exit
        end if
      end do

    end if

    !deallocate
    if (allocated(ext)) deallocate(ext)
    if (allocated(perm)) deallocate(perm)
    if (allocated(pol)) deallocate(pol)
    if (allocated(photon_id)) deallocate(photon_id)
  contains

  subroutine charge_conj(x)
    ! determine charge conjugate of array x(3:)
    implicit none
    type(extparticle), intent(inout)     :: x(:)
    integer :: i

    do i=1,size(x)
      if (.not. x(i)%is_initial) then
        select case(x(i)%id)
        case(0, 21, 22, 23, 25, 35, 36)
          x(i)%id  = x(i)%id
        case default
          x(i)%id = -x(i)%id
        end select
      end if
    end do
  end subroutine charge_conj

  function nextpid(ext,pids)
      use KIND_TYPES, only: DREALKIND
      implicit none
      type(extparticle), intent(in) :: ext(:)
      integer, intent(in) :: pids(:)
      integer nextpid, i, j
      nextpid=0
      do i=1,size(ext)
        do j=1,size(pids)
          if (abs(ext(i)%id) == abs(pids(j))) nextpid=nextpid+1
        end do
      end do
  end function nextpid

  subroutine check_photon_id(ext,photon_id)
    implicit none
    type(extparticle), intent(inout) :: ext(:)
    integer, intent(out) :: photon_id(:)
    integer i
    photon_id=0
    do i=1,size(ext)
      if (abs(ext(i)%id) == 2002) then
        if (ext(i)%id == 2002) then ! on-shell photon
          if (ext(i)%is_initial) then
            photon_id(i)=1  ! on-shell initial-state photon
          else
            photon_id(i)=2  ! on-shell final-state photon
          end if
        else
          if (ext(i)%is_initial) then
            photon_id(i)=-1  ! off-shell initial-state photon
          else
            photon_id(i)=-2  ! off-shell final-state photon
          end if
        end if
        ext(i)%id = 22
      end if
    end do
  end subroutine check_photon_id

  subroutine normal_order(ext, perm, extid, proc)
      use KIND_TYPES, only: DREALKIND
      implicit none
      type(extparticle), intent(in) :: ext(:)
      integer, intent(out) :: perm(:)
      integer, intent(out) :: extid(:)
      character(len=*), intent(inout) :: proc
      integer :: i,j, normal(35), pos
      character(len=3) :: normalc(35)
      logical :: do_normal_order = .true.

      ! define normal ordering and corresponding characters

      ! SM
      normal(  1:10) = [ 12  ,-12  , 14  ,-14  , 16  ,-16  , 11  ,-11  , 13  ,-13  ]
      normalc( 1:10) = ["ne ","nex","nm ","nmx","nl ","nlx","e  ","ex ","m  ","mx "]

      normal( 11:20) = [ 15  ,-15  ,  2  , -2  ,  4  , -4  ,  6  , -6  ,  1  , -1  ]
      normalc(11:20) = ["l  ","lx ","u  ","ux ","c  ","cx ","t  ","tx ","d  ","dx "]

      normal( 21:30) = [  3  , -3  ,  5  , -5  , 25  , 35  , 36  , 37  ,-37  , 22  ]
      normalc(21:30) = ["s  ","sx ","b  ","bx ","h  ","h0 ","a0 ","hp ","hpx","a  "]

      normal( 31:35) = [ 23  , -24 , 24  , 21  ,  0  ]
      normalc(31:35) = ["z  ","w  ","wx ","g  ","g  "]

      perm = 0

      if(do_normal_order) then
        ! normal order, build string and store permutation
        pos = 1
        do i = 1, size(normal)
          do j = 1, size(ext)
            if (ext(j)%id == normal(i)) then
              proc = trim(proc) // trim(normalc(i))
              perm(j) = pos
              extid(j) = ext(j)%id
              pos = pos + 1
            end if
          end do
        end do
        if (pos-1 /= size(ext)) then
          proc = ""
        end if
      else
        ! build string, no normal ordering, no permutation
        do j = 1, size(ext)
          do i = 1, size(normal)
            if (ext(j)%id == normal(i)) then
              proc = trim(proc) // trim(normalc(i))
              cycle
            end if
          end do
          extid(j) = ext(j)%id
          perm(j) = j
        end do
      end if

  end subroutine normal_order

  function loop_over_libraries(proc, amptype, n_in, perm, pol, extid, process_in, photon_id)
  use ol_parameters_decl_/**/DREALKIND, only: check_collection
  ! loop over library types
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/DREALKIND, only: coupling_qcd, coupling_ew, OLMode
    implicit none
    character(len=max_parameter_length), intent(in) :: proc
    integer, intent(in) :: amptype, n_in
    integer, intent(in), optional :: perm(:)
    integer, intent(in), optional :: pol(:)
    integer, intent(in), optional :: extid(:)
    character(len=*), intent(in), optional :: process_in
    integer, intent(in), optional :: photon_id(:)
    integer loop_over_libraries
    integer librarytype, loop_olmode, check

    loop_over_libraries = -1
    librarytype = 0
    LoopLibrarytype: do
      if (OLMode == -1) then
        LoopOLMode: do loop_olmode = 3, 0, -1
          check = check_process(proc, amptype, librarytype, n_in, &
                  olmode=loop_olmode, perm_in=perm, pol=pol, extid=extid, process_string=process_in, photon_id=photon_id)
          if (error > 1) return
          if (check > 0) exit LoopOLMode
        end do LoopOLMode
      else
        check = check_process(proc, amptype, librarytype, n_in, &
                olmode=OLMode, perm_in=perm, pol=pol, extid=extid, process_string=process_in, photon_id=photon_id)
        if (error > 1) return
      end if
      if (check > 0) then ! found & registered
        loop_over_libraries = check
        exit Looplibrarytype
      else if (check == 0) then ! look in next library type
        librarytype = librarytype + 1
      else
        exit Looplibrarytype
      end if
    end do Looplibrarytype

    if (check == -1 .and. check_collection) then  ! not found --> check collections
      check = check_process(proc, 999, librarytype, n_in, perm_in=perm, pol=pol, extid=extid, process_string=process_in)
      if (error > 1) return
    end if

  end function loop_over_libraries

  end function register_process_string


  function register_process_id(ext, amptype, n_in_in)
    ! ext: array with format [in_1, .. , in_n_in, out_1, .. , out_n_out]
    ! amptype: integer 1,2,3,4,11,12
    ! (optional) n_in_in: number of initial state particles, default=2
    ! return (integer) process id to be used in evaluate_process
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: ext(:)
    integer, intent(in) :: amptype
    integer, optional, intent(in) :: n_in_in
    integer :: register_process_id
    character(len=max_parameter_length) :: process
    integer :: n_in, i
    if (present(n_in_in)) then
      n_in = n_in_in
    else
      n_in = 2
    end if
    process = to_string(ext(1:n_in),.false., " ")
    process = trim(process) //  " -> " // to_string(ext(n_in+1:),.false., " ")
    register_process_id = register_process_string(process, amptype)
  end function register_process_id


  function check_process(proc_in, amptype, librarytype, n_in, olmode, perm_in, pol, extid, process_string, photon_id)
  ! 1: found, 0: not found, -1: abort
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/DREALKIND, only: &
      & install_path, rMB, rMC, &
      & allowed_libs, tmp_dir, auto_preset, apicheck
    use ol_generic, only: to_string, to_lowercase, count_substring, string_to_integerlist
    use ol_loop_parameters_decl_/**/DREALKIND, only: stability_mode,a_switch,a_switch_rescue,redlib_qp, &
    use_bubble_vertex
    use ol_version, only: process_api
    implicit none
    character(len=max_parameter_length), intent(in) :: proc_in
    integer, intent(in) :: amptype, librarytype, n_in
    integer, intent(in), optional :: olmode
    integer, intent(in), optional :: perm_in(:)
    integer, intent(in), optional :: pol(:)
    integer, intent(in), optional :: extid(:)
    character(len=*), intent(in), optional :: process_string
    integer, intent(in), optional :: photon_id(:)
    integer, allocatable :: perm(:)
    integer, allocatable :: select_pol(:)
    integer check_process
    integer :: lib_content
    integer, save :: info_files_read = 0
    integer :: readok, ierrg
    integer :: i, j, p, p_unmapped
    integer, save :: max_out_length = 35
    logical :: found, found_condmap
    logical :: is_already_loaded, only_loaded
    character(len=4) :: loops_specification
    character(len=4) :: lib_specification
    character(len=max_parameter_length) :: proc, libfilename, libhandle, libname
    character(len=max_parameter_length) :: map_libname
    character(len=max_parameter_length) :: procunmapped, proclastmapped
    character(len=max_parameter_length) :: mapping_str
    character(len=max_parameter_length) :: outstring

    check_process = -1
    found = .false.
    check_process = 0
    map_libname = ' '
    mapping_str = ' '
    if (present(perm_in)) then
      allocate(perm(size(perm_in)))
      perm = perm_in
    end if

    if (.not. check_proclib_exists()) then
      check_process = -2
      return
    end if

    ! read all info files
    if (info_files_read < 1) then
      call readAllInfoFiles()
      if (error /= 0)  then
        check_process = -2
        return
      end if
      info_files_read = 1
    end if

    only_loaded = .false.
    ! set loops_specification
    select case (amptype)
      case (1,2,3,4) ! tree-like
        loops_specification = "t"
        if (librarytype == 0) then
          lib_specification = "lt"
          only_loaded = .true.
        else if (librarytype == 1) then
          lib_specification = "t"
        else if (librarytype == 2) then
          lib_specification = "lt"
        else if (librarytype == 3) then
          lib_specification = "lpt"
        else if (librarytype == 4) then
          lib_specification = "lst"
        else if (librarytype == 5) then
          lib_specification = "lpst"
        else
          check_process = -1
          return
        end if
      case (11) ! loop
        loops_specification = "l"
        if (librarytype == 0) then
          lib_specification = "lt"
        else if (librarytype == 1) then
          lib_specification = "l"
        else if (librarytype == 2) then
          lib_specification = "lp"
        else if (librarytype == 3) then
          lib_specification = "lst"
        else if (librarytype == 4) then
          lib_specification = "lpt"
        else if (librarytype == 5) then
          lib_specification = "lpst"
        else
          check_process = -1
          return
        end if
      case (12) ! loop-induced
        loops_specification = "s"
        if (librarytype == 0) then
          lib_specification = "ls"
        else if (librarytype == 1) then
          lib_specification = "lst"
        else if (librarytype == 2) then
          lib_specification = "lps"
        else if (librarytype == 3) then
          lib_specification = "lpst"
        else
          check_process = -1
          return
        end if
      case (999) ! check libraries
        lib_specification = "lib"
        loops_specification = ""
        if (info_files_read < 2) then
          call readAllInfoFiles(.true.)
            if (error /= 0)  then
              check_process = -2
              call ol_msg("Error: no process libraries installed.")
              return
            end if
            info_files_read = 2
        end if
      case default
        call ol_msg("register_process: amplitude type not supported: " // to_string(amptype))
        check_process = -2
    end select

    ! find process
      proc = proc_in
      p = 0
      p_unmapped = 0
      procunmapped = proc
      InfoLoop: do
        p = p+1
        if (p > size(process_infos)) then
          if ( len_trim(map_libname) /= 0 ) then
            map_libname = ' '
            mapping_str = ' '
            proc = procunmapped
            p = p_unmapped
            if (allocated(perm)) then
              perm = perm_in
            end if
            cycle
          else
            exit
          end if
        end if

         !process loader
        if (index(proc_in,"_") == 0) then

          libname = trim(process_infos(p)%LIBNAME)

          !check if library is "allowed"
          if (len_trim(allowed_libs) /= 0 &
            .and. index(allowed_libs, " " // trim(libname) // " ") == 0) cycle InfoLoop

          !check if correct mappings
          if (len_trim(map_libname) /= 0 &
            .and. trim(map_libname) /= trim(libname)) cycle InfoLoop
          p_unmapped = p

          ! check correct process
          if ( trim(proc) /= trim(process_infos(p)%PROC) &
            & .or. trim(lib_specification) /= trim(process_infos(p)%LTYPE) &
            & .or. index(trim(process_infos(p)%TYPE), trim(loops_specification)) == 0 &
            & ) cycle InfoLoop

          !follow mapping
          if(len_trim(process_infos(p)%MAP) /= 0) then
            !check for conditional mappings
            call check_parameters_condmap(p, found_condmap)
            if (found_condmap) then
              proclastmapped = proc
              proc = trim(process_infos(p)%MAP)
              call ol_msg(2, "Following info-file mapping: " // trim(proclastmapped) // &
                          " --> " // trim(process_infos(p)%MAP) // "[" // trim(process_infos(p)%MAPPERM) // "].")
              !map permutation
              if(len_trim(process_infos(p)%MAPPERM) /= 0) then
                if (allocated(perm)) then
                  call map_permutation(perm,string_to_integerlist(process_infos(p)%MAPPERM))
                  if (error > 1) return
                end if
              end if
              mapping_str = " (mapped from " // trim(procunmapped) // ")"
              map_libname = libname
              p = 0
              cycle InfoLoop
            else
              call ol_msg(2, "Not following mapping " // trim(proc) // " --> " // trim(process_infos(p)%MAP) // ".")
              cycle InfoLoop
            end if
          end if

          ! get library filename
          if (amptype == 999) then
            libfilename = "collection " // trim(process_infos(p)%LIBNAME)
          else
            libfilename = 'libopenloops_' // trim(process_infos(p)%LIBNAME) // '_' // &
                          & trim(process_infos(p)%LTYPE) // '.' // dynlib_extension
          end if

          if (allocated(select_pol)) deallocate(select_pol)
          allocate(select_pol(size(pol)))
          select_pol = pol

          ! check OLmode
          if (present(OLMode)) then
            if (process_infos(p)%OLMode /= OLMode) cycle InfoLoop
          end if

          ! check parameters
          call check_parameters(p,amptype,found)

        ! direct library loader
        else if (index(proc_in,"_") > 0) then
          libhandle = proc_in
          libname = proc_in(:index(proc_in(:index(proc_in,"_",.true.)-1),"_",.true.)-1)
          proc = proc_in(index(proc_in(:index(proc_in,"_",.true.)-1),"_",.true.)+1:index(proc_in,"_",.true.)-1)
          if ( trim(proc) /= trim(process_infos(p)%PROC) .or. &
            &  trim(libname) /= trim(process_infos(p)%LIBNAME) .or. &
            &  trim(lib_specification) /= trim(process_infos(p)%LTYPE) .or. &
            &  index(trim(process_infos(p)%TYPE), trim(loops_specification)) == 0 .or. &
            &  proc_in(index(proc_in,"_",.true.)+1:) /=  trim(process_infos(p)%ID) &
            & ) cycle InfoLoop
          if (apicheck .and. process_api /= process_infos(p)%API) cycle InfoLoop
          libfilename = 'libopenloops_' // trim(libname) // '_' // &
                           & trim(lib_specification) // '.' // dynlib_extension
          found = .true.
          exit InfoLoop
        else
          call ol_error("register_process: process format not supported.")
          check_process = -2
          return
        end if


        !if required, check if library is already loaded
        if (only_loaded) then
          is_already_loaded = .false.
          if (allocated(loaded_libs)) then
            do j = 1, size(loaded_libs)
              if (trim(loaded_libs(j)%LIBNAME) == trim(libname)  &
                & .and. index(trim(loaded_libs(j)%TYPE), trim(loops_specification)) > 0 &
                &  ) then
                  is_already_loaded = .true.
              end if
            end do
          end if
          if (.not. is_already_loaded) found = .false.
        end if

        ! found correct library
        if (found) then
          call ol_msg(2, "Parameters do match info-file for process " // trim(proc) // " in library " // trim(libfilename))
          if (amptype == 999) then
            if (present(process_string)) then
              call ol_msg("Library for " // trim(process_string) // " not installed but available in: " // trim(libname))
            else
              call ol_msg("Library for " // trim(proc) // " not installed but available in: " // trim(libname))
            end if
            call ol_msg("Note: this library can be downloaded and installed via")
            call ol_msg("$ cd " // trim(install_path))
            call ol_msg("$ ./openloops libinstall " // trim(libname))
            check_process = 1
            return
          end if
          libhandle = trim(to_lowercase(libname)) // "_" // trim(proc) // "_" // trim(process_infos(p)%ID)
          exit
        else
          call ol_msg(2,"Parameters do not match info-file for process " //  &
                      trim(proc) // "_" // trim(process_infos(p)%ID) // " in library " // trim(libname))
          check_process = 0
        end if

      end do InfoLoop

    if (found) then
        libfilename = trim(install_path) // '/proclib/' // libfilename
        lib_content = 0
        do i = 1, len(loops_flags)
          if (index(trim(process_infos(p)%TYPE), loops_flags(i:i)) > 0) lib_content = ibset(lib_content, i-1)
        end do

        if(process_infos(p)%OLMode .le. 0) call set_if_modified(use_bubble_vertex, 0)

        if(auto_preset) then
          if (process_infos(p)%OLMode .ge. 2) then
            call set_if_modified(a_switch,1)
            call set_if_modified(stability_mode,11)
          else
            if (amptype == 12) then
              call set_if_modified(a_switch,1)
              call set_if_modified(a_switch_rescue,7)
              call set_if_modified(stability_mode,21)
            else
               if(process_infos(p)%MODEL=="heft") then
                 call set_if_modified(a_switch,5)
                 call set_if_modified(redlib_qp,5)
                 call set_if_modified(stability_mode,14)
               else
                 call set_if_modified(a_switch,1)
                 call set_if_modified(a_switch_rescue,7)
                 call set_if_modified(redlib_qp,5)
                 call set_if_modified(stability_mode,23)
               end if
            end if
          end if
        end if

        !register
        check_process = register_process_lib(libfilename, libhandle, lib_content, amptype, &
                              & n_in, perm=perm, pol=select_pol, extid=extid, photon_id=photon_id, &
                              & qcd_powers=process_infos(p)%qcdorder)
        if (error > 1) then
          call ol_error("register_process_lib failed")
          check_process = -2
          return
        end if

        if (present(process_string)) then
          outstring = "Library loaded: " //  trim(process_string)
        else
          outstring = "Library loaded: " //  trim(proc)
        end if
        outstring = adjustl(outstring)
        if (len_trim(outstring) > max_out_length) max_out_length = len_trim(outstring)
        outstring = outstring(1:max_out_length) // " @" // &
          & " tree=" // trim(to_string(process_infos(p)%EWorder(0)))  // "," // &
          & trim(to_string(process_infos(p)%QCDorder(0)))
        if (amptype > 10) then
          outstring = trim(outstring) // &
          & " loop=" // trim(to_string(process_infos(p)%EWorder(0)+process_infos(p)%EWorder(1))) // "," // &
          & trim(to_string(process_infos(p)%QCDorder(0)+process_infos(p)%QCDorder(1)))
        end if
          outstring = trim(outstring) // " (EW,QCD) >  " // trim(libhandle)
        if (allocated(perm)) then
          outstring = trim(outstring) //  trim(to_string(perm,.true.))
        end if
        outstring = trim(outstring) // " (id=" // trim(to_string(check_process)) // ")"
        outstring = trim(outstring) // trim(mapping_str)
        call ol_msg(1,outstring)

        !add to list of loaded libraries
        call add_loaded_library(process_infos(p))

    end if


    if (allocated(perm)) deallocate(perm)
    if (allocated(select_pol)) deallocate(select_pol)

  contains

    subroutine map_permutation(perm, map)
    !map permutation
      implicit none
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: map(:)
      integer  :: perm_tmp(size(perm))
      integer :: i, x
      if (size(perm) /= size(map)) then
        call ol_fatal("error in map_permutation")
        return
      end if
      do i = 1, size(map)
        perm_tmp(i) = map(perm(i))
      end do
      perm = perm_tmp
    end subroutine map_permutation



  end function check_process



  function check_proclib_exists()
  ! checks that proclib folder exists within set install_path
    use ol_parameters_decl_/**/DREALKIND, only: install_path
    implicit none
    logical check_proclib_exists
    logical proclib_exists
#ifdef USE_GFORTRAN
    inquire(file=trim(install_path)//"/proclib/.", exist=proclib_exists)
#endif
#ifdef USE_IFORT
    inquire(directory=trim(install_path)//"/proclib", exist=proclib_exists)
#endif
    if (.not. proclib_exists) then
      call ol_fatal("register_process: proclib folder not found, check install_path or install libraries.")
      check_proclib_exists = .false.
      return
    else
      check_proclib_exists = .true.
      return
    end if
  end function check_proclib_exists


  subroutine check_parameters(p,amptype,found)
    use ol_parameters_decl_/**/DREALKIND, only: &
      & rME, rMM, rML, rMU, rMD, rMC, rMS, rMB, rMT, &
      & rYE, rYM, rYL, rYU, rYD, rYC, rYS, rYB, rYT, &
      & leadingcolour, coupling_qcd, coupling_ew, &
      & loop_order_ew, loop_order_qcd, &
      & approximation, CKMORDER, model, allowed_libs, apicheck
    use ol_loop_parameters_decl_/**/DREALKIND , only: nf, nc
    use ol_version, only: process_api
    implicit none
    integer, intent(in) :: p
    integer, intent(in) :: amptype
    logical, intent(out) :: found

    found = .false.
    if (allocated(process_infos)) then
      if (size(process_infos) < p) then
        call ol_error(1,"check_parameters: process not available")
        return
      end if
      found = .true.

      if (apicheck) call check(process_infos(p)%API == process_api, found, "API version not ok.")
      if (trim(model) /= "heft") then
        call check(process_infos(p)%EWorder(0) == coupling_EW(0) .or. coupling_EW(0) == -1, found, "EW tree coupling NOT ok.")
        call check(amptype == 1 .or. process_infos(p)%eworder(1) == coupling_ew(1) &
                  .or. coupling_ew(1) == -1, found, "EW loop not ok.")
      end if
      call check(process_infos(p)%qcdorder(0) == coupling_qcd(0) .or. coupling_qcd(0) == -1, found, "qcd tree coupling not ok.")
      call check(amptype == 1  .or. process_infos(p)%qcdorder(1) == coupling_qcd(1) &
                  .or. coupling_qcd(1) == -1, found, "qcd loop not ok.")
      call check(loop_order_ew == -1 .or. process_infos(p)%eworder(0)+process_infos(p)%eworder(1) == loop_order_ew, &
                     found, "absolute ew loop not ok.")
       call check(loop_order_qcd == -1 .or. process_infos(p)%qcdorder(0)+process_infos(p)%qcdorder(1) == loop_order_qcd, &
                     found, "absolute qcd loop not ok.")
      call check(process_infos(p)%LeadingColour == leadingcolour, found, "LeadingColour OK.")
      call check(process_infos(p)%NC == nc, found, "nc NOT ok.")
      call check(process_infos(p)%NF == nf, found, "nf NOT ok.")
      call check(process_infos(p)%CKMorder == CKMORDER, found, "CKM NOT ok.")
      call check(index(process_infos(p)%MODEL, trim(model)) == 1, found, "model NOT ok.")
      call check((process_infos(p)%ME /= "0" .and. rME /= 0) .or. rME == 0, found, "mass ME NOT ok.")
      call check((process_infos(p)%MM /= "0" .and. rMM /= 0) .or. rMM == 0, found, "mass MM NOT ok.")
      call check((process_infos(p)%ML /= "0" .and. rML /= 0) .or. rML == 0, found, "mass ML NOT ok.")
      call check((process_infos(p)%MU /= "0" .and. rMU /= 0) .or. rMU == 0, found, "mass MU NOT ok.")
      call check((process_infos(p)%MD /= "0" .and. rMD /= 0) .or. rMD == 0, found, "mass MD NOT ok.")
      call check((process_infos(p)%MS /= "0" .and. rMS /= 0) .or. rMS == 0, found, "mass MS NOT ok.")
      call check((process_infos(p)%MC /= "0" .and. rMC /= 0) .or. rMC == 0, found, "mass MC NOT ok.")
      call check((process_infos(p)%MB /= "0" .and. rMB /= 0) .or. rMB == 0, found, "mass MB NOT ok.")
      call check(rME == rYE  .or. process_infos(p)%YE == 1, found, "YukE /= ME NOT ok.")
      call check(rMM == rYM  .or. process_infos(p)%YM == 1, found, "YukM /= MM NOT ok.")
      call check(rML == rYL  .or. process_infos(p)%YL == 1, found, "YukL /= ML NOT ok.")
      call check(rMU == rYU  .or. process_infos(p)%YU == 1, found, "YukU /= MU NOT ok.")
      call check(rMD == rYD  .or. process_infos(p)%YD == 1, found, "YukD /= MD NOT ok.")
      call check(rMS == rYS  .or. process_infos(p)%YS == 1, found, "YukS /= MS NOT ok.")
      call check(rMC == rYC  .or. process_infos(p)%YC == 1, found, "YukC /= MC NOT ok.")
      call check(rMB == rYB  .or. process_infos(p)%YB == 1, found, "YukB /= MB NOT ok.")
      call check(rMT == rYT  .or. process_infos(p)%YT == 1, found, "YukT /= YT NOT ok.")
      call check(amptype == 1 .or. amptype > 10 .or. process_infos(p)%CC /= "0", found, "CC NOT ok.")
      call check(trim(process_infos(p)%APPROX) == trim(approximation), found, "APPROX NOT ok.")
    end if

  end subroutine check_parameters

  subroutine check_parameters_condmap(p,found)
    use ol_parameters_decl_/**/DREALKIND, only: &
      & rME, rMM, rML, rMU, rMD, rMC, rMS, rMB
    implicit none
    integer, intent(in) :: p
    logical, intent(out) :: found

    found = .false.
    if (allocated(process_infos)) then
      if (size(process_infos) < p) then
        call ol_error(1,"check_parameters_mapping: process not available")
        return
      end if
      found = .true.
      call check(.not. (process_infos(p)%ME == "0" .and. rME /= 0), found, "mass ME NOT ok.")
      call check(.not. (process_infos(p)%MM == "0" .and. rMM /= 0), found, "mass MM NOT ok.")
      call check(.not. (process_infos(p)%ML == "0" .and. rML /= 0), found, "mass ML NOT ok.")
      call check(.not. (process_infos(p)%MU == "0" .and. rMU /= 0), found, "mass MU NOT ok.")
      call check(.not. (process_infos(p)%MD == "0" .and. rMD /= 0), found, "mass MD NOT ok.")
      call check(.not. (process_infos(p)%MS == "0" .and. rMS /= 0), found, "mass MS NOT ok.")
      call check(.not. (process_infos(p)%MC == "0" .and. rMC /= 0), found, "mass MC NOT ok.")
      call check(.not. (process_infos(p)%MB == "0" .and. rMB /= 0), found, "mass MB NOT ok.")
    end if

  end subroutine check_parameters_condmap

  subroutine check(test,found,message)
    implicit none
    logical, intent(in) :: test
    logical, intent(inout) :: found
    character(len=*), intent(in) :: message
    if (.not. test) then
      found = .false.
      call ol_msg(3,"Library does not match: " // trim(message))
    end if
  end subroutine

  subroutine readAllInfoFiles(load_channel_lib)
    use ol_parameters_decl_/**/DREALKIND, only: install_path
    use iso_fortran_env, only: iostat_end
    use ol_cwrappers, only: opendir, readdir, closedir
    implicit none
    logical, optional, intent(in) :: load_channel_lib
    integer :: readok
    integer, parameter :: gf_info = 994
    integer :: counter,API
    character(len=500) :: infofilename
    character(len=500) :: infoline
    character(len=5) :: info_file_suffix = 'info'
    logical :: iqopen
    type(processinfos) infos
    type(processinfos), allocatable :: process_infos_bak(:)

    if (present(load_channel_lib)) then
      call ol_msg(1, "Requested library not installed. Checking collection...")
      if (load_channel_lib) then
        info_file_suffix = 'rinfo'
      end if
    end if

    ! open proclib folder
    readok = opendir(trim(install_path) // '/proclib')
    if (readok /= 0) then
      call ol_error('opening proclib directory failed. Check install_path.')
      return
    end if

    ProclibDirLoop: do
      readok = readdir(infofilename)
      if (readok /= 0) then
        call ol_error("reading proclib directory content failed.")
        exit
      end if
      if (len(trim(infofilename)) == 0) exit
      if (index(trim(infofilename),"."//trim(info_file_suffix)) == 0) then
        cycle
      else
        infofilename = trim(install_path) // "/proclib/" // trim(infofilename)
      end if

      inquire(gf_info, opened=iqopen)
      if(iqopen) close(unit=gf_info)
      open(gf_info, file=trim(infofilename), status = "old", iostat=readok)
      if (readok /= 0) then
        call ol_error("in readAllInfoFiles can't open file: " // trim(infofilename) )
        exit
      end if
      counter = 0
      API = 1
      InfoFileLoop: do
        read (gf_info, '(A)', iostat=readok )  infoline
        if (readok /= 0) then ! EOF -> exit
          if (readok == iostat_end) then
            exit InfoFileLoop
          else
            call ol_error("in redAllInfoFiles error reading file: " // trim(infofilename) )
            exit ProclibDirLoop
          end if
        end if

        ! strip empty lines
        if (len_trim(infoline) == 0) then
          cycle InfoFileLoop
        end if

        infoline = adjustl(infoline)
        ! strip possible comment: start with #
        if (infoline(1:1) == "#") then
          cycle InfoFileLoop
        end if

        ! strip lines starting with OPTIONS=... (deprecated)
        if (infoline(1:8) == "options ") then
          call readInfoInt(infoline, 'API', API)
          cycle InfoFileLoop
        end if

        counter = counter+1
        ! strip first line of collection files
        if (info_file_suffix == 'rinfo' .and. counter == 1) then
          cycle InfoFileLoop
        end if

        ! determine library type from name of info file
        if (info_file_suffix == 'rinfo') then
          infos%LTYPE = "lib"
        else
          infos%LTYPE = infofilename(index(infofilename,'_',.true.)+1:index(infofilename,'.info')-1)
        end if

        call readAllInfos(infoline, infos, API)
        if (error > 1) then
          call ol_error("reading infofile line: " // trim(infoline))
          exit ProclibDirLoop
        end if

        ! Add to array of infos
        if (.not. allocated(process_infos)) then
          allocate(process_infos(1))
        else
          allocate(process_infos_bak(size(process_infos)))
          process_infos_bak = process_infos
          deallocate(process_infos)
          allocate(process_infos(size(process_infos_bak)+1))
          process_infos(1:size(process_infos_bak)) = process_infos_bak
          deallocate(process_infos_bak)
        end if
        process_infos(size(process_infos)) = infos
      end do InfoFileLoop
    end do ProclibDirLoop

    ! close directory handle
    call closedir()

    ! close file handle
    inquire(gf_info, opened=iqopen)
    if(iqopen) close(unit=gf_info)

    if (.not. allocated(process_infos)) then
      call ol_error("no processes installed!")
    end if

    contains

    subroutine readAllInfos(lineinfo, infos, api)
      implicit none
      character(len=*), intent(in) :: lineinfo
      integer,          intent(in) :: api
      type(processinfos), intent(inout) :: infos
      integer :: ccount

      call readInfoCol(lineinfo, 1, infos%LIBNAME)
      call readInfoCol(lineinfo, 2, infos%PROC)

      if (index(lineinfo, 'map') > 0) then
        if (index(lineinfo, ' map') > 0) then
          call readInfo(lineinfo, 'map', infos%MAP)
        else if (index(lineinfo, ' condmap') > 0) then
          call readInfo(lineinfo, 'condmap', infos%MAP)
        else
          call ol_fatal("info-file mapping not supported!")
          return
        end if
        infos%MAPPERM =  infos%MAP(index(infos%MAP, '[')+1:index(infos%MAP, ']')-1)
        if (len_trim(infos%MAPPERM) /= 0) then
          infos%MAP = infos%MAP(1:index(infos%MAP,'[')-1)
        end if
        infos%ID = '0'
      else
        call readInfoCol(lineinfo, 3, infos%ID)
        call readInfoCoupling(lineinfo, 'QCD', infos%QCDorder)
        call readInfoCoupling(lineinfo, 'EW', infos%EWorder)
        infos%MAP = ' '
        infos%MAPPERM = ' '
      end if
      call readInfo(lineinfo, 'ME', infos%ME)
      call readInfo(lineinfo, 'MM', infos%MM)
      call readInfo(lineinfo, 'ML', infos%ML)
      call readInfo(lineinfo, 'MU', infos%MU)
      call readInfo(lineinfo, 'MD', infos%MD)
      call readInfo(lineinfo, 'MS', infos%MS)
      call readInfo(lineinfo, 'MC', infos%MC)
      call readInfo(lineinfo, 'MB', infos%MB)
      call readInfoInt(lineinfo, 'YukE', infos%YE)
      call readInfoInt(lineinfo, 'YukM', infos%YM)
      call readInfoInt(lineinfo, 'YukL', infos%YL)
      call readInfoInt(lineinfo, 'YukU', infos%YU)
      call readInfoInt(lineinfo, 'YukD', infos%YD)
      call readInfoInt(lineinfo, 'YukS', infos%YS)
      call readInfoInt(lineinfo, 'YukC', infos%YC)
      call readInfoInt(lineinfo, 'YukB', infos%YB)
      call readInfoInt(lineinfo, 'YukT', infos%YT)
      call readInfo(lineinfo, 'APPROX', infos%APPROX)
      call readInfoInt(lineinfo, 'CKMORDER', infos%CKMorder)
      call readInfoInt(lineinfo, 'nc', infos%NC)
      call readInfoInt(lineinfo, 'nf', infos%NF)
      call readInfoInt(lineinfo, 'LeadingColour', infos%LeadingColour)
      call readInfoInt(lineinfo, 'POLSEL', infos%POLSEL)
      call readInfo(lineinfo, 'CC', infos%CC)
      call readInfo(lineinfo, 'MODEL', infos%Model)
      call readInfoInt(lineinfo, 'OLMode', infos%OLMode)
      infos%API = api
      if (len_trim(infos%Model) == 0) then
        infos%Model = "sm"
      end if
      call readInfo(lineinfo, 'Type', infos%TYPE)
      if (trim(infos%TYPE) == "") then
        infos%TYPE = infos%LTYPE
      end if
    end subroutine readAllinfos

    subroutine readInfo(lineinfo, var, res)
      implicit none
      character(len=*), intent(in) :: lineinfo
      character(len=*), intent(in) :: var
      character(len=*), intent(out) :: res
      if (index(lineinfo, ' '//var//'=') /= 0) then
        res = lineinfo(index(lineinfo, var//'=')+len_trim(var)+1: &
              & index(lineinfo, var//'=')+index(lineinfo(index(lineinfo, var//'='):),' ')-1 )
      else
        res = ""
      end if
    end subroutine readInfo

    subroutine readInfoCol(lineinfo, col, res)
      implicit none
      character(len=*), intent(in) :: lineinfo
      integer, intent(in)  :: col
      character(len=*), intent(out) :: res
      integer sstart, send, i
      sstart = 1
      do i=1,col-1
        sstart = sstart + index(lineinfo(sstart:), " ")
      end do
      send = sstart + index(lineinfo(sstart:), " ") - 2
      res = trim(lineinfo(sstart:send))
    end subroutine readInfoCol

    subroutine readInfoColInt(lineinfo, col, res)
      use ol_generic, only: to_int
      implicit none
      character(len=*), intent(in) :: lineinfo
      integer, intent(in)  :: col
      integer, intent(out) :: res
      character(len=max_parameter_length) :: resc
      integer sstart, send
      call readInfoCol(lineinfo,col, resc)
      res = to_int(trim(resc))
      if (res == -huge(res)) then
        call ol_msg(1, "Warning: problem reading info line: " // trim(lineinfo))
      end if
    end subroutine readInfoColInt

    subroutine readInfoInt(lineinfo, var, res)
      use ol_generic, only: to_int
      implicit none
      character(len=*), intent(in) :: lineinfo
      character(len=*), intent(in) :: var
      character(len=10) :: restemp
      integer, intent(out) :: res
      if (index(lineinfo, var//'=') /= 0) then
        restemp = lineinfo(index(lineinfo, var//'=')+len_trim(var)+1: &
            & index(lineinfo, var//'=')+index(lineinfo(index(lineinfo, var//'='):),' ')-1)
      else
        restemp = "0"
      end if
      res = to_int(restemp)
      if (res == -huge(res)) then
        call ol_msg(1, "Warning: problem reading info line: " // trim(lineinfo))
      end if
    end subroutine readInfoInt

    subroutine readInfoCoupling(lineinfo, var, res)
      use ol_generic, only: to_int
      implicit none
      character(len=*), intent(in) :: lineinfo
      character(len=*), intent(in) :: var
      integer ::  stat(0:1)
      character(len=max_parameter_length) :: restempc
      integer, intent(out) :: res(0:1)
      if (index(lineinfo, var//'=') /= 0) then
        call readInfo(lineinfo,var,restempc)
        res(0) = to_int(trim(restempc(1:index(restempc,",")-1)))
        res(1) = to_int(trim(restempc(index(restempc,",")+1:)))
      else
        res(0) = 0
        res(1) = 0
      end if
      if (any(res == -huge(res))) then
        call ol_msg(1,"Warning: problem reading info line: " // trim(lineinfo))
      end if
    end subroutine readInfoCoupling

  end subroutine readAllInfoFiles


   subroutine add_loaded_library(infos)
   ! add infos of loaded library to the list of already loaded libraries.
    implicit none
    type(processinfos), intent(in) :: infos
    type(processinfos), allocatable :: loaded_libs_bak(:)
    integer :: i,loaded_libs_last
    if (allocated(loaded_libs)) then
      do i = 1, size(loaded_libs)
        if ( trim(infos%LIBNAME) == trim(loaded_libs(i)%LIBNAME) .and.  &
          &   trim(infos%TYPE) == trim(loaded_libs(i)%TYPE) ) return
      end do
      loaded_libs_last = size(loaded_libs)
      allocate(loaded_libs_bak(loaded_libs_last))
      loaded_libs_bak = loaded_libs
      deallocate(loaded_libs)
      allocate(loaded_libs(loaded_libs_last+1))
      loaded_libs(1:loaded_libs_last) = loaded_libs_bak
      deallocate(loaded_libs_bak)
    else
     loaded_libs_last = 0
     allocate(loaded_libs(1))
    end if
    loaded_libs(loaded_libs_last+1) = infos
   end subroutine


  subroutine flavour_mapping(ext,switch_in)
  ! ext: PDG coded integer array
  ! switch: 0: no mapping, 1: quark & lepton mapping, 2: only lepton mapping, 3: only quark mapping
  ! ==Lepton & quark flavour mapping==
  !   Concept (taken from old Sherpa interface):
  !   (1) Given a final state, determine the four (anti)lepton/neutrino
  !       multiplicities in the given process:
  !       a_nubar, a_nu, a_lbar, a_l
  !   (2) Compute the discriminant N[i] as
  !       N[i] = Ngen - i + Ngen*(a_nubar + a_nu*Nmax + a_lbar*Nmax^2 + a_l*Nmax^3)
  !       where Nmax should be chosen such that a_...<=Nmax.
  !       In practice one can safely set Nmax=10 and it will work for any
  !       process with <= 20 final-state leptons.
  !       It is also convenient to set Ngen=10, although i runs only from 1 to 3.
  !   (3) Reassign the lepton generations with a permutation
  !       p1 -> 1, p2 -> 2, p3 -> 3  such that  N[p1] > N[p2] > N[p3] */
    use ol_generic, only: to_string
    use ol_parameters_decl_/**/DREALKIND, only: rMC, rYC, rMM, rYM, rML, rYL
    implicit none
    type(extparticle), intent(inout) :: ext(:)
    integer, optional, intent(in)    :: switch_in
    type(extparticle), allocatable :: new_ext(:)
    integer :: switch
    integer :: i, j
    integer :: Ngen, Nlgen, Nqgen, Nmax
    integer :: l_gen, nu_gen, l_gen_new, nu_gen_new
    integer :: d_gen, u_gen, d_gen_new, u_gen_new
    integer :: a_nu, a_nubar, a_l, a_lbar
    integer :: a_u, a_ubar, a_d, a_dbar
    integer :: Nl(3,2), Nq(2,2)
    integer :: perm(size(ext))

    if (present(switch_in)) then
      switch = switch_in
    else
      switch = 0
    end if

    if (switch < 0 .or. switch > 3) then
      call ol_error("flavour_mapping: only options for switch=0,1,2,3")
    end if

    Ngen=10
    Nmax=10

    call ol_msg(3,"Flavour mapping. Original (all ingoing) process: " // to_string([(ext(j)%id, j=1,size(ext))]) )

    if (rML == 0 .and. rYL == 0 .and. rMM == 0 .and. rYM == 0) then
      Nlgen = 3
    else if ((rML /= 0. .or. rYL /= 0 ) .and. rMM == 0 .and. rYM == 0) then
      Nlgen = 2
    else
      Nlgen = 1
    end if

    if (rMC == 0 .and. rYC == 0) then
      Nqgen = 2
    else
      Nqgen = 1
    end if

    allocate(new_ext(size(ext)))
    new_ext = ext

    if (switch == 1 .or. switch == 2) then
      !lepton flavour mapping
      do i = 1, Nlgen
        l_gen=9+2*i;
        nu_gen=10+2*i;

        a_nu=count_integer(ext,nu_gen);
        a_nubar=count_integer(ext,-nu_gen);
        a_l=count_integer(ext,l_gen);
        a_lbar=count_integer(ext,-l_gen);

        Nl(i,1)=Ngen-i+Ngen*(a_nubar+a_nu*Nmax+a_lbar*Nmax*Nmax+a_l*Nmax*Nmax*Nmax)
        Nl(i,2)=i
      end do

      call sort_pair(Nl,Nlgen)

      do i = 1, Nlgen
        l_gen=9+2*Nl(i,2);
        l_gen_new=9+2*i;
        nu_gen=10+2*Nl(i,2);
        nu_gen_new=10+2*i;

        do j = 1, size(ext)
          if (abs(ext(j)%id)==nu_gen) new_ext(j)%id=sign(nu_gen_new,ext(j)%id)
          if (abs(ext(j)%id)==l_gen)  new_ext(j)%id=sign(l_gen_new,ext(j)%id)
        end do
      end do
      ext = new_ext
    end if

    if (switch == 1 .or. switch == 3) then
      !quark flavour mapping
      do i = 1, Nqgen
        d_gen=2*i-1;
        u_gen=2*i;

        a_d=count_integer(ext,d_gen);
        a_dbar=count_integer(ext,-d_gen);
        a_u=count_integer(ext,u_gen);
        a_ubar=count_integer(ext,-u_gen);

        Nq(i,1)=Ngen-(i+1) + Ngen*(a_ubar+a_u*Nmax+a_dbar*Nmax*Nmax+a_d*Nmax*Nmax*Nmax)
        Nq(i,2)=i
      end do

      call sort_pair(Nq,Nqgen)

      do i = 1, Nqgen
        d_gen=2*Nq(i,2)-1;
        d_gen_new=2*i-1;
        u_gen=2*Nq(i,2);
        u_gen_new=2*i;

        do j = 1, size(ext)
          if (abs(ext(j)%id)==u_gen) new_ext(j)%id=sign(u_gen_new,ext(j)%id)
          if (abs(ext(j)%id)==d_gen) new_ext(j)%id=sign(d_gen_new,ext(j)%id)
        end do
      end do
      ext = new_ext
    end if

    deallocate(new_ext)

    call ol_msg(3, "Flavour mapping. Mapped (all ingoing) process:   " // to_string([(ext(j)%id, j=1,size(ext))]))

    contains

    subroutine sort_pair(a,n)
    ! Simple insertion sort. Sorting descending on the first component of a 2-tuple of length n
      integer, intent(in) :: n
      integer, intent(inout), dimension(:,:) :: a
      integer :: temp(2)
      integer :: i, j

      do i = 2, n
        j = i - 1
        temp(1) = a(i,1)
        temp(2) = a(i,2)
        do while (j>=1 .and. a(j,1)<temp(1))
          a(j+1,1) = a(j,1)
          a(j+1,2) = a(j,2)
          j = j - 1
          if (j==0) exit
        end do
        a(j+1,1) = temp(1)
        a(j+1,2) = temp(2)
      end do
    end subroutine sort_pair

    function count_integer(list, j)
    ! count frequency of integer j in array ilist
      implicit none
      type(extparticle), intent(in) :: list(:)
      integer, intent(in) :: j
      integer :: count_integer
      integer i
      count_integer = 0
      do i = 1, size(list)
        if (list(i)%id == j) count_integer = count_integer+1
      end do
    end function count_integer

  end subroutine


  function write_shop_list(ext, proc)
    use ol_parameters_decl_/**/DREALKIND, only: shopping_list, order_ew, order_qcd, &
      & rMU, rMD, rMC, rMS, rMB, rME, rMM, rML
    use ol_generic, only: to_string
    implicit none
    type(extparticle), intent(in) :: ext(:)
    character(len=max_parameter_length), intent(in) :: proc
    integer write_shop_list
    character(len=500) :: output
    character(len=max_parameter_length), allocatable :: shopped_processes_bak(:)
    integer, save :: id = 1
    integer ::  readok
    integer :: i, already_shopped
    integer :: oqcd, oew
    logical :: set_masses = .false.
    logical :: iqopen

    write_shop_list = -1

    if ( .not. (shopping_list_open) ) then
      inquire(fh_shopping, opened=iqopen)
      if(iqopen) close(unit=fh_shopping)
      !open shopping list
      open(fh_shopping, file=trim(shopping_list), status = "unknown", position='append', iostat=readok)
      if (readok /= 0) then
        call ol_msg("Error opening shopping list " // trim(shopping_list))
        return
      end if
      ! write header
      write(fh_shopping,'(A)') ""
      if (order_ew /= -1) then
        oqcd = size(ext)-2-order_ew
        write(fh_shopping,'(A)') "SelectCoupling = (Exponent[#, gQCD] == " // trim(to_string(oqcd)) // " + 2 * #2 &);"
        write(fh_shopping,'(A)') "SelectInterference = {eQED -> " // trim(to_string(order_ew*2)) // "};"
        write(fh_shopping,'(A)') "UnitaryGauge = True;"
        write(fh_shopping,'(A)') ""
      else if (order_qcd /= -1) then
        oew = size(ext)-2-order_qcd
        write(fh_shopping,'(A)') "SelectCoupling = (Exponent[#, eQED] == " // trim(to_string(oew)) // " + 2 * #2 " &
          & // " ||  Exponent[#1, eQED] == " // trim(to_string(oew+2))  // " - 2 * #2 &);"
        write(fh_shopping,'(A)') "SelectInterference = {gQCD -> " // trim(to_string(order_qcd*2)) // "};"
        write(fh_shopping,'(A)') "UnitaryGauge = False;"
        write(fh_shopping,'(A)') ""
      end if
      shopping_list_open = .true.
    end if

    if (allocated(shopped_processes)) then
      do i = 1, size(shopped_processes)
        if ( trim(proc) == trim(shopped_processes(i)) ) then
          call ol_msg(2, "Not written to shopping list. Already shopped as process " &
                    & // trim(to_string(i)) // ": " //trim(proc) )
          write_shop_list = i
          return
        end if
      end do
      already_shopped = size(shopped_processes)
      allocate(shopped_processes_bak(already_shopped))
      shopped_processes_bak = shopped_processes
      deallocate(shopped_processes)
      allocate(shopped_processes(already_shopped+1))
      shopped_processes(1:already_shopped) = shopped_processes_bak
      deallocate(shopped_processes_bak)
    else
      already_shopped = 0
      allocate(shopped_processes(1))
    end if
    shopped_processes(already_shopped+1) = trim(proc)

    ! process name and id
    output = "(* " //trim(proc)// " *) AddProcess[FeynArtsProcess -> "
    ! inital state
    output = trim(output) // " {" // trim(PDGtoFA(ext(1)%id)) // ", " // trim(PDGtoFA(ext(2)%id)) //  "} -> {"
    ! final state
    do i=3,size(ext)
      output = trim(output) // trim(PDGtoFA(ext(i)%id))
      if (i /= size(ext)) output = trim(output) // ", "
    end do
    output = trim(output) // "}"

    ! massive fermions
    set_masses = .false.
    if (rMU /= 0 .or. rMD /= 0 .or. rMC /= 0 .or. rMS /= 0 .or. rMB /= 0 .or. &
      & rME /= 0 .or. rMM /= 0 .or. rML /= 0) then
      output = trim(output) // ", SetParameters -> JoinOptions[{"
      if (rMU /= 0) then
        output = trim(output) // "MU -> MU"
        set_masses = .true.
      end if
      if (rMD /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "MD -> MD"
        set_masses = .true.
      end if
      if (rMC /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "MC -> MC"
        set_masses = .true.
      end if
      if (rMS /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "MS -> MS"
        set_masses = .true.
      end if
      if (rMB /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "MB -> MB"
        set_masses = .true.
      end if
      if (rME /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "ME -> ME"
        set_masses = .true.
      end if
      if (rMM /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "MM -> MM"
        set_masses = .true.
      end if
      if (rML /= 0) then
        if (set_masses) output = trim(output) // ","
        output = trim(output) // "ML -> ML"
        set_masses = .true.
      end if
      output = trim(output) // "}]"
    end if

    output = trim(output) // "];"

    call ol_msg(1," Write to shopping list "// trim(shopping_list)  //": " // trim(output))

    ! write process to shopping list
    write(fh_shopping,'(A)') trim(output)

    write_shop_list = id
    id = id+1
  end function write_shop_list



  function PDGtoFA(pdg)
  ! PDG number scheme -> FeynArts naming scheme
    implicit none
    integer, intent(in) :: pdg
    character(len=10) :: PDGtoFA

    if (pdg < 0 .and. pdg /= -24) then
      PDGtoFA = "-"
    else if (pdg == 24) then
      PDGtoFA = "-"
    else
      PDGtoFA = ""
    end if

    select case (abs(pdg))
      case  (1)
        PDGtoFA = trim(PDGtoFA) // "F[4,{1}]"
      case  (2)
        PDGtoFA = trim(PDGtoFA) // "F[3,{1}]"
      case  (3)
        PDGtoFA = trim(PDGtoFA) // "F[4,{2}]"
      case  (4)
        PDGtoFA = trim(PDGtoFA) // "F[3,{2}]"
      case  (5)
        PDGtoFA = trim(PDGtoFA) // "F[4,{3}]"
      case  (6)
        PDGtoFA = trim(PDGtoFA) // "F[3,{3}]"
      case  (11)
        PDGtoFA = trim(PDGtoFA) // "F[2,{1}]"
      case  (12)
        PDGtoFA = trim(PDGtoFA) // "F[1,{1}]"
      case  (13)
        PDGtoFA = trim(PDGtoFA) // "F[2,{2}]"
      case  (14)
        PDGtoFA = trim(PDGtoFA) // "F[1,{2}]"
      case  (15)
        PDGtoFA = trim(PDGtoFA) // "F[2,{3}]"
      case  (16)
        PDGtoFA = trim(PDGtoFA) // "F[1,{3}]"
      case  (21,9)
        PDGtoFA = "V[5]"
      case  (22)
        PDGtoFA = "V[1]"
      case  (23)
        PDGtoFA = "V[2]"
      case  (24)
        PDGtoFA = trim(PDGtoFA) // "V[3]"
      case  (25)
        PDGtoFA = "S[1]"
      case default
        call ol_msg("Error: only SM particles are allowed!")
        PDGtoFA = "?"
    end select

  end function PDGtoFA


  function ID_to_extparticle(id_in)
    use KIND_TYPES, only: DREALKIND
    use ol_generic, only: to_int, to_lowercase
    use ol_parameters_decl_/**/DREALKIND, only: use_me_cache
  ! MadGraph naming scheme -> PDG
    implicit none
    character(len=*), intent(in) :: id_in
    character(len=len(id_in)) :: id
    type(extparticle) :: ID_to_extparticle

    if (index(id_in, "(") > 0 .and. index(id_in, ")") > 0) then
      ID_to_extparticle%pol = to_int(id_in(index(id_in,"(")+1:index(id_in,")")-1))
      if (ID_to_extparticle%pol /= 0 .and. abs(ID_to_extparticle%pol) /= 1 .and. ID_to_extparticle%pol /= 2) then
        call ol_error("polarization of external particles has to be: 0 (unpol.), -1 (left), 1 (right), 2 (long.)")
      end if
      id = id_in(1:index(id_in,"(")-1)
      ! deactive me_cache for polarized amplitudes
      if (use_me_cache == 1) then
        use_me_cache = 0
        call ol_msg(2,"Matrix element cache deactivated (not available for polarized amplitudes).")
      end if
    else
      ID_to_extparticle%pol = 0
      id = id_in
    end if

    select case (trim(to_lowercase(id)))
      case  ('d')
        ID_to_extparticle%id = 1
      case  ('d~')
        ID_to_extparticle%id = -1
      case  ('u')
        ID_to_extparticle%id = 2
      case  ('u~')
        ID_to_extparticle%id = -2
      case  ('s')
        ID_to_extparticle%id = 3
      case  ('s~')
        ID_to_extparticle%id = -3
      case  ('c')
        ID_to_extparticle%id = 4
      case  ('c~')
        ID_to_extparticle%id = -4
      case  ('b')
        ID_to_extparticle%id = 5
      case  ('b~')
        ID_to_extparticle%id = -5
      case  ('t')
        ID_to_extparticle%id = 6
      case  ('t~')
        ID_to_extparticle%id = -6
      case  ('e-')
        ID_to_extparticle%id = 11
      case  ('e+')
        ID_to_extparticle%id = -11
      case  ('ve', 'ne', 'nu_e')
        ID_to_extparticle%id = 12
      case  ('ve~', 'ne~', 'nu_e~')
        ID_to_extparticle%id = -12
      case  ('mu-', 'm-')
        ID_to_extparticle%id = 13
      case  ('mu+', 'm+')
        ID_to_extparticle%id = -13
      case  ('vm', 'vmu', 'nm', 'nmu', 'nu_mu')
        ID_to_extparticle%id = 14
      case  ('vm~', 'vmu~', 'nm~', 'nmu~', 'nu_mu~')
        ID_to_extparticle%id = -14
      case  ('ta-', 'tau-', 'l-')
        ID_to_extparticle%id = 15
      case  ('ta+', 'tau+', 'l+')
        ID_to_extparticle%id = -15
      case  ('vt', 'vta', 'vtau', 'nt', 'nta', 'ntau', 'nu_tau')
        ID_to_extparticle%id = 16
      case  ('vt~', 'vta~', 'vtau~', 'nt~', 'nta~', 'ntau~', 'nu_tau~')
        ID_to_extparticle%id = -16
      case  ('g')
        ID_to_extparticle%id = 21
      case  ('a') ! on-shell photon
        ID_to_extparticle%id = 22
      case  ('aon') ! off-shell photon
        ID_to_extparticle%id = 2002
      case  ('aoff') ! off-shell photon
        ID_to_extparticle%id = -2002
      case  ('z')
        ID_to_extparticle%id = 23
      case  ('w+')
        ID_to_extparticle%id = 24
      case  ('w-')
        ID_to_extparticle%id = -24
      case  ('h')
        ID_to_extparticle%id = 25
      case  ('h1')
        ID_to_extparticle%id = 25
      case  ('h2', 'h0')
        ID_to_extparticle%id = 35
      case  ('h3', 'a0')
        ID_to_extparticle%id = 36
      case  ('h-', 'hpx')
        ID_to_extparticle%id = -37
      case  ('h+', 'hp')
        ID_to_extparticle%id = 37
      case default
        ID_to_extparticle%id = to_int(trim(id))
        if (ID_to_extparticle%id == -huge(ID_to_extparticle%id)) then
          call ol_error('unrecognised particle id: ' // trim(id))
        end if
    end select
  end function ID_to_extparticle


  function ewcharge(id)
    use KIND_TYPES, only: DREALKIND
    use ol_generic, only: to_string
  ! PDG -> charge
    implicit none
    integer, intent(in) :: id
    real(DREALKIND) :: ewcharge

    select case (id)
      case  (1,3,5)
        ewcharge = -1./3
      case  (-1,-3,-5)
        ewcharge = +1./3
      case  (2,4,6)
        ewcharge = +2./3
      case  (-2,-4,-6)
        ewcharge = -2./3
      case  (11,13,15,-24,-37)
        ewcharge = -1
      case  (-11,-13,-15,24,37)
        ewcharge = +1
      case  (12,14,16,-12,-14,-16,21,22,23,2002,-2002,25,35,36)
       ewcharge = 0
      case default
          call ol_error('unrecognised particle id: ' // to_string(id))
          ewcharge = 0
    end select
  end function ewcharge


  function process_to_extparticlelist(c_in,is_initial)
  ! convert a comma/space/slash separated string of MadGraph or PDG ids into an array of PDG integers
    implicit none
    character(len=*), intent(in) :: c_in
    logical, optional, intent(in) :: is_initial
    character(len(c_in)+1) :: c
    type(extparticle), allocatable :: process_to_extparticlelist(:)
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

    allocate(process_to_extparticlelist(n))

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
        process_to_extparticlelist(n) = ID_to_extparticle(c(pos1+1:i-1))
        if (present(is_initial)) process_to_extparticlelist(n)%is_initial = is_initial
        pos1 = i
        last_seperator = .true.
      else
        last_seperator = .false.
      end if
    end do
  end function process_to_extparticlelist


  subroutine ol_printparameter(filename)
    ! Write parameters to a file.
    ! [in] filename
    use ol_parameters_init_/**/DREALKIND, only: parameters_write
    implicit none
    character(len=*), intent(in) :: filename
    call parameters_write(filename)
  end subroutine ol_printparameter


  subroutine ol_printparameter_c(filename) bind(c,name="ol_printparameter")
    ! C wrapper to ol_printparameter
    ! [in] filename as C string
    use ol_iso_c_utilities, only: c_f_string
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: filename
    character(len=max_parameter_length) :: f_filename
    call c_f_string(filename, f_filename, max_parameter_length)
    call ol_printparameter(trim(f_filename))
  end subroutine ol_printparameter_c


  function register_process_c(process, amptype) bind(c,name="ol_register_process")
    use ol_iso_c_utilities, only: c_f_string
    use ol_parameters_decl_/**/DREALKIND,  only: max_parameter_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: process
    integer(c_int), value :: amptype
    integer(c_int) :: register_process_c
    character(len=max_parameter_length) :: f_process
    integer :: f_amptype
    f_amptype = amptype
    call c_f_string(process, f_process, max_parameter_length)
    register_process_c = register_process(f_process, f_amptype)
  end function register_process_c


  subroutine stop_invalid_id(id)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    if (id <= 0 .or. id > last_process_id) then
      call ol_fatal("Error: no registered process with id " // to_string(id))
      return
    end if
  end subroutine stop_invalid_id

  pure function amplitudetype(id)
    ! [in] id: a process id
    ! return amptype of type integer
    implicit none
    integer, intent(in) :: id
    integer amplitudetype
    ! call stop_invalid_id(id) not possible here,
    ! because 'print' and 'stop' are not allowed in pure functions.
    if (id <= 0 .or. id > last_process_id) then
      amplitudetype = 0
    else
      amplitudetype = process_handles(id)%amplitude_type
    end if
  end function amplitudetype


  function amplitudetype_c(id) bind(c,name="ol_amplitudetype")
    ! [in] id: a process id
    ! return amptype of type integer
    implicit none
    integer(c_int), value :: id
    integer(c_int) :: amplitudetype_c
    call stop_invalid_id(int(id))
    if (error > 1) return
    amplitudetype_c = process_handles(int(id))%amplitude_type
  end function amplitudetype_c


  function library_content_c(id) bind(c,name='ol_library_content')
    implicit none
    integer(c_int), value :: id
    integer(c_int) :: library_content_c
    call stop_invalid_id(int(id))
    if (error > 1) return
    library_content_c = process_handles(int(id))%content
  end function library_content_c


  pure function n_external(id)
    implicit none
    integer, intent(in) :: id
    integer :: n_external
    ! call stop_invalid_id(id) not possible here,
    ! because 'print' and 'stop' are not allowed in pure functions.
    if (id <= 0 .or. id > last_process_id) then
      n_external = 0
    else
      n_external = process_handles(id)%n_particles
    end if
  end function n_external


  function n_external_c(id) bind(c,name="ol_n_external")
    implicit none
    integer(c_int), value :: id
    integer(c_int) :: n_external_c
    call stop_invalid_id(int(id))
    if (error > 1) return
    n_external_c = process_handles(int(id))%n_particles
  end function n_external_c


  subroutine phase_space_point(id, sqrt_s, psp)
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: sqrt_s
    real(DREALKIND), intent(out) :: psp(:,:)
    type(process_handle) :: subprocess
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    call subprocess%set_permutation(subprocess%permutation)
    n_scatt = subprocess%n_in
    call subprocess%rambo(sqrt_s, psp)
  end subroutine phase_space_point


  subroutine phase_space_point_c(id, sqrt_s, pp) bind(c,name="ol_phase_space_point")
    implicit none
    integer(c_int), value :: id
    real(c_double), value :: sqrt_s
    real(c_double), intent(out) :: pp(5*n_external(int(id)))
    type(process_handle) :: subprocess
    integer :: i
    real(DREALKIND) :: f_sqrt_s
    real(DREALKIND) :: f_psp(0:3,n_external(int(id)))
    ! call stop_invalid_id(id) not needed here
    i = id
    f_sqrt_s = sqrt_s
    subprocess = process_handles(i)
    call phase_space_point(i, f_sqrt_s, f_psp)
    do i = 1, subprocess%n_particles
      pp(5*(i-1)+1:5*(i-1)+4) = f_psp(0:3,i)
      pp(5*i) = subprocess%masses(i)
    end do
  end subroutine phase_space_point_c


  subroutine tree_colbasis_dim(id, ncolb, colelemsz, nhel)
    ! for process with id 'id' return
    ! ncolb = number of tree colour basis elements;
    ! colelemsz = number of colour colour indices in a colour basis element
    ! nhel = number of helicity configuration, including those which vanish
    implicit none
    integer, intent(in) :: id
    integer, intent(out) :: ncolb, colelemsz, nhel
    integer :: extcols(n_external(id)), ncoupl, maxpows, ncolext
    call stop_invalid_id(id)
    if (error > 1) return
    if (.not. associated(process_handles(id)%tree_colbasis_dim)) then
      call ol_msg("Error: colour basis information is not available")
      call ol_fatal("       for process " // process_handles(id)%process_name)
      return
    end if
    call process_handles(id)%tree_colbasis_dim(extcols, ncolb, ncoupl, maxpows, nhel)
    ncolext = count(extcols /= 0)
    colelemsz = ncolext/2 + ncolext - 1
  end subroutine tree_colbasis_dim


  subroutine tree_colbasis_dim_c(id, ncolb, colelemsz, nhel) bind(c,name="ol_tree_colbasis_dim")
    implicit none
    integer(c_int), value :: id
    integer(c_int), intent(out) :: ncolb, colelemsz, nhel
    integer :: f_ncolb, f_colelemsz, f_nhel
    ! call stop_invalid_id(id) not needed here
    call tree_colbasis_dim(int(id), f_ncolb, f_colelemsz, f_nhel)
    ncolb = f_ncolb
    colelemsz = f_colelemsz
    nhel = f_nhel
  end subroutine tree_colbasis_dim_c


  pure function get_tree_colbasis_dim(id)
    ! number of tree colour basis elements; used to declare array sizes
    implicit none
    integer, intent(in) :: id
    integer :: get_tree_colbasis_dim
    integer :: extcols(n_external(id)), ncoupl, maxpows, nhel
    ! call stop_invalid_id(id) not possible here,
    ! because 'print' and 'stop' are not allowed in pure functions.
    if (id <= 0 .or. id > last_process_id) then
      get_tree_colbasis_dim = 0
    else
      call process_handles(id)%tree_colbasis_dim(extcols, get_tree_colbasis_dim, ncoupl, maxpows, nhel)
    end if
  end function get_tree_colbasis_dim


  pure function tree_colbasis_elemsize(id)
    ! number of coloured external particles; used to declare array sizes
    implicit none
    integer, intent(in) :: id
    integer :: tree_colbasis_elemsize
    integer :: extcols(n_external(id)), ncolb, ncoupl, maxpows, nhel, ncolext
    ! call stop_invalid_id(id) not possible here,
    ! because 'print' and 'stop' are not allowed in pure functions.
    if (id <= 0 .or. id > last_process_id) then
      tree_colbasis_elemsize = 0
    else
      call process_handles(id)%tree_colbasis_dim(extcols, ncolb, ncoupl, maxpows, nhel)
      ncolext = count(extcols /= 0)
      tree_colbasis_elemsize = ncolext/2+ncolext-1
      if (ncolext == 0) tree_colbasis_elemsize = 0
    end if
  end function tree_colbasis_elemsize


  pure function get_nhel(id)
    ! number of helicity configurations (all, not just non-vanishing); used to declare array sizes
    implicit none
    integer, intent(in) :: id
    integer :: get_nhel
    integer :: extcols(n_external(id)), ncolb, ncoupl, maxpows
    ! call stop_invalid_id(id) not possible here,
    ! because 'print' and 'stop' are not allowed in pure functions.
    if (id <= 0 .or. id > last_process_id) then
      get_nhel = 0
    else
      call process_handles(id)%tree_colbasis_dim(extcols, ncolb, ncoupl, maxpows, get_nhel)
    end if
  end function get_nhel


  subroutine tree_colbasis(id, basis, needed)
    use ol_generic, only: compositions2, nth_permutation
    implicit none
    integer, intent(in) :: id
    integer, intent(out) :: basis(:,:)
    integer, intent(out) :: needed(:,:)
    integer :: extcols(n_external(id)), ncolb, ncoupl, maxpows, nhel
    integer, allocatable :: pbasis(:,:), selected_powers(:,:), perm(:), compos(:,:), compo(:), basiselem(:)
    integer :: ncolext, i, j, k, m, needij
    integer :: ncol2ext(n_external(id)), invextperm(n_external(id))
    logical :: powok
    call stop_invalid_id(id)
    if (error > 1) return
    call process_handles(id)%tree_colbasis_dim(extcols, ncolb, ncoupl, maxpows, nhel)
    ! number of coloured external particles
    ncolext = 0
    do i = 1, size(extcols)
      if (extcols(i) /= 0) then
        ncolext = ncolext + 1
        ncol2ext(ncolext) = i
      end if
    end do
    do i = 1, size(invextperm)
      invextperm(process_handles(id)%permutation(i)) = i
    end do
    allocate(pbasis(ncoupl+2,ncolb))
    allocate(selected_powers(maxpows,ncoupl))
    allocate(perm(ncolext))
    call compositions2(compos, ncolext)
    allocate(compo(size(compos,1)))
    allocate(basiselem(ncolext/2+ncolext-1))
    call process_handles(id)%tree_colbasis(pbasis, selected_powers)
    do i = 1, ncolb
      do j = i, ncolb
        needij = 1
        do k = 1, ncoupl
          powok = .false.
          do m = 1, maxpows
            if (pbasis(2+k,i) + pbasis(2+k,j) == selected_powers(m,k)) then
              powok = .true.
              exit
            end if
          end do
          if (.not. powok) then
            needij = 0
            exit
          end if
        end do
        needed(i,j) = needij
        needed(j,i) = needij
      end do
      ! TODO:
      ! - apply crossing
      compo = compos(:,pbasis(1,i))
      perm = nth_permutation([(k, k=1, ncolext)], pbasis(2,i))
      basiselem = 0
      m = 1
      do j = 1, count(compo > 0)
        do k = 1, compo(j)
          basiselem(m+j-1) = invextperm(ncol2ext(perm(m)))
          m = m + 1
        end do
      end do
      basis(:,i) = basiselem
    end do
    deallocate(pbasis)
    deallocate(selected_powers)
    deallocate(perm)
    deallocate(compos)
    deallocate(compo)
    deallocate(basiselem)
  end subroutine tree_colbasis


  subroutine tree_colbasis_c(id, basis, needed) bind(c,name="ol_tree_colbasis")
    implicit none
    integer(c_int), value :: id
    integer(c_int), intent(out) :: basis(tree_colbasis_elemsize(id),get_tree_colbasis_dim(id))
    integer(c_int), intent(out) :: needed(get_tree_colbasis_dim(id),get_tree_colbasis_dim(id))
    integer :: f_basis(tree_colbasis_elemsize(id),get_tree_colbasis_dim(id))
    integer :: f_needed(get_tree_colbasis_dim(id),get_tree_colbasis_dim(id))
    ! call stop_invalid_id(id) not needed here
    call tree_colbasis(int(id), f_basis, f_needed)
    basis = f_basis
    needed = f_needed
  end subroutine tree_colbasis_c


  subroutine tree_colourflow(id, flowbasis)
    ! Convert the trace colour basis for tree amplitudes from tree_colbasis()
    ! to a colour flow basis and return it as a rank 3 integer array
    ! flowbasis(2,nexternal,basissize). For each basis element and each
    ! external particle it contains a pair (i,j) with the incoming fundamental
    ! colour i and outgoing colour j. A zero entry means uncoloured.
    ! Incoming quark / outgoing antiquark: (i,0)
    ! Incoming antiquark / outgoing quark: (0,i)
    ! In T(a1,a2,...,an)_ij, j connects to the gluon 'an' and to an incoming quark.
    implicit none
    integer, intent(in) :: id
    integer, intent(out) :: flowbasis(2,n_external(id),get_tree_colbasis_dim(id))
    integer :: tracebasis(tree_colbasis_elemsize(id),get_tree_colbasis_dim(id))
    integer :: neededdummy(get_tree_colbasis_dim(id),get_tree_colbasis_dim(id))
    integer :: ncolb, colelemsz, i, k, lastcol, startmon, endmon
    integer :: powersof2(n_external(id))

    ncolb = size(tracebasis,2)
    colelemsz = size(tracebasis,1)
    flowbasis = 0
    call tree_colbasis(id, tracebasis, neededdummy)

    do k = 1, ncolb ! for each basis element
      do lastcol = colelemsz, 1, -1
        ! size of trace basis element with trailing zeroes dropped
        if (tracebasis(lastcol,k) /= 0) exit
      end do
      if (lastcol == 0) then
        ! no coloured external particles
        return
      end if
      startmon = 1
      do ! for each colour monomial
        endmon = locateval(tracebasis(:lastcol,k), 0, startmon) - 1
        if (endmon <= 0) endmon = lastcol
        if (any(process_handles(id)%extid(tracebasis(endmon,k)) == pdgadjoint)) then
          ! if the last index is adjoint (gluon),
          ! the monomial is a trace
          call tracetoflow(tracebasis(startmon:endmon,k), flowbasis(:,:,k))
        else
          ! the monomial is a chain
          call chaintoflow(tracebasis(startmon:endmon,k), flowbasis(:,:,k))
        end if
        startmon = endmon + 2
        if (startmon > lastcol) exit
      end do
    end do

  contains

    function locateval(arr, val, startpos)
      ! Return the first position (starting from startpos if present) in the
      ! array arr which holds the value val, or -1 if val is not in arr(startpos:).
      implicit none
      integer :: locateval
      integer, intent(in) :: arr(:), val
      integer, intent(in), optional :: startpos
      integer :: sp, k
      sp = 1
      if (present(startpos)) sp = startpos
      if (sp < 1) sp = 1
      locateval = -1
      do k = sp, size(arr)
        if (arr(k) == val) then
          locateval = k
          return
        end if
      end do
    end function locateval

    subroutine chaintoflow(colmon, flowbasisel)
      ! Add information from the colour chain colmon
      ! to the colour flow basis element flowbasisel.
      implicit none
      integer, intent(in) :: colmon(:)
      integer, intent(inout) :: flowbasisel(:,:)
      integer :: k, ng, colin, colout
      ! chain T(a1,a2,...,an)_ij
      ! number of gluons
      ng = size(colmon) - 2
      colin = colmon(ng+1)
      flowbasisel(:,colmon(ng+1)) = [0,colin] ! anti-quark (0,i)
      do k = 1, ng
        ! (i,a1)(a1,a2)...(an-1,an)
        colout = colmon(k)
        flowbasisel(:,colout) = [colin,colout]
        colin = colout
      end do
      flowbasisel(:,colmon(ng+2)) = [colin,0] ! quark (an,0), j->an
    end subroutine chaintoflow

    subroutine tracetoflow(colmon, flowbasisel)
      ! Add information from the colour trace colmon
      ! to the colour flow basis element flowbasisel.
      implicit none
      integer, intent(in) :: colmon(:)
      integer, intent(inout) :: flowbasisel(:,:)
      integer :: k, ng, colin, colout
      ! trace T(a1,a2,...,an)
      ! number of gluons
      ng = size(colmon)
      colin = colmon(ng)
      do k = 1, ng
        ! (a1,an)(a2,a1)...(an,an-1)
        colout = colmon(k)
        flowbasisel(:,colout) = [colin,colout]
        colin = colout
      end do
    end subroutine tracetoflow

  end subroutine tree_colourflow


  subroutine tree_colourflow_c(id, flowbasis) bind(c,name="ol_tree_colourflow")
    ! C interface to tree_colourflow(). Return flowbasis as
    ! int flowbasis[ncolb][nex][2];
    implicit none
    integer(c_int), value :: id
    integer(c_int), intent(out) :: flowbasis(2,n_external(int(id)),get_tree_colbasis_dim(int(id)))
    integer :: f_flowbasis(2,n_external(int(id)),get_tree_colbasis_dim(int(id)))
    call tree_colourflow(int(id), f_flowbasis)
    flowbasis = f_flowbasis
  end subroutine tree_colourflow_c


  subroutine start() bind(c,name="ol_start")
    use ol_parameters_decl_/**/DREALKIND,  only: &
      & write_params_at_start, stability_logdir_not_created, stability_log, stability_logdir, &
      & write_psp, ti_monitor
    use ol_parameters_init_/**/DREALKIND, only: parameters_write
    use ol_cwrappers, only: mkdir
    implicit none
    integer :: mkdirerr
    call parameters_flush()
    if (stability_logdir_not_created .and. (stability_log > 0 .or. &
      & write_psp > 0 .or. ti_monitor > 0)) then
      stability_logdir_not_created = .false.
      mkdirerr = mkdir(stability_logdir)
    end if
    if (write_params_at_start) call parameters_write()
  end subroutine


  subroutine finish() bind(c,name="ol_finish")
    implicit none
    call cleanup()
    call unregister_processes()
    if (shopping_list_open) close(fh_shopping)
  end subroutine finish


  subroutine evaluate_tree(id, psp, res)
   ! Tree matrix element.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: squared tree matrix element
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res
    real(DREALKIND) :: m2cc(0:n_external(id)*(n_external(id)+1)/2+1)
    real(DREALKIND) :: resmunu(4,4)
    type(process_handle) :: subprocess
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal("evaluate: tree routine not available for process " // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call tree_parameters_flush()
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call subprocess%tree(psp, m2cc, 0, &
      & [0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND], &
      & 1, [0], resmunu)
    res = m2cc(0)
  end subroutine evaluate_tree


  subroutine evaluate_tree_c(id, pp, res) bind(c,name="ol_evaluate_tree")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_tree(f_id, f_pp(0:3,:), f_res)
    res = f_res
  end subroutine evaluate_tree_c


  subroutine evaluate_tree_colvect(id, psp, amp, nhel)
    ! Tree amplitude as colour vectors for each helicity configuration.
    ! [in] id: process id as set by register_process
    ! [in] psp: phase space point
    ! [out] amp: amp(:,h) is the colour vector for helicity configuration h
    ! [out] nhel: number of non-zero helicity configurations,
    !       amp(:,nhel+1:) contains no information
    use ol_ew_renormalisation_/**/REALKIND, only: photon_factors
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    complex(DREALKIND), intent(out) :: amp(:,:)
    integer, intent(out) :: nhel
    real(DREALKIND) :: res, bornphotonfactor
    call evaluate_tree(id, psp, res) ! fill colour vector cache
    call process_handles(id)%tree_colvect(amp, nhel)
    call photon_factors(process_handles(id)%photon_id, 0, bornphotonfactor)
    amp = amp * sqrt(bornphotonfactor)
  end subroutine evaluate_tree_colvect


  subroutine evaluate_tree_colvect_c(id, pp, amp, nhel) bind(c,name="ol_evaluate_tree_colvect")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: amp(2*get_tree_colbasis_dim(id),get_nhel(id))
    integer(c_int), intent(out) :: nhel
    real(c_double) :: res
    complex(DREALKIND) :: f_amp(get_tree_colbasis_dim(id),get_nhel(id))
    integer :: f_nhel, k, h
    call evaluate_tree_c(id, pp, res) ! fill colour vector cache
    call process_handles(int(id))%tree_colvect(f_amp, f_nhel)
    do h = 1, f_nhel
      do k = 1, size(f_amp,1)
        amp(2*k-1,h) = real(f_amp(k,h))
        amp(2*k,h) = aimag(f_amp(k,h))
      end do
    end do
    nhel = f_nhel
  end subroutine evaluate_tree_colvect_c


  subroutine evaluate_tree_colvect2(id, psp, m2arr)
    ! Helicity summed squared tree colour vector.
    ! Apart from the missing colour factor these are the squared matrix elements
    ! for each colour flow.
    ! [in] id: process id as set by register_process
    ! [in] psp: phase space point
    ! [out] m2arr: array of squared matrix elements per colour flow
    use ol_ew_renormalisation_/**/REALKIND, only: photon_factors
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2arr(:)
    integer :: nhel
    real(DREALKIND) :: res, bornphotonfactor
    complex(DREALKIND) :: amp(get_tree_colbasis_dim(id),get_nhel(id))
    call evaluate_tree(id, psp, res) ! fill colour vector cache
    call process_handles(id)%tree_colvect(amp, nhel)
    m2arr = sum(real(amp(:,1:nhel)*conjg(amp(:,1:nhel))), 2)
    call photon_factors(process_handles(id)%photon_id, 0, bornphotonfactor)
    m2arr = m2arr * bornphotonfactor
  end subroutine evaluate_tree_colvect2


  subroutine evaluate_tree_colvect2_c(id, pp, m2arr) bind(c,name="ol_evaluate_tree_colvect2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2arr(get_tree_colbasis_dim(id))
    integer :: nhel
    real(c_double) :: res
    complex(DREALKIND) :: amp(get_tree_colbasis_dim(id),get_nhel(id))
    call evaluate_tree_c(id, pp, res) ! fill colour vector cache
    call process_handles(int(id))%tree_colvect(amp, nhel)
    m2arr = sum(real(amp(:,1:nhel)*conjg(amp(:,1:nhel))), 2)
  end subroutine evaluate_tree_colvect2_c


  subroutine evaluate_cc(id, psp, tree, cc, ewcc)
   ! Independent color correlated tree matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] tree: Born matrix element
   ! [out] cc(n_external*(n_external-1)/2): array with the indepenent color correlated
   !       tree amplitudes C_ij = <M|T_iT_j|M>
   !       cc(i+j(j-1)/2+1) = C_ij with 0 <= i < j <= n_external-1
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: tree, cc(:), ewcc
    type(process_handle) :: subprocess
    real(DREALKIND) :: m2cc(0:n_external(id)*(n_external(id)+1)/2+1) ! keep +1 for compatibility
    real(DREALKIND) :: resmunu(4,4)
    integer  :: n_cc, i, j
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal('evaluate: cc routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2
    call tree_parameters_flush()
    call subprocess%tree(psp, m2cc, 0, &
      & [0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND], &
      & n_cc, [(i, i = 0, n_cc)], resmunu)
    tree = m2cc(0)
    ewcc = m2cc(n_cc+1)
    do j = 1, subprocess%n_particles - 1
      do i = 0, j - 1
        cc(i+j*(j-1)/2+1) = m2cc((j+1)*j/2+i+1)
      end do
    end do
  end subroutine evaluate_cc


  subroutine evaluate_cc_c(id, pp, tree, cc, ewcc) bind(c,name="ol_evaluate_cc")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: tree, cc(rval_size(n_external(id),2)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_tree, f_cc(rval_size(n_external(id),2)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_cc(f_id, f_pp(0:3,:), f_tree, f_cc, f_ewcc)
    tree = f_tree
    cc = f_cc
    ewcc = f_ewcc
  end subroutine evaluate_cc_c


  subroutine evaluate_cc2(id, psp, res, cc, ewcc)
   ! Independent color correlated loop^2 matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: loop^2 matrix element
   ! [out] cc(n_external*(n_external-1)/2): array with the indepenent color correlated
   !       loop^2 amplitudes C_ij = <M1|T_iT_j|M1>
   !       cc(i+j(j-1)/2+1) = C_ij with 0 <= i < j <= n_external-1
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, cc(:), ewcc
    type(process_handle) :: subprocess
    integer  :: n_cc, i, j
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(correlator) :: Lcc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 2)) then
      call ol_fatal('evaluate: cc2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2
    call parameters_flush()
    Lcc%type=1
    allocate(Lcc%rescc(0:n_external(id)*(n_external(id)+1)/2+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lcc, acc)
    res = Lcc%rescc(0)
    ewcc = Lcc%rescc(n_cc+1)
    do j = 1, subprocess%n_particles - 1
      do i = 0, j - 1
        cc(i+j*(j-1)/2+1) = Lcc%rescc((j+1)*j/2+i+1)
      end do
    end do
    deallocate(Lcc%rescc)
  end subroutine evaluate_cc2


  subroutine evaluate_cc2_c(id, pp, res, cc, ewcc) bind(c,name="ol_evaluate_cc2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, cc(rval_size(n_external(id),2)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_cc(rval_size(n_external(id),2)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_cc2(f_id, f_pp(0:3,:), f_res, f_cc, f_ewcc)
    res = f_res
    cc = f_cc
    ewcc = f_ewcc
  end subroutine evaluate_cc2_c


  subroutine evaluate_loopcc(id, psp, m2l0, m2l1, cc, ewcc)
   ! Independent color correlated born X loop matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] m2l0: born matrix element
   ! [out] m2l1: born X loop matrix element
   ! [out] cc(n_external*(n_external-1)/2): array with the indepenent color correlated
   !       born X loop C_ij = <M0|T_iT_j|M1>
   !       cc(i+j(j-1)/2+1) = C_ij with 0 <= i < j <= n_external-1
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), cc(:), ewcc
    type(process_handle) :: subprocess
    integer  :: n_cc, i, j
    real(DREALKIND) :: ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(correlator) :: Lcc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loopcc routine not available for process ' // trim(to_string(id)))
      return
    end if
    if (.not. associated(subprocess%loopcr)) then
      call ol_fatal('evaluate: loopcc routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2
    call parameters_flush()
    Lcc%type=11
    allocate(Lcc%rescc(0:n_external(id)*(n_external(id)+1)/2+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lcc, acc)
    ewcc = Lcc%rescc(n_cc+1)
    do j = 1, subprocess%n_particles - 1
      do i = 0, j - 1
        cc(i+j*(j-1)/2+1) = Lcc%rescc((j+1)*j/2+i+1)
      end do
    end do
    deallocate(Lcc%rescc)
  end subroutine evaluate_loopcc


  subroutine evaluate_loopcc_c(id, pp, m2l0, m2l1, cc, ewcc) bind(c,name="ol_evaluate_loopcc")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), cc(rval_size(n_external(id),2)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_cc(rval_size(n_external(id),2)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loopcc(f_id, f_pp(0:3,:), f_m2l0, f_m2l1, f_cc, f_ewcc)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    cc = f_cc
    ewcc = f_ewcc
  end subroutine evaluate_loopcc_c


  subroutine evaluate_ccmatrix(id, psp, tree, ccij, ewcc)
   ! Color correlated tree matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] tree: squared born matrix element
   ! [out] cc(n_external:n_external): array with the color correlated
   !       tree amplitudes C_ij = <M|T_iT_j|M>
   !       cc(i,j) = C_ij with i,j = 1 <= n_external
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: tree, ccij(:,:), ewcc
    type(process_handle) :: subprocess
    real(DREALKIND) :: m2cc(0:n_external(id)*(n_external(id)+1)/2+1)
    integer  :: n_cc, i, j
    real(DREALKIND) :: resmunu(4,4)
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal('evaluate: ccmatrix routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2+1
    call tree_parameters_flush()
    call subprocess%tree(psp, m2cc, 0, &
      & [0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND], &
      & n_cc, [(i, i = 0, n_cc)], resmunu)
    do i = 1, subprocess%n_particles
      do j = 1, i
        ccij(i,j) = m2cc(i*(i-1)/2+j)
        if (i /= j) ccij(j,i) = ccij(i,j)
      end do
    end do
    tree = m2cc(0)
    ewcc = m2cc(n_cc)
  end subroutine evaluate_ccmatrix


  subroutine evaluate_ccmatrix_c(id, pp, tree, ccij, ewcc) bind(c,name="ol_evaluate_ccmatrix")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: tree, ccij(n_external(id)*n_external(id)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_tree, f_ccij(n_external(id),n_external(id)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_ccmatrix(f_id, f_pp(0:3,:), f_tree, f_ccij, f_ewcc)
    tree = f_tree
    ccij = reshape(f_ccij,(/n_external(id)*n_external(id)/))
    ewcc = f_ewcc
  end subroutine evaluate_ccmatrix_c


  subroutine evaluate_ccmatrix2(id, psp, res, ccij, ewcc)
   ! Color correlated loop^2 matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: loop^2 matrix element
   ! [out] cc(n_external:n_external): array with the color correlated
   !       loop^2 amplitudes C_ij = <M|T_iT_j|M>
   !       cc(i,j) = C_ij with i,j = 1 <= n_external
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, ccij(:,:), ewcc
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(process_handle) :: subprocess
    integer  :: n_cc, i, j
    type(correlator) :: Lcc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 2)) then
      call ol_fatal('evaluate: ccmatrix2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2+1
    call parameters_flush()
    Lcc%type=1
    allocate(Lcc%rescc(0:n_external(id)*(n_external(id)+1)/2+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lcc, acc)
    do i = 1, subprocess%n_particles
      do j = 1, i
        ccij(i,j) = Lcc%rescc(i*(i-1)/2+j)
        if (i /= j) ccij(j,i) = ccij(i,j)
      end do
    end do
    res = Lcc%rescc(0)
    ewcc = Lcc%rescc(n_cc)
    deallocate(Lcc%rescc)
  end subroutine evaluate_ccmatrix2


  subroutine evaluate_ccmatrix2_c(id, pp, res, ccij, ewcc) bind(c,name="ol_evaluate_ccmatrix2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, ccij(n_external(id)*n_external(id)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_ccij(n_external(id),n_external(id)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_ccmatrix2(f_id, f_pp(0:3,:), f_res, f_ccij, f_ewcc)
    res = f_res
    ccij = reshape(f_ccij,(/n_external(id)*n_external(id)/))
    ewcc = f_ewcc
  end subroutine evaluate_ccmatrix2_c


  subroutine evaluate_loopccmatrix(id, psp, m2l0, m2l1, ccij, ewcc)
   ! Color correlated tree matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] tree: squared born matrix element
   ! [out] cc(n_external:n_external): array with the color correlated
   !       tree amplitudes C_ij = <M|T_iT_j|M>
   !       cc(i,j) = C_ij with i,j = 1 <= n_external
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), ccij(:,:), ewcc
    type(process_handle) :: subprocess
    real(DREALKIND) :: m2cc(0:n_external(id)*(n_external(id)+1)/2+1)
    real(DREALKIND) :: ir1(0:2), m2l2(0:4), ir2(0:4), acc
    integer  :: n_cc, i, j
    type(correlator) :: Lcc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loopccmatrix routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    n_cc = subprocess%n_particles*(subprocess%n_particles+1)/2
    call parameters_flush()
    Lcc%type=11
    allocate(Lcc%rescc(0:n_cc+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lcc, acc)
    ewcc = Lcc%rescc(n_cc+1)
    do i = 1, subprocess%n_particles
      do j = 1, i
        ccij(i,j) = Lcc%rescc(i*(i-1)/2+j)
        if (i /= j) ccij(j,i) = ccij(i,j)
      end do
    end do
    deallocate(Lcc%rescc)
  end subroutine evaluate_loopccmatrix


  subroutine evaluate_loopccmatrix_c(id, pp, m2l0, m2l1, ccij, ewcc) bind(c,name="ol_evaluate_loopccmatrix")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), ccij(n_external(id)*n_external(id)), ewcc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_ccij(n_external(id),n_external(id)), f_ewcc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loopccmatrix(f_id, f_pp(0:3,:), f_m2l0, f_m2l1, f_ccij, f_ewcc)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    ccij = reshape(f_ccij,(/n_external(id)*n_external(id)/))
    ewcc = f_ewcc
  end subroutine evaluate_loopccmatrix_c


  subroutine evaluate_ccewmatrix(id, psp, tree, ccij)
   ! Charge correlated tree matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] tree: squared born matrix element
   ! [out] cc(n_external:n_external): array with the charge correlated
   !       tree amplitudes C_ij = <M|Q_iQ_j|M>
   !       cc(i,j) = C_ij with i,j = 1 <= n_external
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: tree, ccij(:,:)
    type(process_handle) :: subprocess
    integer  :: i, j
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal('evaluate: ccewmatrix routine not available for process ' // trim(to_string(id)))
      return
    end if
    call evaluate_tree(id,psp,tree)
    do i = 1, subprocess%n_particles
      do j = 1, i
        ccij(i,j) = tree*ewcharge(subprocess%extid(i))*ewcharge(subprocess%extid(j))
        if (i /= j) ccij(j,i) = ccij(i,j)
      end do
    end do
  end subroutine evaluate_ccewmatrix


  subroutine evaluate_ccewmatrix_c(id, pp, tree, ccij) bind(c,name="ol_evaluate_ccewmatrix")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: tree, ccij(n_external(id)*n_external(id))
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_tree, f_ccij(n_external(id),n_external(id))
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_ccewmatrix(f_id, f_pp(0:3,:), f_tree, f_ccij)
    tree = f_tree
    ccij = reshape(f_ccij,(/n_external(id)*n_external(id)/))
  end subroutine evaluate_ccewmatrix_c


  subroutine evaluate_ccewmatrix2(id, psp, res, ccij)
   ! Charge correlated loop-squared matrix elements.
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: squared loop matrix element
   ! [out] cc(n_external:n_external): array with the charge correlated
   !       tree amplitudes C_ij = <M|Q_iQ_j|M>
   !       cc(i,j) = C_ij with i,j = 1 <= n_external
   ! [out] ewcc: charge correlation for EW i-operator
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, ccij(:,:)
    real(DREALKIND) :: acc
    type(process_handle) :: subprocess
    integer  :: i, j
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 2)) then
      call ol_fatal('evaluate: ccewmatrix routine not available for process ' // trim(to_string(id)))
      return
    end if
    call evaluate_loop2(id,psp,res,acc)
    do i = 1, subprocess%n_particles
      do j = 1, i
        ccij(i,j) = res*ewcharge(subprocess%extid(i))*ewcharge(subprocess%extid(j))
        if (i /= j) ccij(j,i) = ccij(i,j)
      end do
    end do
  end subroutine evaluate_ccewmatrix2


  subroutine evaluate_ccewmatrix2_c(id, pp, res, ccij) bind(c,name="ol_evaluate_ccewmatrix2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, ccij(n_external(id)*n_external(id))
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_ccij(n_external(id),n_external(id))
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_ccewmatrix2(f_id, f_pp(0:3,:), f_res, f_ccij)
    res = f_res
    ccij = reshape(f_ccij,(/n_external(id)*n_external(id)/))
  end subroutine evaluate_ccewmatrix2_c


  subroutine evaluate_sc(id, psp, emitter, polvect, res)
    ! Spin correlated matrix elements.
    ! [in] id: process id as set by register_process
    ! [in] psp: phase space point
    ! [in] int emitter: emitter
    ! [in] polvect: polarisation vector
    ! [out] res(n_external): array with results for each spectator j,
    !       res(j) = 1/mom^2 * <emitter,mu|mom^mu cc_ij mom^nu|j,nu>
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:), polvect(4)
    real(DREALKIND), intent(out) :: res(:)
    type(process_handle) :: subprocess
    integer :: j, extcombs(n_external(id)), nextcombs
    real(DREALKIND) :: m2sc(0:n_external(id)*(n_external(id)+1)/2+1)
    real(DREALKIND) :: resmunu(4,4)
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal('evaluate: sc routine not available for process ' // trim(to_string(id)))
      return
    end if
    if (subprocess%extid(emitter) == 21) then
      do j = 1, subprocess%n_particles
        if (j <= emitter) then
          extcombs(j) = emitter*(emitter-1)/2 + j
        else
          extcombs(j) = j*(j-1)/2 + emitter
        end if
      end do
      nextcombs = subprocess%n_particles
    else if (subprocess%extid(emitter) == 22) then    ! no color insertions for photon emitter
      extcombs = 0
      nextcombs = 1
    else
      res = 0
      return
    end if
    n_scatt = subprocess%n_in
    call tree_parameters_flush()
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call subprocess%tree(psp, m2sc, emitter, polvect, nextcombs, extcombs, resmunu)
    do j = 1, subprocess%n_particles
      res(j) = m2sc(extcombs(j))
    end do
  end subroutine evaluate_sc


  subroutine evaluate_sc_c(id, pp, emitter, polvect, res) bind(c,name="ol_evaluate_sc")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id)), polvect(4)
    real(c_double), intent(out) :: res(n_external(id))
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id)), f_polvect(4)
    real(DREALKIND) :: f_res(n_external(id))
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    f_polvect = polvect
    call evaluate_sc(f_id, f_pp(0:3,:), f_emitter, f_polvect, f_res)
    res = f_res
  end subroutine evaluate_sc_c


  subroutine evaluate_sc2(id, psp, emitter, polvect, res)
    ! Spin correlated loop^2 matrix elements.
    ! [in] id: process id as set by register_process
    ! [in] psp: phase space point
    ! [in] int emitter: emitter
    ! [in] polvect: polarisation vector
    ! [out] res(n_external): array with results for each spectator j,
    !       res(j) = 1/mom^2 * <emitter,mu|mom^mu cc_ij mom^nu|j,nu>
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:), polvect(4)
    real(DREALKIND), intent(out) :: res(:)
    type(process_handle) :: subprocess
    integer :: j, extcombs(n_external(id)), nextcombs
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(correlator) :: Lsc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 2)) then
      call ol_fatal('evaluate: sc2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    if (subprocess%extid(emitter) == 21) then ! color insertions for gluon emitter
      do j = 1, subprocess%n_particles
        if (j <= emitter) then
          extcombs(j) = emitter*(emitter-1)/2 + j
        else
          extcombs(j) = j*(j-1)/2 + emitter
        end if
      end do
      nextcombs = subprocess%n_particles
    else if (subprocess%extid(emitter) == 22) then    ! no color insertions for photon emitter
      extcombs = 0
      nextcombs = 1
    else
      res = 0
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    Lsc%type=2
    Lsc%emitter=emitter
    Lsc%mom=polvect
    Lsc%nextcombs=nextcombs
    allocate(Lsc%extcombs(nextcombs))
    Lsc%extcombs = extcombs
    allocate(Lsc%rescc(0:n_external(id)*(n_external(id)+1)/2+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lsc, acc)
    do j = 1, subprocess%n_particles
      res(j) = Lsc%rescc(extcombs(j))
    end do
    deallocate(Lsc%extcombs)
    deallocate(Lsc%rescc)
  end subroutine evaluate_sc2


  subroutine evaluate_sc2_c(id, pp, emitter, polvect, res) bind(c,name="ol_evaluate_sc2")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id)), polvect(4)
    real(c_double), intent(out) :: res(n_external(id))
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id)), f_polvect(4)
    real(DREALKIND) :: f_res(n_external(id))
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    f_polvect = polvect
    call evaluate_sc2(f_id, f_pp(0:3,:), f_emitter, f_polvect, f_res)
    res = f_res
   end subroutine evaluate_sc2_c


  subroutine evaluate_loopsc(id, psp, emitter, polvect, res)
    ! Spin correlated boonXloop matrix elements.
    ! [in] id: process id as set by register_process
    ! [in] psp: phase space point
    ! [in] int emitter: emitter
    ! [in] polvect: polarisation vector
    ! [out] res(n_external): array with results for each spectator j,
    !       res(j) = 1/mom^2 * <emitter,mu|mom^mu cc_ij mom^nu|j,nu>
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:), polvect(4)
    real(DREALKIND), intent(out) :: res(:)
    type(process_handle) :: subprocess
    integer :: j, extcombs(n_external(id)), nextcombs
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(correlator) :: Lsc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loopsc routine not available for process ' // trim(to_string(id)))
      return
    end if
    if (subprocess%extid(emitter) == 21) then ! color insertions for gluon emitter
      do j = 1, subprocess%n_particles
        if (j <= emitter) then
          extcombs(j) = emitter*(emitter-1)/2 + j
        else
          extcombs(j) = j*(j-1)/2 + emitter
        end if
      end do
      nextcombs = subprocess%n_particles
    else if (subprocess%extid(emitter) == 22) then    ! no color insertions for photon emitter
      extcombs = 0
      nextcombs = 1
    else
      res = 0
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    Lsc%type=12
    Lsc%emitter=emitter
    Lsc%mom=polvect
    Lsc%nextcombs=nextcombs
    allocate(Lsc%extcombs(nextcombs))
    Lsc%extcombs = extcombs
    allocate(Lsc%rescc(0:n_external(id)*(n_external(id)+1)/2+1))
    call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lsc, acc)
    do j = 1, subprocess%n_particles
      res(j) = Lsc%rescc(extcombs(j))
    end do
    deallocate(Lsc%extcombs)
    deallocate(Lsc%rescc)
  end subroutine evaluate_loopsc


  subroutine evaluate_loopsc_c(id, pp, emitter, polvect, res) bind(c,name="ol_evaluate_loopsc")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id)), polvect(4)
    real(c_double), intent(out) :: res(n_external(id))
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id)), f_polvect(4)
    real(DREALKIND) :: f_res(n_external(id))
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    f_polvect = polvect
    call evaluate_loopsc(f_id, f_pp(0:3,:), f_emitter, f_polvect, f_res)
    res = f_res
   end subroutine evaluate_loopsc_c



  subroutine evaluate_scpowheg(id, psp, emitter, res, resmunu)
   ! Spin correlated tree matrix elements in POWHEG convention
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: squared born matrix element
   ! [out] res(4:4): array with the spin correlated born matrix element B^(mu,nu)
   ! B^(mu,nu) = sum_l,k M(l)M(k) epsilon_l^mu* epsilon_k^nu
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, resmunu(4,4)
    type(process_handle) :: subprocess
    real(DREALKIND) :: m2cc(0:n_external(id)*(n_external(id)+1)/2+1)
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 0)) then
      call ol_fatal('evaluate: scpowheg routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call tree_parameters_flush()
    call subprocess%tree(psp, m2cc, -emitter, &
      & [0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND, 0._/**/DREALKIND], &
      &  1, [0], resmunu)
    res = m2cc(0)
  end subroutine evaluate_scpowheg


  subroutine evaluate_scpowheg_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_scpowheg")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_scpowheg(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_scpowheg_c


  subroutine evaluate_stensor_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_stensor")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_scpowheg(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_stensor_c


  subroutine evaluate_sctensor(id, psp, emitter, res, resmunu)
   ! Colour-Spin correlated tree tensors
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [in] emitter: j
   ! [out] res: squared born matrix element
   ! [out] res(N:4:4): array with the spin correlated born matrix element B^(jk)(mu,nu)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, resmunu(n_external(id),4,4)
    call stop_invalid_id(id)
    call ol_fatal('evaluate: sctensor routine not available for process ' // trim(to_string(id)))
  end subroutine evaluate_sctensor


  subroutine evaluate_sctensor_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_sctensor")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(n_external(id)*16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_sctensor(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/n_external(id)*4*4/))
  end subroutine evaluate_sctensor_c


  subroutine evaluate_sctensor2(id, psp, emitter, res, resmunu)
   ! Colour-Spin correlated tree tensors
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [in] emitter: j
   ! [out] res: squared born matrix element
   ! [out] res(N:4:4): array with the spin correlated born matrix element B^(jk)(mu,nu)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, resmunu(n_external(id),4,4)
    call stop_invalid_id(id)
    call ol_fatal('evaluate: sctensor2 routine not available for process ' // trim(to_string(id)))
  end subroutine evaluate_sctensor2


  subroutine evaluate_sctensor2_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_sctensor2")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(n_external(id)*16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_sctensor2(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/n_external(id)*4*4/))
  end subroutine evaluate_sctensor2_c


  subroutine evaluate_loopsctensor(id, psp, emitter, m2l0, m2l1, resmunu)
   ! Colour-Spin correlated tree tensors
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [in] emitter: j
   ! [out] res: squared born matrix element
   ! [out] res(N:4:4): array with the spin correlated born matrix element B^(jk)(mu,nu)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), resmunu(n_external(id),4,4)
    call stop_invalid_id(id)
    call ol_fatal('evaluate: loopsctensor routine not available for process ' // trim(to_string(id)))
  end subroutine evaluate_loopsctensor


  subroutine evaluate_loopsctensor_c(id, pp, emitter, m2l0, m2l1, resmunu) bind(c,name="ol_evaluate_loopsctensor")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), resmunu(n_external(id)*16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_loopsctensor(f_id, f_pp(0:3,:), f_emitter, f_m2l0, f_m2l1, f_resmunu)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    resmunu = reshape(f_resmunu,(/n_external(id)*4*4/))
  end subroutine evaluate_loopsctensor_c


  subroutine evaluate_scpowheg2(id, psp, emitter, res, resmunu)
   ! Spin correlated loop-squared matrix elements in POWHEG convention
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: squared loop-squared matrix element
   ! [out] res(4:4): array with the spin correlated loop^2 matrix element L^(mu,nu)
   ! L^(mu,nu) = sum_l,k Mloop(l)Mloop(k) epsilon_l^mu* epsilon_k^nu
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: res, resmunu(4,4)
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(process_handle) :: subprocess
    type(correlator) :: Lmunu
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 2)) then
      call ol_fatal('evaluate: scpowheg2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    if (subprocess%extid(emitter) /= 21 .and. subprocess%extid(emitter) /= 22 .and. subprocess%extid(emitter) /= 0) then
      Lmunu%type=0
      call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lmunu, acc)
      res = m2l2(0)
      resmunu = 0
      return
    else
      Lmunu%type=3
      Lmunu%emitter=emitter
      call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lmunu, acc)
      res = m2l2(0)
      resmunu = Lmunu%resmunu
    end if
  end subroutine evaluate_scpowheg2


  subroutine evaluate_scpowheg2_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_scpowheg2")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_scpowheg2(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_scpowheg2_c


  subroutine evaluate_stensor2_c(id, pp, emitter, res, resmunu) bind(c,name="ol_evaluate_stensor2")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_scpowheg2(f_id, f_pp(0:3,:), f_emitter, f_res, f_resmunu)
    res = f_res
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_stensor2_c


  subroutine evaluate_loopscpowheg(id, psp, emitter, m2l0, m2l1, resmunu)
   ! Spin correlated loop x Born matrix elements in POWHEG convention
   ! [in] id: process id as set by register_process
   ! [in] psp: phase space point
   ! [out] res: squared loop-squared matrix element
   ! [out] res(4:4): array with the spin correlated loop^2 matrix element L^(mu,nu)
   ! L^(mu,nu) = sum_l,k Mloop(l)Mborn(k) epsilon_l^mu* epsilon_k^nu
    use ol_generic, only: to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    implicit none
    integer, intent(in) :: id, emitter
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), resmunu(4,4)
    real(DREALKIND) :: ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(process_handle) :: subprocess
    type(correlator) :: Lmunu
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loopscpowheg routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    if (subprocess%extid(emitter) /= 21 .and. subprocess%extid(emitter) /= 22 .and. subprocess%extid(emitter) /= 0) then
      Lmunu%type=0
      call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lmunu, acc)
      resmunu = 0
      return
    else
      Lmunu%type=13
      Lmunu%emitter=emitter
      call evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, Lmunu, acc)
      resmunu = Lmunu%resmunu
    end if
  end subroutine evaluate_loopscpowheg


  subroutine evaluate_loopscpowheg_c(id, pp, emitter, m2l0, m2l1, resmunu) bind(c,name="ol_evaluate_loopscpowheg")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_loopscpowheg(f_id, f_pp(0:3,:), f_emitter, f_m2l0, f_m2l1, f_resmunu)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_loopscpowheg_c


  subroutine evaluate_loopstensor_c(id, pp, emitter, m2l0, m2l1, resmunu) bind(c,name="ol_evaluate_loopstensor")
    implicit none
    integer(c_int), value :: id, emitter
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), resmunu(16)
    integer :: f_id, f_emitter
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_resmunu(4,4)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    f_emitter = emitter
    call evaluate_loopscpowheg(f_id, f_pp(0:3,:), f_emitter, f_m2l0, f_m2l1, f_resmunu)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    resmunu = reshape(f_resmunu,(/4*4/))
  end subroutine evaluate_loopstensor_c


  recursive subroutine evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
    use ol_stability
    use ol_generic, only: to_string
    use ol_parameters_decl_/**/DREALKIND, only: add_associated_ew, alpha_QCD
    use ol_loop_parameters_decl_/**/DREALKIND, only: IR_is_on, loop_parameters_status, muren_unscaled, CT_is_on, R2_is_on, TP_is_on
    use ol_loop_parameters_decl_/**/DREALKIND, only: stability_mode,a_switch,a_switch_rescue,redlib_qp,use_bubble_vertex
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4)
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: m2l0ew, m2l1ew(0:2), ir1ew(0:2), m2l2ew(0:4), ir2ew(0:4), accew
    real(DREALKIND) :: last_m2ct, new_m2ct, muren_bak
    integer :: IR_is_on_bak
    type(process_handle)  :: subprocess, subprocess_replace
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1) .and. subprocess%replace_loop==0) then
      call ol_fatal('evaluate: loop routine not available for process ' // trim(to_string(id)))
      return
    end if
    if (subprocess%stability_mode/=-1)  call set_if_modified(stability_mode,subprocess%stability_mode)
    if (subprocess%a_switch/=-1)        call set_if_modified(a_switch,subprocess%a_switch)
    if (subprocess%a_switch_rescue/=-1) call set_if_modified(a_switch_rescue,subprocess%a_switch_rescue)
    if (subprocess%redlib_qp/=-1)       call set_if_modified(redlib_qp,subprocess%redlib_qp)
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    if (IR_is_on == 5) then
      ! Return Born (from tree amplitude) as vamp (for debug)
      call evaluate_tree(id,psp,m2l0)
      m2l1 = m2l0
      m2l0 = 1.
      return
    else if (IR_is_on == 6) then
      ! Return zero (for debug)
      call evaluate_tree(id,psp,m2l0)
      m2l1 = 0
      return
    end if

    if (subprocess%replace_loop > 0) then
      call evaluate_full(subprocess%replace_loop, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
      call evaluate_tree(id,psp,m2l0)
    else
      call subprocess%loop(psp, m2l0, m2l1, ir1, m2l2, ir2)
    end if
    acc = last_relative_deviation

    if (last_from_cache .and. alpha_QCD /= subprocess%last_alpha_QCD) then
      call ol_msg(3, "evaluate_full: same psp, alpha_QCD changed: rescale.")
      m2l0 = m2l0*(alpha_QCD/subprocess%last_alpha_QCD)**subprocess%qcd_powers(1)
      m2l1 = m2l1*(alpha_QCD/subprocess%last_alpha_QCD)**(subprocess%qcd_powers(1)+subprocess%qcd_powers(2))
      ir1  = ir1*(alpha_QCD/subprocess%last_alpha_QCD)**(subprocess%qcd_powers(1)+subprocess%qcd_powers(2))
      m2l2 = m2l2*(alpha_QCD/subprocess%last_alpha_QCD)**subprocess%qcd_powers(1)
      ir2  = ir2*(alpha_QCD/subprocess%last_alpha_QCD)**subprocess%qcd_powers(1)
    end if


    if (last_from_cache .and. muren_unscaled /= subprocess%last_muren) then
      call ol_msg(3, "evaluate_full: same psp, muren changed: shift counterterm.")
      muren_bak = muren_unscaled
      call set_parameter("muren",subprocess%last_muren)
      call subprocess%ct(psp, m2l0, last_m2ct)
      call set_parameter("muren",muren_bak)
      call subprocess%ct(psp, m2l0, new_m2ct)
      m2l1(0) = m2l1(0)-last_m2ct+new_m2ct
    end if


    ! add associated one-loop ew
    if (add_associated_ew > 0 .and. subprocess%associated_ew > 0) then
      IR_is_on_bak = IR_is_on
      IR_is_on = 2
      call set_parameter("ew_renorm", 1)
      call evaluate_loop(subprocess%associated_ew, psp, m2l0ew, m2l1ew, accew)
      m2l1 = m2l1+m2l1ew
      acc = max(acc, accew)
      IR_is_on = IR_is_on_bak
      call set_parameter("ew_renorm", 0)
    else if (add_associated_ew > 0 .and. subprocess%associated_ew <= 0) then
      call ol_msg(2, "evaluate_full: associated EW loop library not loaded -> only QCD used.")
    end if

    ! add associated born ew
    if (add_associated_ew > 1 .and. subprocess%associated_born_1 > 0) then
      call evaluate_tree(subprocess%associated_born_1,psp,m2l0ew)
      m2l1(0) = m2l1(0)+m2l0ew
    else if (add_associated_ew > 1 .and. subprocess%associated_born_1 <= 0) then
      call ol_msg(2, "evaluate_full: associated EW born library not loaded -> only V used.")
    end if

    if (add_associated_ew > 2 .and. subprocess%associated_born_2 > 0) then
      call evaluate_tree(subprocess%associated_born_2,psp,m2l0ew)
      m2l1(0) = m2l1(0)+m2l0ew
    else if (add_associated_ew > 2 .and. subprocess%associated_born_2 <= 0) then
      call ol_msg(2, "evaluate_full: associated EW born library not loaded -> only V used.")
    end if

    if (IR_is_on == 3) then
      m2l1 = ir1  ! Return I-Operator as vamp (for debug)
    else if (IR_is_on == 4) then
      m2l1 = m2l0 ! Return Born as vamp (for debug)
    end if
    if (.not. last_from_cache) then
      process_handles(id)%last_psp = psp
      process_handles(id)%last_perm = subprocess%permutation
      process_handles(id)%last_pol = subprocess%pol
      process_handles(id)%last_zero = (m2l1(0) == 0 .and. m2l2(0) == 0)
      process_handles(id)%loop_parameters_status = loop_parameters_status
      process_handles(id)%last_muren = muren_unscaled
      process_handles(id)%last_alpha_QCD = alpha_QCD
    end if
    call ol_msg(5,"evaluate_full: " // trim(to_string(m2l0)) // " " // trim(to_string(m2l1(0))))
  end subroutine evaluate_full


  subroutine evaluate_full_c(id, pp, m2l0, m2l1, ir1, m2l2, ir2, acc) bind(c,name="ol_evaluate_full")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2), f_ir1(0:2), f_m2l2(0:4), f_ir2(0:4)
    real(DREALKIND) :: f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_full(f_id, f_pp(0:3,:), f_m2l0, f_m2l1, f_ir1, f_m2l2, f_ir2, f_acc)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    ir1  = f_ir1
    m2l2 = f_m2l2
    ir2  = f_ir2
    acc  = f_acc
  end subroutine evaluate_full_c


  subroutine evaluate_fullcr(id, psp, m2l0, m2l1, ir1, m2l2, ir2, cr, acc)
    use ol_stability
    use ol_generic, only: to_string, to_string
    use ol_data_types_/**/DREALKIND, only: correlator
    use ol_loop_parameters_decl_/**/DREALKIND, only: loop_parameters_status
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4)
    real(DREALKIND), intent(out) :: acc
    type(correlator) :: cr
    type(process_handle)  :: subprocess
    logical has_loopcc
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loop routine not available for process ' // trim(to_string(id)))
      return
    end if
    m2l1=0
    ir1=0
    m2l2=0
    ir2=0
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    if (cr%type > 10) call subprocess%loopcc(.true., has_loopcc)
    if (any(subprocess%last_psp /= psp) .or. &
      & any(subprocess%last_perm /= subprocess%permutation) .or. &
      & any(subprocess%last_pol /= subprocess%pol) .or. &
      & subprocess%loop_parameters_status /= loop_parameters_status .or. &
        .not. has_loopcc &
      ) then
      call ol_msg(3, "me-cache for correlators: " // trim(subprocess%process_name) &
                 &  // trim(to_string(subprocess%permutation,.true.)) // ' reevaluate' )
      call evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)  ! fill colour/helicity vector cache
      if (m2l1(0) == 0 .and. m2l2(0) == 0) then
        cr%rescc = 0
        cr%resmunu = 0
        return
      end if
    else
      call ol_msg(3, "me-cache for correlators: " // trim(subprocess%process_name) &
                 &  // trim(to_string(subprocess%permutation,.true.)) // ' taken from the cache' )
      if (subprocess%last_zero) then
        cr%rescc = 0
        cr%resmunu = 0
        return ! return for unstable points
      end if
    end if
    if (cr%type < 1) then
      return
    else if (cr%type < 14) then
      call subprocess%loopcr(psp, m2l0, m2l1(0), m2l2(0), &
            cr%type, cr%emitter, cr%nextcombs, cr%extcombs, cr%mom, cr%rescc, cr%resmunu)
    else
      call ol_fatal("evaluate_fullcr: correlator type not available.")
    end if
    acc = last_relative_deviation
  end subroutine evaluate_fullcr


  subroutine evaluate_loop(id, psp, m2l0, m2l1, acc)
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2)
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: ir1(0:2), m2l2(0:4), ir2(0:4)
    call evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
  end subroutine evaluate_loop


  subroutine evaluate_loop_c(id, pp, m2l0, m2l1, acc) bind(c,name="ol_evaluate_loop")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2)
    real(DREALKIND) :: f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loop(f_id, f_pp(0:3,:), f_m2l0, f_m2l1, f_acc)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    acc  = f_acc
  end subroutine evaluate_loop_c


  subroutine evaluate_loopbare(id, psp, m2l0, m2l1, acc)
    use ol_loop_parameters_decl_/**/DREALKIND, only: CT_is_on, R2_is_on, IR_is_on, polecheck_is
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2l1(0:2)
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: ir1(0:2), m2l2(0:4), ir2(0:4), m2ct(0:2)
    integer :: CT_is_on_bak, R2_is_on_bak, IR_is_on_bak
    CT_is_on_bak = CT_is_on
    R2_is_on_bak = R2_is_on
    IR_is_on_bak = IR_is_on
    call set_parameter("ct_on", 1)
    call set_parameter("ir_on", 1)
    call set_parameter("r2_on", 0)
    call evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
    call evaluate_loopct(id, psp, m2l0, m2ct)
    m2l1 = m2l1 - m2ct
    call set_parameter("ct_on", CT_is_on_bak)
    call set_parameter("ir_on", IR_is_on_bak)
    call set_parameter("r2_on", R2_is_on_bak)
  end subroutine evaluate_loopbare


  subroutine evaluate_loopbare_c(id, pp, m2l0, m2l1, acc) bind(c,name="ol_evaluate_loopbare")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2l1(0:2), acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2l1(0:2)
    real(DREALKIND) :: f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loopbare(f_id, f_pp(0:3,:), f_m2l0, f_m2l1, f_acc)
    m2l0 = f_m2l0
    m2l1 = f_m2l1
    acc  = f_acc
  end subroutine evaluate_loopbare_c


  subroutine evaluate_loop2(id, psp, res, acc)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: res
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4)
    if (.not. btest(process_handles(id)%content, 2)) then
      call ol_fatal('evaluate: loop^2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    call evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
    res = m2l2(0)
  end subroutine evaluate_loop2


  subroutine evaluate_loop2_c(id, pp, res, acc) bind(c,name="ol_evaluate_loop2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loop2(f_id, f_pp(0:3,:), f_res, f_acc)
    res = f_res
    acc  = f_acc
  end subroutine evaluate_loop2_c


  subroutine evaluate_iop2(id, psp, res, resir, acc)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: res, resir(0:2)
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: m2l0, m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4)
    if (.not. btest(process_handles(id)%content, 2)) then
      call ol_fatal('evaluate: loop^2 routine not available for process ' // trim(to_string(id)))
      return
    end if
    call evaluate_full(id, psp, m2l0, m2l1, ir1, m2l2, ir2, acc)
    res = m2l2(0)
    resir(0) = ir2(0)
    resir(1) = ir2(1)
    resir(2) = ir2(2)
  end subroutine evaluate_iop2


  subroutine evaluate_iop2_c(id, pp, res, resir, acc) bind(c,name="ol_evaluate_iop2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res, resir(0:2), acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res, f_resir(0:2), f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_iop2(f_id, f_pp(0:3,:), f_res, f_resir, f_acc)
    res = f_res
    resir = f_resir
    acc  = f_acc
  end subroutine evaluate_iop2_c


  subroutine evaluate_associated(id, psp, level, m2l0)
    use ol_stability
    use ol_generic, only: to_string
    use ol_parameters_decl_/**/DREALKIND, only: add_associated_ew
    use ol_loop_parameters_decl_/**/DREALKIND, only: IR_is_on
    implicit none
    integer, intent(in) :: id, level
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0
    real(DREALKIND) :: m2l0ew, m2l1ew(0:2), ir1ew(0:2), m2l2ew(0:4), ir2ew(0:4), accew
    integer :: IR_is_on_bak
    type(process_handle)  :: subprocess, subprocessew
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: loop routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()

    m2l0 = 0
    if (level .eq. 1 .and. subprocess%associated_ew > 0) then ! add associated one-loop ew
      IR_is_on_bak = IR_is_on
      IR_is_on = 2
      call set_parameter("ew_renorm", 1)
      call evaluate_loop(subprocess%associated_ew, psp, m2l0ew, m2l1ew, accew)
      m2l0 = m2l1ew(0)
      IR_is_on = IR_is_on_bak
      call set_parameter("ew_renorm", 0)
    else if (level .eq. 1 .and. subprocess%associated_ew <= 0) then
      call ol_msg(2, "evaluate_associated: associated EW loop library not loaded -> only QCD used.")
    else if (level .eq. 2 .and. subprocess%associated_born_1 > 0) then ! add associated born ew
      call evaluate_tree(subprocess%associated_born_1,psp,m2l0ew)
      m2l0 = m2l0ew
    else if (level .eq. 2 .and. subprocess%associated_born_1 <= 0) then
      call ol_msg(2, "evaluate_associated: associated EW born library not loaded -> only V used.")
    else if (level .eq. 3 .and. subprocess%associated_born_2 > 0) then
      call evaluate_tree(subprocess%associated_born_2,psp,m2l0ew)
      m2l0 = m2l0ew
    else if (level .eq. 3 .and. subprocess%associated_born_2 <= 0) then
      call ol_msg(2, "evaluate_associated: associated EW born library not loaded -> only V used.")
    end if

    call ol_msg(5,"evaluate_associated: " // trim(to_string(m2l0)) )
  end subroutine evaluate_associated


  subroutine evaluate_associated_c(id, pp, level, res) bind(c,name="ol_evaluate_associated")
    implicit none
    integer(c_int), value :: id, level
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: res
    integer :: f_id, f_level
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_res
    f_id = id
    f_level = level
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_associated(f_id, f_pp(0:3,:), f_level, f_res)
    res = f_res
  end subroutine evaluate_associated_c


  subroutine evaluate_ct(id, psp, m2l0, m2ct)
    use ol_parameters_decl_/**/DREALKIND, only: add_associated_ew
    use ol_loop_parameters_decl_/**/DREALKIND, only: CT_is_on, R2_is_on, TP_is_on, do_ew_renorm,use_bubble_vertex
    use ol_init, only: set_if_modified, parameters_flush
    use ol_stability
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    integer :: CT_on_bak, R2_on_bak, TP_on_bak, use_bubble_vertex_bak
    real(DREALKIND), intent(out) :: m2l0, m2ct
    real(DREALKIND) :: m2l0ew, m2ctew
    type(process_handle)  :: subprocess, subprocessew
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: ct routine not available for process ' // trim(to_string(id)))
      return
    end if
    ! For compatibility: signal ctamp2 that this is the (new) interface (CT_is_on=2),
    ! where CT_is_on is set here, not in ctamp2base.
    ! If the old interface is used with a new process, the process will thus notice
    ! that it has to set CT_is_on itself.
    CT_on_bak = CT_is_on
    R2_on_bak = R2_is_on
    TP_on_bak = TP_is_on
    use_bubble_vertex_bak = use_bubble_vertex
    call set_if_modified(CT_is_on, 2)
    call set_if_modified(R2_is_on, 0)
    call set_if_modified(TP_is_on, 0)
    call set_if_modified(use_bubble_vertex, 0)
    ! call parameters_flush() is done in ctamp2base,
    ! because only there we know if ew renormalisation needs to be activated.
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call subprocess%ct(psp, m2l0, m2ct)
    if (add_associated_ew == 1 .and. subprocess%associated_ew > 0) then
      subprocessew = process_handles(subprocess%associated_ew)
      if (.not. btest(subprocessew%content, 1)) then
        call ol_fatal('evaluate: loop routine not available for associated process ' // trim(to_string(subprocess%associated_ew)))
        call set_if_modified(CT_is_on, CT_on_bak)
        call set_if_modified(R2_is_on, R2_on_bak)
        call set_if_modified(TP_is_on, TP_on_bak)
        call set_if_modified(use_bubble_vertex, use_bubble_vertex_bak)
        call parameters_flush()
        return
      end if
      n_scatt = subprocess%n_in
      call subprocessew%set_permutation(subprocessew%permutation)
      if (any(subprocessew%photon_id /= 0)) call subprocessew%set_photons(subprocessew%photon_id)
      call set_if_modified(do_ew_renorm, 1)
      call parameters_flush()
      call subprocessew%ct(psp, m2l0ew, m2ctew)
      m2ct = m2ct + m2ctew
      call set_if_modified(do_ew_renorm, 0)
    else if (add_associated_ew == 1 .and. subprocess%associated_ew <= 0) then
      call ol_error("evaluate_ct: associated EW library not loaded -> only QCD used.")
    end if
    call set_if_modified(CT_is_on, CT_on_bak)
    call set_if_modified(R2_is_on, R2_on_bak)
    call set_if_modified(TP_is_on, TP_on_bak)
    call set_if_modified(use_bubble_vertex, use_bubble_vertex_bak)
    call parameters_flush()
  end subroutine evaluate_ct


  subroutine evaluate_ct_c(id, pp, m2l0, m2ct) bind(c,name="ol_evaluate_ct")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2ct
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2ct
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_ct(f_id, f_pp(0:3,:), f_m2l0, f_m2ct)
    m2l0 = f_m2l0
    m2ct = f_m2ct
  end subroutine evaluate_ct_c


  subroutine evaluate_loopct(id, psp, m2l0, m2ct)
    use ol_loop_parameters_decl_/**/DREALKIND, only: CT_is_on, IR_is_on
    use ol_loop_parameters_decl_/**/DREALKIND, only: de1_UV, de1_IR, de2_i_IR
    use ol_stability
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2ct(0:2)
    real(DREALKIND) :: m2ct_000, m2ct_110, polescale
    real(DREALKIND) :: poleUV1bak, poleIR1bak, poleIR2bak
    integer CT_on_bak
    call stop_invalid_id(id)
    if (error > 1) return
    ! remember original parameters
    poleUV1bak = de1_UV
    poleIR1bak = de1_IR
    poleIR2bak = de2_i_IR
    CT_on_bak = CT_is_on
    call get_parameter("polescale", polescale)

    ! two calls to CT amplitude with different UV & IR poles
    call set_parameter("ct_on", 1)
    call set_parameter("pole_uv", polescale)
    call set_parameter("pole_ir1", polescale)
    call set_parameter("pole_ir2", 0)
    call evaluate_ct(id, psp, m2l0, m2ct_110)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", 0)
    call evaluate_ct(id, psp, m2l0, m2ct_000)

    m2ct(0)=m2ct_000
    m2ct(1)=(m2ct_110-m2ct_000)/polescale
    m2ct(2)=0

    ! restore original parameters
    call set_parameter("pole_uv", poleUV1bak)
    call set_parameter("pole_ir1", poleIR1bak)
    call set_parameter("pole_ir2", poleIR2bak)
    call set_parameter("ct_on", CT_on_bak)
  end subroutine evaluate_loopct


  subroutine evaluate_loopct_c(id, pp, m2l0, m2ct) bind(c,name="ol_evaluate_loopct")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2ct(0:2)
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2ct(0:2)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_loopct(f_id, f_pp(0:3,:), f_m2l0, f_m2ct)
    m2l0 = f_m2l0
    m2ct = f_m2ct
  end subroutine evaluate_loopct_c


  subroutine evaluate_iop(id, psp, m2l0, m2iop)
    use ol_generic, only: to_string
    use ol_loop_parameters_decl_/**/DREALKIND, only: IR_is_on
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2iop(0:2)
    real(DREALKIND) :: m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4), acc
    type(process_handle)  :: subprocess
    integer :: IR_is_on_bak
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1) .or. subprocess%amplitude_type < 10) then
      call ol_fatal("evaluate: iop routine not available for process " // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call tree_parameters_flush()
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    IR_is_on_bak = IR_is_on
    call set_parameter("ir_on", 1)
    call subprocess%iop(psp, m2l0, m2iop)
    call set_parameter("ir_on", IR_is_on_bak)
  end subroutine evaluate_iop


  subroutine evaluate_iop_c(id, pp, m2l0, m2iop) bind(c,name="ol_evaluate_iop")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2iop(0:2)
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2iop(0:2)
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_iop(f_id, f_pp(0:3,:), f_m2l0, f_m2iop)
    m2l0 = f_m2l0
    m2iop = f_m2iop
  end subroutine evaluate_iop_c


  subroutine evaluate_schsf(id, psp, m2l0, m2schsf) bind(c,name="ol_evaluate_schsf")
    use ol_stability
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2schsf
    type(process_handle)  :: subprocess
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: Spin-correlated hard scattering factor not available for process ' // trim(to_string(id)))
      return
    end if
    if(.not. associated(subprocess%schsf)) then
      call ol_fatal('evaluate: Spin-correlated hard scattering factor not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call subprocess%schsf(psp, m2l0, m2schsf)
  end subroutine evaluate_schsf



  subroutine evaluate_r2(id, psp, m2l0, m2ct)
    use ol_parameters_decl_/**/DREALKIND, only: add_associated_ew
    use ol_loop_parameters_decl_/**/DREALKIND, only: CT_is_on, R2_is_on, TP_is_on, do_ew_renorm, use_bubble_vertex
    use ol_init, only: set_if_modified, parameters_flush
    use ol_stability
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    integer :: CT_on_bak, R2_on_bak, TP_on_bak, use_bubble_vertex_bak
    real(DREALKIND), intent(out) :: m2l0, m2ct
    real(DREALKIND) :: m2l0ew, m2ctew
    type(process_handle)  :: subprocess, subprocessew
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1)) then
      call ol_fatal('evaluate: ct routine not available for process ' // trim(to_string(id)))
      return
    end if
    ! For compatibility: signal ctamp2 that this is the (new) interface (CT_is_on=-1),
    ! where CT_is_on is set here, not in ctamp2.
    CT_on_bak = CT_is_on
    R2_on_bak = R2_is_on
    TP_on_bak = TP_is_on
    use_bubble_vertex_bak = use_bubble_vertex
    call set_if_modified(CT_is_on, -1)
    call set_if_modified(R2_is_on, 1)
    call set_if_modified(TP_is_on, 0)
    call set_if_modified(use_bubble_vertex, 0)
    ! call parameters_flush() is done in ctamp2base,
    ! because only there we know if ew renormalisation needs to be activated.
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    call subprocess%pol_init(subprocess%pol)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call subprocess%ct(psp, m2l0, m2ct)
    call set_if_modified(CT_is_on, CT_on_bak)
    call set_if_modified(R2_is_on, R2_on_bak)
    call set_if_modified(TP_is_on, TP_on_bak)
    call set_if_modified(use_bubble_vertex, use_bubble_vertex_bak)
    call parameters_flush()
  end subroutine evaluate_r2


  subroutine evaluate_r2_c(id, pp, m2l0, m2ct) bind(c,name="ol_evaluate_r2")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2ct
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2ct
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_r2(f_id, f_pp(0:3,:), f_m2l0, f_m2ct)
    m2l0 = f_m2l0
    m2ct = f_m2ct
  end subroutine evaluate_r2_c


  subroutine evaluate_pt(id, psp, m2l0, m2pt, m2l1)
    use ol_stability
    use ol_generic, only: to_string
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2pt, m2l1
    type(process_handle)  :: subprocess
    call stop_invalid_id(id)
    if (error > 1) return
    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 3)) then
      call ol_fatal('evaluate: pt routine not available for process ' // trim(to_string(id)))
      return
    end if
    n_scatt = subprocess%n_in
    call subprocess%set_permutation(subprocess%permutation)
    if (any(subprocess%photon_id /= 0)) call subprocess%set_photons(subprocess%photon_id)
    call parameters_flush()
    call subprocess%pt(psp, m2l0, m2pt, m2l1)
  end subroutine evaluate_pt


  subroutine evaluate_pt_c(id, pp, m2l0, m2pt, m2l1) bind(c,name="ol_evaluate_pt")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2pt, m2l1
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2pt, f_m2l1
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_pt(f_id, f_pp(0:3,:), f_m2l0, f_m2pt, f_m2l1)
    m2l0 = f_m2l0
    m2pt = f_m2pt
    m2l1 = f_m2l1
  end subroutine evaluate_pt_c


  subroutine evaluate_poles(id, psp, m2l0, m2bare, m2ct, m2ir, m2sum)
    use ol_init, only: set_if_modified, parameters_flush
    use ol_generic, only: to_string
    use ol_loop_parameters_decl_/**/DREALKIND, only: CT_is_on, IR_is_on
    use ol_loop_parameters_decl_/**/DREALKIND, only: de1_UV, de1_IR, de2_i_IR
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(out) :: m2l0, m2bare(0:2), m2ct(0:2), m2ir(0:2), m2sum(0:2)
    real(DREALKIND) :: m2l1(0:2), ir1(0:2), m2l2(0:4), ir2(0:4)
    real(DREALKIND) :: m2l1_100(0:2), m2l1_010(0:2), m2l1_001(0:2), m2l1_000(0:2)
    real(DREALKIND) :: m2ct_100, m2ct_010, m2ct_000
    real(DREALKIND) :: acc, polescale
    real(DREALKIND) :: poleUV1bak, poleIR1bak, poleIR2bak
    integer polecheck, CT_on_bak, IR_on_bak
    call stop_invalid_id(id)
    if (error > 1) return
    ! remember original parameters
    poleUV1bak = de1_UV
    poleIR1bak = de1_IR
    poleIR2bak = de2_i_IR
    CT_on_bak = CT_is_on
    IR_on_bak = IR_is_on
    call get_parameter("polescale", polescale)
    call get_parameter("polecheck", polecheck)
    if (polecheck > 0) then
      call ol_msg("evaluate_poles: truepoles>0 not allowed.")
      return
    end if

    ! four calls to bare amplitude with different UV & IR poles
    call set_parameter("ct_on", 0)
    call set_parameter("pole_uv", polescale)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", 0)
    call set_parameter("ir_on", 0)
    call evaluate_full(id, psp, m2l0, m2l1_100, ir1, m2l2, ir2, acc)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", polescale)
    call set_parameter("pole_ir2", 0)
    call evaluate_full(id, psp, m2l0, m2l1_010, ir1, m2l2, ir2, acc)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", polescale)
    call evaluate_full(id, psp, m2l0, m2l1_001, ir1, m2l2, ir2, acc)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", 0)
    call set_parameter("ir_on", 1)
    call evaluate_full(id, psp, m2l0, m2l1_000, ir1, m2l2, ir2, acc)

    m2bare(0)=(m2l1_100(0)-m2l1_000(0))/polescale
    m2bare(1)=(m2l1_010(0)-m2l1_000(0))/polescale
    m2bare(2)=(m2l1_001(0)-m2l1_000(0))/polescale

    ! three calls to CT amplitude with different UV & IR poles
    call set_parameter("ct_on", 1)
    call set_parameter("pole_uv", polescale)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", 0)
    call evaluate_ct(id, psp, m2l0, m2ct_100)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", polescale)
    call set_parameter("pole_ir2", 0)
    call evaluate_ct(id, psp, m2l0, m2ct_010)
    call set_parameter("pole_uv", 0)
    call set_parameter("pole_ir1", 0)
    call set_parameter("pole_ir2", 0)
    call evaluate_ct(id, psp, m2l0, m2ct_000)

    m2ct(0)=(m2ct_100-m2ct_000)/polescale
    m2ct(1)=(m2ct_010-m2ct_000)/polescale
    m2ct(2)=0

    m2ir(0)=0
    m2ir(1)=ir1(1)
    m2ir(2)=ir1(2)

    m2sum=m2bare+m2ct+m2ir

    ! restore original parameters
    call set_parameter("pole_uv", poleUV1bak)
    call set_parameter("pole_ir1", poleIR1bak)
    call set_parameter("pole_ir2", poleIR2bak)
    call set_parameter("ct_on", CT_on_bak)
    call set_parameter("ir_on", IR_on_bak)
  end subroutine evaluate_poles


  subroutine evaluate_poles_c(id, pp, m2l0, m2bare, m2ct, m2ir, m2sum) bind(c,name="ol_evaluate_poles")
    implicit none
    integer(c_int), value :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(out) :: m2l0, m2bare(0:2), m2ct(0:2), m2ir(0:2), m2sum(0:2)
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_m2l0, f_m2bare(0:2), f_m2ct(0:2), f_m2ir(0:2), f_m2sum(0:2)
    real(DREALKIND) :: f_acc
    f_id = id
    call stop_invalid_id(f_id) ! needed because of reshape
    if (error > 1) return
    f_pp = reshape(pp, [5,process_handles(id)%n_particles])
    call evaluate_poles(f_id, f_pp(0:3,:), f_m2l0, f_m2bare, f_m2ct, f_m2ir, f_m2sum)
    m2l0 = f_m2l0
    m2bare = f_m2bare
    m2ct   = f_m2ct
    m2ir   = f_m2ir
    m2sum  = f_m2sum
  end subroutine evaluate_poles_c



  subroutine check_poles(id, psp_in)
    use ol_generic, only: to_string
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in), optional  :: psp_in(:,:)
    real(DREALKIND) :: psp(0:3,n_external(int(id)))
    real(DREALKIND) :: m2l0, m2bare(0:2), m2ct(0:2), m2ir(0:2), m2sum(0:2)
    real(DREALKIND) :: energy = 2000
    character(len=150) :: output
    type(process_handle)  :: subprocess
    integer :: polecheck

    subprocess = process_handles(id)
    if (.not. btest(subprocess%content, 1) .and. subprocess%replace_loop==0) then
      return
    end if

    call get_parameter("polecheck", polecheck)
    if (polecheck > 0) then
      call ol_msg("check_poles: truepoles>0 not allowed.")
      return
    end if

    if(present(psp_in)) then
      psp = psp_in
    else
      call phase_space_point(id, energy, psp)
    end if

    call evaluate_poles(id, psp, m2l0, m2bare, m2ct, m2ir, m2sum)
2000 format(4e25.16, 4e25.16, 4e25.16)
    call ol_msg("== UV and IR Poles for process: " // trim(to_string(id)) // " ==")
    write(output, 2000) m2l0
    call ol_msg("Wtree       = " // output)
    call ol_msg("                 UVeps1                   IReps1                   IReps2")
    write(output, 2000) m2bare(0), m2bare(1), m2bare(2)
    call ol_msg("Wbare       = "// output)
    write(output, 2000) m2ct(0), m2ct(1), m2ct(2)
    call ol_msg("Wct         = "// output)
    write(output, 2000) m2ir(0), m2ir(1), m2ir(2)
    call ol_msg("WI-operator = "// output)
    write(output, 2000) m2sum(0), m2sum(1), m2sum(2)
    call ol_msg("W1-loop     = "// output)


  end subroutine check_poles


end module openloops
