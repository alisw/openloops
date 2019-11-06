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


module openloops_blha
!
! This is the BLHA interface of OpenLoops.
! Load just this module to use OpenLoops via the BLHA.
!
  use KIND_TYPES, only: DREALKIND
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_double, c_null_char
  use ol_init, only: set_parameter, set_init_error_fatal
  use ol_parameters_decl_/**/DREALKIND, only: max_parameter_length
  use ol_debug, only: ol_error, ol_msg, ol_fatal
  use openloops, only: ol_printparameter
  use ol_data_types_/**/DREALKIND, only: correlator
  implicit none
  private
  ! BLHA interface
  public :: olp_setparameter, olp_evalsubprocess, olp_evalsubprocess2
  public :: olp_info, olp_printparameter, olp_start, olp_start_line

  ! if 0: write OLP answer only to file; 1: also to stdout; 2: only to stdout
  integer :: stdout_contract = 0

  character (len=max_parameter_length), allocatable :: blha_answer(:)

  type(correlator) corr

  type flag
    integer :: InterfaceVersion ! 1=BLHA1, 2=BLHA2
    integer :: CorrectionType ! 1=QCD, 2=EW
    integer :: AmplitudeType ! 1=Tree, 2=ccTree, 3=scTree, 4=scTree_polvect, 11=Loop, 12=LoopInduced
    integer :: PoleCheck ! 0/1
    character(max_parameter_length) :: answer_file_name ! Answer file as string
  end type flag

  type(flag) flags


  interface olp_printparameter
    module procedure ol_printparameter
  end interface olp_printparameter


  contains

  subroutine olp_evalsubprocess2(id, psp, mu, rval, acc)
    ! Fortran BLHA-like olp_evalsubprocess2.
    ! Evaluate matrix elements for the process corresponding
    ! to process_handles(id) (including crossing).
    use openloops, only: amplitudetype, rval_size, n_external, &
                   evaluate_tree, evaluate_loop, evaluate_cc, evaluate_sc, &
                   evaluate_loop2, evaluate_cc2, evaluate_sc2
    implicit none
    integer, intent(in) :: id
    real(DREALKIND), intent(in) :: psp(:,:)
    real(DREALKIND), intent(in) :: mu
    real(DREALKIND), intent(out) :: rval(:)
    real(DREALKIND), intent(out) :: acc
    real(DREALKIND) :: m2l0, m2l1(0:2)
    rval = 0
    acc = 0
    call set_parameter("muren", mu)
    select case (amplitudetype(id))
      case (1) ! Tree
        call evaluate_tree(id, psp, rval(1))
      case (2) ! ccTree
        call evaluate_cc(id, psp, m2l0, rval(1:rval_size(n_external(id),2)), m2l1(0))
      case (3) ! scTree
        call ol_fatal("[OpenLoops] Error: spin correlations in BLHA notation are not implemented")
        return
      case (4) ! scTree_polvect
        call evaluate_sc2(id, psp, corr%emitter, corr%mom, rval(1:rval_size(n_external(id),4)))
      case (11) ! Loop
        call evaluate_loop(id, psp, m2l0, m2l1, acc)
        rval(1:4) = [m2l1(2), m2l1(1), m2l1(0), m2l0]
      case (12) ! LoopInduced
        call evaluate_loop2(id, psp, rval(1), acc)
      case (22) ! ccLoop^2
        call evaluate_cc2(id, psp, m2l0, rval(1:rval_size(n_external(id),2)), m2l1(0))
      case (24) ! scLoop^2_polvect
        call evaluate_sc2(id, psp, corr%emitter, corr%mom, rval(1:rval_size(n_external(id),4)))
      case default
        call ol_fatal("invalid amplitude type")
        return
    end select
  end subroutine olp_evalsubprocess2

  subroutine olp_scpolvec(emitter, mom)
    ! sets the emitter and polarization vector for the spin correlators
    implicit none
    integer, intent(in) :: emitter
    real(DREALKIND), intent(in) :: mom
      corr%emitter = emitter
      corr%mom = mom
  end subroutine olp_scpolvec

  subroutine olp_evalsubprocess(id, psp, mu, alpha_s, rval)
    ! Fortran BLHA-like olp_evalsubprocess routine.
    ! Wrapper to olp_evalsubprocess2 with alpha_s argument and without accuracy parameter
    implicit none
    integer, intent(in)  :: id
    real(DREALKIND), intent(in)  :: psp(:,:)
    real(DREALKIND), intent(in)  :: mu, alpha_s
    real(DREALKIND), intent(out) :: rval(:)
    real(DREALKIND)              :: acc
    call set_parameter("alphas", alpha_s)
    call olp_evalsubprocess2(id, psp, mu, rval, acc)
  end subroutine olp_evalsubprocess


  subroutine olp_setparameter(para, re, im, ierr)
    ! Fortran BLHA-like olp_setparameter routine.
    ! Must by flushed by parameters_flush() or start() to take effect
    ! (automatically done by all evaluate() variants).
    ! [in]  para: parameter name
    ! [in]  re: double precision value
    ! [in]  im: double precision value, ignored
    ! [out] ierr: 0=error, 1=ok, 2=ignored
    implicit none
    character(len=*), intent(in) :: para
    real(DREALKIND), intent(in) :: re, im
    integer, intent(out) :: ierr
    call set_parameter(para, re, ierr)
    ierr = ierr + 1
    if (ierr == 3) ierr = 0
  end subroutine olp_setparameter


  subroutine olp_setparameter_c(para, re, im, ierr) bind(c,name="OLP_SetParameter")
    ! BLHA OLP_SetParameter routine.
    ! C wrapper to olp_setparameter
    ! [in]  para: parameter name as a C string
    ! [in]  re: c_double value
    ! [in]  im: c_double value, ignored
    ! [out] ierr: 0=error, 1=ok, 2=ignored
    use ol_iso_c_utilities, only: c_f_string
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: para
    real(c_double), intent(in) :: re, im
    integer(c_int), intent(out) :: ierr
    character(len=max_parameter_name_length) :: f_para
    real(DREALKIND) :: f_re, f_im
    integer :: f_ierr
    f_re = re
    f_im = im
    call c_f_string(para, f_para, max_parameter_name_length)
    call OLP_SetParameter(trim(f_para), f_re, f_im, f_ierr)
    ierr = f_ierr
  end subroutine olp_setparameter_c


  subroutine olp_info(olp_name, olp_version, message)
    ! Fortran BLHA-like olp_info routine.
    ! [out] olp_name: "OpenLoops"
    ! [out] olp_version: OpenLoops version or svn revision (if unversioned) as a C string
    ! [out] message: TODO
    use ol_version, only: version, revision
    implicit none
    character(len=14), intent(out) :: olp_name, olp_version
    character(len=254), intent(out) :: message
    if (len(trim(version)) /= 0) then
      olp_version = version(9:) ! cut off "version "
    else
      olp_version = "SVN r" // revision
    end if
    olp_name = "OpenLoops"
    message = "blabla"
  end subroutine olp_info


  subroutine olp_evalsubprocess_c(id, pp, mu, alpha_s, rval) bind(c,name="OLP_EvalSubProcess")
    ! BLHA OLP_EvalSubProcess routine (version 1).
    ! C wrapper to olp_evalsubprocess
    ! [in]  C int* id: process id as set by register_process
    ! [in]  C double* pp: phase space point as a one-dimensional array
    !                     of size 5*n (E1,p1x,p1y,p1z,m1,E2,p2x,p2y,p2z,m2,...)
    !                     in 2 -> n-2 notation (p1+p2=p3+..+pn)
    ! [in]  C double* mu: renormalisation scale
    ! [in]  C double* alpha_s: strong coupling constant
    ! [out] C double* rval: array [loop_ir2, loop_ir1, loop_finite, born]
    use openloops, only: n_external, amplitudetype, rval_size, stop_invalid_id
    implicit none
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(in) :: mu, alpha_s
    real(c_double), intent(out) :: rval(rval_size(n_external(id), amplitudetype(id)))
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_mu, f_alpha_s
    real(DREALKIND) :: f_rval(rval_size(n_external(id), amplitudetype(id)))
    f_id = id
    call stop_invalid_id(f_id)
    f_pp = reshape(pp, [5,n_external(id)])
    f_mu = mu
    f_alpha_s = alpha_s
    call olp_evalsubprocess(f_id, f_pp(0:3,:), f_mu, f_alpha_s, f_rval)
    rval = f_rval
  end subroutine olp_evalsubprocess_c


  subroutine olp_evalsubprocess2_c(id, pp, mu, rval, acc) bind(c,name="OLP_EvalSubProcess2")
    ! BLHA OLP_EvalSubProcess2 routine (version 2).
    ! C wrapper to olp_evalsubprocess2
    ! [in]  C int* id: process id as set by register_process
    ! [in]  C double* pp: phase space point as a one-dimensional array
    !                     of size 5*n (E1,p1x,p1y,p1z,E2,p2x,p2y,p2z...)
    !                     in 2 -> n-2 notation (p1+p2=p3+..+pn)
    ! [in]  C double* mu: renormalisation scale
    ! [out] C double* rval: array with results (depending on amplitude_type)
    ! [out] C int* acc: accuracy, not available, always returns 0 (good)
    use openloops, only: n_external, amplitudetype, rval_size, stop_invalid_id
    implicit none
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: pp(5*n_external(id))
    real(c_double), intent(in) :: mu
    real(c_double), intent(out) :: rval(rval_size(n_external(id), amplitudetype(id)))
    real(c_double), intent(out) :: acc
    integer :: f_id
    real(DREALKIND) :: f_pp(0:4,n_external(id))
    real(DREALKIND) :: f_mu
    real(DREALKIND) :: f_rval(rval_size(n_external(id), amplitudetype(id)))
    real(DREALKIND) :: f_acc
    f_id = id
    call stop_invalid_id(f_id)
    f_pp = reshape(pp, [5,n_external(id)])
    f_mu = mu
    call olp_evalsubprocess2(f_id, f_pp(0:3,:), f_mu, f_rval, f_acc)
    rval = f_rval
    acc = f_acc
  end subroutine olp_evalsubprocess2_c


  subroutine olp_scpolvec_c(emitter, mom) bind(c,name="olp_scpolvec")
    ! sets the emitter and polarization vector for the spin correlators
    implicit none
    integer(c_int), intent(in) :: emitter
    real(c_double), intent(in) :: mom
    integer          :: f_emitter
    real(DREALKIND)  :: f_mom
    f_emitter = emitter
    f_mom = mom
    call olp_scpolvec(f_emitter, f_mom)
  end subroutine olp_scpolvec_c


  subroutine olp_info_c(olp_name, olp_version, message) bind(c,name="OLP_Info")
    ! BLHA OLP_Info routine.
    ! C wrapper to olp_info
    ! [out] olp_name: "OpenLoops" as a C string
    ! [out] olp_version: OpenLoops version or svn revision (if unversioned) as a C string
    ! [out] message as a C string
    implicit none
    character(kind=c_char), intent(out) :: olp_name(15), olp_version(15), message(255)
    character(len=14) :: f_olp_name, f_olp_version
    character(len=254) :: f_message
    integer :: len_olp_name, len_olp_version, len_message, i
    call olp_info(f_olp_name, f_olp_version, f_message)
    len_olp_name = len(trim(f_olp_name))
    len_olp_version = len(trim(f_olp_version))
    len_message = len(trim(f_message))
    do i = 1, len_olp_name
      olp_name(i) = f_olp_name(i:i)
    end do
    olp_name(len_olp_name+1) = c_null_char
    do i = 1, len_olp_version
      olp_version(i) = f_olp_version(i:i)
    end do
    olp_version(len_olp_version+1) = c_null_char
    do i = 1, len_message
      message(i) = f_message(i:i)
    end do
    message(len_message+1) = c_null_char
  end subroutine olp_info_c


  subroutine olp_printparameter_c(filename) bind(c,name="OLP_PrintParameter")
    ! BLHA OLP_PrintParameter routine.
    ! C wrapper to olp_printparameter
    ! [in] filename as C string
    use ol_iso_c_utilities, only: c_f_string
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: filename
    character(len=max_parameter_length) :: f_filename
    call c_f_string(filename, f_filename, max_parameter_length)
    call olp_printparameter(trim(f_filename))
  end subroutine olp_printparameter_c


  subroutine olp_start_c(contract_file_name, ierr) bind(C,name="OLP_Start")
    ! BLHA OLP_Start routine.
    ! C wrapper to olp_start
    ! [in] filename as C string
    use ol_iso_c_utilities, only: c_f_string
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: contract_file_name
    integer(c_int), intent(out) :: ierr
    character(len=max_parameter_length) :: f_contract_file_name
    integer :: f_ierr
    call c_f_string(contract_file_name, f_contract_file_name, max_parameter_length)
    call olp_start(trim(f_contract_file_name), f_ierr)
    ierr = f_ierr
  end subroutine olp_start_c


  subroutine olp_start_line_c(contract_line, answer_line, ierr) bind(C,name="OLP_StartLine")
    ! OLP_StartLine routine (not officially defined in BLHA standard)
    ! C wrapper to olp_start_line
    ! [in]  contract line as C string
    ! [out] answer line as C string
    ! [out] ierr: 1=ok, 0=error
    use ol_iso_c_utilities, only: c_f_string
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: contract_line
    character(kind=c_char), dimension(*), intent(out) :: answer_line
    integer(c_int), intent(out) :: ierr
    character(len=max_parameter_length) :: f_contract_line
    character(len=max_parameter_length) :: f_answer_line
    integer :: i, f_ierr
    call c_f_string(contract_line, f_contract_line, max_parameter_length)
    call olp_start_line(trim(f_contract_line), f_answer_line, f_ierr)
    do i = 1, len(trim(f_answer_line))
      answer_line(i) = f_answer_line(i:i)
    end do
    answer_line(i) = char(0)
    ierr = f_ierr
  end subroutine olp_start_line_c


  subroutine olp_start(contract_file_name, ierr)
    ! Fortran BLHA-like olp_start routine.
    ! Implements contract file handling: reads in BLHA contract file line-by-line
    ! and passes each line on to olp_start_line
    ! file is written.
    ! [in] filename
    ! [out] ierr: 1=ok, 0=error
    use iso_fortran_env, only: iostat_end
#ifdef USE_IFORT
    use ifport, only: system
#endif
    use ol_parameters_decl_/**/DREALKIND, only: tmp_dir
    implicit none
    integer, parameter :: cf = 993
    integer :: readok, commentpos

    character(len=*), intent(in) :: contract_file_name
    integer, intent(out) :: ierr
    integer ierrr
    character (len=max_parameter_length) :: linein
    character (len=max_parameter_length) :: lineout

    !set defaults
    flags%InterfaceVersion = 0
    flags%CorrectionType = 1
    flags%AmplitudeType = 11
    flags%answer_file_name = trim(contract_file_name)

    ierr = 1
    call ol_msg(1,"BLHA interface for OpenLoops invoked.")
    call ol_msg(1,"Reading contract file: " // contract_file_name)

    ! open contract files
    open(cf , file=contract_file_name, status = "old", iostat=readok )

    if (readok /= 0) then
      call ol_error(1, "can't find contract file " // trim(contract_file_name))
      ierr = 0
      return
    end if

    ! read contract file
    ReadLoop: do

      read (cf, "(A)", iostat=readok )  linein

      if (readok /= 0) then ! EOF -> exit
        if (readok == iostat_end) then
          exit ReadLoop
        else
          call ol_error(1, "error reading contract file")
          ierr = 0
          return
        end if
      end if

      lineout = linein ! remember original line
      commentpos = index(linein, "#")
      if (commentpos /= 0) then
        ! strip off comments
        linein = linein(1:commentpos-1)
      end if

      if (len_trim(linein) == 0) then
        ! ignore empty / comment only lines
        call add_blha_answer(trim(lineout))
        cycle ReadLoop
      end if

      call olp_start_line(linein, lineout, ierrr)
      call add_blha_answer(trim(lineout))
      if (ierrr == 1) then
        ierr = ierr*1
      else
        call ol_error(1,"error while reading contract file at: ")
        call ol_msg(1,"> " // trim(lineout))
        ierr = 0
      end if

    end do ReadLoop

    ! close contract files
    close (cf)

    ! Write answer file
    if (stdout_contract < 2) then
      call ol_msg(1,"Writing contract file to: " // trim(flags%answer_file_name))
      call write_blha_answer(ierrr)
      if (ierrr /= 1) then
        call ol_msg("Error: can't write answer file " // trim(flags%answer_file_name))
        ierr = 0
      end if
    end if

    if (ierr /= 1) then
    ! exit with error if something went wrong, i.e. not all required flags were set.
      call ol_error("Error reading/understanding BLHA contract file!")
    else
      ! everything looks fine -> exit without errors
      call ol_msg(1,"OLP_start: done. ")
    end if

  end subroutine olp_start



  subroutine olp_start_line(line, lineout, ierr)
  ! Handling of individual lines of BLHA contract file:
  ! - Set input flags & input parameters
  ! - Specify process and register process handle
  ! [in] line
  ! [out] lineout, ierr: 1=ok, 0=error
    use iso_fortran_env, only: iostat_end
    use ol_generic, only: to_lowercase
    use openloops, only: register_process
    character (len=*), intent(in) :: line
    character (len=*), intent(out) :: lineout
    integer, intent(out) :: ierr

    character(len=max_parameter_length) :: key
    character(len=max_parameter_length) :: val, param
    character(len=max_parameter_length) :: tmp
    character(len=max_parameter_length) :: inp
    character(len=max_parameter_length) :: outp
    character(len=max_parameter_length) :: librarypath
    character :: coupling_order(2)
    real(DREALKIND) :: paramInput
    integer :: ierrparam
    integer :: tmpi
    integer :: readok, ierrg
    integer :: libid

    ierr = 1

    ! everything but subprocesses should be either general initialisation and/or parameters
    if (index(line, '->')==0) then

      ! split: line -> (key, value)
      key = to_lowercase(adjustl(line(1:index(line, " "))))
      val = adjustl(line(index(line, " "):len(line)))

      ! Initialize
      select case (trim(key))
        ! required keys
        case ("interfaceversion")
          if (trim(to_lowercase(val)) == "blha1") then
            lineout = trim(line) // "      | OK"
            flags%InterfaceVersion = 1
          else if (trim(to_lowercase(val)) == "blha2") then
            lineout = trim(line) // "      | Ok"
            flags%InterfaceVersion = 2
          else
            lineout = trim(line) // "      | Error: unsupported flag. Currently only 'BLHA1' and 'BLHA2' are supported."
            ierr = 0
          end if
        case ("model")
          call set_parameter("model", trim(to_lowercase(val)), ierrparam)
          if (ierrparam == 0) then
            lineout = trim(line) // "      | OK"
          else
            lineout = trim(line) // "      | Error: unsupported model. Currently only 'SM', 'HEFT', supported"
            ierr = 0
          end if
        case ("correctiontype")
          if (trim(to_lowercase(val)) == "qcd") then
            lineout = trim(line) // "      | OK"
            flags%CorrectionType = 1
            call set_parameter("coupling_qcd_1", 1, ierrparam)
            call set_parameter("coupling_ew_1", 0, ierrparam)
          else if (trim(to_lowercase(val)) == "ew") then
            lineout = trim(line) // "      | OK"
            flags%CorrectionType = 2
            call set_parameter("coupling_qcd_1", 0, ierrparam)
            call set_parameter("coupling_ew_1", 1, ierrparam)
          else
            lineout = trim(line) // "      | Error: unsupported flag. Currently only 'QCD' and 'EW' supported."
            ierr = 0
          end if
        case ("irregularisation", "irregularization")
          if (trim(to_lowercase(val)) == "cdr") then
            lineout = trim(line) // "      | OK"
          else
            lineout = trim(line) // "      | Error: unsupported flag. Currently only 'CDR' supported"
            ierr = 0
          end if
        !
        ! optional keys
        case ("operationmode")
          if (flags%InterfaceVersion == 1) then ! only BLHA1
            if (trim(to_lowercase(val)) == "leadingcolour" ) then
              call set_parameter("leading_colour", 1, ierrparam)
              if ( ierrparam == 1) then
                lineout = trim(line) // "      | OK"
              else
                lineout = trim(line) // "      | Error"
              end if
            else
              lineout = trim(line) // "      | Error: unsupported flag. Currently only 'LeadingColour' supported"
              ierr = 0
            end if
          else
            lineout = trim(line) // "      | Error: only supported in BLHA1 interface version"
            ierr = 0
          end if
        case ("amplitudetype")            !only BLHA2
          select case (trim(to_lowercase(val)))
            case ("tree")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 1
            case ("cctree")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 11
            case ("sctree")
              lineout = trim(line) // "      | Error: not implemented; use sctree_polvect instead"
              flags%AmplitudeType = 3
              ierr = 0
            case ("sctree_polvect")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 11
            case ("loop")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 11
            case ("loopinduced")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 12
            case ("ccloopinduced")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 22
            case ("scloopinduced_polvect")
              lineout = trim(line) // "      | OK"
              flags%AmplitudeType = 24
            case default
              lineout = trim(line) // "      | Error: unsupported flag. Supported: Tree, ccTree, scTree, Loop, LoopInduced."
              ierr = 0
          end select
        case ("ewrenormalisationscheme")  !BLHA1 keyword
          if (trim(to_lowercase(val)) == "alpha0") then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 0)
          else if (trim(to_lowercase(val)) == "gmu") then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 1)
          else if (trim(to_lowercase(val)) == "alphamz") then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 2)
          else
            lineout = trim(line) // "      | Error: unsupported flag. Currently only 'alpha0' and 'Gmu' supported"
            ierr = 0
          end if
        case ("ewscheme")                 !BLAH2 keyword
          if (trim(to_lowercase(val)) == "alpha0" ) then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 0)
          else if (trim(to_lowercase(val)) == "gmu" ) then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 1)
          else if (trim(to_lowercase(val)) == "alphaMZ" ) then
            lineout = trim(line) // "      | OK"
            call set_parameter("ew_scheme", 2)
          else
            lineout = trim(line) // "      | Error: unsupported flag. Currently only 'alpha0' and 'Gmu' supported"
            ierr = 0
          end if
        case ("extra")
          if (index(to_lowercase(val), "openloopspath ") == 1) then
            val = adjustl(val(15:))
            if (len_trim(val) == 0) then
              lineout = trim(line) // "      | Error: path is missing"
              ierr = 0
            else
              if (len_trim(librarypath) == 0) then
                call set_parameter("install_path", trim(val), ierrparam)
                lineout = trim(line) // "      | OK"
              else
                lineout = trim(line) // "      | OK # Ignored: environment variable 'OpenLoopsPath' is preferred."
              end if
            end if
          else if (index(to_lowercase(val), "answerfile ") == 1) then
            val = adjustl(val(12:))
            if (len_trim(val) == 0) then
              lineout = trim(line) // "      | Error: filename is missing"
              ierr = 0
            else
              flags%answer_file_name = trim(val)
              lineout = trim(line) // "      | OK"
            end if
          else if (index(to_lowercase(val), "openloopspolecheck ") == 1) then
            val = adjustl(val(20:))
            if (len_trim(val) == 0) then
              lineout = trim(line) // "      | Error: switch is missing"
              ierr = 0
            else
              if ( trim(val) == "0") then
                flags%PoleCheck = 0
                call set_parameter("polecheck", 0, ierrparam)
                lineout = trim(line) // "      | OK"
              else if (trim(val) == "1") then
                flags%PoleCheck = 1
                call set_parameter("polecheck", 1, ierrparam)
                lineout = trim(line) // "      | OK"
              else
                lineout = trim(line) // "      | Error: switch not supported."
              end if
            end if
          else if (index(to_lowercase(val), "openloopsallowedlibs ") == 1) then
            val = adjustl(val(22:))
            if (len_trim(val) == 0) then
              lineout = trim(line) // "      | Error: space/comma seperated list is missing"
              ierr = 0
            else
              call set_parameter("allowed_libs", trim(val), ierrparam)
              lineout = trim(line) // "      | OK"
            end if
          else if (index(to_lowercase(val), "openloopsapproximation ") == 1) then
            val = adjustl(val(24:))
            if (len_trim(val) == 0) then
              lineout = trim(line) // "      | Error: specification of approximation missing"
              ierr = 0
            else
              call set_parameter("approximation", trim(val), ierrparam)
              lineout = trim(line) // "      | OK"
            end if
          else    ! try to pass to set_parameter directly
              param = to_lowercase(val(1:index(val, " ")-1))
              val   = to_lowercase(adjustl(val(index(val, " ")+1:)))
              call set_init_error_fatal(0)
              call set_parameter(trim(param), trim(val), ierrparam)
              call set_init_error_fatal(2)
              if (ierrparam == 0) then
                lineout = trim(line) // "      | OK"
              else
                lineout = trim(line) // "      | Error: parameter not supported"
                ierr = 0
              end if
          end if
        !
        !
        case ("couplingpower")       ! CouplingPower can change for different subprocesses. Always use latest.
          ! split: value -> (coupling, order)
          tmp = adjustl(val(index(val, " "):len(val)))
          val = adjustl(val(1:index(val, " ")))
          select case(trim(to_lowercase(val)))
            case("qcd")
              read(tmp(1:1),'(I1)') tmpi
              call set_parameter("coupling_qcd_0", tmpi, ierrparam)
              lineout = trim(line) // "     | OK"
            case("qed")
              read(tmp(1:1),'(I1)') tmpi
              call set_parameter("coupling_ew_0", tmpi, ierrparam)
              lineout = trim(line) // "     | OK"
            case default
              lineout = trim(line) // "     | Error: unsupported coupling type."
              ierr = 0
          end select

        case default !if not listed before let's see if its a parameter
          read(val,*) paramInput
          call set_parameter(key, paramInput, ierrparam)
          if ( ierrparam == 0 ) then
            lineout = trim(line) // "     | OK"
          else if ( ierrparam == 2) then
            lineout = trim(line) // "     | Error: unknown option. ignored."
          else
            lineout = trim(line) //  "     | Error: something went wrong setting this parameter."
            ierr=0
          end if
      end select

    ! and here should be the subprocesses
    else if (index(line, '->') > 0) then
      if (flags%AmplitudeType /= 0) then
        libid = register_process(line, flags%AmplitudeType)
          if (libid <= 0) then
            lineout = trim(line) // "     | Process not found"
            ierr = 0
          else
            write(lineout,'(A,A,I4)') trim(line),  " | 1 " , libid
          end if
      else
        call ol_error(1,"amplitude type not specified!")
        lineout = trim(line) // "     |  Error: amplitude type not specified!"
        ierr=0
        return
      end if
    end if

    if (stdout_contract > 0) call ol_msg(trim(lineout))

  end subroutine olp_start_line


  subroutine add_blha_answer(line)
  ! Appends line to array of BLHA answer 'blha_answer'
  ! [out] line
    implicit none
    character (len=*), intent(in) :: line
    character (len=max_parameter_length), allocatable :: blha_answer_bak(:)

    if (.not. allocated(blha_answer)) then
      allocate(blha_answer(1))
    else
      allocate(blha_answer_bak(size(blha_answer)))
      blha_answer_bak = blha_answer
      deallocate(blha_answer)
      allocate(blha_answer(size(blha_answer_bak)+1))
      blha_answer(1:size(blha_answer_bak)) = blha_answer_bak
      deallocate(blha_answer_bak)
    end if
    blha_answer(size(blha_answer)) = trim(line)
  end subroutine add_blha_answer


  subroutine write_blha_answer(ierr)
  ! Write BLHA answer in blha_answer to disk
  ! [out] ierr: 1=ok, 0=error
    implicit none
    integer, intent(out) :: ierr
    integer :: i, readok
    integer, parameter :: cfo = 992

    open(cfo, file=trim(flags%answer_file_name) , status = "REPLACE", iostat=readok)

    if (readok /= 0) then
      call ol_msg("Error: can't open output file " // trim(flags%answer_file_name))
      ierr = 0
      return
    end if

    do i=1, size(blha_answer)
      write(cfo,'(A)') trim(blha_answer(i))
    end do

    ! close/deallocate output file/array
    close (cfo)
    deallocate(blha_answer)
    ierr = 1
  end subroutine write_blha_answer


end module openloops_blha
