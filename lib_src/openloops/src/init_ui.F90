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


module ol_init
  use KIND_TYPES, only: DREALKIND, QREALKIND
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int
  use ol_iso_c_utilities, only: c_f_string
  use ol_debug, only: ol_fatal, ol_msg, ol_error, set_verbose, get_verbose, do_not_stop
  implicit none
  private
  public :: set_init_error_fatal
  public :: set_parameter, get_parameter, parameters_flush, tree_parameters_flush
  public :: register_cleanup, cleanup
  public :: set_if_modified, set_if_modified_muren

  logical, save :: setparameter_tree_was_called = .true.
  logical, save :: setparameter_loop_was_called = .true.
  logical, save :: setparameter_muren_was_called = .true.
  logical, save :: setparameter_alphaQCD_was_called = .true.
  logical, save :: forwarded_init = .false.
  integer, save :: error = 0
  integer, save :: init_error_fatal = 2

  type cleanup_routine
    procedure(), pointer, nopass :: clean => null()
  end type cleanup_routine
  type(cleanup_routine), allocatable, save :: cleanup_routines(:)
  integer, save :: n_cleanup_routines = 0

  interface set_if_modified
    module procedure set_if_modified_bool, set_if_modified_int, set_if_modified_double, &
                     set_if_modified_cmplx, set_if_modified_cmplxcmplx, &
                     set_if_modified_string, set_if_modified_quad, set_if_modified_quadcmplx, &
                     set_if_modified_quadcmplxcmplx
  end interface set_if_modified

  interface set_parameter
    module procedure setparameter_int, setparameter_string,    &
                     setparameter_single, setparameter_double, &
                     setparameter_dcomplex, setparameter_bool
  end interface set_parameter

  interface get_parameter
    module procedure getparameter_int, getparameter_double
  end interface get_parameter

  contains

  subroutine set_init_error_fatal(flag)
    ! flag = 0: ignore, 1: warn, 2 (default): stop
    implicit none
    integer, intent(in) :: flag
    init_error_fatal = flag
  end subroutine set_init_error_fatal

  subroutine set_init_error_fatal_c(flag) bind(c,name="ol_set_init_error_fatal")
    implicit none
    integer(c_int), value :: flag
    init_error_fatal = flag
  end subroutine set_init_error_fatal_c

  subroutine set_if_modified_bool(current, new)
    implicit none
    logical, intent(inout) :: current
    logical, intent(in) :: new
    if (current .neqv. new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_bool

  subroutine set_if_modified_int(current, new)
    implicit none
    integer, intent(inout) :: current
    integer, intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_int

  subroutine set_if_modified_double(current, new)
    implicit none
    real(DREALKIND), intent(inout) :: current
    real(DREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_double

  subroutine set_if_modified_muren(current, new)
    implicit none
    real(DREALKIND), intent(inout) :: current
    real(DREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_muren_was_called = .true.
    end if
  end subroutine set_if_modified_muren

  subroutine set_if_modified_alphaQCD(current, new)
    implicit none
    real(DREALKIND), intent(inout) :: current
    real(DREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_alphaQCD_was_called = .true.
    end if
  end subroutine set_if_modified_alphaQCD

  subroutine set_if_modified_quad(current, new)
    implicit none
    real(QREALKIND), intent(inout) :: current
    real(QREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_quad

  subroutine set_if_modified_quadcmplx(current, new)
    implicit none
    complex(QREALKIND), intent(inout) :: current
    real(QREALKIND),    intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_quadcmplx

  subroutine set_if_modified_quadcmplxcmplx(current, new)
    implicit none
    complex(QREALKIND), intent(inout) :: current
    complex(QREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_quadcmplxcmplx

  subroutine set_if_modified_cmplx(current, new)
    implicit none
    complex(DREALKIND), intent(inout) :: current
    real(DREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_cmplx

  subroutine set_if_modified_cmplxcmplx(current, new)
    implicit none
    complex(DREALKIND), intent(inout) :: current
    complex(DREALKIND), intent(in) :: new
    if (current /= new) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_cmplxcmplx

  subroutine set_if_modified_string(current, new)
    implicit none
    character(*), intent(inout) :: current
    character(*), intent(in) :: new
    if (trim(current) /= trim(new)) then
      current = new
      setparameter_tree_was_called = .true.
      setparameter_loop_was_called = .true.
    end if
  end subroutine set_if_modified_string

  subroutine setparameter_int(param, val, err)
    ! Set an OpenLoops integer parameter.
    ! Must be flushed by parameters_flush() to take effect
    ! [in]  param: parameter name
    ! [in]  val: integer value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase, to_string
!    use nll_ppvj, only: nll_order, nll_ll, nll_add
    implicit none
    character(*), intent(in) :: param
    integer, intent(in)  :: val
    integer, intent(out), optional :: err

    error = 0

    call ol_msg(5, "setparameter_int: " // trim(param) // " " // to_string(val))

    select case (to_lowercase(param))

      case ("redlib1")
        if (expert_mode) then
          call set_if_modified(a_switch, val)
          if (val == 1 .or. val == 7) then
            call set_if_modified(se_integral_switch, val)
          else
            call set_if_modified(se_integral_switch, 3)
          end if
          auto_preset = .false.
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("redlib2")
        if (expert_mode) then
          call set_if_modified(a_switch_rescue, val)
          auto_preset = .false.
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("redlib3", "redlib_qp")
        if (expert_mode) then
          call set_if_modified(redlib_qp, val)
          auto_preset = .false.
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("redlib", "redlibs")
        if (expert_mode) then
          ! set libraries equal
          call set_if_modified(a_switch, val)
          call set_if_modified(a_switch_rescue, val)
          call set_if_modified(redlib_qp, val)
          if (val == 1 .or. val == 7) then
            call set_if_modified(se_integral_switch, val)
          else
            call set_if_modified(se_integral_switch, 3)
          end if
          auto_preset = .false.
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("stability_mode")
        if (expert_mode) then
          call set_if_modified(stability_mode, val)
          auto_preset = .false.
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("hp_mode")
        call set_if_modified(hp_mode, val)
        ! hp disabled
        if (val == 0) then
          call set_if_modified(hp_switch, 0)
        ! hard region
        else if (val == 1) then
          call set_if_modified(hp_switch, 1)
          call set_if_modified(hp_gamma_trig, .false.)
          call set_if_modified(hp_ir_trig, .false.)
          call set_if_modified(hp_irtri_trig, .false.)
        ! IR regions
        else if (val == 2) then
          call set_if_modified(hp_switch, 1)
          call set_if_modified(hp_gamma_trig, .true.)
          call set_if_modified(hp_ir_trig, .true.)
          call set_if_modified(hp_irtri_trig, .false.)
        end if
      case ("hp_switch")
        if (expert_mode) then
          call set_if_modified(hp_switch, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("use_bubble_vertex")
        if (expert_mode) then
          call set_if_modified(use_bubble_vertex, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("bubble_vertex")
        if (expert_mode) then
          call set_if_modified(bubble_vertex, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("sync_qp_kinematics")
        if (expert_mode) then
          call set_if_modified(sync_qp_kinematics, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("hp_fake_trig")
        if (expert_mode) then
          call set_if_modified(hp_fake_trig, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("hel_mem_opt","hel_mem_opt_switch")
        if (expert_mode) then
          if (val == 0) then
            call set_if_modified(hel_mem_opt_switch, .false.)
          else
            call set_if_modified(hel_mem_opt_switch, .true.)
          end if
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("deviation_mode")
        if (expert_mode) then
          call set_if_modified(deviation_mode, val)
          if (val /= 1 .and. val /= 2) then
            call ol_error(1,"unrecognised deviation_mode:" // to_string(val))
          end if
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("scaling_mode")
        if (expert_mode) then
          call set_if_modified(scaling_mode, val)
          if (val /= 1 .and. val /= 3) then
            call ol_msg(1, "unrecognised scaling_mode:" // to_string(val))
          end if
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("write_psp", "write_points")
        write_psp = val
      case ("write_parameters", "parameters_write")
        if (val == 0) write_params_at_start = .false.
        if (val == 1) write_params_at_start = .true.
      case ("check_poles")
        if (val == 1) then
          do_pole_checks = .true.
        else
          do_pole_checks = .false.
        end if
      case ("ti_monitor", "timonitor")
        ti_monitor = val
      case ("nf", "n_quarks")
        call set_if_modified(nf, val)
      case ("nq_nondecoupled", "minnf_alphasrun", "nf_alphasrun")
        call set_if_modified(nq_nondecoupl, val)
      case ("nc", "ncolours", "ncolors")
        ! affects only renormalisation, r2, and ir-subtraction
        call set_if_modified(nc, val)
      case ("freeyuk_on")
        if (val == 0) then
          yuk_from_mass = .true.
        else if (val == 1) then
          yuk_from_mass = .false.
        else
          call ol_error("freeyuk_on not available:" // trim(to_string(val)))
        end if
      case ("coupling_qcd_0", "coupling_qcd_t", "coupling_qcd_tree")
        coupling_qcd(0) = val
      case ("coupling_qcd_1", "coupling_qcd_l", "coupling_qcd_loop")
        coupling_qcd(1) = val
      case ("coupling_ew_0", "coupling_ew_t", "coupling_ew_tree")
        coupling_ew(0) = val
      case ("coupling_ew_1", "coupling_ew_l", "coupling_ew_loop")
        coupling_ew(1) = val
        call set_if_modified(do_ew_renorm, 1)
      case ("add_associated_ew")
        add_associated_ew = val
      case ("order_ew", "oew")
        coupling_ew(0) = val
        coupling_ew(1) = 0
        coupling_qcd(0) = -1
        coupling_qcd(1) = -1
      case ("order_qcd", "oqcd")
        coupling_ew(0) = -1
        coupling_ew(1) = -1
        coupling_qcd(0) = val
        coupling_qcd(1) = 0
        call set_if_modified(do_ew_renorm, 1)
      case ("loop_order_ew", "loew")
        coupling_ew = -1
        coupling_qcd = -1
        loop_order_ew = val
      case ("loop_order_qcd", "loqcd")
        coupling_ew = -1
        coupling_qcd = -1
        loop_order_qcd = val
      case ("olmode")
        OLmode = val
      case ("ct_on")
        call set_if_modified(ct_is_on, val)
      case ("r2_on")
        call set_if_modified(r2_is_on, val)
      case ("iop_on")
        if (val==1) then
          call set_if_modified(ir_is_on, 2)
        else
          call set_if_modified(ir_is_on, 1)
        end if
      case ("ir_on")
        call set_if_modified(ir_is_on, val)
      case ("tp_on")
        call set_if_modified(tp_is_on, val)
      case ("ckmorder")
        call set_if_modified(CKMORDER, val)
        if (val > 0) then
          call set_if_modified(model, "sm_ckm")
          if (flavour_mapping_on /= 0) then
            flavour_mapping_on = 2
          end if
        end if
      case ("ioperator_mode")
        call set_if_modified(ioperator_mode, val)
      case ("polecheck","truepoles","truepoles_on")
        if (val == 1) then
          call set_if_modified(coli_cache_use, 0)
          call set_if_modified(polecheck_is, val)
        else if (val == 2) then
          call set_if_modified(coli_cache_use, 0)
          call set_if_modified(polecheck_is, 1)
          call set_if_modified(do_ew_renorm, 1)
          call set_if_modified(ir_is_on, 2)
        end if
      case ("fermion_loops")
        call set_if_modified(swf, val)
      case ("nonfermion_loops")
        call set_if_modified(swb, val)
      case ("polenorm")
        call set_if_modified(norm_swi, val)
      case ("leading_colour")
        call set_if_modified(leadingcolour, val)
      case ("stability_log")
        stability_log = val
      case ("last_switch")
        call set_if_modified(l_switch, val)
      case ("check_ward_tree")
        call set_if_modified(ward_tree, val)
      case ("check_ward_loop")
        call set_if_modified(ward_loop, val)
      case ("out_symmetry")
        call set_if_modified(out_symmetry_on, val)
      case ("flavour_mapping")
        if (val < 0 .or. val > 3) then
          call ol_error("flavour_mapping not supported: " // to_string(val))
        else
          flavour_mapping_on = val
        end if
      case ("max_point", "maxpoint")
        call set_if_modified(maxpoint, val)
      case ("max_rank", "maxrank")
        call set_if_modified(maxrank, val)
      case ("me_cache")
        use_me_cache = val
      case ("ew_renorm")
        call set_if_modified(do_ew_renorm, val)
        if (use_bubble_vertex .ne. 0 .and. val .eq. 0) then
          call set_if_modified(bubble_vertex, 1)
        else
          call set_if_modified(bubble_vertex, 0)
        end if
      case ("qcd_renorm")
        call set_if_modified(do_qcd_renorm, val)
      case ("se_integral_switch")
        call set_if_modified(se_integral_switch, val)
      case ("ew_scheme", "ewscheme")
        if (val < 0 .or. val > 4) then
          call ol_error(1,"unrecognised ew_scheme: " // to_string(val))
        else
          call set_if_modified(ew_scheme, val)
          call set_if_modified(ew_renorm_scheme, abs(val))
        end if
      case ("ew_renorm_scheme")
        if (val /= 0 .and. val /= 1 .and. val /= 2) then
          call ol_error(1,"unrecognised ew_renorm_scheme:" // to_string(val))
        else
          call set_if_modified(ew_renorm_scheme, val)
        end if
      case ("select_ew")
        if (val == 0) then
          call set_if_modified(qed_on, .true.)
          call set_if_modified(weak_on, .true.)
        else if (val == 1) then
          call set_if_modified(qed_on, .true.)
          call set_if_modified(weak_on, .false.)
        else if (val == 2) then
          call set_if_modified(qed_on, .false.)
          call set_if_modified(weak_on, .true.)
        else
          call ol_error("select_ew not available:" // trim(to_string(val)))
        end if
      case ("onshell_photons_lsz", "onshell_photon_lsz")
        if (val == 1) then
          call set_if_modified(onshell_photons_lsz, .true.)
        else if (val == 0) then
          call set_if_modified(onshell_photons_lsz, .false.)
        else
          call ol_error(1,"unrecognised " // trim(param) // "=" // to_string(val))
        end if
      case ("offshell_photons_dimreg", "offshell_photon_dimreg", "offshell_photons_lsz", "offshell_photon_lsz")
        if (val == 1) then
          call set_if_modified(offshell_photons_lsz, .true.)
        else if (val == 0) then
          call set_if_modified(offshell_photons_lsz, .false.)
        else
          call ol_error(1,"unrecognised " // trim(param) // "=" // to_string(val))
        end if
      case ("all_photons_dimreg", "all_photon_dimreg", "delta_alphamz_dimreg")
        if (val == 1) then
          call set_if_modified(delta_alphamz_dimreg, .true.)
        else if (val == 0) then
          call set_if_modified(delta_alphamz_dimreg, .false.)
        else
          call ol_error(1,"unrecognised " // trim(param) // "=" // to_string(val))
        end if
      case ("complex_mass_scheme", "use_cms", "cms")
        if (val == 0 .or. val == 1 .or. val == 2) then
          call set_if_modified(cms_on, val)
        else
          call ol_error(1,"unrecognised " // trim(param) // "=" // to_string(val))
        end if
      case ("wf_v_select")
        if (val == 1 .or. val == 2 .or. val == 3) then
          call set_if_modified(wf_v_select, val)
        else
          call ol_error(1,"unrecognised " // trim(param) // "=" // to_string(val))
        end if
      case ("cll_tenred")
        call set_if_modified(cll_tenred, val)
      case ("cll_channels")
        call set_if_modified(cll_channels, val)
      case ("cll_log")
        call set_if_modified(cll_log, val)
      case ("olo_verbose")
        if (olo_verbose /= val) reset_olo = .true.
        call set_if_modified(olo_verbose, val)
      case ("olo_unit")
        if (olo_outunit /= val) reset_olo = .true.
        call set_if_modified(olo_outunit, val)
      case ("cuttools_idig")
        if (oppidig /= val) cuttools_not_init = .true.
        call set_if_modified(oppidig, val)
      case ("cuttools_scaloop")
        if (oppscaloop /= val) cuttools_not_init = .true.
        call set_if_modified(oppscaloop, val)
      case ("dd_red_mode")
        if (dd_red_mode /= val) dd_not_init = .true.
        call set_if_modified(dd_red_mode, val)
      case ("use_coli_cache")
        call set_if_modified(coli_cache_use, val)
      case ("no_collier_stop")
        if (val == 1) call set_if_modified(no_collier_stop, .true.)
      case ("ol_params_verbose", "parameters_verbose")
        parameters_verbose = val
      case ("verbose")
        call set_verbose(val)
      case ("do_not_stop")
        if (val == 1) then
          do_not_stop = .true.
        else
          do_not_stop = .false.
        end if
      case ("no_splash", "nosplash")
        if (val == 0) then
          nosplash = .false.
        else
          nosplash = .true.
        end if
      case ("apicheck")
        if (val == 0) then
          apicheck = .false.
        else
          apicheck = .true.
        end if
      case ("splash")
        if (val == 0) then
          nosplash = .true.
        else
          nosplash = .false.
        end if
      case ("check_collection")
        if (val == 1) then
          check_collection = .true.
        else
          check_collection = .false.
        end if
      case ("auto_preset")
        if (val == 1) then
          auto_preset = .true.
        elseif (val == 0) then
          auto_preset = .false.
        else
          call ol_error("auto_preset not available:" // trim(to_string(val)))
        end if
      case ("expert_mode")
         if (val == 1) then
          expert_mode = .true.
          apicheck = .false.
        elseif (val == 0) then
          expert_mode = .false.
          apicheck = .true.
        else
          call ol_error("expert_mode not available:" // trim(to_string(val)))
        end if
      case ("preset")
          call ol_msg("preset ignored (deprecated for OpenLoops >= 2.0).")
      case ("hp_alloc_mode")
        call set_if_modified(hp_alloc_mode, val)
      case default
        error = 1
        if (.not. forwarded_init) then
          forwarded_init = .true.
          call setparameter_double(param, real(val, DREALKIND))
          forwarded_init = .false.
          if (error == 1 .and. init_error_fatal == 1) then
            call ol_error(1, "ol_setparameter_int ignored unknown parameter '" // trim(param) // "'")
          end if
        end if

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("unknown parameter '" // trim(param) // "' in ol_setparameter_int")
        return
      end if
    end if
  end subroutine setparameter_int

  subroutine setparameter_bool(param, val, err)
    ! Set an OpenLoops integer parameter.
    ! Must be flushed by parameters_flush() to take effect
    ! [in]  param: parameter name
    ! [in]  val: logical value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase, to_string
    implicit none
    character(*), intent(in) :: param
    logical, intent(in)  :: val
    integer, intent(out), optional :: err

    error = 0

    call ol_msg(5, "setparameter_bool: " // trim(param) // " " // to_string(val))

    select case (to_lowercase(param))

      case ("hp_gamma_trig")
        call set_if_modified(hp_gamma_trig, val)
      case ("hp_automerge")
        call set_if_modified(hp_automerge, val)
      case ("hp_irtri_trig")
        call set_if_modified(hp_irtri_trig, val)
      case ("hp_ir_trig")
        call set_if_modified(hp_ir_trig, val)
      case ("ir_hacks")
        call set_if_modified(ir_hacks, val)

      case default
        error = 1

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("unknown parameter '" // trim(param) // "' in ol_setparameter_bool")
        return
      end if
    end if
  end subroutine setparameter_bool

  subroutine getparameter_int(param, val, err)
    ! Get an OpenLoops integer parameter.
    ! [in]  param: parameter name
    ! [out] val: integer value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase
    implicit none
    character(*), intent(in) :: param
    integer, intent(out)  :: val
    integer, intent(out), optional :: err

    error = 0

    select case (to_lowercase(param))

      case ("redlib1")
        val = a_switch
      case ("redlib2")
        val = a_switch_rescue
      case ("redlib3", "redlib_qp")
        val = redlib_qp
      case ("stability_mode")
        val = stability_mode
      case ("deviation_mode")
        val = deviation_mode
      case ("scaling_mode")
        val = scaling_mode
      case ("nf", "n_quarks")
        val = nf
      case ("nq_nondecoupled", "minnf_alphasrun", "nf_alphasrun")
        val = nq_nondecoupl
      case ("nc", "ncolours", "ncolors")
        val = nc
      case ("coupling_qcd_0", "coupling_qcd_t")
        val = coupling_qcd(0)
      case ("coupling_qcd_1", "coupling_qcd_l")
        val = coupling_qcd(1)
      case ("coupling_ew_0", "coupling_ew_t")
        val = coupling_ew(0)
      case ("coupling_ew_1", "coupling_ew_l")
        val = coupling_ew(1)
      case ("order_ew")
        val = order_ew
      case ("order_qcd")
        val = order_qcd
      case ("ct_on")
        val = ct_is_on
      case ("r2_on")
        val = r2_is_on
      case ("ir_on")
        val = ir_is_on
      case ("tp_on")
        val = tp_is_on
      case ("ckmorder")
        val = CKMORDER
      case ("ioperator_mode")
        val = ioperator_mode
      case ("polecheck","truepoles","truepoles_on")
        val = polecheck_is
      case ("fermion_loops")
        val = swf
      case ("nonfermion_loops")
        val = swb
      case ("polenorm")
        val = norm_swi
      case ("leading_colour")
        val = leadingcolour
      case ("stability_log")
        val = stability_log
      case ("last_switch")
        val = l_switch
      case ("check_ward_tree")
        val = ward_tree
      case ("check_ward_loop")
        val = ward_loop
      case ("out_symmetry")
        val = out_symmetry_on
      case ("flavour_mapping")
        val = flavour_mapping_on
      case ("max_point")
        val = maxpoint
      case ("max_rank")
        val = maxrank
      case ("me_cache")
        val = use_me_cache
      case ("ew_renorm")
        val = do_ew_renorm
      case ("se_integral_switch")
        val = se_integral_switch
      case ("ew_scheme")
        val = ew_scheme
      case ("ew_renorm_scheme")
        val = ew_renorm_scheme
      case ("complex_mass_scheme", "use_cms", "cms")
        val = cms_on
      case ("cll_tenred")
        val = cll_tenred
      case ("cll_channels")
        val = cll_channels
      case ("cuttools_idig")
        val = oppidig
      case ("cuttools_scaloop")
        val = oppscaloop
      case ("dd_red_mode")
        val = dd_red_mode
      case ("use_coli_cache")
        val = coli_cache_use
      case ("ol_params_verbose")
        val = parameters_verbose
      case ("hp_mode")
        val = hp_mode
      case ("hp_switch")
        val = hp_switch
      case ("bubble_vertex")
        val = bubble_vertex
      case ("sync_qp_kinematics")
        val = sync_qp_kinematics
      case ("hp_fake_trig")
        val = hp_fake_trig
      case ("verbose")
        call get_verbose(val)
      case("do_not_stop")
        if (do_not_stop) then
          val = 1
        else
          val = 0
        end if
      case ("welcome_length")
        val = welcome_length

      case default
        error = 1
        if (init_error_fatal == 1) then
          call ol_error(1,"getparameter_int ignored unknown parameter '" // trim(param) // "'")
        end if

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("unknown parameter '" // trim(param) // "' in ol_getparameter_int")
        return
      end if
    end if
  end subroutine getparameter_int



  subroutine setparameter_double(param, val, err)
    ! Set an OpenLoops double precision parameter.
    ! Must be flushed by parameters_flush() to take effect.
    ! Calls are passed to ol_setparameter_int() if param does not match and val==int(val)
    ! [in]  param: parameter name
    ! [in]  val: double precision value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase, to_string
    implicit none
    character(*), intent(in) :: param
    real(DREALKIND), intent(in) :: val
    integer, intent(out), optional :: err

    error = 0

    call ol_msg(5, "setparameter_double: " // trim(param) // " " // to_string(val))

    select case (to_lowercase(param))

      case ("mu", "renscale")
        if (mureg /= val) reset_mureg = .true.
        call set_if_modified(mureg_unscaled, val)
        call set_if_modified(muren_unscaled, val)
      case ("muren")
        call set_if_modified_muren(muren_unscaled, val)
      case ("mureg")
        if (mureg /= val) reset_mureg = .true.
        call set_if_modified(mureg_unscaled, val)
      case ("alphas", "alpha_s", "alpha_qcd")
        call set_if_modified_alphaQCD(alpha_QCD, val)
      case ("alpha", "alpha_qed")
        if (ew_scheme == 0 .or. ew_scheme == 2 .or. ew_scheme == -1 .or. ew_scheme == 4) then
          call set_if_modified(alpha_QED_input, val)
        else
           call ol_msg("WARNING: " // trim(param) // " ignored in ew_scheme=" // trim(to_string(ew_scheme)) // ".")
        end if
      case ("alpha_qed_mz")
        call set_if_modified(alpha_QED_MZ, val)
      case ("alpha_qed_0")
        call set_if_modified(alpha_QED_0, val)
      case ("gmu")
        call set_if_modified(Gmu_unscaled, val)
      case ("sw2eff", "sw2", "sintheta2eff")
        call set_if_modified(sw2_input, val)
      case ("scalefactor")
        if (scalefactor == 0) then
          call ol_error("scalefactor == 0 not supported!")
          return
        else
          if (scalefactor /= val) reset_scalefactor = .true.
          call set_if_modified(scalefactor, val)
        end if
      case ("rescalefactor")
        call set_if_modified(rescalefactor, val)
      case ("mass(1)", "d_mass", "rmd")
        call set_if_modified(rMD_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYD_unscaled, val)
      case ("width(1)", "d_width", "wmd")
        call set_if_modified(wMD_unscaled, val)
      case ("yuk(1)", "d_yuk", "yukmass(1)")
        call set_if_modified(rYD_unscaled, val)
      case ("mass(2)", "u_mass", "rmu")
        call set_if_modified(rMU_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYU_unscaled, val)
      case ("width(2)", "u_width", "wmu")
        call set_if_modified(wMU_unscaled, val)
      case ("yuk(2)", "u_yuk", "yukmass(2)")
        call set_if_modified(rYU_unscaled, val)
      case ("mass(3)", "s_mass", "rms")
        call set_if_modified(rMS_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYS_unscaled, val)
      case ("width(3)", "s_width", "wms")
        call set_if_modified(wMS_unscaled, val)
      case ("yuk(3)", "s_yuk", "yukmass(3)")
        call set_if_modified(rYS_unscaled, val)
      case ("mass(4)", "c_mass", "rmc")
        call set_if_modified(rMC_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYC_unscaled, val)
      case ("width(4)", "c_width", "wmc")
        call set_if_modified(wMC_unscaled, val)
      case ("lambdam(4)", "c_lambdam")
        call set_if_modified(LambdaMC2_unscaled, val**2)
        if (yuk_from_mass) call set_if_modified(LambdaYC2_unscaled, val**2)
      case ("yuk(4)", "c_yuk", "yukmass(4)")
        call set_if_modified(rYC_unscaled, val)
      case ("yukw(4)", "c_yukw", "yukwidth(4)")
        call set_if_modified(wYC_unscaled, val)
      case ("lambday(4)", "c_lambday")
        call set_if_modified(LambdaYC2_unscaled, val**2)
      case ("mass(5)", "b_mass", "rmb", "mb")
        call set_if_modified(rMB_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYB_unscaled, val)
      case ("width(5)", "b_width", "wmb")
        call set_if_modified(wMB_unscaled, val)
        if (yuk_from_mass) call set_if_modified(wYB_unscaled, val)
      case ("lambdam(5)", "b_lambdam")
        call set_if_modified(LambdaMB2_unscaled, val**2)
        if (yuk_from_mass) call set_if_modified(LambdaYB2_unscaled, val**2)
      case ("yuk(5)", "b_yuk", "yb", "yukmass(5)")
        call set_if_modified(rYB_unscaled, val)
      case ("yukw(5)", "b_yukw", "yukwidth(5)")
        call set_if_modified(wYB_unscaled, val)
      case ("lambday(5)", "b_lambday")
        call set_if_modified(LambdaYB2_unscaled, val**2)
      case ("mass(6)", "t_mass", "rmt", "mt")
        call set_if_modified(rMT_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYT_unscaled, val)
      case ("width(6)", "t_width", "wmt")
        call set_if_modified(wMT_unscaled, val)
        if (yuk_from_mass) call set_if_modified(wYT_unscaled, val)
      case ("lambdam(6)", "t_lambdam")
        call set_if_modified(LambdaMT2_unscaled, val**2)
        if (yuk_from_mass) call set_if_modified(LambdaYT2_unscaled, val**2)
      case ("yuk(6)", "t_yuk", "yt", "yukmass(6)")
        call set_if_modified(rYT_unscaled, val)
      case ("yukw(6)", "t_yukw", "yukwidth(6)")
        call set_if_modified(wYT_unscaled, val)
      case ("lambday(6)", "t_lambday")
        call set_if_modified(LambdaYT2_unscaled, val**2)
      case ("mass(11)", "e_mass", "rme")
        call set_if_modified(rME_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYE_unscaled, val)
      case ("width(11)", "e_width", "wme")
        call set_if_modified(wME_unscaled, val)
      case ("yuk(11)", "e_yuk")
        call set_if_modified(rYE_unscaled, val)
      case ("mass(13)", "mu_mass", "rmm")
        call set_if_modified(rMM_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYM_unscaled, val)
      case ("width(13)", "mu_width", "wmm")
        call set_if_modified(wMM_unscaled, val)
      case ("yuk(13)", "m_yuk", "mu_yuk")
        call set_if_modified(rYM_unscaled, val)
      case ("mass(15)", "tau_mass", "rml")
        call set_if_modified(rML_unscaled, val)
        if (yuk_from_mass) call set_if_modified(rYL_unscaled, val)
      case ("width(15)", "tau_width", "wml")
        call set_if_modified(wML_unscaled, val)
      case ("yuk(15)", "l_yuk", "tau_yuk")
        call set_if_modified(rYL_unscaled, val)
      case ("mass(23)", "z_mass", "rmz", "mz")
        call set_if_modified(rMZ_unscaled, val)
      case ("width(23)", "z_width", "wmz")
        call set_if_modified(wMZ_unscaled, val)
      case ("mass(24)", "w_mass", "rmw", "mw")
        call set_if_modified(rMW_unscaled, val)
      case ("width(24)", "w_width", "wmw")
        call set_if_modified(wMW_unscaled, val)
      case ("mass(25)", "h_mass", "rmh", "mh")
        call set_if_modified(rMH_unscaled, val)
      case ("width(25)", "h_width", "wmh")
        call set_if_modified(wMH_unscaled, val)
      case("x_width", "wx")
        if (trim(model) /= "sm_vaux") then
          call ol_msg("Warning: x_width can only be used with model sm_vaux")
        end if
        call set_if_modified(wMX_unscaled, val)
      case("y_width", "wy")
        if (trim(model) /= "sm_vaux") then
          call ol_msg("Warning: y_width can only be used with model sm_vaux")
        end if
        call set_if_modified(wMY_unscaled, val)
      case("mass(36)", "rma0", "ma0")
        call set_if_modified(rMA0_unscaled, val)
      case("width(36)", "wma0")
        call set_if_modified(wMA0_unscaled, val)
      case("mass(35)", "rmhh", "mhh")
        call set_if_modified(rMHH_unscaled, val)
      case("width(35)", "wmhh")
        call set_if_modified(wMHH_unscaled, val)
      case("mass(37)", "rmhp", "mhp")
        call set_if_modified(rMHp_unscaled, val)
      case("width(37)", "wmhp")
        call set_if_modified(wMHp_unscaled, val)
      case("tanb", "tan_b")
        call set_if_modified(thdmTB, val)
      case("sinba", "sin_ba")
        call set_if_modified(thdmSBA, val)
      case("lambda5")
        call set_if_modified(thdmL5, val)
      case("hqq_right")
        call set_if_modified(gH(1), val)
      case("hqq_left")
        call set_if_modified(gH(2), val)
      case("vckmdu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdu, val)
      case("vckmsu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMsu, val)
      case("vckmbu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbu, val)
      case("vckmdc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdc, val)
      case("vckmsc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMsc, val)
      case("vckmbc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbc, val)
      case("vckmdt")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdt, val)
      case("vckmst")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMst, val)
      case("vckmbt")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbt, val)
      case("vckmidu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdu, VCKMdu+CI*val)
      case("vckmisu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMsu, VCKMsu+CI*val)
      case("vckmibu")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbu, VCKMbu+CI*val)
      case("vckmidc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdc, VCKMdc+CI*val)
      case("vckmisc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMsc, VCKMsc+CI*val)
      case("vckmibc")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbc, VCKMbc+CI*val)
      case("vckmidt")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdt, VCKMdt+CI*val)
      case("vckmist")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMst, VCKMst+CI*val)
      case("vckmibt")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMbt, VCKMbt+CI*val)

      case("vckmdiag")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdu, val)
        call set_if_modified(VCKMsu, rZERO)
        call set_if_modified(VCKMbu, rZERO)
        call set_if_modified(VCKMdc, rZERO)
        call set_if_modified(VCKMsc, val)
        call set_if_modified(VCKMbc, rZERO)
        call set_if_modified(VCKMdt, rZERO)
        call set_if_modified(VCKMbt, val)
      case("vckmall")
        if (trim(model) /= "sm_ckm" .or. ckmorder == 0) then
          call ol_msg("Warning: non-diagonal CKM matrix can only be used with ckmorder /= 0")
        end if
        call set_if_modified(VCKMdu, val)
        call set_if_modified(VCKMsu, val)
        call set_if_modified(VCKMbu, val)
        call set_if_modified(VCKMdc, val)
        call set_if_modified(VCKMsc, val)
        call set_if_modified(VCKMbc, val)
        call set_if_modified(VCKMdt, val)
        call set_if_modified(VCKMbt, val)
      case("theta_cabi", "theta_cabibbo", "cabi")
        call set_if_modified(ThetaCabi, val)
      case("hpovev")
        call set_if_modified(HPOvev_unscaled, val)
      case("kapww", "hpokapww")
        call set_if_modified(HPOkapWW, val)
      case("kapzz", "hpokapzz")
        call set_if_modified(HPOkapZZ, val)
      case("epsww", "hpoepsww")
        call set_if_modified(HPOepsWW, val)
      case("aepsww", "hpoaepsww")
        call set_if_modified(HPOaepsWW, val)
      case("epszz", "hpoepszz")
        call set_if_modified(HPOepsZZ, val)
      case("aepszz", "hpoaepszz")
        call set_if_modified(HPOaepsZZ, val)
      case("epsza", "hpoepsza")
        call set_if_modified(HPOepsZA, val)
      case("aepsza", "hpoaepsza")
        call set_if_modified(HPOaepsZA, val)
      case("epsaa", "hpoepsaa")
        call set_if_modified(HPOepsAA, val)
      case("aepsaa", "hpoaepsaa")
        call set_if_modified(HPOaepsAA, val)
      case("epszner", "hpoepszner")
        call set_if_modified(HPOepsZnn(1,1), val)
      case("epsznmr", "hpoepsznmr")
        call set_if_modified(HPOepsZnn(2,1), val)
      case("epsznlr", "hpoepsznlr")
        call set_if_modified(HPOepsZnn(3,1), val)
      case("epszer", "hpoepszer")
        call set_if_modified(HPOepsZll(1,1), val)
      case("epszmr", "hpoepszmr")
        call set_if_modified(HPOepsZll(2,1), val)
      case("epszlr", "hpoepszlr")
        call set_if_modified(HPOepsZll(3,1), val)
      case("epsznel", "hpoepsznel")
        call set_if_modified(HPOepsZnn(1,2), val)
      case("epsznml", "hpoepsznml")
        call set_if_modified(HPOepsZnn(2,2), val)
      case("epsznll", "hpoepsznll")
        call set_if_modified(HPOepsZnn(3,2), val)
      case("epszel", "hpoepszel")
        call set_if_modified(HPOepsZll(1,2), val)
      case("epszml", "hpoepszml")
        call set_if_modified(HPOepsZll(2,2), val)
      case("epszll", "hpoepszll")
        call set_if_modified(HPOepsZll(3,2), val)
      case("epszdr", "hpoepszdr")
        call set_if_modified(HPOepsZdd(1,1), val)
      case("epszsr", "hpoepszsr")
        call set_if_modified(HPOepsZdd(2,1), val)
      case("epszbr", "hpoepszbr")
        call set_if_modified(HPOepsZdd(3,1), val)
      case("epszur", "hpoepszur")
        call set_if_modified(HPOepsZuu(1,1), val)
      case("epszcr", "hpoepszcr")
        call set_if_modified(HPOepsZuu(2,1), val)
      case("epsztr", "hpoepsztr")
        call set_if_modified(HPOepsZuu(3,1), val)
      case("epszdl", "hpoepszdl")
        call set_if_modified(HPOepsZdd(1,2), val)
      case("epszsl", "hpoepszsl")
        call set_if_modified(HPOepsZdd(2,2), val)
      case("epszbl", "hpoepszbl")
        call set_if_modified(HPOepsZdd(3,2), val)
      case("epszul", "hpoepszul")
        call set_if_modified(HPOepsZuu(1,2), val)
      case("epszcl", "hpoepszcl")
        call set_if_modified(HPOepsZuu(2,2), val)
      case("epsztl", "hpoepsztl")
        call set_if_modified(HPOepsZuu(3,2), val)
      case("epswlne", "hpoepswlne")
        call set_if_modified(HPOepsWln(1), val)
      case("epswlnm", "hpoepswlnm")
        call set_if_modified(HPOepsWln(2), val)
      case("epswlnl", "hpoepswlnl")
        call set_if_modified(HPOepsWln(3), val)
      case("epswdu", "hpoepswdu")
        call set_if_modified(HPOepsWqq(1), val)
      case("epswsc", "hpoepswsc")
        call set_if_modified(HPOepsWqq(2), val)
      case("epswbt", "hpoepswbt")
        call set_if_modified(HPOepsWqq(3), val)
      case ("fact_uv")
        call set_if_modified(x_uv, 1/val)
      case ("fact_ir")
        call set_if_modified(x_ir, 1/val)
      case ("pole_uv", "pole_uv1")
        call set_if_modified(de1_uv, val)
      case ("pole_ir1")
        call set_if_modified(de1_ir, val)
      case ("pole_ir2")
        call set_if_modified(de2_i_ir, val)
      case ("polescale")
        call set_if_modified(polescale, val)

      case ("stability_triggerratio")
        trigeff_targ = val
      case ("stability_unstable")
        call set_if_modified(abscorr_unst, val)
      case ("stability_kill")
        call set_if_modified(ratcorr_bad, val)
      case ("stability_kill2")
        call set_if_modified(ratcorr_bad_L2, val)

      case ("cll_pvthr")
        call set_if_modified(cll_pvthr, val)
      case ("cll_accthr")
        call set_if_modified(cll_accthr, val)
      case ("cll_mode3thr")
        call set_if_modified(cll_mode3thr, val)
      case ("ti_os_thresh")
        call set_if_modified(ti_os_thresh, val)
      case ("cuttools_rootsvalue")
        if (opprootsvalue /= val) cuttools_not_init = .true.
        call set_if_modified(opprootsvalue_unscaled, val)
      case ("cuttools_limitvalue")
        if (opplimitvalue /= val) cuttools_not_init = .true.
        call set_if_modified(opplimitvalue, val)
      case ("opp_threshold")
        if (oppthrs /= val) reset_olo = .true.
        call set_if_modified(oppthrs, val)
      case ("dd_c_threshold")
        if (c_pv_threshold /= val) dd_not_init = .true.
        call set_if_modified(c_pv_threshold, val)
      case ("dd_d_threshold")
        if (d_pv_threshold /= val) dd_not_init = .true.
        call set_if_modified(d_pv_threshold, val)
      case("psp_tolerance")
        psp_tolerance = val

      case ("lambda_hhh", "kappa_hhh")
        call set_if_modified(lambdaHHH, val)
      case ("lambda_hhhh", "kappa_hhhh")
        call set_if_modified(lambdaHHHH, val)
      case ("lambda_hww", "kappa_hww")
        call set_if_modified(lambdaHWW, val)
      case ("lambda_hzz", "kappa_hzz")
        call set_if_modified(lambdaHZZ, val)

      ! Hybrid precision mode thresholds
      case ("hp_step_thres")
        if (expert_mode) then
          call set_if_modified(hp_step_thres, val)
        else
          call ol_msg(0, trim(param) // " ignored. Only available in expert mode.")
        end if
      case ("hp_loopacc")
        call set_if_modified(hp_loopacc, val)
        call set_if_modified(hp_err_thres, max(0.,16.-hp_loopacc))
      case default

        error = 1
        if (.not. forwarded_init) then
          if (val == int(val)) then
            forwarded_init = .true.
            call setparameter_int(param, int(val))
            forwarded_init = .false.
          end if
          if (error == 1 .and. init_error_fatal == 1) then
            call ol_error(1, "ol_setparameter_double ignored unknown parameter '" // trim(param) // "'")
          end if
        end if

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("unknown parameter '" // trim(param) // "' in ol_setparameter_double")
        return
      end if
    end if
  end subroutine setparameter_double



  subroutine setparameter_single(param, val, err)
    ! single precision wrapper for setparameter_double()
    implicit none
    character(*), intent(in) :: param
    real, intent(in) :: val
    integer, intent(out), optional :: err
    call setparameter_double(param, real(val, DREALKIND), err)
  end subroutine setparameter_single



  subroutine setparameter_dcomplex(param, val, err)
    ! double precision complex wrapper for setparameter_double()
    implicit none
    character(*), intent(in) :: param
    complex(DREALKIND), intent(in) :: val
    integer, intent(out), optional :: err
    call setparameter_double(param, real(val, DREALKIND), err)
    if (aimag(val) /= 0) then
      if (present(err)) then
        call ol_error(1, "non-vanishing imaginary part in real parameter")
        err = 1
      else
        call ol_fatal("non-vanishing imaginary part in real parameter")
        return
      end if
    end if
  end subroutine setparameter_dcomplex



  subroutine getparameter_double(param, val, err)
    ! Get an OpenLoops double precision parameter.
    ! [in]  param: parameter name
    ! [out] val: double precision value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase
    implicit none
    character(*), intent(in) :: param
    real(DREALKIND), intent(out) :: val
    integer, intent(out), optional :: err

    error = 0

    select case (to_lowercase(param))

      case ("mu", "renscale")
        val = mureg
      case ("alphas", "alpha_s", "alpha_qcd")
        val = alpha_QCD
      case ("alpha", "alpha_qed", "alpha_qed_mz")
        val = alpha_QED
      case ("alpha_qed_0")
        val = alpha_qed_0
      case ("gmu")
        val = Gmu
      case ("scalefactor")
        val = scalefactor
      case ("rescalefactor")
        val = rescalefactor
      case ("mass(1)", "d_mass", "rmd")
        val = rMD
      case ("width(1)", "d_width", "wmd")
        val = wMD
      case ("yuk(1)", "d_yuk")
        val = rYD
      case ("mass(2)", "u_mass", "rmu")
        val = rMU
      case ("width(2)", "u_width", "wmu")
        val = wMU
      case ("yuk(2)", "u_yuk")
        val = rYU
      case ("mass(3)", "s_mass", "rms")
        val = rMS
      case ("width(3)", "s_width", "wms")
        val = wMS
      case ("yuk(3)", "s_yuk")
        val = rYS
      case ("mass(4)", "c_mass", "rmc")
        val = rMC
      case ("width(4)", "c_width", "wmc")
        val = wMC
      case ("yuk(4)", "c_yuk")
        val = rYC
      case ("mass(5)", "b_mass", "rmb", "mb")
        val = rMB
      case ("width(5)", "b_width", "wmb")
        val = wMB
      case ("yuk(5)", "b_yuk", "yb")
        val = rYB
      case ("mass(6)", "t_mass", "rmt", "mt")
        val = rMT
      case ("width(6)", "t_width", "wmt")
        val = wMT
      case ("yuk(6)", "t_yuk", "yt")
        val = rYT
      case ("mass(11)", "e_mass", "rme")
        val = rME
      case ("width(11)", "e_width", "wme")
        val = wME
      case ("yuk(11)", "e_yuk")
        val = rYE
      case ("mass(13)", "mu_mass", "rmm")
        val = rMM
      case ("width(13)", "mu_width", "wmm")
        val = wMM
      case ("yuk(13)", "m_yuk")
        val = rYM
      case ("mass(15)", "tau_mass", "rml")
        val = rML
      case ("width(15)", "tau_width", "wml")
        val = wML
      case ("yuk(15)", "l_yuk")
        val = rYL
      case ("mass(23)", "z_mass", "rmz", "mz")
        val = rMZ
      case ("width(23)", "z_width", "wmz")
        val = wMZ
      case ("mass(24)", "w_mass", "rmw", "mw")
        val = rMW
      case ("width(24)", "w_width", "wmw")
        val = wMW
      case ("mass(25)", "h_mass", "rmh", "mh")
        val = rMH
      case ("width(25)", "h_width", "wmh")
        val = wMH
      case ("fact_uv")
        val = 1/x_uv
      case ("fact_ir")
        val = 1/x_ir
      case ("pole_uv","pole_uv1")
        val = de1_uv
      case ("pole_ir1")
        val = de1_ir
      case ("pole_ir2")
        val = de2_i_ir
      case ("polescale")
        val = polescale

      case ("stability_triggerratio")
        val = trigeff_targ
      case ("stability_unstable")
        val = abscorr_unst
      case ("stability_kill")
        val = ratcorr_bad
      case ("stability_kill2")
        val = ratcorr_bad_L2

      case ("cll_pvthr")
        val = cll_pvthr
      case ("cll_accthr")
        val = cll_accthr
      case ("cll_mode3thr")
        val = cll_mode3thr
      case ("ti_os_thresh")
        val = ti_os_thresh
      case ("cuttools_rootsvalue")
        val = opprootsvalue_unscaled
      case ("cuttools_limitvalue")
        val = opplimitvalue
      case ("opp_threshold")
        val = oppthrs
      case ("dd_c_threshold")
        val = c_pv_threshold
      case ("dd_d_threshold")
        val = d_pv_threshold
      case ("psp_tolerance")
        val = psp_tolerance

      case ("lambda_hhh")
        val = lambdaHHH
      case ("lambda_hww")
        val = lambdaHWW
      case ("lambda_hzz")
        val = lambdaHZZ

      case default
        error = 1
        if (init_error_fatal == 1) then
          call ol_error(1,"getparameter_double ignored unknown parameter '" // trim(param) // "'")
        end if

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("error: unknown parameter '" // trim(param) // "' in ol_getparameter_double")
        return
      end if
    end if
  end subroutine getparameter_double



  subroutine setparameter_string(param, val, err)
    ! Set an OpenLoops string parameter.
    ! Calls are passed to set_parameter_double() if param does not match and val represents a number
    ! In this case, must be flushed by parameters_flush() to take effect.
    ! [in]  param: parameter name
    ! [in]  val: string value
    ! sets error flag: 0=ok, 1=ignored, 2=error(unused)
    use ol_parameters_decl_/**/DREALKIND
    use ol_loop_parameters_decl_/**/DREALKIND
    use ol_generic, only: to_lowercase, to_string
    implicit none
    character(*), intent(in) :: param
    character(*), intent(in) :: val
    integer, intent(out), optional :: err
    real(DREALKIND) :: real_parameter
    integer :: i

    error = 0

    call ol_msg(5, "setparameter_string: " // trim(param)  // " " // trim(val))

    if (len(val) > max_parameter_length - 80) then
      call ol_fatal("ol_setparameter_string: " // trim(param) // " value must not exceed " // &
             & trim(to_string(max_parameter_length)) // " characters. If necessary, increase " // &
             & "the limit by setting max_string_length in your openloops.cfg and recompile.")
      return
    end if

    select case (to_lowercase(param))

      case ("install_path")
        install_path = val
      case ("stability_logdir")
        if (stability_logdir /= val) then
          stability_logdir_not_created = .true.
          stability_logdir = val
        end if
      case ("tmp_dir")
        tmp_dir = val
      case ("allowed_libs", "allowed_libraries", "allowedlibs", "allowedlibraries")
        if (len(val) > max_parameter_length-2) then
          ! needs a leading and a trailing space
          call ol_fatal("ol_setparameter_string: " // trim(param) // " value must not exceed " // &
                 & trim(to_string(max_parameter_length-2)) // " characters")
          return
        end if
        allowed_libs = val
        do i = 1, max_parameter_length
          if (allowed_libs(i:i) == ",") allowed_libs(i:i) = " "
        end do
        ! leading and trailing space(s) are needed as boundaries
        allowed_libs = " " // adjustl(allowed_libs)
      case ("approximation", "approx")
        approximation = val
        if (trim(val) == "vbf") then
          partial_normal_order = .true.
        else
          partial_normal_order = .false.
        end if
      case ("shopping_list", "shopping_card", "ol_shopping")
        if (trim(val) /= "1") then
          shopping_list = val
        end if
        write_shopping_list = .true.
      case ("model")
        select case (to_lowercase(trim(val)))
          case ("sm", "smdiag", "sm_yuksel")
            call set_if_modified(model, "sm")
            call set_if_modified(nf, 6)
          case ("smckm", "smnondiag", "sm_ckm")
            call set_if_modified(model, "sm_ckm")
            call set_if_modified(nf, 6)
            if (flavour_mapping_on /= 0) then
              call set_if_modified(flavour_mapping_on, 2)
            end if
          case ("sm_vaux")
            call set_if_modified(model, "sm_vaux")
            call set_if_modified(nf, 6)
            call set_if_modified(cms_on, 0)
          case ("heft", "sm+ehc")
            call set_if_modified(model, "heft")
            call set_if_modified(nf, 5)
          case ("2hdm1", "thdm1", "2hdmi", "thdmi")
            call set_if_modified(model, "2hdm")
            call set_if_modified(thdm_type, 1)
            call set_if_modified(nf, 6)
          case ("2hdm", "thdm", "2hdm2", "thdm2", "2hdmii", "thdmii")
            call set_if_modified(model, "2hdm")
            call set_if_modified(thdm_type, 2)
            call set_if_modified(nf, 6)
          case ("hpoprodmfv_ufo", "hpoprodmfv_ufo_fixed", "higgspo")
            call set_if_modified(model, "higgspo")
            call set_if_modified(nf, 6)
        case default
          call ol_error(1, "unknown model: " // trim(val) // ", model set to: " // trim(model))
        end select


      case default

        error = 1
        ! if the string can be converted to a real number, foward it to ol_setparameter_double();
        ! note that integers will be forwarded to ol_setparameter_int() automatically
        read(val,*,iostat=error) real_parameter
        if (error == 0) then
          call setparameter_double(param, real_parameter)
        else
          error = 1
        end if
        if (error == 1 .and. init_error_fatal == 1) then
          call ol_error(1,"ol_setparameter_string ignored unknown parameter '" // trim(param) // "'")
        end if

    end select

    if (present(err)) then
      err = error
    else
      if (init_error_fatal == 2 .and. error /= 0) then
        call ol_fatal("unknown parameter '" // trim(param) // "' in ol_setparameter_string")
        return
      end if
    end if
  end subroutine setparameter_string



  subroutine parameters_flush() bind(c,name="ol_parameters_flush")
    use ol_parameters_init_/**/DREALKIND, only: parameters_init, loop_parameters_init
    use ol_parameters_init_/**/QREALKIND, only: parameters_init_qp=>parameters_init, &
                                                loop_parameters_init_qp=>loop_parameters_init
    use ol_parameters_init_/**/REALKIND, only: qcd_parameters_init
    implicit none
    if (setparameter_tree_was_called .or. setparameter_loop_was_called) then
      call parameters_init()
      call loop_parameters_init()
      setparameter_tree_was_called = .false.
      setparameter_loop_was_called = .false.
    end if
    if (setparameter_alphaQCD_was_called) then
      call qcd_parameters_init()
      setparameter_alphaQCD_was_called = .false.
    end if
    if (setparameter_muren_was_called) then
      call qcd_parameters_init(.true.)
      setparameter_muren_was_called = .false.
    end if
  end subroutine parameters_flush


  subroutine tree_parameters_flush() bind(c,name="ol_tree_parameters_flush")
    use ol_parameters_init_/**/DREALKIND, only: parameters_init
    use ol_parameters_init_/**/REALKIND, only: qcd_parameters_init
    implicit none
    if (setparameter_tree_was_called) then
      call parameters_init()
      setparameter_tree_was_called = .false.
    end if
    if (setparameter_alphaQCD_was_called) then
      call qcd_parameters_init()
      setparameter_alphaQCD_was_called = .false.
    end if
  end subroutine tree_parameters_flush


  subroutine register_cleanup(sub)
    implicit none
    procedure() :: sub
    type(cleanup_routine), allocatable, save :: cleanup_routines_tmp(:)
    if (.not. allocated(cleanup_routines)) then
      allocate(cleanup_routines(1))
    end if
    if (n_cleanup_routines == size(cleanup_routines)) then
      allocate(cleanup_routines_tmp(n_cleanup_routines))
      cleanup_routines_tmp = cleanup_routines
      deallocate(cleanup_routines)
      allocate(cleanup_routines(2*n_cleanup_routines))
      cleanup_routines(1:n_cleanup_routines) = cleanup_routines_tmp
      deallocate(cleanup_routines_tmp)
    end if
    n_cleanup_routines = n_cleanup_routines + 1
    cleanup_routines(n_cleanup_routines)%clean => sub
  end subroutine register_cleanup

  subroutine cleanup()
    implicit none
    integer :: k
    do k = 1, n_cleanup_routines
      call cleanup_routines(k)%clean()
    end do
    if (allocated(cleanup_routines)) deallocate(cleanup_routines)
    n_cleanup_routines = 0
  end subroutine cleanup


  ! ============ !
  ! C interfaces !
  ! ============ !

  subroutine setparameter_int_c(param, val) bind(c,name="ol_setparameter_int")
    ! Convert null terminated C string to Fortran string, then call set_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    integer(c_int), value :: val
    character(max_parameter_name_length) :: f_param
    integer :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    f_val = val
    call set_parameter(trim(f_param), f_val)
  end subroutine setparameter_int_c

  subroutine setparameter_bool_c(param, val) bind(c,name="ol_setparameter_bool")
    ! Convert null terminated C string to Fortran string, then call set_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    integer(c_int), value :: val
    character(max_parameter_name_length) :: f_param
    integer :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    f_val = val
    if (f_val .eq. 0) then
      call set_parameter(trim(f_param), .false.)
    else
      call set_parameter(trim(f_param), .true.)
    end if
  end subroutine setparameter_bool_c

  subroutine getparameter_int_c(param, val) bind(c,name="ol_getparameter_int")
    ! Convert null terminated C string to Fortran string, then call get_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    integer(c_int), intent(out) :: val
    character(max_parameter_name_length) :: f_param
    integer :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    call get_parameter(trim(f_param), f_val)
    val = f_val
  end subroutine getparameter_int_c

  subroutine setparameter_double_c(param, val) bind(c,name="ol_setparameter_double")
    ! Convert null terminated C string to Fortran string, then call set_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    real(c_double), value :: val
    character(max_parameter_name_length) :: f_param
    real(DREALKIND) :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    f_val = val
    call set_parameter(trim(f_param), f_val)
  end subroutine setparameter_double_c

  subroutine getparameter_double_c(param, val) bind(c,name="ol_getparameter_double")
    ! Convert null terminated C string to Fortran string, then call get_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    real(c_double), intent(out) :: val
    character(max_parameter_name_length) :: f_param
    real(DREALKIND) :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    call get_parameter(trim(f_param), f_val)
    val = f_val
  end subroutine getparameter_double_c

  subroutine setparameter_string_c(param, val) bind(c,name="ol_setparameter_string")
    ! Convert null terminated C strings to Fortran strings, then call set_parameter()
    use ol_parameters_decl_/**/DREALKIND, only: max_parameter_name_length, max_parameter_length
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: param
    character(kind=c_char), dimension(*), intent(in) :: val
    character(max_parameter_name_length) :: f_param
    character(max_parameter_length) :: f_val
    call c_f_string(param, f_param, max_parameter_name_length)
    call c_f_string(val, f_val, max_parameter_length)
    call set_parameter(trim(f_param), trim(f_val))
  end subroutine setparameter_string_c

end module ol_init
