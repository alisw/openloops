
! Example  program how to use the native interface of OpenLoops.
! It calculates the Tree and loop matrix element of the process
! d dbar -> Z u ubar for a random phase-space point.

program main
  use openloops
  implicit none
  integer :: id, error, k
  real(8) :: m2_tree, m2_loop(0:2), acc
  real(8) :: p_ex(0:3,5)
  real(8) :: mZ = 91.2
  real(8) :: mu = 100, alpha_s = 0.1, energy=1000

  ! coupling order alpha_ew^1, implies QCD correction for loop process
  call set_parameter("order_ew", 1)

  ! set Z mass
  call set_parameter("mass(23)", mZ)

  ! Increase verbosity level to list loaded libraries
  call set_parameter("verbose", 1)

  ! register one-loop amplitude for process d dbar -> Z u ubar
  ! The "ppzjj" process library must be installed before via
  ! $ ./scons auto=ppzjj
  !
  ! second argument of register_process:
  ! 1 for tree-like matrix elements (tree, color and spin correlations),
  ! 11 for loop, 12 for loop^2
  id = register_process("1 -1 -> 23 2 -2", 11)

  ! start
  call start()

  if (id > 0) then
    ! generate a random phase space point with rambo
    call phase_space_point(id, energy, p_ex)

    ! set strong coupling
    call set_parameter("alpha_s", alpha_s)
    ! set renormalisation scale
    call set_parameter("mu", mu)

    print *
    print *, "Tree and loop matrix element of the process"
    print *, "d dbar -> Z u ubar"
    print *, "for the phase-space point"
    do k = 1, size(p_ex,2)
      print *, 'P[', int(k,1), '] =', p_ex(:,k)
    end do

    ! evaluate tree matrix element
    call evaluate_tree(id, p_ex, m2_tree)
    ! print tree result
    print *
    print *, "evaluate_tree"
    print *, "Tree:       ", m2_tree

    ! evaluate tree and loop matrix elements
    call evaluate_loop(id, p_ex, m2_tree, m2_loop(0:2), acc)
    ! print loop result
    print *
    print *, "evaluate_loop"
    print *, "Tree:       ", m2_tree
    print *, "Loop ep^0:  ", m2_loop(0)
    print *, "Loop ep^-1: ", m2_loop(1)
    print *, "Loop ep^-2: ", m2_loop(2)
    print *, "accuracy:   ", acc
    print *
  end if

  ! finish
  call finish()

end program main
