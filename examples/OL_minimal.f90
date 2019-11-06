
! Minimal example how to use the native interface of OpenLoops.
! It calculates the Tree and loop matrix element of the process
! u dbar -> W+ g for a random phase-space point.

program main
  use openloops
  implicit none
  integer :: id
  real(8) :: psp(0:3,5), m2l0, m2l1(0:2), acc, sqrt_s = 1000

  call set_parameter("order_ew", 1) ! coupling order
  id = register_process("1 -1 -> 23 2 -2", 11) ! register loop process d dbar --> Z u ubar
  call start() ! initialise
  call set_parameter("mass(23)", 91.2);
  call set_parameter("alpha_s", 0.1)
  call set_parameter("mu", 100)
  call phase_space_point(id, sqrt_s, psp) ! obtain a random phase-space point from Rambo
  call evaluate_loop(id, psp, m2l0, m2l1, acc) ! calculate tree and loop matrix elements
  call finish()

  print *, "  d dbar --> Z u ubar for a RAMBO phase space point"
  print *, "  tree                     loop ep^0                  loop ep^-1               loop ep^-2                 accuracy"
  print *, m2l0, m2l1, acc
end program main
