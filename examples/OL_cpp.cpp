
// C++ example program for the native interface of OpenLoops.
// It calculates the tree and loop matrix element of the process
// d dbar -> Z u ubar for a random phase-space point.

#include <iostream>
#include "openloops.h"

int main() {

  double sqrts = 1000., mu = 100., mZ = 91.2, alphas = 0.1;
  double m2_tree, m2_loop[3], acc;

  // coupling order alpha_ew^1, implies QCD correction for loop process
  ol_setparameter_int("order_ew", 1);

  // Set parameter: Z mass
  ol_setparameter_double("mass(23)", mZ);

  // Increase verbosity level to list loaded libraries
  ol_setparameter_int("verbose", 1);

  //
  // register one-loop amplitude for process d dbar -> Z u ubar.
  // The "ppzjj" process library must be installed before via
  // $ ./scons auto=ppzjj
  //
  // second argument of ol_register_process:
  // 1 for tree-like matrix elements (tree, color and spin correlations),
  // 11 for loop, 12 for loop^2
  int id = ol_register_process("1 -1 -> 23 2 -2", 11);

  // Initialize OpenLoops
  ol_start();


  if (id > 0) {
    // Set parameter: strong coupling
    ol_setparameter_double("alpha_s", alphas);
    // Set parameter: renormalization scale
    ol_setparameter_double("mu", mu);

    // Obtain a random phase-space point in the format pp[5*N] from Rambo
    double pp[5*ol_n_external(id)];
    ol_phase_space_point(id, sqrts, pp);
    std::cout.precision(15);
    std::cout << std::endl;
    std::cout << "Tree and loop matrix element of the process" << std::endl;
    std::cout << "d dbar -> Z u ubar" << std::endl;
    std::cout << "for the phase-space point" << std::endl;
    for (int k = 0; k < 5; k++) {
      std::cout << "P[" << k+1 << "] = " << pp[5*k] << "  " << pp[5*k+1]
                << "  " << pp[5*k+2] << "  " << pp[5*k+3] << std::endl;
    }

    // evaluate tree matrix element
    ol_evaluate_tree(id, pp, &m2_tree);

    // print tree result
    std::cout << std::endl;
    std::cout << "ol_evaluate_tree" << std::endl;
    std::cout << "Tree:       " << m2_tree << std::endl;

    // evaluate loop matrix element (which also returns the tree)
    ol_evaluate_loop(id, pp, &m2_tree, m2_loop, &acc);

    // print loop result
    std::cout << std::endl;
    std::cout << "ol_evaluate_loop" << std::endl;
    std::cout << "Tree:       " << m2_tree << std::endl;
    std::cout << "Loop ep^0:  " << m2_loop[0] << std::endl;
    std::cout << "Loop ep^-1: " << m2_loop[1] << std::endl;
    std::cout << "Loop ep^-2: " << m2_loop[2] << std::endl;
    std::cout << "Accuracy:   " << acc << std::endl;
    std::cout << std::endl;
  }


  ol_finish();

  return 0;
}
