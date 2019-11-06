
#include <iostream>
#include "openloops.h"

int main() {

  int ierr = 0;
  double sqrts = 1000.;
  double mu = 100.;
  double result[4];
  double acc;

  // OLP_Start: OpenLoops reads contract file, registers processes and writes answer file
  OLP_Start("OLE_order.lh", &ierr);

  int id = 1;

  if (ierr == 1) {

    //obtain random phase-space point
    int n_ex = ol_n_external(id);
    double pp[5*n_ex];
    ol_phase_space_point(id, sqrts, pp);

    // evaluate process
    OLP_EvalSubProcess2(&id, pp, &mu, result, &acc);

    std::cout.precision(15);
    std::cout << std::endl;
    std::cout << "Tree and loop matrix element of the process" << std::endl;
    std::cout << "d dbar -> Z u ubar" << std::endl;
    std::cout << "for the phase-space point" << std::endl;
    for (int k = 0; k < 5; k++) {
      std::cout << "P[" << k+1 << "] = " << pp[5*k] << "  " << pp[5*k+1]
                << "  " << pp[5*k+2] << "  " << pp[5*k+3] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Tree:      " << result[3] << std::endl;
    std::cout << "Loop e^0:  " << result[2] << std::endl;
    std::cout << "Loop e^-1: " << result[1] << std::endl;
    std::cout << "Loop e^-2: " << result[0] << std::endl;
    std::cout << "Accuracy:  " << acc << std::endl;
    std::cout << std::endl;

  }

  return 0;
}
