/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Swendsen-Wang Cluster Algorithm for Potts Model

#include <alps/parapack/parapack.h>
#include "potts.hpp"

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION(potts_worker::program());
PARAPACK_SET_COPYRIGHT(potts_worker::copyright());
PARAPACK_REGISTER_WORKER(potts_worker, "potts");
PARAPACK_REGISTER_EVALUATOR(potts_evaluator, "potts");
