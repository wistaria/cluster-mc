/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <alps/parapack/parapack.h>
#include "potts.hpp"
#include "infinite.hpp"
#include "percolation_bond.hpp"
#include "percolation_site.hpp"

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION("Cluster-MC: Cluster Algorithm Monte Carlo Methods");
PARAPACK_SET_COPYRIGHT("Cluster-MC: Cluster Algorithm Monte Carlo Methods\n  Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>");

PARAPACK_REGISTER_WORKER(potts_worker, "potts");
PARAPACK_REGISTER_EVALUATOR(potts_evaluator, "potts");

PARAPACK_REGISTER_WORKER(infinite_worker, "infinite");
PARAPACK_REGISTER_EVALUATOR(infinite_evaluator, "infinite");

PARAPACK_REGISTER_WORKER(percolation_bond_worker, "bond percolation");

PARAPACK_REGISTER_WORKER(percolation_site_worker, "site percolation");
