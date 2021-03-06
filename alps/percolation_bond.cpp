/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Bond Percolation Problem

#include <alps/parapack/parapack.h>
#include "percolation_bond.hpp"

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION(percolation_bond_worker::program());
PARAPACK_SET_COPYRIGHT(percolation_bond_worker::copyright());
PARAPACK_REGISTER_WORKER(percolation_bond_worker, "bond percolation");
