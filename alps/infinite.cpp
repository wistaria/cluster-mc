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

// O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model

#include <alps/parapack/parapack.h>
#include "infinite.hpp"

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION(infinite_worker::program());
PARAPACK_SET_COPYRIGHT(infinite_worker::copyright());
PARAPACK_REGISTER_WORKER(infinite_worker, "infinite");
PARAPACK_REGISTER_EVALUATOR(infinite_evaluator, "infinite");
