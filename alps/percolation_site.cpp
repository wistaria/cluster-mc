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

// Site Percolation Problem

#include <alps/parapack/parapack.h>
#include "percolation_site.hpp"

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION(percolation_site_worker::program());
PARAPACK_SET_COPYRIGHT(percolation_site_worker::copyright());
PARAPACK_REGISTER_WORKER(percolation_site_worker, "site percolation");
