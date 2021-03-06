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

// Bond Percolation Problem on Square Lattice

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <standards/accumulator.hpp>
#include <standards/power.hpp>
#include <standards/timer.hpp>
#include <lattice/graph.hpp>
#include <cluster/union_find.hpp>
#include "percolation_options.hpp"

using standards::power2;

int main(int argc, char* argv[]) {
  std::cout << "Bond Percolation Problem on Square Lattice\n";
  options p(argc, argv, 0.5);
  if (!p.valid) std::exit(127);

  // square lattice
  auto lattice = lattice::graph::simple(2, p.length);

  // random number generators
  std::mt19937 eng(p.seed);
  std::uniform_real_distribution<> r_uniform01;

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());

  // observables
  standards::accumulator num_clusters("Number of Clusters"), strength("Strength of Largest Cluster"),
    cluster_size("Cluster Size");

  standards::timer tm;
  for (unsigned int mcs = 0; mcs < p.sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int b = 0; b < lattice.num_bonds(); ++b)
      if (r_uniform01(eng) < p.probability) unify(fragments, lattice.source(b), lattice.target(b));

    // accumulate cluster properties
    int nc = 0;
    double wmax = 0, mag2 = 0;
    for (auto& f : fragments) {
      if (f.is_root()) {
        ++nc;
        double w = f.weight();
        wmax = std::max(wmax, w);
        mag2 += power2(w);
      }
    }

    num_clusters << (double)nc;
    strength << wmax / lattice.num_sites();
    cluster_size << (mag2 - power2(wmax)) / lattice.num_sites();
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << p.sweeps / elapsed << " MCS/sec\n";
  std::cout << num_clusters << std::endl
            << strength << std::endl
            << cluster_size << std::endl;
}
