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

// Site Percolation Problem on Square Lattice

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <cmath>
#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <math/power.hpp>
#include <lattice/square.hpp>
#include <stat/accumulator.hpp>
#include <cluster/union_find.hpp>
#include "percolation_options.hpp"

using math::power2;

int main(int argc, char* argv[]) {
  std::cout << "Site Percolation Problem on Square Lattice\n";
  options p(argc, argv, 0.592746);
  if (!p.valid) std::exit(127);

  // square lattice
  lattice::square lattice(p.length);

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // configuration
  std::vector<bool> occupied(lattice.num_sites());

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());

  // observables
  stat::accumulator num_clusters("Number of Clusters"), strength("Strength of Largest Cluster"),
    cluster_size("Cluster Size");

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < p.sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int s = 0; s < lattice.num_sites(); ++s) occupied[s] = (uniform_01() < p.probability);
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      int s0 = lattice.source(b);
      int s1 = lattice.target(b);
      if (occupied[s0] && occupied[s1]) unify(fragments, s0, s1);
    }

    // accumulate cluster properties
    int nc = 0;
    double wmax = 0, mag2 = 0;
    BOOST_FOREACH(fragment_t& f, fragments) {
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
