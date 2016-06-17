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

#define ALPS_INDEP_SOURCE

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include <cluster/observable.hpp>
#include <cluster/power.hpp>
#include <cluster/union_find.hpp>
#include <cluster/square_lattice.hpp>
#include "percolation_options.hpp"

using cluster::power2;

int main(int argc, char* argv[]) {
  std::cout << "Bond Percolation Problem on Square Lattice\n";
  options p(argc, argv, 0.5);
  if (!p.valid) std::exit(127);
  double probability = p.probability;
  unsigned int sweeps = p.sweeps;

  // sqaure lattice
  cluster::square_lattice lattice(p.length);
  unsigned int num_sites = lattice.num_sites();
  unsigned int num_bonds = lattice.num_bonds();

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(num_sites);

  // observables
  cluster::observable num_clusters, strength, cluster_size;

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int b = 0; b < num_bonds; ++b)
      if (uniform_01() < probability) unify(fragments, lattice.source(b), lattice.target(b));

    // assign cluster id & accumulate cluster properties
    int nc = 0;
    double wmax = 0, mag2 = 0;
    for (int f = 0; f < fragments.size(); ++f)
      if (fragments[f].is_root()) {
        ++nc;
        double w = fragments[f].weight();
        wmax = std::max(wmax, w);
        mag2 += power2(w);
      }
    num_clusters << (double)nc;
    strength << wmax / num_sites;
    cluster_size << (mag2 - power2(wmax)) / num_sites;
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time       = " << elapsed << " sec\n"
            << "Speed = " << sweeps / elapsed << " MCS/sec\n";
  std::cout << "Number of Clusters = "
            << num_clusters.mean() << " +- " << num_clusters.error() << std::endl
            << "Strength of Largest Cluster = "
            << strength.mean() << " +- " << strength.error() << std::endl
            << "Average Cluster Size = "
            << cluster_size.mean() << " +- " << cluster_size.error() << std::endl;
}
