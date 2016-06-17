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

// Swendsen-Wang Cluster Algorithm for Square-Lattice Potts Model

#define ALPS_INDEP_SOURCE

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <cluster/observable.hpp>
#include <cluster/power.hpp>
#include <cluster/union_find.hpp>
#include <cluster/square_lattice.hpp>
#include "potts_options.hpp"

using cluster::power2;
using cluster::power4;

int main(int argc, char* argv[]) {
  std::cout << "Swendsen-Wang Cluster Algorithm for Square Lattice Potts Model\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  unsigned int q = p.q;
  double prob = 1 - std::exp(-1 / p.temperature);
  unsigned int sweeps = p.sweeps;
  unsigned int therm = p.therm;

  // sqaure lattice
  cluster::square_lattice lattice(p.length);
  unsigned int num_sites = lattice.num_sites();

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // spin configuration
  std::vector<int> spins(num_sites, 0 /* all zero state */);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(num_sites);
  std::vector<int> flip(num_sites);

  // observables
  cluster::observable num_clusters, energy, magnetization2, magnetization4;

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      if (spins[lattice.source(b)] == spins[lattice.target(b)] && uniform_01() < prob)
        unify(fragments, lattice.source(b), lattice.target(b));
    }

    // assign cluster id & accumulate cluster properties
    int nc = 0;
    double mag2 = 0, mag4 = 0;
    for (int f = 0; f < num_sites; ++f) {
      if (fragments[f].is_root()) {
        fragments[f].set_id(nc++);
        double w = fragments[f].weight();
        mag2 += power2(w);
        mag4 += power4(w);
      }
    }
    for (int f = 0; f < num_sites; ++f) fragments[f].set_id(cluster_id(fragments, f));

    // determine whether clusters are flipped or not
    for (int c = 0; c < nc; ++c) flip[c] = static_cast<int>(q * uniform_01());

    // flip spins
    for (int s = 0; s < num_sites; ++s)
      spins[s] = (spins[s] + flip[fragments[s].id()]) % q;

    if (mcs >= therm) {
      num_clusters << (double)nc;
      double ne = 0;
      for (int b = 0; b < lattice.num_bonds(); ++b)
        ne += (spins[lattice.source(b)] == spins[lattice.target(b)] ? 1.0 : 0.0);
      energy << -ne / num_sites;
      magnetization2 << mag2;
      double fc = 2.0 / (q - 1);
      magnetization4 << ((1+fc) * power2(mag2) - fc * mag4);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time       = " << elapsed << " sec\n"
            << "Speed = " << (therm + sweeps) / elapsed << " MCS/sec\n";
  std::cout << "Number of Clusters = "
            << num_clusters.mean() << " +- " << num_clusters.error() << std::endl
            << "Energy per Site    = "
            << energy.mean() << " +- " << energy.error() << std::endl
            << "Order Parameter^2  = "
            << magnetization2.mean() << " +- " << magnetization2.error() << std::endl
            << "Order Parameter^4  = "
            << magnetization4.mean() << " +- " << magnetization4.error() << std::endl
            << "Binder Ratio of Order Parameter = "
            << power2(magnetization2.mean()) / magnetization4.mean() << std::endl;
}
