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

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <algorithm>
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
#include "potts_options.hpp"

using math::power2;
using math::power4;

int main(int argc, char* argv[]) {
  std::cout << "Swendsen-Wang Cluster Algorithm for Square Lattice Potts Model\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  unsigned int q = p.q;
  double prob = 1 - std::exp(-1 / p.temperature);

  // square lattice
  lattice::square lattice(p.length);

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // spin configuration
  std::vector<int> spins(lattice.num_sites(), 0 /* all zero state */);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());
  std::vector<int> flip(lattice.num_sites());

  // observables
  stat::accumulator num_clusters("Number of Clusters"), energy("Energy Density"),
    magnetization2("Order Parameter^2"), magnetization4("Order Parameter^4");

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < p.therm + p.sweeps; ++mcs) {
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
    BOOST_FOREACH(fragment_t& f, fragments) {
      if (f.is_root()) {
        f.set_id(nc++);
        double w = f.weight();
        mag2 += power2(w);
        mag4 += power4(w);
      }
    }
    BOOST_FOREACH(fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));

    // flip spins
    for (int c = 0; c < nc; ++c) flip[c] = static_cast<int>(q * uniform_01());
    for (int s = 0; s < lattice.num_sites(); ++s)
      spins[s] = (spins[s] + flip[fragments[s].id()]) % q;

    double ene = 0;
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      ene -= (spins[lattice.source(b)] == spins[lattice.target(b)] ? 1.0 : 0.0);
    }

    if (mcs >= p.therm) {
      num_clusters << (double)nc;
      energy << ene / lattice.num_sites();
      magnetization2 << mag2;
      double fc = 2.0 / (q - 1);
      magnetization4 << ((1+fc) * power2(mag2) - fc * mag4);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (p.therm + p.sweeps) / elapsed << " MCS/sec\n";
  std::cout << num_clusters << std::endl
            << energy << std::endl
            << magnetization2 << std::endl
            << magnetization4 << std::endl
            << "Binder Ratio of Order Parameter = "
            << power2(magnetization2.mean()) / magnetization4.mean() << std::endl;
}
