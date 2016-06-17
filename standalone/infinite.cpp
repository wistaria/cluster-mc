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

// O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model

#define ALPS_INDEP_SOURCE

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cluster/observable.hpp>
#include <cluster/power.hpp>
#include <cluster/union_find.hpp>
#include <cluster/fully_connected_lattice.hpp>
#include "infinite_options.hpp"

using cluster::power2;
using cluster::power4;

int main(int argc, char* argv[]) {
  std::cout << "O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);

  // fully-connected lattice
  cluster::fully_connected_lattice lattice(p.num_sites);

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(lattice.num_sites() / p.temperature));

  // spin configuration
  std::vector<int> spins(lattice.num_sites(), 1);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());
  std::vector<bool> flip(lattice.num_sites());

  // observables
  cluster::observable num_clusters("Number of Clusters"),
    magnetization_unimp("Magnetization (unimproved)"),
    magnetization2_unimp("Magnetization^2 (unimproved)"),
    magnetization4_unimp("Magnetization^4 (unimproved)"),
    magnetization2("Magnetization^2"), magnetization4("Magnetization^4");

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < p.therm + p.sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (double t = r_time(); t < 1; t += r_time()) {
      int s0 = lattice.num_sites() * uniform_01();
      int s1 = lattice.num_sites() * uniform_01();
      if (spins[s0] == spins[s1]) unify(fragments, s0, s1);
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
    for (int c = 0; c < nc; ++c) flip[c] = (uniform_01() < 0.5);
    for (int s = 0; s < lattice.num_sites(); ++s)
      if (flip[fragments[s].id()]) spins[s] ^= 1;

    double mu = 0;
    for (int s = 0; s < lattice.num_sites(); ++s) mu += 2 * spins[s] - 1;

    if (mcs >= p.therm) {
      num_clusters << (double)nc;
      magnetization_unimp << mu;
      magnetization2_unimp << power2(mu);
      magnetization4_unimp << power4(mu);
      magnetization2 << mag2;
      magnetization4 << (3 * power2(mag2) - 2 * mag4);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (p.therm + p.sweeps) / elapsed << " MCS/sec\n";
  std::cout << num_clusters << std::endl
            << magnetization_unimp << std::endl
            << magnetization2_unimp << std::endl
            << magnetization4_unimp << std::endl
            << magnetization2 << std::endl
            << magnetization4 << std::endl
            << "Binder Ratio of Magnetization = "
            << power2(magnetization2.mean()) / magnetization4.mean() << std::endl;
}
