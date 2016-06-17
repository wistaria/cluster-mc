/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 2014-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model

#define ALPS_INDEP_SOURCE

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cluster/observable.hpp>
#include <cluster/power.hpp>
#include <cluster/union_find.hpp>
#include "infinite_options.hpp"

using cluster::power2;
using cluster::power4;

int main(int argc, char* argv[]) {
  std::cout << "O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  unsigned int num_sites = p.num_sites;
  unsigned int sweeps = p.sweeps;
  unsigned int therm = p.therm;

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(num_sites / p.temperature));

  // spin configuration
  std::vector<int> spins(num_sites, 1);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(num_sites);
  std::vector<bool> to_flip(num_sites);

  // observables
  cluster::observable num_clusters, magnetization2, magnetization4;

  boost::timer tm;
  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (double t = r_time(); t < 1; t += r_time()) {
      int s0 = num_sites * uniform_01();
      int s1 = num_sites * uniform_01();
      if (spins[s0] == spins[s1]) unify(fragments, s0, s1);
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
    for (int f = 0; f < num_sites; ++f)
      fragments[f].set_id(cluster_id(fragments, f));

    // flip spins
    for (int c = 0; c < nc; ++c) to_flip[c] = (uniform_01() < 0.5);
    for (int s = 0; s < num_sites; ++s)
      if (to_flip[fragments[s].id()]) spins[s] ^= 1;

    if (mcs >= therm) {
      num_clusters << (double)nc;
      magnetization2 << mag2;
      magnetization4 << (3 * power2(mag2) - 2 * mag4);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time       = " << elapsed << " sec\n"
            << "Speed = " << (therm + sweeps) / elapsed << " MCS/sec\n";
  std::cout << "Number of Clusters = "
            << num_clusters.mean() << " +- " << num_clusters.error() << std::endl
            << "Magnetization^2 = "
            << magnetization2.mean() << " +- " << magnetization2.error() << std::endl
            << "Magnetization^4 = "
            << magnetization4.mean() << " +- " << magnetization4.error() << std::endl
            << "Binder Ratio of Magnetization = "
            << power2(magnetization2.mean()) / magnetization4.mean() << std::endl;
}
