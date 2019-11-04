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

// Site Percolation Problem

#include <alps/parapack/worker.h>
#include <algorithm>
#include <vector>
#include <standards/power.hpp>
#include <cluster/union_find.hpp>

using math::power2;

class percolation_site_worker : public alps::parapack::mc_worker {
private:
  typedef alps::parapack::mc_worker super_type;
  typedef alps::graph_helper<>::bond_descriptor bond_descriptor;
  typedef cluster::union_find::node fragment_t;

public:
  percolation_site_worker(alps::Parameters const& params) :
    super_type(params), lattice(params), probability(alps::evaluate("PROBABILITY", params)),
    mcs(params), occupied(lattice.num_sites()), fragments(lattice.num_sites()) {
  }
  virtual ~percolation_site_worker() {}

  static std::string program() {
    return "Site Percolation Problem";
  }
  static std::string copyright() {
    return program() + "\n  Copyright (C) 2014-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>";
  }
      
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Number of Sites")
        << alps::SimpleRealObservable("Occupation Probability")
        << alps::RealObservable("Number of Clusters")
        << alps::RealObservable("Strength of Largest Cluster")
        << alps::RealObservable("Cluster Size");
  }

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs;

    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int s = 0; s < lattice.num_sites(); ++s) occupied[s] = (uniform_01() < probability);
    BOOST_FOREACH(bond_descriptor b, lattice.bonds()) {
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

    obs["Number of Sites"] << (double)lattice.num_sites();
    obs["Occupation Probability"] << probability;
    obs["Number of Clusters"] << (double)nc;
    obs["Strength of Largest Cluster"] << wmax / lattice.num_sites();
    obs["Cluster Size"] << (mag2 - power2(wmax)) / lattice.num_sites();
  }

  void save(alps::ODump& dp) const { dp << mcs; }
  void load(alps::IDump& dp) { dp >> mcs; }

private:
  alps::graph_helper<> lattice;
  double probability; // occupation probability
  alps::mc_steps mcs;
  std::vector<bool> occupied;
  std::vector<fragment_t> fragments;
};
