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

#include <alps/parapack/worker.h>
#include <algorithm>
#include <vector>
#include <cluster/power.hpp>
#include <cluster/union_find.hpp>

using cluster::power2;
using cluster::power4;

class infinite_worker : public alps::parapack::mc_worker {
private:
  typedef alps::parapack::mc_worker super_type;
  typedef cluster::union_find::node fragment_t;

public:
  infinite_worker(alps::Parameters const& params) :
    super_type(params), temperature(alps::evaluate("T", params)),
    num_sites(alps::evaluate("N", params)),
    r_time(engine(), boost::exponential_distribution<>(num_sites / temperature)),
    mcs(params), spins(num_sites, 1), fragments(num_sites), to_flip(num_sites) {
  }
  virtual ~infinite_worker() {}

  static std::string program() {
    return "O(N) Swendsen-Wang Cluster Algorithm for Infinite Range Ising Model";
  }
  static std::string copyright() {
    return program() + "\n  Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>";
  }
      
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::RealObservable("Number of Clusters")
        << alps::RealObservable("Magnetization (unimproved)")
        << alps::RealObservable("Magnetization^2 (unimproved)")
        << alps::RealObservable("Magnetization^4 (unimproved)")
        << alps::RealObservable("Magnetization^2")
        << alps::RealObservable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs;

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

    double mu = 0;
    for (int s = 0; s < num_sites; ++s) mu += 2 * spins[s] - 1;

    obs["Number of Clusters"] << (double)nc;
    obs["Magnetization (unimproved)"] << mu;
    obs["Magnetization^2 (unimproved)"] << power2(mu);
    obs["Magnetization^4 (unimproved)"] << power4(mu);
    obs["Magnetization^2"] << mag2;
    obs["Magnetization^4"] << (3 * power2(mag2) - 2 * mag4);
  }

  void save(alps::ODump& dp) const { dp << mcs << spins; }
  void load(alps::IDump& dp) { dp >> mcs >> spins; }

private:
  double temperature; // temperature
  int num_sites; // number of lattice sites
  boost::variate_generator<engine_type&, boost::exponential_distribution<> > r_time;
  alps::mc_steps mcs;
  std::vector<int> spins; // spin configuration
  std::vector<fragment_t> fragments;
  std::vector<bool> to_flip;
};

class infinite_evaluator : public alps::parapack::simple_evaluator {
public:
  infinite_evaluator(alps::Parameters const&) {}
  void evaluate(alps::ObservableSet& obs) const {
    {
      alps::RealObsevaluator m2 = obs["Magnetization^2 (unimproved)"];
      alps::RealObsevaluator m4 = obs["Magnetization^4 (unimproved)"];
      alps::RealObsevaluator binder("Binder Ratio of Magnetization (unimproved)");
      binder = power2(m2) / m4;
      obs.addObservable(binder);
    }
    {
      alps::RealObsevaluator m2 = obs["Magnetization^2"];
      alps::RealObsevaluator m4 = obs["Magnetization^4"];
      alps::RealObsevaluator binder("Binder Ratio of Magnetization");
      binder = power2(m2) / m4;
      obs.addObservable(binder);
    }
  }
};
