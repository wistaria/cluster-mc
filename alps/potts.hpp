/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Swendsen-Wang Cluster Algorithm for Potts Model

#include <alps/parapack/worker.h>
#include <algorithm>
#include <vector>
#include <math/power.hpp>
#include <cluster/union_find.hpp>

using math::power2;
using math::power4;

class potts_worker : public alps::parapack::mc_worker {
private:
  typedef alps::parapack::mc_worker super_type;
  typedef alps::graph_helper<>::bond_descriptor bond_descriptor;
  typedef cluster::union_find::node fragment_t;

public:
  potts_worker(alps::Parameters const& params) :
    super_type(params), lattice(params), q(alps::evaluate("Q", params)),
    temperature(alps::evaluate("T", params)), mcs(params), spins(lattice.num_sites(), 0),
    fragments(lattice.num_sites()), flip(lattice.num_sites()) {
  }
  virtual ~potts_worker() {}

  static std::string program() {
    return "Swendsen-Wang Cluster Algorithm for Potts Model";
  }
  static std::string copyright() {
    return program() + "\n  Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>";
  }
      
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Number of Sites")
        << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::RealObservable("Number of Clusters")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Order Parameter Density^2")
        << alps::RealObservable("Order Parameter Density^4");
  }

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs;
    double prob = 1 - std::exp(-1 / temperature);

    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    BOOST_FOREACH(bond_descriptor b, lattice.bonds()) {
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
    BOOST_FOREACH(bond_descriptor b, lattice.bonds()) {
      ene -= (spins[lattice.source(b)] == spins[lattice.target(b)] ? 1.0 : 0.0);
    }

    obs["Number of Sites"] << (double)lattice.num_sites();
    obs["Temperature"] << temperature;
    obs["Inverse Temperature"] << 1 / temperature;
    obs["Number of Clusters"] << (double)nc;
    obs["Energy"] << ene;
    obs["Energy Density"] << ene / lattice.num_sites();
    obs["Energy^2"] << ene * ene;
    obs["Order Parameter Density^2"] << mag2;
    double fc = 2.0 / (q - 1);
    obs["Order Parameter Density^4"] << ((1 + fc) * power2(mag2) - fc * mag4);
  }

  void save(alps::ODump& dp) const { dp << mcs << spins; }
  void load(alps::IDump& dp) { dp >> mcs >> spins; }

private:
  alps::graph_helper<> lattice;
  int q; // number of states
  double temperature; // temperature
  alps::mc_steps mcs;
  std::vector<int> spins; // spin configuration
  std::vector<fragment_t> fragments;
  std::vector<int> flip;
};

class potts_evaluator : public alps::parapack::simple_evaluator {
public:
  potts_evaluator(alps::Parameters const&) {}
  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
    if (obs.has("Order Parameter Density^2") &&
        obs.has("Order Parameter Density^4")) {
      alps::RealObsevaluator m2 = obs["Order Parameter Density^2"];
      alps::RealObsevaluator m4 = obs["Order Parameter Density^4"];
      alps::RealObsevaluator binder("Binder Ratio of Order Parameter");
      binder = power2(m2) / m4;
      obs.addObservable(binder);
    }
  }
};
