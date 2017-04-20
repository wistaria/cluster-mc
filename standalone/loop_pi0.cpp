/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Loop Algorithm for Spin-1/2 Antiferromagnetic Heisenberg Chain
// [continuous time path integral; using std::list<> for operator string]

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/utility.hpp> // for boost::prior
#include <iostream>
#include <list>
#include <vector>
#include <math/power.hpp>
#include <stat/accumulator.hpp>
#include <cluster/union_find.hpp>
#include "loop_options.hpp"

using math::power2;

enum operator_type { diagonal, offdiagonal };

struct local_operator_t {
  local_operator_t() {}
  local_operator_t(int b, double t) : type(diagonal), bond(b), time(t) {}
  void flip() { type = (type == diagonal ? offdiagonal : diagonal); }
  operator_type type;
  unsigned int bond;
  unsigned int upper_loop, lower_loop;
  double time;
};

struct cluster_t {
  cluster_t(bool t = false) : to_flip(t), size(0), mag(0), length(0) {}
  bool to_flip;
  int size;
  int mag;
  double length;
};

typedef cluster::union_find::node fragment_t;

// lattice helper functions (returns site index at left/right end of a bond)
inline int left(int /* L */, int b) { return b; }
inline int right(int L, int b) { return (b == L-1) ? 0 : b+1; }

int main(int argc, char* argv[]) {
  std::cout << "Loop Algorithm for Spin-1/2 Antiferromagnetic Heisenberg Chain\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  const unsigned int nsites = p.length;
  const unsigned int nbonds = nsites;
  const unsigned int sweeps = p.sweeps;
  const unsigned int therm = p.therm;
  const double beta = 1. / p.temperature;

  // random number generators
  boost::mt19937 eng(p.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    random(eng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(beta * nbonds / 2));

  // list of operators
  std::list<local_operator_t> operators;

  // spin configuration at t = 0 (1 for down and 0 for up)
  std::vector<int> spins(nsites, 0); // initialized with all up

  // cluster information
  std::vector<fragment_t> fragments;
  std::vector<unsigned int> current(nsites); // id of fragments at current time
  std::vector<cluster_t> clusters;

  // oservables
  stat::accumulator energy("Energy Density"), smag("Staggered Magnetizetion^2"),
    ssus("Staggered Susceptibility"), usus("Uniform Susceptibility");

  //
  // Monte Carlo steps
  //

  boost::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {
    // initialize cluster information (setup s cluster fragments)
    fragments.resize(nsites);
    std::fill(fragments.begin(), fragments.end(), fragment_t());
    for (unsigned int s = 0; s < nsites; ++s) current[s] = s;

    double t = r_time();
    for (std::list<local_operator_t>::iterator oi = operators.begin();
         t < 1 || oi != operators.end();) {

      // diagonal update
      if (oi == operators.end() || t < oi->time) {
        unsigned int b = nbonds * random();
        if (spins[left(nbonds, b)] != spins[right(nbonds, b)]) {
          oi = operators.insert(oi, local_operator_t(b, t));
          t += r_time();
        } else {
          t += r_time();
          continue;
        }
      } else {
        if (oi->type == diagonal) {
          operators.erase(boost::prior(++oi));
          continue;
        }
      }

      // cluster generation
      unsigned int s0 = left(nbonds, oi->bond);
      unsigned int s1 = right(nbonds, oi->bond);
      oi->lower_loop = unify(fragments, current[s0], current[s1]);
      oi->upper_loop = current[s0] = current[s1] = add(fragments);
      if (oi->type == offdiagonal) {
        spins[s0] ^= 1;
        spins[s1] ^= 1;
      }
      ++oi;
    }

    // connect bottom and top cluster fragments
    for (int s = 0; s < nsites; ++s) unify(fragments, s, current[s]);

    //
    // cluster flip
    //

    // assign cluster id & determine if clusters are to be flipped
    int nc = 0;
    BOOST_FOREACH(fragment_t& f, fragments) { if (f.is_root()) f.set_id(nc++); }
    clusters.resize(nc);
    BOOST_FOREACH(fragment_t& f, fragments) { f.set_id(cluster_id(fragments, f)); }
    for (int c = 0; c < nc; ++c) clusters[c] = cluster_t(random() < 0.5);

    // 'flip' operators & do improved measurements
    for (std::list<local_operator_t>::iterator oi = operators.begin();
         oi != operators.end(); ++oi) {
      int id_l = fragments[oi->lower_loop].id();
      int id_u = fragments[oi->upper_loop].id();
      clusters[id_l].length += 2 * oi->time;
      clusters[id_u].length -= 2 * oi->time;
      if (clusters[id_l].to_flip ^ clusters[id_u].to_flip) oi->flip();
    }

    // flip spins & do improved measurements
    for (unsigned int s = 0; s < nsites; ++s) {
      int id = fragments[s].id();
      clusters[id].size += 1;
      clusters[id].mag += 1 - 2 * spins[s];
      clusters[id].length += 1;
      if (clusters[id].to_flip) spins[s] ^= 1;
    }

    if (mcs < therm) continue;

    //
    // measurements
    //

    // accumurate loop length and magnetization
    double s2 = 0;
    double m2 = 0;
    double l2 = 0;
    for (std::vector<cluster_t>::const_iterator pi = clusters.begin();
         pi != clusters.end(); ++pi) {
      s2 += power2(pi->size);
      m2 += power2(pi->mag);
      l2 += power2(pi->length);
    }

    energy << (0.25 * nbonds - operators.size() / beta) / nsites;
    smag << 0.25 * s2 / nsites;
    usus << 0.25 * beta * m2 / nsites;
    ssus << 0.25 * beta * l2 / nsites;
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (therm + sweeps) / elapsed << " MCS/sec\n";
  std::cout << energy << std::endl
            << smag<< std::endl
            << usus << std::endl
            << ssus << std::endl;
}
