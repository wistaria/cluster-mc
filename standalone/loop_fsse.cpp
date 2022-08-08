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
// [fixed-length SSE]

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <algorithm> // for std::swap
#include <iostream>
#include <random>
#include <vector>
#include <standards/accumulator.hpp>
#include <standards/power.hpp>
#include <standards/timer.hpp>
#include <cluster/union_find.hpp>
#include "loop_options.hpp"

using standards::power2;

enum operator_type { diagonal, offdiagonal, identity };

struct local_operator_t {
  local_operator_t() : type(identity) {}
  local_operator_t(int b) : type(diagonal), bond(b) {}
  void flip() { type = (type == diagonal ? offdiagonal : diagonal); }
  operator_type type;
  unsigned int bond;
  unsigned int upper_loop, lower_loop;
};

struct cluster_t {
  cluster_t(bool t = false) : to_flip(t), size(0), mag(0), length(0) {}
  bool to_flip;
  int size;
  int mag;
  int length;
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
  const double lb2 = nbonds * beta / 2;

  // random number generators
  std::mt19937 eng(p.seed);
  std::uniform_real_distribution<> r_uniform01;

  // vector of operators
  std::vector<local_operator_t> operators(nbonds);
  unsigned int nop = 0; // number of non-identity operators

  // spin configuration at t = 0 (1 for down and 0 for up)
  std::vector<int> spins(nsites, 0); // initialized with all up

  // cluster information
  std::vector<fragment_t> fragments;
  std::vector<unsigned int> current(nsites); // id of fragments at current time
  std::vector<cluster_t> clusters;

  // observables
  standards::accumulator energy("Energy Density"), smag("Staggered Magnetizetion^2"),
    ssus("Staggered Susceptibility"), usus("Uniform Susceptibility");

  //
  // Monte Carlo steps
  //

  standards::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {
    // adjust length of operator string
    if (nop > 0.8 * operators.size()) {
      std::vector<local_operator_t> operators_new(2 * operators.size());
      std::vector<local_operator_t>::iterator itr_new = operators_new.begin();
      for (auto itr = operators.begin(); itr!= operators.end(); ++itr, itr_new += 2) *itr_new = *itr;
      std::swap(operators, operators_new);
    }

    // initialize cluster information (setup s cluster fragments)
    fragments.resize(nsites);
    std::fill(fragments.begin(), fragments.end(), fragment_t());
    for (unsigned int s = 0; s < nsites; ++s) current[s] = s;

    for (auto oi = operators.begin(); oi != operators.end(); ++oi) {

      // diagonal update
      if (oi->type == identity) {
        unsigned int b = nbonds * r_uniform01(eng);
        if (spins[left(nbonds, b)] != spins[right(nbonds, b)] &&
            (operators.size() - nop) * r_uniform01(eng) < lb2) {
          *oi = local_operator_t(b);
          ++nop;
        } else {
          continue;
        }
      } else {
        if (oi->type == diagonal &&
            lb2 * r_uniform01(eng) < operators.size() - nop + 1) {
          oi->type = identity;
          --nop;
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
    }

    // connect bottom and top cluster fragments
    for (int s = 0; s < nsites; ++s) unify(fragments, s, current[s]);

    //
    // cluster flip
    //

    // assign cluster id & determine if clusters are to be flipped
    int nc = 0;
    for (auto& f : fragments) { if (f.is_root()) f.set_id(nc++); }
    clusters.resize(nc);
    for (auto& f : fragments) { f.set_id(cluster_id(fragments, f)); }
    for (int c = 0; c < nc; ++c) clusters[c] = cluster_t(r_uniform01(eng) < 0.5);

    // 'flip' operators & do improved measurements
    unsigned int t = 0;
    for (auto oi = operators.begin(); oi != operators.end(); ++oi) {
      if (oi->type == identity) continue;
      int id_l = fragments[oi->lower_loop].id();
      int id_u = fragments[oi->upper_loop].id();
      clusters[id_l].length += 2 * t;
      clusters[id_u].length -= 2 * t;
      if (clusters[id_l].to_flip ^ clusters[id_u].to_flip) oi->flip();
      ++t;
    }

    // flip spins & do improved measurements
    for (unsigned int s = 0; s < nsites; ++s) {
      int id = fragments[s].id();
      clusters[id].size += 1;
      clusters[id].mag += 1 - 2 * spins[s];
      clusters[id].length += nop;
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
    for (auto pi = clusters.begin(); pi != clusters.end(); ++pi) {
      s2 += power2(pi->size);
      m2 += power2(pi->mag);
      l2 += power2(pi->length);
    }

    energy << (0.25 * nbonds - nop / beta) / nsites;
    smag << 0.25 * s2 / nsites;
    usus << 0.25 * beta * m2 / nsites;
    ssus << 0.25 * beta * ((nop ? l2 / nop : 0) + s2) / (nop + 1) / nsites;
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (therm + sweeps) / elapsed << " MCS/sec\n";
  std::cout << energy << std::endl
            << smag<< std::endl
            << usus << std::endl
            << ssus << std::endl;
}
