/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <lattice/graph.hpp>
#include "cluster/union_find.hpp"

int main(int argc, char* argv[]) {
  std::size_t num_threads = omp_get_max_threads();
  int seed = 12345;
  int length = 64;
  double p = 0.5;

  // square lattice
  auto lattice = lattice::graph::simple(2, length);

  // generate random bonds
  std::mt19937 eng(seed);
  std::uniform_real_distribution<> r_uniform01;
  std::vector<int> bonds(lattice.num_bonds());
  for (auto& bond : bonds) bond = (r_uniform01(eng) < p);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());

  //
  // serial execution
  //
  
  std::fill(fragments.begin(), fragments.end(), fragment_t());

  for (std::size_t b = 0; b < bonds.size(); ++b)
    if (bonds[b]) unify(fragments, lattice.source(b), lattice.target(b));
  
  // accumulate cluster properties
  std::size_t nc_s = 0;
  std::size_t wmax_s = 0;
  std::size_t w2_s = 0;
  for (auto& f : fragments) {
    if (f.is_root()) {
      ++nc_s;
      std::size_t w = f.weight();
      wmax_s = std::max(wmax_s, w);
      w2_s += w * w;
    }
  }

  //
  // openmp execution
  //
  
  std::fill(fragments.begin(), fragments.end(), fragment_t());

  #pragma omp parallel for
  for (std::size_t b = 0; b < bonds.size(); ++b)
    if (bonds[b]) unify(fragments, lattice.source(b), lattice.target(b));
  
  // accumulate cluster properties
  std::size_t nc_p = 0;
  std::size_t wmax_p = 0;
  std::size_t w2_p = 0;
  for (auto& f : fragments) {
    if (f.is_root()) {
      ++nc_p;
      std::size_t w = f.weight();
      wmax_p = std::max(wmax_p, w);
      w2_p += w * w;
    }
  }

  std::clog << "number of threads = 1, " << num_threads << std::endl
            << "number of clusters = " << nc_s << ", " << nc_p << std::endl
            << "largest cluster size = " << wmax_s << ", " << wmax_p << std::endl
            << "sum of square of cluster size = " << w2_s << ", " << w2_p << std::endl;

  if (nc_s != nc_p || wmax_s != wmax_p || w2_s != w2_p) {
    std::cerr << "result mismatch\n";
    return 127;
  }
  return 0;
}
