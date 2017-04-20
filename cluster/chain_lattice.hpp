/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef CLUSTER_CHAIN_LATTICE_H
#define CLUSTER_CHAIN_LATTICE_H

namespace cluster {

class chain_lattice {
public:
  chain_lattice(unsigned int L) : length_(L) {}
  unsigned int num_sites() const { return length_; }
  unsigned int num_bonds() const { return num_sites(); }
  unsigned int source(unsigned int b) const { return b; }
  unsigned int target(unsigned int b) const { return (b == length_-1) ? 0 : b+1; }
  double phase(unsigned int s) const { return (s & 1) ? 1.0 : -1.0; }

  // dummy functions
  unsigned int num_plqs() const { return 0; }
  unsigned int plq2bond0(unsigned int p) const { return 0; }
  unsigned int plq2bond1(unsigned int p) const { return 0; }

private:
  unsigned int length_;
};

} // end namespace cluster

#endif // CLUSTER_CHAIN_LATTICE_HPP
