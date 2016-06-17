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

#ifndef CLUSTER_FULLY_CONNECTED_LATTICE_HPP
#define CLUSTER_FULLY_CONNECTED_LATTICE_HPP

#include <boost/tuple/tuple.hpp>
#include <vector>

namespace cluster {

class fully_connected_lattice {
public:
  fully_connected_lattice(unsigned int N) : num_sites_(N) {}
  unsigned int num_sites() const { return num_sites_; }
  unsigned int num_bonds() const { return num_sites_ * (num_sites_ - 1) / 2; }
private:
  unsigned int num_sites_;
};

} // end namespace cluster

#endif // CLUSTER_FULLY_CONNECTED_LATTICE_HPP
