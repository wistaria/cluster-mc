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

#ifndef CLUSTER_SQUARE_LATTICE_HPP
#define CLUSTER_SQUARE_LATTICE_HPP

#include <vector>

namespace cluster {

class square_lattice {
public:
  square_lattice(unsigned int L) : length_x_(L), length_y_(L) { init(); }
  square_lattice(unsigned int Lx, unsigned int Ly) : length_x_(Lx), length_y_(Ly) { init(); }
  void init() {
    source_.resize(2 * length_x_ * length_y_);
    target_.resize(2 * length_x_ * length_y_);
    site_phase_.resize(length_x_ * length_y_);
    bond_phase_.resize(2 * length_x_ * length_y_);
    unsigned int n = length_x_ * length_y_ * 2;
    for (unsigned int s = 0; s < length_x_ * length_y_; ++s) {
      int ix = s % length_x_;
      int iy = s / length_y_;
      site_phase_[s] = 2 * ((ix + iy) % 2) - 1;
    }
    for (unsigned int b = 0; b < n; ++b) {
      unsigned int s = b / 2;
      unsigned int t;
      if (b % 2 == 0) {
        t = (s + 1) % length_x_ + (s / length_x_) * length_x_; // target right
        bond_phase_[b] = 2.0 * ((b / 2) % 2) - 1.0;
      } else {
        t = (s + length_x_) % (length_x_ * length_y_); // target below
        bond_phase_[b] = 2.0 * ((b / length_x_ / 2) % 2) - 1.0;
      }
      source_[b] = s;
      target_[b] = t;
    }
  }
  unsigned int get_length_x() const { return length_x_; }
  unsigned int get_length_y() const { return length_y_; }
  unsigned int num_sites() const { return length_x_ * length_y_; }
  unsigned int num_bonds() const { return 2 * num_sites(); }
  unsigned int source(unsigned int b) const { return source_[b]; }
  unsigned int target(unsigned int b) const { return target_[b]; }
  double site_phase(unsigned int s) const { return site_phase_[s]; }
  double bond_phase(unsigned int b) const { return bond_phase_[b]; }
private:
  unsigned int length_x_, length_y_;
  std::vector<unsigned int> source_;
  std::vector<unsigned int> target_;
  std::vector<double> site_phase_;
  std::vector<double> bond_phase_;
};

} // end namespace cluster

#endif // CLUSTER_SQUARE_LATTICE_HPP
