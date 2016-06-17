/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef CLUSTER_OBSERVABLE_HPP
#define CLUSTER_OBSERVABLE_HPP

#include <cmath> // for std::sqrt
#include <string>

namespace cluster {

class observable {
public:
  observable(std::string const& name) : name_(name), count_(0), sum_(0), esq_(0) {}
  void operator<<(double x) { sum_ += x; esq_ += x * x; ++count_; }
  double mean() const { return (count_ > 0) ? (sum_ / count_) : 0.; }
  double error() const {
    return (count_ > 1) ?
      std::sqrt((esq_ / count_ - mean() * mean()) / (count_ - 1)) : 0.;
  }
  std::ostream& write(std::ostream& os) const {
    os << name_ << " = " << this->mean() << " +- " << this->error();
    return os;
  }
private:
  std::string name_;
  unsigned int count_;
  double sum_, esq_;
};

std::ostream& operator<<(std::ostream& os, observable const& obs) {
  return obs.write(os);
}
  
} // end namespace cluster

#endif // CLUSTER_OBSERVABLE_HPP
