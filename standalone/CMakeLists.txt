#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

foreach(name infinite percolation_bond percolation_site)
  add_executable(${name}_sa ${name}.cpp)
  add_alps_test(${name}_sa)
endforeach(name)
