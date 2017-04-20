#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(PROGS infinite ising potts percolation_bond percolation_site loop_pi0 loop_pi1 loop_fsse loop_vsse)
foreach(name ${PROGS})
  set(target_name ${name}_standalone)
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  add_iotest(${name})
endforeach(name)
