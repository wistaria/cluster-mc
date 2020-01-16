#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(OPENMP_FOUND)
  set(PROGS parallel parallel_noweight)
  foreach(name ${PROGS})
    add_executable(${name} ${name}.cpp)
    add_test(${name} ${name})
  endforeach(name)
endif(OPENMP_FOUND)
