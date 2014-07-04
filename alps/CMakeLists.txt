#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# O(N) Swendsen-Wang Algorithm for Infinite Range Ising Model
add_executable(infinite infinite.C)
target_link_libraries(infinite ${ALPS_LIBRARIES})
add_alps_test(infinite)
