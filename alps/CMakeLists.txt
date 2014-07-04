#
# Cluster-MC: Cluster Algorithms Monte Carlo Methods
#
# Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#

# O(N) Algorithm for Infinite Range Ising Model
add_executable(infinite infinite.C)
target_link_libraries(infinite ${ALPS_LIBRARIES})
add_alps_test(infinite)
