# Cluster-MC: Cluster Algorithm Monte Carlo Methods

Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

##Algorithms


* infinite: O(N) Swendsen-Wang Cluster Algorithm for Infinite Ragnge Ising Model
    * Hamiltonian
        H = -\frac{1}{N} \sum_{i<j} \sigma_i \sigma_j
    * paramters
        * N : Number of Sites
        * T : Temperature
    * observables
        * Number of Clusters: average number of clusters
        * Magnetization (unimproved): \sum_i \sigma_i
        * Magnetization^2 (unimproved): (\sum_i \sigma_i)^2
        * Magnetization^4 (unimproved): (\sum_i \sigma_i)^4
        * Magnetization^2: improved estimator for Magnetization^2
        * Magnetization^4: improved estimator for Magnetization^4
        * Binder Ratio of Magnetization (unimproved): Q = (m^2)^2/(m^4)
        * Binder Ratio of Magnetization: improved estimator for Q
