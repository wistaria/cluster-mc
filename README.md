Cluster-MC: Cluster Algorithm Monte Carlo Methods
=================================================

Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

Install
-------

* Prerequisites
    * ALPS Library

* CMake Options
    * ALPS_ROOT_DIR: path to ALPS library
    * BUILD_ALPS_APPLICATIONS: (default ON)
    * BUILD_STANDALONE_APPLICATIONS: (default ON)

Directories
-----------

* alps: contains applications that use the ALPS Libraries
    (lattice, model, scheudler, etc)

* standalone: contains standalone applications

* loop: contains common header files

Algorithms
----------

* infinite: O(N) Swendsen-Wang Cluster Algorithm for Infinite Ragnge Ising Model
