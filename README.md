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

* Configure & make

  cmake -DALPS_ROOT_DIR /where_to_alps
  make

Directories
-----------

* alps: contains applications that use the ALPS Libraries
    (lattice, model, scheudler, etc)

* standalone: contains standalone applications

* loop: contains common header files

Algorithms
----------

* infinite: O(N) Swendsen-Wang Cluster Algorithm for Infinite Ragnge Ising Model

References
----------

* Synge Todo, Kiyoshi Kato, Cluster Algorithms for General-S Quantum Spin System
s, Phys. Rev. Lett. 87, 047203 (2001).

* B. Bauer, L. D. Carr, A. Feiguin, J. Freire, S. Fuchs, L. Gamper, J. Gukelberger, E. Gull, S. Guertler, A. Hehn, R. Igarashi, S.V. Isakov, D. Koop, P.N. Ma, P. Mates, H. Matsuo, O. Parcollet, G. Pawlowski, J.D. Picon, L. Pollet, E. Santos, V.W. Scarola, U. Schollwoeck, C. Silva, B. Surer, S. Todo, S. Trebst, M. Troyer, M.L. Wall, P. Werner, S. Wessel, The ALPS project release 2.0: Open source software for strongly correlated systems, J. Stat. Mech. P05001 (2011).

* Synge Todo, Loop Algorithm, in Strongly Correlated Systems: Numerical Methods (Springer Series in Solid-State Sciences), ed. A. Avella, F. Mancini, pp. 153-184 (Springer-Verlag, Berlin, 2013).
