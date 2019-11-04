# Cluster-MC: Cluster Algorithm Monte Carlo Methods

Copyright (C) 1997-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

## Install

* Prerequisites
    * Eigen3
    * ALPS Library (optional)
* CMake Options
    * ALPS_ROOT_DIR: path to ALPS library
    * BUILD_ALPS_APPLICATIONS: (default ON)
    * BUILD_STANDALONE_APPLICATIONS: (default ON)
* Build standalone programs only
    ```
    mkdir build
    cd build
    cmake ..
    make
    ```
* Build ALPS and standalone programs
    ```
    mkdir build
    cd build
    cmake -DALPS_ROOT_DIR=/where_to_alps ..
      (or if environment variable ALPS_ROOT is set, just type "cmake ..")
    make
    ```

## Directories

* alps: contains applications that use the ALPS Libraries (lattice, model, scheudler, etc)
* standalone: contains standalone applications
* cluster: contains common header files
* tool/standards: from standards library https://github.com/todo-group/standards
* tool/lattice: from lattice library https://github.com/todo-group/lattice

## Algorithms

* infinite: O(N) Swendsen-Wang Cluster Algorithm for Infinite Ragnge Ising Model
* ising: Swendsen-Wang Cluster Algorithm for Ising Model (standalone version only)
* potts: Swendsen-Wang Cluster Algorithm for Potts Model
* loop_*: Loop Algorithm for Spin-1/2 Antiferromagnetic Heisenberg Chain (standalone version only)
   * loop_pi0: continuous time path integral; using std::list<> for operator string
   * loop_pi1: continuous time path integral; using std::vector<> for operator string
   * loop_fsse: fixed-length SSE
   * loop_vsse: variable-length SSE

## Release Note

* Release 0.3
    * ising: Swendsen-Wang Cluster Algorithm for Ising Model (standalone version)
    * loop_*: Loop Algorithm for Spin-1/2 Antiferromagnetic Heisenberg Chain (standalone version)
* Release 0.2
    * potts: Swendsen-Wang Cluster Algorithm for Potts Model
    * standalone version of applications
    * site and bond percolation problems
* Release 0.1 (2014/07/05)
    * initial version
    * infinite: O(N) Swendsen-Wang Cluster Algorithm for Infinite Ragnge Ising Model

## License

* Distributed under Boost Software License - Version 1.0 - August 17th, 2003

    Permission is hereby granted, free of charge, to any person or organization obtaining a copy of the software and accompanying documentation covered by this license (the "Software") to use, reproduce, display, distribute, execute, and transmit the Software, and to prepare derivative works of the Software, and to permit third-parties to whom the Software is furnished to do so, all subject to the following:
    
    The copyright notices in the Software and this entire statement, including the above license grant, this restriction and the following disclaimer, must be included in all copies of the Software, in whole or in part, and all derivative works of the Software, unless such copies or derivative works are solely in the form of machine-executable object code generated by a source language processor.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
## References

* _Cluster Algorithms for General-S Quantum Spin System,_ S. Todo, K. Kato, Phys. Rev. Lett. 87, 047203 (2001).
* _The ALPS project release 2.0: Open source software for strongly correlated systems,_ B. Bauer, L. D. Carr, A. Feiguin, J. Freire, S. Fuchs, L. Gamper, J. Gukelberger, E. Gull, S. Guertler, A. Hehn, R. Igarashi, S.V. Isakov, D. Koop, P.N. Ma, P. Mates, H. Matsuo, O. Parcollet, G. Pawlowski, J.D. Picon, L. Pollet, E. Santos, V.W. Scarola, U. Schollwoeck, C. Silva, B. Surer, S. Todo, S. Trebst, M. Troyer, M.L. Wall, P. Werner, S. Wessel, J. Stat. Mech. P05001 (2011).
* _Loop Algorithm,_ S. Todo, in _Strongly Correlated Systems: Numerical Methods_ (Springer Series in Solid-State Sciences), ed. A. Avella, F. Mancini, pp. 153-184 (Springer-Verlag, Berlin, 2013).
* ALPS: [http://alps.comp-phys.org/](http://alps.comp-phys.org/).
* ALPS/looper: [http://exa.phys.s.u-tokyo.ac.jp/alps-looper](http://exa.phys.s.u-tokyo.ac.jp/alps-looper)
* For the variable length SSE, see _Optimized broad-histogram ensembles for the simulation of quantum systems,_ S. Wessel, N. Stoop, E. Gull, S. Trebst and M. Troyer, J. Stat. Mech. P12005 (2007).
