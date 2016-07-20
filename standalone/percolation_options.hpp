/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/lexical_cast.hpp>
#include <iostream>

struct options {
  unsigned int seed, length;
  double probability;
  unsigned int sweeps;
  bool valid;

  options(unsigned int argc, char *argv[], double default_probability, bool print = true) :
    seed(29833), length(256), probability(default_probability), sweeps(1 << 8),
    valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 's' :
          if (++i == argc) { usage(print); return; }
          seed = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'l' :
          if (++i == argc) { usage(print); return; }
          length = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'p' :
          if (++i == argc) { usage(print); return; }
          probability = boost::lexical_cast<double>(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(print); return; }
          sweeps = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'h' :
          usage(print, std::cout); return;
        default :
          usage(print); return;
        }
        break;
      default :
        usage(print); return;
      }
    }
    if (length == 0 || probability < 0 || probability > 1 || sweeps == 0) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "System Linear Size     = " << length << std::endl
                << "Occupation Probability = " << probability << std::endl
                << "Monte Carlo Steps      = " << sweeps << std::endl;;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -l int    System Linear Size\n"
         << "  -p double Occupation Probability\n"
         << "  -m int    Monte Carlo Steps\n"
         << "  -h        this help\n";
    valid = false;
  }
};
