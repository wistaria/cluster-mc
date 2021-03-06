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

#include <cstdlib>
#include <iostream>

struct options {
  unsigned int seed, q, length;
  double temperature;
  unsigned int sweeps, therm;
  bool valid;

  options(unsigned int argc, char *argv[], bool print = true) :
    // default parameters
    seed(29833), length(8), temperature(2.27), sweeps(1 << 16), therm(sweeps >> 3),
    valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 's' :
          if (++i == argc) { usage(print); return; }
          seed = std::atoi(argv[i]); break;
        case 'l' :
          if (++i == argc) { usage(print); return; }
          length = std::atoi(argv[i]); break;
        case 't' :
          if (++i == argc) { usage(print); return; }
          temperature = std::atof(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(print); return; }
          sweeps = std::atoi(argv[i]);
          therm = sweeps >> 3; break;
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
    if (length == 0 || temperature <= 0. || sweeps == 0) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "System Linear Size     = " << length << std::endl
                << "Temperature            = " << temperature << std::endl
                << "MCS for Thermalization = " << therm << std::endl
                << "MCS for Measurement    = " << sweeps << std::endl;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -l int    System Linear Size\n"
         << "  -t double Temperature\n"
         << "  -m int    MCS for Measurement\n"
         << "  -h        this help\n";
    valid = false;
  }
};
