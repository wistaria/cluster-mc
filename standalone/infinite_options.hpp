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
  unsigned int seed, num_sites;
  double temperature;
  unsigned int sweeps, therm;
  bool valid;
  options(unsigned int argc, char *argv[], bool print = true) :
    seed(29833), num_sites(128), temperature(1.0), sweeps(1 << 16), therm(sweeps >> 3),
    valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 's' :
          if (++i == argc) { usage(print); return; }
          seed = std::atoi(argv[i]); break;
        case 'n' :
          if (++i == argc) { usage(print); return; }
          num_sites = std::atoi(argv[i]); break;
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
    if (num_sites == 0 || temperature <= 0 || sweeps == 0) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "Number of Sites        = " << num_sites << std::endl
                << "Temperature            = " << temperature << std::endl
                << "MCS for Thermalization = " << therm << std::endl
                << "MCS for Measurement    = " << sweeps << std::endl;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -n int    Number of Sites\n"
         << "  -t double Temperature\n"
         << "  -m int    MCS for Measurement\n"
         << "  -h        this help\n\n";
    valid = false;
  }
};
