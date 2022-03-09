#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include <utility>
#include <vector>

#ifndef SIGMAOPTS_H
#define SIGMAOPTS_H

class SIGMAINDICES {
  // sigma indexes w, n1 and n2: <psi_n1|sigma(w)|psi_n2>
 public:
  double w;
  std::vector<std::pair<int, int>> n12;  // n1 and n2
  SIGMAINDICES(double _w, int _n1, int _n2) {
    w = _w;
    std::pair<int, int> new_n12(_n1, _n2);
    n12.push_back(new_n12);
  }

  void add_mtrxel(int _n1, int _n2) {
    std::pair<int, int> new_n12(_n1, _n2);
    n12.push_back(new_n12);
  }
  
};

#endif