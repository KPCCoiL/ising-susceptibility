#include <iostream>
#include <numeric>
#include <sstream>
#include <fstream>
#include <string>
#include "ising.hpp"

int main(int argc, char** argv) {
  using namespace std::literals;
  auto get_param = [argc, argv](int index, auto def) {
    if (argc <= index) {
      return def;
    }
    std::stringstream ss;
    ss << argv[index];
    decltype(def) res;
    ss >> res;
    return res;
  };
  auto outfile = get_param(1, "result.txt"s);
  long long steps = get_param(2, 100000000ll);
  int repetition = get_param(3, 1);
  double J = get_param(4, 1.0);
  const double T_c = 2.26918531421 * J;
  std::ofstream result(outfile);

  constexpr int MaxL = 30;
  std::vector<int> Ls(MaxL - 1);
  for (int L = 2; L <= MaxL; L++) {
    Ls[L - 2] = L;
  }
  
  const int N = Ls.size();
  std::vector<double> chis(N), stds(N);
  
  #pragma omp parallel for
  for (int i = 0; i < N; i++) {
    IsingModel ising(Ls[i], Ls[i], J);
    std::vector<double> res(repetition);
    for (int j = 0; j < repetition; j++) {
      res[j] = susceptibility(ising, T_c, steps);
    }
    double mean = std::accumulate(std::begin(res), std::end(res), 0.0) / repetition;
    double var = 0;
    for (auto chi : res) {
      auto diff = chi - mean;
      var += diff * diff;
    }
    var /= repetition;
    chis[i] = mean;
    stds[i] = std::sqrt(var);
  }
  for (int i = 0; i < N; i++) {
    result << Ls[i] << ' ' << chis[i] << ' ' << stds[i] << '\n';
  }
}
