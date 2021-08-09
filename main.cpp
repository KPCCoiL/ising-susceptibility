#include <iostream>
#include <sstream>
#include "ising.hpp"

int main(int argc, char** argv) {
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
  double T_from = get_param(1, 1.0);
  double T_to = get_param(2, 10.0);
  double T_by = get_param(3, 1.0);
  int L = get_param(4, 4);
  int steps = get_param(5, 100000);
  double J = get_param(6, 1.0);
  IsingModel ising(L, L, J);
  for (double T = T_from; T < T_to; T += T_by) {
    IsingQuantities qtys = expectation_values(ising, T, steps);
    std::cout << T << ' '
              << qtys.energy << ' '
              << qtys.magnetization_square << ' '
              << qtys.specific_heat << '\n';
  }
}
