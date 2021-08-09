#pragma once

#include <vector>
#include <random>
#include <utility>

class IsingModel {
public:
  using Site = std::size_t;
  IsingModel(std::size_t num_spin, double J, std::vector<std::vector<Site>> const& neighbors);
  IsingModel(std::size_t rows, std::size_t columns, double J);

  void flip_spin(std::size_t site);

  double magnetization_square() const;
  inline double energy() const {
    return -interaction * same_direction;
  }
  inline std::size_t count_spins() const {
    return num_spin;
  }

private:
  std::size_t num_spin;
  std::vector<bool> spin_up;
  double interaction;
  int same_direction;
  std::vector<std::vector<Site>> neighbors;
};

struct IsingQuantities {
  double energy, magnetization_square, specific_heat;
};

IsingQuantities expectation_values(IsingModel& ising, double temperature, int steps);
