#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <numeric>
#include "ising.hpp"

IsingModel::IsingModel(std::size_t num_spin, double J,
                       std::vector<std::vector<IsingModel::Site>> const& neighbors)
  : num_spin(num_spin), spin_up(num_spin), interaction(J), neighbors(neighbors) {
  same_direction = 0;
  for (int i = 0; i < num_spin; i++) {
    same_direction += neighbors[i].size();
  }
  same_direction /= 2;
}

std::vector<std::vector<IsingModel::Site>> construct_lattice(std::size_t rows, std::size_t columns) {
  std::vector neighbors(rows * columns, std::vector<IsingModel::Site>());
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      int pos = i * columns + j;
      std::array<std::pair<int, int>, 4> drs = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
      for (auto [dx, dy] : drs) {
        int newx = (i + dx + rows) % rows,
          newy = (j + dy + columns) % columns,
          newpos = newx * columns + newy;
        neighbors[pos].push_back(newpos);
      }
    }
  }
  return neighbors;
}

IsingModel::IsingModel(std::size_t rows, std::size_t columns, double J)
  : IsingModel(rows * columns, J, construct_lattice(rows, columns)) {}

void IsingModel::flip_spin(std::size_t site) {
  spin_up[site].flip();
  for (auto neighbor : neighbors[site]) {
    same_direction += spin_up[site] == spin_up[neighbor] ? 2 : -2;
  }
}

double IsingModel::magnetization_square() const {
  int M = 0;
  for (int i = 0; i < num_spin; i++) {
    M += spin_up[i] ? 1 : -1;
  }
  return static_cast<double>(M * M) / static_cast<double>(num_spin * num_spin);
}

IsingQuantities expectation_values(IsingModel& ising, double temperature, int steps) {
  // constexpr double Boltzmann = 1.386503e-23;
  constexpr double Boltzmann = 1;
  double const beta = 1 / (Boltzmann * temperature);
  int spins = ising.count_spins();
  std::mt19937 rng(std::random_device{}());
  std::bernoulli_distribution coin;
  for (int i = 0; i < spins; i++) {
    if (coin(rng)) {
      ising.flip_spin(i);
    }
  }

  int discard = steps / 10;
  double msquare_sum = 0;
  std::vector<double> energys;
  double energy = ising.energy();
  std::uniform_int_distribution<> choose_spin(0, spins - 1);
  int acceptance = 0;

  for (int step = 0; step < steps; step++) {
    int site = choose_spin(rng);
    ising.flip_spin(site);
    double new_energy = ising.energy(),
      accept = std::min(1., std::exp(-beta * (new_energy - energy)));
    if (std::bernoulli_distribution(accept)(rng)) {
      energy = new_energy;
      acceptance++;
    }
    else {
      ising.flip_spin(site);
    }

    if (step >= discard) {
      energys.push_back(ising.energy());
      msquare_sum += ising.magnetization_square();
    }
  }

  int measured_steps = steps - discard;
  std::clog << (static_cast<double>(acceptance) / steps) << std::endl;
  double mean_energy = std::accumulate(std::begin(energys), std::end(energys), 0.0) / measured_steps;
  double energy_variance = 0;
  for (auto e : energys) {
    double diff = e - mean_energy;
    energy_variance += diff * diff;
  }
  energy_variance /= measured_steps;
  return {
    mean_energy,
    msquare_sum / measured_steps,
    energy_variance / (spins * temperature * temperature)
  };
}
