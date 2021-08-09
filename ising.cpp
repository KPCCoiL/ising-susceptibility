#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <numeric>
#include "ising.hpp"

IsingModel::IsingModel(std::size_t num_spin, double J,
                       std::vector<std::vector<IsingModel::Site>> const& neighbors)
  : num_spin(num_spin), spin_up(num_spin), total_spin(-num_spin), interaction(J), neighbors(neighbors) {
  same_direction = 0;
  for (std::size_t i = 0; i < num_spin; i++) {
    same_direction += neighbors[i].size();
  }
  same_direction /= 2;
}

std::vector<std::vector<IsingModel::Site>> construct_lattice(std::size_t rows, std::size_t columns) {
  std::vector neighbors(rows * columns, std::vector<IsingModel::Site>());
  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = 0; j < columns; j++) {
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
  total_spin += spin_up[site] ? 2 : -2;
  for (auto neighbor : neighbors[site]) {
    same_direction += spin_up[site] == spin_up[neighbor] ? 2 : -2;
  }
}

int IsingModel::magnetization() const {
  return total_spin;
}

double susceptibility(IsingModel& ising, double temperature, long long steps) {
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

  long long discard = steps / 10;
  long long sum_M = 0, sum_Msquare = 0;
  double energy = ising.energy();
  std::uniform_int_distribution<> choose_spin(0, spins - 1);
  int acceptance = 0;

  for (long long step = 0; step < steps; step++) {
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
      auto M = ising.magnetization();
      sum_M += M;
      sum_Msquare += M * M;
    }
  }

  int measured_steps = steps - discard;
  std::clog << (static_cast<double>(acceptance) / steps) << std::endl;
  double mean_M = static_cast<double>(sum_M) / static_cast<double>(measured_steps),
    mean_Msquare = static_cast<double>(sum_Msquare) / static_cast<double>(measured_steps),
    magnetization_variance = mean_Msquare - mean_M * mean_M;
  return magnetization_variance * beta / ising.count_spins();
}
