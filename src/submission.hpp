#pragma once

#include <cstddef>
#include <vector>
#pragma GCC optimize("Ofast")

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid
{
private:
  std::size_t rows_;
  std::size_t cols_;
  std::size_t maxrow_;
  std::size_t maxcol_;

public:
  Grid(std::size_t rows, std::size_t cols);

  std::vector<std::vector<double>> array;

  double &operator()(std::size_t i, std::size_t j);
  double operator()(std::size_t i, std::size_t j) const;

  int rows() const {
    return rows_;
  }
  int cols() const {
    return cols_;
  }
  int maxrow() const {
    return maxrow_;
  }
  int maxcol() const {
    return maxcol_;
  }

  void setrowcol(std::size_t rows, std::size_t cols) {
    rows_ = rows;
    cols_ = cols;
    maxrow_ = rows-1;
    maxcol_ = cols-1;
  }
};

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid &old_grid, Grid &new_grid);

Grid::Grid(std::size_t rows, std::size_t cols) {
  setrowcol(rows, cols);
  
  std::vector<std::vector<double>> data(rows, std::vector<double>(cols, 0));
  this->array = data;
}

inline double& Grid::operator()(std::size_t i, std::size_t j) {
  return this->array[i][j];
}

inline double Grid::operator()(std::size_t i, std::size_t j) const {
  return this->array[i][j];
}

inline void apply_stencil(const Grid &old_grid, Grid &new_grid) {

  #pragma omp simd
  for (int i = 1; i < old_grid.maxrow(); i++) {
    for (int j = 1; j < old_grid.maxcol(); j++) {
      new_grid.array[i][j] = 0.5*(old_grid.array[i][j]) +0.125*(old_grid.array[i-1][j] + old_grid.array[i+1][j] + old_grid.array[i][j-1] + old_grid.array[i][j+1]);
    }
  }

  for (int i = 0; i < old_grid.rows(); i++) {
    new_grid.array[i][0] = old_grid.array[i][0];
    new_grid.array[i][old_grid.maxcol()] = old_grid.array[i][old_grid.maxcol()];
  }
  for (int i = 1; i < old_grid.maxcol(); i++) {
    new_grid.array[0][i] = old_grid.array[0][i];
    new_grid.array[old_grid.maxrow()][i] = old_grid.array[old_grid.maxrow()][i];
  }
}