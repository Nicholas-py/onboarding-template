#pragma once

#include <cstddef>
#include <vector>

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid
{
private:

public:
  Grid(std::size_t rows, std::size_t cols);
  std::size_t rows_;
  std::size_t cols_;
  std::vector<std::vector<double>> array;

  double &operator()(std::size_t i, std::size_t j);
  double operator()(std::size_t i, std::size_t j) const;
};

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid &old_grid, Grid &new_grid);

Grid::Grid(std::size_t rows, std::size_t cols) {
  this->cols_ = cols;
  this->rows_ = rows;
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


  int rows = old_grid.rows_;
  int cols = old_grid.cols_;
  std::vector<std::vector<double>> oldarr = old_grid.array;
  std::vector<std::vector<double>> newarr = new_grid.array;
  for (int i = 1; i < rows-1; i++) {
    for (int j = 1; j < cols-1; j++) {
      newarr[i][j] = 0.5*(oldarr[i][j]) + 0.125*(oldarr[i-1][j] + oldarr[i+1][j] + oldarr[i][j-1] + oldarr[i][j+1]);
    }
  }
  for (int i = 0; i < rows; i++) {
    newarr[i][0] = oldarr[i][0];
    newarr[i][cols-1] = oldarr[i][cols-1];
  }
  for (int i = 1; i < cols-1; i++) {
    newarr[0][i] = oldarr[0][i];
    newarr[rows-1][i] = oldarr[rows-1][i];

  }
  new_grid.array = newarr;

}