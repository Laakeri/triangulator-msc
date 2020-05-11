#pragma once

#include <vector>
#include <cstdlib>

namespace triangulator {
// Interface
template<typename T>
class Matrix {
 public:
  Matrix();
  explicit Matrix(size_t rows, size_t columns);
  void Resize(size_t rows, size_t columns);
  T* operator[](size_t pos);
  const T* operator[](size_t pos) const;
  size_t Rows() const;
  size_t Columns() const;
 private:
  size_t rows_, columns_;
  std::vector<T> data_;
};

// Implementation
template<typename T>
Matrix<T>::Matrix(size_t rows, size_t columns) : rows_(rows), columns_(columns), data_(rows * columns) { }

template<typename T>
Matrix<T>::Matrix() : Matrix(0, 0) {}

template<typename T>
void Matrix<T>::Resize(size_t rows, size_t columns) {
  rows_ = rows;
  columns_ = columns;
  data_.resize(rows * columns);
}

template<typename T>
T* Matrix<T>::operator[](size_t pos) {
  return data_.data() + pos * columns_;
}

template<typename T>
const T* Matrix<T>::operator[](size_t pos) const {
  return data_.data() + pos * columns_;
}

template<typename T>
size_t Matrix<T>::Rows() const {
  return rows_;
}

template<typename T>
size_t Matrix<T>::Columns() const {
  return columns_;
}

template<typename T>
class VVector {
 public:
  VVector();
  explicit VVector(const std::vector<std::vector<T>>& data);
  size_t Size(size_t i) const;
  T* operator[](size_t i);
  const T* operator[](size_t i) const;
 private:
  std::vector<T> data_;
  std::vector<size_t> pos_;
};

template<typename T>
VVector<T>::VVector(const std::vector<std::vector<T>>& data) {
  pos_.resize(data.size() + 1);
  for (int i = 0; i < (int)data.size(); i++) {
    pos_[i+1] = pos_[i] + data[i].size();
  }
  data_.resize(pos_.back());
  for (int i = 0; i < (int)data.size(); i++) {
    for (int j = 0; j < (int)data[i].size(); j++) {
      data_[pos_[i] + j] = data[i][j];
    }
  }
}

template<typename T>
VVector<T>::VVector() : VVector(std::vector<std::vector<T>>({})) {}

template<typename T>
size_t VVector<T>::Size(size_t i) const {
  return pos_[i+1] - pos_[i];
}

template<typename T>
T* VVector<T>::operator[](size_t i) {
  return data_.data() + pos_[i];
}

template<typename T>
const T* VVector<T>::operator[](size_t i) const {
  return data_.data() + pos_[i];
}


} // namespace triangulator
