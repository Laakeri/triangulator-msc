#pragma once

#include <vector>

#include "utils.hpp"

namespace triangulator {
// Interface
template<typename T>
class StaticSet {
public:
  StaticSet();
  explicit StaticSet(const std::vector<T>& values);
  explicit StaticSet(const std::vector<std::pair<T, T> >& values);
  void Init(const std::vector<T>& values);
  void Init(const std::vector<std::pair<T, T> >& values);
  int Rank(T value) const;
  T Kth(int k) const;
  int Size() const;
  std::vector<T> Values() const;
private:
  std::vector<T> values_;
};


// Implementation
template<typename T>
StaticSet<T>::StaticSet(const std::vector<T>& values) {
  Init(values);
}

template<typename T>
StaticSet<T>::StaticSet() : StaticSet(std::vector<T>()) { }

template<typename T>
StaticSet<T>::StaticSet(const std::vector<std::pair<T, T> >& values) {
  Init(values);
}

template<typename T>
void StaticSet<T>::Init(const std::vector<T>& values) {
  values_ = values;
  utils::SortAndDedup(values_);
}

template<typename T>
void StaticSet<T>::Init(const std::vector<std::pair<T, T> >& values) {
  values_.clear();
  for (const std::pair<T, T>& value : values) {
    values_.push_back(value.first);
    values_.push_back(value.second);
  }
  utils::SortAndDedup(values_);
}

template<typename T>
int StaticSet<T>::Rank(T value) const {
  return std::lower_bound(values_.begin(), values_.end(), value) - values_.begin();
}

template<typename T>
T StaticSet<T>::Kth(int k) const {
  return values_[k];
}

template<typename T>
int StaticSet<T>::Size() const {
  return values_.size();
}

template<typename T>
std::vector<T> StaticSet<T>::Values() const {
  return values_;
}
} // namespace triangulator
