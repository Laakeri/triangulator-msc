#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <queue>

#include "bitset.hpp"

namespace triangulator {
// Interface

namespace utils {
template<typename T>
bool IsSorted(const std::vector<T>& vec);

template<typename T>
void SortAndDedup(std::vector<T>& vec);

template<typename T>
void InitZero(std::vector<T>& vec, size_t size);

template<typename T>
std::vector<std::pair<T, T> > CompleteEdges(const std::vector<T>& vertices);

template<typename T>
std::vector<T> PermInverse(const std::vector<T>& perm);

template<typename... Args>
void Warning(Args... message);

template<typename... Args>
void ErrorDie(Args... message);

// Is a subset of b? (or a == b)
bool Subsumes(const std::vector<int>& a, const std::vector<int>& b);
Bitset ToBitset(const std::vector<int>& a, int n);

template<typename T>
bool BS(const std::vector<T>& a, const T& x);

template<typename T>
int BSFind(const std::vector<T>& a, const T& x);

template<typename T>
std::vector<T> Union(const std::vector<T>& a, const std::vector<T>& b);

template<typename T>
T MaxInSub(const std::vector<std::vector<T>>& a, const std::vector<int>& idx);
} // namespace utils

class Log {
 public:
  template<typename T, typename... Args>
  static void P(T first_message, Args... message);
  
  template<typename T>
  static void Write(int lvl, T message);

  template<typename T, typename... Args>
  static void Write(int lvl, T first_message, Args... message);

  static void SetLogLevel(int lvl);
 private:
  template<typename T>
  static void WriteImpl(std::vector<T> message);

  template<typename T>
  static void WriteImpl(std::pair<T, T> message);
  
  template<typename T>
  static void WriteImpl(T message);
  static int log_level_;
};

class Timer {
 private:
  bool timing;
  std::chrono::duration<double> elapsedTime;
  std::chrono::time_point<std::chrono::steady_clock> startTime;
 public:
  Timer();
  void start();
  void stop();
  std::chrono::duration<double> getTime();
};

class UniqQue {
 private:
  std::queue<int> q_;
  std::vector<int> inq_;
  int n_;
 public:
  explicit UniqQue(int n);
  void Add(int x);
  void Add(const std::vector<int>& xs);
  int Pop();
  bool Empty() const;
};

// Implementation

namespace utils {
template<typename T>
bool IsSorted(const std::vector<T>& vec) {
  for (int i = 1; i < (int)vec.size(); i++) {
    if (vec[i] < vec[i-1]) return false;
  }
  return true;
}

template<typename T>
void SortAndDedup(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

template<typename T>
void InitZero(std::vector<T>& vec, size_t size) {
  vec.resize(size);
  std::fill(vec.begin(), vec.begin() + size, 0);
}

template<typename T>
T MaxInSub(const std::vector<std::vector<T>>& a, const std::vector<int>& idx) {
  for (int i : idx) {
    assert(i < (int)a.size());
  }
  T res = 0;
  for (int i = 0; i < (int)idx.size(); i++) {
    for (int ii = i+1; ii < (int)idx.size(); ii++) {
      res = std::max(res, a[idx[i]][idx[ii]]);
    }
  }
  return res;
}

template<typename T>
std::vector<T> Union(const std::vector<T>& a, const std::vector<T>& b) {
  std::vector<T> ret;
  ret.reserve(a.size() + b.size());
  for (const T& x : a) {
    ret.push_back(x);
  }
  for (const T& x : b) {
    ret.push_back(x);
  }
  SortAndDedup(ret);
  return ret;
}

template<typename T>
void Append(std::vector<T>& a, const std::vector<T>& b) {
  a.reserve(a.size() + b.size());
  for (const T& x : b) {
    a.push_back(x);
  }
}

template<typename T>
bool BS(const std::vector<T>& a, const T& x) {
  return std::binary_search(a.begin(), a.end(), x);
}

template<typename T>
int BSFind(const std::vector<T>& a, const T& x) {
  auto it = std::lower_bound(a.begin(), a.end(), x);
  if (it == a.end()) return -1;
  int i = it - a.begin();
  if (a[i] == x) return i;
  else return -1;
}

template<typename T>
std::vector<std::pair<T, T> > CompleteEdges(const std::vector<T>& vertices) {
  if (vertices.size() < 2) return { };
  std::vector<std::pair<T, T> > edges;
  edges.reserve(vertices.size() * (vertices.size() - 1) / 2); // care with unsigned
  for (int i = 0; i < vertices.size(); i++) {
    for (int ii = i + 1; ii < vertices.size(); ii++) {
      edges.push_back({vertices[i], vertices[ii]});
    }
  }
  return edges;
}

template<typename T>
std::vector<T> PermInverse(const std::vector<T>& perm) {
  std::vector<T> ret(perm.size());
  for (int i = 0; i < (int)perm.size(); i++) {
    ret[perm[i]] = i;
  }
  return ret;
}

template<typename... Args>
void Warning(Args... message) {
  Log::Write(5, "[WARNING] ", message...);
}
// TODO: this is bad idea because assert shows line number
template<typename... Args>
void ErrorDie(Args... message) {
  Log::Write(1, "\033[0;31m[FATAL ERROR]\033[0m ", message...);
  std::exit(EXIT_FAILURE);
}
} // namespace utils
template<typename T>
void Log::WriteImpl(T message) {
  std::cerr<<message;
}
template<>
inline void Log::WriteImpl(std::vector<char> message) {
  std::cerr<<"{ ";
  for (char e : message) {
    std::cerr<<(int)e<<" ";
  }
  std::cerr<<"}";
}
template<>
inline void Log::WriteImpl(std::vector<int> message) {
  std::cerr<<"{";
  for (int i=0;i<(int)message.size();i++){
    std::cerr<<message[i];
    if (i+1<(int)message.size()) {
      std::cerr<<", ";
    }
  }
  std::cerr<<"}";
}
template<typename T>
inline void Log::WriteImpl(std::pair<T, T> message) {
  std::cerr<<"{";
  WriteImpl(message.first);
  std::cerr<<", ";
  WriteImpl(message.second);
  std::cerr<<"}";
}
template<typename T>
inline void Log::WriteImpl(std::vector<T> message) {
  std::cerr<<"{";
  for (int i=0;i<(int)message.size();i++){
    WriteImpl(message[i]);
    if (i+1<(int)message.size()) {
      std::cerr<<", ";
    }
  }
  std::cerr<<"}";
}
template<>
inline void Log::WriteImpl(Bitset b) {
  std::cerr<<"{";
  for (size_t i=0;i<b.chunks_*BITS;i++){
    std::cerr<<b.Get(i);
  }
  std::cerr<<"}";
}
template<typename T>
void Log::Write(int lvl, T message) {
  if (lvl > log_level_) return;
  WriteImpl(message);
  std::cerr<<std::endl;
}
template<typename T, typename... Args>
void Log::Write(int lvl, T first_message, Args... message) {
  if (lvl > log_level_) return;
  WriteImpl(first_message);
  Write(lvl, message...);
}
template<typename T, typename... Args>
void Log::P(T first_message, Args... message) {
  Write(0, first_message, message...);
}
} // namespace triangulator
