#include "utils.hpp"
#include "bitset.hpp"

#include <chrono>
#include <cassert>

namespace triangulator {
namespace utils {
// Is a subset of b? (or a == b)
bool Subsumes(const std::vector<int>& a, const std::vector<int>& b) {
  if (a.size() > b.size()) return false;
  int i2 = 0;
  for (int x : a) {
    while (i2 < (int)b.size() && b[i2] < x) {
      i2++;
    }
    if (i2 >= (int)b.size() || b[i2] > x) return false;
  }
  return true;
}

Bitset ToBitset(const std::vector<int>& a, int n) {
  assert(n>=0);
  Bitset bs(n);
  for (int x : a) {
    assert(x>=0&&x<n);
    bs.SetTrue(x);
  }
  return bs;
}

} // namespace utils

int Log::log_level_ = 10000;
void Log::SetLogLevel(int lvl) {
  log_level_ = lvl;
}

Timer::Timer() {
  timing = false;
  elapsedTime = std::chrono::duration<double>(std::chrono::duration_values<double>::zero());
}

void Timer::start() {
  if (timing) return;
  timing = true;
  startTime = std::chrono::steady_clock::now();
}

void Timer::stop() {
  if (!timing) return;
  timing = false;
  std::chrono::time_point<std::chrono::steady_clock> endTime = std::chrono::steady_clock::now();
  elapsedTime += (endTime - startTime);
}

std::chrono::duration<double> Timer::getTime() {
  if (timing) {
    stop();
    std::chrono::duration<double> ret = elapsedTime;
    start();
    return ret;
  }
  else {
    return elapsedTime;
  }
}

UniqQue::UniqQue(int n) : inq_(n), n_(n) {}
void UniqQue::Add(int x) {
  assert(x>=0 && x<n_);
  if (!inq_[x]) {
    q_.push(x);
    inq_[x] = true;
  }
}
void UniqQue::Add(const std::vector<int>& xs) {
  for (int x : xs) {
    Add(x);
  }
}
int UniqQue::Pop() {
  assert(q_.size()>0);
  int r = q_.front();
  assert(inq_[r]);
  inq_[r] = false;
  q_.pop();
  return r;
}
bool UniqQue::Empty() const {
  return q_.empty();
}
} // namespace triangulator
