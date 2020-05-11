#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <random>
#include <limits>
#include <cassert>

#define BITS 64

namespace triangulator {
class Bitset {
 public:
  uint64_t* data_;
  size_t chunks_;
  Bitset() {
    chunks_ = 0;
    data_ = nullptr;
  }
  explicit Bitset(size_t size) {
    chunks_ = (size + BITS - 1) / BITS;
    data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
    for (size_t i=0;i<chunks_;i++){
      data_[i] = 0;
    }
  }
  ~Bitset() {
    std::free(data_);
  }
  Bitset(const Bitset& other) {
    chunks_ = other.chunks_;
    data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
    for (size_t i=0;i<chunks_;i++){
      data_[i] = other.data_[i];
    }
  }
  Bitset& operator=(const Bitset& other) {
    if (this != &other) {
      if (chunks_ != other.chunks_) {
        std::free(data_);
        chunks_ = other.chunks_;
        data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
      }
      for (size_t i=0;i<chunks_;i++){
        data_[i] = other.data_[i];
      }
    }
    return *this;
  }
  Bitset(Bitset&& other) {
    data_ = other.data_;
    chunks_ = other.chunks_;
    other.data_ = nullptr;
    other.chunks_ = 0;
  }
  Bitset& operator=(Bitset&& other) {
    if (this != &other) {
      std::free(data_);
      data_ = other.data_;
      chunks_ = other.chunks_;
      other.data_ = nullptr;
      other.chunks_ = 0;
    }
    return *this;
  }
  bool operator<(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i]<other.data_[i]) return true;
      else if(data_[i]>other.data_[i]) return false;
    }
    return false;
  }
  bool operator==(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i] != other.data_[i]) return false;
    }
    return true;
  }
  bool operator!=(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i] != other.data_[i]) return true;
    }
    return false;
  }
  Bitset operator|(const Bitset& other) const {
    Bitset ret(BITS*chunks_);
    for (size_t i=0;i<chunks_;i++){
      ret.data_[i] = data_[i] | other.data_[i];
    }
    return ret;
  }
  Bitset operator&(const Bitset& other) const {
    Bitset ret(BITS*chunks_);
    for (size_t i=0;i<chunks_;i++){
      ret.data_[i] = data_[i] & other.data_[i];
    }
    return ret;
  }
  Bitset operator~() const {
    Bitset ret(BITS*chunks_);
    for (size_t i=0;i<chunks_;i++){
      ret.data_[i] = (~data_[i]);
    }
    return ret;
  }
  void CopyFrom(const Bitset& other) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = other.data_[i];
    }
  }
  void Set(size_t i, bool v) {
    if (v) {
      data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
    } else {
      data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
    }
  }
  void SetTrue(size_t i) {
    data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  void SetFalse(size_t i) {
    data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
  }
  void SetTrue(const std::vector<size_t>& v) {
    for (size_t x : v) {
      SetTrue(x);
    }
  }
  void SetTrue(const std::vector<int>& v) {
    for (int x : v) {
      SetTrue(x);
    }
  }
  void SetFalse(const std::vector<int>& v) {
    for (int x : v) {
      SetFalse(x);
    }
  }
  void FillTrue() {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = ~0;
    }
  }
  void FillUpTo(size_t n) {
    for (size_t i=0;i<chunks_;i++){
      if ((i+1)*BITS <= n) {
        data_[i] = ~0;
      } else if (i*BITS < n) {
        for (size_t j=i*BITS;j<n;j++){
          SetTrue(j);
        }
      } else {
        return;
      }
    }
  }
  bool Get(size_t i) const {
    return data_[i/BITS] & ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  void Clear() {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = 0;
    }
  }
  bool IsEmpty() const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i]) return false;
    }
    return true;
  }
  uint64_t Chunk(size_t i) const {
    return data_[i];
  }
  uint64_t& Chunk(size_t i) {
    return data_[i];
  }
  size_t Chunks() const {
    return chunks_;
  }
  void operator |= (const Bitset& rhs) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] |= rhs.data_[i];
    }
  }
  void operator &= (const Bitset& rhs) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] &= rhs.data_[i];
    }
  }
  void TurnOff(const Bitset& rhs) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] &= (~rhs.data_[i]);
    }
  }
  void InvertAnd(const Bitset& rhs) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = (~data_[i]) & rhs.data_[i];
    }
  }
  void SetNeg(const Bitset& rhs) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = ~rhs.data_[i];
    }
  }
  void SetNegAnd(const Bitset& rhs1, const Bitset& rhs2) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = (~rhs1.data_[i]) & rhs2.data_[i];
    }
  }
  void SetAnd(const Bitset& rhs1, const Bitset& rhs2) {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = rhs1.data_[i] & rhs2.data_[i];
    }
  }
  bool Subsumes(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if ((data_[i] | other.data_[i]) != data_[i]) return false;
    }
    return true;
  }
  std::vector<int> Elements() const {
    std::vector<int> ret;
    for (size_t i=0;i<chunks_;i++){
      uint64_t td = data_[i];
      while (td) {
        ret.push_back(i*BITS + __builtin_ctzll(td));
        td &= ~-td;
      }
    }
    return ret;
  }
  int Popcount() const {
    int cnt = 0;
    for (size_t i=0;i<chunks_;i++) {
      cnt += __builtin_popcountll(data_[i]);
    }
    return cnt;
  }
  bool Intersects(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i] & other.data_[i]) return true;
    }
    return false;
  }

  class BitsetIterator {
   private:
    const Bitset* const bitset_;
    size_t pos_;
    uint64_t tb_;
   public:
    BitsetIterator(const Bitset* const bitset, size_t pos, uint64_t tb) : bitset_(bitset), pos_(pos), tb_(tb) { }
    bool operator!=(const BitsetIterator& other) const {
      return pos_ != other.pos_ || tb_ != other.tb_;
    }
    const BitsetIterator& operator++() {
      tb_ &= ~-tb_;
      while (tb_ == 0 && pos_ < bitset_->chunks_) {
        pos_++;
        if (pos_ < bitset_->chunks_) {
          tb_ = bitset_->data_[pos_];
        }
      }
      return *this;
    }
    int operator*() const {
      return pos_*BITS + __builtin_ctzll(tb_);
    }
  };

  BitsetIterator begin() const {
    size_t pos = 0;
    while (pos < chunks_ && data_[pos] == 0) {
      pos++;
    }
    if (pos < chunks_) {
      return BitsetIterator(this, pos, data_[pos]);
    } else {
      return BitsetIterator(this, pos, 0);
    }
  }
  BitsetIterator end() const {
    return BitsetIterator(this, chunks_, 0);
  }
};

class BitsetSet {
 public:
  BitsetSet() {}
  BitsetSet(size_t len, size_t capacity, size_t load_factor) {
    chunks_ = (len + BITS - 1) / BITS;
    load_factor_ = load_factor;
    assert(chunks_ > 0);
    assert(load_factor_ >= 2);
    capacity_ = NextPrime((capacity + 1) * load_factor_);
    container_.resize(capacity_ * chunks_);
  }
  bool Insert(const Bitset& bitset) {
    assert(!bitset.IsEmpty());
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) break;
      else if (Equal(IndToPtr(ind, container_), bitset.data_)) return false;
      else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
    Copy(bitset.data_, IndToPtr(ind, container_));
    elements_++;
    if (elements_ * load_factor_ > capacity_) {
      Resize();
      assert(elements_ * load_factor_ < capacity_);
    }
    return true;
  }
  bool Contains(const Bitset& bitset) const {
    assert(!bitset.IsEmpty());
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) return false;
      else if (Equal(IndToPtr(ind, container_), bitset.data_)) return true;
      else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
  }
  bool Inited() const {
    return capacity_ > 0;
  }
  size_t ContainerSize() const {
    return container_.size();
  }
 private:
  size_t chunks_ = 0;
  size_t elements_ = 0;
  size_t load_factor_ = 0;
  size_t capacity_ = 0;
  std::vector<uint64_t> container_;

  bool IsPrime(size_t n) const {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n%2 == 0) return false;
    for (size_t i=3;i*i<=n;i+=2) {
      if (n%i == 0) return false;
    }
    return true;
  }

  size_t NextPrime(size_t n) const {
    while (!IsPrime(n)) {
      n++;
    }
    return n;
  }

  size_t Hash(const uint64_t* const data, size_t mod) const {
    uint64_t h = 0;
    for (size_t i=0;i<chunks_;i++){
      h = h*65599 + data[i];
    }
    h %= (uint64_t)mod;
    return h;
  }

  bool Equal(const uint64_t* const data1, const uint64_t* const data2) const {
    for (size_t i=0;i<chunks_;i++){
      if (data1[i] != data2[i]) return false;
    }
    return true;
  }

  bool Zero(const uint64_t* const data) const {
    for (size_t i=0;i<chunks_;i++){
      if (data[i]) return false;
    }
    return true;
  }

  const uint64_t* IndToPtr(size_t ind, const std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks_);
  }

  uint64_t* IndToPtr(size_t ind, std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks_);
  }

  void Copy(const uint64_t* const from, uint64_t* const to) {
    for (size_t i=0;i<chunks_;i++){
      to[i] = from[i];
    }
  }

  void Resize() {
    size_t new_capacity = NextPrime(capacity_ * 2);
    std::vector<uint64_t> new_container(new_capacity * chunks_);
    for (size_t i = 0; i < capacity_; i++) {
      if (Zero(IndToPtr(i, container_))) continue;
      size_t ind = Hash(IndToPtr(i, container_), new_capacity);
      while (!Zero(IndToPtr(ind, new_container))) {
        ind++;
        if (ind == new_capacity) {
          ind = 0;
        }
      }
      Copy(IndToPtr(i, container_), IndToPtr(ind, new_container));
    }
    container_ = std::move(new_container);
    capacity_ = new_capacity;
  }
};
} // namespace triangulator