#pragma once

#include <vector>
#include <ostream>
#include <set>

#include "utils.hpp"
#include "staticset.hpp"
#include "matrix.hpp"
#include "bitset.hpp"

namespace triangulator {

typedef std::pair<int, int> Edge;

class Graph {
public:
  explicit Graph(int n);
  explicit Graph(std::vector<Edge> edges);
  void AddEdge(int v, int u);
  void AddEdge(Edge e);
  void AddEdges(const std::vector<Edge>& edges);

  void RemoveEdge(int v, int u);
  void RemoveEdgesBetween(int v, const std::vector<int>& vs);
  
  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  std::vector<Edge> Edges() const;
  std::vector<int> Vertices() const;

  int Degeneracy() const;
  
  bool IsConnected() const;
  bool IsConnectedOrIsolated() const;
  
  Bitset Neighbors(const Bitset& vs) const;
  const std::vector<int>& Neighbors(int v) const;
  std::vector<std::vector<int> > Components(const std::vector<int>& separator) const;
  std::vector<std::vector<int> > NComponents(const std::vector<int>& separator) const;
  std::vector<Bitset> NComponents(const Bitset& bs) const;
  std::vector<int> Neighbors(const std::vector<int>& vs) const;
  std::vector<int> FindComponentAndMark(int v, std::vector<char>& block) const;
  std::vector<Edge> EdgesIn(const std::vector<int>& vs) const;
  
  bool IsClique(const std::vector<int>& clique) const;
  bool IsAlmostClique(const std::vector<int>& clq) const;
  bool IsClique(Bitset bs) const;

  std::vector<Edge> FillEdges(const std::vector<int>& clq) const;
  std::vector<Edge> FillEdges(const Graph& other) const;
  std::vector<Edge> FillEdges(Bitset bs) const;
  void FillBS(Bitset bs);
  int FillSize(Bitset bs) const;
  
  int MapBack(int v) const;
  std::vector<int> MapBack(std::vector<int> vs) const;
  Edge MapBack(Edge e) const;
  std::vector<Edge> MapBack(std::vector<Edge> es) const;
  std::pair<int, int> MapBack(int v, int u) const;
  int MapInto(int v) const;
  std::vector<int> MapInto(std::vector<int> vs) const;
  std::set<std::pair<int, int>> MapInto(std::set<std::pair<int, int>> vs) const;
  
  void InheritMap(const Graph& parent);

  std::vector<std::vector<int>> CompNeighs(const std::vector<int>& block) const;
  std::vector<int> CompNeigh(const std::vector<int>& block, int v) const;
  std::vector<Bitset> CompNeighsBit(const Bitset& block) const;

  Bitset AnotherComp(int x, const Bitset& minsep) const;

  bool IsMinsep(const std::vector<int>& separator) const;
  bool IsMinsep(const std::vector<int>& separator, int a, int b) const;
  bool IsMinsep(const Bitset& separator) const;
  bool HasNFullComponents(const Bitset& separator, int n) const;
  bool HasNFullComponents(const std::vector<int>& separator, int n) const;
  
  void Dfs2(int v, Bitset& sep, Bitset& vis, std::vector<int>& f) const;
  bool IsFull(int v, Bitset sep, Bitset vis) const;
  bool IsFull2(int v, Bitset sep, Bitset& vis) const;
  std::vector<Bitset> BitComps(Bitset vis) const;
  void Dfs22(int v, Bitset& sep, Bitset& vis, std::vector<int>& f, const Bitset& good) const;
  void Dfs2Bit(Bitset& vis, Bitset& ne) const;

  std::vector<Bitset> adj_mat2_;

  std::vector<int> Distances(const std::vector<int>& start) const;
  std::vector<std::vector<int>> DistanceMatrix() const;

  std::vector<Bitset> AllMinseps() const;

  int MaximalIS(const Bitset& vs) const;

  Graph(const Graph& rhs) = default;
  Graph& operator=(const Graph& rhs) = default;
  
private:
  int n_, m_;
  StaticSet<int> vertex_map_;
  std::vector<std::vector<int> > adj_list_;
  Matrix<char> adj_mat_;
  void Dfs(int v, std::vector<char>& blocked, std::vector<int>& component) const;
};
} // namespace triangulator
