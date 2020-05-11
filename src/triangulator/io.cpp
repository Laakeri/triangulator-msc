#include "io.hpp"

#include <vector>
#include <string>
#include <istream>
#include <sstream>
#include <algorithm>
#include <cassert>

#include "graph.hpp"
#include "hypergraph.hpp"
#include "utils.hpp"
#include "phyl_mat.hpp"
#include "matrix.hpp"

namespace triangulator {

int NumTokens(std::string s) {
  std::stringstream ss;
  ss<<s;
  int num = 0;
  while (ss>>s) {
    num++;
  }
  return num;
}

std::vector<std::string> GetTokens(std::string s) {
  std::stringstream ss;
  ss<<s;
  std::vector<std::string> tokens;
  while (ss>>s) {
    tokens.push_back(s);
  }
  return tokens;
}

Graph Io::ReadGraph(std::istream& in) {
  std::vector<std::pair<std::string, std::string> > edges;
  std::string input_line;
  bool format_detected = false;
  dimacs_ = false;
  bool pace_ = false;
  int line_num = 0;
  while (std::getline(in, input_line)) {
    line_num++;
    assert(input_line.size() > 0);
    if (input_line.size() >= 2 && input_line.substr(0, 2) == "p " && NumTokens(input_line) == 3) {
      if (format_detected) {
        utils::Warning("Detected another begin of a file format in line ", line_num, " lines before that are ignored");
      }
      Log::Write(10, "Dimacs graph format detected");
      edges.clear();
      dimacs_ = true;
      format_detected = true;
    } else if (input_line.size() >= 2 && input_line.substr(0, 5) == "p tw " && NumTokens(input_line) == 4) {
      if (format_detected) {
        utils::Warning("Detected another begin of a file format in line ", line_num, " lines before that are ignored");
      }
      Log::Write(10, "PACE graph format detected");
      edges.clear();
      pace_ = true;
      format_detected = true;
    } else if (pace_ && input_line.size() >= 2 && NumTokens(input_line) == 2) {
      auto tokens = GetTokens(input_line);
      edges.push_back({tokens[0], tokens[1]});
    } else if (dimacs_ && input_line.size() >= 2 && input_line.substr(0, 2) == "e " && NumTokens(input_line) == 3) {
      auto tokens = GetTokens(input_line);
      assert(tokens[0] == "e");
      edges.push_back({tokens[1], tokens[2]});
    } else if (input_line.size() >= 2 && NumTokens(input_line) == 2) {
      auto tokens = GetTokens(input_line);
      edges.push_back({tokens[0], tokens[1]});
    }
  }
  vertex_map_.Init(edges);
  Graph graph(vertex_map_.Size());
  for (auto edge : edges) {
    graph.AddEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
  return graph;
}

HyperGraph Io::ReadHyperGraph(std::istream& in) {
  std::map<std::string, int> vertices;
  std::vector<std::vector<int>> edges;
  std::string tmp;
  while (getline(in, tmp)) {
    std::string line;
    int open = 0;
    int close = 0;
    for (char& c : tmp) {
      if (c == ',') {
        line += " , ";
      } else if (c == '(') {
        line += " ( ";
        open++;
      }
      else if (c == ')') {
        line += " ) ";
        close++;
      } else {
        line += c;
      }
    }
    assert(open == close);
    std::stringstream ss;
    ss<<line;
    for (int i=0;i<open;i++){
      ss>>tmp;
      ss>>tmp;
      assert(tmp == "(");
      std::vector<int> edge;
      while (ss>>tmp) {
        if (!vertices.count(tmp)) {
          int ti = vertices.size();
          vertices[tmp] = ti;
        }
        edge.push_back(vertices[tmp]);
        ss>>tmp;
        if (tmp == ")") break;
        else assert(tmp == ",");
      }
      utils::SortAndDedup(edge);
      edges.push_back(edge);
      assert(tmp == ")");
      ss>>tmp;
      assert(tmp == "," || tmp == ".");
    }
  }
  HyperGraph hg(edges);
  return hg;
}

PhylMat Io::ReadPhylMat(std::istream& in) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  int n, m;
  in>>n>>m;
  assert(n > 0 && m > 0);
  Matrix<int> mat(n, m);
  for (int i = 0; i < n; i++) {
    for (int ii = 0; ii < m; ii++) {
      std::string t;
      in>>t;
      if (t == "?") {
        mat[i][ii] = -1;
      } else {
        try {
          mat[i][ii] = stoi(t);
        } catch (const std::exception& e) {
          assert(0);
        }
      }
    }
  }
  PhylMat pm(mat);
  return pm;
}

std::string Io::MapBack(int v) {
  return vertex_map_.Kth(v);
}

std::pair<std::vector<int64_t>, Graph> Io::ReadBayesNet(std::istream& in) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  std::string tmp;
  bool net = false;
  std::map<std::string, int> nodes;
  std::vector<int64_t> domain;
  std::vector<std::vector<int>> inc;
  while (in>>tmp){
    if (tmp == "net") {
      assert(!net);
      net = true;
      in>>tmp;
      assert(tmp == "{");
      in>>tmp;
      assert(tmp == "}");
      continue;
    }
    if (!net) continue;
    if (tmp == "node") {
      in>>tmp;
      assert(!nodes.count(tmp));
      int ni = nodes.size();
      nodes[tmp] = ni;
      in>>tmp;
      assert(tmp == "{");
      in>>tmp;
      assert(tmp == "states");
      in>>tmp;
      assert(tmp == "=");
      in>>tmp;
      assert(tmp == "(");
      int dm = 0;
      while (in>>tmp){
        if (tmp == ");") break;
        assert(tmp[0] == '\"' && tmp.back() == '\"');
        dm++;
      }
      assert(dm > 0);
      domain.push_back(dm);
      inc.push_back({});
      assert(nodes.size() == domain.size() && nodes.size() == inc.size());
      in>>tmp;
      assert(tmp == "}");
      continue;
    }
    if (tmp == "potential") {
      in>>tmp;
      assert(tmp == "(");
      in>>tmp;
      assert(nodes.count(tmp));
      int tv = nodes[tmp];
      in>>tmp;
      if (tmp == "|") {
        while (in>>tmp) {
          if (tmp == ")") break;
          assert(nodes.count(tmp));
          inc[tv].push_back(nodes[tmp]);
        }
      } else {
        assert(tmp == ")");
      }
      in>>tmp;
      assert(tmp == "{");
      in>>tmp;
      assert(tmp == "data");
      while (in>>tmp){
        if (tmp == "}") break;
      }
      continue;
    }
    assert(0);
  }
  assert(nodes.size() == domain.size() && nodes.size() == inc.size());
  Graph moral(nodes.size());
  for (int i=0;i<moral.n();i++){
    for (int a : inc[i]) {
      moral.AddEdge(a, i);
      for (int b : inc[i]) {
        if (a!=b) moral.AddEdge(a, b);
      }
    }
  }
  assert(moral.n() > 0);
  return {domain, moral};
}

} // namespace triangulator
