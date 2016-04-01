#ifndef _GRAPH_UTIL_H_
#define _GRAPH_UTIL_H_

#include "datastructure.h"

// If delta is -ve then residual bandwidth is reduced. If delta is +ve then
// residual bandwidth is increased.
void UpdateResidualBandwidth(matrix_t<long>& res_bw_matrix, const path_t& path,
                             long delta) {
  path_t::const_iterator path_it;
  for (path_it = path.begin(); path_it != path.end(); ++path_it) {
    int u = path_it->first, v = path_it->second;
    res_bw_matrix.matrix[u][v] += delta;
    res_bw_matrix.matrix[v][u] += delta;
  }
}

void UpdateUtilMatrix(const Graph* graph, matrix_t<double>& util_matrix,
                      const matrix_t<long>& res_bw_matrix, const path_t& path) {
  path_t::const_iterator path_it;
  for (path_it = path.begin(); path_it != path.end(); ++path_it) {
    int u = path_it->first, v = path_it->second;
    long b_uv = graph->GetEdgeBandwidth(u, v);
    long res_uv = res_bw_matrix.matrix[u][v];
    util_matrix.matrix[u][v] = util_matrix.matrix[v][u] =
        1.0 - static_cast<double>(res_uv) / static_cast<double>(b_uv);
  }
}

int GetNumBottleneckLinksOnPath(const matrix_t<double>& util_matrix,
                                const path_t& path, double util_threshold) {
  path_t::const_iterator path_it;
  int n_bottlenecks = 0;
  for (path_it = path.begin(); path_it != path.end(); ++path_it) {
    if (util_matrix.matrix[path_it->first][path_it->second] > util_threshold)
      ++n_bottlenecks;
  }
  return n_bottlenecks;
}

struct dijkstra_node {
  int node_id;
  long cost;
  dijkstra_node(int id, long c) : node_id(id), cost(c) {}
  bool operator<(const dijkstra_node& dnode) const { return cost < dnode.cost; }
};

// Find the least cost path between src and dest having at least bw bandwidth
// and avoiding edges in to_avoid. Cost of a link is computed as the product of
// per unit bandwidth cost and the used bandwidth on the link. Hence, the least
// cost path will contain as less bottleneck links as possible.
unique_ptr<path_t> ConstrainedVLinkEmbed(const Graph* phys_topology,
                                         const matrix_t<double>& util_matrix,
                                         const matrix_t<long>& res_bw_matrix,
                                         int src, int dest, long bw,
                                         const std::set<edge_t>& to_avoid) {
  std::vector<long> d(phys_topology->node_count(), INF);
  std::vector<int> pre(phys_topology->node_count(), NIL);
  std::vector<bool> finished(phys_topology->node_count(), false);
  std::priority_queue<dijkstra_node> pq;

  d[src] = 0;
  pq.push(dijkstra_node(src, 0.0));

  while (!pq.empty()) {
    dijkstra_node dnode = pq.top();
    pq.pop();
    int u = dnode.node_id;
    if (finished[u]) {
      continue;
    }
    finished[u] = true;
    for (int i = 0; i < phys_topology->adj_list()->at(u).size(); ++i) {
      const edge_endpoint& end_point = phys_topology->adj_list()->at(u)[i];
      int v = end_point.node_id;
      long link_res_bw = res_bw_matrix.matrix[u][v];
      if ((to_avoid.find(edge_t(u, v)) != to_avoid.end()) ||
          (to_avoid.find(edge_t(v, u)) != to_avoid.end()))
        continue;
      if (link_res_bw < bw) continue;
      long cost =
          end_point.cost *
          (bw + (static_cast<long>(util_matrix.matrix[u][v] *
                                   static_cast<double>(end_point.bandwidth))));
      if (d[v] > d[u] + cost) {
        d[v] = d[u] + cost;
        pre[v] = u;
        pq.push(dijkstra_node(v, d[v]));
      }
    }
  }
  unique_ptr<path_t> ret_path(new path_t());
  int current = dest;
  while (pre[current] != NIL) {
    edge_t new_edge(pre[current], current);
    if (pre[current] > current) {
      std::swap(new_edge.first, new_edge.second);
    }
    ret_path->push_back(new_edge);
    current = pre[current];
  }
  std::reverse(ret_path->begin(), ret_path->end());
  return boost::move(ret_path);
}

#endif  // _GRAPH_UTIL_H_
