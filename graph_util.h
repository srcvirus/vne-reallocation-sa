#ifndef _GRAPH_UTIL_H_
#define _GRAPH_UTIL_H_

#include "datastructure.h"

// If erase is true then remove channel (path.first) from res_channel_matrix,
// otherwise, add channel (path.first) to res_channel_matrix.
void UpdateResidualChannel(matrix_t<std::vector<int> >& res_channel_matrix,
                           const dwdm_path_t& path, bool erase) {
  path_t::const_iterator path_it;
  int ch_id = path.first;
  for (path_it = path.second.begin(); path_it != path.second.end(); ++path_it) {
    int u = path_it->first, v = path_it->second;
    std::vector<int>::iterator it;
    if (erase) {
      it = std::find(res_channel_matrix.matrix[u][v].begin(),
                     res_channel_matrix.matrix[u][v].end(), ch_id);
      res_channel_matrix.matrix[u][v].erase(it);
      it = std::find(res_channel_matrix.matrix[v][u].begin(),
                     res_channel_matrix.matrix[v][u].end(), ch_id);
      res_channel_matrix.matrix[v][u].erase(it);
    } else {
      res_channel_matrix.matrix[u][v].insert(
          std::upper_bound(res_channel_matrix.matrix[u][v].begin(),
                           res_channel_matrix.matrix[u][v].end(), ch_id),
          ch_id);
      res_channel_matrix.matrix[v][u].insert(
          std::upper_bound(res_channel_matrix.matrix[v][u].begin(),
                           res_channel_matrix.matrix[v][u].end(), ch_id),
          ch_id);
    }
  }
}

void UpdateUtilMatrix(const Graph* graph, matrix_t<double>& util_matrix,
                      const matrix_t<std::vector<int> >& res_channel_matrix,
                      const path_t& path) {
  path_t::const_iterator path_it;
  for (path_it = path.begin(); path_it != path.end(); ++path_it) {
    int u = path_it->first, v = path_it->second;
    long ch_uv = graph->GetEdgeNumChannels(u, v);
    long res_ch = res_channel_matrix.matrix[u][v].size();
    util_matrix.matrix[u][v] = util_matrix.matrix[v][u] =
        1.0 - static_cast<double>(res_ch) / static_cast<double>(ch_uv);
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

// Find the least cost path between src and dest and allocate the first channel
// that is available. Avoid the edges listed in to_avoid. Cost of a link is
// computed as the product of per channel cost and the used channels on the
// link. Hence, the least cost path will contain as less bottleneck links as
// possible.
unique_ptr<dwdm_path_t> ConstrainedDWDMVLinkEmbed(
    const Graph* phys_topology, const matrix_t<double>& util_matrix,
    const matrix_t<std::vector<int> >& res_channel_matrix, int src, int dest,
    int required_channels, int max_channels, const std::set<edge_t>& to_avoid) {
  int ch_lo = 0, ch_hi = max_channels - 1;

  unique_ptr<dwdm_path_t> ret_path(new dwdm_path_t(NIL, path_t()));
  for (int ch = ch_lo; ch <= ch_hi; ++ch) {
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
        bool has_valid_channel =
            (std::find(res_channel_matrix.matrix[u][v].begin(),
                       res_channel_matrix.matrix[u][v].end(),
                       ch) != res_channel_matrix.matrix[u][v].end());
        if (!has_valid_channel) continue;
        int link_res_channels = res_channel_matrix.matrix[u][v].size();
        if ((to_avoid.find(edge_t(u, v)) != to_avoid.end()) ||
            (to_avoid.find(edge_t(v, u)) != to_avoid.end()))
          continue;
        if (link_res_channels < required_channels) continue;
        long cost =
            end_point.cost *
            (required_channels +
             (static_cast<long>(util_matrix.matrix[u][v] *
                                static_cast<double>(end_point.channels))));
        if (d[v] > d[u] + cost) {
          d[v] = d[u] + cost;
          pre[v] = u;
          pq.push(dijkstra_node(v, d[v]));
        }
      }
    }
    if (d[dest] != INF) {
      ret_path->second.clear();
      ret_path->first = ch;
      int current = dest;
      while (pre[current] != NIL) {
        edge_t new_edge(pre[current], current);
        if (pre[current] > current) {
          std::swap(new_edge.first, new_edge.second);
        }
        ret_path->second.push_back(new_edge);
        current = pre[current];
      }
      std::reverse(ret_path->second.begin(), ret_path->second.end());
    }
  }
  return boost::move(ret_path);
}

#endif  // _GRAPH_UTIL_H_
