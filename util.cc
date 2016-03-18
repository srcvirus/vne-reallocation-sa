#include "util.h"

#include <algorithm>
#include <map>
#include <stdarg.h>
#include <stdio.h>
#include <string>

void PrintDebugMessage(const char* location, const char* fmt_string, ...) {
  va_list args;
  va_start(args, fmt_string);
  std::string str = location;
  str += fmt_string;
  vprintf(str.c_str(), args);
  fflush(stdout);
  va_end(args);
}

template <class T>
double GetMean(const std::vector<T>& data) {
  T sum = T(0);
  const size_t kNumElements = data.size();
  for (int i = 0; i < data.size(); ++i) {
    sum += data[i];
  }
  return sum / static_cast<T>(kNumElements);
}

template <class T>
T GetNthPercentile(const std::vector<T>& data, int n) {
  std::vector<T> temp_data_buffer = data;
  sort(temp_data_buffer.begin(), temp_data_buffer.end());
  const size_t kNumElements = data.size();
  int rank = n * kNumElements;
  if (rank % 100) {
    rank = (rank / 100) + 1;
  } else
    rank /= 100;
  --rank;
  return temp_data_buffer[rank];
}

template <class T>
std::vector<std::pair<T, double> > GetCDF(const std::vector<T>& data) {
  int precision = 1;
  std::vector<T> temp_data_buffer = data;
  if (typeid(temp_data_buffer[0]) == typeid(double)||
      typeid(temp_data_buffer[0]) == typeid(float)) {
    precision = 1000;
  }
  std::map<int, int> cdf;
  for (int i = 0; i < temp_data_buffer.size(); ++i) {
    int bucket_index = temp_data_buffer[i] * precision;
    if (cdf[bucket_index])
      cdf[bucket_index]++;
    else
      cdf[bucket_index] = 1;
  }
  std::map<int, int>::iterator prev = cdf.begin(), current = cdf.begin();
  current++;
  for (; current != cdf.end(); current++, prev++) {
    current->second += prev->second;
  }
  int total = temp_data_buffer.size();
  std::vector<std::pair<T, double> > ret;
  for (current = cdf.begin(); current != cdf.end(); ++current) {
    T first = static_cast<T>(current->first) / static_cast<T>(precision);
    double second =
        static_cast<double>(current->second) / static_cast<double>(total);
    ret.push_back(std::make_pair(first, second));
  }
  return ret;
}

unique_ptr<ReverseEmbedding> GetInverseEmbedding(
    const Graph* phys_topology, const ptr_vector<VNEmbedding>& embeddings,
    int num_vns) {
  unique_ptr<ReverseEmbedding> r_embedding(new ReverseEmbedding());
  r_embedding->rnode_map.resize(phys_topology->node_count());
  for (int i = 0; i < num_vns; ++i) {
    // Generate physical_node --> set_of_vnodes mapping.
    for (int vnode = 0; vnode < embeddings[i].node_map.size(); ++vnode) {
      // DEBUG("i = %d, vnode = %d, embeddings[i].node_map[vnode] = %d\n", i,
      //      vnode, embeddings[i].node_map[vnode]);
      r_embedding->rnode_map[embeddings[i].node_map[vnode]]
          .insert(vnode_t(i, vnode));
    }

    // Generate physical_link --> set of virtual links mapping.
    edge_map_t::const_iterator emap_it;
    for (emap_it = embeddings[i].edge_map.begin();
         emap_it != embeddings[i].edge_map.end(); ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const std::vector<edge_t>& plinks = emap_it->second;
      std::vector<edge_t>::const_iterator plink_it;
      for (plink_it = plinks.begin(); plink_it != plinks.end(); ++plink_it) {
        if (r_embedding->redge_map.find(*plink_it) ==
            r_embedding->redge_map.end()) {
          r_embedding->redge_map.insert(
              std::make_pair(*plink_it, vedge_set_t()));
        }
        r_embedding->redge_map[*plink_it].insert(vn_edge_t(i, vlink));
      }
    }
  }
  return boost::move(r_embedding);
}

/*
long CostFunction(const Graph* phys_topology, int num_vns,
                  const ptr_vector<Graph>& virt_topologies,
                  const ptr_vector<VNEmbedding>& embeddings,
                  const ReverseEmbedding* r_embedding,
                  const VNRParameters& vnr_params) {
  double bw_cost = 0.0;
  double bl_cost = 0.0;

  // Bandwidth consumption cost.
  for (int i = 0; i < num_vns; ++i) {
    edge_map_t::const_iterator emap_it;
    for (emap_it = embeddings[i].edge_map.begin();
         emap_it != embeddings[i].edge_map.end(); ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      path_t::const_iterator plink_it;
      for (plink_it = plinks.begin(); plink_it != plinks.end(); ++plink_it) {
        bw_cost += static_cast<double>(
            phys_topology->GetEdgeCost(plink_it->first, plink_it->second) *
            virt_topologies[i].GetEdgeBandwidth(vlink.first, vlink.second));
      }
    }
  }
  bw_cost *= vnr_params.alpha;

  // Bottleneck Link cost.
  reverse_edge_map_t::const_iterator remap_it;
  int n_bottlenecks = 0;
  for (remap_it = r_embedding->redge_map.begin();
       remap_it != r_embedding->redge_map.end(); ++remap_it) {
    vedge_set_t::const_iterator vlink_it;
    const edge_t& plink = remap_it->first;
    double util = 0.0;
    for (vlink_it = remap_it->second.begin();
         vlink_it != remap_it->second.end(); ++vlink_it) {
      int vn_id = vlink_it->first;
      const edge_t& vlink = vlink_it->second;
      util += static_cast<double>(
          virt_topologies[vn_id].GetEdgeBandwidth(vlink.first, vlink.second));
    }
    util /= static_cast<double>(
        phys_topology->GetEdgeBandwidth(plink.first, plink.second));
    if (util >= vnr_params.util_threshold) ++n_bottlenecks;
  }
  bl_cost = vnr_params.beta * static_cast<double>(n_bottlenecks);

  double cost = bw_cost + bl_cost;
  return cost;
}
*/

unsigned long GetBandwidthCost(
    const Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings) {
  unsigned long cost = 0;
  for (int i = 0; i < virt_topologies.size(); ++i) {
    const VNEmbedding& embedding = vn_embeddings[i];
    const Graph& virt_topology = virt_topologies[i];
    edge_map_t::const_iterator emap_it;
    for (emap_it = embedding.edge_map.begin();
         emap_it != embedding.edge_map.end(); ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      path_t::const_iterator plink_it;
      for (plink_it = plinks.begin(); plink_it != plinks.end(); ++plink_it) {
        const edge_t& e = *plink_it;
        cost += phys_topology->GetEdgeCost(e.first, e.second) *
                virt_topology.GetEdgeBandwidth(vlink.first, vlink.second);
      }
    }
  }
  return cost;
}

int GetNumBottleneckLinks(const Graph* phys_topology,
                          const boost::ptr_vector<Graph>& virt_topologies,
                          const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                          const VNRParameters* vnr_param) {
  int num_bottlenecks = 0;
  matrix_t<double> util_matrix(phys_topology->node_count(),
                               phys_topology->node_count(), 0.0);
  for (int i = 0; i < virt_topologies.size(); ++i) {
    const Graph& virt_topology = virt_topologies[i];
    const VNEmbedding& embedding = vn_embeddings[i];
    edge_map_t::const_iterator emap_it;
    for (emap_it = embedding.edge_map.begin();
         emap_it != embedding.edge_map.end(); ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      path_t::const_iterator path_it;
      for (path_it = plinks.begin(); path_it != plinks.end(); ++path_it) {
        const edge_t& e = *path_it;
        long b_mn = virt_topology.GetEdgeBandwidth(vlink.first, vlink.second);
        util_matrix.matrix[e.first][e.second] += b_mn;
        util_matrix.matrix[e.second][e.first] += b_mn;
      }
    }
  }
  for (int u = 0; u < phys_topology->node_count(); ++u) {
    std::vector<edge_endpoint>::const_iterator end_point_it;
    const std::vector<edge_endpoint>& u_neighbors =
        phys_topology->adj_list()->at(u);
    for (end_point_it = u_neighbors.begin(); end_point_it != u_neighbors.end();
         ++end_point_it) {
      const edge_endpoint& end_point = *end_point_it;
      int v = end_point.node_id;
      long b_uv = phys_topology->GetEdgeBandwidth(u, v);

      util_matrix.matrix[u][v] /= static_cast<double>(b_uv);
      if (u < v && util_matrix.matrix[u][v] > vnr_param->util_threshold) {
        ++num_bottlenecks;
      }
    }
  }
  return num_bottlenecks;
}

double GetMaxPLinkUtilization(
    const Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings) {
  double max_util = 0.0;
  matrix_t<double> util_matrix(phys_topology->node_count(),
                               phys_topology->node_count(), 0.0);
  for (int i = 0; i < virt_topologies.size(); ++i) {
    const Graph& virt_topology = virt_topologies[i];
    const VNEmbedding& embedding = vn_embeddings[i];
    edge_map_t::const_iterator emap_it;
    for (emap_it = embedding.edge_map.begin();
         emap_it != embedding.edge_map.end(); ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      path_t::const_iterator path_it;
      for (path_it = plinks.begin(); path_it != plinks.end(); ++path_it) {
        const edge_t& e = *path_it;
        long b_mn = virt_topology.GetEdgeBandwidth(vlink.first, vlink.second);
        util_matrix.matrix[e.first][e.second] += b_mn;
        util_matrix.matrix[e.second][e.first] += b_mn;
      }
    }
  }
  for (int u = 0; u < phys_topology->node_count(); ++u) {
    std::vector<edge_endpoint>::const_iterator end_point_it;
    const std::vector<edge_endpoint>& u_neighbors =
        phys_topology->adj_list()->at(u);
    for (end_point_it = u_neighbors.begin(); end_point_it != u_neighbors.end();
         ++end_point_it) {
      const edge_endpoint& end_point = *end_point_it;
      int v = end_point.node_id;
      long b_uv = phys_topology->GetEdgeBandwidth(u, v);
      util_matrix.matrix[u][v] /= static_cast<double>(b_uv);
      if (u < v && util_matrix.matrix[u][v] > max_util) {
        max_util = util_matrix.matrix[u][v];
      }
    }
  }
  return max_util;
}

double CostFunction(const Graph* phys_topology,
                    const boost::ptr_vector<Graph>& virt_topologies,
                    const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                    const VNRParameters* vnr_param) {
  long bw_cost =
      GetBandwidthCost(phys_topology, virt_topologies, vn_embeddings);
  int num_bottlenecks = GetNumBottleneckLinks(phys_topology, virt_topologies,
                                              vn_embeddings, vnr_param);
  double cost = vnr_param->alpha * static_cast<double>(bw_cost) +
                vnr_param->beta * static_cast<double>(num_bottlenecks);
  return cost;
}

void ComputePhysicalNetworkCapacity(
    Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings) {
  for (int i = 0; i < virt_topologies.size(); ++i) {
    const Graph& vn = virt_topologies[i];
    const VNEmbedding& vne = vn_embeddings[i];
    edge_map_t::const_iterator emap_it;
    for (emap_it = vne.edge_map.begin(); emap_it != vne.edge_map.end();
         ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      int m = vlink.first, n = vlink.second;
      long b_mn = vn.GetEdgeBandwidth(m, n);
      path_t::const_iterator plink_it;
      for (plink_it = plinks.begin(); plink_it != plinks.end(); ++plink_it) {
        int u = plink_it->first, v = plink_it->second;
        long cur_bw = phys_topology->GetEdgeBandwidth(u, v);
        phys_topology->SetEdgeBandwidth(u, v, cur_bw + b_mn);
        phys_topology->SetEdgeBandwidth(v, u, cur_bw + b_mn);
      }
    }
  }
}

void WriteSolutionToFile(const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                         const std::string& vnr_directory) {
  for (int vn_index = 0; vn_index < vn_embeddings.size(); ++vn_index) {
    std::string emap_file = vnr_directory + "/vn" +
                            boost::lexical_cast<std::string>(vn_index) +
                            ".edge_remap.txt";
    std::string nmap_file = vnr_directory + "/vn" +
                            boost::lexical_cast<std::string>(vn_index) +
                            ".node_remap.txt";
    FILE* nmap_f = fopen(nmap_file.c_str(), "w");
    const VNEmbedding& vne = vn_embeddings[vn_index];
    for (int i = 0; i < vne.node_map.size(); ++i) {
      fprintf(nmap_f, "VN %d: Virtual node %d --> physical node %d\n", vn_index,
              i, vne.node_map[i]);
    }
    fclose(nmap_f);

    FILE* emap_f = fopen(emap_file.c_str(), "w");
    edge_map_t::const_iterator emap_it;
    for (emap_it = vne.edge_map.begin(); emap_it != vne.edge_map.end();
         ++emap_it) {
      const edge_t& vlink = emap_it->first;
      const path_t& plinks = emap_it->second;
      path_t::const_iterator pit;
      for (pit = plinks.begin(); pit != plinks.end(); ++pit) {
        fprintf(emap_f,
                "VN %d: Virtual edge (%d, %d) --> physical edge (%d, %d)\n",
                vn_index, vlink.first, vlink.second, pit->first, pit->second);
      }
    }
    fclose(emap_f);
  }
}
