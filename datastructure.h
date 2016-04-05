#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <boost/heap/fibonacci_heap.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/move/unique_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <list>
#include <set>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <map>
#include <memory>
#include <stdlib.h>

#define INF 99999999
#define MAXN 1000
#define NIL -1

using boost::movelib::unique_ptr;
using boost::ptr_vector;
using boost::heap::fibonacci_heap;

// Type definition for convenience.
typedef std::pair<int, int> edge_t;
typedef std::pair<int, edge_t> dwdm_edge_t;
typedef std::pair<int, int> vnode_t;
typedef std::pair<int, edge_t> vn_edge_t;
typedef std::vector<edge_t> path_t;
typedef std::pair<int, path_t> dwdm_path_t;
typedef std::set<edge_t> edge_set_t;
typedef std::set<dwdm_edge_t> dwdm_edge_set_t;
typedef std::set<vn_edge_t> vedge_set_t;
typedef std::map<edge_t, path_t> edge_map_t;
typedef std::map<edge_t, dwdm_path_t> dwdm_edge_map_t;
typedef std::map<dwdm_edge_t, vedge_set_t> dwdm_redge_map_t;
typedef std::map<edge_t, vedge_set_t> reverse_edge_map_t;

inline edge_t ConstructEdge(int u, int v) {
  if (u > v) std::swap(u, v);
  return edge_t(u, v);
}

inline dwdm_edge_t ConstructDWDMEdge(int u, int v, int w) {
  if (u > v) std::swap(u, v);
  return dwdm_edge_t(w, edge_t(u, v));
}

template <typename T>
struct matrix_t {
  std::vector<std::vector<T> > matrix;
  matrix_t() {}
  matrix_t(int rows, int columns, T fill_value = T())
      : matrix(rows, std::vector<T>(columns, fill_value)) {}
};

struct bneck_edge_element_t {
  edge_t bneck_edge;
  double util;

  bneck_edge_element_t() {}
  bneck_edge_element_t(const edge_t& be, const double u = 0.0)
      : bneck_edge(be), util(u) {}
  bool operator<(const bneck_edge_element_t& x) const { return util < x.util; }
};

struct edge_plength_set_element_t {
  vn_edge_t vn_edge;
  int path_len;

  edge_plength_set_element_t() {}
  edge_plength_set_element_t(const vn_edge_t& ve, const int pl = 0)
      : vn_edge(ve), path_len(pl) {}

  bool operator<(const edge_plength_set_element_t& x) const {
    return path_len < x.path_len;
  }
};

struct node_bw_set_element_t {
  vnode_t vnode;
  long bw_usage;

  node_bw_set_element_t() {}
  node_bw_set_element_t(const vnode_t& vn, const long bwu = 0)
      : vnode(vn), bw_usage(bwu) {}
  bool operator<(const node_bw_set_element_t& x) const {
    return bw_usage < x.bw_usage;
  }
};

// An entry in an adjacent list. An entry contains the node_id of the endpoint.
// The entry contains #channels, delay and cost of the corresponding edge.
struct edge_endpoint {
  int node_id;
  long channels;
  int delay;
  int cost;
  edge_endpoint() : node_id(NIL), channels(0), delay(INF), cost(INF) {}
  edge_endpoint(int node_id, long ch, int delay, int cost)
      : node_id(node_id), channels(ch), delay(delay), cost(cost) {}
  std::string GetDebugString() const {
    return "ndoe_id = " + boost::lexical_cast<std::string>(node_id) +
           ", channels = " + boost::lexical_cast<std::string>(channels) +
           ", delay = " + boost::lexical_cast<std::string>(delay) +
           ", cost = " + boost::lexical_cast<std::string>(cost);
  }
};

// Graph data structure. Represents an undirected graph. Internally, the graph
// is represented as an adjacency list. As a speed hack, a function
// Graph::Matrixize() is provided. When invoked, this function converts the
// internal adjacency list to a adjacency matrix and then can lookup/update
// link bandwidth in constant time. Be advised to use Graph::Matrixize() only
// after the structure of the graph is stable, i.e., no link/node
// addition/deletion
// is required.
class Graph {
 public:
  Graph() {
    adj_list_ = unique_ptr<std::vector<std::vector<edge_endpoint> > >(
        new std::vector<std::vector<edge_endpoint> >);
    node_count_ = edge_count_ = max_channels_ = 0;
    is_matrixized_ = false;
  }

  Graph(const Graph& g) {
    this->node_count_ = g.node_count_;
    this->edge_count_ = g.edge_count_;
    this->max_channels_ = g.max_channels_;
    this->adj_list_ = unique_ptr<std::vector<std::vector<edge_endpoint> > >(
        new std::vector<std::vector<edge_endpoint> >(*g.adj_list_.get()));
    this->is_matrixized_ = g.is_matrixized_;
    this->adj_matrix_ = g.adj_matrix_;
  }

  // Accessor methods.
  int node_count() const { return node_count_; }
  int edge_count() const { return edge_count_; }
  int max_channels() const { return max_channels_; }
  const std::vector<std::vector<edge_endpoint> >* adj_list() const {
    return static_cast<const std::vector<std::vector<edge_endpoint> >*>(
        adj_list_.get());
  }

  // Populate the internal adjacency matrix data structure with the contents of
  // adjacency list. This is added for some efficiency gain while reading link
  // meta-data, e.g., bandwidth, cost.
  void Matrixize() {
    adj_matrix_.resize(node_count_);
    for (int u = 0; u < node_count_; ++u) adj_matrix_[u].resize(node_count_);
    for (int u = 0; u < node_count_; ++u) {
      const std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
      std::vector<edge_endpoint>::const_iterator end_point_it;
      for (end_point_it = neighbors.begin(); end_point_it != neighbors.end();
           ++end_point_it) {
        int v = end_point_it->node_id;
        adj_matrix_[u][v] = *end_point_it;
      }
    }
    is_matrixized_ = true;
  }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::AddEdge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  void AddEdge(int u, int v, long channels, int delay, int cost) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    const std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
    std::vector<edge_endpoint>::const_iterator end_point_it;
    for (end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id == v) {
        return;
      }
    }
    adj_list_->at(u).push_back(edge_endpoint(v, channels, delay, cost));
    adj_list_->at(v).push_back(edge_endpoint(u, channels, delay, cost));
    ++edge_count_;
    node_count_ = adj_list_->size();
    if (max_channels_ < channels) max_channels_ = channels;
  }

  int GetEdgeCost(int u, int v) const {
    if (is_matrixized_) return adj_matrix_[u][v].cost;
    const std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
    std::vector<edge_endpoint>::const_iterator end_point_it;
    for (end_point_it = neighbors.begin(); end_point_it < neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id == v) return end_point_it->cost;
    }
  }

  long GetEdgeNumChannels(int u, int v) const {
    if (is_matrixized_) return adj_matrix_[u][v].channels;
    const std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
    std::vector<edge_endpoint>::const_iterator end_point_it;
    for (end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id == v) return end_point_it->channels;
    }
  }

  void SetEdgeNumChannels(int u, int v, long channels) {
    if (is_matrixized_) adj_matrix_[u][v].channels = channels;
    std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
    std::vector<edge_endpoint>::iterator end_point_it;
    for (end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id == v) {
        end_point_it->channels = channels;
        break;
      }
    }
  }

  long GetTotalNodeChannels(int u) const {
    const std::vector<edge_endpoint>& neighbors = adj_list_->at(u);
    std::vector<edge_endpoint>::const_iterator end_point_it;
    long total_ch = 0;
    for (end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      total_ch += end_point_it->channels;
    }
    return total_ch;
  }

  int GetNodeDegree(int u) const { return adj_list_->at(u).size(); }

  std::string GetDebugString() const {
    std::string ret_string =
        "node_count = " + boost::lexical_cast<std::string>(node_count_);
    ret_string += ", edge_count = " +
                  boost::lexical_cast<std::string>(edge_count_) + "\n";
    for (int i = 0; i < node_count_; ++i) {
      const std::vector<edge_endpoint>& neighbors = adj_list_->at(i);
      ret_string += boost::lexical_cast<std::string>(i) + " --> ";
      for (int i = 0; i < neighbors.size(); ++i) {
        const edge_endpoint& neighbor = neighbors[i];
        ret_string += " (" + neighbor.GetDebugString() + ")";
      }
      ret_string += "\n";
    }
    return ret_string;
  }

  virtual ~Graph() { adj_list_.reset(); }

 private:
  unique_ptr<std::vector<std::vector<edge_endpoint> > > adj_list_;
  std::vector<std::vector<edge_endpoint> > adj_matrix_;
  int node_count_, edge_count_;
  int max_channels_;
  bool is_matrixized_;
};

struct VNEmbedding {
  std::vector<int> node_map;
  dwdm_edge_map_t edge_map;
  VNEmbedding() {}
  VNEmbedding(const std::vector<int>& nmap, const dwdm_edge_map_t& emap)
      : node_map(nmap), edge_map(emap) {}
  virtual ~VNEmbedding() {}
};

// Contains reverse embedding, i.e., physical --> virtual embedding.
// rnode_map contains a mapping from a physical node to the set of virtual
// nodes  mapped on that physical node. redge_map contains mapping from a
// physical
// link to the set of virtual links mapped on a path containing that physical
// link.
struct ReverseEmbedding {
  std::vector<std::set<vnode_t> > rnode_map;
  reverse_edge_map_t redge_map;
  ReverseEmbedding() {}
  ReverseEmbedding(const std::vector<std::set<vnode_t> >& rnmap,
                   const reverse_edge_map_t& remap)
      : rnode_map(rnmap), redge_map(remap) {}
  virtual ~ReverseEmbedding() {}
};

struct VNRParameters {
  // Weights of the cost components.
  double alpha, beta;
  // Threshold to determine bottlenec links.
  double util_threshold;
};

struct SASolution {
  // Immutable pointer to the topology.
  const Graph* physical_topology;

  // Immutable pointer vector the set of virtual topologies.
  const ptr_vector<Graph>& virt_topologies;

  // Immutable pointer vector to the set of location constraints for the virtual
  // topologies.
  const ptr_vector<std::vector<std::vector<int> > >& location_constraints;

  // Optimization parameters.
  const VNRParameters& vnr_params;

  // Number of VNs.
  int num_vns;

  // Cost of this solution.
  double cost;

  // Embedding of each virtual network from virt_topologies.
  ptr_vector<VNEmbedding> vn_embeddings;

  // Reverse embedding of the virtual networks, i.e., physical --> virtual
  // mapping.
  unique_ptr<ReverseEmbedding> r_embedding;

  // Physical network utilization matrix.
  matrix_t<double> util_matrix;

  // Physical network residual number of channels matrix.
  matrix_t<std::vector<int> > res_channel_matrix;

  // A heap of bottleneck physical links.
  fibonacci_heap<bneck_edge_element_t> bottleneck_edges;

  // A heap of virtual links sorted by their mapped physical path lengths.
  fibonacci_heap<edge_plength_set_element_t> vlinks_by_plength;

  // A heap of virtual nodes sorted by the cumulative number of channel
  // consumption
  // of their incident virtual links.
  fibonacci_heap<node_bw_set_element_t> node_bw_usage;

  // Bookkeeping datastructure to keep track of the heap element from each of
  // the aforementioned heaps.
  std::map<edge_t, fibonacci_heap<bneck_edge_element_t>::handle_type>
      bneck_heap_handlers;
  std::map<vn_edge_t, fibonacci_heap<edge_plength_set_element_t>::handle_type>
      vlp_heap_handlers;
  std::map<vnode_t, fibonacci_heap<node_bw_set_element_t>::handle_type>
      node_bw_heap_handlers;

  // Overloaded assignment opereator.
  SASolution& operator=(const SASolution& sa_sol) {
    util_matrix = sa_sol.util_matrix;
    res_channel_matrix = sa_sol.res_channel_matrix;
    r_embedding = unique_ptr<ReverseEmbedding>(new ReverseEmbedding(
        sa_sol.r_embedding->rnode_map, sa_sol.r_embedding->redge_map));
    cost = sa_sol.cost;
    vn_embeddings.clear();
    for (int i = 0; i < num_vns; ++i) {
      vn_embeddings.push_back(new VNEmbedding(
          sa_sol.vn_embeddings[i].node_map, sa_sol.vn_embeddings[i].edge_map));
    }
    bottleneck_edges = sa_sol.bottleneck_edges;
    // Handle the handlers here.
    for (fibonacci_heap<bneck_edge_element_t>::iterator be_it =
             bottleneck_edges.begin();
         be_it != bottleneck_edges.end(); ++be_it) {
      bneck_heap_handlers[be_it->bneck_edge] =
          fibonacci_heap<bneck_edge_element_t>::s_handle_from_iterator(be_it);
    }

    vlinks_by_plength = sa_sol.vlinks_by_plength;
    // Handle the handlers here.
    for (fibonacci_heap<edge_plength_set_element_t>::iterator ep_it =
             vlinks_by_plength.begin();
         ep_it != vlinks_by_plength.end(); ++ep_it) {
      vlp_heap_handlers[ep_it->vn_edge] =
          fibonacci_heap<edge_plength_set_element_t>::s_handle_from_iterator(
              ep_it);
    }

    node_bw_usage = sa_sol.node_bw_usage;
    // Handle the handlers here.
    for (fibonacci_heap<node_bw_set_element_t>::iterator nb_it =
             node_bw_usage.begin();
         nb_it != node_bw_usage.end(); ++nb_it) {
      node_bw_heap_handlers[nb_it->vnode] =
          fibonacci_heap<node_bw_set_element_t>::s_handle_from_iterator(nb_it);
    }
    return *this;
  }

  // Overloaded copy constructor.
  SASolution(const SASolution& sa_sol)
      : physical_topology(sa_sol.physical_topology),
        virt_topologies(sa_sol.virt_topologies),
        location_constraints(sa_sol.location_constraints),
        vnr_params(sa_sol.vnr_params),
        num_vns(sa_sol.num_vns),
        util_matrix(sa_sol.util_matrix),
        res_channel_matrix(sa_sol.res_channel_matrix),
        cost(sa_sol.cost) {
    // Populate reverse embedding.
    r_embedding = unique_ptr<ReverseEmbedding>(new ReverseEmbedding(
        sa_sol.r_embedding->rnode_map, sa_sol.r_embedding->redge_map));

    for (int i = 0; i < num_vns; ++i) {
      vn_embeddings.push_back(new VNEmbedding(
          sa_sol.vn_embeddings[i].node_map, sa_sol.vn_embeddings[i].edge_map));
    }
    bottleneck_edges = sa_sol.bottleneck_edges;
    // Handle the handlers here.
    for (fibonacci_heap<bneck_edge_element_t>::iterator be_it =
             bottleneck_edges.begin();
         be_it != bottleneck_edges.end(); ++be_it) {
      bneck_heap_handlers[be_it->bneck_edge] =
          fibonacci_heap<bneck_edge_element_t>::s_handle_from_iterator(be_it);
    }

    vlinks_by_plength = sa_sol.vlinks_by_plength;
    // Handle the handlers here.
    for (fibonacci_heap<edge_plength_set_element_t>::iterator ep_it =
             vlinks_by_plength.begin();
         ep_it != vlinks_by_plength.end(); ++ep_it) {
      vlp_heap_handlers[ep_it->vn_edge] =
          fibonacci_heap<edge_plength_set_element_t>::s_handle_from_iterator(
              ep_it);
    }

    node_bw_usage = sa_sol.node_bw_usage;
    // Handle the handlers here.
    for (fibonacci_heap<node_bw_set_element_t>::iterator nb_it =
             node_bw_usage.begin();
         nb_it != node_bw_usage.end(); ++nb_it) {
      node_bw_heap_handlers[nb_it->vnode] =
          fibonacci_heap<node_bw_set_element_t>::s_handle_from_iterator(nb_it);
    }
  }

  SASolution(const Graph* pt, const ptr_vector<Graph>& vts,
             const ptr_vector<std::vector<std::vector<int> > >& lc,
             const VNRParameters& params, const ptr_vector<VNEmbedding>& vnes,
             const ReverseEmbedding* re, int nvns)
      : physical_topology(pt),
        virt_topologies(vts),
        location_constraints(lc),
        vnr_params(params),
        num_vns(nvns),
        util_matrix(pt->node_count(), pt->node_count(), 0.0),
        res_channel_matrix(pt->node_count(), pt->node_count()) {
    r_embedding = unique_ptr<ReverseEmbedding>(
        new ReverseEmbedding(re->rnode_map, re->redge_map));

    // Initialize residual channel matrix.
    for (int u = 0; u < physical_topology->node_count(); ++u) {
      std::vector<edge_endpoint>::const_iterator end_point_it;
      const std::vector<edge_endpoint>& u_neighbors =
          physical_topology->adj_list()->at(u);
      for (end_point_it = u_neighbors.begin();
           end_point_it != u_neighbors.end(); ++end_point_it) {
        int v = end_point_it->node_id;
        for (int w = 0; w < end_point_it->channels; ++w) {
          res_channel_matrix.matrix[u][v].push_back(w);
        }
      }
    }

    // Populate utilization and residual channel matrix.
    for (int i = 0; i < num_vns; ++i) {
      vn_embeddings.push_back(
          new VNEmbedding(vnes[i].node_map, vnes[i].edge_map));
      dwdm_edge_map_t::const_iterator emap_it;
      for (emap_it = vn_embeddings[i].edge_map.begin();
           emap_it != vn_embeddings[i].edge_map.end(); ++emap_it) {
        const edge_t& vlink = emap_it->first;
        int ch_id = emap_it->second.first;
        const path_t& plinks = emap_it->second.second;
        // Populate the set of vlinks sorted by mapped physical path length.
        vlp_heap_handlers[vn_edge_t(i, vlink)] = vlinks_by_plength.push(
            edge_plength_set_element_t(vn_edge_t(i, vlink), plinks.size()));
        // Populate the list of vnodes sorted according to total channel
        // consumption.
        vnode_t nb_key(i, vlink.first);
        if (node_bw_heap_handlers.find(nb_key) == node_bw_heap_handlers.end()) {
          long total_cost = 0;
          std::vector<edge_endpoint>::const_iterator vend_point_it;
          const std::vector<edge_endpoint>& neighbors =
              virt_topologies[i].adj_list()->at(vlink.first);
          for (vend_point_it = neighbors.begin();
               vend_point_it != neighbors.end(); ++vend_point_it) {
            long ch = vend_point_it->channels;
            const dwdm_path_t& mapped_path = vn_embeddings[i].edge_map[vlink];
            path_t::const_iterator path_it;
            for (path_it = mapped_path.second.begin();
                 path_it != mapped_path.second.end(); ++path_it) {
              total_cost += ch * physical_topology->GetEdgeCost(
                                     path_it->first, path_it->second);
            }
          }
          node_bw_set_element_t nb_value(nb_key, total_cost);
          node_bw_heap_handlers[nb_key] = node_bw_usage.push(nb_value);
        }
        nb_key.second = vlink.second;
        if (node_bw_heap_handlers.find(nb_key) == node_bw_heap_handlers.end()) {
          long total_cost = 0;
          std::vector<edge_endpoint>::const_iterator vend_point_it;
          const std::vector<edge_endpoint>& neighbors =
              virt_topologies[i].adj_list()->at(vlink.second);
          for (vend_point_it = neighbors.begin();
               vend_point_it != neighbors.end(); ++vend_point_it) {
            long ch = vend_point_it->channels;
            path_t::const_iterator path_it;
            const dwdm_path_t& mapped_path = vn_embeddings[i].edge_map[vlink];
            for (path_it = mapped_path.second.begin();
                 path_it != mapped_path.second.end(); ++path_it) {
              total_cost += ch * physical_topology->GetEdgeCost(
                                     path_it->first, path_it->second);
            }
          }
          node_bw_set_element_t nb_value(nb_key, total_cost);
          node_bw_heap_handlers[nb_key] = node_bw_usage.push(nb_value);
        }
        path_t::const_iterator plink_it;
        for (plink_it = plinks.begin(); plink_it != plinks.end(); ++plink_it) {
          int u = plink_it->first, v = plink_it->second;
          // Populate residual channel matrix.
          std::vector<int>::iterator it =
              std::find(res_channel_matrix.matrix[u][v].begin(),
                        res_channel_matrix.matrix[u][v].end(), ch_id);
          if (it == res_channel_matrix.matrix[u][v].end()) {
            printf("!!!!BANG!!!! Culprit: %d %d %d %d %d %d\n", i, u, v,
                   vlink.first, vlink.second, ch_id);
          } else {
            res_channel_matrix.matrix[u][v].erase(it);
            it = std::find(res_channel_matrix.matrix[v][u].begin(),
                           res_channel_matrix.matrix[v][u].end(), ch_id);
            if (it != res_channel_matrix.matrix[v][u].end()) {
              res_channel_matrix.matrix[v][u].erase(it);
            }
          }
        }
      }

      // Populate utilization matrix.
      for (int u = 0; u < physical_topology->node_count(); ++u) {
        std::vector<edge_endpoint>::const_iterator end_point_it;
        const std::vector<edge_endpoint>& u_neighbors =
            physical_topology->adj_list()->at(u);
        for (end_point_it = u_neighbors.begin();
             end_point_it != u_neighbors.end(); ++end_point_it) {
          int v = end_point_it->node_id;
          long res_ch = res_channel_matrix.matrix[u][v].size();
          long ch_uv = physical_topology->GetEdgeNumChannels(u, v);
          util_matrix.matrix[u][v] = util_matrix.matrix[v][u] =
              1.0 - static_cast<double>(res_ch) / static_cast<double>(ch_uv);
          // Populate the set of bottleneck physical links.
          if (util_matrix.matrix[u][v] > vnr_params.util_threshold) {
            if (bneck_heap_handlers.find(ConstructEdge(u, v)) ==
                bneck_heap_handlers.end()) {
              bneck_heap_handlers[ConstructEdge(u, v)] =
                  bottleneck_edges.push(bneck_edge_element_t(
                      ConstructEdge(u, v), util_matrix.matrix[u][v]));
            }
          }
        }
      }
    }
  }
};

#endif  // MIDDLEBOX_PLACEMENT_SRC_DATASTRUCTURE_H_
