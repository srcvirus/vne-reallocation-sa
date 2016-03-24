#ifndef _SIMULATED_ANNEALING_H_
#define _SIMULATED_ANNEALING_H_

#include "datastructure.h"
#include "graph_util.h"
#include "sa_util.h"
#include "util.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <math.h>

double GetCostDelta(SASolution* sa_sol_ptr,
                    const std::vector<vn_edge_t>& vlinks,
                    const std::vector<dwdm_path_t>& old_embedded_paths,
                    const std::vector<dwdm_path_t>& new_embedded_paths) {
  double cost_delta = 0.0;
  int n_bottleneck_delta = 0;
  for (int vl_index = 0; vl_index < vlinks.size(); ++vl_index) {
    const vn_edge_t& vlink = vlinks[vl_index];
    int vn_index = vlink.first;
    int m = vlink.second.first, n = vlink.second.second;
    long required_channels = 
      sa_sol_ptr->virt_topologies[vn_index].GetEdgeNumChannels(m, n);
    path_t::const_iterator path_it;
    for (path_it = old_embedded_paths[vl_index].second.begin();
         path_it != old_embedded_paths[vl_index].second.end(); ++path_it) {
      int u = path_it->first, v = path_it->second;
      // Subtract the cost of old embedding first.
      cost_delta -= (required_channels * 
                      sa_sol_ptr->physical_topology->GetEdgeCost(u, v));
      if (sa_sol_ptr->bneck_heap_handlers.find(*path_it) !=
          sa_sol_ptr->bneck_heap_handlers.end()) {
        DEBUG("(%d, %d) was a bottleneck\n", u, v);
        if (sa_sol_ptr->util_matrix.matrix[u][v] <=
            sa_sol_ptr->vnr_params.util_threshold) {
          --n_bottleneck_delta;
          sa_sol_ptr->bottleneck_edges.erase(
              sa_sol_ptr->bneck_heap_handlers[*path_it]);
          sa_sol_ptr->bneck_heap_handlers.erase(*path_it);
          DEBUG("Old bottleneck (%d, %d) eleminated\n", u, v);
        }
      } else if (sa_sol_ptr->util_matrix.matrix[u][v] >
                 sa_sol_ptr->vnr_params.util_threshold) {
        sa_sol_ptr->bneck_heap_handlers[*path_it] =
            sa_sol_ptr->bottleneck_edges.push(bneck_edge_element_t(
                *path_it, sa_sol_ptr->util_matrix.matrix[u][v]));
        ++n_bottleneck_delta;
        DEBUG("New bottleneck created (%d, %d) with util = %lf\n", u, v,
              sa_sol_ptr->util_matrix.matrix[u][v]);
      }
    }

    for (path_it = new_embedded_paths[vl_index].second.begin();
         path_it != new_embedded_paths[vl_index].second.end(); ++path_it) {
      int u = path_it->first, v = path_it->second;
      cost_delta += required_channels * 
                      sa_sol_ptr->physical_topology->GetEdgeCost(u, v);
      if (sa_sol_ptr->bneck_heap_handlers.find(*path_it) !=
          sa_sol_ptr->bneck_heap_handlers.end()) {
        if (sa_sol_ptr->util_matrix.matrix[u][v] <=
            sa_sol_ptr->vnr_params.util_threshold) {
          --n_bottleneck_delta;
          sa_sol_ptr->bottleneck_edges.erase(
              sa_sol_ptr->bneck_heap_handlers[*path_it]);
          sa_sol_ptr->bneck_heap_handlers.erase(*path_it);
        }
      } else if (sa_sol_ptr->util_matrix.matrix[u][v] >
                 sa_sol_ptr->vnr_params.util_threshold) {
        sa_sol_ptr->bneck_heap_handlers[*path_it] =
            sa_sol_ptr->bottleneck_edges.push(bneck_edge_element_t(
                *path_it, sa_sol_ptr->util_matrix.matrix[u][v]));
        ++n_bottleneck_delta;
        DEBUG("new bottleneck created (%d, %d) with util = %lf\n", u, v,
              sa_sol_ptr->util_matrix.matrix[u][v]);
      }
    }
  }

  cost_delta *= sa_sol_ptr->vnr_params.alpha;
  cost_delta += (n_bottleneck_delta * sa_sol_ptr->vnr_params.beta);
  return cost_delta;
}

bool ReallocateBottleneckPLink(const SASolution* prev_sol,
                               SASolution* new_sol,
                               boost::random::mt19937& generator) {
  int total_bottlenecks = new_sol->bottleneck_edges.size();
  int nTopTenPercent = total_bottlenecks;
  // static_cast<int>(ceil(static_cast<double>(total_bottlenecks) / 10.0));
  boost::random::uniform_int_distribution<> dist(0, nTopTenPercent - 1);
  int plink_rank = dist(generator);
  // DEBUG("Removing %d-th top bottlneck link out of %d\n", plink_rank,
  //      new_sol->bottleneck_edges.size());

  // Get the plink_rank-th most utilized bottleneck link from the PQ.
  fibonacci_heap<bneck_edge_element_t>::ordered_iterator plink_it =
      new_sol->bottleneck_edges.ordered_begin();
  while (--plink_rank > 0) {
    ++plink_it;
  }
  int bneck_u = plink_it->bneck_edge.first,
      bneck_v = plink_it->bneck_edge.second;
  DEBUG("Reallocating bneck plink = (%d, %d) with util %lf\n", bneck_u, bneck_v,
        new_sol->util_matrix.matrix[bneck_u][bneck_v]);
  bool reallocation_possible = false;

  while (new_sol->util_matrix.matrix[bneck_u][bneck_v] >
         new_sol->vnr_params.util_threshold) {
    // Find the VN on the selected plink which in turn has an embedding with the
    // most number of bottleneck plinks.
    const vedge_set_t& mapped_vlinks =
        new_sol->r_embedding->redge_map[plink_it->bneck_edge];
    vedge_set_t::const_iterator vlink_it;
    vn_edge_t selected_vlink(NIL, edge_t(NIL, NIL));
    int max_bottlenecks = 0;
    for (vlink_it = mapped_vlinks.begin(); vlink_it != mapped_vlinks.end();
         ++vlink_it) {
      int n_bottlenecks = GetNumBottleneckLinksOnPath(
          new_sol->util_matrix, new_sol->vn_embeddings[vlink_it->first]
                                    .edge_map.find(vlink_it->second)
                                    ->second.second,
          new_sol->vnr_params.util_threshold);
      DEBUG("vn = %d, vlink (%d, %d), n_bottlenecks = %d\n", vlink_it->first,
            vlink_it->second.first, vlink_it->second.second, n_bottlenecks);
      if (n_bottlenecks > max_bottlenecks) {
        max_bottlenecks = n_bottlenecks;
        selected_vlink = *vlink_it;
      }
    }

    if (selected_vlink.first != NIL) {
      int vn_index = selected_vlink.first;
      const edge_t& vlink = selected_vlink.second;
      int m = vlink.first, n = vlink.second;
      DEBUG("Reallocating vlink (%d, %d) from vn %d\n", m, n, vn_index);
      const Graph& virt_topology = new_sol->virt_topologies[vn_index];
      VNEmbedding& embedding = new_sol->vn_embeddings[vn_index];
      int ph_src = embedding.node_map[m];
      int ph_dest = embedding.node_map[n];
      int required_channels = virt_topology.GetEdgeNumChannels(m, n);
      std::set<edge_t> forbidden;
      forbidden.insert(plink_it->bneck_edge);
      dwdm_path_t& old_path = embedding.edge_map.find(vlink)->second;

      // Release channel from the old path first.
      // revisit
      UpdateResidualChannel(new_sol->res_channel_matrix, old_path, false);
      UpdateUtilMatrix(new_sol->physical_topology, new_sol->util_matrix,
                       new_sol->res_channel_matrix, old_path.second);

      // Find a new embedding for the vlink.
      unique_ptr<dwdm_path_t> embedded_path(ConstrainedDWDMVLinkEmbed(
          new_sol->physical_topology, new_sol->util_matrix,
          new_sol->res_channel_matrix, ph_src, ph_dest, required_channels, 
          new_sol->physical_topology->max_channels(), forbidden).release());
      if (embedded_path->second.size() <= 0) {
        DEBUG("Cannot find a path to embed (%d, %d) from vn %d\n", m, n,
              vn_index);
        // Since new path could not be found, restore bandwidth on the old path.
        // revisit
        UpdateResidualChannel(new_sol->res_channel_matrix, old_path, true);
        UpdateUtilMatrix(new_sol->physical_topology, new_sol->util_matrix,
                         new_sol->res_channel_matrix, old_path.second);
        break;
      }
      DEBUG("New path of length %d found\n", embedded_path->second.size());
      for (path_t::const_iterator pit = embedded_path->second.begin();
           pit != embedded_path->second.end(); ++pit) {
        DEBUG("new_path: (%d %d)\n", pit->first, pit->second);
      }

      // path_t& old_path = embedding.edge_map.find(vlink)->second;
      DEBUG("Old path had length %d\n", old_path.second.size());
      for (path_t::const_iterator pit = old_path.second.begin(); pit != old_path.second.end();
           ++pit) {
        DEBUG("old_path: (%d %d)\n", pit->first, pit->second);
      }

      UpdateReverseLinkEmbedding(new_sol->r_embedding.get(), selected_vlink,
                                 old_path.second, embedded_path->second);
      // UpdateResidualChannel(new_sol->res_channel_matrix, old_path, bw);
      // UpdateUtilMatrix(new_sol->physical_topology,
      //    new_sol->util_matrix, new_sol->res_channel_matrix, old_path);

      // Update residual bandwidth and utilization along the new path.
      // revisit
      UpdateResidualChannel(new_sol->res_channel_matrix, *embedded_path, true);
      UpdateUtilMatrix(new_sol->physical_topology, new_sol->util_matrix,
                       new_sol->res_channel_matrix, embedded_path->second);

      std::vector<vn_edge_t> vlinks_changed(1, selected_vlink);
      std::vector<dwdm_path_t> old_paths(1, old_path);
      std::vector<dwdm_path_t> new_paths(1, *embedded_path);
      double cost_delta =
          GetCostDelta(new_sol, vlinks_changed, old_paths, new_paths);
      UpdateReverseLinkEmbedding(new_sol->r_embedding.get(), selected_vlink,
                                 old_path.second, embedded_path->second);
      new_sol->vlinks_by_plength.update(
          new_sol->vlp_heap_handlers[selected_vlink],
          edge_plength_set_element_t(selected_vlink, embedded_path->second.size()));
      old_path.first = embedded_path->first;
      old_path.second.assign(embedded_path->second.begin(), embedded_path->second.end());
      DEBUG("Old cost = %lf, cost_delta = %lf\n", new_sol->cost, cost_delta);
      new_sol->cost += cost_delta;
      reallocation_possible = true;
      DEBUG("After vlink reallocation, util of plink (%d, %d) is %lf\n",
            bneck_u, bneck_v, new_sol->util_matrix.matrix[bneck_u][bneck_v]);
      DEBUG(
          "After vlink reallocation, embedded vlinks on plink (%d, %d) are:\n",
          bneck_u, bneck_v);
      const vedge_set_t& vet =
          new_sol->r_embedding->redge_map.find(edge_t(bneck_u, bneck_v))
              ->second;
      vedge_set_t::const_iterator vit;
      for (vit = vet.begin(); vit != vet.end(); ++vit) {
        DEBUG("Vlink (%d, %d) from vn %d\n", vit->second.first,
              vit->second.second, vit->first);
      }
    } else {
      break;
    }
  }
  return reallocation_possible;
}

bool PerformVLinkMigration(const SASolution* current_solution,
                           SASolution* neighbor_solution) {
  bool reallocation_feasible = false;
  fibonacci_heap<edge_plength_set_element_t>::ordered_iterator vlp_it =
      neighbor_solution->vlinks_by_plength.ordered_begin();
  std::vector<edge_plength_set_element_t> changed_vlps;
  int total_vlinks = neighbor_solution->vlinks_by_plength.size();
  int num_migration_candidates = total_vlinks / 10;
  int num_migration_attempt = 0;
  num_migration_candidates = 1;
  while (num_migration_attempt++ < num_migration_candidates) {
    edge_plength_set_element_t vlp = *vlp_it;
    ++vlp_it;
    //  edge_plength_set_element_t vlp =
    //    neighbor_solution->vlinks_by_plength.top();
    // neighbor_solution->vlinks_by_plength.pop();
    int vn_index = vlp.vn_edge.first;
    int vlink_u = vlp.vn_edge.second.first, vlink_v = vlp.vn_edge.second.second;
    DEBUG("Migrating (%d, %d) from vn %d\n", vlink_u, vlink_v, vn_index);
    int ph_src = neighbor_solution->vn_embeddings[vn_index].node_map[vlink_u];
    int ph_dest = neighbor_solution->vn_embeddings[vn_index].node_map[vlink_v];
    long required_channels = neighbor_solution->virt_topologies[vn_index]
                  .GetEdgeNumChannels(vlink_u, vlink_v);
    dwdm_path_t& old_path =
        neighbor_solution->vn_embeddings[vn_index].edge_map[vlp.vn_edge.second];
    // Release residual bandwidth on the old embedded path.
    // revisit
    UpdateResidualChannel(neighbor_solution->res_channel_matrix, old_path, false);
    UpdateUtilMatrix(neighbor_solution->physical_topology,
                     neighbor_solution->util_matrix,
                     neighbor_solution->res_channel_matrix, old_path.second);
    // Try to embed to plink on a new path.
    unique_ptr<dwdm_path_t> embedded_path(ConstrainedDWDMVLinkEmbed(
        neighbor_solution->physical_topology, neighbor_solution->util_matrix,
        neighbor_solution->res_channel_matrix, ph_src, ph_dest, required_channels,
        neighbor_solution->physical_topology->max_channels(), std::set<edge_t>()).release());

    if (embedded_path->second.size() > 0) {
      // path_t& old_path = neighbor_solution->vn_embeddings[vn_index]
      //                       .edge_map[vlp.vn_edge.second];
      // UpdateResidualChannel(neighbor_solution->res_channel_matrix, old_path,
      // bw);
      // UpdateUtilMatrix(neighbor_solution->physical_topology,
      //    neighbor_solution->util_matrix, neighbor_solution->res_channel_matrix,
      // old_path);

      // Update residual bandwidth and utilization along the newly found path.
      DEBUG("New path of length %d found\n", embedded_path->second.size());
      for (path_t::const_iterator pit = embedded_path->second.begin();
           pit != embedded_path->second.end(); ++pit) {
        DEBUG("new_path: (%d %d) : %d\n", pit->first, pit->second, 
            embedded_path->first);
      }

      DEBUG("Old path had length %d\n", old_path.second.size());
      for (path_t::const_iterator pit = old_path.second.begin(); 
           pit != old_path.second.end(); ++pit) {
        DEBUG("old_path: (%d %d) : %d\n", pit->first, pit->second,
            old_path.first);
      }
      // revisit
      UpdateResidualChannel(neighbor_solution->res_channel_matrix, *embedded_path, true);
      UpdateUtilMatrix(neighbor_solution->physical_topology,
                       neighbor_solution->util_matrix,
                       neighbor_solution->res_channel_matrix, embedded_path->second);
      std::vector<vn_edge_t> vlinks_changed(1, vlp.vn_edge);
      std::vector<dwdm_path_t> old_paths(1, old_path);
      std::vector<dwdm_path_t> new_paths(1, *embedded_path);
      double cost_delta =
          GetCostDelta(neighbor_solution, vlinks_changed, old_paths, new_paths);
      UpdateReverseLinkEmbedding(neighbor_solution->r_embedding.get(),
                                 vlp.vn_edge, old_path.second, embedded_path->second);
      old_path.first = embedded_path->first;
      old_path.second.assign(embedded_path->second.begin(), 
                             embedded_path->second.end());
      neighbor_solution->cost += cost_delta;
      reallocation_feasible = true;
      vlp.path_len = old_path.second.size();
      // Update bookkeeping information.
      changed_vlps.push_back(vlp);
      // neighbor_solution->vlp_heap_handlers[vlp.vn_edge] =
      //  neighbor_solution->vlinks_by_plength.push(
      //    edge_plength_set_element_t(vlp.vn_edge, embedded_path->size()));
    } else {
      // Since no new embedding could be found, restore bandwidth on the old
      // path.
      // revisit
      UpdateResidualChannel(neighbor_solution->res_channel_matrix, old_path, true);
      UpdateUtilMatrix(neighbor_solution->physical_topology,
                       neighbor_solution->util_matrix,
                       neighbor_solution->res_channel_matrix, old_path.second);
    }
  }

  for (int i = 0; i < changed_vlps.size(); ++i) {
    edge_plength_set_element_t& vlp = changed_vlps[i];
    neighbor_solution->vlinks_by_plength.update(
        neighbor_solution->vlp_heap_handlers[vlp.vn_edge], vlp);
  }
  return reallocation_feasible;
}

bool PerformNodeMigration(const SASolution* current_solution,
                          SASolution* neighbor_solution,
                          boost::random::mt19937& generator) {
  bool reallocation_feasible = false;
  int num_total_vnodes = neighbor_solution->node_bw_usage.size();
  boost::random::uniform_int_distribution<> dist(0, num_total_vnodes - 1);
  int vnode_to_move = dist(generator);
  fibonacci_heap<node_bw_set_element_t>::ordered_iterator nbit =
      neighbor_solution->node_bw_usage.ordered_begin();
  while (--vnode_to_move > 0) ++nbit;
  if (nbit == neighbor_solution->node_bw_usage.ordered_end()) 
    return reallocation_feasible;
  node_bw_set_element_t nbw = *nbit;
  neighbor_solution->node_bw_usage.erase(
      neighbor_solution->node_bw_heap_handlers[nbw.vnode]);
  neighbor_solution->node_bw_heap_handlers.erase(nbw.vnode);
  int vn_index = nbw.vnode.first;
  int vnode_u = nbw.vnode.second;
  int vnode_u_map =
      neighbor_solution->vn_embeddings[vn_index].node_map[vnode_u];
  DEBUG("Reallocating vnode %d of vn %d from pnode %d\n", vnode_u, vn_index,
        vnode_u_map);
  int best_u_map = NIL;
  long min_total_path_cost = INF;
  long prev_cost = nbw.bw_usage;

  for (int i = 0;
       i < neighbor_solution->location_constraints[vn_index][vnode_u].size();
       ++i) {
    int temp_u_map =
        neighbor_solution->location_constraints[vn_index][vnode_u][i];
    if (temp_u_map == vnode_u_map) continue;
    std::vector<edge_endpoint>::const_iterator vend_point_it;
    ptr_vector<dwdm_path_t> embedded_paths;
    bool embedding_feasible = true;
    long total_path_cost = 0;
    const std::vector<edge_endpoint>& neighbors =
        neighbor_solution->virt_topologies[vn_index].adj_list()->at(vnode_u);
    for (vend_point_it = neighbors.begin(); vend_point_it != neighbors.end();
         ++vend_point_it) {
      int vnode_v = vend_point_it->node_id;
      int required_channels = vend_point_it->channels;
      int vnode_v_map =
          neighbor_solution->vn_embeddings[vn_index].node_map[vnode_v];
      embedded_paths.push_back(ConstrainedDWDMVLinkEmbed(
          neighbor_solution->physical_topology, neighbor_solution->util_matrix,
          neighbor_solution->res_channel_matrix, temp_u_map, vnode_v_map, 
          required_channels, neighbor_solution->physical_topology->max_channels(),
          std::set<edge_t>()).release());
      if (embedded_paths.back().second.size() <= 0) {
        embedding_feasible = false;
        break;
      }
      const dwdm_path_t& recent_path = embedded_paths.back();
      total_path_cost += GetEmbeddedPathCost(
          vn_edge_t(vn_index, ConstructEdge(vnode_u, vnode_v)), 
          recent_path.second, neighbor_solution);
    }
    if (!embedding_feasible) continue;
    if (total_path_cost < min_total_path_cost) {
      min_total_path_cost = total_path_cost;
      best_u_map = temp_u_map;
    }
  }
  if (best_u_map != NIL) {
    DEBUG("Reallocating vnode %d of vn %d from pnode %d to pnode %d\n", vnode_u,
          vn_index, vnode_u_map, best_u_map);
    const std::vector<edge_endpoint>& neighbors =
        neighbor_solution->virt_topologies[vn_index].adj_list()->at(vnode_u);
    std::vector<edge_endpoint>::const_iterator vend_point_it;
    std::vector<dwdm_path_t> old_paths, new_paths;
    std::vector<vn_edge_t> changed_vlinks;
    for (vend_point_it = neighbors.begin(); vend_point_it != neighbors.end();
         ++vend_point_it) {
      int vnode_v = vend_point_it->node_id;
      long required_channels = vend_point_it->channels;
      DEBUG("Moving vlink (%d, %d) of vn %d\n", vnode_u, vnode_v, vn_index);
      int vnode_v_map =
          neighbor_solution->vn_embeddings[vn_index].node_map[vnode_v];
      DEBUG("Finding a path between pnodes %d --> %d\n", best_u_map,
            vnode_v_map);
      dwdm_path_t& old_path = neighbor_solution->vn_embeddings[vn_index]
                             .edge_map[ConstructEdge(vnode_u, vnode_v)];
      // revisit
      UpdateResidualChannel(neighbor_solution->res_channel_matrix, old_path, false);
      UpdateUtilMatrix(neighbor_solution->physical_topology,
                       neighbor_solution->util_matrix,
                       neighbor_solution->res_channel_matrix, old_path.second);
      unique_ptr<dwdm_path_t> embedded_path(ConstrainedDWDMVLinkEmbed(
          neighbor_solution->physical_topology, neighbor_solution->util_matrix,
          neighbor_solution->res_channel_matrix, best_u_map, vnode_v_map,
          required_channels, neighbor_solution->physical_topology->max_channels(), 
          std::set<edge_t>()).release());
      vn_edge_t vl(vn_index, ConstructEdge(vnode_u, vnode_v));
      changed_vlinks.push_back(vl);
      for (int x = 0; x < old_path.second.size(); ++x) {
        DEBUG("[old_path]: (%d, %d)\n", old_path.second[x].first, old_path.second[x].second);
      }
      for (int x = 0; x < embedded_path->second.size(); ++x) {
        DEBUG("[new_path]: (%d, %d)\n", embedded_path->second[x].first,
              embedded_path->second[x].second);
      }
      old_paths.push_back(old_path);
      new_paths.push_back(*embedded_path);
      // UpdateResidualChannel(neighbor_solution->res_channel_matrix, old_path,
      //                        bw);
      // UpdateUtilMatrix(neighbor_solution->physical_topology,
      //    neighbor_solution->util_matrix, neighbor_solution->res_channel_matrix,
      // old_path);
      // revisit
      UpdateResidualChannel(neighbor_solution->res_channel_matrix, *embedded_path, true);
      UpdateUtilMatrix(neighbor_solution->physical_topology,
                       neighbor_solution->util_matrix,
                       neighbor_solution->res_channel_matrix, embedded_path->second);
      UpdateReverseLinkEmbedding(neighbor_solution->r_embedding.get(), vl,
                                 old_path.second, embedded_path->second);
      old_path.first = embedded_path->first;
      old_path.second.assign(embedded_path->second.begin(), 
                             embedded_path->second.end());
      neighbor_solution->vlinks_by_plength.update(
          neighbor_solution->vlp_heap_handlers[vl],
          edge_plength_set_element_t(vl, embedded_path->second.size()));
    }
    int old_nmap = neighbor_solution->vn_embeddings[vn_index].node_map[vnode_u];
    neighbor_solution->vn_embeddings[vn_index].node_map[vnode_u] = best_u_map;
    UpdateReverseNodeEmbedding(neighbor_solution->r_embedding.get(),
                               vnode_t(vn_index, vnode_u), old_nmap,
                               best_u_map);
    double cost_delta =
        GetCostDelta(neighbor_solution, changed_vlinks, old_paths, new_paths);
    DEBUG("Old cost = %lf, cost_delta = %lf\n", neighbor_solution->cost,
          cost_delta);
    neighbor_solution->cost += cost_delta;
    // Update bookkeeping information.
    neighbor_solution->node_bw_heap_handlers[nbw.vnode] =
        neighbor_solution->node_bw_usage.push(
            node_bw_set_element_t(nbw.vnode, min_total_path_cost));
    reallocation_feasible = true;
  }
  return reallocation_feasible;
}

unique_ptr<SASolution> GenerateNeighbor(
    const SASolution& current_solution, 
    boost::random::mt19937& generator) {
  unique_ptr<SASolution> neighbor_solution(new SASolution(current_solution));
  enum {
    ALTER_BNECK = 0,
    NODE_MIGRATE,
    LINK_MIGRATE
  };
  int n_possible_operations = 3;
  boost::random::uniform_int_distribution<> dist(
      0, n_possible_operations - 1);
  int current_operation = dist(generator);

  bool result;
  int total_bottlenecks = 0;
  // current_operation = ALTER_BNECK;
  switch (current_operation) {
    case ALTER_BNECK:
      total_bottlenecks = neighbor_solution->bottleneck_edges.size();
      if (total_bottlenecks > 0) {
        DEBUG("Performing Bneck link configuration\n");
        result = ReallocateBottleneckPLink(&current_solution,
                                           neighbor_solution.get(),
                                           generator);
        break;
      }
    case NODE_MIGRATE:
      DEBUG("Performing Node migration\n");
      result = PerformNodeMigration(&current_solution, neighbor_solution.get(), generator);
      break;
    case LINK_MIGRATE:
      DEBUG("Perform link migration\n");
      result =
          PerformVLinkMigration(&current_solution, neighbor_solution.get());
      break;
  }
  /*
  DEBUG("ALTER_BNECK\n");
  result = ReallocateBottleneckPLink(
      &current_solution, neighbor_solution.get());
  DEBUG("NODE_MIGRATE\n");
  result = PerformNodeMigration(
      &current_solution, neighbor_solution.get());
  DEBUG("LINK_MIGRATE\n");
  result = PerformVLinkMigration(
      &current_solution, neighbor_solution.get());
  */
  return boost::move(neighbor_solution);
}
#endif  // _SIMULATED_ANNEALING_H_
