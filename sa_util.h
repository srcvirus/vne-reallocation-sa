#ifndef _SA_UTIL_H_
#define _SA_UTIL_H_

inline void UpdateReverseLinkEmbedding(ReverseEmbedding* re,
                                       const vn_edge_t& vlink,
                                       const path_t& old_embedding,
                                       const path_t& new_embedding) {
  for (int i = 0; i < old_embedding.size(); ++i) {
    re->redge_map[old_embedding[i]].erase(vlink);
  }
  for (int i = 0; i < new_embedding.size(); ++i) {
    re->redge_map[new_embedding[i]].insert(vlink);
  }
}

inline void UpdateReverseNodeEmbedding(ReverseEmbedding* re,
                                       const vnode_t& vnode, int old_pnode,
                                       int new_pnode) {
  re->rnode_map[old_pnode].erase(vnode);
  re->rnode_map[new_pnode].insert(vnode);
}

inline long GetEmbeddedPathCost(const vn_edge_t& vlink,
                                const path_t& embedded_path,
                                const SASolution* sa_sol) {
  long total_cost = 0;
  int vn_index = vlink.first;
  const edge_t& vedge = vlink.second;
  long ch = sa_sol->virt_topologies[vn_index]
                .GetEdgeNumChannels(vedge.first, vedge.second);
  for (int i = 0; i < embedded_path.size(); ++i) {
    total_cost += ch * sa_sol->physical_topology->GetEdgeCost(
                           embedded_path[i].first, embedded_path[i].second);
  }
  return total_cost;
}
#endif  // _SA_UTIL_H_
