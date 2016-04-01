#ifndef IO_H_
#define IO_H_

#include "datastructure.h"
#include "util.h"

#include <boost/move/unique_ptr.hpp>
#include <boost/move/utility.hpp>
#include <map>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>

using boost::movelib::unique_ptr;
typedef std::vector<std::vector<std::string> > csv_vector_t;
typedef unique_ptr<csv_vector_t> csv_vector_ptr_t;

unique_ptr<std::map<std::string, std::string> > ParseArgs(int argc,
                                                          char* argv[]) {
  unique_ptr<std::map<std::string, std::string> > arg_map(
      new std::map<std::string, std::string>());
  for (int i = 1; i < argc; ++i) {
    char* key = strtok(argv[i], "=");
    char* value = strtok(NULL, "=");
    DEBUG(" [%s] => [%s]\n", key, value);
    arg_map->insert(std::make_pair(key, value));
  }
  return boost::move(arg_map);
}

csv_vector_ptr_t ReadCSVFile(const char* filename) {
  DEBUG("[Parsing %s]\n", filename);
  FILE* file_ptr = fopen(filename, "r");
  if (!file_ptr) {
    DEBUG("Invalid file %s\n", filename);
    return csv_vector_ptr_t(NULL);
  }
  const static int kBufferSize = 1024;
  char line_buffer[kBufferSize];
  csv_vector_ptr_t ret_vector(new csv_vector_t());
  std::vector<std::string> current_line;
  int row_number = 0;
  while (fgets(line_buffer, kBufferSize, file_ptr)) {
    DEBUG("Read %d characters\n", strlen(line_buffer));
    if (strlen(line_buffer) <= 0) continue;
    if (line_buffer[0] == '\n' || line_buffer[0] == '\r') continue;
    current_line.clear();
    char* token = strtok(line_buffer, ",\n\r");
    current_line.push_back(token);
    while ((token = strtok(NULL, ",\n"))) {
      current_line.push_back(token);
    }
    ret_vector->push_back(current_line);
  }
  fclose(file_ptr);
  DEBUG("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
  return boost::move(ret_vector);
}

unique_ptr<Graph> InitializeTopologyFromFile(const char* filename) {
  int node_count = 0, edge_count = 0;
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return unique_ptr<Graph>(NULL);
  }
  unique_ptr<Graph> graph(new Graph());
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);

    // Each line has the following format:
    // LinkID, SourceID, DestinationID, PeerID, Cost, Bandwidth, Delay.
    int u = atoi(row[1].c_str());
    int v = atoi(row[2].c_str());
    int cost = atoi(row[4].c_str());
    long bw = atol(row[5].c_str());
    int delay = atoi(row[6].c_str());

    DEBUG("Line[%d]: u = %d, v = %d, cost = %d, bw = %ld, delay = %d\n", i, u,
          v, cost, bw, delay);
    graph->AddEdge(u, v, bw, delay, cost);
  }
  return boost::move(graph);
}

unique_ptr<std::vector<std::vector<int> > > InitializeVNLocationsFromFile(
    const char* filename, int num_virtual_nodes) {
  DEBUG("Parsing %s\n", filename);
  unique_ptr<std::vector<std::vector<int> > > ret_vector(
      new std::vector<std::vector<int> >(num_virtual_nodes));
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return unique_ptr<std::vector<std::vector<int> > >(NULL);
  }
  DEBUG("Parsing %s successful\n", filename);
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);
    int vnode_id = atoi(row[0].c_str());
    for (int j = 1; j < row.size(); ++j) {
      ret_vector->at(vnode_id).push_back(atoi(row[j].c_str()));
    }
  }
  return boost::move(ret_vector);
}

unique_ptr<VNEmbedding> InitializeVNEmbeddingFromFile(const char* nmap_file,
                                                      const char* emap_file) {
  DEBUG("Reading node embedding from: %s\n", nmap_file);
  DEBUG("Reading edge embedding from: %s\n", emap_file);
  unique_ptr<VNEmbedding> vn_embedding(new VNEmbedding());
  FILE* nmap = fopen(nmap_file, "r");
  FILE* emap = fopen(emap_file, "r");
  if (!nmap || !emap) {
    return unique_ptr<VNEmbedding>(NULL);
  }
  char buf[256];
  while (fgets(buf, sizeof(buf), nmap) != NULL) {
    int vnode, vnode_map;
    sscanf(buf, "%d %d", &vnode, &vnode_map);
    if (vnode > static_cast<int>(vn_embedding->node_map.size()) - 1)
      vn_embedding->node_map.resize(vnode + 1);
    vn_embedding->node_map[vnode] = vnode_map;
  }
  fclose(nmap);
  while (fgets(buf, sizeof(buf), emap) != NULL) {
    int u, v, m, n;
    sscanf(buf, "%d %d %d %d", &u, &v, &m, &n);
    if (u > v) std::swap(u, v);
    if (m > n) std::swap(m, n);
    edge_t vlink(m, n), plink(u, v);
    if (vn_embedding->edge_map.find(vlink) == vn_embedding->edge_map.end()) {
      vn_embedding->edge_map[vlink] = path_t();
    }
    vn_embedding->edge_map[vlink].push_back(plink);
  }
  fclose(emap);
  DEBUG("Embedding read successfully\n");
  return boost::move(vn_embedding);
}

VNRParameters InitializeParametersFromFile(const char* parameter_file) {
  VNRParameters parameters;
  FILE* param_file = fopen(parameter_file, "r");
  if (!param_file) {
    return parameters;
  }
  char buffer[256];
  const char* prefix[] = {"Goal Utilization", "alpha", "beta"};
  enum {
    UTIL = 0,
    ALPHA,
    BETA
  };

  while (fgets(buffer, sizeof(buffer), param_file) != NULL) {
    for (int i = 0; i < 4; ++i) {
      if (strncmp(buffer, prefix[i], strlen(prefix[i])) == 0) {
        switch (i) {
          case UTIL:
            sscanf(buffer + strlen(prefix[i]) + 2, "%lf",
                   &parameters.util_threshold);
            parameters.util_threshold /= 100.0;
            break;
          case ALPHA:
            sscanf(buffer + strlen(prefix[i]) + 2, "%lf", &parameters.alpha);
            break;
          case BETA:
            sscanf(buffer + strlen(prefix[i]) + 2, "%lf", &parameters.beta);
            break;
          default:
            DEBUG("Invalid parameter specified in %s\n", parameter_file);
        }
        break;
      }
    }
  }
  fclose(param_file);
  return parameters;
}
#endif  // IO_H_
