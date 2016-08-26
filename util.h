#ifndef UTIL_H_
#define UTIL_H_

#include "datastructure.h"

#include <vector>
#include <utility>

#define ONE_GIG 1000000000ULL
#define EPS 1e-8

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__) " "

#ifdef DBG
#define DEBUG(...) PrintDebugMessage(AT, __VA_ARGS__)
#else
#define DEBUG(...)
#endif

// Prints a message formatted using fmt_string. Prepends the location of the
// file before the message.
void PrintDebugMessage(const char* location, const char* fmt_string, ...);

// Returns the mean of the vector of data.
template <typename T>
double GetMean(const std::vector<T>& data);

// Returns the Nth percentile of the vector of data. Nth percentile is
// calculated as the (ceil(N / 100) * data.size() - 1)-th element from an array
// obtained by sorting data.
template <typename T>
T GetNthPercentile(const std::vector<T>& data, int n);

// Returns the Cumalitive Distribution Frequence of the data items sotred in
// data. If the data items sotred in data are of type double, then a precision
// of 3 digits after the decimal points is used.
template <typename T>
std::vector<std::pair<T, double> > GetCDF(const std::vector<T>& data);

// Inverses a set of VN Embedding and represents them in terms of Physical
// network entities.
unique_ptr<ReverseEmbedding> GetInverseEmbedding(
    const Graph* physical_topology, const ptr_vector<VNEmbedding>& embeddings,
    int num_vns);

// Returns the cost of allocating bandwidth for embedding all the vns in
// virt_topology on phys_topology.
unsigned long GetBandwidthCost(
    const Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings);

// Returns the cost of bottleneck physical links in the physical topology.
int GetNumBottleneckLinks(const Graph* phys_topology,
                          const boost::ptr_vector<Graph>& virt_topologies,
                          const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                          const VNRParameters* vnr_param);

// Returns the maximum physical link utilization.
double GetMaxPLinkUtilization(
    const Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings);

// Compute the cost of a set of embeddings using our cost function.
double CostFunction(const Graph* phys_topology,
                    const boost::ptr_vector<Graph>& virt_topologies,
                    const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                    const VNRParameters* vnr_param);

// If the input physical network contains residual bandwidth then compute the
// bandwidth capacity of the physical links by adding the bandwidths of embedded
// virtual links.
void ComputePhysicalNetworkCapacity(
    Graph* phys_topology, const boost::ptr_vector<Graph>& virt_topologies,
    const boost::ptr_vector<VNEmbedding>& vn_embeddings);

// Write VN embeddings to file.
void WriteSolutionToFile(const boost::ptr_vector<VNEmbedding>& vn_embeddings,
                         const std::string& vnr_directory,
                         const std::vector<int>& valid_indices);
#endif  // UTIL_H_
