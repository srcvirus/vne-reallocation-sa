#include "datastructure.h"
#include "io.h"
#include "simulated_annealing.h"
#include "util.h"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

const std::string kUsage =
    "./vne_reallocation "
    "--case_directory=<case_directory>\n"
    "\t[--max_iterations=<max_iterations>]\n"
    "\t[--iterations_per_temperature=<iterations_per_temperature>]";

// Maximum number of main simulated annealing iterations to perform. Can also be
// set by command line option --max_iterations
int max_iterations = 300;

// Number of iterations to perform for a fixed temperature. Can also be set by
// command line option --iterations_per_temperature
int iterations_per_temperature = 150;

// Execution thread for simulated annealing search. It takes as parameter a
// pointer to an initial solution and returns the best solution after the search
// is complete.
void* SimulatedAnnealingThread(void* args) {
  SASolution* initial = reinterpret_cast<SASolution*>(args);
  unique_ptr<SASolution> current_solution(new SASolution(*initial));
  SASolution* best_solution = new SASolution(*initial);
  double best_cost = initial->cost;
  double temparature = 0.95;
  double ro = 0.99;
  int k = 0;
  boost::random::mt19937 seed(0x5414ab);
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(
      seed, dist);
  do {
    DEBUG("Iteration = %d\n", k);
    DEBUG("Current solution cost = %lf\n", current_solution->cost);
    int l = 0;
    do {
      unique_ptr<SASolution> next(
          GenerateNeighbor(*current_solution, seed).release());
      double cost_difference = next->cost - current_solution->cost;
      double energy_value = exp(-cost_difference / temparature);
      double rand_val = random();
      DEBUG("Cost difference = %lf, Energy value = %lf, rand_val = %lf\n",
            cost_difference, energy_value, rand_val);
      if (cost_difference < 0.0) {
        current_solution.reset(next.release());
        DEBUG("Better taken\n");
      } else if (rand_val < energy_value) {
        current_solution.reset(next.release());
        DEBUG("Worse taken\n");
      }
      if (best_cost > current_solution->cost) {
        best_cost = current_solution->cost;
        delete best_solution;
        best_solution = new SASolution(*current_solution);
        DEBUG("Best updated, current best cost = %lf\n", best_cost);
      }
    } while (++l < iterations_per_temperature);
    temparature *= ro;
  } while (++k < max_iterations);
  pthread_exit(reinterpret_cast<void*>(best_solution));
}

int main(int argc, char* argv[]) {
  // Parse the command line args.
  using std::string;
  unique_ptr<std::map<string, string> > arg_map(
      ParseArgs(argc, argv).release());
  string case_directory = "";
  std::map<string, string>::iterator arg_map_it;
  for (arg_map_it = arg_map.get()->begin(); arg_map_it != arg_map.get()->end();
       ++arg_map_it) {
    if (arg_map_it->first == "--case_directory") {
      case_directory = arg_map_it->second;
    } else if (arg_map_it->first == "--max_iterations") {
      max_iterations = atoi(arg_map_it->second.c_str());
    } else if (arg_map_it->first == "--iterations_per_temperature") {
      iterations_per_temperature = atoi(arg_map_it->second.c_str());
    } else {
      printf("Invalid command line option: %s\n", arg_map_it->first.c_str());
      printf("Usage: %s\n", kUsage.c_str());
      return 1;
    }
  }

  // Parse the topologies and embeddings.
  const string kPhysicalTopologyFile = case_directory + "/sn.txt";
  unique_ptr<Graph> physical_topology(
      InitializeTopologyFromFile(kPhysicalTopologyFile.c_str()).release());
  int num_vns = 0;
  boost::ptr_vector<Graph> virt_topologies;
  boost::ptr_vector<std::vector<std::vector<int> > > location_constraints;
  boost::ptr_vector<VNEmbedding> vn_embeddings;

  // Read the VNs and their embeddings.
  while (true) {
    const string kVirtTopologyFile = case_directory + "/vnr/vn" +
                                     boost::lexical_cast<string>(num_vns) +
                                     ".txt";
    const string kVNLocationConstraintFile =
        case_directory + "/vnr/vnloc" + boost::lexical_cast<string>(num_vns) +
        ".txt";
    const string kVLinkEmbeddingFile = kVirtTopologyFile + ".semap";
    const string kVNodeEmbeddingFile = kVirtTopologyFile + ".nmap";
    unique_ptr<Graph> virt_topology(
        InitializeTopologyFromFile(kVirtTopologyFile.c_str()).release());
    if (virt_topology.get() == NULL) {
      break;
    }
    virt_topologies.push_back(virt_topology.release());
    virt_topologies[num_vns].Matrixize();
    DEBUG(virt_topologies[num_vns].GetDebugString().c_str());
    location_constraints.push_back(InitializeVNLocationsFromFile(
        kVNLocationConstraintFile.c_str(),
        virt_topologies[num_vns].node_count()).release());
    vn_embeddings.push_back(InitializeVNEmbeddingFromFile(
        kVNodeEmbeddingFile.c_str(), kVLinkEmbeddingFile.c_str()).release());
    ++num_vns;
  }

  // Compute the physical network capacity from the VN embeddings and residual
  // capacity.
  ComputePhysicalNetworkCapacity(physical_topology.get(), virt_topologies,
                                 vn_embeddings);
  // Convert the internal adjacency list to adjacency matrix for fast lookup.
  physical_topology->Matrixize();

  VNRParameters vnr_parameters = InitializeParametersFromFile(
      (case_directory + "/optimize_para.txt").c_str());
  unique_ptr<ReverseEmbedding> reverse_embedding(GetInverseEmbedding(
      physical_topology.get(), vn_embeddings, num_vns).release());

  // log old valus first.
  double prev_cost = CostFunction(physical_topology.get(), virt_topologies,
                                  vn_embeddings, &vnr_parameters);
  int prev_num_bottlenecks = GetNumBottleneckLinks(
      physical_topology.get(), virt_topologies, vn_embeddings, &vnr_parameters);
  long prev_bw_cost = static_cast<long>(
      (prev_cost - vnr_parameters.beta * prev_num_bottlenecks) /
      vnr_parameters.alpha);
  double prev_max_util = GetMaxPLinkUtilization(physical_topology.get(),
                                                virt_topologies, vn_embeddings);

  FILE* f = fopen((case_directory + "/vnr/prev_cost").c_str(), "w");
  fprintf(f, "%lf\n", prev_cost);
  fclose(f);
  f = fopen((case_directory + "/vnr/prev_bnecks").c_str(), "w");
  fprintf(f, "%d\n", prev_num_bottlenecks);
  fclose(f);
  f = fopen((case_directory + "/vnr/prev_bw_cost").c_str(), "w");
  fprintf(f, "%ld\n", prev_bw_cost);
  fclose(f);
  f = fopen((case_directory + "/vnr/prev_max_plink_util").c_str(), "w");
  fprintf(f, "%lf\n", prev_max_util);
  fclose(f);
  
  // Start a timer.
  std::time_t start = std::time(NULL);

  // Generate the initial solutions and start the solution threads.
  ptr_vector<SASolution> initial_solutions;
  initial_solutions.push_back(new SASolution(
      physical_topology.get(), virt_topologies, location_constraints,
      vnr_parameters, vn_embeddings, reverse_embedding.get(), num_vns));
  initial_solutions[0].cost = CostFunction(
      physical_topology.get(), virt_topologies, vn_embeddings, &vnr_parameters);
  const int kNumCores = sysconf(_SC_NPROCESSORS_ONLN);
  // const int kNumThreads = kNumCores;
  const int kNumThreads = 1;
  std::vector<pthread_t> threads(kNumThreads);
  boost::random::mt19937 seed(0x5414ab);
  for (int i = 0, thread_id = 0, current_core = 0; i < kNumThreads;
       ++i, ++thread_id, current_core = (current_core + 1) % kNumCores) {
    if (i > 0) {
      initial_solutions.push_back(
          GenerateNeighbor(initial_solutions[i - 1], seed).release());
    }
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(current_core, &cpu_set);
    pthread_create(&threads[thread_id], NULL, &SimulatedAnnealingThread,
                   &initial_solutions[i]);
    pthread_setaffinity_np(threads[thread_id], sizeof(cpu_set), &cpu_set);
  }

  // Wait for the threads to finish and collect each solution.
  ptr_vector<SASolution> solutions;
  for (int i = 0; i < kNumThreads; ++i) {
    void* ret_value;
    pthread_join(threads[i], &ret_value);
    solutions.push_back(reinterpret_cast<SASolution*>(ret_value));
  }
  double best_solution_cost = solutions[0].cost;
  SASolution* best_solution = &solutions[0];
  for (int i = 0; i < solutions.size(); ++i) {
    if (solutions[i].cost < best_solution_cost) {
      best_solution_cost = solutions[i].cost;
      best_solution = &solutions[i];
    }
  }

  // Stop the timer.
  std::time_t end = std::time(NULL);
  double duration = std::difftime(end, start);

  // Log the new values.
  int new_bnecks =
      GetNumBottleneckLinks(physical_topology.get(), virt_topologies,
                            best_solution->vn_embeddings, &vnr_parameters);
  printf("B.Necks = %ld\n", best_solution->bottleneck_edges.size());
  printf("Exec time = %lfs\n", duration);
  double new_cost = best_solution->cost;
  long new_bw_cost = static_cast<long>(
      (new_cost - vnr_parameters.beta * new_bnecks) / vnr_parameters.alpha);
  double new_max_util = GetMaxPLinkUtilization(physical_topology.get(),
                                               virt_topologies, vn_embeddings);

  DEBUG("new_bnecks = %d, new_cost = %lf, new_bw_cost = %ld\n", new_bnecks,
        new_cost, new_bw_cost);
  f = fopen((case_directory + "/vnr/sol_time").c_str(), "w");
  fprintf(f, "%lf\n", duration);
  fclose(f);
  f = fopen((case_directory + "/vnr/new_cost").c_str(), "w");
  fprintf(f, "%lf\n", new_cost);
  fclose(f);
  f = fopen((case_directory + "/vnr/new_bnecks").c_str(), "w");
  fprintf(f, "%d\n", new_bnecks);
  fclose(f);
  f = fopen((case_directory + "/vnr/new_bw_cost").c_str(), "w");
  fprintf(f, "%ld\n", new_bw_cost);
  fclose(f);
  f = fopen((case_directory + "/vnr/new_max_plink_util").c_str(), "w");
  fprintf(f, "%lf\n", new_max_util);
  fclose(f);
  f = fopen((case_directory + "/vnr/max_iterations").c_str(), "w");
  fprintf(f, "%d\n", max_iterations);
  fclose(f);
  f = fopen((case_directory + "/vnr/iterations_per_temperature").c_str(), "w");
  fprintf(f, "%d\n", iterations_per_temperature);
  fclose(f);
  printf("Initial Solution Cost = %lf\n", initial_solutions[0].cost);
  printf("Best Solution Cost = %lf\n", best_solution->cost);
  WriteSolutionToFile(best_solution->vn_embeddings, case_directory + "/vnr");
  return 0;
}
